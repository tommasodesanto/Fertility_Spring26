"""Native equivalence check for the population-free finite-boundary policy seed.

This driver performs two bounded checks against a pinned, accepted A0 count-six
forecast.  It first compares the old full boundary evaluation with the new
policy-only constructor at the accepted terminal coordinates.  It then reruns
the accepted forecast from its final coordinates with exactly two root
evaluations and compares the complete dated result with the saved reference.

The script does not fit preferences, advance a historical window, or promote a
source tree.  Any unsupported object, changed pin, missing actual-population
audit, or nonexact numerical comparison is a hard failure.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pickle
import struct
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

for _thread_key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"
):
    os.environ[_thread_key] = "1"

import numpy as np


MAXIMUM_SECONDS = 20 * 60


def _sha(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path | str) -> Any:
    return json.loads(Path(path).read_text())


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_jsonable(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def _tagged_update(digest: Any, tag: str, payload: bytes = b"") -> None:
    encoded = tag.encode("utf-8")
    digest.update(struct.pack(">Q", len(encoded)))
    digest.update(encoded)
    digest.update(struct.pack(">Q", len(payload)))
    digest.update(payload)


def _digest_value(value: Any) -> str:
    """Hash supported numerical containers with dtype, shape, keys, and types."""
    digest = hashlib.sha256()
    active: set[int] = set()

    def visit(item: Any, location: str) -> None:
        if item is None:
            _tagged_update(digest, "none")
        elif isinstance(item, (bool, np.bool_)):
            _tagged_update(digest, "bool", b"1" if bool(item) else b"0")
        elif isinstance(item, np.ndarray):
            if item.dtype.hasobject:
                raise TypeError(f"Object-dtype array is unsupported at {location}")
            _tagged_update(digest, "ndarray-dtype", item.dtype.str.encode())
            _tagged_update(digest, "ndarray-shape", repr(item.shape).encode())
            _tagged_update(digest, "ndarray-bytes", np.ascontiguousarray(item).tobytes())
        elif isinstance(item, np.generic):
            _tagged_update(digest, "numpy-scalar-dtype", item.dtype.str.encode())
            _tagged_update(digest, "numpy-scalar-bytes", item.tobytes())
        elif isinstance(item, float):
            _tagged_update(digest, "float64", struct.pack(">d", item))
        elif isinstance(item, int):
            _tagged_update(digest, "int", str(item).encode())
        elif isinstance(item, str):
            _tagged_update(digest, "str", item.encode())
        elif isinstance(item, bytes):
            _tagged_update(digest, "bytes", item)
        elif isinstance(item, Path):
            _tagged_update(digest, "path", str(item).encode())
        elif isinstance(item, dict):
            identity = id(item)
            if identity in active:
                raise TypeError(f"Cyclic dictionary is unsupported at {location}")
            active.add(identity)
            _tagged_update(digest, "dict", str(len(item)).encode())
            for key in sorted(item, key=lambda candidate: (type(candidate).__name__, repr(candidate))):
                visit(key, location + ".<key>")
                visit(item[key], location + f"[{key!r}]")
            active.remove(identity)
        elif isinstance(item, (tuple, list)):
            identity = id(item)
            if identity in active:
                raise TypeError(f"Cyclic sequence is unsupported at {location}")
            active.add(identity)
            _tagged_update(digest, type(item).__name__, str(len(item)).encode())
            for index, child in enumerate(item):
                visit(child, location + f"[{index}]")
            active.remove(identity)
        elif hasattr(item, "__dict__"):
            identity = id(item)
            if identity in active:
                raise TypeError(f"Cyclic object is unsupported at {location}")
            active.add(identity)
            _tagged_update(
                digest, "object-type", f"{type(item).__module__}.{type(item).__qualname__}".encode()
            )
            visit(vars(item), location + ".__dict__")
            active.remove(identity)
        else:
            raise TypeError(f"Unsupported exact-hash type {type(item)!r} at {location}")

    visit(value, "root")
    return digest.hexdigest()


def _assert_exact(left: Any, right: Any, location: str = "root") -> None:
    """Require recursive numerical identity; tolerate NaN only at identical bits."""
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        if not isinstance(left, np.ndarray) or not isinstance(right, np.ndarray):
            raise AssertionError(f"Type mismatch at {location}: {type(left)} versus {type(right)}")
        if left.dtype != right.dtype or left.shape != right.shape:
            raise AssertionError(
                f"Array metadata mismatch at {location}: {left.dtype}/{left.shape} "
                f"versus {right.dtype}/{right.shape}"
            )
        if left.dtype.hasobject or left.tobytes() != right.tobytes():
            raise AssertionError(f"Array bytes differ at {location}")
        return
    if isinstance(left, np.generic) or isinstance(right, np.generic):
        if type(left) is not type(right) or left.tobytes() != right.tobytes():
            raise AssertionError(f"NumPy scalar differs at {location}")
        return
    if isinstance(left, dict) or isinstance(right, dict):
        if not isinstance(left, dict) or not isinstance(right, dict) or left.keys() != right.keys():
            raise AssertionError(f"Dictionary keys/types differ at {location}")
        for key in left:
            _assert_exact(left[key], right[key], location + f"[{key!r}]")
        return
    if isinstance(left, (tuple, list)) or isinstance(right, (tuple, list)):
        if type(left) is not type(right) or len(left) != len(right):
            raise AssertionError(f"Sequence type/length differs at {location}")
        for index, (a, b) in enumerate(zip(left, right)):
            _assert_exact(a, b, location + f"[{index}]")
        return
    if hasattr(left, "__dict__") or hasattr(right, "__dict__"):
        if not hasattr(left, "__dict__") or not hasattr(right, "__dict__"):
            raise AssertionError(f"Object/container mismatch at {location}")
        _assert_exact(vars(left), vars(right), location + ".__dict__")
        return
    if isinstance(left, float) or isinstance(right, float):
        if not isinstance(left, (float, np.floating)) or not isinstance(right, (float, np.floating)):
            raise AssertionError(f"Float type mismatch at {location}")
        if struct.pack(">d", float(left)) != struct.pack(">d", float(right)):
            raise AssertionError(f"Float bits differ at {location}: {left!r} versus {right!r}")
        return
    if type(left) is not type(right) or left != right:
        raise AssertionError(f"Values differ at {location}: {left!r} versus {right!r}")


def _field_hashes(value: Any) -> dict[str, str]:
    if not hasattr(value, "__dict__"):
        raise TypeError("Policy must expose native fields through __dict__")
    return {name: _digest_value(field) for name, field in sorted(vars(value).items())}


def _require_pinned(pins: dict[str, str], path: Path, label: str) -> None:
    resolved = path.resolve()
    matches = [(Path(candidate).resolve(), digest) for candidate, digest in pins.items()
               if Path(candidate).resolve() == resolved]
    if len(matches) != 1 or matches[0][1] != _sha(resolved):
        raise ValueError(f"{label} is not uniquely SHA256-pinned by the manifest: {resolved}")


def _load_old_boundary(path: Path) -> Any:
    name = "e5f_closed_finite_boundary_old_reference"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load old boundary helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    if not hasattr(module, "boundary_evaluation") or hasattr(module, "boundary_policy"):
        raise ValueError("Old helper must expose boundary_evaluation and predate boundary_policy")
    return module


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--old-helper", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--initial-raw-summary-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    started = time.monotonic()
    deadline = started + MAXIMUM_SECONDS
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summary: dict[str, Any] = {
        "status": "running",
        "schema": "e5f_boundary_policy_seed_native_equivalence_v1",
        "maximum_seconds": MAXIMUM_SECONDS,
        "checks": {},
    }
    _write(out / "summary.json", summary)

    try:
        manifest_path = args.manifest.resolve()
        helper = args.helper.resolve()
        old_helper = args.old_helper.resolve()
        case_dir = args.case_dir.resolve()
        for directory, label in ((helper, "new helper"), (old_helper, "old helper"),
                                 (case_dir, "reference case")):
            if not directory.is_dir():
                raise ValueError(f"Missing {label} directory: {directory}")

        manifest = _read(manifest_path)
        pins = manifest.get("file_sha256")
        if not isinstance(pins, dict) or not pins:
            raise ValueError("New manifest needs a nonempty file_sha256 contract")
        prior_plan_path = Path(manifest.get("prior_plan", "")).resolve()
        if str(prior_plan_path) not in pins:
            raise ValueError("New manifest must pin its complete prior plan")
        plan = _read(prior_plan_path)
        source_root = Path(plan["source_root"]).resolve()

        sys.path[:0] = [str(helper), str(source_root / "code/model/tools"),
                        str(source_root / "code/model")]
        import run_e5f_final_rebated_history as driver
        import e5f_closed_finite_boundary as new_boundary
        import e5f_rebated_surprises as rebated
        from e5f_balanced_terminal import TerminalAuditControls
        from e5f_rebated_initial_bridge import build_rebated_initial_state

        if Path(driver.__file__).resolve() != helper / "run_e5f_final_rebated_history.py":
            raise ImportError("New driver did not load from --helper")
        if Path(new_boundary.__file__).resolve() != helper / "e5f_closed_finite_boundary.py":
            raise ImportError("New boundary module did not load from --helper")
        driver.verify_pins(pins)
        driver.verify_pins(plan["file_sha256"])
        equivalence_path = Path(manifest['kernel_equivalence'])
        _require_pinned(pins, equivalence_path, 'initial/history kernel equivalence')
        pairs = _read(equivalence_path)['pairs']
        if len(pairs) != int(manifest['kernel_equivalence_pairs']):
            raise ValueError('Unexpected initial/history kernel pair count')
        for pair in pairs:
            if _sha(pair['initial']) != pair['sha256'] or _sha(pair['history']) != pair['sha256']:
                raise ValueError('Initial/history model kernel equivalence failed')
        for filename in ("run_e5f_final_rebated_history.py", "e5f_closed_finite_boundary.py"):
            _require_pinned(pins, helper / filename, "new " + filename)

        old_case_manifest_path = case_dir.parents[1] / "history_manifest_corrected_seeded_v1.json"
        old_case_manifest = _read(old_case_manifest_path)
        contract = _read(case_dir / "contract_receipt.json")
        if _sha(old_case_manifest_path) != contract.get("manifest_sha256"):
            raise ValueError("Reference case does not match its original manifest")
        driver.verify_pins(old_case_manifest["file_sha256"])
        old_case_plan_path = Path(old_case_manifest["prior_plan"])
        if str(old_case_plan_path) not in old_case_manifest["file_sha256"]:
            raise ValueError("Reference manifest does not pin its prior plan")
        old_case_plan = _read(old_case_plan_path)
        driver.verify_pins(old_case_plan["file_sha256"])
        if (contract.get("case") != "A0" or contract.get("count") != 6
                or case_dir.name != "A0_6" or contract.get("source_root") != str(source_root)):
            raise ValueError("Reference must be the accepted A0 count-six case on the same source root")

        old_boundary_path = old_helper / "e5f_closed_finite_boundary.py"
        old_boundary_digest = _sha(old_boundary_path)
        seeded_boundary_pins = [digest for path, digest in old_case_manifest["file_sha256"].items()
                                if Path(path).name == "e5f_closed_finite_boundary.py"]
        if old_boundary_digest not in seeded_boundary_pins:
            raise ValueError("Old boundary helper does not match the reference case's pinned boundary")
        old_boundary = _load_old_boundary(old_boundary_path)

        summary_path = Path(manifest["initial_summary"])
        if str(summary_path.resolve()) not in pins:
            raise ValueError("New manifest must pin the initial summary")
        initial_summary = _read(summary_path)
        if initial_summary.get("status") != "verified_rebated_initial_smoke":
            raise ValueError("Initial rebated smoke is not verified")
        checkpoint_entry = initial_summary["checkpoint"]
        checkpoint_path = Path(checkpoint_entry.get("path", checkpoint_entry.get("checkpoint")))
        checkpoint_digest = checkpoint_entry.get(
            "sha256", checkpoint_entry.get("checkpoint_sha256")
        )
        if (_sha(checkpoint_path) != checkpoint_digest
                or checkpoint_digest != contract.get("initial_checkpoint_sha256")):
            raise ValueError("New manifest and old accepted case do not share the pinned checkpoint")

        rebated._runtime()
        with gzip.open(checkpoint_path, "rb") as stream:
            packet = pickle.load(stream)
        raw_summary_path = Path(manifest.get(
            "initial_raw_summary", str(checkpoint_path.parent / "summary.json")
        ))
        if _sha(raw_summary_path) != args.initial_raw_summary_sha256:
            raise ValueError("Initial raw summary differs from its explicit launch pin")
        raw_summary = _read(raw_summary_path)
        old_state = build_rebated_initial_state(
            packet=packet, normalization=raw_summary["normalization"],
            outside_origin_entry_share=plan["outside_origin_entry_share"],
            preference_change_2023=0.0,
        )
        demographics = driver.migration_case(packet["demographic_seed"], "A0")
        inherited = rebated.InheritedState(2007, old_state.initial_state)
        audit = TerminalAuditControls(**plan["terminal_template"]["audit_controls"])
        runtime = rebated._runtime()

        reference_folder = case_dir / "window_2007/trial_00"
        fit_path = reference_folder / "fit.json"
        root_path = reference_folder / "root_receipt.json"
        accepted_path = reference_folder / "accepted_forecast.pkl.gz"
        fit = _read(fit_path)
        root_receipt = _read(root_path)
        realized = _read(case_dir / "realized_fit.json")
        if not realized or realized[0] != fit or Path(fit.get("folder", "")).resolve() != reference_folder:
            raise ValueError("Reference first-window fit does not match its realized-fit ledger")
        psi = float(fit["psi"])
        final = root_receipt.get("final")
        if (root_receipt.get("converged") is not True or final is None
                or final.get("mapping_valid") is not True or root_receipt.get("count") != 6
                or root_receipt.get("case") != "A0" or root_receipt.get("start_year") != 2007
                or float(root_receipt.get("psi")) != psi):
            raise ValueError("Reference first-window root is not an accepted A0 count-six result")
        coordinates = np.asarray(final["prices"], dtype=float)
        driver.unpack_coordinates(coordinates, 6)
        with gzip.open(accepted_path, "rb") as stream:
            accepted = pickle.load(stream)
        old_result = accepted.get("result")
        old_saved_boundary = accepted.get("boundary")
        if (old_result is None or old_saved_boundary is None
                or not np.array_equal(np.asarray(accepted.get("coordinates")), coordinates)
                or driver.clean(old_result.root_receipt) != root_receipt):
            raise ValueError("Reference accepted checkpoint differs from its root receipt")

        prices, pensions, transfers = driver.unpack_coordinates(coordinates, 6)
        parameters = copy.deepcopy(old_state.parameters)
        parameters.psi_child = psi
        boundary_arguments = dict(
            parameters=parameters, grid=old_state.b_grid,
            price=float(prices[-1]), pension=float(pensions[-1]),
            transfer=float(transfers[-1]), deadline_monotonic=deadline,
        )
        old_evaluation = old_boundary.boundary_evaluation(
            **boundary_arguments, g_pre=inherited.households.g_pre,
            supply_rule=old_state.supply_rule, audit_controls=audit,
        )
        new_policy = new_boundary.boundary_policy(**boundary_arguments)
        _assert_exact(vars(old_evaluation.policy), vars(new_policy.policy), "boundary.policy")
        fiscal_fields = ("pension", "property_tax_lump_sum_transfer", "tau_pay", "income")
        for name in fiscal_fields:
            if not hasattr(old_evaluation.parameters, name) or not hasattr(new_policy.parameters, name):
                raise AttributeError("Missing bound fiscal parameter: " + name)
            _assert_exact(
                getattr(old_evaluation.parameters, name), getattr(new_policy.parameters, name),
                "boundary.parameters." + name,
            )
        new_initial_audit = driver.cached_boundary(
            new_policy, inherited.households.g_pre, old_state, audit, runtime
        )
        # cached_boundary deliberately stores the compact root ledger.  Compare
        # every field it exposes, then verify that the old full evaluator has
        # exactly the documented reporting-only additions and that they are
        # internally exact.  A broader schema difference is a hard failure.
        account_extras = {"equal_transfer_period_units", "implied_equal_transfer_period"}
        residual_extras = {"housing_demand", "housing_supply", "housing_absolute",
                           "pension_absolute", "rebate_absolute"}
        if (set(old_evaluation.actual_accounts) - set(new_initial_audit.actual_accounts)
                != account_extras
                or set(new_initial_audit.actual_accounts) - set(old_evaluation.actual_accounts)):
            raise AssertionError("Unexpected old-full/new-cached account schema difference")
        if (set(old_evaluation.residuals) - set(new_initial_audit.residuals)
                != residual_extras
                or set(new_initial_audit.residuals) - set(old_evaluation.residuals)):
            raise AssertionError("Unexpected old-full/new-cached residual schema difference")
        for name, value in new_initial_audit.actual_accounts.items():
            _assert_exact(old_evaluation.actual_accounts[name], value,
                          "initial_boundary.actual_accounts." + name)
        for name, value in new_initial_audit.residuals.items():
            _assert_exact(old_evaluation.residuals[name], value,
                          "initial_boundary.residuals." + name)
        _assert_exact(old_evaluation.actual_accounts["equal_transfer_period_units"],
                      new_policy.parameters.property_tax_lump_sum_transfer,
                      "initial_boundary.equal_transfer_period_units")
        implied = (new_initial_audit.actual_accounts["property_tax_revenue"]
                   / new_initial_audit.actual_accounts["household_heads"])
        _assert_exact(old_evaluation.actual_accounts["implied_equal_transfer_period"], implied,
                      "initial_boundary.implied_equal_transfer_period")
        _assert_exact(old_evaluation.residuals["housing_absolute"],
                      old_evaluation.residuals["housing_demand"]
                      - old_evaluation.residuals["housing_supply"],
                      "initial_boundary.housing_absolute")
        _assert_exact(old_evaluation.residuals["pension_absolute"],
                      new_initial_audit.actual_accounts["pension_budget_residual"],
                      "initial_boundary.pension_absolute")
        _assert_exact(old_evaluation.residuals["rebate_absolute"],
                      new_initial_audit.actual_accounts["property_tax_budget_residual"],
                      "initial_boundary.rebate_absolute")
        _assert_exact(old_evaluation.gates, new_initial_audit.gates, "initial_boundary.gates")
        summary["checks"]["valid_boundary_policy"] = {
            "passed": True,
            "coordinates": [float(prices[-1]), float(pensions[-1]), float(transfers[-1])],
            "policy_sha256": _digest_value(vars(new_policy.policy)),
            "policy_field_sha256": _field_hashes(new_policy.policy),
            "bound_fiscal_sha256": {
                name: _digest_value(getattr(new_policy.parameters, name)) for name in fiscal_fields
            },
            "old_full_and_new_cached_account_fields_exact": True,
            "old_full_and_new_cached_residual_fields_exact": True,
            "documented_compact_schema_verified": True,
            "old_full_and_new_cached_gates_exact": True,
        }
        _write(out / "summary.json", summary)

        controls = dict(plan["history_root_controls"])
        controls.update(manifest.get("root_controls", {}))
        controls.setdefault("transfer_bounds", [1e-10, 10.0])
        controls["max_evaluations"] = 2
        controls["automatic_fiscal_polish"] = False

        native_boundary_policy = new_boundary.boundary_policy
        native_cached_boundary = driver.cached_boundary
        templates: list[dict[str, Any]] = []
        actual_audits: list[dict[str, Any]] = []

        def tracked_policy(**kwargs: Any) -> Any:
            template = native_boundary_policy(**kwargs)
            templates.append({
                "object": template,
                "before": _digest_value(vars(template.policy)),
                "fields_before": _field_hashes(template.policy),
            })
            return template

        def tracked_cached(template: Any, actual_g: np.ndarray, *call_args: Any) -> Any:
            match = next((entry for entry in reversed(templates) if entry["object"] is template), None)
            if match is None:
                raise RuntimeError("Actual boundary audit received an untracked policy template")
            before_actual_audit = _digest_value(vars(template.policy))
            if before_actual_audit != match["before"]:
                raise RuntimeError("Boundary policy mutated during the backward/forward forecast")
            value = native_cached_boundary(template, actual_g, *call_args)
            actual_audits.append({
                "terminal_g_sha256": _digest_value(np.asarray(actual_g)),
                "policy_sha256_before_actual_audit": before_actual_audit,
                "policy_sha256_after_actual_audit": _digest_value(vars(template.policy)),
                "gates": dict(value.gates),
            })
            return value

        new_boundary.boundary_policy = tracked_policy
        driver.cached_boundary = tracked_cached
        try:
            new_result, detail = driver.solve_forecast(
                inherited=inherited, old=old_state, demographics=demographics,
                psi=psi, count=6, initial=coordinates, controls=controls,
                audit=audit, deadline=deadline, folder=out / "solve_forecast", case="A0",
            )
        finally:
            new_boundary.boundary_policy = native_boundary_policy
            driver.cached_boundary = native_cached_boundary

        if time.monotonic() > deadline:
            raise TimeoutError("Native validation exceeded its 20-minute hard budget")
        new_root = new_result.root_receipt
        evaluations = int(new_root.get("evaluations", -1))
        if (new_result.next_state is None or new_root.get("converged") is not True
                or new_root.get("finite_horizon_market_fiscal_converged") is not True
                or evaluations != 2 or len(templates) != evaluations
                or len(actual_audits) != evaluations):
            raise RuntimeError("Two-map native replay did not execute one actual boundary audit per map")
        for index, (template, audit_record) in enumerate(zip(templates, actual_audits)):
            after = _digest_value(vars(template["object"].policy))
            if (after != template["before"]
                    or audit_record["policy_sha256_after_actual_audit"] != template["before"]
                    or not all(audit_record["gates"].values())):
                raise RuntimeError(f"Policy mutation or failed actual boundary gate in map {index + 1}")

        new_final = new_root.get("final")
        if new_final is None or new_final.get("mapping_valid") is not True:
            raise RuntimeError("Native replay lacks a valid final mapping")
        _assert_exact(np.asarray(new_final["prices"]), coordinates, "root.final.coordinates")
        _assert_exact(np.asarray(new_final["residual"]), np.asarray(final["residual"]),
                      "root.final.residual")
        _assert_exact(new_result.path.rows, old_result.path.rows, "path.rows")
        _assert_exact(new_result.path.values, old_result.path.values, "path.values")
        _assert_exact(new_result.path.person_tail.terminal_state.g_pre,
                      old_result.path.person_tail.terminal_state.g_pre,
                      "path.person_tail.terminal_state.g_pre")
        _assert_exact(detail["boundary"].actual_accounts, old_saved_boundary.actual_accounts,
                      "terminal_boundary.actual_accounts")
        _assert_exact(detail["boundary"].gates, old_saved_boundary.gates,
                      "terminal_boundary.gates")
        _assert_exact(detail["boundary"].residuals, old_saved_boundary.residuals,
                      "terminal_boundary.residuals")

        generated_checkpoint = out / "solve_forecast/accepted_forecast.pkl.gz"
        if not generated_checkpoint.is_file():
            raise RuntimeError("Native solve did not write its accepted forecast checkpoint")
        with gzip.open(generated_checkpoint, "rb") as stream:
            generated = pickle.load(stream)
        _assert_exact(np.asarray(generated["coordinates"]), coordinates,
                      "generated_checkpoint.coordinates")
        _assert_exact(generated["result"].path.rows, new_result.path.rows,
                      "generated_checkpoint.path.rows")

        input_paths = [
            manifest_path, prior_plan_path, summary_path, raw_summary_path, checkpoint_path,
            old_case_manifest_path, old_case_plan_path, case_dir / "contract_receipt.json",
            fit_path, root_path, accepted_path, old_boundary_path,
            helper / "e5f_closed_finite_boundary.py",
            helper / "run_e5f_final_rebated_history.py", Path(__file__).resolve(),
        ]
        summary["checks"]["accepted_count_six_replay"] = {
            "passed": True,
            "psi": psi,
            "root_evaluations": evaluations,
            "actual_boundary_audit_calls": len(actual_audits),
            "final_residual_sha256": _digest_value(np.asarray(new_final["residual"])),
            "path_rows_sha256": _digest_value(new_result.path.rows),
            "path_values_sha256": _digest_value(new_result.path.values),
            "terminal_g_pre_sha256": _digest_value(
                new_result.path.person_tail.terminal_state.g_pre
            ),
            "boundary_accounts_sha256": _digest_value(detail["boundary"].actual_accounts),
            "boundary_gates_sha256": _digest_value(detail["boundary"].gates),
            "policy_mutation_checks": [
                {
                    "before_sha256": item["before"],
                    "after_sha256": _digest_value(vars(item["object"].policy)),
                    "field_sha256": item["fields_before"],
                }
                for item in templates
            ],
            "actual_boundary_audits": actual_audits,
            "generated_checkpoint": {
                "path": str(generated_checkpoint), "sha256": _sha(generated_checkpoint)
            },
        }
        summary.update(
            status="passed",
            elapsed_seconds=time.monotonic() - started,
            source_root=str(source_root),
            case="A0",
            count=6,
            input_sha256={str(path): _sha(path) for path in input_paths},
            no_history_advanced=True,
            source_promoted=False,
        )
        _write(out / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"], "elapsed_seconds": summary["elapsed_seconds"],
            "checks": sorted(summary["checks"]),
        }), flush=True)
    except BaseException as error:
        summary.update(
            status="failed", elapsed_seconds=time.monotonic() - started,
            error_type=type(error).__name__, error=str(error),
        )
        _write(out / "summary.json", summary)
        raise


if __name__ == "__main__":
    main()
