"""One-call local stationary-GE workflow; imports do not initialize or solve."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

for _thread_var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = "1"

import numpy as np

from .storage import StoredResult, load_case, publish_latest, reserve_case, save_case

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_ROOT = ROOT / "output/model/local_solution"


def _write_rows(path: Path, rows) -> None:
    rows = list(rows or [])
    if not rows: return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def _read_rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _write_summary(case: Path, outcome: dict) -> None:
    report = Path(outcome["report_directory"])
    target_fit = outcome.get("target_fit") or (_read_rows(report / "target_fit.csv") if (report / "target_fit.csv").is_file() else [])
    parameters = outcome.get("parameters") or (_read_rows(report / "parameters.csv") if (report / "parameters.csv").is_file() else [])
    if len(target_fit) != 14 or len(parameters) != 31:
        raise RuntimeError("A completed GE must provide 14 target-fit rows and 31 parameter rows")
    _write_rows(case / "target_fit.csv", target_fit)
    _write_rows(case / "parameters.csv", parameters)
    closure = outcome.get("closure", {})
    if isinstance(closure, dict):
        closure_name = closure.get("name", closure.get("closure_mode",
                                                         closure.get("closure", "stationary GE")))
        status = closure.get("status", "converged" if closure.get("converged", True) else "not converged")
        closure_items = [(key, value) for key, value in closure.items()
                         if isinstance(value, (str, int, float, bool)) or value is None]
    else:
        closure_name, status, closure_items = str(closure), "converged", []
    lines = ["# Local stationary GE", "", f"Status: {status}",
             f"Price: {float(outcome['price']):.12g}", f"Closure: {closure_name}"]
    if outcome.get("report_directory"):
        report = Path(outcome["report_directory"]).resolve()
        try:
            report_link = report.relative_to(case.resolve()).as_posix()
        except ValueError as exc:
            raise RuntimeError("Native report directory must be inside the private case") from exc
        lines.append(f"Native report directory: [{report_link}]({report_link})")
    lines.append("Native diagnostic figures: [standard_diagnostics/](standard_diagnostics/)")
    for key, value in closure_items:
        if key in {"name", "closure", "status", "converged"}:
            continue
        lines.append(f"{key}: {value}")
    if isinstance(closure, dict):
        if "fixed_h0_population_scale" in closure:
            lines.append("Fixed-H0 interpretation: population scale = "
                         f"{closure['fixed_h0_population_scale']}")
        if "implied_H0_at_population_one" in closure:
            lines.append("Population-one interpretation: implied H0 = "
                         f"{closure['implied_H0_at_population_one']}")
    lines += ["", f"Target-fit rows: {len(target_fit)} ([target_fit.csv](target_fit.csv)).",
              f"Parameter rows: {len(parameters)} ([parameters.csv](parameters.csv))."]
    (case / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def _cached_plots(case: Path) -> None:
    """Run existing cached plotters only after serialization; never calls a solver."""
    import sys
    model_dir = ROOT / "code/model"
    if str(model_dir) not in sys.path:
        sys.path.insert(0, str(model_dir))
    from plot_model_policies import _plot_run
    import plot_model_aggregates as aggregate_plots
    result, _ = load_case(case)
    policy_paths = _plot_run(result, case)
    if len(policy_paths) != 8:
        raise RuntimeError(f"Expected 8 cached policy plots, found {len(policy_paths)}")
    prior_run_directory = aggregate_plots.RUN_DIRECTORY
    aggregate_plots.RUN_DIRECTORY = case
    try:
        aggregate_directory = aggregate_plots.main()
    finally:
        aggregate_plots.RUN_DIRECTORY = prior_run_directory
    if len(list(Path(aggregate_directory).glob("*.png"))) != 7:
        raise RuntimeError("Cached aggregate plotter did not produce all 7 figures")


def _write_explorer_assets(case: Path, result: StoredResult) -> None:
    """Export the native fields consumed by the saved-case browser explorer."""
    required = ("b_grid", "V", "c_pol", "hR_pol", "bp_pol", "c_pol_stay",
                "bp_pol_stay", "tenure_probs", "fert_probs", "fert2_probs",
                "g_beginning_distribution", "g", "g_stay_distribution", "type_values")
    arrays = {}
    for name in required:
        if not hasattr(result.solution, name):
            raise RuntimeError(f"Explorer export is missing native solution field {name}")
        value = np.asarray(getattr(result.solution, name))
        if value.dtype.hasobject:
            raise TypeError(f"Explorer cannot safely export object array {name}")
        arrays[name] = value
    arrays_path = case / "explorer_arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    P = result.P
    timing = str(getattr(P, "purchase_timing", "inherited_only"))
    if timing not in {"inherited_only", "transaction_inside"}:
        raise RuntimeError(f"Unsupported explorer purchase timing label: {timing}")
    spec = {"id": "latest", "label": "Local stationary GE · post-interest timing",
            "arrays": str(arrays_path.resolve()),
            "sha256": hashlib.sha256(arrays_path.read_bytes()).hexdigest(),
            "price": float(result.price), "timing": "inherited_only"}
    config = {"common": {"age_start": int(P.age_start),
                         "period_years": int(P.period_years),
                         "houses": np.asarray(P.H_own, dtype=float).reshape(-1).tolist(),
                         "R_gross": float(P.R_gross),
                         "selling_cost": float(P.psi)},
              "cases": [spec], "report_root": str(case.resolve())}
    (case / "explorer_cases.json").write_text(json.dumps(config, indent=2) + "\n")


def _expose_native_diagnostics(case: Path, report: Path) -> None:
    source = report / "standard_diagnostics"
    figures = list(source.glob("*.png")) if source.is_dir() else []
    if len(figures) != 17:
        raise RuntimeError(f"Expected 17 native standard diagnostics, found {len(figures)}")
    destination = case / "standard_diagnostics"
    if destination.exists():
        raise FileExistsError(f"Refusing to replace existing diagnostics directory: {destination}")
    shutil.copytree(source, destination)


def run_stationary(parameters, external_inputs, native_overrides=None, price_guess=None,
                   budget_seconds=1800, closure="fixed_h0", output_root=None,
                   parameter_file_metadata=None):
    """Solve once, validate/cache its result, then atomically publish ``latest``.

    Core equilibrium code owns all numerical work.  This wrapper deliberately
    has no fallback solve, calibration loop, or economic defaults.
    """
    from .inputs import _base, load_inputs
    from .equilibrium import solve_stationary_ge
    import numba
    numba.set_num_threads(1)

    root = Path(output_root or DEFAULT_OUTPUT_ROOT)
    case = reserve_case(root)
    try:
        external_inputs = dict(external_inputs or {})
        native_overrides = dict(native_overrides or {})
        contract_path = case / "input_contract.json"
        input_contract = {
            "parameters": _jsonable(dict(parameters or {})),
            "external_inputs": _jsonable(external_inputs),
            "native_overrides": _jsonable(native_overrides),
            "price_guess": _jsonable(price_guess),
            "closure": _jsonable(closure),
            "budget_seconds": float(budget_seconds),
            "entry_law": "reference entrant-wealth mapping retained; no entry-law input was edited",
            "fiscal_mapping": {
                "status": "reference earnings/payroll primitives retained; no fiscal mapping change",
                "supplied_primitives": [],
                "edited_primitives": [],
            },
        }
        if parameter_file_metadata is not None:
            input_contract["parameter_file"] = _jsonable(parameter_file_metadata)
        contract_path.write_text(json.dumps(input_contract, indent=2, sort_keys=True,
                                            allow_nan=False) + "\n")
        P, grid = load_inputs(parameters=parameters, external_inputs=external_inputs,
                              native_overrides=native_overrides)
        receipt = getattr(P, "_production_input_receipt", None)
        fiscal_primitives = {"w_hat", "income_age_profile", "tau_pay"}
        fiscal_values = {**external_inputs, **native_overrides}
        supplied_fiscal = sorted(fiscal_primitives & set(fiscal_values))
        base_P, _ = _base()
        edited_fiscal = sorted(key for key in supplied_fiscal
                               if not np.array_equal(np.asarray(fiscal_values[key]),
                                                     np.asarray(getattr(base_P, key))))
        if receipt is not None:
            input_contract["fiscal_mapping"] = _jsonable(receipt)
        elif edited_fiscal:
            input_contract["fiscal_mapping"] = {
                "status": "native fixed-payroll balanced-pension mapping applied",
                "mapping": "native mapping derives disposable income and pension from gross wage, age profile, and payroll tax",
                "supplied_primitives": supplied_fiscal,
                "edited_primitives": edited_fiscal,
            }
        else:
            input_contract["fiscal_mapping"] = {
                "status": "reference earnings/payroll primitives retained; no fiscal mapping change",
                "supplied_primitives": supplied_fiscal,
                "edited_primitives": [],
            }
        contract_path.write_text(json.dumps(input_contract, indent=2, sort_keys=True,
                                            allow_nan=False) + "\n")
        # The core owns its report directory and requires a fresh output path.
        outcome = solve_stationary_ge(P, grid, out=case / "native", price_start=price_guess,
                                      budget_seconds=budget_seconds, closure=closure)
        required = {"solution", "P", "b_grid", "price", "closure", "report_directory"}
        missing = required - set(outcome)
        if missing: raise RuntimeError(f"Production GE return is missing: {sorted(missing)}")
        if isinstance(outcome["closure"], dict):
            outcome["closure"] = dict(outcome["closure"])
            outcome["closure"].setdefault("status", "converged")
        result = StoredResult(outcome["solution"], outcome["P"], outcome["b_grid"], float(outcome["price"]),
                              parameters=parameters, label="local stationary GE")
        result.closure = outcome["closure"]
        result.report_directory = str(outcome["report_directory"])
        _write_summary(case, outcome)
        _write_explorer_assets(case, result)
        save_case(result, case, metadata={"closure": outcome["closure"],
                                          "report_directory": str(outcome["report_directory"]),
                                          "budget_seconds": budget_seconds,
                                          "price_guess": _jsonable(price_guess),
                                          "input_contract_file": "input_contract.json",
                                          "input_contract_sha256": hashlib.sha256(
                                              contract_path.read_bytes()).hexdigest()})
        report = Path(outcome["report_directory"])
        native_plots = list(report.rglob("*.png")) if report.is_dir() else []
        if len(native_plots) != 17:
            raise RuntimeError(f"Expected 17 native diagnostic PNGs, found {len(native_plots)}")
        _expose_native_diagnostics(case, report)
        _cached_plots(case)
        # Reopen once more immediately before publication to catch any damaged
        # archive or metadata written by a later artifact step.
        reopened, _ = load_case(case)
        from .storage import _assert_roundtrip
        _assert_roundtrip(result, reopened)
        publish_latest(case, root)
        return result, case
    except Exception as exc:
        # Private attempts remain inspectable, while the previous latest is untouched.
        (case / "failure.json").write_text(json.dumps({"status": "failed",
                                                       "exception_type": type(exc).__name__,
                                                       "message": str(exc)},
                                                      indent=2) + "\n")
        raise
