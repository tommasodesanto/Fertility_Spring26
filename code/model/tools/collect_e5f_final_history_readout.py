#!/usr/bin/env python3
"""Extract the final 2019--2023 readout from accepted native snapshots.

This collector performs no Bellman solve, path replay, root update, stationary
population construction, or plotting.  It configures the manifest-pinned
sequential runtime before unpickling, checks the accepted root and checkpoint
hashes, and then applies the retained measurement observers to the saved 2019
and 2023 evaluations.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import shutil
import sys
import time

import numpy as np


TOLERANCE = 2e-10


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value):
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def load_gzip_pickle(path):
    with gzip.open(path, "rb") as stream:
        return pickle.load(stream)


def resolve_path(value, base):
    path = Path(value)
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def normalized_pins(manifest, manifest_path):
    pins = manifest.get("file_sha256")
    if not isinstance(pins, dict) or not pins:
        raise ValueError("Manifest requires a nonempty explicit file_sha256 mapping")
    base = manifest_path.parent
    result = {}
    for raw_path, digest in pins.items():
        path = resolve_path(raw_path, base)
        if path in result or not isinstance(digest, str) or len(digest) != 64:
            raise ValueError("Manifest file pins must be unique SHA256 entries")
        result[path] = digest.lower()
    return result


def verify_pins(pins):
    for path, expected in pins.items():
        if not path.is_file():
            raise FileNotFoundError("Pinned input is missing: " + str(path))
        actual = file_sha256(path)
        if actual != expected:
            raise ValueError(f"Pinned input changed: {path}; {actual} != {expected}")


def validate_manifest(manifest, manifest_path):
    required = {"source_root", "accepted_forecast", "first_period_diagnostics",
        "native_2023_snapshot", "root_receipt", "kernel_files", "file_sha256"}
    missing = required - manifest.keys()
    if missing:
        raise ValueError("Manifest omits required fields: " + ", ".join(sorted(missing)))
    pins = normalized_pins(manifest, manifest_path)
    base = manifest_path.parent
    artifacts = {name: resolve_path(manifest[name], base) for name in (
        "accepted_forecast", "first_period_diagnostics", "native_2023_snapshot",
        "root_receipt")}
    kernels = [resolve_path(path, base) for path in manifest["kernel_files"]]
    if not kernels or len(kernels) != len(set(kernels)):
        raise ValueError("kernel_files must be a nonempty unique explicit list")
    for label, path in {**artifacts, **{f"kernel_{i}": p for i, p in enumerate(kernels)}}.items():
        if path not in pins:
            raise ValueError(f"{label} is not covered by the explicit file pins: {path}")
    source_root = resolve_path(manifest["source_root"], base)
    if not (source_root / "code/model/tools").is_dir() or not (source_root / "code/model").is_dir():
        raise ValueError("source_root does not contain the pinned model runtime")
    return pins, artifacts, kernels, source_root


def import_explicit_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError("Cannot load explicit module: " + str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_loaded_kernel(pins, kernels, item, label):
    path = Path(item.__code__.co_filename if hasattr(item, "__code__")
                else item.__file__).resolve()
    if path not in kernels or path not in pins:
        raise ValueError(f"Loaded {label} is not an explicitly pinned kernel: {path}")


def configure_runtime(*, source_root, pins, kernels, manifest, manifest_path):
    """Load and configure the exact runtime before any scientific unpickle."""
    tools = source_root / "code/model/tools"
    model = source_root / "code/model"
    runtime_paths = [resolve_path(path, manifest_path.parent)
        for path in manifest.get("runtime_paths", [manifest_path.parent, tools, model])]
    if not runtime_paths or any(not path.is_dir() for path in runtime_paths):
        raise ValueError("Explicit runtime search paths must be existing directories")
    sys.path[:0] = [str(path) for path in runtime_paths]
    profile_path = resolve_path(
        manifest.get("profile_kernel", tools / "collect_e5f_patch_readout.py"),
        manifest_path.parent)
    if profile_path not in kernels or profile_path not in pins:
        raise ValueError("The retained profile helper must be an explicitly pinned kernel")
    saved_path = list(sys.path)
    profile_module = import_explicit_module("_pinned_e5f_patch_profile", profile_path)
    sys.path[:] = saved_path
    import run_e5f_matched_pf_smoke as primitive
    import run_e5f_transition_calibration as fertility
    import e5f_rebated_surprises as rebated
    import e5f_closed_finite_boundary as closed_boundary
    import run_e5f_final_rebated_history as final_driver
    import e5f_initial_fertility_observer as fertility_stock
    import e5f_initial_housing_observer as housing_wealth
    import e5f_recent_parent_flow_observer as recent_parent
    primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = (
        primitive.pf.transition.advance_sequential_calendar_distribution)
    rebated._runtime()
    loaded = ((profile_module.profile, "profile"), (fertility, "fertility"),
        (rebated, "rebated forecast"), (closed_boundary, "finite boundary"),
        (final_driver, "final history driver"), (primitive, "matched PF runtime"),
        (fertility_stock, "fertility stock observer"),
        (housing_wealth, "housing/wealth observer"),
        (recent_parent, "recent-parent observer"))
    for item, label in loaded:
        require_loaded_kernel(pins, kernels, item, label)
    return dict(profile=profile_module.profile, fertility=fertility,
        fertility_stock=fertility_stock, housing_wealth=housing_wealth,
        recent_parent=recent_parent)


def validate_pack(pack, label):
    if not isinstance(pack, dict):
        raise ValueError(label + " must be a dictionary snapshot")
    required = {"parameters", "b_grid", "evaluation", "shared", "supply_rule"}
    missing = required - pack.keys()
    if missing:
        raise ValueError(label + " omits: " + ", ".join(sorted(missing)))
    grid = np.asarray(pack["b_grid"], dtype=float)
    evaluation = pack["evaluation"]
    if (grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all()
            or np.any(np.diff(grid) <= 0)):
        raise ValueError(label + " has an invalid wealth grid")
    for name in ("g_pre", "g_post_fertility", "g_current"):
        values = np.asarray(getattr(evaluation, name), dtype=float)
        if (values.ndim != 7 or values.shape[0] != len(grid)
                or not np.isfinite(values).all() or np.any(values < 0)):
            raise ValueError(f"{label} has invalid {name}")
    return pack


def row_for_year(rows, year):
    matches = [row for row in rows if int(row["calendar_year"]) == year]
    if len(matches) != 1:
        raise ValueError(f"Accepted path must contain exactly one {year} row")
    return matches[0]


def snapshot_aggregate_check(pack, row, profile_function, fertility_function):
    """Pure aggregate comparison of one native packet to its accepted row."""
    P, grid, evaluation = pack["parameters"], pack["b_grid"], pack["evaluation"]
    profile = profile_function(evaluation, P, grid)
    fertility = fertility_function(evaluation, P)
    totals = profile["totals"]
    heads = float(np.sum(evaluation.g_current))
    quantities = dict(
        housing_demand=float(evaluation.demand_by_loc[0]),
        housing_supply=float(evaluation.supply_by_loc[0]),
        owner_rate=float(np.sum(evaluation.g_current[:, 1:])) / heads,
        asset_price=float(np.asarray(evaluation.policy.price).reshape(-1)[0]),
        birth_children_topcode_adjusted=float(np.sum(
            fertility["birth_flow_topcode_adjusted"])),
        pension_period_units=float(P.pension),
        equal_transfer_period_units=float(P.property_tax_lump_sum_transfer))
    # The historical PF row calls current household mass ``adult_population``;
    # the person-demography row distinguishes it from resident persons and calls
    # it ``household_heads``.  Compare only the producer's date-specific field.
    head_key = "household_heads" if "household_heads" in row else "adult_population"
    if head_key not in row:
        raise ValueError("Accepted row omits its household-head aggregate")
    quantities[head_key] = heads
    gaps = {}
    for key, model_value in quantities.items():
        if key not in row:
            raise ValueError("Accepted row omits required snapshot aggregate: " + key)
        data_value = float(row[key])
        if not np.isfinite([model_value, data_value]).all():
            raise ValueError("Nonfinite snapshot comparison: " + key)
        gaps[key] = model_value - data_value
    gaps["profile_rooms_minus_evaluation_demand"] = (
        float(totals["rooms"]) - quantities["housing_demand"])
    gaps["profile_owners_minus_owner_rate"] = (
        float(totals["owners"]) / float(totals["households"])
        - quantities["owner_rate"])
    maximum = max(abs(value) for value in gaps.values())
    if maximum > TOLERANCE:
        raise RuntimeError(f"Native snapshot differs from accepted row by {maximum:.3e}")
    return dict(calendar_year=int(row["calendar_year"]), gaps=gaps,
                maximum_abs=maximum, profile=profile, fertility=fertility)


def validate_accepted(accepted, receipt):
    if not isinstance(accepted, dict) or not {"result", "boundary", "coordinates"} <= accepted.keys():
        raise ValueError("Accepted forecast checkpoint has the wrong schema")
    result, boundary = accepted["result"], accepted["boundary"]
    if result.path is None or result.next_state is None:
        raise ValueError("Accepted forecast lacks its path or realized next state")
    if clean(result.root_receipt) != receipt:
        raise ValueError("Accepted checkpoint root receipt differs from root_receipt.json")
    if (not receipt.get("converged")
            or not receipt.get("finite_horizon_market_fiscal_converged")
            or receipt.get("final") is None
            or not receipt["final"].get("mapping_valid")):
        raise ValueError("Root receipt does not certify a finite accepted mapping")
    reproduction = receipt.get("final_reproduction_max_abs")
    if reproduction is None or not np.isfinite(reproduction) or reproduction > TOLERANCE:
        raise ValueError("Root receipt lacks exact final reproduction")
    coordinates = np.asarray(accepted["coordinates"], dtype=float)
    final_coordinates = np.asarray(receipt["final"]["prices"], dtype=float)
    if (coordinates.shape != final_coordinates.shape or not np.isfinite(coordinates).all()
            or not np.array_equal(coordinates, final_coordinates)):
        raise ValueError("Accepted coordinates differ from the exact final root")
    required_boundary = {"parameters", "b_grid", "g_pre", "policy", "residuals",
        "actual_accounts", "diagnostics", "gates"}
    missing_boundary = [name for name in required_boundary if not hasattr(boundary, name)]
    if missing_boundary:
        raise ValueError("Accepted finite boundary omits: " + ", ".join(sorted(missing_boundary)))
    diagnostics = boundary.diagnostics
    if (hasattr(boundary, "fixed_point") or not isinstance(diagnostics, dict)
            or diagnostics.get("horizon_status") != "unverified_finite_truncation"
            or diagnostics.get("stationary_population_computed", False)
            or not isinstance(boundary.gates, dict) or not boundary.gates
            or not all(bool(value) for value in boundary.gates.values())
            or receipt.get("boundary_status") != "unverified_finite_truncation"
            or receipt.get("terminal_distance_passed") is not False
            or receipt.get("horizon_verified") is not False
            or receipt.get("production_eligible") is not False):
        raise ValueError("Accepted checkpoint is not the verified closed finite-boundary diagnostic")
    return result


def extract_measurements(pack2019, pack2023, runtime):
    fertility = runtime["fertility"]
    stock = runtime["fertility_stock"]
    housing = runtime["housing_wealth"]
    recent = runtime["recent_parent"]
    profile = runtime["profile"]
    errors = {}
    values = {}

    def measure(name, call):
        try:
            values[name] = call()
        except Exception as exc:  # every missing family must remain visible
            errors[name] = dict(error_type=type(exc).__name__, error=str(exc))

    e19, P19, grid19, shared19 = (pack2019[key] for key in (
        "evaluation", "parameters", "b_grid", "shared"))
    e23, P23, grid23, shared23 = (pack2023[key] for key in (
        "evaluation", "parameters", "b_grid", "shared"))
    branch = None
    try:
        branch = fertility.begin_dated_first_birth_housing_branch(
            e19, P19, grid19, shared19, origin_period=0)
    except Exception as exc:
        errors["first_birth_origin"] = dict(error_type=type(exc).__name__, error=str(exc))
    measure("profile", lambda: profile(e23, P23, grid23))
    measure("fertility", lambda: fertility.period_fertility_diagnostics(e23, P23))
    measure("fertility_stock_timing", lambda: stock.observe_initial_fertility(
        e23, P23, age_projection="uniform_birth_time"))
    measure("housing_wealth", lambda: housing.observe_initial_housing_wealth(
        e23, P23, grid23, shared23, diagnostic_enabled=True,
        age_projection="uniform_within_age_cell",
        diagnostic_allow_family_proxies=True, include_wealth=True,
        include_birth_response=False))
    measure("recent_parent", lambda: recent.observe_recent_parent_flow(
        e23, P23, diagnostic_enabled=True, snapshot=recent.SNAPSHOT,
        age_projection=recent.AGE_PROJECTION,
        diagnostic_allow_residence_proxy=True))
    if branch is not None:
        measure("dated_first_birth_rooms", lambda:
            fertility.finish_dated_first_birth_housing_branch(
                branch, e23, P23, grid23, shared23, destination_period=1))
    elif "first_birth_origin" not in errors:
        errors["first_birth_origin"] = dict(
            error_type="RuntimeError", error="Missing saved 2019 branch")
    required = {"profile", "fertility", "fertility_stock_timing",
        "housing_wealth", "recent_parent", "dated_first_birth_rooms"}
    for name in required - values.keys() - errors.keys():
        errors[name] = dict(error_type="RuntimeError", error="Observer did not return")
    return values, errors


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest_path = args.manifest.resolve()
    if file_sha256(manifest_path) != args.manifest_sha256.lower():
        raise ValueError("Manifest SHA256 mismatch")
    manifest = read_json(manifest_path)
    pins, artifacts, kernels, source_root = validate_manifest(manifest, manifest_path)
    verify_pins(pins)
    runtime = configure_runtime(source_root=source_root, pins=pins, kernels=kernels,
                                manifest=manifest, manifest_path=manifest_path)
    # Scientific objects are unpickled only after the sequential runtime above.
    accepted = load_gzip_pickle(artifacts["accepted_forecast"])
    pack2019 = validate_pack(load_gzip_pickle(
        artifacts["first_period_diagnostics"]), "2019 native snapshot")
    pack2023 = validate_pack(load_gzip_pickle(
        artifacts["native_2023_snapshot"]), "2023 native snapshot")
    receipt = read_json(artifacts["root_receipt"])
    result = validate_accepted(accepted, receipt)
    rows = list(result.path.rows)
    row2019, row2023 = row_for_year(rows, 2019), row_for_year(rows, 2023)
    check2019 = snapshot_aggregate_check(
        pack2019, row2019, runtime["profile"],
        runtime["fertility"].period_fertility_diagnostics)
    check2023 = snapshot_aggregate_check(
        pack2023, row2023, runtime["profile"],
        runtime["fertility"].period_fertility_diagnostics)
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(artifacts["root_receipt"], args.out / "root_receipt.json")
    started = time.monotonic()
    observed, errors = extract_measurements(pack2019, pack2023, runtime)
    source_hashes = {str(path): digest for path, digest in pins.items()}
    historical_status = copy.deepcopy(manifest.get("historical_fit_status", {
        "complete": False,
        "status": "not_independently_verified_by_final-window_extractor"}))
    horizon_status = dict(
        boundary_status=receipt.get("boundary_status"),
        terminal_distance_passed=bool(receipt.get("terminal_distance_passed", False)),
        horizon_verified=bool(receipt.get("horizon_verified", False)))
    common = dict(
        verification_method="native_saved_snapshot_aggregate_match",
        replay_performed=False, replay_maximum_abs=None,
        snapshot_maximum_abs=max(check2019["maximum_abs"], check2023["maximum_abs"]),
        dated_snapshot_checks=[check2019, check2023],
        finite_converged=True, historical_fit_status=historical_status,
        horizon_status=horizon_status, horizon_verified=horizon_status["horizon_verified"],
        production_eligible=False, stationary_reset_used=False,
        accepted_forecast_sha256=pins[artifacts["accepted_forecast"]],
        first_period_diagnostics_sha256=pins[artifacts["first_period_diagnostics"]],
        native_2023_snapshot_sha256=pins[artifacts["native_2023_snapshot"]],
        root_receipt_sha256=pins[artifacts["root_receipt"]],
        manifest_sha256=args.manifest_sha256.lower(), source_sha256=source_hashes)
    measurement = dict(errors=errors, stationary_restart=False,
        origin_snapshot=str(artifacts["first_period_diagnostics"]),
        destination_snapshot=str(artifacts["native_2023_snapshot"]),
        verification_method=common["verification_method"])
    write_json(args.out / "measurement_verification.json", measurement)
    if errors:
        write_json(args.out / "verification.json", dict(
            common, status="FAIL", seconds=time.monotonic() - started,
            error="One or more required saved-state measurements failed",
            measurement_errors=errors))
        raise RuntimeError("Required saved-state measurements failed: "
                           + ", ".join(sorted(errors)))
    model = dict(observed, calendar_year=2023,
        forecast_receipt_sha256=pins[artifacts["root_receipt"]],
        finite_converged=True, horizon_verified=horizon_status["horizon_verified"])
    write_json(args.out / "model_2023.json", model)
    verify_pins(pins)
    if file_sha256(manifest_path) != args.manifest_sha256.lower():
        raise ValueError("Readout manifest changed during extraction")
    verification = dict(common, status="PASS", year=2023,
        seconds=time.monotonic() - started, measurement_errors={})
    write_json(args.out / "verification.json", verification)
    print(json.dumps(dict(status="PASS", year=2023,
        verification_method=common["verification_method"],
        snapshot_maximum_abs=common["snapshot_maximum_abs"])))


if __name__ == "__main__":
    main()
