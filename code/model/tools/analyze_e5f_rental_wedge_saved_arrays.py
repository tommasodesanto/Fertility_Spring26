"""Read-only diagnostics for retained rental-wedge policy arrays.

This reader deliberately imports no model modules and never solves or writes to
the experiment results directory.  The array axis contract is copied from the
factorial/cohort renderer: ``(wealth, tenure, location, age, income,
children_ever_born, child_state)``.  The renter state is tenure index 0.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import numpy as np


CASE_ORDER = (
    "cap6zero",
    "cap10zero",
    "cap10s005",
    "cap10s02",
    "cap10s1",
    "cap10s02_phi1",
)
AXES = (
    "wealth",
    "tenure",
    "location",
    "age",
    "income",
    "children_ever_born",
    "child_state",
)
TENURE_RENTER_INDEX = 0
POLICY_NAMES = (
    "V",
    "c_pol",
    "hR_pol",
    "bp_pol",
    "tenure_choice",
    "tenure_probs",
    "loc_probs",
    "fert_probs",
    "fert_value",
    "fert2_probs",
    "price",
)
STATE_ARRAY_NAMES = ("g_pre", "g_post_fertility", "g_current")
REQUIRED_NPZ_KEYS = POLICY_NAMES + STATE_ARRAY_NAMES + ("births",)
VALUE_PAIRS = (
    ("cap10zero", "cap10s005"),
    ("cap10s005", "cap10s02"),
    ("cap10s02", "cap10s1"),
    ("cap10s1", "cap6zero"),
    ("cap10s02_phi1", "cap10s02"),
)
VIOLATION_TOLERANCE = 1e-7
FINITE_VALUE_FLOOR = -1e9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def resolve_case_root(results: Path) -> tuple[Path, Path]:
    """Return ``(case_root, metadata_root)`` for either accepted CLI layout."""
    results = results.resolve()
    nested = results / "results"
    if nested.is_dir():
        return nested, results
    return results, results.parent


def assert_shape(name: str, array: np.ndarray, shape: tuple[int, ...]) -> None:
    if array.ndim != len(AXES) or array.shape != shape:
        raise ValueError(
            f"{name} shape {array.shape} does not match {AXES} contract {shape}"
        )


def check_case(case_root: Path, label: str) -> dict[str, Any]:
    case = case_root / label
    receipt_path = case / "receipt.json"
    array_path = case / "policy_arrays.npz"
    if not receipt_path.is_file() or not array_path.is_file():
        raise FileNotFoundError(f"missing retained case inputs for {label}")
    receipt = load_json(receipt_path)
    expected_hash = receipt.get("policy_arrays", {}).get("sha256")
    receipt_file_hash = receipt.get("file_hashes", {}).get("policy_arrays.npz")
    if not expected_hash or expected_hash != receipt_file_hash:
        raise ValueError(f"receipt has inconsistent policy-array hashes for {label}")
    actual_hash = sha256(array_path)
    if actual_hash != expected_hash:
        raise ValueError(f"policy-array SHA mismatch for {label}")

    with np.load(array_path, allow_pickle=False) as loaded:
        keys = tuple(sorted(loaded.files))
        missing = sorted(set(REQUIRED_NPZ_KEYS) - set(loaded.files))
        if missing:
            raise ValueError(f"{label} missing renderer keys: {missing}")
        # Copy only arrays needed after this case is closed; this keeps the
        # reader sequential and avoids retaining NPZ file handles.
        arrays = {
            name: np.array(loaded[name], copy=True)
            for name in ("V", "g_pre", "g_current", "hR_pol")
        }
        recorded_shapes = receipt["policy_arrays"].get("shapes", {})

    shape = tuple(int(x) for x in arrays["V"].shape)
    for name in ("V", "g_pre", "g_current", "hR_pol"):
        assert_shape(name, arrays[name], shape)
        if recorded_shapes and list(arrays[name].shape) != recorded_shapes.get(name):
            raise ValueError(f"receipt shape mismatch for {label}:{name}")
    for name in ("g_pre", "g_current"):
        if not np.isfinite(arrays[name]).all():
            raise ValueError(f"nonfinite mass in {label}:{name}")
        if float(arrays[name].min()) < -1e-12:
            raise ValueError(f"negative mass in {label}:{name}")

    renter = arrays["g_current"][:, TENURE_RENTER_INDEX, ...]
    renter_mass = float(renter.sum())
    reported_renter_mass = float(receipt["metrics"]["renter_mass"])
    if abs(renter_mass - reported_renter_mass) > 1e-12:
        raise ValueError(f"renter mass mismatch for {label}")
    rooms = arrays["hR_pol"][:, TENURE_RENTER_INDEX, ...]
    above_six = renter * (rooms > 6.0 + 1e-8)
    return {
        "label": label,
        "status": receipt.get("status"),
        "policy_arrays_sha256": actual_hash,
        "receipt_sha256": sha256(receipt_path),
        "array_keys": list(keys),
        "shape": list(shape),
        "axis_order": list(AXES),
        "g_pre_mass": float(arrays["g_pre"].sum()),
        "g_current_mass": float(arrays["g_current"].sum()),
        "renter_mass": renter_mass,
        "reported_renter_mass": reported_renter_mass,
        "renter_mass_abs_gap": abs(renter_mass - reported_renter_mass),
        "renter_hR_gt6": {
            "numerator_mass": float(above_six.sum()),
            "denominator_mass": renter_mass,
            "share": float(above_six.sum() / renter_mass),
        },
        "_arrays": arrays,
    }


def value_pair(
    higher: str,
    lower: str,
    values: dict[str, np.ndarray],
    g_pre: np.ndarray,
    joint_mask: np.ndarray,
) -> dict[str, Any]:
    # A positive difference is a violation of V(higher) >= V(lower).
    difference = values[lower] - values[higher]
    candidate = np.where(joint_mask, difference, -np.inf)
    flat_index = int(np.argmax(candidate))
    state = [int(x) for x in np.unravel_index(flat_index, candidate.shape)]
    raw_max = float(candidate[tuple(state)])
    max_violation = max(0.0, raw_max)
    v_higher = float(values[higher][tuple(state)])
    v_lower = float(values[lower][tuple(state)])
    relative_gap = float((v_lower - v_higher) / max(abs(v_higher), 1e-300))
    over = joint_mask & (difference > VIOLATION_TOLERANCE)
    occupied = over & (g_pre > 1e-12)
    positive_mass = over & (g_pre > 0.0)
    positive_difference = np.where(joint_mask, np.maximum(difference, 0.0), 0.0)
    return {
        "higher_case": higher,
        "lower_case": lower,
        "max_violation": max_violation,
        "raw_max_difference": raw_max,
        "max_state_index": state,
        "max_state_g_pre_mass": float(g_pre[tuple(state)]),
        "max_state_values": {"higher": v_higher, "lower": v_lower},
        "max_state_relative_gap": relative_gap,
        "count_violation_gt_1e-7": int(over.sum()),
        "g_pre_mass_on_violation_gt_1e-7": float(g_pre[over].sum()),
        "count_violation_gt_1e-7_gpre_positive": int(positive_mass.sum()),
        "count_violation_gt_1e-7_gpre_gt_1e-12": int(occupied.sum()),
        "g_pre_mass_on_violation_gt_1e-7_gpre_gt_1e-12": float(g_pre[occupied].sum()),
        "occupied_weighted_positive_violation_sum": float(
            positive_difference[occupied].dot(g_pre[occupied])
        ),
    }


def analyze(results: Path, output: Path) -> dict[str, Any]:
    case_root, metadata_root = resolve_case_root(results)
    output = output.resolve()
    if output == case_root or case_root in output.parents:
        raise ValueError("output must be separate from the immutable results root")
    output.mkdir(parents=True, exist_ok=True)

    cases: dict[str, dict[str, Any]] = {}
    values: dict[str, np.ndarray] = {}
    common_g_pre: np.ndarray | None = None
    source_identity: dict[str, Any] = {"cases": {}}
    for label in CASE_ORDER:
        checked = check_case(case_root, label)
        arrays = checked.pop("_arrays")
        values[label] = arrays["V"]
        case_g_pre = arrays["g_pre"]
        if common_g_pre is None:
            common_g_pre = case_g_pre
            source_identity["reference_case"] = label
            source_identity["reference_mass"] = float(common_g_pre.sum())
            source_identity["cases"][label] = True
        else:
            same = bool(np.array_equal(common_g_pre, case_g_pre))
            source_identity["cases"][label] = same
            if not same:
                raise ValueError(f"g_pre is not bitwise identical for {label}")
        cases[label] = checked

    assert common_g_pre is not None
    if not all(source_identity["cases"].values()):
        raise ValueError("common g_pre source identity failed")
    finite_all = np.isfinite(common_g_pre) & (common_g_pre >= 0.0)
    for value in values.values():
        finite_all &= np.isfinite(value) & (value > FINITE_VALUE_FLOOR)
    if not finite_all.any():
        raise ValueError("empty common finite value mask")
    pairs = {
        f"{higher}_ge_{lower}": value_pair(
            higher, lower, values, common_g_pre, finite_all
        )
        for higher, lower in VALUE_PAIRS
    }

    source_files = {}
    for name in ("source_manifest.json", "launch_manifest.json", "driver_plan.json"):
        path = metadata_root / name
        if path.is_file():
            source_files[name] = {"path": str(path), "sha256": sha256(path)}
    reader_path = Path(__file__).resolve()
    result = {
        "status": "complete",
        "scope": "read-only retained saved-array diagnostics; no household solve, replay, or plot regeneration",
        "results_root": str(case_root),
        "output_root": str(output),
        "reader": {"path": str(reader_path), "sha256": sha256(reader_path)},
        "source_hashes": source_files,
        "definitions": {
            "axis_order": list(AXES),
            "renter_tenure_index": TENURE_RENTER_INDEX,
            "renter_rooms_cutoff": "> 6 + 1e-8",
            "finite_value_domain": f"all six V arrays finite and > {FINITE_VALUE_FLOOR:g}",
            "joint_mask": "common g_pre finite/nonnegative AND all six V arrays finite and > -1e9",
            "violation": "V(lower) - V(higher)",
            "violation_tolerance": VIOLATION_TOLERANCE,
            "array_renderer_keys": list(POLICY_NAMES),
        },
        "source_identity": source_identity,
        "joint_mask_count": int(finite_all.sum()),
        "cases": cases,
        "value_choice_inequalities": pairs,
        "prior_axis_fix_audit": {
            "preserved_artifact": "results/saved_array_supplement_pre_axis_fix.json",
            "corrected_artifact": "results/saved_array_supplement.json",
            "correction": "previous [...,0,:] selected childless state; renter is [:,0,...] under the renderer axis contract",
        },
    }
    # Arrays are retained only in local variables and never serialized.
    (output / "lead_saved_array_review.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.results, args.output)
    print(json.dumps({"status": result["status"], "output": str(args.output.resolve()), "cases": len(result["cases"]), "joint_mask_count": result["joint_mask_count"]}))


if __name__ == "__main__":
    main()
