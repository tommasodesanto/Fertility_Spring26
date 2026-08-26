#!/usr/bin/env python3
"""Audit horizon stability of coherent-person-law policy transitions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


MARKET_TOLERANCE = 2.0e-4
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5
DEFAULT_LEVEL_RELATIVE_TOLERANCE = 1.0e-3
DEFAULT_EFFECT_ABSOLUTE_TOLERANCE = 0.02


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--reform-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--comparison-periods", type=int, default=20)
    parser.add_argument(
        "--level-relative-tolerance",
        type=float,
        default=DEFAULT_LEVEL_RELATIVE_TOLERANCE,
    )
    parser.add_argument(
        "--effect-absolute-tolerance",
        type=float,
        default=DEFAULT_EFFECT_ABSOLUTE_TOLERANCE,
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_run(directory: Path, expected_case: str) -> dict[str, Any]:
    root = directory.resolve()
    summary_path = root / "summary.json"
    path_file = root / "transition_path.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
        raise RuntimeError(f"Incomplete path packet: {root}")
    if summary.get("case") != expected_case:
        raise RuntimeError(f"Wrong case in {root}")
    accounting = dict(summary.get("accounting_gates") or {})
    if not accounting or not all(accounting.values()):
        raise RuntimeError(f"Accounting gates failed in {root}")
    rows = read_csv(path_file)
    horizon = int(summary["horizon"])
    if len(rows) != horizon:
        raise RuntimeError(f"Path row count differs from horizon in {root}")
    return {
        "case": expected_case,
        "directory": root,
        "summary": summary,
        "summary_path": summary_path,
        "summary_sha256": sha256(summary_path),
        "path_file": path_file,
        "path_sha256": sha256(path_file),
        "horizon": horizon,
        "rows": rows,
    }


def path_gate_audit(run: dict[str, Any]) -> dict[str, Any]:
    summary = run["summary"]
    rows = run["rows"]
    market = np.asarray(
        [abs(float(row["relative_market_residual"])) for row in rows], dtype=float
    )
    fiscal = np.asarray(
        [abs(float(row["government_budget_residual"])) for row in rows],
        dtype=float,
    )
    finite = bool(np.all(np.isfinite(market)) and np.all(np.isfinite(fiscal)))
    maximum_market = float(np.max(market)) if finite else float("inf")
    maximum_fiscal = float(np.max(fiscal)) if finite else float("inf")
    market_pass = finite and maximum_market <= MARKET_TOLERANCE
    fiscal_pass = finite and maximum_fiscal <= FISCAL_ABSOLUTE_TOLERANCE
    recomputed_converged = market_pass and fiscal_pass
    summary_matches = bool(
        np.isclose(
            float(summary["maximum_market_residual"]),
            maximum_market,
            rtol=1e-12,
            atol=1e-14,
        )
        and np.isclose(
            float(summary["maximum_fiscal_residual"]),
            maximum_fiscal,
            rtol=1e-12,
            atol=1e-14,
        )
    )
    flag_consistent = bool(summary.get("path_converged")) == recomputed_converged
    return {
        "case": run["case"],
        "horizon": run["horizon"],
        "all_path_values_finite": finite,
        "maximum_market_residual_recomputed": maximum_market,
        "maximum_fiscal_residual_recomputed": maximum_fiscal,
        "market_tolerance": MARKET_TOLERANCE,
        "fiscal_absolute_tolerance": FISCAL_ABSOLUTE_TOLERANCE,
        "market_pass": market_pass,
        "fiscal_pass": fiscal_pass,
        "summary_matches_saved_path": summary_matches,
        "summary_path_flag_consistent": flag_consistent,
        "path_checks_pass": bool(
            recomputed_converged and summary_matches and flag_consistent
        ),
        "terminal_checks_pass": bool(
            (summary.get("terminal_convergence") or {}).get("all_checks_pass")
        ),
        "iterations": int(summary["iterations"]),
        "elapsed_seconds": float(summary["elapsed_seconds"]),
        "summary_path": str(run["summary_path"]),
        "summary_sha256": run["summary_sha256"],
        "transition_path": str(run["path_file"]),
        "transition_path_sha256": run["path_sha256"],
    }


def annual_birth_rate(row: dict[str, str]) -> float:
    return float(row["birth_children_topcode_adjusted"]) / max(
        4.0 * float(row["household_heads"]), 1e-15
    )


LEVEL_VARIABLES: tuple[tuple[str, Callable[[dict[str, str]], float]], ...] = (
    ("asset_price", lambda row: float(row["asset_price"])),
    ("renter_price", lambda row: float(row["renter_price"])),
    (
        "equal_transfer",
        lambda row: float(row["equal_transfer_period_units"]),
    ),
    ("annual_births_per_head", annual_birth_rate),
    (
        "resident_persons_actual_scale",
        lambda row: float(row["resident_persons_actual_scale"]),
    ),
    ("household_heads", lambda row: float(row["household_heads"])),
    ("owner_rate", lambda row: float(row["owner_rate"])),
)


def percent_change(reform: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    if np.any(np.abs(baseline) <= 1e-15):
        raise RuntimeError("Cannot form a policy percent change around zero")
    return 100.0 * (reform / baseline - 1.0)


def policy_effects(
    baseline: dict[str, Any], reform: dict[str, Any], periods: int
) -> dict[str, np.ndarray]:
    left_rows = baseline["rows"][:periods]
    right_rows = reform["rows"][:periods]
    left_years = np.asarray([int(row["calendar_year"]) for row in left_rows])
    right_years = np.asarray([int(row["calendar_year"]) for row in right_rows])
    if not np.array_equal(left_years, right_years):
        raise RuntimeError("Baseline and reform calendars differ")
    arrays = {}
    for name, getter in LEVEL_VARIABLES:
        left = np.asarray([getter(row) for row in left_rows], dtype=float)
        right = np.asarray([getter(row) for row in right_rows], dtype=float)
        if name == "owner_rate":
            arrays["owner_rate_percentage_point_change"] = 100.0 * (right - left)
        else:
            arrays[f"{name}_percent_change"] = percent_change(right, left)
    arrays["calendar_year"] = left_years
    return arrays


def contract_gate(runs: Sequence[dict[str, Any]]) -> dict[str, Any]:
    reference_hashes = runs[0]["summary"]["source_hashes"]
    checks: dict[str, bool] = {}
    for run in runs:
        checks[
            f"{run['case']}_h{run['horizon']}_source_and_input_fingerprint"
        ] = run["summary"]["source_hashes"] == reference_hashes
    for case in sorted({str(run["case"]) for run in runs}):
        case_runs = [run for run in runs if run["case"] == case]
        reference_terminal = case_runs[0]["summary"]["terminal_root"]
        for run in case_runs:
            checks[f"{case}_h{run['horizon']}_terminal_root"] = (
                run["summary"]["terminal_root"] == reference_terminal
            )
    by_horizon: dict[int, list[dict[str, Any]]] = {}
    for run in runs:
        by_horizon.setdefault(int(run["horizon"]), []).append(run)
    for horizon, pair in by_horizon.items():
        if len(pair) != 2:
            raise RuntimeError(f"Horizon {horizon} does not have both policy cases")
        initial_people = [float(run["rows"][0]["resident_persons"]) for run in pair]
        initial_heads = [float(run["rows"][0]["household_heads"]) for run in pair]
        checks[f"h{horizon}_common_initial_persons"] = bool(
            np.isclose(initial_people[0], initial_people[1], rtol=0.0, atol=2e-12)
        )
        checks[f"h{horizon}_common_initial_heads"] = bool(
            np.isclose(initial_heads[0], initial_heads[1], rtol=0.0, atol=2e-12)
        )
    return {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "source_hashes": reference_hashes,
    }


def level_stability(
    case_runs: Sequence[dict[str, Any]], periods: int, tolerance: float
) -> dict[str, Any]:
    ordered = sorted(case_runs, key=lambda run: int(run["horizon"]))
    reference = ordered[-1]
    records: list[dict[str, Any]] = []
    for candidate in ordered[:-1]:
        common = min(periods, candidate["horizon"], reference["horizon"])
        candidate_years = [
            int(row["calendar_year"]) for row in candidate["rows"][:common]
        ]
        reference_years = [
            int(row["calendar_year"]) for row in reference["rows"][:common]
        ]
        if candidate_years != reference_years:
            raise RuntimeError("Horizon paths do not share calendar dates")
        for name, getter in LEVEL_VARIABLES:
            left = np.asarray(
                [getter(row) for row in candidate["rows"][:common]], dtype=float
            )
            right = np.asarray(
                [getter(row) for row in reference["rows"][:common]], dtype=float
            )
            scale = np.maximum(np.abs(right), 1e-12)
            records.append(
                {
                    "comparison": candidate["case"],
                    "short_horizon": candidate["horizon"],
                    "reference_horizon": reference["horizon"],
                    "comparison_periods": common,
                    "variable": name,
                    "maximum_absolute_gap": float(np.max(np.abs(left - right))),
                    "maximum_relative_gap": float(
                        np.max(np.abs(left - right) / scale)
                    ),
                }
            )
    maximum = max(float(row["maximum_relative_gap"]) for row in records)
    return {
        "status": "passed" if maximum <= tolerance else "failed",
        "relative_tolerance": tolerance,
        "maximum_relative_gap": maximum,
        "rows": records,
    }


def effect_stability(
    baseline_runs: Sequence[dict[str, Any]],
    reform_runs: Sequence[dict[str, Any]],
    periods: int,
    tolerance: float,
) -> dict[str, Any]:
    baseline_by_horizon = {run["horizon"]: run for run in baseline_runs}
    reform_by_horizon = {run["horizon"]: run for run in reform_runs}
    horizons = sorted(baseline_by_horizon)
    reference_horizon = horizons[-1]
    reference_periods = min(periods, reference_horizon)
    reference = policy_effects(
        baseline_by_horizon[reference_horizon],
        reform_by_horizon[reference_horizon],
        reference_periods,
    )
    records: list[dict[str, Any]] = []
    for horizon in horizons[:-1]:
        common = min(periods, horizon, reference_horizon)
        candidate = policy_effects(
            baseline_by_horizon[horizon], reform_by_horizon[horizon], common
        )
        if not np.array_equal(
            candidate["calendar_year"], reference["calendar_year"][:common]
        ):
            raise RuntimeError("Policy-effect calendars differ across horizons")
        for name in sorted(key for key in candidate if key != "calendar_year"):
            gap = np.abs(candidate[name] - reference[name][:common])
            records.append(
                {
                    "comparison": "policy_effect",
                    "short_horizon": horizon,
                    "reference_horizon": reference_horizon,
                    "comparison_periods": common,
                    "variable": name,
                    "maximum_absolute_gap": float(np.max(gap)),
                }
            )
    maximum = max(float(row["maximum_absolute_gap"]) for row in records)
    return {
        "status": "passed" if maximum <= tolerance else "failed",
        "absolute_tolerance_percentage_points": tolerance,
        "maximum_absolute_gap": maximum,
        "rows": records,
    }


def make_figure(
    baseline_runs: Sequence[dict[str, Any]],
    reform_runs: Sequence[dict[str, Any]],
    periods: int,
    output_dir: Path,
) -> None:
    figure, axes = plt.subplots(4, 2, figsize=(11, 12), constrained_layout=True)
    panels = (
        ("asset_price_percent_change", "Asset-price effect"),
        ("renter_price_percent_change", "Rent-equivalent effect"),
        ("equal_transfer_percent_change", "Transfer effect"),
        ("annual_births_per_head_percent_change", "Annual-birth-rate effect"),
        (
            "resident_persons_actual_scale_percent_change",
            "Resident-person effect",
        ),
        ("household_heads_percent_change", "Household-head effect"),
        ("owner_rate_percentage_point_change", "Ownership effect"),
    )
    reform_by_horizon = {run["horizon"]: run for run in reform_runs}
    for baseline in sorted(baseline_runs, key=lambda run: run["horizon"]):
        horizon = baseline["horizon"]
        common = min(periods, horizon)
        effects = policy_effects(
            baseline, reform_by_horizon[horizon], common
        )
        for axis, (name, title) in zip(axes.flat, panels, strict=False):
            axis.plot(
                effects["calendar_year"], effects[name], label=f"H={horizon}"
            )
            axis.set_title(title)
            axis.set_xlabel("Year")
            axis.set_ylabel("Percentage points")
            axis.axhline(0.0, color="0.5", linewidth=1)
            axis.grid(alpha=0.2)
    for axis in axes.flat[: len(panels)]:
        axis.legend(frameon=False)
    axes.flat[-1].axis("off")
    figure.savefig(output_dir / "horizon_policy_effect_stability.png", dpi=200)
    figure.savefig(output_dir / "horizon_policy_effect_stability.pdf")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if int(args.comparison_periods) < 1:
        raise ValueError("--comparison-periods must be positive")
    if not 0.0 < float(args.level_relative_tolerance) < 1.0:
        raise ValueError("--level-relative-tolerance must lie in (0,1)")
    if not 0.0 < float(args.effect_absolute_tolerance) < 1.0:
        raise ValueError("--effect-absolute-tolerance must lie in (0,1)")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    baseline_runs = sorted(
        (
            load_run(path, "rebated-tax1-baseline")
            for path in args.baseline_dirs
        ),
        key=lambda run: run["horizon"],
    )
    reform_runs = sorted(
        (load_run(path, "rebated-tax2-reform") for path in args.reform_dirs),
        key=lambda run: run["horizon"],
    )
    baseline_horizons = [run["horizon"] for run in baseline_runs]
    reform_horizons = [run["horizon"] for run in reform_runs]
    if len(baseline_horizons) < 2:
        raise RuntimeError("Horizon stability requires at least two horizons")
    if baseline_horizons != reform_horizons:
        raise RuntimeError("Baseline and reform horizon sets differ")
    if len(set(baseline_horizons)) != len(baseline_horizons):
        raise RuntimeError("Duplicate horizons were supplied")

    runs = [*baseline_runs, *reform_runs]
    contract = contract_gate(runs)
    path_gates = [path_gate_audit(run) for run in runs]
    baseline_stability = level_stability(
        baseline_runs,
        int(args.comparison_periods),
        float(args.level_relative_tolerance),
    )
    reform_stability = level_stability(
        reform_runs,
        int(args.comparison_periods),
        float(args.level_relative_tolerance),
    )
    policy_stability = effect_stability(
        baseline_runs,
        reform_runs,
        int(args.comparison_periods),
        float(args.effect_absolute_tolerance),
    )
    longest_horizon = baseline_horizons[-1]
    longest_gates = [
        row for row in path_gates if int(row["horizon"]) == longest_horizon
    ]
    accepted = bool(
        contract["status"] == "passed"
        and all(bool(row["path_checks_pass"]) for row in path_gates)
        and all(bool(row["terminal_checks_pass"]) for row in longest_gates)
        and baseline_stability["status"] == "passed"
        and reform_stability["status"] == "passed"
        and policy_stability["status"] == "passed"
    )
    payload = {
        "status": (
            "certified_person_demography_policy_horizon_stability"
            if accepted
            else "failed_person_demography_policy_horizon_stability"
        ),
        "accepted": accepted,
        "promotion_status": "not_promoted",
        "interpretation": (
            "Certification requires all supplied horizons to pass independently "
            "recomputed market and fiscal gates, the longest horizon in each "
            "policy case to pass every terminal-state gate, common source and "
            "initial-state contracts, stable early level paths, and stable early "
            "policy effects."
        ),
        "horizons": baseline_horizons,
        "longest_horizon": longest_horizon,
        "comparison_periods": int(args.comparison_periods),
        "common_contract_gate": contract,
        "path_gates": path_gates,
        "baseline_level_stability": {
            key: value for key, value in baseline_stability.items() if key != "rows"
        },
        "reform_level_stability": {
            key: value for key, value in reform_stability.items() if key != "rows"
        },
        "policy_effect_stability": {
            key: value for key, value in policy_stability.items() if key != "rows"
        },
    }
    stability_rows = [
        *baseline_stability["rows"],
        *reform_stability["rows"],
        *policy_stability["rows"],
    ]
    write_json(output_dir / "certification.json", payload)
    write_csv(output_dir / "path_gates.csv", path_gates)
    write_csv(output_dir / "horizon_stability.csv", stability_rows)
    make_figure(
        baseline_runs,
        reform_runs,
        int(args.comparison_periods),
        output_dir,
    )


if __name__ == "__main__":
    main()
