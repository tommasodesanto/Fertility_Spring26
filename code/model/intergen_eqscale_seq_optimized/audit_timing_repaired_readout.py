#!/usr/bin/env python3
"""Re-measure collected E2/E3/E3b winners after the wealth-timing repair."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any

# Match the deterministic thread contract of run_e1_chain before NumPy/Numba
# heavy imports occur.
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = ROOT / "output/model/eqscale_timing_repair_readout_20260723"
CASES: dict[str, dict[str, Any]] = {
    "e2": {
        "record": ROOT / "output/model/eqscale_seq_optimized_recalibration_20260719/report/results.json",
        "l4": False,
    },
    "e3": {
        "record": ROOT / "output/model/eqscale_seq_l4_recalibration_20260722/report/results.json",
        "l4": True,
    },
    "e3b": {
        "record": ROOT / "output/model/eqscale_seq_l4_recalibration_20260722/report_continuation/results.json",
        "l4": True,
    },
}
WEALTH_ROWS = {
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535",
    "old_total_estate_wealth_to_annual_income_median_7684",
    "old_nonhousing_ge_1x_income_share_6575",
}
P90_P50 = "old_total_wealth_to_annual_income_p90_p50_7684"
P90_P50_LEGACY = "old_total_estate_wealth_to_annual_income_p90_p50_7684"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cases", default="e2,e3,e3b",
        help="Comma-separated subset of e2,e3,e3b (default: all).",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def selected_cases(value: str) -> list[str]:
    names = [name.strip().lower() for name in value.split(",") if name.strip()]
    unknown = sorted(set(names) - set(CASES))
    if unknown or not names:
        raise ValueError(f"--cases must select from {', '.join(CASES)}; got {value!r}")
    return list(dict.fromkeys(names))


def load_winner(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    try:
        winner = payload["winners"]["E1"]
    except (KeyError, TypeError) as exc:
        raise RuntimeError(f"winner record schema differs: {path} has no winners.E1") from exc
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise RuntimeError(f"winner record schema differs: {path} has no winners.E1.theta")
    if not isinstance(winner.get("moments"), dict):
        raise RuntimeError(f"winner record schema differs: {path} has no winners.E1.moments")
    return winner


def solve_twice(case: str, metadata: dict[str, Any], theta: dict[str, float]) -> tuple[dict[str, float], Any, bool]:
    """Return target moments and the first strict solution after an exact gate."""
    if metadata["l4"]:
        os.environ["E3_L4"] = "1"
    else:
        os.environ.pop("E3_L4", None)

    # This module intentionally delays NumPy/model imports.  Calling its loader
    # after the environment assignment makes common_overrides observe the arm.
    from intergen_eqscale_seq_optimized import run_e1_chain as chain
    from intergen_eqscale_seq_optimized.e1_profile import E1_TARGET_SET

    chain.load_runtime()
    args = argparse.Namespace(J=17, Nb=120, max_iter_eq=10, tol_eq=1e-4)
    overrides = chain.common_overrides(args)
    targets, weights = chain.get_target_set(E1_TARGET_SET)
    # This is part of the E-chain target contract, applied immediately after
    # fetching the atomic target registry in run_e1_chain.main().
    targets["aggregate_mean_occupied_rooms_18_85"] = 5.779970481941968
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    if set(targets) != set(weights):
        raise RuntimeError(f"{case}: target and weight keys differ")

    repeated: list[tuple[dict[str, float], Any]] = []
    for _ in range(2):
        sol, parameters, _ = chain.run_model_cp_dt(
            {**overrides, "max_iter_eq": 40, "tol_eq": 2.5e-5, **theta},
            verbose=False,
        )
        all_moments = chain.extract_moments(sol, parameters)
        moments = {name: float(all_moments[name]) for name in targets}
        repeated.append((moments, sol))

    first, second = repeated[0][0], repeated[1][0]
    differing = [name for name in targets if first[name] != second[name]]
    if differing:
        values = "; ".join(
            f"{name}: {first[name]!r} versus {second[name]!r}" for name in differing
        )
        raise RuntimeError(f"{case}: strict repeats are not bitwise-identical: {values}")
    return first, repeated[0][1], True


def nonwealth_difference(rows: list[dict[str, Any]]) -> float:
    eligible = [
        abs(float(row["diff"])) for row in rows
        if row["moment"] not in WEALTH_ROWS
        and "wealth" not in str(row["moment"])
        and "estate" not in str(row["moment"])
    ]
    return max(eligible) if eligible else math.nan


def run_case(case: str) -> dict[str, Any]:
    metadata = CASES[case]
    record_path = Path(metadata["record"])
    winner = load_winner(record_path)
    theta = {name: float(value) for name, value in winner["theta"].items()}
    moments, solution, repeats_bitwise = solve_twice(case, metadata, theta)

    # Re-importing the small helper preserves the chain's table calculation.
    from intergen_eqscale_seq_optimized import run_e1_chain as chain
    from intergen_eqscale_seq_optimized.e1_profile import E1_TARGET_SET

    targets, weights = chain.get_target_set(E1_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = 5.779970481941968
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    fit = chain.target_fit(moments, targets, weights)
    repaired_loss = float(sum(float(row["loss_contribution"]) for row in fit))
    stored_moments = winner["moments"]
    missing = [name for name in targets if name not in stored_moments]
    if missing:
        raise RuntimeError(f"{case}: stored winner lacks target moments: {', '.join(missing)}")
    comparison = [
        {
            "moment": name,
            "stored": float(stored_moments[name]),
            "repaired": float(moments[name]),
            "diff": float(moments[name]) - float(stored_moments[name]),
        }
        for name in targets
    ]

    outdir = OUTPUT_ROOT / case
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(
        outdir / "target_fit_full.csv", fit,
        ["moment", "target", "model", "gap", "weight", "loss_contribution"],
    )
    write_csv(outdir / "stored_vs_repaired.csv", comparison, ["moment", "stored", "repaired", "diff"])
    summary = {
        "case": case,
        "winner_record": str(record_path),
        "output_directory": str(outdir),
        "l4_enabled": bool(metadata["l4"]),
        "target_set": E1_TARGET_SET,
        "repaired_loss": repaired_loss,
        "stored_rank_loss": float(winner.get("rank_loss", math.nan)),
        "repeats_bitwise": repeats_bitwise,
        "market_residual": float(solution.best_max_abs_rel_excess),
        "wealth_moment_timing": str(getattr(solution, "wealth_moment_timing", "")),
        "current_distribution_timing": str(getattr(solution, "current_distribution_timing", "")),
        "p90_p50": float(getattr(solution, P90_P50, math.nan)),
        "p90_p50_legacy_alias": float(getattr(solution, P90_P50_LEGACY, math.nan)),
        "theta": theta,
        "nonwealth_max_abs_diff": nonwealth_difference(comparison),
    }
    write_json(outdir / "summary.json", summary)
    return {**summary, "comparison": comparison}


def write_readme(results: list[dict[str, Any]]) -> None:
    lines = [
        "# E-series timing-repair readout",
        "",
        "Each row re-solves the collected E1 winner twice with the strict evaluator "
        "after the wealth-timing repair. Stored values are cluster records and are "
        "reported for comparison only.",
        "",
        "| Case | Stored loss | Repaired loss | Young renter wealth (stored → repaired) | "
        "Old estate median (stored → repaired) | Old nonhousing share (stored → repaired) | "
        "Old p90/p50 (stored → repaired) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        values = {row["moment"]: row for row in result["comparison"]}
        def shift(name: str) -> str:
            row = values.get(name)
            return "n/a" if row is None else f"{row['stored']:.6g} → {row['repaired']:.6g}"
        lines.append(
            f"| {result['case']} | {result['stored_rank_loss']:.9g} | {result['repaired_loss']:.9g} | "
            f"{shift('young_childless_renter_liquid_wealth_to_annual_gross_income_2535')} | "
            f"{shift('old_total_estate_wealth_to_annual_income_median_7684')} | "
            f"{shift('old_nonhousing_ge_1x_income_share_6575')} | {shift(P90_P50_LEGACY)} |"
        )
    (OUTPUT_ROOT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    cases = selected_cases(args.cases)
    results = []
    for case in cases:
        print(f"running {case}", flush=True)
        results.append(run_case(case))
    write_readme(results)


if __name__ == "__main__":
    main()
