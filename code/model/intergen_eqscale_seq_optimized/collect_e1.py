#!/usr/bin/env python3
"""Pick the strict, exactly repeated E1 winner and write a compact report."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULTS_ROOT = ROOT / "output/model/eqscale_seq_recalibration_20260718/production"
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_recalibration_20260718/report"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def eligible_tight(summary: dict[str, Any]) -> dict[str, Any] | None:
    tight = summary.get("best_tight")
    repeat = summary.get("tight_repeat_check") or {}
    if not (isinstance(tight, dict) and tight.get("strict_converged") and repeat.get("both_strict")
            and float(repeat.get("loss_abs_difference", math.inf)) == 0.0
            and float(repeat.get("max_abs_moment_difference", math.inf)) == 0.0):
        return None
    return tight


def parameter_rows(
    winner: dict[str, Any],
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    """Report every free parameter in the scale used by its search bound."""
    theta = winner.get("theta") or {}
    domain = metadata.get("active_domain") or []
    expected = int(metadata.get("free_parameter_count", -1))
    if len(domain) != expected:
        raise RuntimeError(
            f"active domain has {len(domain)} rows but metadata declares {expected} free parameters"
        )
    rows: list[dict[str, Any]] = []
    for spec in domain:
        name = str(spec["name"])
        theta_name = "beta" if name == "beta_annual" else name
        if theta_name not in theta:
            raise RuntimeError(f"winner theta omits free parameter {theta_name}")
        estimate = float(theta[theta_name])
        if name == "beta_annual":
            estimate = estimate ** 0.25
        lower, upper = float(spec["lower"]), float(spec["upper"])
        span = upper - lower
        lower_distance, upper_distance = estimate - lower, upper - estimate
        distance = min(lower_distance, upper_distance)
        near_bound = bool(distance <= 0.02 * span)
        external_restriction = ">= 0.94" if name == "beta_annual" else ""
        restriction_satisfied = estimate >= 0.94 if name == "beta_annual" else ""
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower_bound": lower,
                "upper_bound": upper,
                "transform": spec.get("transform"),
                "distance_to_nearest_bound": distance,
                "distance_as_share_of_range": distance / span,
                "near_bound_2pct": near_bound,
                "near_bound_side": (
                    "lower" if near_bound and lower_distance <= upper_distance
                    else "upper" if near_bound
                    else ""
                ),
                "external_restriction": external_restriction,
                "external_restriction_satisfied": restriction_satisfied,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    paths = sorted(args.results_root.glob("chain_*/summary.json"))
    if not paths:
        raise RuntimeError(f"no E1 chain summaries under {args.results_root}")
    eligible: list[
        tuple[dict[str, Any], dict[str, Any], dict[str, Any], str]
    ] = []
    chain_rows: list[dict[str, Any]] = []
    for path in paths:
        summary = json.loads(path.read_text())
        meta = summary.get("metadata") or {}
        is_legacy = meta.get("arm") in {"E1", "E3_L4", "E4_SPLIT"} and int(meta.get("free_parameter_count", -1)) == 12 and int(meta.get("target_count", -1)) == 15
        is_e5 = meta.get("arm") == "E5" and int(meta.get("free_parameter_count", -1)) == 10 and int(meta.get("target_count", -1)) == 12
        is_e6 = meta.get("arm") in {"E6A", "E6B", "E6AB"} and int(meta.get("free_parameter_count", -1)) == 10 and int(meta.get("target_count", -1)) == 12
        is_e5_probe = meta.get("arm") == "E5_PROBE_KE" and int(meta.get("free_parameter_count", -1)) == 9 and int(meta.get("target_count", -1)) == 12
        if not (is_legacy or is_e5 or is_e6 or is_e5_probe):
            raise RuntimeError(f"{path} does not satisfy an accepted E1/E3/E4/E5/E6 calibration contract")
        tight = eligible_tight(summary)
        chain_rows.append({"summary_path": str(path), "arm": meta.get("arm"), "seed": meta.get("seed"), "eligible": tight is not None,
                           "search_cases": summary.get("n_cases_completed"), "search_strict": summary.get("n_strict"),
                           "tight_loss": tight.get("rank_loss") if tight else None,
                           "tight_residual": tight.get("market_residual") if tight else None})
        if tight is not None:
            eligible.append(
                (tight, meta, summary.get("tight_repeat_check") or {}, str(path))
            )
    if not eligible:
        raise RuntimeError("no strict, exactly repeated E1 winner")
    winner, winner_metadata, winner_repeat_check, winner_summary_path = min(
        eligible, key=lambda item: float(item[0]["rank_loss"])
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)
    write_csv(args.outdir / "target_fit_full.csv", list(winner["target_fit"]))
    write_csv(args.outdir / "parameter_table_full.csv", parameter_rows(winner, winner_metadata))
    results = {
        "winners": {"E1": winner},
        "winner_metadata": winner_metadata,
        "winner_repeat_check": winner_repeat_check,
        "winner_summary_path": winner_summary_path,
        "eligible_E1_chains": len(eligible),
        "production_chain_count": len(paths),
    }
    (args.outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    (args.outdir / "README.md").write_text("\n".join([
        "# E1 collector",
        "",
        "This collector scans `production/chain_*/summary.json`.",
        "It retains only strict, exactly repeated tight winners.",
        "Exact equality is required for loss and target moments.",
        "The lowest eligible tight-repeat loss is the E1 winner.",
        "`results.json` exposes `winners.E1` in the M5 winner shape.",
        "`results.json` also preserves the selected chain's complete calibration metadata.",
        "`results.json` preserves the exact-repeat check and source summary path.",
        "`target_fit_full.csv` contains every row in the selected target system.",
        "`parameter_table_full.csv` contains every estimated parameter, its search bounds, and a two-percent near-bound flag.",
        "`chain_summary.csv` records eligibility for every discovered chain.",
        "No model solve is performed by this collector.",
    ]) + "\n")


if __name__ == "__main__":
    main()
