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


def main() -> None:
    args = parse_args()
    paths = sorted(args.results_root.glob("chain_*/summary.json"))
    if not paths:
        raise RuntimeError(f"no E1 chain summaries under {args.results_root}")
    eligible: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    for path in paths:
        summary = json.loads(path.read_text())
        meta = summary.get("metadata") or {}
        if meta.get("arm") != "E1" or int(meta.get("free_parameter_count", -1)) != 12 or int(meta.get("target_count", -1)) != 15:
            raise RuntimeError(f"{path} does not satisfy the E1 12-parameter/15-target contract")
        tight = eligible_tight(summary)
        chain_rows.append({"summary_path": str(path), "seed": meta.get("seed"), "eligible": tight is not None,
                           "search_cases": summary.get("n_cases_completed"), "search_strict": summary.get("n_strict"),
                           "tight_loss": tight.get("rank_loss") if tight else None,
                           "tight_residual": tight.get("market_residual") if tight else None})
        if tight is not None:
            eligible.append(tight)
    if not eligible:
        raise RuntimeError("no strict, exactly repeated E1 winner")
    winner = min(eligible, key=lambda row: float(row["rank_loss"]))
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)
    write_csv(args.outdir / "target_fit_full.csv", list(winner["target_fit"]))
    results = {"winners": {"E1": winner}, "eligible_E1_chains": len(eligible), "production_chain_count": len(paths)}
    (args.outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    (args.outdir / "README.md").write_text("\n".join([
        "# E1 collector",
        "",
        "This collector scans `production/chain_*/summary.json`.",
        "It retains only strict, exactly repeated tight winners.",
        "Exact equality is required for loss and target moments.",
        "The lowest eligible tight-repeat loss is the E1 winner.",
        "`results.json` exposes `winners.E1` in the M5 winner shape.",
        "`target_fit_full.csv` contains all 15 unchanged target rows.",
        "`chain_summary.csv` records eligibility for every discovered chain.",
        "No model solve is performed by this collector.",
    ]) + "\n")


if __name__ == "__main__":
    main()
