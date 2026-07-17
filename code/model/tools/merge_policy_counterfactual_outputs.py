#!/usr/bin/env python3
"""Merge per-case policy counterfactual outputs into one presentation packet."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from run_policy_counterfactuals_from_record import write_csv, write_latex_table, write_summary_plot


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    baseline_written = False
    for case_dir in args.case_dir:
        summary = case_dir / "summary.csv"
        if not summary.exists():
            raise FileNotFoundError(summary)
        with summary.open() as handle:
            for row in csv.DictReader(handle):
                if row.get("case") == "baseline_fixed_outside":
                    if baseline_written:
                        continue
                    baseline_written = True
                rows.append(row)
    write_csv(args.outdir / "summary.csv", rows)
    write_latex_table(args.outdir / "summary_table.tex", rows)
    write_summary_plot(args.outdir / "summary_effects.png", rows)
    print(f"Wrote merged policy packet to {args.outdir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, action="append", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    main()
