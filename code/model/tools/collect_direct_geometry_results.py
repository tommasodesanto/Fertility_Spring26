#!/usr/bin/env python3
"""Collect direct no-inversion calibration JSON outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


SUMMARY_FIELDS = [
    "job_id",
    "eval_id",
    "loss",
    "solve_elapsed_sec",
    "tfr",
    "childless_rate",
    "own_rate",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "implied_total_population",
    "scale_factor",
    "inv_pop_share_C",
    "inv_rent_ratio_C_over_P",
    "beta",
    "psi_child",
    "h_bar_jump",
    "h_bar_n",
    "c_bar_n",
    "kappa_fert",
    "chi",
    "kappa_loc",
    "mu_move",
    "theta0",
    "theta_n",
    "h_bar_0",
    "E_C",
    "r_bar_C",
    "outside_value",
    "outside_entry_flow",
    "renewal_retention",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect direct calibration best records.")
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--csv", type=Path, default=None)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    records = []
    for best_path in sorted(args.results_dir.glob("job_*/best.json")):
        try:
            rec = json.loads(best_path.read_text())
        except json.JSONDecodeError:
            continue
        rec["_path"] = str(best_path)
        records.append(rec)

    records.sort(key=lambda rec: float(rec.get("loss", math.inf)))
    csv_path = args.csv or args.results_dir / "direct_geometry_summary.csv"
    json_path = args.json or args.results_dir / "direct_geometry_best.json"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(csv_path, records)
    json_path.write_text(json.dumps(records, indent=2, sort_keys=True))

    print(f"collected {len(records)} jobs")
    if records:
        best = records[0]
        moments = best.get("moments", {})
        print(
            "best: "
            f"job={best.get('job_id')} eval={best.get('eval_id')} loss={float(best.get('loss')):.6g} "
            f"TFR={float(moments.get('tfr', math.nan)):.3f} "
            f"own={float(moments.get('own_rate', math.nan)):.3f} "
            f"N={float(moments.get('implied_total_population', math.nan)):.3f} "
            f"popC={float(moments.get('inv_pop_share_C', math.nan)):.3f} "
            f"rr={float(moments.get('inv_rent_ratio_C_over_P', math.nan)):.3f}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for rec in records:
            row = flatten_record(rec)
            writer.writerow({key: row.get(key, "") for key in SUMMARY_FIELDS})


def flatten_record(rec: dict[str, Any]) -> dict[str, Any]:
    out = {
        "job_id": rec.get("job_id"),
        "eval_id": rec.get("eval_id"),
        "loss": rec.get("loss"),
        "solve_elapsed_sec": rec.get("solve_elapsed_sec"),
    }
    out.update(rec.get("moments", {}))
    out.update(rec.get("parameters", {}))
    return out


if __name__ == "__main__":
    main()
