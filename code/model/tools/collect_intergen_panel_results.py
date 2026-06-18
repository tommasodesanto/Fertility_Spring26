#!/usr/bin/env python3
"""Collect intergen_housing_fertility local-panel task outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


MOMENT_FIELDS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_family_gap",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "young_liquid_wealth_to_income",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
    "old_nonhousing_wealth_to_income_6575",
    "old_nonhousing_wealth_to_income_median_6575",
    "old_total_wealth_to_income_6575",
    "old_total_wealth_to_income_median_6575",
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "old_parent_childless_nonhousing_wealth_to_income_median_gap_6575",
    "old_parent_childless_total_wealth_to_income_gap_6575",
    "old_parent_childless_total_wealth_to_income_median_gap_6575",
    "liquid_wealth_to_income",
    "housing_user_cost_share",
    "prime_childless_renter_median_rooms",
    "prime_childless_owner_median_rooms",
    "prime30_55_childless_renter_mean_rooms",
    "prime30_55_childless_owner_mean_rooms",
    "prime30_55_childless_renter_share_rooms_ge6",
    "prime30_55_childless_owner_share_rooms_ge6",
    "mean_age_first_birth",
    "parity_share_0",
    "parity_share_1",
    "parity_share_2plus",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    results_dir = args.results_dir
    outdir = args.outdir or results_dir
    outdir.mkdir(parents=True, exist_ok=True)

    task_dirs = sorted(p for p in results_dir.glob("task_*") if p.is_dir())
    records: list[dict[str, Any]] = []
    task_summaries: list[dict[str, Any]] = []
    for task_dir in task_dirs:
        task_id = task_dir.name.removeprefix("task_")
        summary = read_json(task_dir / "summary.json")
        if summary:
            meta = dict(summary.get("metadata", {}))
            task_summaries.append(
                {
                    "task_id": task_id,
                    "n_cases_completed": summary.get("n_cases_completed"),
                    "n_cases_submitted": summary.get("n_cases_submitted"),
                    "elapsed_sec": summary.get("elapsed_sec"),
                    "stopped_by_time_budget": summary.get("stopped_by_time_budget"),
                    "target_set": meta.get("rank_target_set"),
                    "seed": meta.get("seed"),
                    "include_anchors": meta.get("include_anchors"),
                    "best_rank_loss": nested(summary, "best", "rank_loss"),
                    "best_case": nested(summary, "best", "case"),
                    "best_label": nested(summary, "best", "label"),
                }
            )
        for record in read_jsonl(task_dir / "cases.jsonl"):
            record = dict(record)
            record["task_id"] = task_id
            records.append(record)

    ok_records = [r for r in records if math.isfinite(as_float(r.get("rank_loss"), math.inf))]
    ranked = sorted(ok_records, key=lambda r: as_float(r.get("rank_loss"), math.inf))
    best = ranked[0] if ranked else None

    summary_out = {
        "results_dir": str(results_dir),
        "n_tasks_found": len(task_dirs),
        "n_task_summaries": len(task_summaries),
        "n_case_records": len(records),
        "n_ok_records": len(ok_records),
        "best": compact_record(best) if best else None,
        "case_elapsed_sec": summarize([as_float(r.get("elapsed_sec"), math.nan) for r in ok_records]),
        "markov_income_solve_time_sec": summarize(
            [as_float(dict(r.get("timings", {})).get("markov_income_solve_time"), math.nan) for r in ok_records]
        ),
        "iterations_completed": summarize(
            [as_float(dict(r.get("timings", {})).get("iterations_completed"), math.nan) for r in ok_records]
        ),
        "task_summaries": task_summaries,
    }
    (outdir / "collected_summary.json").write_text(json.dumps(summary_out, indent=2, sort_keys=True))
    write_ranked_cases(outdir / "collected_ranked_cases.csv", ranked[: max(1, int(args.top_n))])
    write_task_summary(outdir / "collected_task_summary.csv", task_summaries)

    print(f"tasks={len(task_dirs)} cases={len(records)} ok={len(ok_records)}")
    if best:
        print(
            "best "
            f"task={best.get('task_id')} case={best.get('case')} "
            f"label={best.get('label')} rank_loss={as_float(best.get('rank_loss'), math.nan):.6g}"
        )
    print(f"wrote {outdir / 'collected_summary.json'}")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def nested(obj: dict[str, Any], *keys: str) -> Any:
    cur: Any = obj
    for key in keys:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def as_float(value: Any, default: float) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out


def summarize(values: list[float]) -> dict[str, float | int | None]:
    vals = sorted(v for v in values if math.isfinite(v))
    if not vals:
        return {"n": 0, "min": None, "p25": None, "p50": None, "p75": None, "max": None, "mean": None}
    n = len(vals)
    return {
        "n": n,
        "min": vals[0],
        "p25": vals[n // 4],
        "p50": vals[n // 2],
        "p75": vals[(3 * n) // 4],
        "max": vals[-1],
        "mean": sum(vals) / n,
    }


def compact_record(record: dict[str, Any]) -> dict[str, Any]:
    moments = dict(record.get("moments", {}))
    timings = dict(record.get("timings", {}))
    return {
        "task_id": record.get("task_id"),
        "case": record.get("case"),
        "label": record.get("label"),
        "rank_loss": record.get("rank_loss"),
        "full_old_nonlocation_loss": record.get("full_old_nonlocation_loss"),
        "market_residual": record.get("market_residual"),
        "elapsed_sec": record.get("elapsed_sec"),
        "markov_income_solve_time": timings.get("markov_income_solve_time"),
        "iterations_completed": timings.get("iterations_completed"),
        "moments": {name: moments.get(name) for name in MOMENT_FIELDS if name in moments},
        "theta": record.get("theta"),
    }


def write_ranked_cases(path: Path, records: list[dict[str, Any]]) -> None:
    fieldnames = [
        "task_id",
        "case",
        "label",
        "rank_loss",
        "full_old_nonlocation_loss",
        "market_residual",
        "elapsed_sec",
        "markov_income_solve_time",
        "iterations_completed",
        *MOMENT_FIELDS,
        "theta_json",
    ]
    rows = []
    for record in records:
        moments = dict(record.get("moments", {}))
        timings = dict(record.get("timings", {}))
        row = {
            "task_id": record.get("task_id"),
            "case": record.get("case"),
            "label": record.get("label"),
            "rank_loss": record.get("rank_loss"),
            "full_old_nonlocation_loss": record.get("full_old_nonlocation_loss"),
            "market_residual": record.get("market_residual"),
            "elapsed_sec": record.get("elapsed_sec"),
            "markov_income_solve_time": timings.get("markov_income_solve_time"),
            "iterations_completed": timings.get("iterations_completed"),
            "theta_json": json.dumps(record.get("theta", {}), sort_keys=True),
        }
        for name in MOMENT_FIELDS:
            row[name] = moments.get(name)
        rows.append(row)
    write_csv(path, fieldnames, rows)


def write_task_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "task_id",
        "n_cases_completed",
        "n_cases_submitted",
        "elapsed_sec",
        "stopped_by_time_budget",
        "target_set",
        "seed",
        "include_anchors",
        "best_rank_loss",
        "best_case",
        "best_label",
    ]
    write_csv(path, fieldnames, rows)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name) for name in fieldnames})


if __name__ == "__main__":
    main()
