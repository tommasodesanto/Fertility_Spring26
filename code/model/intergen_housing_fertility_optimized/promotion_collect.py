#!/usr/bin/env python3
"""Collect the isolated promotion battery and adjudicate every declared gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .m5_profile import m5_target_system
from .promotion_contract import CALIBRATION_CASES, DEMAND_PRICES, ROOT_CASES, bound_cases


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def load_records(folder: Path, expected: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index in range(1, expected + 1):
        path = folder / f"task_{index}.json"
        if path.exists():
            records.append(json.loads(path.read_text()))
    return records


def read_json(path: Path) -> dict[str, Any] | None:
    return json.loads(path.read_text()) if path.exists() else None


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bound_gate(records: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for record in records:
        rows.append(
            {
                "index": record["index"],
                "case": record["case"]["name"],
                "production_status": record["production"]["status"],
                "optimized_status": record["optimized"]["status"],
                "pass": bool(record["comparison"]["pass"]),
            }
        )
    return {
        "expected": len(bound_cases()),
        "completed": len(records),
        "passed": sum(row["pass"] for row in rows),
        "rows": rows,
        "pass": bool(len(records) == len(bound_cases()) and all(row["pass"] for row in rows)),
    }


def demand_gate(records: list[dict[str, Any]]) -> dict[str, Any]:
    ordered = sorted(records, key=lambda record: float(record["price"]))
    feasible = [record for record in ordered if record["production"]["status"] == "ok"]
    sign_changes = []
    monotonicity_violations = []
    for left, right in zip(feasible, feasible[1:]):
        left_excess = float(left["production"]["excess_demand"])
        right_excess = float(right["production"]["excess_demand"])
        if left_excess * right_excess <= 0.0:
            sign_changes.append([float(left["price"]), float(right["price"])])
        if right_excess > left_excess + 1e-10:
            monotonicity_violations.append(
                {
                    "left_price": float(left["price"]),
                    "right_price": float(right["price"]),
                    "left_excess": left_excess,
                    "right_excess": right_excess,
                }
            )
    parity_pass = all(bool(record["comparison"]["pass"]) for record in ordered)
    endpoints_present = bool(
        ordered
        and float(ordered[0]["price"]) == min(DEMAND_PRICES)
        and float(ordered[-1]["price"]) == max(DEMAND_PRICES)
    )
    return {
        "expected": len(DEMAND_PRICES),
        "completed": len(records),
        "feasible_points": len(feasible),
        "full_configured_endpoints_present": endpoints_present,
        "sign_change_brackets": sign_changes,
        "monotonicity_violations": monotonicity_violations,
        "parity_pass": parity_pass,
        "pass": bool(
            len(records) == len(DEMAND_PRICES)
            and endpoints_present
            and parity_pass
            and sign_changes
        ),
    }


def root_gate(records: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for record in records:
        production = record["production"]
        optimized = record["optimized"]
        refine = optimized.get("timings", {}).get("scalar_market_refine", {})
        rows.append(
            {
                "case": record["case"]["name"],
                "production_status": production["status"],
                "optimized_status": optimized["status"],
                "production_seconds": production.get("elapsed_seconds"),
                "optimized_seconds": optimized.get("elapsed_seconds"),
                "speedup": record["comparison"].get("optimized_speedup"),
                "optimized_residual": optimized.get("market_residual"),
                "bracket_found": refine.get("bracket_found"),
                "directional_fallback_used": refine.get("directional_fallback_used"),
                "pass": bool(record["comparison"]["pass"]),
            }
        )
    speedups = [float(row["speedup"]) for row in rows if row["speedup"] is not None]
    return {
        "expected": len(ROOT_CASES),
        "completed": len(records),
        "median_speedup": float(np.median(speedups)) if speedups else None,
        "minimum_speedup": min(speedups) if speedups else None,
        "actual_directional_fallback_count": sum(
            row["directional_fallback_used"] is True for row in rows
        ),
        "synthetic_wrong-direction_fallback_test": (
            "test_wrong_direction_uses_bilateral_fallback_and_finds_root"
        ),
        "rows": rows,
        "pass": bool(len(records) == len(ROOT_CASES) and all(row["pass"] for row in rows)),
    }


def repeat_gate(records: list[dict[str, Any]]) -> dict[str, Any]:
    if len(records) != 2:
        return {"expected": 2, "completed": len(records), "pass": False}
    left = records[0]["optimized"]
    right = records[1]["optimized"]
    moment_max = max(
        abs(float(left["moments"][name]) - float(right["moments"][name]))
        for name in left["moments"]
    )
    comparisons = {
        "price_abs_difference": abs(float(left["price"][0]) - float(right["price"][0])),
        "loss_abs_difference": abs(float(left["loss"]) - float(right["loss"])),
        "mass_abs_difference": abs(float(left["mass"]) - float(right["mass"])),
        "residual_abs_difference": abs(
            float(left["market_residual"]) - float(right["market_residual"])
        ),
        "moment_max_abs_difference": moment_max,
    }
    return {
        "expected": 2,
        "completed": 2,
        **comparisons,
        "pass": bool(left["status"] == right["status"] == "ok" and not any(comparisons.values())),
    }


def calibration_gate(
    production: dict[str, Any] | None,
    optimized: dict[str, Any] | None,
) -> dict[str, Any]:
    if production is None or optimized is None:
        return {"pass": False, "reason": "missing_package_result"}
    production_by_name = {
        record["case"]["name"]: record["result"] for record in production["records"]
    }
    optimized_by_name = {
        record["case"]["name"]: record["result"] for record in optimized["records"]
    }
    rows = []
    for case in CALIBRATION_CASES:
        name = case["name"]
        left = production_by_name.get(name)
        right = optimized_by_name.get(name)
        statuses_equal = bool(
            left is not None
            and right is not None
            and left["status"] == right["status"]
            and left["status"] in {"ok", "infeasible"}
        )
        rows.append(
            {
                "case": name,
                "production_status": left["status"] if left else "missing",
                "optimized_status": right["status"] if right else "missing",
                "production_loss": left.get("loss") if left else None,
                "optimized_loss": right.get("loss") if right else None,
                "loss_abs_difference": (
                    abs(float(left["loss"]) - float(right["loss"]))
                    if left and right and left["status"] == right["status"] == "ok"
                    else None
                ),
                "production_seconds": left.get("elapsed_seconds") if left else None,
                "optimized_seconds": right.get("elapsed_seconds") if right else None,
                "speedup": (
                    float(left["elapsed_seconds"]) / float(right["elapsed_seconds"])
                    if left and right and float(right.get("elapsed_seconds", 0.0)) > 0.0
                    else None
                ),
                "statuses_equal": statuses_equal,
            }
        )
    speedups = [float(row["speedup"]) for row in rows if row["speedup"] is not None]
    left_best = production.get("best")
    right_best = optimized.get("best")
    best_equal = bool(
        left_best
        and right_best
        and left_best["case"]["name"] == right_best["case"]["name"]
    )
    return {
        "status": "numerical_throughput_smoke_not_a_research_calibration",
        "expected_cases": len(CALIBRATION_CASES),
        "production_completed": len(production["records"]),
        "optimized_completed": len(optimized["records"]),
        "production_total_seconds": production.get("elapsed_seconds"),
        "optimized_total_seconds": optimized.get("elapsed_seconds"),
        "total_throughput_speedup": (
            float(production["elapsed_seconds"]) / float(optimized["elapsed_seconds"])
        ),
        "median_case_speedup": float(np.median(speedups)) if speedups else None,
        "same_selected_case": best_equal,
        "production_selected_case": left_best["case"]["name"] if left_best else None,
        "optimized_selected_case": right_best["case"]["name"] if right_best else None,
        "rows": rows,
        "pass": bool(
            len(production["records"]) == len(CALIBRATION_CASES)
            and len(optimized["records"]) == len(CALIBRATION_CASES)
            and all(row["statuses_equal"] for row in rows)
            and best_equal
        ),
    }


def target_fit_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    system = m5_target_system()
    rows = []
    for name, target, weight in zip(
        system.moment_names, system.target_values, system.weights
    ):
        model = float(result["moments"][name])
        gap = model - target
        rows.append(
            {
                "moment": name,
                "target": target,
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap**2,
            }
        )
    return rows


def write_report(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        "# Optimized package promotion battery",
        "",
        "This is a numerical promotion test, not a new model specification or paper calibration.",
        "",
        "| Gate | Result |",
        "|---|---:|",
    ]
    for name in ("bound_parity", "demand_bracket", "root_cases", "strict_repeats", "diagnostics", "calibration_throughput"):
        lines.append(f"| {name} | {'PASS' if summary[name]['pass'] else 'FAIL'} |")
    lines.extend(
        [
            "",
            f"Overall: **{'PASS' if summary['overall_pass'] else 'FAIL'}**.",
            "",
            "The calibration-throughput arm fixes the M5 target system and compares a predeclared",
            "ten-case objective surface. It is not an identified SMM estimate and must not replace M5.",
            "Complete machine-readable rows and selected-case target-fit tables are adjacent to this report.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    output = args.output or args.run_root / "collected"
    output.mkdir(parents=True, exist_ok=True)
    bounds = load_records(args.run_root / "bounds", len(bound_cases()))
    demand = load_records(args.run_root / "demand", len(DEMAND_PRICES))
    roots = load_records(args.run_root / "roots", len(ROOT_CASES))
    repeats = load_records(args.run_root / "repeats", 2)
    production_calibration = read_json(args.run_root / "calibration" / "production.json")
    optimized_calibration = read_json(args.run_root / "calibration" / "optimized.json")
    diagnostics_payload = read_json(args.run_root / "diagnostics" / "comparison.json")
    diagnostics_gate = (
        diagnostics_payload["comparison"]
        if diagnostics_payload is not None
        else {"pass": False, "reason": "missing"}
    )
    summary = {
        "bound_parity": bound_gate(bounds),
        "demand_bracket": demand_gate(demand),
        "root_cases": root_gate(roots),
        "strict_repeats": repeat_gate(repeats),
        "diagnostics": diagnostics_gate,
        "calibration_throughput": calibration_gate(
            production_calibration, optimized_calibration
        ),
    }
    summary["overall_pass"] = all(gate["pass"] for gate in summary.values())
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_csv(output / "bound_parity.csv", summary["bound_parity"].pop("rows"))
    write_csv(output / "root_cases.csv", summary["root_cases"].pop("rows"))
    write_csv(
        output / "calibration_throughput.csv",
        summary["calibration_throughput"].pop("rows", []),
    )
    if production_calibration and production_calibration.get("best"):
        write_csv(
            output / "production_selected_target_fit.csv",
            target_fit_rows(production_calibration["best"]["result"]),
        )
    if optimized_calibration and optimized_calibration.get("best"):
        write_csv(
            output / "optimized_selected_target_fit.csv",
            target_fit_rows(optimized_calibration["best"]["result"]),
        )
    # Rewrite the row-free compact summary after emitting the sidecar tables.
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_report(output / "REPORT.md", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    raise SystemExit(0 if summary["overall_pass"] else 1)


if __name__ == "__main__":
    main()
