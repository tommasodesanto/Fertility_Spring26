#!/usr/bin/env python3
"""Synthetic packet tests for the policy mechanism collector."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import collect_e5f_post2023_policy_mechanisms as collector
import run_e5f_post2023_policy_mechanisms as policy


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def make_case(root: Path, case: str) -> None:
    case_dir = root / case
    case_dir.mkdir(parents=True)
    rows = []
    for year in (2023, 2027):
        shift = 0.0 if case == "baseline" else 0.01
        rows.append(
            {
                "policy_case": case,
                "policy_label": policy.POLICIES[case].label,
                "calendar_year": year,
                "asset_price": 1.0 - shift,
                "housing_user_cost": 0.08 + shift,
                "owner_rate": 0.6 + shift,
                "dependent_child_owner_rate": 0.55 + shift,
                "topcode_adjusted_births_per_adult": 0.08 + shift,
                "population_index_2023": 1.0,
                "housing_demand_per_adult": 5.0 + shift,
                "relative_market_residual": 1e-8,
                "mass_accounting_residual": 1e-12,
                "feasibility_frontier_projection_mass": 0.0,
                "annual_property_tax_rate": policy.POLICIES[
                    case
                ].annual_property_tax_rate,
                "property_tax_revenue_discarded": 0.1,
                "property_tax_lump_sum_transfer": 0.0,
                "purchase_grant": False,
            }
        )
    path = case_dir / "policy_path.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        case_dir / "summary.json",
        {
            "status": "complete_post2023_policy_mechanism_case",
            "policy": {
                "case": case,
                "supply_intercept_multiplier": policy.POLICIES[case].supply_multiplier,
                "dependent_child_ltv95": policy.POLICIES[case].dependent_child_ltv95,
                "annual_property_tax_rate": policy.POLICIES[
                    case
                ].annual_property_tax_rate,
                "fiscal_treatment": (
                    "property-tax revenue is discarded; no rebate and no purchase grant"
                ),
            },
            "post_2023_periods": 1,
            "input_hashes": {"contract": "same"},
            "reconstruction": {
                "matched_2023_distribution_sha256": "distribution",
                "matched_2023_queue_sha256": "queue",
            },
        },
    )
    write_json(
        case_dir / "manifest.json",
        {
            "schema": "e5f_post2023_policy_mechanism_manifest_v1",
            "artifacts": {"policy_path.csv": sha(path)},
        },
    )


def test_full_collection_and_corruption_rejection() -> None:
    with tempfile.TemporaryDirectory() as tmp_raw:
        tmp = Path(tmp_raw)
        root = tmp / "cases"
        for case in collector.CASE_ORDER:
            make_case(root, case)
        out = tmp / "report"
        old_argv = sys.argv
        try:
            sys.argv = [
                "collector",
                "--case-root",
                str(root),
                "--output-dir",
                str(out),
                "--expected-post-2023-periods",
                "1",
            ]
            collector.main()
        finally:
            sys.argv = old_argv
        assert (out / "summary.json").is_file()
        summary = json.loads((out / "summary.json").read_text())
        assert summary["status"] == "complete_post2023_policy_mechanism_comparison"
        assert summary["terminal_year"] == 2027

        bad = root / "combined" / "policy_path.csv"
        bad.write_text(bad.read_text() + "\n")
        try:
            collector.validate_case(root / "combined", "combined")
        except RuntimeError as error:
            assert "hash mismatch" in str(error)
        else:
            raise AssertionError("Collector accepted a corrupted policy path")


if __name__ == "__main__":
    test_full_collection_and_corruption_rejection()
    print("POST2023_POLICY_COLLECTOR_TESTS_PASS tests=1")
