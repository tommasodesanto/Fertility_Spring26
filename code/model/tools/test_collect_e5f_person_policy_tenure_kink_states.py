#!/usr/bin/env python3
"""Focused merge and state-comparison tests for the tenure-kink collector."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from argparse import ArgumentTypeError
from pathlib import Path

import collect_e5f_person_policy_tenure_kink_states as collector


def row(year: int, wealth: int, choice: int, mass: float) -> dict[str, str]:
    payload = {name: "0" for name in collector.STATE_KEY_FIELDS}
    payload.update(
        calendar_year=str(year),
        wealth_index=str(wealth),
        stored_tenure=str(choice),
        decision_mass=str(mass),
        decision_mass_share=str(mass),
    )
    return payload


def write(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class MergeTests(unittest.TestCase):
    def test_expected_years_parser_accepts_new_kink_window(self) -> None:
        self.assertEqual(
            collector.parse_expected_years("2355,2359,2363"),
            (2355, 2359, 2363),
        )

    def test_expected_years_parser_rejects_duplicates_and_reordering(self) -> None:
        for value in ("2355,2355", "2359,2355", ""):
            with self.subTest(value=value), self.assertRaises(ArgumentTypeError):
                collector.parse_expected_years(value)

    def test_audit_packet_honors_custom_expected_years(self) -> None:
        expected_years = (2355, 2359, 2363)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (root / "transition_path.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "period",
                        "calendar_year",
                        "feasibility_frontier_projection_mass",
                    ),
                )
                writer.writeheader()
                for period, year in enumerate(expected_years):
                    writer.writerow(
                        {
                            "period": period,
                            "calendar_year": year,
                            "feasibility_frontier_projection_mass": 0.0,
                        }
                    )
            summary = {
                "status": "complete_unpromoted_person_demography_policy_path",
                "history_reproduction_status": "passed",
                "accounting_gates": {"synthetic_gate": True},
            }
            collector.write_json(root / "summary.json", summary)
            validations = []
            for year in expected_years:
                table = root / f"tenure_value_states_{year}.csv.gz"
                table.write_bytes(f"synthetic state table {year}\n".encode())
                validations.append(
                    {
                        "calendar_year": year,
                        "choice_mismatch_count": 0,
                        "maximum_compiled_value_error": 0.0,
                        "maximum_one_market_location_probability_gap": 0.0,
                        "state_csv": table.name,
                        "state_csv_sha256": collector.sha256(table),
                    }
                )
            collector.write_json(
                root / "tenure_diagnostic_summary.json",
                {
                    "status": "complete_unpromoted_tenure_value_diagnostic",
                    "initial_path_sha256": "synthetic-seed",
                    "production_summary_sha256": collector.sha256(root / "summary.json"),
                    "diagnostic_years": list(expected_years),
                    "validations": validations,
                },
            )
            _, _, tables, path_rows = collector.audit_packet(
                root,
                expected_initial_path_sha256="synthetic-seed",
                expected_years=expected_years,
            )
            self.assertEqual(tuple(sorted(tables)), expected_years)
            self.assertEqual(tuple(sorted(path_rows)), expected_years)
            with self.assertRaisesRegex(RuntimeError, "Diagnostic years differ"):
                collector.audit_packet(
                    root,
                    expected_initial_path_sha256="synthetic-seed",
                )

    def test_union_merge_preserves_missing_states(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pre = root / "pre.csv.gz"
            post = root / "post.csv.gz"
            write(pre, [row(2351, 0, 0, 0.2), row(2351, 2, 1, 0.3)])
            write(post, [row(2351, 1, 0, 0.4), row(2351, 2, 0, 0.3)])
            merged = list(collector.merge_state_rows(pre, post))
        self.assertEqual(len(merged), 3)
        self.assertIsNone(merged[0][2])
        self.assertIsNone(merged[1][1])
        self.assertEqual(merged[2][1]["stored_tenure"], "1")
        self.assertEqual(merged[2][2]["stored_tenure"], "0")

    def test_state_key_includes_every_declared_dimension(self) -> None:
        base = row(2351, 3, 0, 0.1)
        key = collector.state_key(base)
        self.assertEqual(len(key), len(collector.STATE_KEY_FIELDS))
        self.assertEqual(key[0], 2351)
        self.assertEqual(key[1], 3)

    def test_feasibility_gate_is_computed_from_saved_path_rows(self) -> None:
        clean = {
            2023: {"feasibility_frontier_projection_mass": "0"},
            2027: {"feasibility_frontier_projection_mass": "0.0"},
        }
        projected = {
            **clean,
            2031: {"feasibility_frontier_projection_mass": "2e-9"},
        }
        self.assertEqual(collector.maximum_path_feasibility_projection(clean), 0.0)
        self.assertEqual(
            collector.maximum_path_feasibility_projection(projected), 2.0e-9
        )


if __name__ == "__main__":
    unittest.main()
