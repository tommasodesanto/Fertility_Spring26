#!/usr/bin/env python3
"""Focused merge and state-comparison tests for the tenure-kink collector."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
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
