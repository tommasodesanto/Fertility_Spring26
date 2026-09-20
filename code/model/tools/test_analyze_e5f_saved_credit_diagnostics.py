#!/usr/bin/env python3
"""Meaningful pure checks for the saved-credit diagnostic."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import analyze_e5f_saved_credit_diagnostics as diag


class SavedCreditDiagnosticTests(unittest.TestCase):
    def test_weighted_mass_mean_and_endpoints(self) -> None:
        total, mean, endpoints = diag.weighted_mass_stats(np.array([0.25, 0.75]), np.array([-1.0, 3.0]))
        self.assertAlmostEqual(total, 1.0)
        self.assertAlmostEqual(mean, 2.0)
        self.assertEqual(endpoints, (-1.0, 3.0))

    def test_cohort_group_rows_uses_age_specific_mass_and_rejects_locations(self) -> None:
        g = np.zeros((1, 3, 1, 1, 1, 1, 1))
        g[0, :, 0, 0, 0, 0, 0] = [0.2, 0.3, 0.5]
        stacked = g[None, ...]
        rows = diag.cohort_group_rows(stacked, np.array([-1.0]), np.array([0.0, 0.0]),
                                      np.array([1.0, 1.0]), np.array([18.0]), 0.0)
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(sum(row["mass"] for row in rows), 1.0)
        self.assertEqual(rows[0]["grid_floor"], -1.0)
        bad = np.zeros((1, 3, 2, 1, 1, 1, 1))
        with self.assertRaises(ValueError):
            diag.group_rows(bad, np.array([-1.0]), np.array([0.0, 0.0]),
                            np.array([1.0, 1.0]), np.array([18.0]), 0.0)

    def test_identity_report_detects_plateau_and_state_difference(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            a = root / "lambda1.npz"
            b = root / "lambda5.npz"
            c = root / "lambda025.npz"
            values = {key: np.array([1.0, 2.0]) for key in diag.POLICY_KEYS}
            values["bp_pol"] = np.array([0.0, 1.0])
            values["tenure_choice"] = np.array([True, False])
            values["price"] = np.array([np.inf, np.inf])
            np.savez(a, **values)
            np.savez(b, **values)
            values["bp_pol"] = np.array([0.0, 1.5])
            np.savez(c, **values)
            equal = diag.identity_report(a, b)
            self.assertTrue(all(row.get("identical_atol", False) for row in equal["rows"]))
            unequal = diag.identity_report(a, c)
            bp = next(row for row in unequal["rows"] if row["array"] == "bp_pol")
            self.assertEqual(unequal["status"], "complete")
            self.assertFalse(unequal["identical_all_atol"])
            self.assertEqual(bp["different_entries"], 1)
            self.assertGreater(bp["max_abs_diff"], 0.0)
            missing = root / "missing.npz"
            np.savez(missing, g_pre=np.array([1.0]))
            self.assertEqual(diag.identity_report(a, missing)["status"], "failed")

    def test_flow_report_recomputes_first_and_continuation_sums(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "cohort_by_age.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["age_years", "explicit_birth_flow", "exact_first_birth_flow"])
                writer.writeheader()
                writer.writerow({"age_years": 18, "explicit_birth_flow": 0.6, "exact_first_birth_flow": 0.5})
                writer.writerow({"age_years": 22, "explicit_birth_flow": 0.7, "exact_first_birth_flow": 0.4})
            result = diag.flow_report(path, {"cohort": {"cumulative_explicit_births_per_initial_household": 1.3, "first_births_per_initial_household": 0.9}})
            self.assertEqual(result["status"], "complete")
            self.assertAlmostEqual(result["explicit_sum"], 1.3)
            self.assertAlmostEqual(result["first_sum"], 0.9)
            self.assertAlmostEqual(result["continuation_sum"], 0.4)
            self.assertAlmostEqual(result["explicit_gap"], 0.0)
            bad = diag.flow_report(path, {"cohort": {"cumulative_explicit_births_per_initial_household": 1.2,
                                                       "first_births_per_initial_household": 0.9}})
            self.assertEqual(bad["status"], "failed")

    def test_select_rows_keeps_only_phi08_cap6_doses(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "summary.json"
            rows = []
            for lam in (0.0, 0.25, 1.0, 5.0):
                rows.append({"label": f"arm_{lam}", "phi": 0.8, "lambda": lam, "rental_cap": 6.0,
                             "status": "completed", "contract": {"checkpoint": "checkpoint", "source_root": "frozen", "source_manifest": "manifest", "checkpoint_sha256": "hash"}})
            rows.append({"label": "wrong", "phi": 1.0, "lambda": 5.0, "rental_cap": 6.0, "status": "completed", "contract": {}})
            path.write_text(json.dumps({"status": "complete", "cases": rows}))
            selected, info = diag.select_rows("refit_new_income", path, False)
            self.assertEqual([row["lambda"] for row in selected], [0.0, 0.25, 1.0, 5.0])
            self.assertEqual(info["errors"], [])


if __name__ == "__main__":
    unittest.main()
