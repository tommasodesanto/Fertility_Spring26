"""Zero-lifecycle orchestration tests for the closed-price search."""
from __future__ import annotations

import tempfile
import time
import unittest
import csv
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import phase_b_ge as ge


class Budget:
    def __init__(self, count=9):
        self.remaining_lifecycle = count
        self.deadline_epoch = time.time() + 2400
        self.stage_deadline_seconds = 300
        self.claims = []

    def claim_lifecycle(self, label):
        self.remaining_lifecycle -= 1
        self.claims.append(label)


class PhaseBTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.writes = {}
        self.context = dict(out=self.out, q_ref=1., P=SimpleNamespace(psi_child=.135),
                            b_grid=np.array([0.]), fp=SimpleNamespace(write=self.write),
                            prepared=object(), manifest={}, objective={}, runtime=object(),
                            reference={})

    def write(self, path, value):
        self.writes[str(path)] = value

    def _solve(self, context, d_bar, q, budget, label, stage_dir):
        budget.claim_lifecycle(label)
        return dict(price=np.asarray([q]), sol=SimpleNamespace(x=np.asarray([d_bar, q])),
                    sd=SimpleNamespace(y=np.asarray([1.])), case_deadline_epoch=time.time() + 300)

    def _observe(self, context, live, label, *, final=False):
        q = float(live["price"][0])
        if final:
            directory = self.out / "phase_b_ge" / label
            directory.mkdir(parents=True, exist_ok=True)
            for filename, count, field in (("target_fit.csv", 14, "moment"),
                                           ("parameters.csv", 31, "parameter")):
                with (directory / filename).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=[field, "value"])
                    writer.writeheader()
                    writer.writerows({field: str(i), "value": "1"} for i in range(count))
        return dict(price=q, renewal_residual=(q - 1.02) * .01,
                    population_scale=1.03, final=final)

    def test_root_repeat_reuses_phase_a_stage_and_respects_cap(self):
        budget = Budget(9)
        selected = dict(selected_d_bar=.14, selected_live=self._solve(
            self.context, .14, 1., budget, "phase_a_selected_qref", self.out))
        with patch.object(ge, "solve_fixed_price", self._solve), patch.object(ge, "observe_price", self._observe):
            result = ge.run_phase_b(self.context, selected, budget)
        self.assertEqual(result["status"], "passed")
        self.assertAlmostEqual(result["selected_price"], 1.02)
        self.assertEqual(budget.claims.count("phase_a_selected_qref"), 1)
        self.assertEqual(budget.claims[-1], "selected_repeat")
        self.assertLessEqual(len(budget.claims), 10)
        self.assertEqual(result["selected_d_bar"], .14)

    def test_unbracketed_stops_with_no_repeat(self):
        budget = Budget(9)
        selected = dict(selected_d_bar=.14, selected_live=self._solve(
            self.context, .14, 1., budget, "phase_a_selected_qref", self.out))
        def same_sign(context, live, label, *, final=False):
            return dict(price=float(live["price"][0]), renewal_residual=.01,
                        population_scale=1.)
        with patch.object(ge, "solve_fixed_price", self._solve), patch.object(ge, "observe_price", same_sign):
            with self.assertRaisesRegex(RuntimeError, "unbracketed"):
                ge.run_phase_b(self.context, selected, budget)
        self.assertEqual(budget.claims, ["phase_a_selected_qref", "lower_085", "upper_115"])

    def test_no_budget_fails_before_new_solve(self):
        budget = Budget(1)
        with self.assertRaisesRegex(RuntimeError, "No GE solve"):
            ge.run_phase_b(self.context, dict(selected_d_bar=.14), budget)
        self.assertEqual(budget.claims, [])

    def test_zero_lifecycle_smoke_exact_loop(self):
        budget = Budget(10)
        result = ge.smoke_phase_b(self.context, budget)
        self.assertEqual(result["lifecycle_solves"], 0)
        self.assertEqual(result["status"], "passed_mock_zero_lifecycle")
        self.assertEqual(budget.claims, [])


if __name__ == "__main__":
    unittest.main()
