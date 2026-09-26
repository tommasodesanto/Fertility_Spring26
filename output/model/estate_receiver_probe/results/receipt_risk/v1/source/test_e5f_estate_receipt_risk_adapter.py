"""Focused, unrun tests for the experimental estate receipt-risk adapter.

Run on Torch only:
``python -m unittest code/model/tools/test_e5f_estate_receipt_risk_adapter.py``.
These are adapter tests, not a native-solver or GE certification.
"""

import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from e5f_estate_receipt_risk_adapter import (
    _backward,
    _forward,
    configure,
    normalized_ledger,
    pooled_profile,
)


MODEL_AGES = list(range(18, 83, 4))


def _published_rows():
    """Complete 25--80 by three-group fixture with hand-checkable moments."""
    rows = []
    for age in range(25, 81):
        for group, probability, amount in (
            ("bottom50", 0.10, float(age)),
            ("middle40", 0.20, float(2 * age)),
            ("top10", 0.30, float(3 * age)),
        ):
            rows.append({"age": age, "income_group": group,
                         "probability_4y": probability,
                         "mean_amount_4y": probability * amount,
                         "conditional_amount_4y": amount})
    return rows


def _write_profile(directory, rows):
    path = Path(directory) / "published.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("age", "income_group", "probability_4y",
                                                    "mean_amount_4y", "conditional_amount_4y"))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _parameters():
    return SimpleNamespace(J=len(MODEL_AGES), age_start=18.0, da=4.0,
                           estate_probe_case="net_valuation", estate_probe_transfer=0.0,
                           income_transfer=7.25, earnings_level=3.5, preference_marker="fixed")


class EstateReceiptRiskAdapterTests(unittest.TestCase):
    def test_pooled_profile_uses_exact_groups_and_explicit_tail_zeroes(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = pooled_profile(_write_profile(directory, _published_rows()), MODEL_AGES)
        self.assertEqual([row["age"] for row in profile], [float(age) for age in MODEL_AGES])
        self.assertEqual([row["published_age_supported"] for row in profile],
                         [False, False] + [True] * 14 + [False])
        row = next(row for row in profile if row["age"] == 26.0)
        expected_p = .5 * .10 + .4 * .20 + .1 * .30
        expected_mu = .5 * (.10 * 26.) + .4 * (.20 * 52.) + .1 * (.30 * 78.)
        self.assertEqual(row["probability"], expected_p)
        self.assertEqual(row["relative_mean"], expected_mu)
        self.assertEqual(row["relative_positive_amount"], expected_mu / expected_p)
        self.assertEqual(row["probability"] * row["relative_positive_amount"], row["relative_mean"])
        for age in (18., 22., 82.):
            tail = next(row for row in profile if row["age"] == age)
            self.assertEqual((tail["probability"], tail["relative_mean"],
                              tail["relative_positive_amount"]), (0., 0., 0.))

    def test_pooled_profile_rejects_duplicate_incomplete_and_invalid_rows(self):
        rows = _published_rows()
        invalid_sets = [
            rows + [dict(rows[0])],
            rows[:-1],
            [dict(row, probability_4y=1.1) if row is rows[0] else row for row in rows],
            [dict(row, probability_4y=float("nan")) if row is rows[0] else row for row in rows],
            [dict(row, mean_amount_4y=float("nan")) if row is rows[0] else row for row in rows],
            [dict(row, conditional_amount_4y=float("nan")) if row is rows[0] else row for row in rows],
            [dict(row, mean_amount_4y=999.) if row is rows[0] else row for row in rows],
            [dict(row, age=25.5) if row is rows[0] else row for row in rows],
        ]
        with tempfile.TemporaryDirectory() as directory:
            for bad_rows in invalid_sets:
                with self.assertRaises(ValueError):
                    pooled_profile(_write_profile(directory, bad_rows), MODEL_AGES)

    def test_configure_preserves_means_base_inputs_and_no_receipt_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = pooled_profile(_write_profile(directory, _published_rows()), MODEL_AGES)
        base = _parameters()
        conditional = configure(base, "conditional_mean", profile, scale=2.5, grid=[0., 2., 5., 9., 20.])
        lottery = configure(base, "receipt_lottery", profile, scale=2.5, grid=[0., 2., 5., 9., 20.])
        control = configure(base, "no_receipt", profile, scale=2.5, grid=[0., 2., 5., 9., 20.])
        self.assertEqual(base.__dict__, _parameters().__dict__)
        for name in ("income_transfer", "earnings_level", "preference_marker", "estate_probe_case",
                     "estate_probe_transfer"):
            self.assertEqual(getattr(conditional, name), getattr(base, name))
            self.assertEqual(getattr(lottery, name), getattr(base, name))
            self.assertEqual(getattr(control, name), getattr(base, name))
        for j, row in enumerate(profile):
            conditional_mean = np.dot(conditional._estate_receipt_jump_plans[j].probabilities,
                                      conditional._estate_receipt_jump_plans[j].amounts)
            lottery_mean = np.dot(lottery._estate_receipt_jump_plans[j].probabilities,
                                  lottery._estate_receipt_jump_plans[j].amounts)
            self.assertAlmostEqual(conditional_mean, lottery_mean)
            self.assertAlmostEqual(lottery_mean, 2.5 * row["relative_mean"])
            self.assertTrue(np.isfinite(lottery_mean))
            self.assertGreaterEqual(lottery_mean, 0.)
        self.assertTrue(all(plan.is_identity for plan in control._estate_receipt_jump_plans))
        self.assertTrue(control._estate_receipt_jump_plans[0].is_identity)  # exact entrant control
        self.assertEqual(conditional.estate_receipt_risk_profile, profile)

    def test_configure_rejects_age_misalignment_and_invalid_scale(self):
        profile = [dict(age=float(age), probability=.2, relative_mean=.5,
                        relative_positive_amount=2.5) for age in MODEL_AGES]
        profile[3]["age"] += .5
        with self.assertRaises(ValueError):
            configure(_parameters(), "receipt_lottery", profile, 1., [0., 2., 5.])
        aligned = [dict(age=float(age), probability=.2, relative_mean=.5,
                        relative_positive_amount=2.5) for age in MODEL_AGES]
        with self.assertRaises(ValueError):
            configure(_parameters(), "receipt_lottery", aligned, -1., [0., 2., 5.])
        with self.assertRaises(ValueError):
            configure(_parameters(), "receipt_lottery", aligned, float("nan"), [0., 2., 5.])

    def test_backward_forward_are_transposes_and_unrecorded_transport_is_native_clean(self):
        profile = [dict(age=float(age), probability=(.5 if age == 26 else 0.),
                        relative_mean=(.5 if age == 26 else 0.),
                        relative_positive_amount=(1. if age == 26 else 0.)) for age in MODEL_AGES]
        configured = configure(_parameters(), "receipt_lottery", profile, 1., [0., 1., 2., 4.])
        values = np.arange(8., dtype=float).reshape(4, 2)
        mass = np.array([[.2, .1], [.3, .2], [.1, .1], [0., 0.]])
        backward = _backward(configured, 2, values)
        forward = _forward(configured, 2, mass, record=False)
        self.assertAlmostEqual(float(np.vdot(mass, backward)), float(np.vdot(forward, values)))
        self.assertEqual(configured._estate_receipt_forward_ledger, [])

    def test_normalized_ledger_scales_full_age_mass_and_rejects_bad_ledgers(self):
        raw_mass = np.arange(1., len(MODEL_AGES) + 1.)
        ledger = [dict(age_index=j, transported_mass=float(mass), expected_receipt_flow=float(2 * mass),
                       clipped_wealth_loss=0., wealth_before=float(3 * mass), wealth_after=float(5 * mass))
                  for j, mass in enumerate(raw_mass)]
        solution = SimpleNamespace(g=(2. * raw_mass).reshape(1, 1, 1, len(raw_mass), 1, 1, 1))
        configured = _parameters()
        configured._estate_receipt_forward_ledger = ledger
        result = normalized_ledger(configured, solution)
        self.assertEqual(result["normalization_scale"], 2.)
        self.assertEqual(result["paid_period"], float(np.sum(4. * raw_mass)))
        self.assertEqual([row["transported_mass"] for row in result["by_age"]], list(2. * raw_mass))
        configured._estate_receipt_forward_ledger = ledger[:-1]
        with self.assertRaises(RuntimeError):
            normalized_ledger(configured, solution)
        configured._estate_receipt_forward_ledger = [dict(row, age_index=0) for row in ledger]
        with self.assertRaises(RuntimeError):
            normalized_ledger(configured, solution)
        configured._estate_receipt_forward_ledger = ledger
        wrong_age_mass = np.array(solution.g, copy=True)
        wrong_age_mass[:, :, :, 4, :, :, :] += 1.
        with self.assertRaises(RuntimeError):
            normalized_ledger(configured, SimpleNamespace(g=wrong_age_mass))


if __name__ == "__main__":
    unittest.main()
