from __future__ import annotations

import unittest

import numpy as np

from demographic_transition import CohortState, advance_person_state_block


class FourYearPersonBridgeTest(unittest.TestCase):
    def test_block_spreads_births_and_closes_every_annual_identity(self) -> None:
        persons = np.zeros((2, 8))
        headship = np.zeros_like(persons)
        headship[:, 4:] = 0.5
        state = CohortState(2025, persons, persons)
        years = range(2026, 2030)
        survival = {year: np.ones_like(persons) for year in years}
        migration = {year: np.zeros_like(persons) for year in years}
        shares = {year: np.array([0.51, 0.49]) for year in years}
        result, ledger = advance_person_state_block(
            state,
            model_period_births=40.0,
            period_years=4,
            birth_sex_shares=shares,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            topcoded_last_age=False,
        )
        self.assertEqual(result.year, 2029)
        self.assertAlmostEqual(float(np.sum(result.persons)), 40.0)
        self.assertAlmostEqual(float(np.sum(result.heads)), 0.0)
        self.assertAlmostEqual(
            sum(float(np.sum(row.births)) for row in ledger.annual_ledgers), 40.0
        )
        self.assertLessEqual(ledger.person_identity_max_abs, 1e-12)
        self.assertLessEqual(ledger.head_identity_max_abs, 1e-12)

    def test_block_rejects_incomplete_calendar_inputs(self) -> None:
        state = CohortState(2025, np.zeros((1, 3)), np.zeros((1, 3)))
        with self.assertRaisesRegex(KeyError, "2027"):
            advance_person_state_block(
                state,
                model_period_births=1.0,
                period_years=2,
                birth_sex_shares={2026: np.array([1.0])},
                survival={2026: np.ones((1, 3))},
                net_migration={2026: np.zeros((1, 3))},
                headship_rates=np.zeros((1, 3)),
            )


if __name__ == "__main__":
    unittest.main()
