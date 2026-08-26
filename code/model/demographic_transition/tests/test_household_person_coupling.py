from __future__ import annotations

import unittest

import numpy as np

from demographic_transition import CohortState, advance_household_person_block


class HouseholdPersonCouplingTest(unittest.TestCase):
    def test_person_heads_replace_raw_household_age_masses(self) -> None:
        persons = np.zeros((2, 10))
        persons[:, 2:6] = 10.0
        heads = np.zeros_like(persons)
        heads[:, 2:6] = 5.0
        state = CohortState(2023, persons, heads)
        years = range(2024, 2028)
        survival = {year: np.ones_like(persons) for year in years}
        migration = {year: np.zeros_like(persons) for year in years}
        shares = {year: np.array([0.5, 0.5]) for year in years}
        headship = np.zeros_like(persons)
        headship[:, 2:10] = 0.5

        raw = np.zeros((2, 1, 1, 2, 1))
        raw[:, :, :, 1, :] = 2.0
        template = np.ones(raw[:, :, :, 0, :].shape)
        coupled, next_state, ledger = advance_household_person_block(
            raw,
            state,
            model_period_births=0.0,
            period_years=4,
            birth_sex_shares=shares,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            model_age_start=2,
            model_age_cell_width=4,
            number_of_model_age_cells=2,
            empty_age_templates={0: template},
            topcoded_last_age=False,
        )

        expected = np.array(
            [
                np.sum(next_state.heads[:, 2:6]),
                np.sum(next_state.heads[:, 6:10]),
            ]
        )
        observed = coupled.sum(axis=(0, 1, 2, 4))
        np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0.0)
        self.assertAlmostEqual(ledger.household_person_head_gap, 0.0)
        self.assertAlmostEqual(ledger.coupled_household_mass, expected.sum())


if __name__ == "__main__":
    unittest.main()
