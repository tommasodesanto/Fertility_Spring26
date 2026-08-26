from __future__ import annotations

import unittest

import numpy as np

from demographic_transition.person_cohort_law import (
    CohortState,
    advance_one_year,
    cohort_survivors,
    infer_net_migration_path,
    simulate_path,
    solve_demographic_steady_state,
)


class PersonCohortLawTest(unittest.TestCase):
    def test_survival_uses_an_open_last_age_cell(self) -> None:
        persons = np.array([[10.0, 20.0, 30.0, 40.0]])
        births = np.array([5.0])
        survival = np.full_like(persons, 0.5)
        survivors, newborn_deaths, existing_deaths = cohort_survivors(
            persons, births, survival
        )
        np.testing.assert_allclose(survivors, [[2.5, 5.0, 10.0, 35.0]])
        np.testing.assert_allclose(newborn_deaths, [2.5])
        np.testing.assert_allclose(existing_deaths, [50.0])

    def test_one_year_person_and_head_identities_close(self) -> None:
        persons = np.array(
            [[12.0, 10.0, 8.0, 6.0], [11.0, 9.0, 7.0, 5.0]]
        )
        headship_0 = np.array(
            [[0.0, 0.20, 0.55, 0.40], [0.0, 0.15, 0.45, 0.35]]
        )
        state = CohortState(2025, persons, persons * headship_0)
        headship_1 = np.array(
            [[0.0, 0.25, 0.50, 0.30], [0.0, 0.20, 0.50, 0.30]]
        )
        survival = np.full_like(persons, 0.98)
        migration = np.array(
            [[1.0, -0.5, 0.25, 0.0], [0.5, 0.25, -0.2, 0.1]]
        )
        next_state, ledger = advance_one_year(
            state,
            births=np.array([4.0, 3.8]),
            survival=survival,
            net_migration=migration,
            headship_rates=headship_1,
        )
        np.testing.assert_allclose(next_state.heads, next_state.persons * headship_1)
        self.assertLessEqual(abs(ledger.person_identity_residual), 1e-12)
        self.assertLessEqual(abs(ledger.head_identity_residual), 1e-12)
        self.assertLessEqual(ledger.headship_mapping_residual, 1e-12)
        self.assertGreater(float(np.sum(ledger.new_heads_from_nonheads)), 0.0)
        self.assertGreater(float(np.sum(ledger.head_dissolutions)), 0.0)

    def test_inferred_migration_exactly_reproduces_target_path(self) -> None:
        target = {
            2025: np.array([[10.0, 9.0, 8.0], [9.5, 8.5, 7.5]]),
            2026: np.array([[4.2, 10.0, 16.0], [4.0, 9.7, 15.2]]),
            2027: np.array([[4.3, 4.4, 25.0], [4.1, 4.2, 23.5]]),
        }
        births = {2026: np.array([4.0, 3.8]), 2027: np.array([4.1, 3.9])}
        survival = {
            2026: np.full((2, 3), 0.98),
            2027: np.full((2, 3), 0.98),
        }
        migration = infer_net_migration_path(target, births, survival)
        headship = np.array([[0.0, 0.2, 0.5], [0.0, 0.15, 0.45]])
        initial = CohortState(2025, target[2025], target[2025] * headship)
        states, ledgers = simulate_path(
            initial,
            births=births,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            final_year=2027,
        )
        for year in (2025, 2026, 2027):
            np.testing.assert_allclose(states[year].persons, target[year])
            np.testing.assert_allclose(states[year].heads, target[year] * headship)
        self.assertTrue(
            all(abs(row.person_identity_residual) <= 1e-12 for row in ledgers.values())
        )

    def test_policy_births_affect_persons_before_headship(self) -> None:
        persons = np.zeros((2, 22))
        headship = np.zeros_like(persons)
        headship[:, 18:] = 0.25
        survival = {year: np.ones_like(persons) for year in range(2026, 2047)}
        migration = {year: np.zeros_like(persons) for year in range(2026, 2047)}
        baseline_births = {
            year: np.array([10.0, 10.0]) for year in range(2026, 2047)
        }
        reform_births = {
            year: np.array([11.0, 11.0]) for year in range(2026, 2047)
        }
        initial = CohortState(2025, persons, persons)
        baseline, _ = simulate_path(
            initial,
            births=baseline_births,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            final_year=2046,
            topcoded_last_age=False,
        )
        reform, _ = simulate_path(
            initial,
            births=reform_births,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            final_year=2046,
            topcoded_last_age=False,
        )
        self.assertAlmostEqual(
            float(np.sum(reform[2026].persons - baseline[2026].persons)), 2.0
        )
        self.assertAlmostEqual(
            float(np.sum(reform[2026].heads - baseline[2026].heads)), 0.0
        )
        self.assertGreater(
            float(np.sum(reform[2044].heads - baseline[2044].heads)), 0.0
        )

    def test_demographic_steady_state_matches_closed_form_scalar_case(self) -> None:
        survival = np.array([[0.5]])
        headship = np.array([[1.0]])
        migration = np.array([[2.0]])
        result = solve_demographic_steady_state(
            annual_births_per_head=0.25,
            birth_sex_shares=np.array([1.0]),
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
        )
        # P = .5 P + .5 * .25 P + 2, hence P = 16/3.
        self.assertAlmostEqual(float(result.persons[0, 0]), 16.0 / 3.0)
        self.assertAlmostEqual(result.annual_births, 4.0 / 3.0)
        self.assertAlmostEqual(result.birth_head_multiplier, 1.0)
        self.assertAlmostEqual(result.renewal_ratio, 0.25)
        self.assertLessEqual(result.population_fixed_point_residual, 1e-11)

    def test_demographic_steady_state_rejects_explosive_renewal(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "renewal_ratio"):
            solve_demographic_steady_state(
                annual_births_per_head=1.1,
                birth_sex_shares=np.array([1.0]),
                survival=np.array([[0.5]]),
                net_migration=np.array([[1.0]]),
                headship_rates=np.array([[1.0]]),
            )


if __name__ == "__main__":
    unittest.main()
