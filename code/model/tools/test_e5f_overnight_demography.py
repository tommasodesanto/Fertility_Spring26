from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from e5f_overnight_demography import (
    FormationUnitContract,
    OutsideEntryContract,
    advance_surviving_maturation_distribution,
    assert_zero_post_origin_migration,
    zero_post_origin_migration,
)


def _distribution() -> np.ndarray:
    # wealth, tenure, location, age, income, children-ever-born, children-at-home
    g = np.zeros((1, 1, 1, 3, 1, 4, 4))
    g[0, 0, 0, 0, 0, 2, 2] = 10.0
    g[0, 0, 0, 1, 0, 1, 1] = 4.0
    g[0, 0, 0, 2, 0, 3, 3] = 2.0
    return g


class OvernightDemographyTests(unittest.TestCase):
    def setUp(self) -> None:
        Pi_child = np.zeros((4, 4, 4))
        for n in range(4):
            for m in range(n + 1):
                Pi_child[m, max(m - 1, 0), n] = 1.0
        self.P = SimpleNamespace(
            I=1,
            use_age_survival=True,
            survival_probs=np.array([0.8, 0.5]),
            use_stochastic_aging=True,
            child_state_mode="independent_count",
            Pi_child=Pi_child,
        )
        self.evaluation = SimpleNamespace(g_post_fertility=_distribution())
        self.contract = FormationUnitContract(
            mature_persons_per_model_child_unit=1.1,
            household_heads_per_mature_person=0.5,
            classification="externally_fixed",
            source="synthetic test contract",
            model_child_unit_definition="capped children currently at home",
            mature_person_definition="surviving dependents leaving home",
            household_head_definition="new lifecycle decision unit",
        )

    @staticmethod
    def _advance(evaluation, entry, P, b_grid, shared):
        del entry, b_grid, shared
        current = evaluation.g_post_fertility
        nxt = np.zeros_like(current)
        mature = np.zeros(P.I)
        deaths = 0.0
        for age, survival in enumerate(P.survival_probs):
            cohort = current[:, :, :, age, :, :, :]
            deaths += (1.0 - survival) * np.sum(cohort)
            survivors = survival * cohort
            # Deterministic stand-in: one child matures in every surviving family.
            for n in range(1, current.shape[5]):
                for m in range(1, min(n, current.shape[6] - 1) + 1):
                    mass = survivors[:, :, :, :, n, m]
                    nxt[:, :, :, age + 1, :, n, m - 1] += mass
                    mature += np.sum(mass, axis=(0, 1, 3))
            nxt[:, :, :, age + 1, :, 0, 0] += survivors[:, :, :, :, 0, 0]
        deaths += np.sum(current[:, :, :, -1, :, :, :])
        return nxt, mature, float(deaths), 0.0

    @staticmethod
    def _entrants(entry, P, b_grid):
        del b_grid
        out = np.zeros((1, 1, P.I, 1, 4, 4))
        out[0, 0, :, 0, 0, 0] = entry
        return out

    def test_b0_joint_death_surviving_maturation_and_identities(self) -> None:
        nxt, entry, ledger = advance_surviving_maturation_distribution(
            self.evaluation,
            self.P,
            np.array([0.0]),
            None,
            advance_distribution=self._advance,
            entrant_cohort=self._entrants,
            formation_contract=self.contract,
            zero_migration=True,
            production=True,
        )
        # Survivors: age-0 8 families and age-1 2 families; one modeled child
        # matures per surviving family. Terminal-age dependents die jointly.
        self.assertAlmostEqual(ledger.surviving_maturations_model_child_units, 10.0)
        self.assertAlmostEqual(ledger.domestic_formation_entries, 5.5)
        self.assertAlmostEqual(entry[0], 5.5)
        self.assertEqual(ledger.outside_migration_entries, 0.0)
        self.assertAlmostEqual(ledger.child_identity_residual, 0.0)
        self.assertAlmostEqual(ledger.household_identity_residual, 0.0)
        self.assertAlmostEqual(np.sum(nxt), ledger.next_households)

    def test_production_rejects_diagnostic_conversion(self) -> None:
        diagnostic = FormationUnitContract(
            **{**self.contract.__dict__, "classification": "diagnostic"}
        )
        with self.assertRaisesRegex(ValueError, "cannot be promoted"):
            advance_surviving_maturation_distribution(
                self.evaluation,
                self.P,
                np.array([0.0]),
                None,
                advance_distribution=self._advance,
                entrant_cohort=self._entrants,
                formation_contract=diagnostic,
                production=True,
            )

    def test_diagnostic_conversion_is_explicit_and_labelled(self) -> None:
        diagnostic = FormationUnitContract(
            **{**self.contract.__dict__, "classification": "diagnostic"}
        )
        _, entry, ledger = advance_surviving_maturation_distribution(
            self.evaluation,
            self.P,
            np.array([0.0]),
            None,
            advance_distribution=self._advance,
            entrant_cohort=self._entrants,
            formation_contract=diagnostic,
            production=False,
        )
        self.assertGreater(entry[0], 0.0)
        self.assertEqual(ledger.formation_contract_classification, "diagnostic")

    def test_b0_rejects_offsetting_outside_entries(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly zero"):
            advance_surviving_maturation_distribution(
                self.evaluation,
                self.P,
                np.array([0.0]),
                None,
                advance_distribution=self._advance,
                entrant_cohort=self._entrants,
                formation_contract=self.contract,
                outside_entries_by_location=np.array([1.0]),
                zero_migration=True,
            )

    def test_b_plus_rejects_uncontracted_mapped_entries(self) -> None:
        with self.assertRaisesRegex(ValueError, "state-allocation contract"):
            advance_surviving_maturation_distribution(
                self.evaluation,
                self.P,
                np.array([0.0]),
                None,
                advance_distribution=self._advance,
                entrant_cohort=self._entrants,
                formation_contract=self.contract,
                outside_entries_by_location=np.array([1.0]),
                zero_migration=False,
            )
        migration_contract = OutsideEntryContract(
            classification="externally_fixed",
            source="synthetic mapped flow",
            input_unit_definition="net migrant household heads",
            household_state_allocation_definition="already allocated to model location",
        )
        _, entries, ledger = advance_surviving_maturation_distribution(
            self.evaluation,
            self.P,
            np.array([0.0]),
            None,
            advance_distribution=self._advance,
            entrant_cohort=self._entrants,
            formation_contract=self.contract,
            outside_entries_by_location=np.array([1.0]),
            outside_entry_contract=migration_contract,
            zero_migration=False,
        )
        self.assertAlmostEqual(entries[0], 6.5)
        self.assertAlmostEqual(ledger.outside_migration_entries, 1.0)

    def test_zero_switch_preserves_history_and_zeros_every_forecast_cell(self) -> None:
        migration = {2022: np.array([[1.0, -1.0]]), 2023: np.array([[2.0, -2.0]])}
        outside = {2022: np.array([0.2]), 2027: np.array([0.3])}
        zero_migration, zero_outside = zero_post_origin_migration(migration, outside)
        np.testing.assert_array_equal(zero_migration[2022], migration[2022])
        np.testing.assert_array_equal(zero_migration[2023], np.zeros((1, 2)))
        np.testing.assert_array_equal(zero_outside[2027], np.zeros(1))
        assert_zero_post_origin_migration(zero_migration, zero_outside)
        with self.assertRaisesRegex(ValueError, "at least one cell"):
            assert_zero_post_origin_migration(migration, zero_outside)


if __name__ == "__main__":
    unittest.main()
