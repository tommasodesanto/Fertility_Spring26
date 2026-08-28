#!/usr/bin/env python3
"""Focused tests for the isolated tenure-complementarity oracle."""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_person_policy_tenure_complementarity_oracle as oracle


class ProbabilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.choices = np.zeros((2, 3, 1, 2, 1, 1, 1), dtype=np.int16)
        self.kwargs = dict(
            wealth_index=0,
            origin_tenure=0,
            market=0,
            age_index=0,
            child_state=0,
            renter_tenure=0,
            owner_tenure=2,
        )

    def test_default_one_hot_probabilities_reproduce_choices(self) -> None:
        probabilities = oracle.deterministic_tenure_probabilities(self.choices)
        np.testing.assert_array_equal(np.argmax(probabilities, axis=-1), self.choices)
        np.testing.assert_array_equal(np.sum(probabilities, axis=-1), 1.0)

    def test_only_declared_atom_is_mixed(self) -> None:
        probabilities, record = oracle.probabilities_with_atom_share(
            self.choices, owner_share=0.25, **self.kwargs
        )
        self.assertEqual(record["strictly_mixed_probability_rows"], 1)
        np.testing.assert_array_equal(probabilities[0, 0, 0, 0, 0, 0, 0], [0.75, 0.0, 0.25])
        np.testing.assert_array_equal(probabilities[1, 0, 0, 0, 0, 0, 0], [1.0, 0.0, 0.0])

    def test_endpoint_shares_are_pure_renter_and_owner(self) -> None:
        renter, _ = oracle.probabilities_with_atom_share(
            self.choices, owner_share=0.0, **self.kwargs
        )
        owner, _ = oracle.probabilities_with_atom_share(
            self.choices, owner_share=1.0, **self.kwargs
        )
        self.assertEqual(np.argmax(renter[0, 0, 0, 0, 0, 0, 0]), 0)
        self.assertEqual(np.argmax(owner[0, 0, 0, 0, 0, 0, 0]), 2)

    def test_complementarity_conditions(self) -> None:
        self.assertEqual(oracle.complementarity_violation(0.0, np.array([-1.0])), 0.0)
        self.assertEqual(oracle.complementarity_violation(1.0, np.array([1.0])), 0.0)
        self.assertEqual(oracle.complementarity_violation(0.5, np.array([0.0])), 0.0)
        self.assertEqual(oracle.complementarity_violation(0.5, np.array([2e-8])), 2e-8)

    def test_mixing_atom_accepts_live_renter_owner_margin(self) -> None:
        options = np.zeros((2, 3, 3))
        options[..., 0] = 1.0
        options[..., 1] = 0.5
        options[..., 2] = 1.0 + 2e-8
        choices = np.full((2, 3), 2, dtype=np.int16)
        gaps = oracle.validate_mixing_atom_options(
            options,
            choices,
            renter_tenure=0,
            owner_tenure=2,
        )
        np.testing.assert_allclose(gaps, 2e-8, rtol=0.0, atol=2e-16)

    def test_mixing_atom_rejects_infeasible_branch(self) -> None:
        options = np.zeros((2, 3, 3))
        options[..., 0] = 1.0
        options[..., 2] = 1.0
        options[1, 2, 2] = oracle.kink.NEGATIVE_INFINITY
        choices = np.zeros((2, 3), dtype=np.int16)
        with self.assertRaisesRegex(RuntimeError, "infeasible renter/owner branch"):
            oracle.validate_mixing_atom_options(
                options,
                choices,
                renter_tenure=0,
                owner_tenure=2,
            )

    def test_mixing_atom_rejects_third_tenure_argmax(self) -> None:
        options = np.ones((2, 3, 3))
        choices = np.zeros((2, 3), dtype=np.int16)
        choices[0, 1] = 1
        with self.assertRaisesRegex(RuntimeError, "third tenure"):
            oracle.validate_mixing_atom_options(
                options,
                choices,
                renter_tenure=0,
                owner_tenure=2,
            )

    def test_current_housing_and_tax_are_linear_in_atom_share(self) -> None:
        g_post = np.zeros((2, 3, 1, 2, 1, 1, 1))
        g_post[0, 0, 0, 0, 0, 0, 0] = 1.0
        loc_probs = np.ones((2, 3, 1, 1, 2, 1, 1, 1))
        lmm_idx = np.zeros((1, 3, 2), dtype=int)
        lmm_wt = np.zeros((1, 3, 2))
        tmx_idx = np.zeros((1, 3, 3, 1, 1, 2), dtype=int)
        tmx_wt = np.zeros((1, 3, 3, 1, 1, 2))
        hR_pol = np.zeros_like(g_post)
        hR_pol[:, 0, ...] = 1.0
        P = SimpleNamespace(
            I=1,
            n_house=2,
            H_own=np.array([2.0, 4.0]),
            tau_H=0.02,
            use_numba_scatter=False,
        )

        outcomes = {}
        for share in (0.0, 0.3, 1.0):
            probabilities, _ = oracle.probabilities_with_atom_share(
                self.choices, owner_share=share, **self.kwargs
            )
            current = calendar.model.realize_current_cross_section(
                g_post,
                loc_probs,
                self.choices,
                probabilities,
                lmm_idx,
                lmm_wt,
                tmx_idx,
                tmx_wt,
            )
            demand = float(np.sum(calendar.housing_demand_by_location(current, hR_pol, P)))
            tax = calendar.model.property_tax_revenue_from_distribution(
                current, hR_pol, np.array([0.5]), P
            )
            outcomes[share] = (current, demand, tax)
        np.testing.assert_allclose(
            outcomes[0.3][0], 0.7 * outcomes[0.0][0] + 0.3 * outcomes[1.0][0]
        )
        self.assertAlmostEqual(outcomes[0.3][1], 0.7 * outcomes[0.0][1] + 0.3 * outcomes[1.0][1])
        self.assertAlmostEqual(outcomes[0.3][2], 0.7 * outcomes[0.0][2] + 0.3 * outcomes[1.0][2])

    def test_forward_distribution_uses_the_same_atom_share(self) -> None:
        cohort = np.zeros((2, 3, 1, 1, 1, 1))
        cohort[0, 0, 0, 0, 0, 0] = 1.0
        loc_probs = np.ones((2, 3, 1, 1, 2, 1, 1, 1))
        lmm_idx = np.zeros((1, 3, 2), dtype=int)
        lmm_wt = np.zeros((1, 3, 2))
        tmx_idx = np.zeros((1, 3, 3, 1, 1, 2), dtype=int)
        tmx_wt = np.zeros((1, 3, 3, 1, 1, 2))
        bp_pol = np.zeros_like(self.choices, dtype=float)
        P = SimpleNamespace(
            I=1,
            n_house=2,
            n_parity=1,
            n_child_states=1,
            n_child_stages=0,
            use_numba_scatter=False,
            readiness_gate_enabled=False,
        )
        shared = SimpleNamespace(nc=1)
        outputs = {}
        for share in (0.0, 0.3, 1.0):
            probabilities, _ = oracle.probabilities_with_atom_share(
                self.choices, owner_share=share, **self.kwargs
            )
            outputs[share] = calendar.model.advance_cohort_one_period_markov_income(
                cohort,
                0,
                loc_probs,
                self.choices,
                probabilities,
                bp_pol,
                P,
                np.array([0.0, 1.0]),
                shared,
                lmm_idx,
                lmm_wt,
                tmx_idx,
                tmx_wt,
                False,
                None,
                np.ones((1, 1)),
            )
        np.testing.assert_allclose(
            outputs[0.3], 0.7 * outputs[0.0] + 0.3 * outputs[1.0]
        )


if __name__ == "__main__":
    unittest.main()
