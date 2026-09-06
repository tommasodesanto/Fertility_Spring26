"""Small plumbing checks for dated joint-nested policy snapshots."""
from __future__ import annotations

from types import SimpleNamespace as NS
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import run_dynamic_population_transition as calendar


class JointNestedIntegrationTests(TestCase):
    def setUp(self) -> None:
        self.P = NS(
            joint_nested_choice=True, I=1, J=1, n_house=1, n_parity=2,
            n_child_states=2, H_own=np.array([1.0]), H0=np.array([1.0]),
            user_cost_rate=1.0, r_bar=np.array([1.0]), xi_supply=np.array([1.0]),
        )
        self.g = np.zeros((1, 2, 1, 1, 1, 2, 2))
        self.g[0, 0, 0, 0, 0, 0, 0] = 1.0
        self.maps = calendar.TransitionMaps(None, None, None, None)
        self.joint = object()
        self.policy = calendar.PolicyBundle(
            V=np.ones_like(self.g), c_pol=np.ones_like(self.g),
            hR_pol=np.zeros_like(self.g), bp_pol=np.zeros_like(self.g),
            tenure_choice=np.zeros_like(self.g, dtype=int),
            tenure_probs=np.zeros(self.g.shape + (2,)), loc_probs=np.ones_like(self.g),
            fert_probs=np.zeros_like(self.g), fert_value=np.zeros(self.g.shape[:-2]),
            price=np.array([1.0]), maps=self.maps, joint_choice=self.joint,
        )

    def test_policy_bundle_owns_joint_choice(self) -> None:
        self.assertIs(self.policy.joint_choice, self.joint)

    def test_evaluation_uses_distribution_specific_kernel(self) -> None:
        effective = np.full(self.g.shape + (2,), 0.5)
        births = np.array([[0.25, 0.0]])

        def factor(g_pre, policy, P, *, mode="natural"):
            self.assertEqual(mode, "natural")
            self.assertIs(policy.joint_choice, self.joint)
            return g_pre.copy(), effective, births, births, births

        with (
            patch.object(calendar, "factor_joint_distribution", side_effect=factor),
            patch.object(calendar, "gate_pre_fertility_distribution", side_effect=lambda g, *x: (g, 0.0)),
            patch.object(calendar.model, "realize_current_cross_section", side_effect=lambda g, *x, **k: g),
        ):
            evaluation = calendar.evaluate_period(
                np.array([1.0]), self.g, self.P, np.array([0.0]), NS(),
                calendar.SolveCounter(), supplied_policy=self.policy,
            )
        self.assertIs(evaluation.policy.joint_choice, self.joint)
        self.assertIs(evaluation.policy.tenure_probs, effective)
        self.assertAlmostEqual(evaluation.births, 0.25)


if __name__ == "__main__":
    main()
