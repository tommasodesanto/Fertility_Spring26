"""Pure regression tests for dated PF choice ownership; no model solves.

The Bellman boundary is stubbed with its observed allocation/replacement
contract. Actual joint allocation and population factorization are exercised.
"""
from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import run_e5f_perfect_foresight_transition as pf
from intergen_eqscale_seq_optimized import joint_nested


SHAPE = (2, 2, 1, 1, 1, 2, 2)


def parameters(nested=True):
    return SimpleNamespace(
        I=1, J=1, n_house=1, n_parity=2, n_child_states=2,
        sequential_births=True, child_maturation_mode="independent",
        joint_nested_choice=nested, fertility_nest_choice=nested,
        two_shock_choice=False, tenure_choice_kappa=.005,
        kappa_fert=.02, kappa_fert_continuation=.03,
        A_f_start=1, A_f_end=1, age_start=30., da=4.,
        fecundity_omega1=.6, fecundity_omega2=0.,
    )


def bellman_outputs(value):
    # Arrays returned directly by the real solver are allocated per call.
    arrays = [np.full(SHAPE, value) for _ in range(9)]
    return (*arrays, None)


def choice_object(P, attempt):
    joint = joint_nested.allocate(SHAPE, P)
    joint.products[...] = np.arange(2)
    joint.wait_probabilities[..., 0] = .75
    joint.wait_probabilities[..., 1] = .25
    joint.probabilities[..., 0, 0] = (1 - attempt) * .75
    joint.probabilities[..., 1, 0] = (1 - attempt) * .25
    joint.probabilities[..., 0, 1] = attempt * .25
    joint.probabilities[..., 1, 1] = attempt * .75
    joint.failure_probabilities[..., 0] = attempt * .8
    joint.failure_probabilities[..., 1] = attempt * .2
    return joint


class DatedChoiceOwnershipTests(unittest.TestCase):
    def setUp(self):
        self.maps = object()
        self.map_patch = patch.object(pf.calendar, "build_transition_maps", return_value=self.maps)
        self.map_patch.start()
        self.addCleanup(self.map_patch.stop)

    def test_nested_policy_retains_original_operator_after_next_date(self):
        P = parameters()
        joint_objects = []
        continuation = np.zeros(SHAPE)

        def boundary(rent, price, parameters, grid, shared, *, continuation_V):
            self.assertIs(continuation_V, continuation)
            # The real solver allocates a new joint and fert2 array on every
            # call and assigns both attributes only after solving its ages.
            joint = choice_object(parameters, .6 if not joint_objects else .2)
            joint_objects.append(joint)
            parameters._joint_choice = joint
            parameters._fert2_probs = np.full((2, 2), len(joint_objects) / 10)
            return bellman_outputs(float(len(joint_objects)))

        with patch.object(pf.calendar.model, "solve_bellman_full_markov_income", side_effect=boundary):
            first = pf.solve_date_policy(price=1., rent=.1, P=P, b_grid=np.arange(2.),
                                         shared=SimpleNamespace(), continuation_V=continuation)
            inherited = np.zeros(SHAPE)
            inherited[0, 0, 0, 0, 0, 0, 0] = 1.
            before = pf.calendar.factor_joint_distribution(inherited, first, P)
            second = pf.solve_date_policy(price=1.1, rent=.12, P=P, b_grid=np.arange(2.),
                                          shared=SimpleNamespace(), continuation_V=continuation)
        self.assertIs(first.joint_choice, joint_objects[0])
        self.assertIs(second.joint_choice, joint_objects[1])
        self.assertIs(P._joint_choice, second.joint_choice)
        for name in ("probabilities", "products", "wait_probabilities", "failure_probabilities"):
            self.assertFalse(np.shares_memory(getattr(first.joint_choice, name),
                                             getattr(second.joint_choice, name)))
        # Mutating the later cache must not alter the earlier dated operator.
        P._joint_choice.probabilities.fill(0.)
        P._joint_choice.failure_probabilities.fill(0.)
        P._fert2_probs.fill(.9)
        after = pf.calendar.factor_joint_distribution(inherited, first, P)
        for expected, actual in zip(before, after):
            np.testing.assert_array_equal(actual, expected)
        post, effective, births, attempts, risk = after
        childless = (0, 0, 0, 0, 0, 0, 0)
        parent = (0, 0, 0, 0, 0, 1, 1)
        self.assertAlmostEqual(post[childless], .76)
        self.assertAlmostEqual(post[parent], .24)
        np.testing.assert_allclose(effective[childless], np.array([.588, .172]) / .76)
        np.testing.assert_allclose(effective[parent], [.25, .75])
        self.assertAlmostEqual(births.sum(), .24)
        self.assertAlmostEqual(attempts.sum(), .6)
        self.assertAlmostEqual(risk.sum(), 1.)
        self.assertAlmostEqual(post.sum(), inherited.sum())
        np.testing.assert_array_equal(first.fert2_probs, np.full((2, 2), .1))
        np.testing.assert_array_equal(first.V, np.ones(SHAPE))
        np.testing.assert_array_equal(first.price, [1.])

    def test_sequential_continuation_probabilities_are_owned_snapshots(self):
        P = parameters(nested=False)
        cache = np.array([[.8, .2], [.3, .7]])
        P._fert2_probs = cache
        P._joint_choice = None
        first = pf.policy_from_objects(bellman_outputs(1.), 1., P, np.arange(2.), SimpleNamespace())
        expected = cache.copy()
        # Even in-place reuse of this mutable parameter cache cannot corrupt
        # continuation-birth probabilities belonging to an earlier date.
        cache[:] = [[.1, .9], [.6, .4]]
        second = pf.policy_from_objects(bellman_outputs(2.), 2., P, np.arange(2.), SimpleNamespace())
        P._fert2_probs = np.zeros_like(cache)
        np.testing.assert_array_equal(first.fert2_probs, expected)
        np.testing.assert_array_equal(second.fert2_probs, cache)
        self.assertFalse(np.shares_memory(first.fert2_probs, cache))
        self.assertFalse(np.shares_memory(second.fert2_probs, cache))
        self.assertIsNone(first.joint_choice)
        self.assertIsNone(second.joint_choice)

    def test_policy_without_optional_choice_caches_preserves_legacy_default(self):
        P = parameters(nested=False)
        policy = pf.policy_from_objects(bellman_outputs(1.), 1., P, np.arange(2.), SimpleNamespace())
        self.assertIsNone(policy.fert2_probs)
        self.assertIsNone(policy.joint_choice)
        self.assertIs(policy.maps, self.maps)


if __name__ == "__main__":
    unittest.main()
