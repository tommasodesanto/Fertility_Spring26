"""Real solver type and tiny allocation-kernel tests; no lifecycle or GE solve.

Run with NUMBA_DISABLE_JIT=1 for the inexpensive local check. Running these same
tests with enabled Numba separately checks compiled kernels; that is not claimed
by the pure local receipt.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]

import e5f_parenthood_utility as adapter
from intergen_eqscale_seq_optimized import kernels, solver
from intergen_eqscale_seq_optimized.e5f_income_entry_profile import E5F_INCOME_ENTRY_DOMAIN
from intergen_eqscale_seq_optimized.parameters import configure_child_state_process, setup_parameters


def old_parameters(*, jump=0.45, slope=0.35, psi=0.0):
    P = setup_parameters()
    overrides = {
        "n_parity": 4, "A_m": 18.0, "child_state_mode": "independent_count",
        "use_stochastic_aging": True, "sequential_births": True,
        "preference_spec": "eqscale", "eqscale_form": "power",
        "alpha_cons": 0.733, "psi_child": psi, "child_room_floor": True,
        "c_bar_0": 0.0, "hbar_first_child_jump": jump,
        "hbar_child_rooms": slope, "beta": 0.9952765791 ** 4,
    }
    # Populate the actual P, then use its actual child-state constructor. The
    # unrelated generic override coercer requires Python 3.10; local /usr/bin
    # Python is 3.9. No utility/type/kernel routine is replaced or mocked here.
    for name, value in overrides.items():
        setattr(P, name, value)
    P.phi = np.full(4, 0.8)
    configure_child_state_process(P)
    return P


def new_parameters(**kwargs):
    return adapter.initialize_parenthood_utility(old_parameters(**kwargs))


def flat(shared, name):
    return np.asarray(getattr(shared, name)).reshape(-1)


def renter_block(P, shared, resources, rent=0.7):
    grid = np.array([0.0, 1.0, 2.0])
    zeros = np.zeros((len(grid), shared.nc))
    return kernels.full_renter_block_kernel(
        Rv1d=np.asarray(resources), Rvt1d=np.asarray(resources), Vc_flat=zeros,
        bp_prev=zeros, has_prev=False, b_grid=grid,
        cb_v=flat(shared, "cb_flat"), hb_v=flat(shared, "hb_flat"),
        psi_v=flat(shared, "psi_flat"), gb_v=flat(shared, "gb_flat"),
        alpha_v=flat(shared, "alpha_flat"), esc_v=flat(shared, "escale_flat"),
        ri=rent, hR_max=P.hR_max, c_min=1e-8, c_bar_0=P.c_bar_0,
        h_bar_0=P.h_bar_0, alpha=P.alpha_cons, oms=1.0-P.sigma,
        beta=P.beta, s_next=1.0, D_next=0.0,
        gs_alpha1=(3.0-np.sqrt(5.0))/2.0,
        gs_alpha2=(np.sqrt(5.0)-1.0)/2.0, gs_tol=1e-10,
        exhaustive_saving=1,
    )


def owner_block(P, shared, resources, housing):
    grid = np.array([0.0, 1.0, 2.0])
    zeros = np.zeros((len(grid), shared.nc))
    return kernels.full_owner_block_kernel(
        Rv1d=np.asarray(resources), Rvt1d=np.asarray(resources), Vco_flat=zeros,
        bp_prev=zeros, has_prev=False, b_grid=grid,
        cb_v=flat(shared, "cb_flat"), hb_v=flat(shared, "hb_flat"),
        psi_v=flat(shared, "psi_flat"), gb_v=flat(shared, "gb_flat"),
        alpha_v=flat(shared, "alpha_flat"), esc_v=flat(shared, "escale_flat"),
        bf_v=np.zeros(shared.nc), oc=0.2, hsv=housing,
        owner_h_bar_scale=P.owner_h_bar_scale, owner_service_premium=P.chi,
        c_min=1e-8, alpha=P.alpha_cons, oms=1.0-P.sigma,
        beta=P.beta, s_next=1.0, D_next=0.0,
        gs_alpha1=(3.0-np.sqrt(5.0))/2.0,
        gs_alpha2=(np.sqrt(5.0)-1.0)/2.0, gs_tol=1e-10,
        strict_hbar_feasibility=1, exhaustive_saving=1,
    )


class ParenthoodUtilityTests(unittest.TestCase):
    def test_initial_first_child_mapping_and_no_unrelated_primitive_rebuild(self):
        old = old_parameters(psi=0.19)
        # A deliberately externally bound fiscal array must survive; generic
        # apply_overrides would recompute it from its other income parameters.
        old.income = np.asarray(old.income) + 123.0
        old.pension = 2.345
        snapshot = copy.deepcopy(vars(old))
        prior = solver.precompute_shared(old, np.array([0.0, 1.0]))
        P = adapter.initialize_parenthood_utility(old)
        shared = solver.precompute_shared(P, np.array([0.0, 1.0]))
        self.assertEqual(adapter.initial_parenthood_requirement(old), 0.8)
        self.assertEqual(P.hbar_first_child_jump, 0.8)
        self.assertEqual(P.hbar_child_rooms, 0.0)
        np.testing.assert_array_equal(shared.h_bar[:, 1], prior.h_bar[:, 1])
        changed = {"hbar_child_rooms", "hbar_first_child_jump", "c_bar_n"}
        self.assertEqual(set(vars(P)), set(snapshot))
        for name, value in snapshot.items():
            np.testing.assert_equal(getattr(old, name), value, err_msg=name)
            if name not in changed:
                np.testing.assert_equal(getattr(P, name), value, err_msg=name)
        self.assertIsNot(P.income, old.income)
        self.assertIsNot(P.Pi_child, old.Pi_child)

    def test_actual_type_arrays_and_scale_for_all_current_children(self):
        for psi in (0.0, 0.19):
            P = new_parameters(psi=psi)
            shared = solver.precompute_shared(P, np.array([0.0, 1.0]))
            for parity in range(4):
                for m in range(parity + 1):
                    with self.subTest(psi=psi, parity=parity, m=m):
                        col = parity + 4*m
                        self.assertEqual(shared.h_bar[parity, m], 0.8 if m else 0.0)
                        self.assertEqual(shared.c_bar[parity, m], 0.0)
                        self.assertEqual(shared.psi_v[parity, m], psi*m)
                        self.assertEqual(shared.alpha_flat[0, col], P.alpha_cons)
                        expected = ((2.0+0.7*m)/2.0)**0.7
                        self.assertAlmostEqual(shared.escale_flat[0, col], expected, places=15)
            if psi == 0.0:
                # Legacy triple deduplication really DOES collapse parents;
                # kernels must retain full per-state equivalence scales.
                self.assertEqual(shared.n_types, 2)
                self.assertEqual(len(set(shared.type_map[[7, 11, 15]])), 1)
                self.assertEqual(len(set(shared.escale_flat[0, [7, 11, 15]])), 3)

    def test_actual_maturation_transition_reaches_zero_floor_and_reward(self):
        P = new_parameters(psi=0.19)
        shared = solver.precompute_shared(P, np.array([0.0, 1.0]))
        for parity in (1, 2, 3):
            self.assertGreater(P.Pi_child[parity, parity, 0], 0.0)
            self.assertEqual(shared.h_bar[parity, 0], 0.0)
            self.assertEqual(shared.psi_v[parity, 0], 0.0)
            self.assertEqual(shared.escale_flat[0, parity], 1.0)

    def test_annual_beta_power_four_and_preserved_fiscal_arrays(self):
        P = new_parameters()
        P.income = np.asarray(P.income) + 123.0
        P.pension = 2.345
        old_beta = P.beta
        annual = 0.98
        new = adapter.bind_parenthood_utility(P, {"beta_annual": annual, "H0": 6.0,
                                                 "kappa_fert": 0.1})
        self.assertEqual(new.beta, annual**4)
        self.assertEqual(new.rho, 1.0/new.beta-1.0)
        self.assertEqual(new.rho_hat, new.rho)
        self.assertEqual(new.kappa_fert, new.eps_fert)
        self.assertEqual(P.beta, old_beta)
        self.assertEqual(new.pension, P.pension)
        self.assertEqual(np.asarray(new.H0).shape, np.asarray(P.H0).shape)
        np.testing.assert_array_equal(new.income, P.income)

    def test_reload_and_candidate_slope_changes_fail_exactly(self):
        P = new_parameters()
        reloaded = SimpleNamespace(**copy.deepcopy(vars(P)))
        adapter.validate_parenthood_utility(reloaded)
        self.assertEqual(adapter.bind_parenthood_utility(reloaded).hbar_first_child_jump, 0.8)
        for slope in (1e-300, 0.2, -1e-300, np.nan):
            broken = copy.deepcopy(P)
            broken.hbar_child_rooms = slope
            with self.subTest(slope=slope), self.assertRaises(ValueError):
                adapter.bind_parenthood_utility(broken, {"h_P": 0.9})
        for candidate in ({"hbar_child_rooms": 0.0}, {"hbar_child_rooms": 0.2},
                          {"hbar_first_child_jump": 0.8}, {"psi_child": 0.2},
                          {"sigma": 2.0}, {"alpha_cons": 0.733}):
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                adapter.validate_parenthood_candidate(candidate)

    def test_fixed_utility_and_lifecycle_corruption_fails(self):
        for name, value in {"sigma": 1.5, "alpha_cons": 0.70, "delta_alpha": 0.01,
                            "delta_alpha_jump": 0.01, "c_bar_0": 0.1, "c_bar_n": 0.1,
                            "eqscale_form": "sqrt", "preference_spec": "stone_geary",
                            "child_room_floor": False, "period_years": 1.0,
                            "child_state_mode": "shared_clock"}.items():
            P = new_parameters()
            setattr(P, name, value)
            with self.subTest(field=name), self.assertRaises(ValueError):
                adapter.validate_parenthood_utility(P)

    def test_nine_coordinates_preserve_existing_bounds_and_transforms(self):
        domain = adapter.PARENTHOOD_SEARCH_DOMAIN
        self.assertEqual(adapter.PARENTHOOD_SEARCH_NAMES, (
            "beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
            "theta0", "theta1", "first_birth_fixed_cost", "h_P"))
        inherited = {row[0]: row for row in E5F_INCOME_ENTRY_DOMAIN}
        for row in domain[:-1]:
            self.assertEqual(row, inherited[row[0]])
        self.assertEqual(domain[-1], ("h_P", 0.1, 2.3, "log"))
        candidate = {name: (lo+hi)/2 for name, lo, hi, _ in domain}
        adapter.validate_parenthood_candidate(candidate, require_complete=True)
        with self.assertRaises(ValueError):
            adapter.validate_parenthood_candidate({"h_P": 0.8}, require_complete=True)
        for candidate in ({"h_P": 2.31}, {"beta_annual": np.nan}, {"H0": [1, 2]}):
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                adapter.validate_parenthood_candidate(candidate)
        self.assertEqual(json.loads(json.dumps(adapter.parenthood_utility_metadata()))["free_parameter_count"], 9)

    def test_h_p_is_translated_and_infeasible_in_place_bind_is_atomic(self):
        P = new_parameters()
        new = adapter.bind_parenthood_utility(P, {"h_P": 0.9})
        self.assertEqual(new.hbar_first_child_jump, 0.9)
        self.assertFalse(hasattr(new, "h_P"))
        self.assertEqual(P.hbar_first_child_jump, 0.8)
        self.assertIs(adapter.bind_parenthood_utility(P, {"h_P": 0.9}, copy_parameters=False), P)
        P.hR_max = 1.0
        with self.assertRaises(ValueError):
            adapter.bind_parenthood_utility(P, {"h_P": 1.2, "theta0": 0.9}, copy_parameters=False)
        self.assertEqual(P.hbar_first_child_jump, 0.9)

    def test_actual_renter_allocation_foc_crra_scale_and_cap(self):
        rent = 0.7
        resources = np.array([2.0, 8.0, 30.0])
        for psi in (0.0, 0.19):
            P = new_parameters(psi=psi)
            shared = solver.precompute_shared(P, np.array([0.0, 1.0, 2.0]))
            value, saving, consumption, housing = renter_block(P, shared, resources, rent)
            for m in range(4):
                col, floor = 3+4*m, 0.8 if m else 0.0
                scale = ((2.0+0.7*m)/2.0)**0.7
                X = resources-saving[:, col]
                unconstrained = (1.0-P.alpha_cons)*X/rent+P.alpha_cons*floor
                np.testing.assert_allclose(housing[:, col], np.minimum(unconstrained, P.hR_max), rtol=0, atol=1e-12)
                np.testing.assert_allclose(consumption[:, col]+rent*housing[:, col], X, rtol=0, atol=1e-12)
                np.testing.assert_array_equal(saving[:, col], 0.0)
                composite = consumption[:, col]**P.alpha_cons*(housing[:, col]-floor)**(1.0-P.alpha_cons)
                np.testing.assert_allclose(value[:, col], -scale/composite+psi*m, rtol=0, atol=1e-12)
                for row in (0, 1):
                    mrs = (1.0-P.alpha_cons)/P.alpha_cons*consumption[row, col]/(housing[row, col]-floor)
                    self.assertAlmostEqual(mrs, rent, places=12)
                self.assertGreater(unconstrained[-1], P.hR_max)
            # At fixed resources, current parents have identical housing, even
            # though their effective material utilities differ.
            np.testing.assert_allclose(housing[:, 7], housing[:, 11], rtol=0, atol=1e-12)
            np.testing.assert_allclose(housing[:, 11], housing[:, 15], rtol=0, atol=1e-12)

    def test_actual_renter_affordability_boundary_is_infeasible(self):
        P = new_parameters()
        shared = solver.precompute_shared(P, np.array([0.0, 1.0, 2.0]))
        floor_cost = 0.7*P.hbar_first_child_jump
        value, _, _, _ = renter_block(P, shared, np.array([floor_cost-0.01, floor_cost, floor_cost+0.01]))
        self.assertTrue(np.all(value[:2, [7, 11, 15]] == -1e10))
        self.assertTrue(np.all(value[-1, [7, 11, 15]] > -1e9))

    def test_actual_owner_kernel_retains_scale_and_strict_floor_boundary(self):
        resources = np.array([1.0, 5.0, 12.0])
        for psi in (0.0, 0.19):
            P = new_parameters(psi=psi)
            shared = solver.precompute_shared(P, np.array([0.0, 1.0, 2.0]))
            value, saving, consumption = owner_block(P, shared, resources, housing=2.0)
            for m in range(4):
                col, floor = 3+4*m, 0.8 if m else 0.0
                scale = ((2.0+0.7*m)/2.0)**0.7
                composite = consumption[:, col]**P.alpha_cons*(P.chi*(2.0-floor))**(1.0-P.alpha_cons)
                np.testing.assert_allclose(value[:, col], -scale/composite+psi*m, rtol=0, atol=1e-12)
                np.testing.assert_array_equal(saving[:, col], 0.0)
                np.testing.assert_allclose(consumption[:, col]+saving[:, col]+0.2, resources, rtol=0, atol=1e-12)
            invalid, _, _ = owner_block(P, shared, resources, housing=0.8)
            self.assertTrue(np.all(invalid[:, [7, 11, 15]] == -1e10))
            self.assertTrue(np.all(invalid[:, 3] > -1e9))


if __name__ == "__main__":
    unittest.main()
