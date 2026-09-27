"""Native array checks for the proposed utility; run on Torch."""
import copy
import math
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from intergen_eqscale_seq_optimized import solver
from intergen_eqscale_seq_optimized.child_preferences import apply_child_preferences


def parameters():
    return SimpleNamespace(n_parity=4, n_child_states=4, n_child_stages=4, n_house=1,
        child_state_mode='independent_count', n_child_count=4, alpha_cons=.733, sigma=2.,
        c_bar_0=0., c_bar_n=0., h_bar_0=0., h_bar_n=0., h_bar_jump=0.,
        psi_child=.06, eqscale_form='power', child_room_floor=False, hbar_child_rooms=0.,
        hbar_first_child_jump=0., preference_spec='eqscale', child_housing_spec='linear_only',
        delta_alpha_jump=.2, delta_alpha=0., gamma_e=0.,
        child_benefit_curvature=.14, compensated_child_housing_shares=True,
        utility_reference_rent=.1104659270)


def shared(P):
    with patch.object(solver, 'independent_child_maturation_active', return_value=True), \
         patch.object(solver, 'readiness_gate_active', return_value=False), \
         patch.object(solver, 'has_birth_dp_grant', return_value=False), \
         patch.object(solver, 'get_phi_state_matrix', return_value=np.ones((4,4))), \
         patch.object(solver, 'get_phi_choice_tensor', return_value=np.ones((4,4,2))), \
         patch.object(solver, 'get_birth_entry_grant_tensor', return_value=np.zeros((4,4,2))):
        return solver.precompute_shared(P, np.array([0.,1.]))


class ChildPreferencesTests(unittest.TestCase):
    def test_default_off_all_arrays_identical_to_unmodified_precompute(self):
        P = parameters()
        del P.child_benefit_curvature
        del P.compensated_child_housing_shares
        with patch.object(solver, 'apply_child_preferences'):
            old = shared(P)
        new = shared(P)
        self.assertEqual(set(vars(old)), set(vars(new)))
        for key, value in vars(old).items():
            np.testing.assert_array_equal(value, getattr(new,key), err_msg=key)

    def test_valid_states_and_compressed_types_use_children_at_home(self):
        P = parameters(); S = shared(P)
        for n in range(4):
            for m in range(4):
                index = n + 4*m
                active = 0 < m <= n
                expected = P.psi_child * float(m)**.86 if active else 0.
                self.assertEqual(S.psi_v[n,m], expected)
                self.assertAlmostEqual(S.alpha_flat[0,index], .533 if active else .733, places=15)
                self.assertEqual(S.hb_flat[0,index], 0.)
                if not active:
                    self.assertEqual(S.escale_flat[0,index], 1.)
        np.testing.assert_array_equal(S.type_psi[S.type_map], S.psi_flat.ravel())

    def test_compensated_expenditure_and_one_child_benefit(self):
        for loading in (0., .1, .2):
            P = parameters(); P.delta_alpha_jump = loading; S = shared(P)
            alpha = P.alpha_cons-loading
            r = P.utility_reference_rent
            K = lambda a: a**a*((1-a)/r)**(1-a)
            A = K(P.alpha_cons)/K(alpha)
            for m in (1,2,3):
                e = ((2+.7*m)/2)**.7
                self.assertAlmostEqual(S.escale_flat[0,m+4*m], e**(P.sigma-1)*A**(1-P.sigma), places=14)
            self.assertEqual(S.psi_v[1,1], P.psi_child)
            self.assertAlmostEqual(A*K(alpha), K(P.alpha_cons), places=14)

    def test_curvature_alone_retains_floor_and_negative_trial_intercept(self):
        P = parameters(); P.compensated_child_housing_shares=False
        P.child_room_floor=True; P.hbar_first_child_jump=1.89; P.psi_child=-.05
        S=shared(P)
        self.assertEqual(S.hb_flat[0,5],1.89)
        self.assertEqual(S.psi_v[3,3],P.psi_child*3.**.86)

    def test_invalid_curvature_and_mixed_adapter_are_rejected(self):
        for value in (-.01,1.,float('nan')):
            P=parameters(); P.child_benefit_curvature=value
            with self.assertRaises(ValueError): shared(P)
        P=parameters(); P.utility_comparison_arm='shares_concave'
        with self.assertRaisesRegex(ValueError,'legacy comparison'): shared(P)

    def test_invalid_compensation_contract_is_rejected(self):
        for key,value in [('utility_reference_rent',0.),('utility_reference_rent',float('nan')),
                          ('delta_alpha',.1),('delta_alpha_jump',.8),('child_room_floor',True),
                          ('hbar_first_child_jump',1.),('eqscale_form','linear')]:
            P=parameters(); setattr(P,key,value)
            with self.assertRaises(ValueError, msg=key): shared(P)


if __name__ == '__main__':
    unittest.main()
