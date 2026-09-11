"""Independent arithmetic and routing checks; no household solve."""
from types import SimpleNamespace
import unittest
import numpy as np
from e5f_stationary_paygo import (stationary_age_income_mass, bind_initial_balanced_pension,
                                certify_initial_pension, rebase_initial_supply)


def parameters():
    return SimpleNamespace(I=1, J=4, J_R=2, income_type_transition='markov',
        z_grid=np.array([.5, 1.5]), z_weights=np.array([1., 0.]),
        Pi_z=np.array([[0., 1.], [1., 0.]]), use_age_survival=True,
        survival_probs=np.array([1., .5, .5, 0.]),
        scale_flows_to_period=True, period_years=4., w_hat=np.array([2.]),
        income_age_profile=np.array([1., 2., 1., 1.]), tau_pay=.2, pension=1.,
        retirement_income_z_scale=.5, H0=np.array([3.]), xi_supply=np.array([1.75]),
        r_bar=np.array([.2]), user_cost_rate=.1)


class InitialPAYGOTests(unittest.TestCase):
    def test_nonstationary_earnings_composition_and_age_survival(self):
        expected=np.array([[[1.,0.],[0.,1.],[.5,0.],[0.,.25]]])/2.75
        np.testing.assert_allclose(stationary_age_income_mass(parameters()), expected,
                                   atol=1e-15, rtol=0)

    def test_actual_budget_and_period_income_from_hand_calculation(self):
        raw=parameters()
        P, _=bind_initial_balanced_pension(raw, payroll_tax=.2)
        # Before common mass normalization: payroll=4*(2*.5+2*2*1.5)=28;
        # retiree exposure=.5*.75+.25*1.25=.6875, revenue=5.6.
        self.assertAlmostEqual(P.pension, 5.6/.6875, places=13)
        np.testing.assert_allclose(P.income[0], [6.4,12.8,P.pension,P.pension])
        self.assertEqual(raw.pension,1.)
        g=np.array([[[[[[[1.]],[[0.]]],[[[0.]],[[1.]]],[[[.5]],[[0.]]],[[[0.]],[[.25]]]]]]])
        self.assertEqual(g.shape,(1,1,1,4,2,1,1))
        certificate=certify_initial_pension(g,P,marginal_tolerance=1e-10,fiscal_tolerance=1e-6)
        self.assertTrue(certificate['fiscal_gate'])
        # A changed age distribution must reject this stationary shortcut.
        g[0,0,0,3,1,0,0] += .05
        with self.assertRaisesRegex(RuntimeError,'certification failed'):
            certify_initial_pension(g,P,marginal_tolerance=1e-10,fiscal_tolerance=1e-6)

    def test_invalid_shortcut_and_zero_budget_fail(self):
        P=parameters(); P.I=2
        with self.assertRaises(ValueError): stationary_age_income_mass(P)
        P=parameters(); P.Pi_z[0,0]=.2
        with self.assertRaises(ValueError): stationary_age_income_mass(P)
        with self.assertRaises(ValueError): bind_initial_balanced_pension(parameters(),payroll_tax=0.)

    def test_supply_anchor_preserved_only_at_selected_price(self):
        P=parameters()
        revised, receipt=rebase_initial_supply(P,asset_prices=np.array([4.]),elasticity=.63)
        self.assertAlmostEqual(revised.H0[0],3.*2.**(1.75-.63))
        self.assertAlmostEqual(revised.H0[0]*2.**.63,3.*2.**1.75)
        self.assertEqual(P.xi_supply[0],1.75)
        self.assertEqual(receipt['new_elasticity'],[.63])


if __name__=='__main__': unittest.main()
