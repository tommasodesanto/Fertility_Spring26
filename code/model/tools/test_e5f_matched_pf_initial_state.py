"""Pure contracts and explicit-parameter routing tests; no numerical solve."""
from types import SimpleNamespace
import time
import unittest
from unittest.mock import patch

import numpy as np
import e5f_matched_pf_initial_state as initial


class InitialStateContracts(unittest.TestCase):
    def fiscal(self):
        return initial.CalibrationFiscalContract(.01, .04, 0., 'retained_calibration_unrebated')

    def summary(self):
        return dict(old_psi_child=.3, best_candidate={'new_psi_child': -.1},
            housing_supply=dict(retained_date0_asset_price=.625,
                date0_normalized_housing_stock=6.4, transition_housing_supply_elasticity=.63,
                retained_housing_supply_elasticity=1.75))

    def test_renormalization_preserves_estimated_change_not_absolute_endpoint(self):
        years, path, delta = initial.historical_preference_path(.5, self.summary())
        np.testing.assert_array_equal(years, [2007, 2011, 2015, 2019, 2023])
        np.testing.assert_allclose(path, [.5, .4, .3, .2, .1], atol=1e-15, rtol=0)
        self.assertAlmostEqual(delta, -.4)

    def test_fiscal_rate_units_and_no_silent_rebate(self):
        P = SimpleNamespace(period_years=4, tau_H=.04, property_tax_lump_sum_transfer=0.)
        self.fiscal().validate(P)
        P.tau_H = .01
        with self.assertRaisesRegex(ValueError, 'parameters differ'):
            self.fiscal().validate(P)
        P.tau_H = .04
        P.property_tax_lump_sum_transfer = .02
        with self.assertRaisesRegex(ValueError, 'parameters differ'):
            self.fiscal().validate(P)

    def test_rebased_supply_packet_is_rejected(self):
        rule = SimpleNamespace(mode='static-elastic', initial_price=.625, initial_stock=6.4, elasticity=.63)
        initial.verify_selected_supply_anchor(rule, self.summary())
        rule.initial_price = .9
        with self.assertRaisesRegex(ValueError, 'initial_price'):
            initial.verify_selected_supply_anchor(rule, self.summary())

    def test_adapter_preserves_provided_parameters_grid_and_old_supply(self):
        P = SimpleNamespace(period_years=4, tau_H=.04, property_tax_lump_sum_transfer=0.,
            Nb=120, J=17, I=1, beta=.97, psi_child=-.1, H0=np.array([13.7]),
            kappa_fert=2., kappa_fert_continuation=1.8, xi_supply=np.array([1.75]),
            tenure_choice_kappa=.005, tol_eq=2.5e-5, normalize_transition_mass_roundoff=True)
        grid = np.linspace(0, 10, 120)
        rule = SimpleNamespace(mode='static-elastic', initial_price=.625, initial_stock=6.4, elasticity=.63)
        calls = []
        class EndOfRoutingTest(Exception):
            pass
        def ge(price, actual, actual_grid, verbose):
            calls.append((price, actual, actual_grid))
            self.assertEqual(actual.beta, .97)
            np.testing.assert_array_equal(actual.H0, P.H0)
            np.testing.assert_array_equal(actual.xi_supply, [1.75])
            np.testing.assert_array_equal(actual_grid, grid)
            self.assertFalse(actual.joint_nested_choice)
            self.assertTrue(actual.exhaustive_saving_control)
            self.assertEqual(actual.psi_child, .4)
            return SimpleNamespace(converged=True, timings={'strict_converged': True, 'best_eq_error': 0.}), actual, price
        model = SimpleNamespace(make_grid=lambda _: grid.copy(), solve_markov_income_equilibrium=ge)
        chain = SimpleNamespace(extract_moments=lambda *_: {'tfr': 2.1})
        def normalization(adapter, overrides, **kwargs):
            self.assertEqual(overrides, {})
            self.assertEqual(kwargs['initial_psi'], .3)
            self.assertTrue(kwargs['normalize'])
            adapter.run_model_cp_dt({'psi_child': .4})
            with self.assertRaisesRegex(ValueError, 'unauthorized primitive'):
                adapter.run_model_cp_dt({'psi_child': .4, 'H0': 99.})
            raise EndOfRoutingTest()
        with patch.object(initial.primitive, 'verify_selected_parameters'), \
             patch.object(initial.pf.transition, 'configure_sequential_model', return_value=(chain, model)), \
             patch.object(initial.calibration, 'solve_old_steady_state', side_effect=normalization):
            with self.assertRaises(EndOfRoutingTest):
                initial.initialize_normalized_old_state(parameters=P, b_grid=grid,
                    selected_summary=self.summary(), selected_supply_rule=rule, arm='sequential',
                    fiscal_contract=self.fiscal(), completed_fertility_tolerance=5e-4,
                    max_stationary_solves=2, deadline_monotonic=time.monotonic() + 30.)
        self.assertEqual(len(calls), 1)
        self.assertEqual(P.psi_child, -.1)
        self.assertFalse(hasattr(P, 'joint_nested_choice'))


if __name__ == '__main__':
    unittest.main()
