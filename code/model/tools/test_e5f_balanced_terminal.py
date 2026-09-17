"""Pure fake-evaluator and routing tests; no household model or Numba work."""
from __future__ import annotations

import copy
from dataclasses import replace
import inspect
from pathlib import Path
import sys
import time
from types import SimpleNamespace as NS
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import e5f_balanced_terminal as balanced
from e5f_social_security import fiscal_accounts


class NoCopyTensor:
    def __deepcopy__(self, memo):
        raise AssertionError('Root attempted to copy an endpoint tensor')


class BalancedTerminalTest(unittest.TestCase):
    def setUp(self):
        self.P = NS(I=1, J=2, J_R=1, Nb=2, period_years=4.,
            scale_flows_to_period=True, tau_H=.04, tau_pay=.179,
            pension=7., property_tax_lump_sum_transfer=0., psi_child=.2,
            user_cost_rate=.1, w_hat=np.array([2.]),
            income_age_profile=np.array([1., 0.]), z_grid=np.array([1.]),
            retirement_income_z_scale=0., joint_nested_choice=False,
            fertility_nest_choice=False, two_shock_choice=False,
            exhaustive_saving_control=True, H0=np.array([3.]),
            r_bar=np.array([.2]), xi_supply=np.array([.63]),
            child_state_mode='independent_count', n_parity=4, n_child_states=4,
            use_stochastic_aging=True, alpha_cons=.733,
            hbar_child_rooms=0., hbar_first_child_jump=.8, hR_max=8.,
            preference_spec='eqscale', eqscale_form='power', child_room_floor=True,
            sigma=2., delta_alpha=0., delta_alpha_jump=0., c_bar_0=0., c_bar_n=0.,
            beta=.9952765791**4)
        self.grid = np.array([0., 1.])
        self.demographics = NS(start_year=2023, last_empirical_year=2100)
        self.supply = NS(mode='static-elastic', initial_price=2.,
                         initial_stock=3., elasticity=.63)
        self.supply.quantity = lambda prices: np.array([self.supply.initial_stock
            * (float(prices[0]) / self.supply.initial_price) ** self.supply.elasticity])
        self.controls = balanced.EndpointControls(20, .5, 1e-9, 1e-10,
                                                  1e-8, 2e-4, 2.5e-5, 2e-9)
        self.audit = balanced.TerminalAuditControls(5e-9, 1e-6, 1e-12, 1e-12, 1e-7)
        self.args = dict(parameters=self.P, b_grid=self.grid,
            demographic_primitives=self.demographics, supply_rule=self.supply,
            controls=self.controls, audit_controls=self.audit,
            start_price=2., start_pension_period=1., price_bounds=(.5, 5.),
            pension_bounds=(.1, 5.), fiscal_tolerance=1e-6,
            market_slope=1., fiscal_slope=1., max_log_step=.5, damping=1.,
            max_evaluations=20, deadline_monotonic=time.monotonic()+30,
            max_condition_number=1e8, worsening_factor=2.,
            final_reproduction_tolerance=1e-10, callback=None)
        self.calls = []

    def fake(self, *, asset_price, pension_period, **kwargs):
        self.calls.append((asset_price, pension_period))
        P = copy.deepcopy(self.P)
        P.pension = pension_period
        market = float(np.log(2. / asset_price))
        account = dict(payroll_tax_revenue=1., pension_outlays=pension_period,
            pension_budget_residual=1.-pension_period,
            scaled_pension_budget_residual=(1.-pension_period)/max(1., pension_period))
        endpoint = NS(mapping_valid=True,
            residuals=dict(housing_relative=market), gates=dict(inner=True),
            fixed_point=NoCopyTensor())
        return balanced.BalancedTerminalEndpoint(P, self.grid,
            NS(price=np.array([asset_price]), tensor=NoCopyTensor()),
            endpoint, account, dict(check='fake'), dict(households=True))

    def solve(self, fake=None, **changes):
        with mock.patch.object(balanced, '_evaluate_terminal_trial', side_effect=fake or self.fake):
            return balanced.solve_balanced_terminal(**dict(self.args, **changes))

    def test_exact_root_still_executes_fresh_replay_without_copying_tensors(self):
        result = self.solve()
        self.assertEqual(self.calls, [(2., 1.), (2., 1.)])
        self.assertTrue(result.production_eligible)
        self.assertEqual(result.root_receipt['evaluations'], 2)
        self.assertEqual(result.root_receipt['returned_endpoint'], 'fresh_final')
        self.assertEqual(result.root_receipt['final']['payload']['trial'], 2)
        self.assertIsInstance(result.endpoint.fixed_point, NoCopyTensor)

    def test_joint_root_changes_both_unknowns_and_replays_same_parameters(self):
        result = self.solve(start_price=1.9, start_pension_period=.8)
        self.assertTrue(result.production_eligible)
        self.assertGreater(len(self.calls), 2)
        self.assertAlmostEqual(self.calls[-1][0], 2., delta=4e-4)
        self.assertAlmostEqual(self.calls[-1][1], 1., delta=1e-6)
        np.testing.assert_array_equal(result.root_receipt['best']['prices'],
                                      result.root_receipt['final']['prices'])
        np.testing.assert_array_equal(result.root_receipt['best']['fiscal_values'],
                                      result.root_receipt['final']['fiscal_values'])

    def test_replay_accounting_drift_blocks_eligibility(self):
        def drift(**kw):
            trial = self.fake(**kw)
            if len(self.calls) == 2:
                trial.social_security['pension_outlays'] += 5e-7
            return trial
        result = self.solve(fake=drift)
        self.assertFalse(result.production_eligible)
        self.assertTrue(result.root_receipt['gates']['social_security'])
        self.assertFalse(result.root_receipt['gates']['fiscal_replay'])

    def test_failed_household_mapping_is_not_market_equilibrium(self):
        def invalid(**kw):
            trial = self.fake(**kw)
            trial.gates['households'] = False
            return trial
        result = self.solve(fake=invalid)
        self.assertFalse(result.production_eligible)
        self.assertEqual(len(self.calls), 1)
        self.assertIsNone(result.endpoint)

    def test_failed_final_mapping_cannot_be_replaced_by_valid_best(self):
        def invalid_replay(**kw):
            trial = self.fake(**kw)
            if len(self.calls) == 2:
                trial.gates['households'] = False
            return trial
        result = self.solve(fake=invalid_replay)
        self.assertFalse(result.production_eligible)
        self.assertFalse(result.endpoint.mapping_valid)

    def test_budget_is_bounded_and_unconverged_endpoint_is_diagnostic(self):
        result = self.solve(start_price=1., max_evaluations=2)
        self.assertEqual(len(self.calls), 2)
        self.assertFalse(result.production_eligible)
        self.assertIsNotNone(result.endpoint)

    def test_timeout_preserves_best_without_claiming_fresh_replay(self):
        def timeout(**kw):
            if self.calls:
                raise TimeoutError('caller watchdog')
            return self.fake(**kw)
        result = self.solve(fake=timeout)
        self.assertFalse(result.production_eligible)
        self.assertFalse(result.root_receipt['fresh_endpoint_matches_final'])
        self.assertEqual(result.root_receipt['returned_endpoint'], 'best_diagnostic')
        self.assertIsNotNone(result.endpoint)

    def test_no_unannounced_defaults_or_diagnostic_population_substitutions(self):
        for name, parameter in inspect.signature(balanced.solve_balanced_terminal).parameters.items():
            if name != 'initial_jacobian':
                self.assertIs(parameter.default, inspect.Parameter.empty, name)
        cases = [dict(price_bounds=(3., 4.)), dict(pension_bounds=(2., 4.)),
            dict(fiscal_tolerance=1e-4), dict(controls=replace(self.controls, market_tolerance=.01)),
            dict(controls=replace(self.controls, maximum_inner_iterations=0)),
            dict(controls=replace(self.controls, accounting_absolute_tolerance=float('nan'))),
            dict(audit_controls=replace(self.audit, value_drop_tolerance=.1)),
            dict(demographic_primitives=NS(start_year=2007, last_empirical_year=2100))]
        for change in cases:
            with self.subTest(change=change), mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
                with self.assertRaises(ValueError):
                    balanced.solve_balanced_terminal(**dict(self.args, **change))
                trial.assert_not_called()

    def test_wrong_fiscal_or_choice_contract_rejected_before_households(self):
        for name, value in [('tau_H', .08), ('tau_pay', .18),
                            ('property_tax_lump_sum_transfer', .1),
                            ('period_years', 1.), ('joint_nested_choice', True)]:
            P = copy.deepcopy(self.P)
            setattr(P, name, value)
            with self.subTest(name=name), mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
                with self.assertRaises(ValueError):
                    balanced.solve_balanced_terminal(**dict(self.args, parameters=P))
                trial.assert_not_called()

    def test_approved_utility_and_exhaustive_saving_are_required_before_solving(self):
        for name, value in [('exhaustive_saving_control', False),
                            ('hbar_child_rooms', .1), ('sigma', 1.5),
                            ('eqscale_form', 'linear'), ('alpha_cons', .7),
                            ('c_bar_n', .01), ('child_state_mode', 'legacy')]:
            P = copy.deepcopy(self.P)
            setattr(P, name, value)
            with self.subTest(name=name), mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
                with self.assertRaises(ValueError):
                    balanced.solve_balanced_terminal(**dict(self.args, parameters=P))
                trial.assert_not_called()

    def test_supply_elasticity_and_parameter_curve_must_match_without_rebasing(self):
        for name, value in [('H0', np.array([3.1])), ('r_bar', np.array([.21])),
                            ('xi_supply', np.array([1.75]))]:
            P = copy.deepcopy(self.P)
            setattr(P, name, value)
            with self.subTest(name=name), mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
                with self.assertRaisesRegex(ValueError, 'supply'):
                    balanced.solve_balanced_terminal(**dict(self.args, parameters=P))
                trial.assert_not_called()
        for name, value in [('elasticity', 1.75), ('initial_stock', 3.1)]:
            old = getattr(self.supply, name)
            setattr(self.supply, name, value)
            with self.subTest(name=name), mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
                with self.assertRaisesRegex(ValueError, 'supply'):
                    balanced.solve_balanced_terminal(**self.args)
                trial.assert_not_called()
                self.assertEqual(getattr(self.supply, name), value)
            setattr(self.supply, name, old)

    def test_supply_guard_checks_two_anchor_prices(self):
        quantity = self.supply.quantity
        probed = []
        def wrong_slope(prices):
            probed.append(float(prices[0]))
            value = quantity(prices)
            return value if len(probed) == 1 else value * 1.001
        self.supply.quantity = wrong_slope
        with mock.patch.object(balanced, '_evaluate_terminal_trial') as trial:
            with self.assertRaisesRegex(ValueError, 'supply law differs'):
                balanced.solve_balanced_terminal(**self.args)
            trial.assert_not_called()
        self.assertEqual(probed, [2., 2.2])

    def test_zero_budget_rule(self):
        for revenue, outlays, expected in [(0., 0., 0.), (0., 1e-15, -1.)]:
            def zero(**kw):
                trial = self.fake(**kw)
                trial.social_security.update(payroll_tax_revenue=revenue, pension_outlays=outlays)
                return trial
            result = self.solve(fake=zero, max_evaluations=2)
            payload = result.root_receipt['best']['payload']
            self.assertEqual(payload['scaled_pension_budget_residual'], expected)
            self.assertEqual(result.production_eligible, expected == 0.)

    def test_actual_trial_binds_before_households_preserves_inputs_and_accounts_actual_mass(self):
        events = []
        shape = (2, 1, 1, 2, 1, 1, 1)
        seed = np.ones(shape)
        actual = seed.copy()
        actual[:, :, :, 0] *= 3.
        policy = NS(price=np.array([2.]))
        endpoint = NS(mapping_valid=True, gates=dict(inner=True),
            fixed_point=NS(g_pre=actual), residuals=dict(housing_relative=0.,
                housing_demand=3., housing_supply=3., tax_revenue=.4, renewal_ratio=.8))
        shared = object()
        def precompute(P, grid):
            events.append('shared')
            self.assertEqual(P.pension, 1.5)
            self.assertEqual(P.tau_pay, .179)
            np.testing.assert_allclose(P.income, [[4.*(1.-.179)*2., 1.5]])
            return shared
        def household(price, P, grid, SD):
            events.append('household')
            self.assertIs(SD, shared)
            self.assertEqual(P.social_security_income_units, 'period')
            return NS()
        def reconstruct(*args):
            events.append('seed')
            return seed, dict(stationary_post_fertility_nesting_l1=0.,
                stationary_post_fertility_nesting_max_abs=0., stationary_feasibility_projection_mass=0.)
        def inner(**kw):
            events.append('terminal')
            self.assertIs(kw['initial_g_pre'], seed)
            self.assertIs(kw['demographic_primitives'], self.demographics)
            self.assertIs(kw['supply_rule'], self.supply)
            self.assertEqual(kw['fiscal_regime'], 'fixed_transfer')
            self.assertEqual(kw['transfer'], 0.)
            self.assertEqual(kw['parameters'].pension, 1.5)
            return endpoint
        current = NS(g_current=actual, demand_by_loc=np.array([3.]), supply_by_loc=np.array([3.]))
        calendar = NS(policy_from_solution=mock.Mock(return_value=policy),
            reconstruct_stationary_pre_fertility=reconstruct,
            evaluate_period=mock.Mock(return_value=current), SolveCounter=lambda: None)
        model = NS(precompute_shared=precompute, solve_markov_income_at_prices=household)
        with mock.patch.object(balanced, '_runtime', return_value=(model, calendar, object())), \
             mock.patch.object(balanced, 'evaluate_endpoint', side_effect=inner), \
             mock.patch.object(balanced, '_household_checks', return_value=({}, {'audit': True})):
            result = balanced._evaluate_terminal_trial(parameters=self.P, b_grid=self.grid,
                demographic_primitives=self.demographics, supply_rule=self.supply,
                controls=self.controls, audit_controls=self.audit, asset_price=2., pension_period=1.5)
        self.assertEqual(events, ['shared', 'household', 'seed', 'terminal'])
        self.assertEqual(self.P.pension, 7.)
        self.assertFalse(hasattr(self.P, 'income'))
        expected = fiscal_accounts(actual, result.parameters)
        self.assertEqual(result.social_security, expected)
        self.assertNotEqual(result.social_security['payroll_tax_revenue'],
                            fiscal_accounts(seed, result.parameters)['payroll_tax_revenue'])
        self.assertEqual(result.social_security['pension_outlays'], 3.)
        self.assertTrue(result.mapping_valid)
        self.assertIs(calendar.evaluate_period.call_args.kwargs['supply_rule'], self.supply)

    def test_household_checks_reject_value_probability_and_finite_failures(self):
        mass = np.ones((2, 1, 1, 2, 1, 1, 1))
        policy = NS(V=np.broadcast_to(np.array([0., 1.]).reshape(2, 1, 1, 1, 1, 1, 1), mass.shape).copy(),
                    fert_probs=np.array([.5, .5]))
        ev = NS(policy=policy, g_pre=mass, g_post_fertility=mass,
                g_current=mass, feasibility_projection_mass=0.)
        primitive = NS(dated_budget=lambda *a: dict(budget_excess_mass=0.),
                       policy_arrays=lambda p: dict(V=p.V, fert_probs=p.fert_probs))
        def check():
            return balanced._household_checks(ev, self.P, None, self.grid, .2, primitive, self.audit)[1]
        self.assertTrue(all(check().values()))
        policy.V[1] = -1.
        self.assertFalse(check()['occupied_value_monotonicity'])
        policy.fert_probs[0] = 1.1
        self.assertFalse(check()['probability_bounds'])
        policy.fert_probs[0] = np.nan
        self.assertFalse(check()['finite_policy_arrays'])
        ev.feasibility_projection_mass = 2e-6
        self.assertFalse(check()['feasibility_projection'])


if __name__ == '__main__':
    unittest.main()
