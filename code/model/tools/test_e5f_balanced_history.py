"""Pure complete-path adapter tests with fake household/demographic kernels."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import time
from types import SimpleNamespace as NS
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import e5f_balanced_history as history
from e5f_social_security import bind_social_security_income, fiscal_accounts
from e5f_parenthood_utility import PARENTHOOD_SEARCH_NAMES
from dataclasses import replace
import test_e5f_balanced_terminal as terminal_tests

NoCopyTensor = terminal_tests.NoCopyTensor


class BalancedHistoryTest(unittest.TestCase):
    def setUp(self):
        fixture = terminal_tests.BalancedTerminalTest()
        fixture.setUp()
        self.P, self.grid, self.supply = fixture.P, fixture.grid, fixture.supply
        for name in PARENTHOOD_SEARCH_NAMES:
            key = 'beta' if name == 'beta_annual' else 'hbar_first_child_jump' if name == 'h_P' else name
            if not hasattr(self.P, key): setattr(self.P, key, .5)
        self.P.age_start, self.P.da = 18, 4
        self.P.q, self.P.delta, self.P.R_gross = .04, .02, 1.04
        self.P.psi_child = .3
        bind_social_security_income(self.P, pension_period=2., payroll_tax=.179)
        self.g = np.zeros((2, 1, 1, 2, 1, 1, 1))
        self.g[0, 0, 0, 0, 0, 0, 0] = 1.
        self.g[0, 0, 0, 1, 0, 0, 0] = .716
        people = NS(year=2023, persons=np.array([3.]), heads=np.array([1.716]))
        self.demography = NS(start_year=2023, last_empirical_year=2100,
            scale_model_units_per_person=.1, initial_person_state=people,
            headship_rates=np.array([.5]), net_migration={2100: np.array([.01])})
        self.old = NS(parameters=self.P, b_grid=self.grid, supply_rule=self.supply,
            years=np.arange(2007, 2024, 4), psi_path=np.linspace(.3, .2, 5),
            initial_state=NS(g_pre=self.g.copy()),
            historical_conditioning=NS(initial_mass=self.g.sum()),
            diagnostics=dict(schema='e5f_approved_parenthood_initial_state_v1', arm='sequential',
                normalization=dict(target=2.1, completed_fertility=2.1, psi_child=.3),
                verified_solution_fertility=2.1, birth_to_entry_conversion=1/2.1,
                stationary_pension=dict(marginal_gate=True, fiscal_gate=True), preference_change_2023=-.1))
        Q = copy.deepcopy(self.P); Q.psi_child = .2
        self.terminal = NS(parameters=Q, b_grid=self.grid, policy=NS(price=np.array([2.]), V=NoCopyTensor()),
            fixed_point=NS(g_pre=self.g.copy(), persons=people), mapping_valid=True,
            endpoint=NS(contract={'supply': {n:getattr(self.supply,n) for n in
                ('mode','initial_price','initial_stock','elasticity')}}))
        self.terminal_receipt = dict(schema='e5f_balanced_terminal_v1', converged=True,
            endpoint_production_eligible=True, fresh_endpoint_matches_final=True, gates={'all':True},
            final=dict(prices=[2.], fiscal_values=[2.]))
        self.kw = dict(old_state=self.old, terminal=self.terminal,
            terminal_root_receipt=self.terminal_receipt, demographic_primitives=self.demography,
            terminal_demographic_primitives=copy.deepcopy(self.demography), count=6,
            initial_prices=np.full(6,2.), initial_pensions=np.full(6,2.),
            price_bounds=(1.,4.), pension_bounds=(.1,4.), market_tolerance=2e-4,
            audit_controls=fixture.audit,
            fiscal_tolerance=1e-6, market_slope=1., fiscal_slope=1., max_log_step=.2,
            damping=1., max_evaluations=8, deadline_monotonic=time.monotonic()+30,
            max_condition_number=1e10, worsening_factor=1.5,
            final_reproduction_tolerance=2e-10, callback=None)
        self.calls, self.paths, self.budget_calls = [], [], []
        self.joined = NS(person_pf=NS(PersonPFState=lambda g,p: NS(g_pre=g,persons=p)),
            pf=NS(rents_from_asset_prices=lambda p,*args: np.asarray(p)*.1),
            evaluate_history_and_person_tail=self.fake_path,
            check_smoke_gates=self.check_gates)
        self.primitive = NS(dated_budget=self.budget,
            policy_arrays=lambda p: dict(V=p.V, tenure_probs=p.tenure_probs))
        self.policy_fault = None
        self.checks = NS(terminal_convergence_diagnostics=lambda *args,**kwargs:
            dict(all_checks_pass=False, metrics={'resident_persons_relative_gap':.1}))
        self.rent_domain = NS(project_price_path_to_positive_rents=lambda p,**kw: (np.asarray(p),{}))

    def check_gates(self, path, expected_years):
        self.assertEqual([r['calendar_year'] for r in path.rows], expected_years)
        self.assertEqual(path.bellman_solves,2*len(expected_years))
        return {'all': {'passed': True}}

    def budget(self, evaluation, P, shared, grid, rent):
        self.budget_calls.append((P.pension,P.tau_pay,rent))
        return dict(violations=0, budget_excess_mass=0.)

    def fake_path(self, **kw):
        self.calls.append(kw)
        rows=[]
        for i,(p,b) in enumerate(zip(kw['prices'],kw['pension_path'])):
            P=copy.deepcopy(kw['base_parameters']); P.psi_child=kw['psi_path'][i]
            bind_social_security_income(P,pension_period=float(b),payroll_tax=float(kw['payroll_tax_path'][i]))
            values=np.zeros(self.g.shape);values[1]=1.
            probabilities=np.full(self.g.shape,.5)
            if self.policy_fault=='probability': probabilities[0]=1.1
            elif self.policy_fault=='value_drop': values[1]=-1.
            elif self.policy_fault=='nonfinite': values[1]=np.nan
            e=NS(g_pre=self.g.copy(),g_post_fertility=self.g.copy(),g_current=self.g.copy(),
                feasibility_projection_mass=0.,policy=NS(V=values,tenure_probs=probabilities))
            kw['observer'](i,e,P,kw['b_grid'],None)
            rows.append(dict(calendar_year=int(kw['years'][i]), housing_supply=3.,
                housing_demand=3.*(1+np.log(2./p)),**fiscal_accounts(e.g_current,P)))
        path=NS(rows=rows,bellman_solves=2*len(rows),person_tail=NS(),values=NoCopyTensor())
        self.paths.append(path)
        return path

    def solve(self, **changes):
        with mock.patch.object(history,'_runtime',return_value=(self.joined,self.primitive,self.checks,self.rent_domain)):
            return history.solve_balanced_history(**dict(self.kw,**changes))

    def test_complete_path_routes_pensions_payroll_zero_rebate_and_fresh_final(self):
        seen=[]
        result=self.solve(observer=lambda i,e,P,g,s:seen.append((i,P.pension,P.tau_pay)))
        self.assertEqual(len(self.calls),2)
        self.assertIs(result.path,self.paths[-1])
        self.assertTrue(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertFalse(result.root_receipt['horizon_verified'])
        self.assertFalse(result.root_receipt['historical_equilibrium_certified'])
        self.assertFalse(result.root_receipt['terminal_distance_passed'])
        self.assertEqual(result.root_receipt['final']['payload']['trial'],2)
        audits=result.root_receipt['final']['payload']['dated_household_audits']
        self.assertEqual(len(audits),6)
        self.assertTrue(all(all(row['gates'].values()) for row in audits))
        self.assertEqual(audits[0]['diagnostics']['occupied_negative_steps'],0)
        self.assertEqual(len(seen),12)
        for kw in self.calls:
            np.testing.assert_array_equal(kw['pension_path'],np.full(6,2.))
            np.testing.assert_array_equal(kw['payroll_tax_path'],np.full(6,.179))
            np.testing.assert_array_equal(kw['transfer_path'],np.zeros(6))
            np.testing.assert_allclose(kw['psi_path'],[.3,.275,.25,.225,.2,.2])
            self.assertIs(kw['terminal_V'],self.terminal.policy.V)
            self.assertIs(kw['initial_state'],self.old.initial_state)
            self.assertIs(kw['demographic_primitives'],self.demography)
        self.assertEqual(self.budget_calls,[(2.,.179,.2)]*12)

    def test_both_prices_and_pensions_change_and_replay_exactly(self):
        result=self.solve(initial_prices=np.full(6,1.9),initial_pensions=np.full(6,1.9))
        self.assertTrue(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertGreater(len(self.calls),2)
        np.testing.assert_array_equal(self.calls[-1]['prices'],self.calls[-2]['prices'])
        np.testing.assert_array_equal(self.calls[-1]['pension_path'],self.calls[-2]['pension_path'])

    def test_wrong_initial_normalization_or_pension_certificate_rejected(self):
        for key,value in [('target',2.),('completed_fertility',2.2),('psi_child',.4)]:
            old=copy.deepcopy(self.old);old.diagnostics['normalization'][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'normalization'):
                self.solve(old_state=old)
        old=copy.deepcopy(self.old);old.diagnostics['stationary_pension']['fiscal_gate']=False
        with self.assertRaises(ValueError):self.solve(old_state=old)

    def test_terminal_psi_receipt_coordinates_and_science_rejected(self):
        for variant in ('psi','pension','structural','receipt'):
            terminal=copy.copy(self.terminal);terminal.parameters=copy.deepcopy(self.terminal.parameters)
            receipt=copy.deepcopy(self.terminal_receipt)
            if variant=='psi':terminal.parameters.psi_child=.21
            elif variant=='pension':receipt['final']['fiscal_values']=[2.1]
            elif variant=='structural':terminal.parameters.chi+=.1
            else:receipt['converged']=False
            with self.subTest(variant=variant),self.assertRaises(ValueError):
                self.solve(terminal=terminal,terminal_root_receipt=receipt)

    def test_changed_demography_initial_mass_and_reanchored_supply_rejected(self):
        d=copy.deepcopy(self.demography);d.net_migration[2100]*=2
        with self.assertRaisesRegex(ValueError,'demographic'):self.solve(demographic_primitives=d)
        old=copy.deepcopy(self.old);old.historical_conditioning.initial_mass*=2
        with self.assertRaisesRegex(ValueError,'2007 scale'):self.solve(old_state=old)
        old=copy.copy(self.old);old.supply_rule=copy.copy(self.supply);old.supply_rule.initial_stock=4.
        with self.assertRaisesRegex(ValueError,'supply'):self.solve(old_state=old)

    def test_incomplete_dates_missing_guess_and_budget_rejected_before_path(self):
        for changes in (dict(count=5),dict(count=101),dict(initial_pensions=[2.]),
                        dict(max_evaluations=9),dict(fiscal_tolerance=1e-3)):
            with self.subTest(changes=changes),self.assertRaises(ValueError):self.solve(**changes)
        self.assertEqual(self.calls,[])

    def test_saved_fiscal_row_mismatch_or_failed_joined_gate_is_not_accepted(self):
        original=self.fake_path
        def bad(**kw):
            path=original(**kw);path.rows[0]['pension_outlays']+=.1;return path
        self.joined.evaluate_history_and_person_tail=bad
        with self.assertRaisesRegex(RuntimeError,'fiscal row'):self.solve()
        self.joined.evaluate_history_and_person_tail=original
        self.joined.check_smoke_gates=mock.Mock(side_effect=RuntimeError('failed person identity'))
        with self.assertRaisesRegex(RuntimeError,'person identity'):self.solve()

    def test_fresh_fiscal_drift_fails_even_inside_market_and_fiscal_gates(self):
        original=self.fake_path
        def drift(**kw):
            if self.calls:self.g[0,0,0,1,0,0,0]+=.0000001
            return original(**kw)
        self.joined.evaluate_history_and_person_tail=drift
        result=self.solve()
        self.assertFalse(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertFalse(result.root_receipt['gates']['fiscal_replay'])
        self.assertIs(result.path,self.paths[-1])

    def test_budget_exhaustion_returns_diagnostic_path_without_certification(self):
        result=self.solve(initial_prices=np.full(6,1.1),initial_pensions=np.full(6,.2),max_evaluations=2)
        self.assertEqual(len(self.calls),2)
        self.assertFalse(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertFalse(result.root_receipt['historical_equilibrium_certified'])

    def test_dated_probability_value_and_finite_policy_failures_rejected(self):
        for fault,gate in [('probability','probability_bounds'),('value_drop','occupied_value_monotonicity'),
                           ('nonfinite','finite_policy_arrays')]:
            self.policy_fault=fault
            with self.subTest(fault=fault),self.assertRaisesRegex(RuntimeError,gate):
                self.solve()

    def test_explicit_audit_controls_cannot_relax_retained_checks(self):
        with self.assertRaisesRegex(ValueError,'TerminalAuditControls'):
            self.solve(audit_controls=None)
        for name,value in [('probability_tolerance',1e-8),('value_drop_tolerance',1e-4),
                           ('occupied_mass_tolerance',1e-8),('feasibility_projection_tolerance',1e-3),
                           ('reconstruction_tolerance',1e-6)]:
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,name):
                self.solve(audit_controls=replace(self.kw['audit_controls'],**{name:value}))
        self.assertEqual(self.calls,[])

    def test_real_joined_composition_slices_pension_path_for_both_passes(self):
        import run_e5f_matched_pf_history as joined
        people_heads=np.zeros((2,101));people_heads[0,18]=1.;people_heads[0,22]=.716
        people=joined.person_pf.CohortState(2023,people_heads*2,people_heads)
        primitives=NS(initial_person_state=people)
        # The actual joined composition calls dataclasses.replace on conditioning.
        from dataclasses import make_dataclass
        Condition=make_dataclass('Condition',[('start_year',int),('observer',object),('validate',object)])
        condition=Condition(2007,None,mock.Mock())
        values=[np.zeros(self.g.shape) for _ in range(3)]
        prefix=NS(terminal_state=NS(g_pre=self.g),values=[np.zeros(self.g.shape)]*5,
            bellman_solves=8,rows=[dict(period=i,calendar_year=2007+4*i) for i in range(4)])
        tail=NS(values=values,bellman_solves=2,
            rows=[dict(period=i,calendar_year=2023+4*i) for i in range(2)])
        benefits=np.arange(6)*.1+2.;taxes=np.full(6,.179)
        with mock.patch.object(joined.pf,'backward_value_path',return_value=(values,2)) as backward, \
             mock.patch.object(joined.pf,'evaluate_path_at_prices',return_value=prefix) as historical, \
             mock.patch.object(joined.pf,'rents_from_asset_prices',return_value=np.ones(2)), \
             mock.patch.object(joined.person_pf,'evaluate_path_at_prices_person_demography',return_value=tail) as forward:
            result=joined.evaluate_history_and_person_tail(years=np.arange(2007,2031,4),
                prices=np.ones(6),psi_path=np.ones(6),transfer_path=np.zeros(6),terminal_price=1.,
                terminal_V=values[-1],base_parameters=self.P,b_grid=self.grid,
                initial_state=self.old.initial_state,historical_conditioning=condition,
                initial_2023_persons=people,demographic_primitives=primitives,supply_rule=self.supply,
                birth_to_entry_conversion=1/2.1,pension_path=benefits,payroll_tax_path=taxes)
        self.assertEqual(result.bellman_solves,12)
        for call in (backward,forward):
            np.testing.assert_array_equal(call.call_args.kwargs['pension_path'],benefits[4:])
            np.testing.assert_array_equal(call.call_args.kwargs['payroll_tax_path'],taxes[4:])
        np.testing.assert_array_equal(historical.call_args.kwargs['pension_path'],benefits[:4])
        np.testing.assert_array_equal(historical.call_args.kwargs['payroll_tax_path'],taxes[:4])
        self.assertIs(historical.call_args.kwargs['terminal_V'],values[0])


if __name__=='__main__':
    unittest.main()
