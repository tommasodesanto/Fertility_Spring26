"""Small non-model tests of information, inherited state and full-root routing."""
from dataclasses import dataclass
from types import SimpleNamespace as NS
import copy
import unittest
from unittest.mock import patch
import numpy as np
import e5f_successive_surprises as s

@dataclass
class Conditioning:
    start_year:int=2007
    initial_mass:float=100.
    next_age_targets:object=None
    outside_flow:float=.2
    retention:float=.8
    observer:object=None

class SequencingTest(unittest.TestCase):
    def initial(self):
        return s.InheritedState(2007,NS(g_pre=np.array([1.,2.]),scheduled_entries=[1,2,3,4],scheduled_raw_entries=[5,6,7,8]))

    def test_reindex_without_reanchoring_mass(self):
        c=s.local_conditioning(Conditioning(),2015)
        self.assertEqual(c.next_age_targets,{1:2019,2:2023})
        self.assertEqual((c.initial_mass,c.outside_flow,c.retention),(100.,.2,.8))
        with self.assertRaises(ValueError):s.local_conditioning(c,2012)

    def test_sequence_carries_all_states_and_never_passes_future_shocks(self):
        seen=[];stored=[];original=self.initial()
        def solve(*,inherited,psi):
            seen.append((copy.deepcopy(inherited),psi))
            state=inherited.households;state.g_pre+=psi
            state.scheduled_entries=state.scheduled_entries[1:]+[psi]
            state.scheduled_raw_entries=state.scheduled_raw_entries[1:]+[2*psi]
            return s.SurpriseResult(None,dict(finite_horizon_market_fiscal_converged=True,terminal_distance_passed=True),
                s.InheritedState(inherited.year+4,state),dict(calendar_year=inherited.year,renter_price=psi,expected_next_asset_price=9.))
        result=s.run_sequence(initial_state=original,dated_preferences=[(2007,.2),(2011,.1)],solve_episode=solve,persist=lambda y,r:stored.append(y))
        np.testing.assert_allclose(seen[1][0].households.g_pre,[1.2,2.2])
        self.assertEqual(seen[1][0].households.scheduled_entries,[2,3,4,.2])
        self.assertEqual(seen[1][0].households.scheduled_raw_entries,[6,7,8,.4])
        np.testing.assert_array_equal(original.households.g_pre,[1,2])
        self.assertEqual(stored,[2007,2011]);self.assertEqual(result['realized_rows'][0]['renter_price'],.2)

    def test_later_surprise_cannot_change_first_realized_observation(self):
        def run(last):
            def solve(*,inherited,psi):
                return s.SurpriseResult(None,dict(finite_horizon_market_fiscal_converged=True,terminal_distance_passed=True),
                    s.InheritedState(inherited.year+4,inherited.households),dict(calendar_year=inherited.year,births=psi))
            return s.run_sequence(initial_state=self.initial(),dated_preferences=[(2007,.2),(2011,last)],solve_episode=solve,persist=lambda *a:None)
        self.assertEqual(run(.1)['realized_rows'][0],run(.9)['realized_rows'][0])

    def test_failure_persisted_and_stops_chain(self):
        calls=[];saved=[]
        def fail(**kw):
            calls.append(kw)
            return s.SurpriseResult(None,dict(finite_horizon_market_fiscal_converged=True,terminal_distance_passed=False),None,None)
        with self.assertRaises(RuntimeError):
            s.run_sequence(initial_state=self.initial(),dated_preferences=[(2007,.2),(2011,.1)],solve_episode=fail,persist=lambda *a:saved.append(a))
        self.assertEqual(len(calls),1);self.assertEqual(len(saved),1)

    def test_solver_exception_is_persisted_before_propagation(self):
        saved=[]
        def fail(**kw):raise RuntimeError('occupied infeasibility')
        with self.assertRaisesRegex(RuntimeError,'occupied infeasibility'):
            s.run_sequence(initial_state=self.initial(),dated_preferences=[(2007,.2)],solve_episode=fail,persist=lambda *a:saved.append(a))
        self.assertEqual(saved[0][1].root_receipt['error_type'],'RuntimeError')
        self.assertIsNone(saved[0][1].next_state)

    def test_first_step_uses_forecast_next_price_value_and_preserves_queues(self):
        initial=self.initial();next_house=NS(g_pre=np.array([3.,4.]),scheduled_entries=[2,3,4,.7],scheduled_raw_entries=[6,7,8,.8])
        row={k:1. for k in ('asset_price','renter_price','housing_demand','owner_rate','birth_children_topcode_adjusted','pension_period_units','payroll_tax_revenue','pension_outlays')}
        path=NS(values=[np.array([5.]),np.array([6.])],rows=[row])
        seen=[]
        def replay(**kw):
            seen.append(kw)
            return NS(values=[path.values[0]],rows=[row],terminal_state=next_house,
                maximum_mass_accounting_error=0.,maximum_policy_reproduction_error=0.,maximum_feasibility_projection_mass=0.)
        runtime=(None,NS(pf=NS(evaluate_path_at_prices=replay),person_pf=NS()),None,None,None)
        old=NS(parameters=NS(),b_grid=np.arange(2),supply_rule=object(),historical_conditioning=Conditioning())
        with patch.object(s,'_runtime',return_value=runtime):
            result=s.first_period_state(inherited=initial,old_state=old,demographics=None,path=path,prices=[2.,3.,99.],pensions=[4.,5.,6.],psi=.2)
        self.assertEqual(seen[0]['terminal_price'],3.);self.assertIs(seen[0]['terminal_V'],path.values[1])
        self.assertEqual(seen[0]['psi_path'],[.2]);self.assertIs(result.households,next_house)
        self.assertEqual(result.year,2011)

class RootRoutingTest(unittest.TestCase):
    def setUp(self):
        import test_e5f_balanced_history as fixture
        import e5f_balanced_history as balanced
        self.f=fixture.BalancedHistoryTest();self.f.setUp();f=self.f
        self.runtime=(balanced,f.joined,f.primitive,NS(terminal_convergence_diagnostics=lambda *a,**kw:dict(checks={'tail':False},all_checks_pass=False)),f.rent_domain)
        self.kw=dict(f.kw);self.kw['root_controls']={k:self.kw.pop(k) for k in ('price_bounds','pension_bounds','market_tolerance','fiscal_tolerance','market_slope','fiscal_slope','max_log_step','damping','max_evaluations','max_condition_number','worsening_factor','final_reproduction_tolerance')}
        self.kw.update(inherited=s.InheritedState(2011,NS(g_pre=f.g.copy())),psi=.2,pension_tail_tolerance=.01)

    def fake_forecast(self,**kw):
        result=self.f.fake_path(years=kw['inherited'].year+4*np.arange(len(kw['prices'])),prices=kw['prices'],
            psi_path=np.full(len(kw['prices']),kw['psi']),pension_path=kw['pensions'],
            payroll_tax_path=np.full(len(kw['prices']),.179),transfer_path=np.zeros(len(kw['prices'])),
            base_parameters=kw['old_state'].parameters,b_grid=kw['old_state'].b_grid,observer=kw['observer'])
        result.person_tail.rows=result.rows
        return result

    def solve(self,**changes):
        with patch.object(s,'_runtime',return_value=self.runtime),patch.object(s,'evaluate_forecast',side_effect=self.fake_forecast):
            return s.solve_surprise(**dict(self.kw,**changes))

    def test_full_price_pension_root_constant_preferences_and_failed_tail(self):
        result=self.solve(initial_prices=np.full(6,1.9),initial_pensions=np.full(6,1.9))
        self.assertTrue(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertGreater(len(self.f.calls),2)
        self.assertIsNone(result.next_state)
        for call in self.f.calls:np.testing.assert_array_equal(call['psi_path'],np.full(6,.2))
        self.assertEqual(self.f.calls[-1]['years'][0],2011)
        self.assertEqual(result.root_receipt['final']['payload']['trial'],len(self.f.calls))

    def test_terminal_cannot_use_future_eventual_preference(self):
        with self.assertRaisesRegex(ValueError,'own constant-preference terminal'):self.solve(psi=.1)

    def test_occupied_feasibility_audit_not_bypassed(self):
        self.f.policy_fault='nonfinite'
        with self.assertRaisesRegex(RuntimeError,'household audit'):self.solve()

    def test_normalization_unchanged(self):
        old=copy.deepcopy(self.f.old);old.diagnostics['normalization']['target']=2.
        with self.assertRaisesRegex(ValueError,'normalization'):self.solve(old_state=old)


    def test_accepted_forecast_advances_only_after_final_replay(self):
        self.runtime[3].terminal_convergence_diagnostics=lambda *a,**kw:dict(checks={'tail':True},all_checks_pass=True)
        next_state=s.InheritedState(2015,NS(g_pre=self.f.g.copy()))
        with patch.object(s,'first_period_state',return_value=next_state) as replay:
            result=self.solve()
        self.assertTrue(result.root_receipt['terminal_distance_passed'])
        self.assertIs(result.next_state,next_state)
        self.assertEqual(replay.call_count,1)
        self.assertEqual(result.realized_row['forecast_vintage_year'],2011)
        self.assertEqual(result.realized_row['expected_next_asset_price'],2.)
        self.assertFalse(result.root_receipt['horizon_verified'])

class ForecastJoinTest(unittest.TestCase):
    def test_all_historical_restart_dates_preserve_continuation_and_single_2023(self):
        for year in (2007,2011,2015,2019,2023):
            with self.subTest(year=year):
                h=(2023-year)//4;n=8;g=np.ones((1,1,1,2,1,1,1));seen=[]
                people=NS(year=2023,persons=np.array([4.]),heads=np.array([2.]));people.validated=lambda:people
                def value(i):return np.full_like(g,i)
                def backward(**kw):
                    seen.append(('backward',kw));self.assertEqual(len(kw['prices']),n-h)
                    return [value(h+i) for i in range(n-h+1)],n-h
                def historical(**kw):
                    seen.append(('history',kw));c=kw['historical_conditioning']
                    self.assertEqual(c.start_year,year);self.assertEqual(c.initial_mass,100.)
                    np.testing.assert_array_equal(kw['terminal_V'],value(h))
                    return NS(rows=[dict(calendar_year=year+4*i,period=i) for i in range(h)],
                        values=[value(i) for i in range(h+1)],terminal_state=NS(g_pre=g),bellman_solves=2*h)
                def tail(**kw):
                    seen.append(('tail',kw));self.assertEqual(kw['initial_state'].persons.year,2023)
                    self.assertEqual(len(kw['precomputed_value_path']),n-h+1)
                    return NS(rows=[dict(calendar_year=2023+4*i,period=i) for i in range(n-h)],
                        values=kw['precomputed_value_path'],bellman_solves=n-h)
                def gates(result,expected_years):
                    self.assertEqual(result.bellman_solves,2*n)
                    self.assertEqual([r['calendar_year'] for r in result.rows],expected_years)
                joined=NS(pf=NS(backward_value_path=backward,rents_from_asset_prices=lambda p,*a:p*.1,evaluate_path_at_prices=historical),
                    person_pf=NS(aggregate_heads_to_model_age_cells=lambda *a,**kw:np.array([1.,1.]),
                        evaluate_path_at_prices_person_demography=tail,PersonPFState=lambda **kw:NS(**kw)),
                    ConditionalHistoryEvaluation=lambda **kw:NS(**kw),check_smoke_gates=gates)
                old=NS(parameters=NS(age_start=18,da=4,J=2),b_grid=np.array([0.]),supply_rule=None,historical_conditioning=Conditioning())
                inherited=s.InheritedState(year,NS(g_pre=g))
                terminal=NS(policy=NS(price=[2.],V=value(n)))
                with patch.object(s,'_runtime',return_value=(None,joined,None,None,None)):
                    result=s.evaluate_forecast(inherited=inherited,old_state=old,demographics=NS(initial_person_state=people),
                        prices=np.arange(n)+2.,pensions=np.ones(n),psi=.18,terminal=terminal)
                self.assertEqual(len(result.values),n+1)
                for i,v in enumerate(result.values):np.testing.assert_array_equal(v,value(i))
                for label,kw in seen:np.testing.assert_array_equal(kw['psi_path'],np.full(len(kw['prices']),.18))
                self.assertEqual(sum(r['calendar_year']==2023 for r in result.rows),1)

class PersistenceTest(unittest.TestCase):
    def test_checkpoint_reloads_both_birth_queues_and_refuses_overwrite(self):
        import tempfile,json
        from pathlib import Path
        state=s.InheritedState(2011,NS(g_pre=np.array([1.,2.]),scheduled_entries=[1.,2.,3.,4.],scheduled_raw_entries=[5.,6.,7.,8.]))
        result=s.SurpriseResult(NS(rows=[dict(calendar_year=2007)]),dict(converged=True),state,dict(calendar_year=2007))
        with tempfile.TemporaryDirectory() as d:
            out=s.persist_episode(d,2007,result,provenance={'fixture':'test'})
            self.assertTrue(json.loads((out/'checkpoint_verification.json').read_text())['reloaded_exactly'])
            with self.assertRaises(FileExistsError):s.persist_episode(d,2007,result,provenance={'fixture':'test'})

    def test_failed_nonfinite_trial_is_saved_without_next_state(self):
        import tempfile,json
        result=s.SurpriseResult(None,dict(converged=False,score=float('inf')),None,None)
        with tempfile.TemporaryDirectory() as d:
            out=s.persist_episode(d,2007,result,provenance={'fixture':'test'})
            self.assertEqual(json.loads((out/'root_receipt.json').read_text())['score'],'inf')
            self.assertFalse((out/'next_state.pkl.gz').exists())

if __name__=='__main__':unittest.main()
