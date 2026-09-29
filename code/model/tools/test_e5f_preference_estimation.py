"""Synthetic fitting and surprise timing tests; no historical shock estimation."""
import copy
import tempfile
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch
import numpy as np
import e5f_preference_shock_fit as fit
import run_e5f_preference_estimation as driver


class ShockEstimationTests(unittest.TestCase):
    def targets(self,values):
        return [dict(moment='tfr_'+str(i),decision_year=2007+4*i,period=i,target=float(v)) for i,v in enumerate(values)]

    def controls(self):
        return dict(driver.draft_plan()['fit'],max_evaluations=35,total_seconds=20,
                    fertility_tolerance=1e-9,damping=1.,max_log_step=.3,reproduction_tolerance=1e-12)

    def test_scalar_fit_recovers_shock_and_requires_fresh_replay(self):
        calls=[]
        def evaluate(psi):
            self.assertIsInstance(psi,float);calls.append(psi)
            return dict(certified=True,model=2.+np.log(psi/.15))
        result=fit.fit_one(evaluate=evaluate,target=1.65,initial_level=.15,bounds=[.02,.3],controls=self.controls())
        self.assertTrue(result['converged']);self.assertAlmostEqual(result['parameter']['estimate'],.15*np.exp(-.35),places=9)
        self.assertEqual(result['equilibrium_evaluations'],len(calls));self.assertEqual(calls[-1],calls[-2])
        self.assertAlmostEqual(result['initial_fertility_derivative_log_psi'],1.)

    def test_four_surprises_fit_sequentially_with_inherited_state(self):
        desired=[.13,.12,.11,.10];state=0.;advances=[];calls=[]
        targets=self.targets([2.+np.log(psi/.15)+.1*i for i,psi in enumerate(desired)])
        def factory(stage):
            self.assertEqual(len(advances),stage)
            inherited=state
            def evaluate(psi):
                self.assertIsInstance(psi,float);self.assertEqual(state,inherited)
                calls.append((stage,psi,inherited))
                return dict(certified=True,model=2.+np.log(psi/.15)+inherited)
            return evaluate
        def advance(stage,result):
            nonlocal state
            self.assertTrue(result['converged']);advances.append(stage);state+=.1
        result=fit.fit_sequence(evaluate_factory=factory,advance=advance,targets=targets,
            initial_level=.15,bounds=[.02,.3],controls=self.controls())
        self.assertTrue(result['converged']);self.assertEqual(advances,[0,1,2,3])
        np.testing.assert_allclose([r['parameter']['estimate'] for r in result['stages']],desired,atol=1e-10)
        self.assertEqual(result['equilibrium_evaluations'],len(calls))

    def test_failed_equilibrium_is_never_scored_or_advanced(self):
        advances=[];visited=[]
        def factory(stage):
            visited.append(stage)
            return lambda psi:dict(certified=False,model=2.)
        with self.assertRaisesRegex(RuntimeError,'no certified equilibrium'):
            fit.fit_sequence(evaluate_factory=factory,advance=lambda *args:advances.append(args),
                targets=self.targets([1.]*4),initial_level=.15,bounds=[.02,.3],controls=self.controls())
        self.assertEqual(visited,[0]);self.assertFalse(advances)

    def test_unmatched_stage_stops_before_next_surprise(self):
        visited=[];advances=[]
        def factory(stage):
            visited.append(stage)
            return lambda psi:dict(certified=True,model=2.+np.log(psi/.15))
        result=fit.fit_sequence(evaluate_factory=factory,advance=lambda *args:advances.append(args),
            targets=self.targets([9.]*4),initial_level=.15,bounds=[.02,.3],
            controls=dict(self.controls(),max_evaluations=5))
        self.assertFalse(result['converged']);self.assertEqual(visited,[0]);self.assertFalse(advances)

    def test_uninformative_target_and_insufficient_budget_stop(self):
        args=dict(evaluate=lambda psi:dict(certified=True,model=2.),target=1.65,
            initial_level=.15,bounds=[.02,.3],controls=self.controls())
        with self.assertRaisesRegex(RuntimeError,'underidentified'):fit.fit_one(**args)
        args['controls']['max_evaluations']=4
        with self.assertRaisesRegex(ValueError,'Budget'):fit.fit_one(**args)

    def test_final_replay_difference_rejects_fit(self):
        count=0
        def evaluate(psi):
            nonlocal count
            count+=1
            return dict(certified=True,model=2.+np.log(psi/.15)+(0.02 if count==5 else 0.))
        result=fit.fit_one(evaluate=evaluate,target=2.,initial_level=.15,bounds=[.02,.3],controls=self.controls())
        self.assertFalse(result['converged'])

    def test_complete_fit_table_retains_validation_and_pending_rows(self):
        targets=self.targets([1.97,1.86,1.75,1.65])
        rows=fit.fit_rows(targets,[2.,1.9,1.8,1.65],'one_permanent')
        self.assertEqual([r['weight'] for r in rows],[0.,0.,0.,1.])
        rows=fit.fit_rows(targets,[1.97,None,None,None],'four_successive')
        self.assertIsNone(rows[1]['model']);self.assertIsNone(rows[1]['loss_contribution'])

    def test_disabled_plans_need_budgets_not_user_supplied_shock_values(self):
        for kind in ('one_permanent','four_successive'):
            p=driver.draft_plan(kind);missing=driver.validate_plan(p)
            self.assertEqual(p['initial_level'],'saved_reference_level')
            self.assertNotIn('shock_values',missing)
            with self.assertRaisesRegex(ValueError,'disabled'):driver.validate_plan(p,launching=True)
        with self.assertRaises(ValueError):driver.draft_plan('four_announced')

    def test_empirical_window_builder_recomputes_annual_rates(self):
        with tempfile.TemporaryDirectory() as directory:
            blocks=Path(directory)/'blocks.csv';annual=Path(directory)/'annual.csv';b=[];a=[]
            for i,year in enumerate((2007,2011,2015,2019)):
                value=2.-.1*i
                b.append(dict(decision_year=year,birth_year_start=year+1,birth_year_end=year+4,period_tfr_arithmetic_mean=value))
                a.extend(dict(year=y,period_tfr_births_per_woman=value,status='verified_published_final',source_url='https://example.org/source')
                         for y in range(year+1,year+5))
            driver.table(blocks,b);driver.table(annual,a)
            result=driver.target_contract(blocks,annual)
            np.testing.assert_allclose([r['target'] for r in result['rows']],[2.,1.9,1.8,1.7])
            a[-1]['year']=a[-2]['year'];driver.table(annual,a)
            with self.assertRaisesRegex(ValueError,'unique annual'):driver.target_contract(blocks,annual)

    def test_native_candidate_sees_only_current_level_and_selects_correct_window(self):
        for kind in ('four_successive','one_permanent'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as directory:
                p=driver.draft_plan(kind);p.update(horizons=[6,8],target_contract={'rows':self.targets([1.]*4)})
                p['budget'].update(total_seconds=30,candidate_seconds=30)
                runner=driver.NativeEstimator(p,Path(directory),{},None,None);calls=[]
                runner.endpoint=lambda psi:(calls.append(('endpoint',psi)) or ({},{'price':1.}))
                def path(psi,terminal,endpoint,H,folder):
                    self.assertIsInstance(psi,float);calls.append(('path',psi,H))
                    row={k:1. for k in ('asset_price','renter_price','adult_population','birth_children','housing_demand','pension_period_units')}
                    return dict(reference_manifest_sha256='same',source_pins={},housing='fixed_stock',
                        shock_contract={'psi':psi},horizon=H,root_and_terminal_pass=True,
                        rows=[row]*H,fertility=[{'period_tfr_topcode_adjusted':1.+i/10} for i in range(H)]),{}
                runner.path=path;reply=runner(.12)
                self.assertTrue(reply['certified']);self.assertEqual(calls,[('endpoint',.12),('path',.12,6),('path',.12,8)])
                self.assertAlmostEqual(reply['model'],1. if kind=='four_successive' else 1.3)
                self.assertTrue((Path(directory)/'best_so_far.json').exists())

    def test_nonconverged_native_candidate_cannot_return_observations(self):
        with tempfile.TemporaryDirectory() as directory:
            p=driver.draft_plan();p['budget'].update(total_seconds=30,candidate_seconds=30)
            runner=driver.NativeEstimator(p,Path(directory),{},None,None)
            def fail(psi):raise driver.CandidateRejected('endpoint failed')
            runner.endpoint=fail;reply=runner(.1)
            self.assertFalse(reply['certified']);self.assertIsNone(reply['model'])

    def test_prefix_uses_own_vintage_boundary_and_preserves_complete_state(self):
        for kind,count in (('four_successive',1),('one_permanent',4)):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as directory:
                p=driver.draft_plan(kind);p['budget'].update(total_seconds=30,candidate_seconds=30,mapping_seconds=20)
                p['target_contract']={'rows':self.targets([1.]*4)}
                rt=NS(rt={'primitive':NS(pf=NS(birth_queue_values=np.asarray))})
                runner=driver.NativeEstimator(p,Path(directory),{},None,rt)
                inherited=NS(g_pre=np.array([2.,3.]),scheduled_entries=np.array([.1,.2]),scheduled_raw_entries=np.array([.3,.4]))
                runner.inherited=inherited
                next_state=NS(g_pre=np.array([1.,4.]),scheduled_entries=np.array([.2,.5]),scheduled_raw_entries=np.array([.4,.6]))
                prices=np.arange(6.)+10;values=[np.array([float(i)]) for i in range(7)]
                rows=[dict(calendar_year=2007+4*i,asset_price=float(prices[i])) for i in range(6)]
                runner.latest=dict(psi=.12,latest=dict(prices=prices,pensions=np.ones(6),result=NS(values=values),record={'rows':rows}),
                    summary={'payload':{'candidate':5,'models':[1.]*count}})
                result=dict(converged=True,root={'final':dict(mapping_valid=True,prices=[.12],payload={'candidate':5})},
                    parameter=dict(estimate=.12,lower=.02,upper=.3,near_bound=False))
                def mapping(*args,**kwargs):
                    self.assertIs(kwargs['initial_state'],inherited);self.assertEqual(kwargs['start_year'],2007)
                    self.assertEqual(args[3]['price'],prices[count]);np.testing.assert_array_equal(args[2]['evaluation'].policy.V,values[count])
                    np.testing.assert_array_equal(args[4],prices[:count]);np.testing.assert_array_equal(args[6],[.12]*count)
                    return NS(values=values[:count+1],terminal_state=next_state),dict(gates={'all':True},rows=rows[:count],
                        fertility=[{'period_tfr_topcode_adjusted':1.}]*count,market_residual=[0.]*count,fiscal_residual=[0.]*count)
                with patch.object(driver.inner,'mapping',side_effect=mapping):runner.advance(0,result,diagnostics=False)
                self.assertEqual(runner.year,2007+4*count);self.assertIsNot(runner.inherited,next_state)
                for name in ('g_pre','scheduled_entries','scheduled_raw_entries'):
                    np.testing.assert_array_equal(getattr(runner.inherited,name),getattr(next_state,name))

    def test_pension_binding_does_not_permit_changed_earnings_or_credit(self):
        from e5f_social_security import bind_social_security_income
        P=NS(I=1,J=3,J_R=2,w_hat=np.array([1.]),income_age_profile=np.array([1.,1.2,0.]),
             period_years=4,scale_flows_to_period=True,pension=.5,tau_pay=.1,psi_child=.15,native_due_stayer_credit=True)
        bind_social_security_income(P)
        Q=copy.deepcopy(P);Q.psi_child=.12;bind_social_security_income(Q,pension_period=.6)
        driver.inner.check_endpoint_primitives(P,Q)
        Q.native_due_stayer_credit=False
        with self.assertRaises(ValueError):driver.inner.check_endpoint_primitives(P,Q)
        Q=copy.deepcopy(P);Q.income[0,0]+=1
        with self.assertRaises(ValueError):driver.inner.check_endpoint_primitives(P,Q)


if __name__=='__main__':unittest.main()
