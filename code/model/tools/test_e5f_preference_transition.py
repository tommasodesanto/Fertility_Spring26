"""Contract and solver-wiring tests; run on Torch, without model imports."""
import copy
import tempfile
from pathlib import Path
import unittest
import numpy as np
import run_e5f_preference_transition as driver


class PreferenceTransitionTests(unittest.TestCase):
    def spec(self,kind):
        p=driver.draft_plan(kind)['shocks']
        p.update(levels=[.1] if kind=='one_permanent' else [.13,.12,.11,.1],
                 interpretation='author_supplied',provenance='synthetic fixture')
        return p

    def complete_plan(self):
        p=driver.draft_plan(); p['shocks']=self.spec('four_announced')
        p['execution_enabled']=True
        p['budget'].update(horizon=6,max_evaluations=4,total_seconds=60,case_seconds=10,
                            observed_mapping_seconds=1,maximum_policy_calls=48)
        p['numerics'].update(initial_prices=[1]*6,initial_pensions=[.2]*6,
            price_bounds=[.1,2],pension_bounds=[.01,1],market_slope=1,fiscal_slope=1,
            max_log_step=.1,damping=.5,max_condition_number=1e8,worsening_factor=2,
            terminal_tolerances={k:.01 for k in driver.TERMINAL_KEYS},raw_queue_relative_tolerance=.01)
        p.update(endpoint={'path':'synthetic'},readiness_receipt={'path':'synthetic'},source_pins={'synthetic':'test'})
        return p

    def test_one_and_four_are_distinct_announced_paths(self):
        self.assertEqual(driver.shock_path(self.spec('one_permanent'),6),[.1]*6)
        self.assertEqual(driver.shock_path(self.spec('four_announced'),6),[.13,.12,.11,.1,.1,.1])

    def test_rejects_surprises_wrong_dates_implicit_or_invalid_levels(self):
        for key,value in [('expectations','successive_surprises'),('years',[2007,2012,2015,2019]),
                          ('levels',None),('levels',[.1,.2]),('levels',[.1,.1,-.1,.1]),
                          ('interpretation',None)]:
            p=self.spec('four_announced'); p[key]=value
            with self.assertRaises(ValueError):driver.shock_path(p,6)
        with self.assertRaises(ValueError):driver.shock_path(self.spec('four_announced'),3)

    def test_drafts_are_inspectable_but_cannot_execute(self):
        for kind in ('one_permanent','four_announced'):
            p=driver.draft_plan(kind)
            self.assertIn('endpoint',driver.validate_plan(p))
            with self.assertRaisesRegex(ValueError,'disabled'):driver.validate_plan(p,launching=True)

    def test_no_credit_or_fiscal_or_entry_change_is_allowed(self):
        for key,value in [('credit','natural_solvency'),('property_rebate',.01),
                          ('outside_entry',.1),('retention',.8),('birth_to_entry_conversion',.5),
                          ('historical_shocks_automatically_imported',True)]:
            p=driver.draft_plan(); p[key]=value
            with self.assertRaises(ValueError):driver.validate_plan(p)

    def test_requires_full_budget_replay_count_and_strict_gates(self):
        p=self.complete_plan(); self.assertEqual(driver.validate_plan(p,launching=True),[])
        for section,key,value in [('budget','maximum_policy_calls',36),('budget','max_evaluations',1),
            ('budget','cache_max_bytes',3*1024**3),('numerics','fiscal_tolerance',1e-4),
            ('numerics','initial_prices',[1]),('numerics','terminal_tolerances',{}),
            ('numerics','raw_queue_relative_tolerance',float('inf'))]:
            q=copy.deepcopy(p); q[section][key]=value
            with self.assertRaises(ValueError):driver.validate_plan(q,launching=True)

    def test_fixed_stock_uses_reference_quantity_not_supply_intercept(self):
        from types import SimpleNamespace
        class Supply:
            def quantity(self,price):return np.array([5.8])
        class Fixed:
            def __init__(self,*args):self.args=args
        packet=dict(supply_rule=Supply(),solution=SimpleNamespace(p_eq=np.array([.8])))
        pf=SimpleNamespace(calendar=SimpleNamespace(HousingSupplyRule=Fixed))
        fixed=driver.supply_rule(packet,pf,'fixed_stock')
        self.assertEqual(fixed.args,('fixed-stock',.8,5.8,0.))
        self.assertIs(driver.supply_rule(packet,pf,'elastic_reference'),packet['supply_rule'])

    def test_horizon_certificate_checks_same_experiment_and_pensions(self):
        row={k:1. for k in ('asset_price','renter_price','adult_population','birth_children',
                            'housing_demand','pension_period_units')}
        a=dict(horizon=4,rows=[row]*4,root_and_terminal_pass=True,reference_manifest_sha256='same',
               source_pins={'same':'source'},housing='fixed_stock',shock_contract=self.spec('four_announced'))
        b=copy.deepcopy(a); b.update(horizon=8,rows=[copy.deepcopy(row) for _ in range(8)])
        self.assertTrue(driver.compare_horizons(a,b,3,1e-3)['passed'])
        b['rows'][0]['pension_period_units']=1.1
        self.assertFalse(driver.compare_horizons(a,b,3,1e-3)['passed'])
        b['shock_contract']['levels'][0]=.14
        with self.assertRaises(ValueError):driver.compare_horizons(a,b,3,1e-3)

    def test_raw_entry_queue_failure_overrides_terminal_pass(self):
        from types import SimpleNamespace as NS
        stationary=NS(scheduled_raw_entries=np.array([.2,.2]))
        pf=NS(stationary_initial_state=lambda *args:stationary,
              terminal_convergence_diagnostics=lambda **kwargs:dict(status='passed',all_checks_pass=True),
              birth_queue_values=lambda queue:queue)
        evaluator=NS(rt={'primitive':NS(pf=pf)})
        terminal=dict(parameters=NS(psi_child=.1),stationary_g_pre=np.ones((1,1,1,1)),
                      evaluation=NS(births=.42))
        result=NS(terminal_state=NS(scheduled_raw_entries=np.array([.2,.3])))
        check=driver.terminal_checks({},evaluator,terminal,dict(price=1.,population_scale=1.),
            result,[.1],dict(terminal_tolerances={},raw_queue_relative_tolerance=.01))
        self.assertFalse(check['all_checks_pass'])
        self.assertEqual(check['status'],'not_converged')
        self.assertFalse(check['raw_queue_pass'])

    def test_jacobian_measures_two_physical_blocks_in_five_calls(self):
        calls=[]
        def evaluate(q,p):
            calls.append((q.copy(),p.copy()))
            return dict(mapping_valid=True,market_residual=-2*np.log(q)+.5*np.log(p),
                        fiscal_residual=.2*np.log(q)-np.log(p))
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'jacobian'
            matrix=driver.measure_jacobian(evaluate,np.ones(3),np.ones(3),1,1e-5,output,{'housing':'fixed_stock'})
            np.testing.assert_allclose(matrix[:3,:3],-2*np.eye(3),atol=1e-10)
            np.testing.assert_allclose(matrix[3:,3:],-np.eye(3),atol=1e-10)
            np.testing.assert_allclose(matrix[:3,3:],.5*np.eye(3),atol=1e-10)
            self.assertEqual(len(calls),5)
            self.assertEqual(driver.read(output/'receipt.json')['residual_units'],'physical_unscaled')


if __name__=='__main__':unittest.main()
