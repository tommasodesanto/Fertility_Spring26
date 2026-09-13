"""Pure joint-root wiring tests. No native model solve or launch."""
from dataclasses import dataclass
import importlib.util
from pathlib import Path
import sys
import tempfile
import time
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_e5f_final_rebated_history as driver


@dataclass(frozen=True)
class Demo:
    net_migration: dict


class DriverTests(unittest.TestCase):
    def test_seed_profiles_accept_legacy_lists_and_named_profiles(self):
        self.assertEqual(driver.initial_seed_step({'seed_steps':[-.02,-.03]},{}),-.02)
        self.assertEqual(driver.initial_seed_step({'seed_steps':{'a':-.04}},{}),-.04)
        self.assertEqual(driver.initial_seed_step({},dict(seed_step=-.05)),-.05)

    def test_migration_zeroes_each_cell_without_mutating_comparison(self):
        demo=Demo({2024:np.array([[2.,-2.]]),2100:np.array([[3.,4.]])})
        zero=driver.migration_case(demo,'A0')
        self.assertTrue(all(np.count_nonzero(a)==0 for a in zero.net_migration.values()))
        np.testing.assert_array_equal(demo.net_migration[2024],[[2.,-2.]])
        self.assertIs(driver.migration_case(demo,'A+'),demo)

    def test_coordinates_include_independent_boundary_in_every_block(self):
        parts=driver.unpack_coordinates(np.arange(1,22,dtype=float),6)
        self.assertEqual([p[-1] for p in parts],[7.,14.,21.])
        with self.assertRaises(ValueError):driver.unpack_coordinates(np.ones(18),6)

    def test_source_hash_drift_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source';p.write_text('first')
            pins={str(p):driver.sha(p)};driver.verify_pins(pins)
            p.write_text('changed')
            with self.assertRaises(ValueError):driver.verify_pins(pins)

    def test_joint_root_uses_actual_terminal_state_and_single_boundary_solve(self):
        count=6;g=np.ones((2,1,1,2,1,1,1));actual_g=3*g
        P=NS(psi_child=.1,pension=1.,property_tax_lump_sum_transfer=.2,tau_pay=.179)
        old=NS(parameters=P,b_grid=np.array([0.,1.]),supply_rule=NS())
        inherited=NS(year=2007,households=NS(g_pre=g))
        calls=[];actual_calls=[]
        pf=NS(rents_from_asset_prices=lambda p,t,P:np.ones(len(p)),solve_date_policy=lambda **kw:None)
        primitive=NS(model=NS(property_tax_revenue_from_distribution=lambda *args:float(g.sum())*.2))
        runtime=(None,NS(pf=pf),primitive,None,
            NS(project_price_path_to_positive_rents=lambda p,**kw:(p,None)))
        def boundary(**kw):
            calls.append(kw['g_pre'].copy())
            return NS(parameters=kw['parameters'],policy=NS(V=np.ones_like(g),price=np.array([kw['price']])))
        def cached(template,carried,*args):
            actual_calls.append(carried.copy())
            return NS(mapping_valid=True,actual_accounts={'household_heads':float(carried.sum())},
                residuals=dict(housing_relative=0.,pension_relative=0.,rebate_relative=0.),gates={'valid':True})
        def forecast(**kw):
            rows=[]
            for i in range(count):
                q=NS(pension=kw['pensions'][i],property_tax_lump_sum_transfer=kw['transfers'][i],tau_pay=.179)
                e=NS(g_current=g,policy=NS(hR_pol=np.ones_like(g),price=np.ones(1)),
                    demand_by_loc=np.ones(1),supply_by_loc=np.ones(1))
                kw['observer'](i,e,q,np.array([0.,1.]),NS())
                rows.append(dict(annual_net_migration_over_period=0.,net_migrant_heads_over_period=0.))
            return NS(rows=rows,person_tail=NS(rows=rows,terminal_state=NS(g_pre=actual_g)))
        def root(**kw):
            x=kw['project'](kw['initial_prices'])
            first=kw['evaluate'](x);second=kw['evaluate'](x)
            self.assertEqual(first['residual'].shape,(21,))
            np.testing.assert_array_equal(first['residual'],second['residual'])
            return dict(converged=True,final=dict(prices=x,payload=second['payload']))
        rebated=NS(_runtime=lambda:runtime,evaluate_forecast=forecast,
            rebated_tax_accounts=lambda **kw:{},dated_residual=lambda **kw:np.zeros(3),
            stack_dated_residuals=lambda rows:np.asarray(rows).T.ravel(),
            first_period_state=lambda **kw:NS(year=2011))
        modules={'e5f_rebated_surprises':rebated,
            'e5f_closed_finite_boundary':NS(boundary_evaluation=boundary),
            'e5f_balanced_terminal':NS(_household_checks=lambda *args:({},dict(valid=True))),
            'e5f_social_security':NS(fiscal_accounts=lambda *args:{}),
            'e5f_matched_pf_path_root':NS(solve_price_path=root),
            'run_e5f_transition_calibration':NS(period_fertility_diagnostics=lambda *args:{'period_tfr_topcode_adjusted':1.7})}
        controls=dict(price_bounds=[.01,5.],pension_bounds=[.01,5.],transfer_bounds=[.01,5.],
            slope=1.,market_tolerance=2e-4,max_log_step=.2,damping=.5,max_evaluations=8,
            max_condition_number=1e8,worsening_factor=2.,final_reproduction_tolerance=2e-10)
        with tempfile.TemporaryDirectory() as d,patch.dict(sys.modules,modules),patch.object(driver,'cached_boundary',cached),patch.object(driver,'checkpoint',lambda *args:None):
            result,detail=driver.solve_forecast(inherited=inherited,old=old,demographics=NS(),
                psi=.1,count=count,initial=np.r_[np.ones(14),np.full(7,.2)],controls=controls,
                audit=NS(),deadline=time.monotonic()+10,folder=d,case='A0')
        self.assertEqual(len(calls),2)  # One per mapping, including fresh replay.
        self.assertEqual(len(actual_calls),2)
        for value in actual_calls:np.testing.assert_array_equal(value,actual_g)
        self.assertTrue(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertFalse(result.root_receipt['horizon_verified'])
        self.assertEqual(result.next_state.year,2011)


if __name__=='__main__':unittest.main()
