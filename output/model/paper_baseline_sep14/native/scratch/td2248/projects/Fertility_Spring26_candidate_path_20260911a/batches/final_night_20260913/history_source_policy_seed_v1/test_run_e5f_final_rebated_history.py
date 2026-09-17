"""Pure joint-root wiring tests. No native model solve or launch."""
from dataclasses import dataclass
import copy
import importlib.util
import json
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
import run_e5f_forecast_jacobian_probe as jacobian_probe


@dataclass(frozen=True)
class Demo:
    net_migration: dict


class DriverTests(unittest.TestCase):
    def resume_fixture(self,d):
        d=Path(d);source=d/'source';current=d/'current';trial=d/'history/window_2007/trial_00'
        source.mkdir();current.mkdir();trial.mkdir(parents=True)
        kernel=source/'kernel.py';kernel.write_text('kernel')
        old_driver=source/'run_e5f_final_rebated_history.py';old_driver.write_text('old')
        new_driver=current/'run_e5f_final_rebated_history.py';new_driver.write_text('new')
        source_manifest=dict(prior_plan='/fixed/plan.json',initial_summary='/fixed/summary.json',
            file_sha256={str(kernel):driver.sha(kernel),str(old_driver):driver.sha(old_driver)},
            policy_reserve_seconds=100)
        source_manifest_path=d/'source_manifest.json';driver.save(source_manifest_path,source_manifest)
        source_root=d/'model_source';initial_digest='initial-checkpoint'
        contract=dict(manifest_sha256=driver.sha(source_manifest_path),case='A0',count=6,
            initial_checkpoint_sha256=initial_digest,source_root=str(source_root))
        contract_path=d/'history/contract_receipt.json';driver.save(contract_path,contract)
        fit=dict(year=2007,psi=.12,target=1.7,model=1.7,gap=0.,folder=str(trial))
        fit_path=trial/'fit.json';driver.save(fit_path,fit)
        prices=np.arange(1.,22.)
        root=dict(converged=True,status='converged',finite_horizon_market_fiscal_converged=True,
            start_year=2007,case='A0',count=6,psi=.12,final_reproduction_max_abs=0.,
            final=dict(prices=prices.tolist(),mapping_valid=True))
        root_path=trial/'root_receipt.json';driver.save(root_path,root)
        state=NS(year=2011,households=NS(g_pre=np.arange(4.).reshape(2,2),
            persons=NS(year=2011,persons=np.array([3.,4.]),heads=np.array([1.,2.]))))
        accepted_path=trial/'accepted_forecast.pkl.gz'
        driver.checkpoint(accepted_path,dict(result=NS(next_state=state,root_receipt=root),
            boundary=NS(),coordinates=prices.copy()))
        state_path=d/'history/realized_state_2011.pkl.gz';driver.checkpoint(state_path,state)
        realized_path=d/'history/realized_fit.json';driver.save(realized_path,[fit])
        artifact=lambda path:dict(path=str(path),sha256=driver.sha(path))
        spec=dict(source_manifest=artifact(source_manifest_path),source_contract=artifact(contract_path),
            realized_fit=artifact(realized_path),last_realized_state=artifact(state_path),
            windows=[dict(year=2007,fit=artifact(fit_path),root_receipt=artifact(root_path),
                accepted_forecast=artifact(accepted_path))])
        current_manifest=copy.deepcopy(source_manifest)
        current_manifest['file_sha256'].pop(str(old_driver));current_manifest['file_sha256'][str(new_driver)]=driver.sha(new_driver)
        current_manifest.update(reuse_forecast_jacobian=True,resume_history=spec,forecast_seconds=200)
        kwargs=dict(current_manifest=current_manifest,
            targets=[dict(decision_year='2007',period_tfr_arithmetic_mean='1.7')],tolerance=.005,
            case='A0',count=6,initial_checkpoint_sha256=initial_digest,source_root=source_root)
        return spec,kwargs,dict(contract_path=contract_path,fit_path=fit_path,
            realized_path=realized_path,state_path=state_path,state=state,fit=fit)

    def test_next_vintage_guess_advances_all_three_blocks(self):
        old=np.arange(1.,22.)
        shifted=driver.shift_forecast_coordinates(old,6).reshape(3,7)
        np.testing.assert_array_equal(shifted,np.array([[2,3,4,5,6,7,7],[9,10,11,12,13,14,14],[16,17,18,19,20,21,21]]))
        np.testing.assert_array_equal(old,np.arange(1.,22.))

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

    def test_optional_initial_coordinate_seed_is_exact_pinned_guess_only(self):
        default=np.full(21,.5)
        unchanged,receipt=driver.initial_coordinate_seed({},'A0',6,default)
        np.testing.assert_array_equal(unchanged,default);self.assertIsNot(unchanged,default)
        self.assertEqual(receipt['mode'],'default')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);root_path=d/'root_receipt.json';seed_path=d/'seed.json'
            prices=np.arange(1.,22.)
            def write_seed(*,seed_case='A0',coordinates=prices,root_case='A0',selection='final'):
                driver.save(root_path,dict(case=root_case,count=6,start_year=2007,
                    final=dict(prices=prices.tolist(),mapping_valid=True),
                    best=dict(prices=(prices+.1).tolist(),mapping_valid=True)))
                driver.save(seed_path,dict(case=seed_case,count=6,start_year=2007,
                    coordinates=np.asarray(coordinates).tolist(),label='numerical_guess_only',
                    source_root_receipt=dict(path=str(root_path),sha256=driver.sha(root_path)),
                    selection=selection))
                return dict(initial_coordinate_seeds={'A0_6':dict(
                    path=str(seed_path),sha256=driver.sha(seed_path))},
                    file_sha256={str(seed_path):driver.sha(seed_path)})
            manifest=write_seed();got,receipt=driver.initial_coordinate_seed(manifest,'A0',6,default)
            np.testing.assert_array_equal(got,prices)
            self.assertEqual(receipt['mode'],'pinned_numerical_guess')
            self.assertEqual(receipt['selection'],'final')
            self.assertEqual(receipt['driver']['sha256'],driver.sha(Path(driver.__file__).resolve()))
            mismatch=write_seed(seed_case='A+')
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(mismatch,'A0',6,default)
            mismatch=write_seed(coordinates=prices+.2)
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(mismatch,'A0',6,default)
            negative=prices.copy();negative[0]=-1.;bad=write_seed(coordinates=negative)
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(bad,'A0',6,default)
            best=write_seed(coordinates=prices+.1,selection='best')
            got,receipt=driver.initial_coordinate_seed(best,'A0',6,default)
            np.testing.assert_array_equal(got,prices+.1);self.assertEqual(receipt['selection'],'best')
            stale=write_seed();stale['initial_coordinate_seeds']['A0_6']['sha256']='0'*64
            stale['file_sha256'][str(seed_path)]='0'*64
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(stale,'A0',6,default)

    def test_initial_coordinate_seed_can_hold_each_source_block_tail(self):
        default=np.full(303,.5)
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);root_path=d/'root_receipt.json';seed_path=d/'seed.json'
            source=np.arange(1.,76.)
            expected=np.concatenate([np.pad(block,(0,76),'edge')
                for block in source.reshape(3,25)])
            def write_seed(*,root_count=24,seed_case='A0',coordinates=expected):
                driver.save(root_path,dict(case='A0',count=root_count,start_year=2007,
                    final=dict(prices=source.tolist(),mapping_valid=True),
                    best=dict(prices=source.tolist(),mapping_valid=True)))
                driver.save(seed_path,dict(case=seed_case,count=100,start_year=2007,
                    coordinates=np.asarray(coordinates).tolist(),label='numerical_guess_only',
                    source_root_receipt=dict(path=str(root_path),sha256=driver.sha(root_path)),
                    selection='final',source_count=24,rule='hold_last'))
                return dict(initial_coordinate_seeds={'A0_100':dict(
                    path=str(seed_path),sha256=driver.sha(seed_path))},
                    file_sha256={str(seed_path):driver.sha(seed_path)})
            got,receipt=driver.initial_coordinate_seed(write_seed(),'A0',100,default)
            np.testing.assert_array_equal(got,expected)
            self.assertEqual(receipt['seed_origin']['source_count'],24)
            self.assertEqual(receipt['seed_origin']['rule'],'hold_last')
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(
                write_seed(root_count=100),'A0',100,default)
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(
                write_seed(seed_case='A+'),'A0',100,default)
            corrupt=expected.copy();corrupt[25+24]=999.
            with self.assertRaises(ValueError):driver.initial_coordinate_seed(
                write_seed(coordinates=corrupt),'A0',100,default)

    def test_reuse_requires_a_finite_verified_same_width_root(self):
        matrix=np.arange(21*21,dtype=float).reshape(21,21)
        receipt=dict(finite_horizon_market_fiscal_converged=True,final_jacobian=matrix)
        result=NS(next_state=NS(),root_receipt=receipt)
        reused=driver.reusable_forecast_jacobian(result,6,True)
        np.testing.assert_array_equal(reused,matrix)
        self.assertIsNot(reused,matrix)
        self.assertIsNone(driver.reusable_forecast_jacobian(result,6,False))
        self.assertIsNone(driver.reusable_forecast_jacobian(
            NS(next_state=None,root_receipt=receipt),6,True))
        self.assertIsNone(driver.reusable_forecast_jacobian(
            NS(next_state=NS(),root_receipt=dict(receipt,final_jacobian=np.eye(3))),6,True))
        self.assertIsNone(driver.reusable_forecast_jacobian(
            NS(next_state=NS(),root_receipt=dict(receipt,fixed_asset_prices=np.ones(7))),6,True))

    def test_probe_seed_requires_accepted_exact_same_track_root(self):
        prices=np.arange(1.,22.);jacobian=np.eye(21)
        receipt=dict(converged=True,status='converged',
            finite_horizon_market_fiscal_converged=True,start_year=2007,case='A0',count=6,
            psi=.13,final_reproduction_max_abs=1e-12,final_damping=.25,
            final=dict(prices=prices,mapping_valid=True),best=dict(prices=prices.copy()),
            final_jacobian=jacobian)
        got_prices,got_jacobian,psi=jacobian_probe.accepted_seed(receipt,'A0',6,2e-10)
        np.testing.assert_array_equal(got_prices,prices)
        np.testing.assert_array_equal(got_jacobian,jacobian)
        self.assertEqual(psi,.13)
        with self.assertRaises(ValueError):
            jacobian_probe.accepted_seed(dict(receipt,case='A+'),'A0',6,2e-10)
        with self.assertRaises(ValueError):
            jacobian_probe.accepted_seed(dict(receipt,final_reproduction_max_abs=1e-4),'A0',6,2e-10)

    def test_fiscal_polish_seed_requires_complete_passing_housing_block(self):
        prices=np.arange(1.,22.);residual=np.r_[np.full(7,1e-5),np.full(7,3e-5),np.full(7,8e-4)]
        receipt=dict(start_year=2007,case='A+',count=6,psi=.13,final_reproduction_max_abs=0.,
            final=dict(prices=prices,residual=residual,mapping_valid=True),
            best=dict(prices=prices.copy(),residual=residual.copy(),mapping_valid=True))
        got,psi,fixed=jacobian_probe.fiscal_polish_seed(receipt,'A+',6,2e-10,2e-4)
        np.testing.assert_array_equal(got,prices);np.testing.assert_array_equal(fixed,prices[:7])
        self.assertEqual(psi,.13)
        bad=copy.deepcopy(receipt);bad['final']['residual'][3]=2e-4
        with self.assertRaises(ValueError):
            jacobian_probe.fiscal_polish_seed(bad,'A+',6,2e-10,2e-4)
        bad=copy.deepcopy(receipt);bad['final']['mapping_valid']=False
        with self.assertRaises(ValueError):
            jacobian_probe.fiscal_polish_seed(bad,'A+',6,2e-10,2e-4)

    def test_fixed_asset_prices_require_bounds_and_positive_pf_rents(self):
        pf=NS(rents_from_asset_prices=lambda prices,terminal,P:np.ones(len(prices)))
        fixed=driver.validated_fixed_asset_prices(np.arange(1.,8.),6,(.5,8.),pf,NS())
        np.testing.assert_array_equal(fixed,np.arange(1.,8.))
        for invalid in (np.ones(6),np.r_[np.ones(6),0.],np.r_[np.ones(6),9.]):
            with self.assertRaises(ValueError):
                driver.validated_fixed_asset_prices(invalid,6,(.5,8.),pf,NS())
        bad_pf=NS(rents_from_asset_prices=lambda prices,terminal,P:np.r_[np.ones(5),0.])
        with self.assertRaises(ValueError):
            driver.validated_fixed_asset_prices(np.arange(1.,8.),6,(.5,8.),bad_pf,NS())

    def test_automatic_polish_switch_predicate_requires_clear_housing_and_open_fiscal(self):
        prices=np.arange(1.,22.);eligible=dict(phase='iterate',mapping_valid=True,evaluation=2,
            prices=prices,residual=np.r_[np.full(7,1e-5),np.zeros(7),np.full(7,3e-4)])
        np.testing.assert_array_equal(driver.fiscal_polish_switch_prices(eligible,7,2e-4,8),prices[:7])
        already=copy.deepcopy(eligible);already['residual'][14:]=1e-5
        self.assertIsNone(driver.fiscal_polish_switch_prices(already,7,2e-4,8))
        bad_market=copy.deepcopy(eligible);bad_market['residual'][0]=2e-4
        self.assertIsNone(driver.fiscal_polish_switch_prices(bad_market,7,2e-4,8))
        late=copy.deepcopy(eligible);late['evaluation']=7
        self.assertIsNone(driver.fiscal_polish_switch_prices(late,7,2e-4,8))

    def test_source_hash_drift_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source';p.write_text('first')
            pins={str(p):driver.sha(p)};driver.verify_pins(pins)
            p.write_text('changed')
            with self.assertRaises(ValueError):driver.verify_pins(pins)

    def test_resume_loads_exact_contiguous_state_and_shifted_coordinates(self):
        with tempfile.TemporaryDirectory() as d:
            spec,kwargs,paths=self.resume_fixture(d)
            resumed=driver.load_resume_history(spec,**kwargs)
            self.assertEqual(resumed['completed_years'],[2007])
            self.assertTrue(driver._same_state(resumed['inherited'],paths['state']))
            np.testing.assert_array_equal(resumed['initial'].reshape(3,7),
                np.array([[2,3,4,5,6,7,7],[9,10,11,12,13,14,14],[16,17,18,19,20,21,21]]))

    def test_resume_rejects_wrong_case_source_hash_year_fit_and_state(self):
        mutations=('case','source','hash','year','fit','state')
        for mutation in mutations:
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as d:
                spec,kwargs,paths=self.resume_fixture(d)
                if mutation in ('case','source'):
                    contract=json.loads(paths['contract_path'].read_text())
                    contract['case']='A+' if mutation=='case' else contract['case']
                    contract['source_root']='/different/source' if mutation=='source' else contract['source_root']
                    driver.save(paths['contract_path'],contract)
                    spec['source_contract']['sha256']=driver.sha(paths['contract_path'])
                elif mutation=='hash':spec['realized_fit']['sha256']='0'*64
                elif mutation=='year':spec['windows'][0]['year']=2011
                elif mutation=='fit':
                    fit=dict(paths['fit'],model=1.71,gap=1.71-1.7)
                    driver.save(paths['fit_path'],fit);driver.save(paths['realized_path'],[fit])
                    spec['windows'][0]['fit']['sha256']=driver.sha(paths['fit_path'])
                    spec['realized_fit']['sha256']=driver.sha(paths['realized_path'])
                else:
                    state=copy.deepcopy(paths['state']);state.households.persons.heads[0]=99.
                    driver.checkpoint(paths['state_path'],state)
                    spec['last_realized_state']['sha256']=driver.sha(paths['state_path'])
                with self.assertRaises(ValueError):driver.load_resume_history(spec,**kwargs)

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
            self.assertNotIn('g_pre',kw)
            calls.append(kw['price'])
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
        fixed_prices=np.linspace(1.1,1.7,7);root_initials=[];adaptive=[]
        def root(**kw):
            if adaptive:
                adaptive.append(len(adaptive))
                if len(adaptive)==2:
                    self.assertEqual(kw['max_evaluations'],8)
                    x=kw['project'](kw['initial_prices']);first=kw['evaluate'](x)
                    kw['callback'](dict(evaluation=1,phase='initial',prices=x,
                        residual=first['residual'],mapping_valid=True,score=.01,elapsed_seconds=.1))
                    trial=x.copy();trial[:7]=fixed_prices;trial=kw['project'](trial)
                    second=kw['evaluate'](trial)
                    kw['callback'](dict(evaluation=2,phase='iterate',prices=trial,
                        residual=np.r_[np.full(7,1e-5),np.zeros(7),np.full(7,3e-4)],
                        mapping_valid=True,score=3e-4,elapsed_seconds=.2))
                    self.fail('Automatic fiscal switch callback did not interrupt the coupled root')
                self.assertEqual(kw['max_evaluations'],6);self.assertIsNone(kw['initial_jacobian'])
                raw=np.asarray(kw['initial_prices']).copy();raw[:7]=4.8;x=kw['project'](raw)
                np.testing.assert_array_equal(x[:7],fixed_prices)
                first=kw['evaluate'](x);second=kw['evaluate'](x)
                history=[dict(evaluation=1,phase='initial',prices=x,residual=first['residual'],
                    mapping_valid=True,score=.1,elapsed_seconds=.1),
                    dict(evaluation=2,phase='final',prices=x,residual=second['residual'],
                    mapping_valid=True,score=.1,elapsed_seconds=.2)]
                return dict(converged=True,status='converged',best=dict(prices=x,residual=first['residual'],
                    mapping_valid=True,payload=first['payload']),final=dict(prices=x,residual=second['residual'],
                    mapping_valid=True,payload=second['payload']),final_reproduction_max_abs=0.,evaluations=2,
                    elapsed_seconds=.2,history=history,final_jacobian=np.eye(21),final_damping=.5)
            root_initials.append(kw['initial_jacobian'])
            expected=np.diag(np.r_[np.full(7,-1.),np.full(14,-200.)])
            np.testing.assert_array_equal(kw['default_jacobian'],expected)
            raw=np.asarray(kw['initial_prices']).copy();raw[:7]=4.5
            x=kw['project'](raw)
            np.testing.assert_array_equal(x[:7],fixed_prices)
            first=kw['evaluate'](x);second=kw['evaluate'](x)
            self.assertEqual(first['residual'].shape,(21,))
            np.testing.assert_array_equal(first['residual'],second['residual'])
            np.testing.assert_array_equal(first['residual'],
                np.r_[np.full(6,.1),0.,np.full(6,.2),0.,np.full(6,.3),0.])
            return dict(converged=True,final=dict(prices=x,payload=second['payload']))
        rebated=NS(_runtime=lambda:runtime,evaluate_forecast=forecast,
            rebated_tax_accounts=lambda **kw:{},dated_residual=lambda **kw:np.array([.1,.2,.3]),
            stack_dated_residuals=lambda rows:np.asarray(rows).T.ravel(),
            first_period_state=lambda **kw:NS(year=2011))
        modules={'e5f_rebated_surprises':rebated,
            'e5f_closed_finite_boundary':NS(boundary_policy=boundary),
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
                audit=NS(),deadline=time.monotonic()+10,folder=Path(d)/'fixed',case='A0',
                initial_jacobian=None,fixed_asset_prices=fixed_prices)
            driver.solve_forecast(inherited=inherited,old=old,demographics=NS(),
                psi=.1,count=count,initial=np.r_[np.ones(14),np.full(7,.2)],controls=controls,
                audit=NS(),deadline=time.monotonic()+10,folder=Path(d)/'learned',case='A0',
                initial_jacobian=np.eye(21),fixed_asset_prices=fixed_prices)
            adaptive.append(0)
            polished,_=driver.solve_forecast(inherited=inherited,old=old,demographics=NS(),
                psi=.1,count=count,initial=np.r_[np.ones(14),np.full(7,.2)],
                controls=dict(controls,automatic_fiscal_polish=True),audit=NS(),
                deadline=time.monotonic()+10,folder=Path(d)/'automatic',case='A0')
        self.assertIsNone(root_initials[0]);np.testing.assert_array_equal(root_initials[1],np.eye(21))
        self.assertEqual(polished.root_receipt['evaluations'],4)
        self.assertTrue(polished.root_receipt['automatic_fiscal_polish']['switched'])
        self.assertEqual(polished.root_receipt['automatic_fiscal_polish']['total_actual_mappings'],4)
        self.assertEqual([row['evaluation'] for row in polished.root_receipt['root_phase_ledger']],[1,2,3,4])
        self.assertEqual([row['root_phase'] for row in polished.root_receipt['root_phase_ledger']],
            ['coupled','coupled','fiscal_polish','fiscal_polish'])
        np.testing.assert_array_equal(polished.root_receipt['fixed_asset_prices'],fixed_prices)
        self.assertEqual(len(calls),8)  # Two mappings in each of four root phases/runs.
        self.assertEqual(len(actual_calls),8)
        for value in actual_calls:np.testing.assert_array_equal(value,actual_g)
        self.assertTrue(result.root_receipt['finite_horizon_market_fiscal_converged'])
        self.assertFalse(result.root_receipt['horizon_verified'])
        self.assertEqual(result.next_state.year,2011)


if __name__=='__main__':unittest.main()
