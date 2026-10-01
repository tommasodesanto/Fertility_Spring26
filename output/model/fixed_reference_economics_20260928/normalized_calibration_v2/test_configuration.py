"""Zero-lifecycle checks of the new normalization and actual 24-chain controller."""
from pathlib import Path
from types import SimpleNamespace as NS
import ast,copy,json,tempfile,time
import numpy as np
import normalized_objective as normal
HERE=Path(__file__).resolve().parent
import run_psi as run

def test_normalization():
    P=NS(H0=np.array([9.]),user_cost_rate=.04,r_bar=np.array([.03]),xi_supply=np.array([.7]),N_target=1.,tol_eq=1e-4)
    sol=NS(total_mass=1.,housing_supply=np.array([20.]),V=np.array([4.,8.]))
    ev=NS(demand_by_loc=np.array([5.7]),g_current=np.array([.4,.6]),supply_by_loc=np.array([20.]),relative_market_residual=1.)
    cal=NS(HousingSupplyRule=lambda *a:a);ctx={'expected_parameters':{'H0':9.},'P':P};live={'sol':sol,'P':P}
    before=sol.V.copy();c,q,e,s=normal.normalize_observation(ctx,live,P,ev,cal,np.array([.72]),True)
    assert q is not P and P.H0[0]==9. and c['expected_parameters']['H0']==q.H0[0]
    np.testing.assert_allclose(q.H0*(q.user_cost_rate*.72/q.r_bar)**q.xi_supply,e.demand_by_loc,rtol=0,atol=1e-15)
    assert sol.aggregate_housing_supply==e.supply_by_loc.sum() and sol.aggregate_housing_excess==e.demand_by_loc.sum()-e.supply_by_loc.sum()
    np.testing.assert_array_equal(sol.V,before)
    ev.demand_by_loc=np.array([1000.])
    normal.normalize_observation(ctx,live,P,ev,cal,np.array([.72]),False)
    try:normal.normalize_observation(ctx,live,P,ev,cal,np.array([.72]),True)
    except normal.H0BoundError:pass
    else:raise AssertionError('Missing accepted-root H0 constraint')
    ev.g_current=np.array([1.01]);ev.demand_by_loc=np.array([5.7])
    try:normal.normalize_observation(ctx,live,P,ev,cal,np.array([.72]),True)
    except RuntimeError:pass
    else:raise AssertionError('Missing unit-mass gate')

def test_source_variants():
    source=(HERE.parent/'utility_calibration_round1_v1/phase_b_pilot.py').read_text();tree=ast.parse(source)
    fn={n.name:ast.get_source_segment(source,n)+'\n' for n in tree.body if isinstance(n,ast.FunctionDef)}
    full=normal.observer_variant(fn['observe_price']);compile(full,'normalized_observer','exec')
    native=(HERE.parent/'utility_floor_round2_v1/runner.py').read_text();nt=ast.parse(native)
    evaluator=next(ast.get_source_segment(native,n)+'\n' for n in nt.body if isinstance(n,ast.FunctionDef) and n.name=='native_evaluator')
    fast,root,_=normal.fast.variants(fn['observe_price'],fn['run_phase_b'],evaluator)
    compile(normal.observer_variant_no_plot(fast),'normalized_fast_observer','exec')
    assert 'population = 1.0' in full and 'H0_derived' in full and 'H0_derived' in normal.observer_variant_no_plot(fast)
    # Original root numerical source kept verbatim before repeat-removal boundary.
    repeat='    repeat, repeat_live = trial(root["price"], "selected_repeat", repeat=True, reason="fresh exact selected-price repeat")'
    assert root.split(repeat)[0].startswith(fn['run_phase_b'].split(repeat)[0])

def test_actual_loops():
    config=run.CONFIG;coords=tuple(config['verified_incumbent_retained']['parameters']);bounds=dict(config['bounds']);bounds['psi_child']=config['psi_bounds']
    seen=set();records=[]
    with tempfile.TemporaryDirectory() as tmp:
      for design in config['nearby_starts']:
        chain=design['chain'];seed=design['parameters'];seen.add(tuple(seed[k] for k in coords))
        for k in coords:assert bounds[k][0]<=seed[k]<=bounds[k][1]
        def evaluate(label,point,end):
            assert end<=time.time()+10800-899
            return dict(status='passed',residual=[(point[k]-seed[k])/(bounds[k][1]-bounds[k][0])+.03 for k in coords],lifecycle_solves=0,report=str(Path(tmp)/label),population=1.,H0_derived=5.8)
        result=run.optimize(Path(tmp)/str(chain),seed,bounds,coords,evaluate,time.time()+10800,toy=True,maxeval=150,chain=chain)
        assert result['objective_calls']<=150 and result['lifecycle_solves']==0 and result['completed_full_ge']>10
        records.append({'chain':chain,'objective_calls':result['objective_calls'],'computed_mock_cases':result['completed_full_ge']})
      stop=run.optimize(Path(tmp)/'reserve',seed,bounds,coords,evaluate,time.time()+899,toy=True,maxeval=150)
      assert stop['completed_full_ge']==0 and stop['search_stop_reason']=='three_hour_actual_start_final_reserve'
      def reject(*a):return dict(status='inadmissible_numerical',reason='Derived normalized H0 outside original [.2,80] bounds',derived_H0_bound_rejection=True,lifecycle_solves=3)
      result=run.optimize(Path(tmp)/'reject',seed,bounds,coords,reject,time.time()+10800,toy=True,maxeval=1)
      assert result['selected'] is None and result['lifecycle_solves']==3
      def fatal(*a):raise RuntimeError('accounting contract failure')
      try:run.optimize(Path(tmp)/'fatal',seed,bounds,coords,fatal,time.time()+10800,toy=True,maxeval=1)
      except RuntimeError as exc:assert str(exc)=='accounting contract failure'
      else:raise AssertionError('Fatal failures swallowed')
    assert len(seen)==24
    return records

def main():
    test_normalization();test_source_variants();records=test_actual_loops()
    for p in HERE.glob('*.py'):compile(p.read_text(),str(p),'exec')
    receipt=dict(status='passed_zero_model_solves',chains=24,distinct_seed_vectors=24,lifecycle_solves=0,normalization_arithmetic=True,accepted_root_H0_bound=True,unit_population_gate=True,full_and_fast_variants=True,budget_reserve_fatal_constraint_gates=True,mock_controller=records)
    (HERE/'configuration_check.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='mock_controller'}))
if __name__=='__main__':main()
