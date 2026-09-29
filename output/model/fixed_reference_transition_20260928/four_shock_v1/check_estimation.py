"""Test-only harness: synthetic fits and unchanged-preference native integration."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'source_estimation_v1'
sys.path[:0]=[str(SOURCE),'/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools']
import run_e5f_preference_estimation as estimator
engine=estimator.inner


def carried_gaps(actual_rows,expected_rows,offset):
    gaps=[]
    engine.require(len(actual_rows)==len(expected_rows),'Carried forecast length differs')
    for actual,expected in zip(actual_rows,expected_rows):
        engine.require(actual['period']+offset==expected['period'],'Wrong within-forecast period offset')
        for key,value in expected.items():
            if key!='period' and isinstance(value,(int,float)):
                gaps.append(abs(float(actual[key])-float(value))/max(1.,abs(float(value))))
    return gaps


def write_drafts(out,pins,contract):
    for kind in ('four_successive','one_permanent'):
        draft=estimator.draft_plan(kind);draft.update(source_pins=pins,target_contract=contract,
            readiness_receipt=dict(path=str(out/'readiness.json'),sha256=engine.sha(out/'readiness.json')))
        engine.write(out/(kind+'_draft_plan.json'),draft)


def recheck_carried(previous,out,pins,contract):
    """Recheck only the saved unchanged-baseline handoff after a test assertion repair."""
    import gzip
    import pickle
    started=time.monotonic();old=HERE/'estimation_tests'/previous
    engine.require(engine.read(old/'plan.json')['estimator_sources']==pins,'Cannot reuse evidence after a model-source change')
    engine.require(engine.read(old/'failure.json')['error']=='Inherited continuation changed the original forecast',
                   'This repair mode handles only the counter-offset assertion')
    engine.write(out/'plan.json',dict(test_only=True,previous=previous,maximum_policy_calls=10,
        total_seconds=300,historical_fit=False,preference_changes=False,estimator_sources=pins,
        harness_sha256=engine.sha(__file__),repair='Compare relative counters with offset; recheck all economic quantities and complete terminal state'))
    names=['test_e5f_preference_estimation','test_e5f_preference_transition','test_e5f_four_shock_acceleration','test_e5f_exact_policy_cache']
    with (out/'tests.log').open('w') as log:
        tested=subprocess.run([sys.executable,'-m','unittest','-v',*names],cwd=SOURCE,
            env=dict(os.environ,PYTHONPATH=os.pathsep.join(sys.path[:2])),stdout=log,stderr=subprocess.STDOUT,timeout=90)
    engine.require(tested.returncode==0,'Pure tests failed')
    m,packet,runtime=engine.load_reference(out/'reference')
    import numpy as np
    endpoint=engine.read(old/'native/endpoints/endpoint_1/receipt.json')
    with gzip.open(engine.pinned(endpoint['checkpoint']),'rb') as stream:terminal=pickle.load(stream)
    saved_state=engine.read(old/'native/accepted_2007/checkpoint.json')
    with gzip.open(engine.pinned(saved_state),'rb') as stream:inherited=pickle.load(stream)
    with gzip.open(old/'native/forecast/latest_state.pkl.gz','rb') as stream:forecast=pickle.load(stream)
    receipt=engine.read(old/'native/forecast/receipt.json')
    root=engine.read(old/'native/forecast/root.json')
    prefix=engine.read(old/'native/accepted_2007/replay.json')
    engine.require(root['converged'] and receipt['root_and_terminal_pass'] and all(prefix['gates'].values()),'Original certified checks missing')
    psi=float(packet['parameters'].psi_child)
    engine.require(np.all(forecast['psi_path']==psi) and endpoint['psi_child']==psi and inherited['year']==2011,
                   'Test must retain the reference preference and inherited 2011 state')
    p=estimator.draft_plan();p['budget'].update(total_seconds=300,candidate_seconds=300,mapping_seconds=200)
    native=estimator.NativeEstimator(p,out/'native',m,packet,runtime);native.deadline=started+300
    carried,record=native.guarded(200,lambda:engine.mapping(packet,runtime,terminal,endpoint,
        forecast['prices'][1:],forecast['pensions'][1:],np.full(5,psi),'fixed_stock',out/'carried',
        p['path']['cache_max_bytes'],initial_state=inherited['households'],start_year=2011,measure_fertility=True))
    gaps=carried_gaps(record['rows'],receipt['rows'][1:],1)
    pf=runtime.rt['primitive'].pf
    g_gap=float(np.max(np.abs(carried.terminal_state.g_pre-forecast['terminal_state'].g_pre)))
    queues={name:float(np.max(np.abs(pf.birth_queue_values(getattr(carried.terminal_state,name))-
        pf.birth_queue_values(getattr(forecast['terminal_state'],name))))) for name in ('scheduled_entries','scheduled_raw_entries')}
    engine.require(all(record['gates'].values()) and max(gaps)<=2e-10 and g_gap<=1e-12 and max(queues.values())<=1e-12,
                   'Carried-state replay changed the original forecast or queues')
    for name,digest in pins.items():engine.require(engine.sha(SOURCE/name)==digest,'Source changed during recheck')
    engine.write(out/'readiness.json',dict(status='PASS',test_only=True,estimator_sources=pins,tests_passed=True,
        native_endpoint_and_fertility_smoke_passed=True,historical_fit=False,preference_changes=False,
        seconds=time.monotonic()-started,maximum_policy_calls=10,original_check_job=previous,
        original_plan=dict(path=str(old/'plan.json'),sha256=engine.sha(old/'plan.json')),
        endpoint=endpoint,terminal=receipt['terminal'],prefix_reproduced=True,
        carried_row_maximum_relative_gap=max(gaps),carried_queue_absolute_gaps=queues,
        carried_distribution_maximum_gap=g_gap,cache=record['cache'],
        execution_enabled=False,production_horizon_verified=False))
    write_drafts(out,pins,contract)


def main(previous=None):
    engine.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only')
    out=HERE/'estimation_tests'/os.environ['SLURM_JOB_ID'];out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();pins={name:engine.sha(SOURCE/name) for name in estimator.SOURCES}
    contract=estimator.target_contract(HERE/'estimation_inputs/empirical_blocks.csv',HERE/'estimation_inputs/annual_fertility_2007_2023.csv')
    if previous is not None:
        engine.require(previous.isdigit(),'Explicit prior test job required')
        return recheck_carried(previous,out,pins,contract)
    engine.write(out/'plan.json',dict(test_only=True,reference_manifest_sha256=engine.MANIFEST_SHA,
        estimator_sources=pins,harness_sha256=engine.sha(__file__),historical_fit=False,
        preference_changes=False,maximum_policy_calls=40,total_seconds=1100,
        workers=1,threads=1,memory_gib=24,no_retry=True,
        cases=['synthetic unit tests','unchanged-psi endogenous endpoint',
               'six-date unchanged baseline root','original-vintage one-period replay','carried-state five-date replay']))
    try:
        names=['test_e5f_preference_estimation','test_e5f_preference_transition',
               'test_e5f_four_shock_acceleration','test_e5f_exact_policy_cache']
        env=dict(os.environ,PYTHONPATH=os.pathsep.join(sys.path[:2]))
        with (out/'tests.log').open('w') as log:
            result=subprocess.run([sys.executable,'-m','unittest','-v',*names],cwd=SOURCE,
                env=env,stdout=log,stderr=subprocess.STDOUT,timeout=90)
        engine.require(result.returncode==0,'Pure tests failed; inspect tests.log')
        engine.write(out/'progress.json',dict(phase='pure_tests_passed',epoch=time.time()))
        m,packet,runtime=engine.load_reference(out/'reference')
        import numpy as np
        psi=float(packet['parameters'].psi_child);q=float(packet['solution'].p_eq[0])
        plan=estimator.draft_plan();plan.update(source_pins=pins,target_contract=contract,horizons=[6,8])
        plan['budget'].update(total_seconds=1100,candidate_seconds=1000,endpoint_seconds=500,
            mapping_seconds=200,path_seconds=400,maximum_policy_calls=40)
        plan['endpoint']['max_evaluations']=2;plan['path']['max_evaluations']=2
        native=estimator.NativeEstimator(plan,out/'native',m,packet,runtime)
        native.deadline=started+1100;native.candidate_deadline=native.deadline
        engine.write(out/'progress.json',dict(phase='unchanged_endpoint',epoch=time.time()))
        terminal,endpoint=native.endpoint(psi)
        engine.require(endpoint['psi_child']==psi and abs(endpoint['price']/q-1)<=1e-12,'Test changed reference preferences or price')
        engine.write(out/'progress.json',dict(phase='unchanged_six_date_root',epoch=time.time()))
        receipt,latest=native.path(psi,terminal,endpoint,6,out/'native/forecast')
        engine.require(np.all(np.array([r['psi_child'] for r in receipt['rows']])==psi),'Test changed preferences')
        observed=float(receipt['fertility'][0]['period_tfr_topcode_adjusted'])
        # Synthetic acceptance at the baseline moment exercises state transfer,
        # without fitting any historical target or calling the shock optimizer.
        native.latest=dict(psi=psi,latest=latest,terminal=terminal,endpoint=endpoint,
            summary={'payload':{'candidate':0,'models':[observed]}})
        synthetic=dict(converged=True,root={'final':dict(mapping_valid=True,prices=[psi],payload={'candidate':0})},
            parameter=dict(estimate=psi,lower=psi*.01,upper=psi*2,near_bound=False))
        engine.write(out/'progress.json',dict(phase='original_forecast_replay',epoch=time.time()))
        first=native.advance(0,synthetic,diagnostics=False)
        engine.require(native.year==2011 and native.inherited is not None,'First period did not advance')
        from types import SimpleNamespace
        boundary=dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=latest['result'].values[6])))
        carried,record=native.guarded(200,lambda:engine.mapping(packet,runtime,boundary,endpoint,
            latest['prices'][1:],latest['pensions'][1:],np.full(5,psi),'fixed_stock',out/'native/carried_five_dates',
            plan['path']['cache_max_bytes'],initial_state=native.inherited,start_year=2011,measure_fertility=True))
        engine.require(all(record['gates'].values()),'Carried-state scientific check failed')
        gaps=carried_gaps(record['rows'],receipt['rows'][1:],1)
        engine.require(max(gaps)<=2e-10,'Inherited continuation changed the original forecast')
        pf=runtime.rt['primitive'].pf
        engine.require(np.allclose(carried.terminal_state.g_pre,latest['result'].terminal_state.g_pre,rtol=0,atol=1e-12),
                       'Inherited distribution differs from six-date forecast')
        queue_gaps={}
        for name in ('scheduled_entries','scheduled_raw_entries'):
            queue_gaps[name]=float(np.max(np.abs(pf.birth_queue_values(getattr(carried.terminal_state,name))-
                pf.birth_queue_values(getattr(latest['result'].terminal_state,name)))))
        engine.require(max(queue_gaps.values())<=1e-12,'Inherited birth queues differ from six-date forecast')
        for name,digest in pins.items():engine.require(engine.sha(SOURCE/name)==digest,'Source changed during test')
        # Six native maps are bounded above by 40 Bellman calls before cache reuse.
        calls=2+2+2*6*2+2+2*5
        engine.require(calls==40,'Test solve-count accounting changed')
        ready=dict(status='PASS',test_only=True,estimator_sources=pins,tests_passed=True,
            native_endpoint_and_fertility_smoke_passed=True,historical_fit=False,preference_changes=False,
            seconds=time.monotonic()-started,maximum_policy_calls=calls,endpoint=endpoint,
            native_fertility=observed,terminal=receipt['terminal'],prefix_reproduced=True,
            carried_row_maximum_relative_gap=max(gaps),carried_queue_absolute_gaps=queue_gaps,
            carried_distribution_maximum_gap=float(np.max(np.abs(carried.terminal_state.g_pre-latest['result'].terminal_state.g_pre))),
            execution_enabled=False,production_horizon_verified=False)
        engine.write(out/'readiness.json',ready)
        write_drafts(out,pins,contract)
    except BaseException as exc:
        engine.write(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),
            seconds=time.monotonic()-started,historical_fit=False,preference_changes=False));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--recheck-carried')
    main(parser.parse_args().recheck_carried)
