"""Bounded preparation: pure tests and two no-shock native maps; no policy path."""
from pathlib import Path
import json
import os
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE/os.environ.get('E5F_TRANSITION_SOURCE_SUBDIR','source')
sys.path.insert(0,str(SOURCE))
import run_e5f_preference_transition as engine


def main():
    engine.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only')
    out=HERE/'runs'/os.environ['SLURM_JOB_ID']
    out.mkdir(parents=True,exist_ok=False)
    sources={name:engine.sha(SOURCE/name) for name in engine.SOURCE_NAMES}
    engine.check_sources(sources)
    planned=dict(schema='block0506_pf_preparation_v1',reference_label=engine.LABEL,
        reference_manifest_sha256=engine.MANIFEST_SHA,source_pins=sources,
        harness_sha256=engine.sha(__file__), maximum_actual_policy_calls=14,total_seconds=900,
        workers=1,threads=1,memory_gib=24,cache_bytes=2*1024**3,
        cases=[dict(name='one_date_cached',horizon=1,housing='elastic_reference'),
               dict(name='six_date_fixed_stock_cached',horizon=6,housing='fixed_stock')],
        preference_shocks=False,production_transition=False,no_retry=True,
        economics='unchanged reference preferences/credit; fixed-stock plumbing at exact baseline supply')
    engine.write(out/'plan.json',planned)
    started=time.monotonic()
    def deadline(*_):raise TimeoutError('900-second preparation budget exhausted')
    signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,900)
    completed=[]
    try:
        tests=['test_e5f_preference_transition','test_e5f_four_shock_acceleration','test_e5f_exact_policy_cache']
        env=dict(os.environ,PYTHONPATH=str(SOURCE))
        with (out/'tests.log').open('w') as log:
            tested=subprocess.run([sys.executable,'-m','unittest','-v',*tests],cwd=SOURCE,env=env,
                stdout=log,stderr=subprocess.STDOUT,timeout=90)
        engine.require(tested.returncode==0,'Pure contract/cache/root tests failed; inspect tests.log')
        engine.write(out/'progress.json',dict(phase='pure_tests_passed',epoch=time.time()))
        m,packet,evaluator=engine.load_reference(out/'reference')
        import numpy as np
        P=packet['parameters'];q=float(packet['solution'].p_eq[0])
        endpoint=dict(price=q,population_scale=1.)
        actual_calls=0
        for case in planned['cases']:
            engine.write(out/'progress.json',dict(phase=case['name'],epoch=time.time()))
            result,record=engine.mapping(packet,evaluator,packet,endpoint,
                np.full(case['horizon'],q),np.full(case['horizon'],P.pension),
                np.full(case['horizon'],P.psi_child),case['housing'],out/case['name'],planned['cache_bytes'])
            actual_calls+=record['cache']['actual_solves']
            engine.require(actual_calls<=planned['maximum_actual_policy_calls'],'Policy-call budget exceeded')
            gates=dict(scientific=all(record['gates'].values()),
                housing=max(abs(x) for x in record['market_residual'])<=2e-4,
                fiscal=max(abs(x) for x in record['fiscal_residual'])<=1e-6,
                retained_credit=P.native_due_stayer_credit and not getattr(P,'native_solvency_credit',False),
                cached=record['cache']['hits']>0)
            if case['horizon']==1:
                old=engine.read(HERE.parent/'preparation_v1/runs/18739319/no_shock_1/rows.json')
                gaps={key:abs(float(record['rows'][0][key])-float(old[0][key])) for key in old[0]
                      if isinstance(old[0][key],(int,float)) and key!='calendar_year'}
                gates['uncached_saved_rows_exact']=all(v==0 for v in gaps.values())
                record['uncached_comparison']=gaps
            else:
                pf=evaluator.rt['primitive'].pf
                initial=pf.stationary_initial_state(packet['stationary_g_pre'],
                    float(packet['stationary_g_pre'][:,:,:,0].sum()),float(packet['evaluation'].births),P,1/2.1)
                record['adjusted_queue_max_change']=float(np.max(np.abs(
                    pf.birth_queue_values(result.terminal_state.scheduled_entries)-pf.birth_queue_values(initial.scheduled_entries))))
                record['raw_queue_max_change']=float(np.max(np.abs(
                    pf.birth_queue_values(result.terminal_state.scheduled_raw_entries)-pf.birth_queue_values(initial.scheduled_raw_entries))))
                gates['entry_lags_covered']=case['horizon']*P.period_years>=20
                record['terminal_check']=engine.terminal_checks(packet,evaluator,packet,endpoint,result,
                    np.full(case['horizon'],P.psi_child),dict(
                        terminal_tolerances={key:1e-6 for key in engine.TERMINAL_KEYS},
                        raw_queue_relative_tolerance=1e-6))
                gates['terminal_population_distribution_both_queues']=record['terminal_check']['all_checks_pass']
            engine.write(out/case['name']/'comparison.json',dict(case=case,gates=gates,record=record))
            engine.require(all(gates.values()),'No-shock/cache gate failed; no retry')
            completed.append(dict(case=case,gates=gates,record=record))
            engine.write(out/'latest_completed.json',completed[-1])
            engine.write(out/'best_so_far.json',dict(completed_cases=[c['case']['name'] for c in completed]))
            del result
        for kind in ('one_permanent','four_announced'):
            plan=engine.draft_plan(kind);plan['source_pins']=sources
            engine.write(out/(kind+'_draft_plan.json'),plan)
        engine.check_sources(sources)
        engine.require(engine.sha(engine.MANIFEST)==engine.MANIFEST_SHA,'Reference manifest changed')
        engine.write(out/'readiness.json',dict(planned,status='PASS',tests_passed=True,native_smoke_passed=True,
            completed=completed,actual_policy_calls=actual_calls,seconds=time.monotonic()-started,
            execution_enabled=False,launch_inputs_pending=['preference levels/provenance','terminal endpoint',
                'production horizon/budget/numerical settings','later explicit execution authorization']))
    except BaseException as exc:
        engine.write(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),completed=completed,
            seconds=time.monotonic()-started,production_transition=False))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)


if __name__=='__main__':main()
