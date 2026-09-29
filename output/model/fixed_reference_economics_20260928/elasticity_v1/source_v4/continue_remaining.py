#!/usr/bin/env python3
"""Budget-preserving continuation after the v2 forecast stop; no model logic here."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

SOURCE = Path('/work/elasticity_source/run_elasticity.py')
EXPECTED_PLAN_SHA = '6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f'
EXPECTED_JOB = '18801318'
EXPECTED_DEADLINE = 1790700222.6872504
ORDER = ('grid_control', 'grid_control_repeat', 'credit', 'credit_repeat')
GROUPS = (('reference_990','credit_990','reference_1010','credit_1010'),
          ('reference_980','credit_980','reference_1020','credit_1020'))

def require(ok, message):
    if not ok: raise RuntimeError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def read(path): return json.loads(Path(path).read_text())

def write(path, value):
    p=Path(path); tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'); tmp.replace(p)

def load_v2():
    sys.path.insert(0,str(SOURCE.parent))
    spec=importlib.util.spec_from_file_location('frozen_elasticity_v2',SOURCE)
    m=importlib.util.module_from_spec(spec); sys.modules[spec.name]=m; spec.loader.exec_module(m)
    return m

def terminal(job, contract):
    state=contract['original_job_terminal']
    require(state==dict(job_id=job,state='COMPLETED',exit_code='0:0'),
            'Host Slurm verification did not certify original job terminal')
    return state

def verify_original(old, plan, v2, contract):
    require(sha(plan)==EXPECTED_PLAN_SHA and contract['v2_plan_sha256']==EXPECTED_PLAN_SHA,
            'Original v2 plan identity differs')
    require(read(old/'launch.json')==contract['original_launch'] and
            sha(old/'launch.json')==contract['original_launch_sha256'], 'Original launch changed')
    launch=read(old/'launch.json')
    require(str(launch['slurm_job'])==EXPECTED_JOB and
            launch['deadline_epoch']==EXPECTED_DEADLINE and
            contract['original_deadline_epoch']==EXPECTED_DEADLINE,
            'Original job/deadline identity differs')
    require(launch['plan_sha256']==sha(plan) and launch['plan']==read(plan) and
            launch['deadline_epoch']==contract['original_deadline_epoch'], 'Original plan/deadline changed')
    require(sha(plan)==contract['v2_plan_sha256'] and sha(SOURCE)==contract['v2_driver_sha256'] and
            sha(__file__)==contract['continuation_driver_sha256'], 'Frozen source/plan changed')
    stop=read(old/'stopped.json')
    require(stop['status']=='stopped_before_next_case' and
            stop['reason']=='Conservative observed wall-time forecast exceeds remaining total budget' and
            stop['completed_cases']==4 and stop['next_case']==GROUPS[0][0],
            'Original did not stop solely for forecast after four q0 controls')
    require(not list(old.rglob('failure.json')), 'Original has model/controller failure')
    state=terminal(launch['slurm_job'],contract)
    latest=read(old/'latest_completed.json')
    records=latest['completed']
    require(latest['lifecycle_solves']==len(records)==4 and
            [r['case'] for r in records]==list(ORDER) and
            records==contract['completed_records'], 'Four completed records differ')
    manifest=read(v2.base.MANIFEST); names=sorted(manifest['standard_diagnostic_names'])
    for row in records:
        case=old/row['case']; receipt_path=case/'receipt.json'; receipt=read(receipt_path)
        require(sha(receipt_path)==row['receipt_sha256'] and receipt['status']=='passed' and
                receipt['plan_sha256']==sha(plan) and receipt['lifecycle_solves']==1 and
                receipt['standard_plot_count']==17, 'Completed q0 receipt changed')
        require(sha(case/'conditional_cohort_state.pkl.gz')==receipt['checkpoint']['sha256'],
                'Completed q0 checkpoint changed')
        require(len(list(csv.DictReader((case/'target_fit.csv').open())))==14 and
                len(list(csv.DictReader((case/'parameters.csv').open())))==31,
                'Completed q0 tables incomplete')
        require(sorted(p.name for p in (case/'standard_diagnostics').glob('*.png'))==names,
                'Completed q0 standard plots incomplete')
        if row['case'].endswith('_repeat'):
            require(receipt['exact_q0_repeat']['status']=='passed' and
                    receipt['exact_q0_repeat']['standard_plots_exact']==17,
                    'Fresh q0 repeat failed')
    return records,state

def forecast(records, cases):
    by_regime={regime:min(600.,1.25*max(r['wall_seconds'] for r in records if r['regime']==regime))
               for regime in ('reference','credit')}
    return sum(by_regime['reference' if name.startswith('reference_') else 'credit'] for name in cases)+60.,by_regime

def partial_one_percent(out, v2, records):
    cases={r['case']:read(out/r['case']/'receipt.json') for r in records}
    if not all(name in cases for name in ORDER+GROUPS[0]): return None
    pins={r['case']:r['receipt_sha256'] for r in records if r['case'] in ORDER+GROUPS[0]}
    require(all(sha(out/name/'receipt.json')==digest for name,digest in pins.items()),
            'Partial input receipt pin changed')
    rows=[]; elast=[]
    outcomes=('births_per_household','first_births','second_births','third_bin_entries',
              'ownership_rate','rooms_per_household')
    for regime,center in (('reference','grid_control'),('credit','credit')):
        for scope in ('impact','cohort'):
            for name in outcomes+(('completed_fertility','mean_first_birth_age') if scope=='cohort' else ()):
                vals=[]
                for factor,case in ((.99,regime+'_990'),(1.,center),(1.01,regime+'_1010')):
                    receipt=cases[case]
                    value=(receipt[name] if name in ('completed_fertility','mean_first_birth_age') else
                           receipt[scope+'_summary'][name])
                    if name in ('first_births','second_births','third_bin_entries'):
                        value/=receipt[scope+'_summary']['household_mass']
                    vals.append(float(value))
                    unit=('births per household' if name in ('births_per_household','first_births',
                          'second_births','third_bin_entries') else 'rooms per household'
                          if name=='rooms_per_household' else 'share' if name=='ownership_rate'
                          else 'children' if name=='completed_fertility' else 'years')
                    rows.append(dict(regime=regime,scope=scope,price_factor=factor,outcome=name,
                                     value=value,unit=unit,prescribed_price=receipt['price'],
                                     mapped_rent=receipt['rent']))
                lo,mid,hi=vals
                def slope(a,b,fa,fb):
                    return (math.log(b)-math.log(a))/(math.log(fb)-math.log(fa)) if a>0 and b>0 else None
                elast.append(dict(regime=regime,scope=scope,outcome=name,step=.01,
                    central_log_elasticity=slope(lo,hi,.99,1.01),
                    lower_one_sided_log_elasticity=slope(lo,mid,.99,1.),
                    upper_one_sided_log_elasticity=slope(mid,hi,1.,1.01)))
    for name,data in (('partial_comparison_1pct.csv',rows),('partial_elasticities_1pct.csv',elast)):
        with (out/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    summary=dict(status='partial_only',reference_label=v2.LABEL,price_range='±1 percent',
        original_deadline_epoch=EXPECTED_DEADLINE,input_plan_sha256=EXPECTED_PLAN_SHA,
        input_receipt_sha256=pins,
        comparison_sha256=sha(out/'partial_comparison_1pct.csv'),
        elasticities_sha256=sha(out/'partial_elasticities_1pct.csv'),
        complete_five_price_elasticities=False,comparison_rows=len(rows),elasticity_rows=len(elast),
        formula='[ln Y(q+) - ln Y(q-)] / [ln q+ - ln q-]; one-sided analogues',
        required_cases=list(ORDER+GROUPS[0]))
    write(out/'partial_1pct.json',summary)
    return summary

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--contract',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True); parser.add_argument('--original',type=Path,required=True)
    parser.add_argument('--plan',type=Path,required=True); args=parser.parse_args()
    out=args.output.resolve(); old=args.original.resolve(); plan=args.plan.resolve()
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only')
    require(not out.exists(),'Continuation output exists')
    contract=read(args.contract); require(contract['schema']=='block0506_elasticity_forecast_continuation_v4',
                                    'Wrong continuation contract')
    v2=load_v2(); p=v2.verify_plan(plan)
    records,state=verify_original(old,plan,v2,contract)
    require(time.time()<contract['original_deadline_epoch'],'Original absolute deadline passed')
    out.mkdir(parents=True)
    try:
        # Child mode requires these exact parent records and inherited launch bytes.
        (out/'launch.json').write_bytes((old/'launch.json').read_bytes())
        for row in records: (out/row['case']).symlink_to(old/row['case'],target_is_directory=True)
        write(out/'latest_completed.json',read(old/'latest_completed.json'))
        write(out/'best_so_far.json',read(old/'best_so_far.json'))
        write(out/'continuation.json',dict(status='started',original_job=state,
            original_deadline_epoch=contract['original_deadline_epoch'],original_solve_count=4,
            original_plan_sha256=sha(plan),continuation_contract_sha256=sha(args.contract)))
        deadline=contract['original_deadline_epoch']
        for group_index,group in enumerate(GROUPS):
            estimate,by_regime=forecast(records,group)
            available=deadline-time.time()
            if estimate>available:
                partial=partial_one_percent(out,v2,records)
                write(out/'stopped.json',dict(status='stopped_before_group',group_index=group_index,
                    next_group=list(group),estimate_seconds=estimate,available_seconds=available,
                    regime_forecast_seconds=by_regime,completed_cases=len(records),
                    all_five_price_elasticities_complete=False,partial_1pct=partial,no_retry=True))
                return
            for name in group:
                require(time.time()<deadline,'Original deadline passed before case')
                factor,regime=next((f,r) for n,f,r in v2.CASES if n==name)
                case_deadline=min(deadline,time.time()+p['case_seconds'])
                command=[sys.executable,str(SOURCE),'--plan',str(plan),'--output',str(out/name),
                         '--case',name,'--factor',repr(factor),'--regime',regime,'--deadline',repr(case_deadline)]
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
                         NUMEXPR_NUM_THREADS='1',NUMBA_NUM_THREADS='1',MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1')
                start=time.time()
                with (out/(name+'.log')).open('w') as log:
                    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                    try:
                        while process.poll() is None:
                            write(out/'progress.json',dict(phase='case_running',case=name,
                                time_epoch=time.time(),original_deadline_epoch=deadline,
                                completed_cases=len(records),case_deadline_epoch=case_deadline))
                            if time.time()>=case_deadline: raise TimeoutError('Original case/deadline: '+name)
                            time.sleep(2)
                        require(process.returncode==0,'Frozen v2 child failed; no retry: '+name)
                    finally:
                        if process.poll() is None:
                            os.killpg(process.pid,signal.SIGTERM)
                            try: process.wait(timeout=3)
                            except subprocess.TimeoutExpired:
                                os.killpg(process.pid,signal.SIGKILL); process.wait()
                receipt_path=out/name/'receipt.json'; receipt=read(receipt_path)
                require(receipt['status']=='passed' and receipt['lifecycle_solves']==1 and
                        receipt['plan_sha256']==sha(plan) and receipt['standard_plot_count']==17,
                        'Frozen v2 case incomplete')
                records.append(dict(case=name,factor=factor,regime=regime,receipt_sha256=sha(receipt_path),
                    abs_market_residual=abs(receipt['cohort_gates']['relative_market_residual']),
                    wall_seconds=time.time()-start))
                write(out/'latest_completed.json',dict(completed=records,lifecycle_solves=len(records),
                    remaining_cases=12-len(records)))
                best=min(records,key=lambda r:r['abs_market_residual'])
                write(out/'best_so_far.json',dict(status='case_completed',**best,
                    criterion='smallest absolute prescribed-price housing excess; no GE certification'))
                require(len(records)<=12,'Original solve cap exceeded')
            if group_index==0: partial_one_percent(out,v2,records)
        require(len(records)==12 and time.time()<deadline,'Complete cases exceed original budget')
        comparison=v2.make_comparison(out,records)
        write(out/'completed.json',dict(status='passed',reference_label=v2.LABEL,
            elapsed_seconds=time.time()-contract['original_launch']['started_epoch'],
            lifecycle_solves=12,
            original_solves=4,continuation_solves=8,original_deadline_epoch=deadline,
            normalization_performed=False,market_clearing_certified=False,
            plan_sha256=sha(plan),**comparison))
    except BaseException as exc:
        if out.exists(): write(out/'failure.json',dict(status='failed',error=str(exc),
            traceback=traceback.format_exc(),retries=0,time_epoch=time.time()))
        raise

if __name__=='__main__': main()
