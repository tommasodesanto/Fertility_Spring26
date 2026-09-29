#!/usr/bin/env python3
"""Two fixed diagnostic evaluations, supervised on Torch; never a resumed search.

Smoke uses the exact subprocess, receipt, validation and process-group timeout
paths with conspicuously synthetic artifacts. Run requires pinned smoke approval.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parents[3]/'code/model/tools'))
import run_e5f_fertility_identification as original
from adapter import apply as apply_adapter

core = original.core
read, write, sha, canon = core.read, core.write, core.sha, core.canon
ARMS = ('retained_start', 'predicted_start')


def require(condition, message):
    if not condition: raise RuntimeError(message)


def verify_fingerprint(p,c,obj):
    require(p['target_weight_fingerprint'] == c['lanes']['primary']['target_weight_fingerprint'] == canon(obj['target_rows']), 'Target/weight fingerprint changed')


def verify(path):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm required')
    require(not sys.flags.optimize, 'Assertions must remain enabled')
    p = read(path)
    require(p['schema'] == 'e5f_numerical_pair_v1', 'Wrong pair schema')
    require(sha(path) == os.environ.get('EXPECTED_NUMERICAL_PAIR_SHA256'), 'Explicit pair contract pin differs')
    for pin in [p[k] for k in ('original_contract','plan','predictions','reference_receipt','reference_fit','reference_parameters','objective','original_source_manifest')] + list(p['files'].values()):
        require(sha(pin['path']) == pin['sha256'], 'Changed pair input: '+pin['path'])
    require(Path(p['files']['runner.py']['path']).resolve() == Path(__file__).resolve(), 'Wrong executing driver')
    require(Path(p['files']['adapter.py']['path']).resolve() == Path(apply_adapter.__code__.co_filename).resolve(), 'Wrong executing adapter')
    c, objs = original.verify(Path(p['original_contract']['path']))
    plan = read(p['plan']['path'])
    require(p['objective'] == c['lanes']['primary']['objective'], 'Objective pin changed')
    require(p['original_source_manifest'] == c['source_manifest'], 'Source manifest pin changed')
    verify_fingerprint(p,c,objs['primary'])
    require(plan['specification_changes'] == [] and p['specification_changes'] == [], 'Economic specification changed')
    require(plan['objective_count'] == plan['workers'] == 2, 'Pair dimensions changed')
    require(plan['objective_cap_seconds'] == 1800 and plan['whole_experiment_cap_seconds'] == 4200, 'Time budgets changed')
    require(plan['maximum_stationary_solves_per_objective'] == c['normalization']['maximum_stationary_solves'] == 23, 'Solve cap changed')
    require(plan['maximum_stationary_solves_total'] == 46, 'Global solve cap changed')
    require(c['normalization']['initial_psi'] == .14281100340255604 and c['normalization']['initial_step'] == plan['initial_bracket_step'] == .005, 'Retained normalization start changed')
    require(plan['predicted_initial_psi_full'] == .13052783066857948, 'Predicted start changed')
    require(set(plan['trial_point']) == set(core.FREE), 'Ten-dimensional proposal required')
    for row in objs['primary']['parameter_restrictions']:
        value = plan['trial_point'][row['parameter']]
        require(math.isfinite(value) and row['lower'] <= value <= row['upper'], 'Proposal violates bound')
    require(plan['paired_consistency'] == dict(failure_action='report sensitivity; do not claim equivalent faster normalization or promote', loss_absolute_difference=.05, price_relative_difference=.0001, psi_absolute_difference=.0001, scored_max_abs_difference_in_working_scale=.01, status='proposed numerical equivalence screens, additional to unchanged scientific gates', validation_absolute_difference_rule='.01*max(1,abs(reference),abs(target))'), 'Paired screens changed')
    return p, c, objs, plan


def arm_contract(c, plan, arm):
    changed = copy.deepcopy(c)
    changed['normalization']['initial_psi'] = c['normalization']['initial_psi'] if arm == ARMS[0] else plan['predicted_initial_psi_full']
    return changed


def identity(p, c, plan):
    # The original contract did not authorize a new normalization guess. This
    # explicit diagnostic identity records both its ancestry and the new adapter.
    return canon(dict(pair=p, point=plan['trial_point'], normalization=c['normalization']))


def expected_request(p, c, plan, arm, contract_sha, folder, deadline, synthetic):
    ac = arm_contract(c, plan, arm)
    ctx = dict(candidate_id=folder.name, stage='initial', contract_sha256=contract_sha,
               source_sha256=p['files']['runner.py']['sha256'], target_sha256=p['objective']['sha256'], point_sha256=canon(plan['trial_point']))
    return dict(id=folder.name, arm=arm, lane='primary', point=plan['trial_point'], context=ctx,
                contract_sha256=contract_sha, scientific_candidate_id=core.identity(ac,'primary',plan['trial_point']),
                diagnostic_identity=identity(p,ac,plan), normalization_inputs=ac['normalization'],
                controller_pid=os.getpid(), deadline_epoch=deadline, graphs=True, synthetic=synthetic)


def extra_receipt(receipt, p, c, plan, req):
    receipt.update(contract_sha256=req['contract_sha256'], scientific_identity=core.science(c),
        scientific_candidate_id=req['scientific_candidate_id'], lane='primary',
        target_weight_fingerprint=p['target_weight_fingerprint'], diagnostic_identity=req['diagnostic_identity'],
        original_contract_sha256=p['original_contract']['sha256'], pair_adapter_sha256=p['files']['adapter.py']['sha256'],
        pair_plan_sha256=p['plan']['sha256'], diagnostic_economic_changes=plan['economic_changes'],
        diagnostic_specification_changes=[], synthetic_fixture=bool(req['synthetic']),
        promotion='Not adopted; final exact repeats not performed')


def fixture(folder, p, c, plan, req):
    """Synthetic values only, explicitly labelled. Never imported by main run."""
    behavior = req['synthetic']
    if behavior == 'timeout':
        child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(300)'])
        write(folder/'grandchild.json', dict(pid=child.pid))
        time.sleep(300)
    if behavior in ('fatal','inadmissible'):
        write(folder/'failure.json', dict(context=req['context'], status=behavior,
            classification='explicit_economic_gate' if behavior=='inadmissible' else 'unknown_or_integrity_failure',
            authenticated=behavior=='inadmissible', error='Injected synthetic '+behavior))
        return 2
    case = folder/'case'; case.mkdir()
    receipt = read(p['reference_receipt']['path'])
    receipt['point'] = req['point']; receipt['normalization_inputs'] = c['normalization']
    receipt['objective_stationary_solves'] = 1
    receipt['normalization']['stationary_solves'] = 1
    receipt['normalization']['stationary_solve_seconds'] = 0.
    receipt['objective_stationary_solve_seconds'] = 0.
    write(case/'stationary_solves.json', [dict(status='completed', psi_child=c['normalization']['initial_psi'], seconds=0., synthetic=True)])
    (case/'initial_state.pkl.gz').write_bytes(b'SYNTHETIC TEST ONLY; NOT A CHECKPOINT\n')
    receipt['case_checkpoint_sha256'] = sha(case/'initial_state.pkl.gz')
    (case/'target_fit.csv').write_bytes(Path(p['reference_fit']['path']).read_bytes())
    params = list(core.table(p['reference_parameters']['path'],'parameter').values())
    for row in params:
        if row['parameter'] in req['point']: row['estimate'] = req['point'][row['parameter']]
    csv_write(case/'parameters.csv',params)
    plots = case/'standard_diagnostics'; plots.mkdir()
    for name in c['standard_diagnostic_names']: (plots/name).write_bytes(b'SYNTHETIC PLOT PLACEHOLDER\n')
    extra_receipt(receipt,p,c,plan,req)
    write(case/'receipt.json',receipt)
    write(folder/'success.json',dict(context=req['context'],receipt_sha256=sha(case/'receipt.json'),loss=receipt['loss'],checkpoint_sha256=receipt['case_checkpoint_sha256']))
    if behavior == 'bad_receipt':
        receipt['diagnostic_identity'] = 'injected_integrity_mismatch'
        write(case/'receipt.json',receipt)
    return 0


def worker(a,p,c,objs,plan):
    req = read(a.request); folder = a.output
    expected = expected_request(p,c,plan,req['arm'],sha(a.contract),folder,req['deadline_epoch'],req['synthetic'])
    expected['controller_pid'] = os.getppid()
    require(req == expected, 'Request identity differs from fixed paired proposal')
    if req['synthetic']:
        require('smoke_v1' in folder.parts, 'Synthetic mode may only write the smoke subtree')
    ac = arm_contract(c,plan,req['arm'])
    folder.mkdir(parents=True,exist_ok=False)
    write(folder/'startup.json',dict(context=req['context'],pid=os.getpid(),parent_pid=os.getppid(),normalization_inputs=ac['normalization']))
    if req['synthetic']: return fixture(folder,p,ac,plan,req)
    evaluator = None; phase='setup'; started=time.monotonic()
    try:
        runtime=core.module('numerical_pair_native_runtime',c['files']['runtime']['path'])
        original_runtime_contract=copy.deepcopy(dict(c,objective=p['objective']))
        evaluator=runtime.setup(original_runtime_contract,objs['primary'],folder)
        # Apply only after the native factory has authenticated unchanged source,
        # original scientific settings, ancestry, objectives and parameter bounds.
        effective=apply_adapter(evaluator,original_runtime_contract,ac['normalization']['initial_psi'])
        require(effective == dict(ac,objective=p['objective']), 'Unexpected adapter difference')
        write(folder/'adapter_receipt.json',dict(original_contract_sha256=p['original_contract']['sha256'],
            pair_contract_sha256=sha(a.contract),adapter_sha256=p['files']['adapter.py']['sha256'],
            original_normalization=c['normalization'],effective_normalization=effective['normalization'],
            economic_specification_changes=[],changed_fields=['normalization.initial_psi'] if req['arm']==ARMS[1] else []))
        phase='objective'
        receipt=evaluator.evaluate(req['point'],folder/'case',req['deadline_epoch'],graphs=True,due=True,fixed_reference=False)
        extra_receipt(receipt,p,ac,plan,req)
        receipt['worker_elapsed_seconds']=time.monotonic()-started
        write(folder/'case/receipt.json',receipt)
        write(folder/'success.json',dict(context=req['context'],receipt_sha256=sha(folder/'case/receipt.json'),loss=receipt['loss'],checkpoint_sha256=receipt['case_checkpoint_sha256']))
        return 0
    except Exception as exc:
        failure=dict(status='fatal',classification='unknown_or_integrity_failure')
        if evaluator is not None:
            try: failure.update(evaluator.classify_failure(exc,req['context']))
            except Exception as err: failure.update(status='fatal',classifier_error=str(err))
        failure.update(context=req['context'],phase=phase,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
        audit=getattr(exc,'audit',getattr(exc,'ledger',None))
        if audit is not None: write(folder/'error_ledger.json',dict(context=req['context'],audit=audit))
        write(folder/'failure.json',failure)
        raise


def validate_extra(folder,p,ac,plan,req):
    r=read(folder/'case/receipt.json')
    for key,value in dict(diagnostic_identity=identity(p,ac,plan),original_contract_sha256=p['original_contract']['sha256'],pair_adapter_sha256=p['files']['adapter.py']['sha256'],pair_plan_sha256=p['plan']['sha256'],synthetic_fixture=bool(req['synthetic'])).items():
        require(r[key] == value, 'Diagnostic receipt mismatch: '+key)
    require(abs(r['normalization']['completed_fertility']-2.1)<=.0005, 'Completed fertility gate')
    require(abs(r['adult_entry_gate']['fertility_gap'])<=.0005, 'Demographic renewal gate')
    require(abs(r['market_residual'])<=2e-4, 'Market audit gate')
    return dict(stationary_solves=r['objective_stationary_solves'],stationary_solve_seconds=r['objective_stationary_solve_seconds'],price=r['price'],psi_child=r['normalization']['psi_child'])


def csv_write(path,rows):
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def dispatch(a,p,c,objs,plan,output,behaviors=None,cap=None):
    """Same finite two-process scheduler for the synthetic smoke and real pair."""
    output.mkdir(parents=True,exist_ok=False)
    helper=core.module('numerical_pair_supervisor',c['files']['recovery_search']['path'])
    started=time.time(); global_end=started+(30 if behaviors else plan['whole_experiment_cap_seconds'])
    cap=cap if cap is not None else plan['objective_cap_seconds']
    records=[]; active={}; fatal=False
    write(output/'clock.json',dict(start=started,end=global_end,case_cap_seconds=cap,workers=2))
    write(output/'latest_completed.json',dict(status='none_completed'))
    write(output/'best_so_far.json',dict(status='none_completed',promotion='disabled'))
    write(output/'arm_failures.json',[])
    def save():
        write(output/'records.json',records)
        write(output/'heartbeat.json',dict(epoch=time.time(),elapsed=time.time()-started,active=list(active),completed=len(records),fatal=fatal))
        write(output/'arm_failures.json',[r for r in records if r['status']!='success'])
        if records: write(output/'latest_completed.json',records[-1])
        success=[r for r in records if r['status']=='success']
        if success: write(output/'best_so_far.json',dict(min(success,key=lambda x:x['loss']),promotion='disabled'))
    try:
        for i,arm in enumerate(ARMS):
            deadline=min(global_end,time.time()+cap)
            folder=output/arm
            req=expected_request(p,c,plan,arm,sha(a.contract),folder,deadline,None if behaviors is None else behaviors[i])
            path=output/(arm+'.request.json');write(path,req)
            env=os.environ.copy();env.update({key:'1' for key in core.THREADS})
            proc=helper.ManagedProcess([sys.executable,__file__,'--mode','worker','--contract',str(a.contract),'--output',str(folder),'--request',str(path)],output/(arm+'.log'),deadline,env)
            active[arm]=(req,proc,folder)
        while active:
            for arm,(req,proc,folder) in list(active.items()):
                code=proc.poll()
                if code is None: continue
                ac=arm_contract(c,plan,arm)
                status,data,error=core.classify(folder,ac,objs,req,proc,code)
                try:
                    verify(a.contract)
                    if status=='success': data.update(validate_extra(folder,p,ac,plan,req))
                except Exception as exc: status,data,error='fatal',{},dict(error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
                row=dict(arm=arm,status=status,error=error,returncode=code,elapsed_seconds=time.time()-proc.started_epoch,request_path=str(output/(arm+'.request.json')),synthetic=bool(req['synthetic']),**data)
                ledger_path=folder/'case/stationary_solves.json'
                ledger=read(ledger_path) if ledger_path.exists() else []
                row.update(stationary_solves_started=len(ledger),stationary_solves_completed=sum(r['status']=='completed' for r in ledger),stationary_solve_seconds_recorded=sum(float(r.get('seconds',0)) for r in ledger))
                if len(ledger)>23: row.update(status='fatal',error='Normalization solve cap violated');status='fatal'
                records.append(row);proc.close();del active[arm]
                if status=='fatal': fatal=True
                save()
            if fatal:
                for arm,(req,proc,folder) in list(active.items()):
                    proc.cancel()
                    ledger_path=folder/'case/stationary_solves.json'
                    ledger=read(ledger_path) if ledger_path.exists() else []
                    records.append(dict(arm=arm,status='cancelled_after_fatal',returncode=proc.process.returncode,error='Sibling unknown/integrity failure; no retry',elapsed_seconds=time.time()-proc.started_epoch,stationary_solves_started=len(ledger),stationary_solves_completed=sum(r['status']=='completed' for r in ledger),stationary_solve_seconds_recorded=sum(float(r.get('seconds',0)) for r in ledger)))
                    del active[arm]
                save();break
            save()
            if active: time.sleep(.2 if behaviors else 5.)
    finally:
        for _,proc,_ in active.values(): proc.close()
    write(output/'dispatch_receipt.json',dict(records=records,elapsed_seconds=time.time()-started,model_solves_started=0 if behaviors else sum(r.get('stationary_solves_started',0) for r in records),model_solves_completed=0 if behaviors else sum(r.get('stationary_solves_completed',0) for r in records),scope='synthetic zero-solve loop' if behaviors else 'all arm ledgers, including failed/censored cases',no_retries=True,complete=len(records)==2))
    return records


def paired_screens(left,right,lf,rf,reference,plan):
    limits=plan['paired_consistency']
    details=[]
    for key,value,limit in [('loss',abs(left['loss']-right['loss']),limits['loss_absolute_difference']),
          ('price_relative',abs(left['price']-right['price'])/max(abs(left['price']),abs(right['price'])),limits['price_relative_difference']),
          ('psi',abs(left['normalization']['psi_child']-right['normalization']['psi_child']),limits['psi_absolute_difference'])]:
        details.append(dict(screen=key,value=value,limit=limit,passed=value<=limit))
    scored=[]
    for moment,row in lf.items():
        if row['role']=='scored': scored.append(math.sqrt(float(row['weight']))*abs(float(row['model'])-float(rf[moment]['model'])))
        elif row['role']=='validation':
            difference=abs(float(row['model'])-float(rf[moment]['model']))
            limit=.01*max(1.,abs(float(reference[moment]['model'])),abs(float(row['target'])))
            details.append(dict(screen='validation:'+moment,value=difference,limit=limit,passed=difference<=limit))
    details.append(dict(screen='scored_max_sqrt_weight_scaled_difference',value=max(scored),limit=limits['scored_max_abs_difference_in_working_scale'],passed=max(scored)<=limits['scored_max_abs_difference_in_working_scale']))
    return dict(passed=all(r['passed'] for r in details),screens=details,failed_action=limits['failure_action'])


def comparison(p,plan,output,records):
    if len(records)!=2 or any(r['status']!='success' for r in records):
        write(output/'paired_comparison.json',dict(status='incomplete_pair',numerical_equivalence_claim=False,records=records));return
    cases=[output/arm/'case' for arm in ARMS]
    receipts=[read(case/'receipt.json') for case in cases]
    fits=[core.table(case/'target_fit.csv','moment') for case in cases]
    ref=core.table(p['reference_fit']['path'],'moment')
    screens=paired_screens(*receipts,*fits,ref,plan)
    prediction=core.table(p['predictions']['path'],'moment')
    rows=[]
    for moment,row in fits[0].items():
        rows.append(dict(moment=moment,role=row['role'],target=float(row['target']),reference=float(ref[moment]['model']),predicted_full=float(prediction[moment]['predicted_full']),predicted_half=float(prediction[moment]['predicted_half']),retained_start=float(row['model']),predicted_start=float(fits[1][moment]['model']),weight=row['weight'],retained_loss=row['loss_contribution'],predicted_start_loss=fits[1][moment]['loss_contribution']))
    csv_write(output/'predictions_vs_actual.csv',rows)
    indexed={r['arm']:r for r in records}
    speed=dict(retained_seconds=indexed[ARMS[0]]['elapsed_seconds'],predicted_seconds=indexed[ARMS[1]]['elapsed_seconds'],retained_stationary_seconds=receipts[0]['objective_stationary_solve_seconds'],predicted_stationary_seconds=receipts[1]['objective_stationary_solve_seconds'],retained_solves=receipts[0]['objective_stationary_solves'],predicted_solves=receipts[1]['objective_stationary_solves'])
    write(output/'paired_comparison.json',dict(status='paired_screens_pass' if screens['passed'] else 'numerical_start_sensitivity',screens=screens,speed=speed,reference_loss=read(p['reference_receipt']['path'])['loss'],predicted_full_loss=plan['predicted_full_loss'],predicted_half_loss=plan['predicted_half_loss'],actual_losses=[r['loss'] for r in receipts],numerical_equivalence_claim=screens['passed'],promotion='disabled; no final exact repeats'))


def smoke(a,p,c,objs,plan):
    output=BASE/'smoke_v1';output.mkdir(exist_ok=False)
    # Full original-source authentication, but no runtime.setup or model solves.
    runtime=core.module('numerical_pair_smoke_native_runtime',c['files']['runtime']['path'])
    runtime.verify_sources(c)
    tests=[]
    for name,behaviors,wanted,cap in [('success_pair',['success','success'],['success','success'],10),('explicit_failure',['inadmissible','success'],['inadmissible','success'],10),('owned_timeout',['timeout','success'],['censored_timeout','success'],3),('fatal_stop',['fatal','timeout'],['fatal','cancelled_after_fatal'],10),('receipt_integrity',['bad_receipt','timeout'],['fatal','cancelled_after_fatal'],10)]:
        records=dispatch(a,p,c,objs,plan,output/name,behaviors=behaviors,cap=cap)
        got={r['arm']:r['status'] for r in records}
        require([got[arm] for arm in ARMS]==wanted,'Smoke dispatch failure: '+name+': '+str(got))
        tests.append(dict(test=name,passed=True,statuses=got))
    timeout_pid=read(output/'owned_timeout/retained_start/grandchild.json')['pid']
    stat=Path('/proc')/str(timeout_pid)/'stat'
    require(not stat.exists() or stat.read_text().split()[2]=='Z','Timeout left live descendant')
    tests.append(dict(test='process_group_descendant_killed',passed=True))
    success=output/'success_pair';records=read(success/'records.json')
    comparison(p,plan,success,records)
    require(read(success/'paired_comparison.json')['screens']['passed'],'Identical pair rejected')
    cases=[success/arm/'case' for arm in ARMS]
    receipts=[read(case/'receipt.json') for case in cases]
    fits=[core.table(case/'target_fit.csv','moment') for case in cases]
    ref=core.table(p['reference_fit']['path'],'moment')
    for field in ('loss','price','psi','scored','validation'):
        rr=copy.deepcopy(receipts);ff=copy.deepcopy(fits)
        if field=='loss':rr[1]['loss']+=.051
        elif field=='price':rr[1]['price']*=1.001
        elif field=='psi':rr[1]['normalization']['psi_child']+=.00011
        else:
            role='scored' if field=='scored' else 'validation'
            key=next(k for k,r in ff[1].items() if r['role']==role)
            ff[1][key]['model']=float(ff[1][key]['model'])+1.
        require(not paired_screens(*rr,*ff,ref,plan)['passed'],'Sensitivity screen failed to reject '+field)
        tests.append(dict(test='paired_screen_'+field,passed=True))
    altered=copy.deepcopy(p);altered['target_weight_fingerprint']='changed'
    try: verify_fingerprint(altered,c,objs['primary'])
    except RuntimeError: pass
    else: raise RuntimeError('Bad target fingerprint was accepted')
    tests.append(dict(test='fingerprint_failfast',passed=True))
    from types import SimpleNamespace
    original_view=copy.deepcopy(dict(c,objective=p['objective']))
    fake=SimpleNamespace(c=original_view)
    adapted=apply_adapter(fake,original_view,plan['predicted_initial_psi_full'])
    require(adapted==dict(arm_contract(c,plan,ARMS[1]),objective=p['objective']),'Adapter result differs')
    require(original_view==dict(c,objective=p['objective']),'Adapter mutated original contract')
    tests.append(dict(test='adapter_only_copied_initial_guess',passed=True))
    write(output/'smoke_receipt.json',dict(status='synthetic_exact_loop_passed',contract_sha256=sha(a.contract),model_solves=0,tests=tests,original_source_manifest_sha256=c['source_manifest']['sha256'],no_main_run=True,slurm_job=os.environ['SLURM_JOB_ID']))
    print(json.dumps(dict(status='synthetic_exact_loop_passed',tests=len(tests),model_solves=0)))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['smoke','run','worker'],required=True)
    parser.add_argument('--contract',type=Path,required=True);parser.add_argument('--output',type=Path)
    parser.add_argument('--request',type=Path);parser.add_argument('--approval',type=Path)
    a=parser.parse_args();a.contract=a.contract.resolve()
    p,c,objs,plan=verify(a.contract)
    if a.mode=='worker':return worker(a,p,c,objs,plan)
    if a.mode=='smoke':smoke(a,p,c,objs,plan);return 0
    require(a.approval is not None,'Reviewed main-run approval required')
    require(sha(a.approval)==os.environ.get('EXPECTED_NUMERICAL_PAIR_APPROVAL_SHA256'),'Approval hash differs')
    approved=read(a.approval)
    require(approved['status']=='approved_bounded_numerical_pair' and approved['pair_contract_sha256']==sha(a.contract),'Wrong launch approval')
    pin=approved['smoke_receipt'];require(sha(pin['path'])==pin['sha256'],'Smoke receipt changed')
    receipt=read(pin['path']);require(receipt['status']=='synthetic_exact_loop_passed' and receipt['contract_sha256']==sha(a.contract) and receipt['model_solves']==0,'Smoke is not valid for this contract')
    require(time.time()<approved['launch_not_after_epoch'],'Launch authorization expired; no late restart')
    output=BASE/'run_v1'
    records=dispatch(a,p,c,objs,plan,output)
    comparison(p,plan,output,records)
    return 1 if any(r['status']=='fatal' for r in records) else 0


if __name__=='__main__':sys.exit(main())
