#!/usr/bin/env python3
"""Torch-side continuation after the bounded initial search has completed.

Schedule this utility afterok:17596347. It reads local receipts, preserves the
original six histories, and prepares six separate refits only after a strictly
improved selected initial state passes its two exact repetitions. Submission is
opt-in. No model solve or scheduler query occurs in this utility.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(value,indent=2,allow_nan=False)+'\n'
    # Replaying this continuation must not overwrite a different selected state
    # or launched manifest. Identical published artifacts may be reused.
    if path.exists():
        if read(path)==value:return
        mutable=path.name in ('history_manifest_refit.json','history_array_manifest_refit.json',
            'continuation_after_search_receipt.json','continuation_after_search_failure.json')
        is_manifest=path.name in ('history_manifest_refit.json','history_array_manifest_refit.json')
        if not mutable or (is_manifest and (path.parent/'histories_refit/submission.json').exists()):
            raise ValueError('Refusing to replace existing continuation artifact: '+str(path))
        path.write_text(text);return
    with path.open('x') as stream:stream.write(text)


def verify(pins):
    if not isinstance(pins,dict) or not pins:raise ValueError('Nonempty immutable source pins required')
    for path,digest in pins.items():
        if sha(path)!=digest:raise ValueError('Source/checkpoint hash mismatch: '+str(path))


def option(command,name):
    if command.count(name)!=1:raise ValueError('Exactly one command option required: '+name)
    index=command.index(name)
    if index+1>=len(command):raise ValueError('Missing option value: '+name)
    return index+1,command[index+1]


def checkpoint_item(item):
    path=item.get('checkpoint',item.get('path'))
    digest=item.get('checkpoint_sha256',item.get('sha256'))
    if not path or not digest:raise ValueError('Selected checkpoint path and SHA256 required')
    verify({path:digest})
    return path,digest


def prepare(batch,deadline_unix,now=None):
    """Validate receipts and publish a separate immutable six-track manifest."""
    batch=Path(batch).resolve();search=batch/'initial_search_joint'
    summary=read(search/'summary.json');contract=read(search/'search_contract.json')
    exact=read(search/'cases/selected_exact_repetitions/result.json')
    smoke=read(batch/'joint_initial_smoke/summary.json')
    if summary.get('status')!='completed_rebated_initial_search':
        raise ValueError('Initial search prerequisite is not completed')
    if (summary.get('selected_exact_repetitions_verified') is not True
            or contract.get('final_exact_repetitions')!=2
            or exact.get('status')!='verified' or exact.get('second_signature_equal') is not True
            or exact.get('proposal',{}).get('repetitions')!=2
            or len(exact.get('accounting',[]))!=2):
        raise ValueError('Two exact selected-point repetitions are not verified')
    loss=float(summary['selected_loss']);smoke_loss=float(smoke['loss'])
    if not math.isfinite(loss) or not math.isfinite(smoke_loss) or loss!=float(exact['loss']):
        raise ValueError('Selected loss is nonfinite or differs from exact repetitions')
    if smoke.get('status')!='verified_rebated_initial_smoke':raise ValueError('Original smoke is not verified')
    source_root=str(Path(contract['source_root']).resolve())
    if str(Path(smoke['source_root']).resolve())!=source_root:
        raise ValueError('Search and original smoke source roots differ')
    if Path(contract['smoke']).resolve() not in ((batch/'joint_initial_smoke').resolve(),(batch/'joint_initial_smoke/summary.json').resolve()):
        raise ValueError('Search did not use this batch original smoke')
    selected_path,selected_sha=checkpoint_item(summary['checkpoint'])
    if summary['checkpoint']!=exact['accounting'][-1]:
        raise ValueError('Published selected checkpoint differs from second exact repetition')
    for receipt in exact['accounting']:checkpoint_item(receipt)
    if len({r.get('repetition') for r in exact['accounting']})!=2:
        raise ValueError('Exact repetitions must have distinct recorded identities')
    result=dict(selected_loss=loss,original_smoke_loss=smoke_loss,source_root=source_root,
        search_job='17596347',deadline_unix=deadline_unix,search_summary_sha256=sha(search/'summary.json'),
        selected_checkpoint_sha256=selected_sha)
    if not loss<smoke_loss:
        return dict(result,status='skipped_no_strict_initial_improvement')
    remaining=math.floor(float(deadline_unix)-(time.time() if now is None else float(now)))
    if not math.isfinite(float(deadline_unix)) or remaining>43200:
        raise ValueError('Fixed shared deadline must be within the original twelve-hour envelope')
    driver_seconds=remaining-60
    if driver_seconds<=3*3600:
        return dict(result,status='skipped_insufficient_shared_time',remaining_seconds=max(0,remaining))
    history_path=batch/'history_manifest_joint.json';array_path=batch/'history_array_manifest_joint.json'
    history=read(history_path);array=read(array_path)
    verify(history['file_sha256']);verify(array['source_pins'])
    if Path(history['initial_summary']).resolve()!=(batch/'joint_initial_smoke/summary.json').resolve():
        raise ValueError('Original history manifest does not use the valid smoke')
    if len(array['stages'])!=6:raise ValueError('Exactly six original history tracks required')
    selected_summary=batch/'selected_initial_summary.json'
    selected=dict(status='verified_rebated_initial_smoke',loss=loss,
        checkpoint=summary['checkpoint'],source_root=source_root,
        selected_exact_repetitions_verified=True,selection_source=str(search/'summary.json'),
        selection_source_sha256=sha(search/'summary.json'))
    write(selected_summary,selected)
    history_new=copy.deepcopy(history);history_new['initial_summary']=str(selected_summary)
    # Explicit old raw-summary overrides must follow the newly selected packet.
    if 'initial_raw_summary' in history_new:
        history_new['initial_raw_summary']=str(Path(selected_path).parent/'summary.json')
    history_new['policy_reserve_seconds']=min(2*3600,int(driver_seconds*.25))
    history_new['file_sha256'][str(selected_summary)]=sha(selected_summary)
    history_new_path=batch/'history_manifest_refit.json';write(history_new_path,history_new)
    array_new=copy.deepcopy(array);array_new['runroot']=str(batch/'histories_refit')
    array_new['max_workers']=6;array_new['horizon_hours']=remaining/3600.
    array_new['source_pins'][str(selected_summary)]=sha(selected_summary)
    array_new['source_pins'][str(history_new_path)]=sha(history_new_path)
    matrix=set()
    for stage in array_new['stages']:
        if stage.get('cpus',1)!=1:raise ValueError('One numerical CPU per history track required')
        command=stage['command'];case=option(command,'--case')[1];count=int(option(command,'--count')[1])
        matrix.add((case,count))
        for flag,value in (('--manifest',str(history_new_path)),
                ('--output',str(batch/'histories_refit'/stage['name'])),('--seconds',str(driver_seconds))):
            index,_=option(command,flag);command[index]=value
        # Queue time counts against the shared deadline. Recompute the driver
        # budget when the array element starts, rather than giving a late start
        # the original submission-time allowance.
        launch_code=("import json,math,pathlib,subprocess,sys,time; "
            "deadline=float(sys.argv[1]); argv=sys.argv[2:]; "
            "seconds=min(int(argv[argv.index('--seconds')+1]),math.floor(deadline-time.time())-60); "
            "out=pathlib.Path(argv[argv.index('--output')+1]); out.mkdir(parents=True,exist_ok=True); "
            "(out/'queue_deadline_receipt.json').write_text(json.dumps({'remaining_driver_seconds':seconds,'deadline_unix':deadline})); "
            "argv[argv.index('--seconds')+1]=str(seconds); "
            "sys.exit(subprocess.call(argv) if seconds>10800 else 0)")
        stage['command']=[sys.executable,'-c',launch_code,str(deadline_unix),*command]
        # The submitter uses stage seconds for Slurm, leaving sixty seconds
        # after the driver's stop while remaining inside the shared deadline.
        stage['seconds']=driver_seconds+60
    if matrix!={(case,count) for case in ('A0','A+') for count in (6,24,100)}:
        raise ValueError('Refit matrix must contain A0/A+ at six, twenty-four and one hundred dates')
    array_new_path=batch/'history_array_manifest_refit.json';write(array_new_path,array_new)
    return dict(result,status='prepared',array_manifest=str(array_new_path),
        array_manifest_sha256=sha(array_new_path),history_manifest=str(history_new_path),
        selected_initial_summary=str(selected_summary),driver_seconds=driver_seconds,
        scheduler_seconds=driver_seconds+60,existing_history_workers=6,
        pilot_workers=2,new_history_workers=6,total_numerical_workers=14,
        maximum_numerical_workers=18,policy_reserve_seconds=history_new['policy_reserve_seconds'])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--deadline-unix',type=float,required=True)
    parser.add_argument('--submit',action='store_true')
    args=parser.parse_args(argv);batch=args.batch.resolve()
    receipt_path=batch/'continuation_after_search_receipt.json'
    existing_submission=batch/'histories_refit/submission.json'
    if existing_submission.exists():
        existing=read(existing_submission)
        print(json.dumps(dict(status='already_submitted',submission=existing),indent=2));return
    if receipt_path.exists():
        previous=read(receipt_path)
        if previous.get('status')=='submitted':
            print(json.dumps(previous,indent=2));return
    try:
        result=prepare(batch,args.deadline_unix)
        if result['status']=='prepared':
            command=[sys.executable,str(batch/'submit_e5f_final_night.py'),result['array_manifest']]
            if args.submit:command.append('--submit')
            # The existing Torch-side submitter validates source pins, emits
            # the local array script and only calls sbatch with --submit.
            completed=subprocess.run(command,check=True,capture_output=True,text=True)
            submission=json.loads(completed.stdout)
            result.update(status='submitted' if args.submit else 'dry_run',
                submitter_command=command,submission=submission)
            if args.submit:
                job=str(submission.get('job',''))
                if not job.isdigit():raise RuntimeError('Submitter returned no numeric job ID; do not retry automatically')
                result['job']=job
        write(receipt_path,result)
        print(json.dumps(result,indent=2))
    except Exception as exc:
        failure=dict(status='failed',error_type=type(exc).__name__,error=str(exc))
        write(batch/'continuation_after_search_failure.json',failure)
        raise


if __name__=='__main__':main()
