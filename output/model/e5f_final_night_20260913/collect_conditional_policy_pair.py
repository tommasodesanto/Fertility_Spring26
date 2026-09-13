"""Assemble an explicitly conditional pair and verify its native input states."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--verifier',type=Path,required=True);ap.add_argument('--deadline',type=float,required=True);args=ap.parse_args()
    read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    run=args.run;spec=read(run/'spec.json');out=run/'conditional_pair'
    out.mkdir(exist_ok=True)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    if time.time()>args.deadline:save('collection_status.json',dict(status='author_cutoff'));return
    selected={};available={}
    for policy in ['baseline_rebate','tax2_rebate']:
        available[policy]=[]
        for seed in ['long','flat']:
            folder=run/(policy+'_'+seed);summary=folder/'summary.json'
            if not summary.exists():continue
            s=read(summary)
            if s.get('status')!='passed':continue
            assert s['spec_sha256']==sha(run/'spec.json') and s['source_history_count']==6 and s['forecast_count']==24 and s['source_history_refitted'] is False
            check=s['results'][policy];assert check['finite_converged'] and check['same_initial_g_pre']
            assert check['root_sha256']==sha(folder/policy/'root_receipt.json')
            available[policy].append(str(folder));selected.setdefault(policy,folder/policy)
    if len(selected)!=2:save('collection_status.json',dict(status='incomplete_pair',available=available));return
    assert not (out/'contract_receipt.json').exists(),'Immutable pair already assembled'
    source=Path(spec['source_case']);base=read(source/'contract_receipt.json')
    contract=dict(case='A0',count=24,conditional_history_count=6,history_refitted=False,
        manifest_sha256=sha(spec['runtime_manifest']),source_history=str(source),source_history_contract_sha256=sha(source/'contract_receipt.json'),
        initial_checkpoint_sha256=base['initial_checkpoint_sha256'],horizon_verified=False,production_eligible=False,
        selected={k:str(v) for k,v in selected.items()},available=available,selection='First available in fixed long-then-flat order; seed is numerical only')
    save('contract_receipt.json',contract)
    shutil.copyfile(source/'realized_fit.json',out/'realized_fit.json')
    shutil.copyfile(source/'finite_history_complete.json',out/'conditioning_history_complete.json')
    (out/'realized_state_2023.pkl.gz').symlink_to(source/'realized_state_2023.pkl.gz')
    (out/'policies').mkdir(exist_ok=True)
    for policy,folder in selected.items():(out/'policies'/policy).symlink_to(folder)
    cmd=['/share/apps/anaconda3/2025.06/bin/python','-B',str(args.verifier),'--case-dir',str(out),
        '--manifest',spec['runtime_manifest'],'--helper',spec['helper'],'--out',str(out/'policy_state_verification.json')]
    result=subprocess.run(cmd,timeout=max(1,min(300,args.deadline-time.time())),capture_output=True,text=True)
    save('collection_status.json',dict(status='passed' if result.returncode==0 else 'native_verification_failed',
        conditional_history_count=6,forecast_count=24,selected=contract['selected'],available=available,
        verifier_sha256=sha(args.verifier),returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    if result.returncode:raise RuntimeError('Native conditional pair verification failed')

if __name__=='__main__':main()
