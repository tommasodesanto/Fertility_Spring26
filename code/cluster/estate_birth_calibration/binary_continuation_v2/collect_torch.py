"""Read-only validation/collection for all 20 continuation task receipts."""
import argparse,hashlib,json
from pathlib import Path
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def collect(stage):
    inv=read(stage/'inventory.json');plan=read(stage/'control/starts.json');starts_sha=sha(stage/'control/starts.json')
    assert len(plan['starts'])==20 and plan['continuation_arm']=='binary'
    rows=[]
    for task in range(20):
        folder=stage/'results'/f'production_binary_chain_{task}';run=folder/'run'
        if not (run/'completed.json').exists():
            rows.append(dict(task=task,status='pending_or_incomplete',heartbeat=read(run/'heartbeat.json') if (run/'heartbeat.json').exists() else None,
              failure=read(run/'failure.json') if (run/'failure.json').exists() else None));continue
        launch=read(folder/'launcher_start.json');terminal=read(folder/'launcher_terminal.json');contract=read(run/'start_contract.json');done=read(run/'completed.json')
        assert launch['mode']=='production' and launch['arm']=='binary' and launch['chain']==task and launch['cpus']==1 and launch['memory_GiB']==24 and launch['wall_seconds']==43200
        assert launch['stage_inventory_sha256']==sha(stage/'inventory.json') and launch['starts_sha256']==starts_sha
        assert terminal['slurm_job_id']==launch['slurm_job_id'] and terminal['slurm_array_job_id']==launch['slurm_array_job_id'] and terminal['exit_code']==0
        assert contract['arm']=='binary' and contract['birth_cap']==1 and contract['chain']==task and contract['seed']==plan['starts'][task] and contract['starts_count']==20
        assert contract['starts_file_sha256']==starts_sha and contract['target_fingerprint']==inv['target_fingerprint'] and contract['weight_fingerprint']==inv['weight_fingerprint']
        assert done['arm']=='binary' and done['birth_cap']==1 and done['chain']==task and done['starts_file_sha256']==starts_sha
        if done['status']=='selected_numerically_verified':
            assert done['selected_postcheck']['status']=='passed' and done['repeat']['status']=='exact_full_ge_repeat_passed'
            assert done['repeat']['experimental_target_fit_exact'] and len(done['target_fit'])==14 and len(done['parameters'])==31 and len(done['repeat']['standard_plot_hashes'])==17
            child=read(run/'native_postcheck/completed.json');search=read(run/'search_completed.json')
            assert child['status']=='full_native_postcheck_passed' and child['search_receipt_sha256']==sha(run/'search_completed.json')
            assert child['starts_file_sha256']==starts_sha==search['starts_file_sha256']
        rows.append(dict(task=task,**done))
    best=min((r for r in rows if r.get('status')=='selected_numerically_verified'),key=lambda r:r['native_loss'],default=None)
    return dict(status='read_only_20_chain_collection',job_id='from_launch_receipts',inventory_sha256=sha(stage/'inventory.json'),starts_sha256=starts_sha,
      target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'],chains=rows,best_verified=best,no_adoption=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not a.out.exists(),'Refusing existing collection';result=collect(a.stage);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(chains=len(result['chains']),verified=sum(r.get('status')=='selected_numerically_verified' for r in result['chains']))))
