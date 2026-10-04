"""Read-only collector for the 5 binary and 10 count-three recovery tasks."""
import argparse,hashlib,json
from pathlib import Path
REL='output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def collect(stage):
    inv=read(stage/'inventory.json');submission=read(stage/'control/production_submission.json');starts=read(stage/'control/starts_receipt.json')
    assert submission['array']=='0-14%15' and submission['cores_max']==15 and submission['cutoff_epoch']==1791122400
    assert submission['stage_inventory_sha256']==sha(stage/'inventory.json') and submission['starts_receipt_sha256']==sha(stage/'control/starts_receipt.json')
    assert starts['status']=='provisional_parent_saved_best_only' and starts['parent_inventory_sha256']==inv['parent_inventory_sha256']
    plans={arm:stage/'source'/REL/f'plan_{arm}.json' for arm in ('binary','count3')}
    for arm in plans:
        assert sha(plans[arm])==starts[f'plan_{arm}_sha256']
    rows=[]
    for task in range(15):
        arm='binary' if task<5 else 'count3';chain=task if task<5 else task-5;cap=1 if arm=='binary' else 3
        plan=read(plans[arm]);plan_sha=sha(plans[arm]);folder=stage/'results'/f'production_{arm}_chain_{chain}';run=folder/'run'
        if not (run/'completed.json').exists():
            rows.append(dict(task=task,arm=arm,chain=chain,status='pending_or_incomplete',
                             heartbeat=read(run/'heartbeat.json') if (run/'heartbeat.json').exists() else None,
                             failure=read(run/'failure.json') if (run/'failure.json').exists() else None,
                             launcher_terminal=read(folder/'launcher_terminal.json') if (folder/'launcher_terminal.json').exists() else None));continue
        contract=read(run/'start_contract.json');done=read(run/'completed.json');launcher=read(folder/'launcher_start.json')
        assert launcher['stage_inventory_sha256']==sha(stage/'inventory.json') and launcher['mode']=='production' and launcher['arm']==arm and launcher['chain']==chain
        assert launcher['cpus']==1 and launcher['memory_GiB']==24 and launcher['deadline_epoch']<=1791122400
        assert contract['arm']==arm and contract['chain']==chain and contract['birth_cap']==cap
        assert contract['starts_file_sha256']==plan_sha and contract['seed']==plan['starts'][chain] and contract['all_starts']==plan['starts']
        assert contract['bounds']==plan['bounds'] and contract['selected_source_sha256']==inv['selected_source_sha256']
        assert contract['provisional_seed_source_sha256']==inv['parent_inventory_sha256']
        assert all(contract[k]==done[k]==inv[k] for k in ('target_fingerprint','weight_fingerprint'))
        assert done['arm']==arm and done['chain']==chain and done['starts_file_sha256']==plan_sha
        if done['status']=='selected_numerically_verified':
            assert done['selected_postcheck']['status']=='passed' and done['selected_postcheck']['experiment_flags']['birth_count_choice_cap']==cap
            assert len(done['target_fit'])==14 and len(done['parameters'])==31 and len(done['repeat']['standard_plot_hashes'])==17
            assert done['repeat']['status']=='exact_full_ge_repeat_passed' and done['repeat']['experimental_target_fit_exact']
            child=read(run/'native_postcheck/completed.json');search=read(run/'search_completed.json')
            assert child['status']=='full_native_postcheck_passed' and child['search_receipt_sha256']==sha(run/'search_completed.json')
            assert child['starts_file_sha256']==search['starts_file_sha256']==plan_sha
        rows.append(dict(task=task,**done))
    best={arm:min((r for r in rows if r['arm']==arm and r['status']=='selected_numerically_verified'),key=lambda r:r['native_loss'],default=None) for arm in plans}
    return dict(status='read_only_recovery_collection',job_id=submission['job_id'],inventory_sha256=sha(stage/'inventory.json'),
                target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'],chains=rows,best_by_arm=best,
                no_adoption=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();assert not args.out.exists(),'Refusing existing collection';result=collect(args.stage)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(chains=len(result['chains']),verified=sum(r['status']=='selected_numerically_verified' for r in result['chains']))))
