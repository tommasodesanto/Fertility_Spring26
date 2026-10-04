"""Fail closed unless each arm passed exact two-call loop and fresh native repeat."""
import hashlib,json,sys
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_continuation_20261004_v1')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def verify(stage):
    inv=read(stage/'inventory.json');fingerprint=sha(stage/'inventory.json')
    for arm,cap in [('binary',1),('count3',3)]:
        launcher=stage/'results'/f'smoke_{arm}_chain_0';run=launcher/'run'
        started=read(launcher/'launcher_start.json');contract=read(run/'start_contract.json');done=read(run/'completed.json')
        child=read(run/'native_postcheck/completed.json');search=read(run/'search_completed.json')
        assert started['stage_inventory_sha256']==fingerprint and started['wall_seconds']==5400
        assert started['deadline_epoch']-started['start_epoch']==5400
        assert contract['birth_cap']==cap and contract['starts_count']==5 and contract['bounds']['beta_annual']==[.93,.99]
        assert contract['starts_file_sha256']==sha(stage/'source/output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1'/f'plan_{arm}.json')
        assert all(contract[k]==done[k]==child[k]==inv[k] for k in ('target_fingerprint','weight_fingerprint'))
        assert done['status']=='selected_numerically_verified' and done['objective_calls']==2
        assert len(read(run/'cases.json'))==2 and read(run/'heartbeat.json')['status']=='completed'
        assert child['status']=='full_native_postcheck_passed' and child['search_receipt_sha256']==sha(run/'search_completed.json')
        assert len(done['target_fit'])==14 and len(done['parameters'])==31
        assert len(done['repeat']['standard_plot_hashes'])==17 and done['repeat']['experimental_target_fit_exact']
        assert done['selected_postcheck']['status']=='passed'
        assert done['smoke_fast_full_comparison']['status']=='search_full_new_target_exact'
        assert search['search_evaluator']=='full_native' and search['selected_evaluator']=='full_native'
    return dict(status='matched_both_arm_smoke_passed',stage_inventory_sha256=fingerprint,target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'])
if __name__=='__main__':print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv)==2 else REMOTE)))
