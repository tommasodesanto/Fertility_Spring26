"""Record or verify the required exact-loop full-native smoke gate."""
import argparse,hashlib,json
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(receipt_path):
    r=read(receipt_path);assert r['status']=='passed_exact_two_call_native_smoke'
    submitted=read(REMOTE/'control/smoke_submission.json')
    assert r['task']==0 and r['objective_calls']==2 and r['native_selected_postcheck']=='passed' and r['exact_repeat']=='passed'
    assert submitted['status']=='submitted' and submitted['job_id']==r['job_id'] and submitted['chain']==r['task']
    assert r['stage_inventory_sha256']==sha(REMOTE/'inventory.json') and r['starts_sha256']==sha(REMOTE/'control/starts.json')
    p=REMOTE/'results/smoke_binary_chain_0'
    assert sha(p/'launcher_terminal.json')==r['launcher_terminal_sha256']
    assert sha(p/'run/completed.json')==r['completed_sha256']
    d=read(p/'run/completed.json')
    assert d['status']=='selected_numerically_verified' and d['objective_calls']==2
    assert d['selected_postcheck']['status']=='passed' and d['repeat']['status']=='exact_full_ge_repeat_passed' and len(d['repeat']['standard_plot_hashes'])==17
    return {'status':'verified','receipt_sha256':sha(receipt_path)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify-existing',type=Path);ap.add_argument('--smoke-dir',type=Path);ap.add_argument('--job-id');a=ap.parse_args()
    if a.verify_existing: print(json.dumps(verify(a.verify_existing)));return
    if not a.smoke_dir or not a.job_id:ap.error('provide --smoke-dir and --job-id after the smoke finishes')
    d=read(a.smoke_dir/'run/completed.json');launch=read(a.smoke_dir/'launcher_start.json');terminal=read(a.smoke_dir/'launcher_terminal.json')
    assert launch['mode']=='smoke' and launch['chain']==0 and launch['slurm_array_job_id']==a.job_id and terminal['slurm_array_job_id']==a.job_id and terminal['exit_code']==0
    submitted=read(REMOTE/'control/smoke_submission.json');assert submitted['job_id']==a.job_id and submitted['chain']==0
    assert d['status']=='selected_numerically_verified' and d['objective_calls']==2 and d['birth_cap']==1
    assert d['selected_postcheck']['status']=='passed' and d['repeat']['status']=='exact_full_ge_repeat_passed' and len(d['repeat']['standard_plot_hashes'])==17 and len(d['target_fit'])==14 and len(d['parameters'])==31
    inv=read(REMOTE/'inventory.json');assert d['target_fingerprint']==inv['target_fingerprint'] and d['weight_fingerprint']==inv['weight_fingerprint']
    r=dict(status='passed_exact_two_call_native_smoke',job_id=a.job_id,task=0,objective_calls=2,native_selected_postcheck='passed',exact_repeat='passed',
      stage_inventory_sha256=sha(REMOTE/'inventory.json'),starts_sha256=sha(REMOTE/'control/starts.json'),launcher_terminal_sha256=sha(a.smoke_dir/'launcher_terminal.json'),completed_sha256=sha(a.smoke_dir/'run/completed.json'),
      target_fingerprint=d['target_fingerprint'],weight_fingerprint=d['weight_fingerprint'],no_production_submission=True)
    path=REMOTE/'control/smoke_gate.json';assert not path.exists(),'Refusing existing smoke gate receipt';path.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(verify(path)))
if __name__=='__main__':main()
