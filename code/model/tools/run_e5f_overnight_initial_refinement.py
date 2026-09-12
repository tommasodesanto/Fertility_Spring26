"""Up to three unchanged capped-beta search batches; no agent/API calls."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, time

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(d,indent=2)+'\n');q.replace(p)
def best_single(folder):
    records=read(Path(folder)/'cases.json')
    return min((r for r in records if r['status']=='verified' and r['proposal'].get('repetitions',1)==1),key=lambda r:r['loss'])
def prepare(template, destination, seed):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=False)
    original=read(template/'plan_capped_beta_099.json')
    for name in ['run_capped_beta.py','run_profile.py','test_run_capped_beta.py','submit_capped_beta.sh']:
        shutil.copy2(template/name,destination/name)
    p=Path(seed['output']);case=p.parent;score=p/'scored_repetition_01/score.json'
    files=[p/'preflight.json',score,p/'summary.json',case/'initial_contract.json',case/'run_contract.json']
    manifest=dict(file_sha256={str(f):sha(f) for f in files},prior_repetitions=1,required_new_exact_seed_repetitions=2,seed_case=seed['case_id'])
    write(destination/'seed_manifest.json',manifest)
    original.update(output_dir=str(destination/'results'),resume_score_path=str(score),resume_score_sha256=sha(score),
        resume_summary_sha256=sha(p/'summary.json'),resume_proposal=seed['proposal'],
        seed_manifest_path=str(destination/'seed_manifest.json'),seed_manifest_sha256=sha(destination/'seed_manifest.json'))
    write(destination/'plan_capped_beta_099.json',original)
    return original

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--template',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--preflight-only',action='store_true');a=ap.parse_args()
    start=time.monotonic();a.output.mkdir(parents=True,exist_ok=True)
    seed=best_single(a.template/'results');initial_loss=seed['loss'];completed=[]
    for number in range(3):
        if time.monotonic()-start>9*3600-10860:break
        folder=a.output/f'round_{number}';plan=prepare(a.template,folder,seed)
        # Original full-loop controller tests and all seed/source pins are retained.
        cmd=[sys.executable,'-B','-m','unittest','discover','-s',str(folder),'-p','test_run_capped_beta.py','-q']
        subprocess.run(cmd,check=True,timeout=180)
        sys.path.insert(0,str(folder));import run_capped_beta as controller
        controller.validate_plan(plan);controller.validate_seed_manifest(plan)
        if a.preflight_only:
            write(a.output/'preflight.json',dict(status='passed',initial_loss=initial_loss,maximum_batches=3,maximum_new_search_cases=270,maximum_repetitions=282,maximum_GE=2256));return
        write(a.output/'latest_stage.json',dict(round=number,status='running',prior_best=seed['loss']))
        with (folder/'controller.log').open('w') as log:
            result=subprocess.run([sys.executable,'-B',str(folder/'run_capped_beta.py'),'--plan',str(folder/'plan_capped_beta_099.json')],stdout=log,stderr=subprocess.STDOUT,timeout=10860)
        summary=read(folder/'results/summary.json') if (folder/'results/summary.json').exists() else {'status':'missing_summary','exit':result.returncode}
        completed.append(summary);write(a.output/'latest_completed.json',dict(round=number,summary=summary))
        if summary.get('best_loss') is not None:write(a.output/'best_so_far.json',dict(loss=summary['best_loss'],folder=str(folder/'results')))
        if result.returncode or not summary.get('selected_exact_repetitions_verified'):break
        next_seed=best_single(folder/'results')
        if next_seed['loss']>=seed['loss']-1e-8:break
        seed=next_seed
    write(a.output/'summary.json',dict(status='finished',initial_loss=initial_loss,batches=completed,elapsed_seconds=time.monotonic()-start,
        targets_changed=False,rooms_target=.7202462623815278,beta_cap=.99,production_eligible=False))
if __name__=='__main__':main()
