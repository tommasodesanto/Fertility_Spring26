"""Isolated, bounded retry after the shared-diagnostics mutation was identified."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--frontier',action='store_true');args=parser.parse_args()
    base=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    old=base/'batches/night_surprises_20260912';new=base/('batches/night_surprise_frontier_20260912' if args.frontier else 'batches/night_surprises_recovery_20260912')
    plan=json.loads((old/'plan.json').read_text())
    plan.update(output_root=str(new/'results'),total_seconds=25200,seed_steps=[-.03,-.045],
        empirical_blocks=str(new/'empirical_blocks.csv'),
        recovery_reason='Deep-copy native terminal diagnostics before adding pension tail checks; avoid mutating shared default tolerances. Original failures preserved.',
        original_plan_sha256=sha(old/'plan.json'))
    if args.frontier:
        plan.update(native_smoke_only=True,total_seconds=5400,policy_reserve_seconds=0,maximum_trials_per_window=1,
            seed_steps=[-.0225,-.0275],diagnostic_question='Can smaller permanent shocks approach the first empirical fertility window while retaining a finite endpoint? Native six-date diagnostics only; no accepted history or policy.')
    oldpins=plan['file_sha256'];plan['file_sha256']={p:h for p,h in oldpins.items() if not Path(p).is_relative_to(old)}
    plan['file_sha256'].update({str(p):sha(p) for p in [*new.glob('*.py'),new/'empirical_blocks.csv']})
    (new/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    env=dict(os.environ,PYTHONPATH=f'{new}:{base}/code/model/tools:{base}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','test_e5f_successive_surprises','-q'],env=env,check=True,timeout=180)
    launch=(old/'submit_e5f_surprise_night.sh').read_text().replace('night_surprises_20260912','night_surprises_recovery_20260912').replace('--array=0-2','--array=0-1').replace('--time=10:00:00','--time=07:10:00').replace('e5f_surprise_fit','e5f_surprise_retry')
    if args.frontier:launch=launch.replace('night_surprises_recovery_20260912','night_surprise_frontier_20260912').replace('--time=07:10:00','--time=01:40:00').replace('--cpus-per-task=3','--cpus-per-task=1').replace('--mem=48G','--mem=16G').replace('e5f_surprise_retry','e5f_native_frontier')
    script=new/'submit_recovery.sh';script.write_text(launch)
    job=subprocess.check_output(['sbatch','--parsable',str(script)],text=True).strip().split(';')[0]
    collector=subprocess.check_output(['sbatch','--parsable','--account=torch_pr_570_general','--job-name=e5f_retry_report',f'--dependency=afterany:{job}:17440306','--cpus-per-task=1','--mem=4G','--time=00:10:00',f'--output={new}/report_%j.out',f'--wrap=module load anaconda3/2025.06; python -B {new}/collect_e5f_surprise_night.py --root {new}'],text=True).strip().split(';')[0]
    receipt=dict(job=job,collector=collector,plan_sha256=sha(new/'plan.json'),tests=19,total_seconds=plan['total_seconds'],seed_steps=plan['seed_steps'])
    (new/'submission.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
