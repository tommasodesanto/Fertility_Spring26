"""Resume two improving matched-horizon roots with their full saved coordinates."""
from pathlib import Path
import hashlib,json,os,subprocess,sys

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    prior=root/'batches/night_matched_horizon_20260912'
    batch=root/'batches/matched_continuation_20260912'
    controls=[]
    for i,n in enumerate((28,56)):
        p=read(prior/f'plan_{i}.json')
        receipt=prior/f'results/arm_{i}/trial_00_2007/forecast_{n}_0/vintage/2007/root_receipt.json'
        d=read(receipt);f=d.get('final') or d['best'];scores=[h['score'] for h in d['history']]
        assert not d['finite_horizon_market_fiscal_converged']
        assert scores[-1]<scores[0]/10 and len(f['prices'])==n
        assert len(d['final_jacobian'])==2*n
        p.update(output_root=str(batch/'results'),total_seconds=14400,
            initialization_receipt=str(receipt),initial_price_multiplier=1.,
            initialization_native_prefix=False,
            diagnostic_question='Resume an improving root from its entire saved price/pension path and Jacobian; unchanged first-window preference and economic contract.')
        p['history_root_controls']['max_evaluations']=8 if n==28 else 6
        p['file_sha256'].update({str(q):sha(q) for q in [receipt,*batch.glob('*.py')]})
        (batch/f'plan_{i}.json').write_text(json.dumps(p,indent=2)+'\n')
        controls.append(dict(dates=n,initial_score=scores[0],saved_score=scores[-1],
            receipt=str(receipt),receipt_sha256=sha(receipt),plan_sha256=sha(batch/f'plan_{i}.json'),
            seconds_per_observed_mapping=814 if n==28 else 1554,
            expected_mapping_capacity_in_budget=16 if n==28 else 8))
    env=dict(os.environ,PYTHONPATH=f'{batch}:{root}/code/model/tools:{root}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','test_e5f_successive_surprises','-q'],env=env,check=True,timeout=180)
    script=f'''#!/bin/bash
#SBATCH --job-name=e5f_matched_continue
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:10:00
#SBATCH --output={batch}/slurm_%A_%a.out
#SBATCH --error={batch}/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export PYTHONPATH="{batch}:{root}/code/model/tools:{root}/code/model"
export NUMBA_CACHE_DIR="{root}/output/cache/numba"
python -B {batch}/run_e5f_successive_surprises_overnight.py --plan {batch}/plan_${{SLURM_ARRAY_TASK_ID}}.json --arm "$SLURM_ARRAY_TASK_ID"
'''
    path=batch/'submit.sh';path.write_text(script)
    subprocess.run(['bash','-n',str(path)],check=True)
    job=subprocess.check_output(['sbatch','--parsable',str(path)],text=True).strip().split(';')[0]
    result=dict(job=job,arms=controls,budget_seconds=14400,
        stop='Finite root closure, unchanged gate failure, or four-hour budget; no automatic repetition of this continuation.',
        smoke='19 loop tests plus pinned native equilibrium smoke; exact full-coordinate equality checked before first root mapping',
        scientific_changes=False,historical_fit_complete=False)
    (batch/'submission.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
