"""Parallel 28/56-date checks of the candidate matching the first short window."""
from pathlib import Path
import copy,hashlib,json,os,subprocess,sys
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a');batch=root/'batches/night_matched_horizon_20260912'
    matched=root/'batches/night_warm_followup_20260912/results/arm_0/trial_00_2007'
    native=matched/'forecast_6_1/vintage/2007/root_receipt.json';assert read(native)['finite_horizon_market_fiscal_converged']
    initial=root/'batches/night_surprises_recovery_20260912/results/arm_1/trial_00_2007/forecast_28_2/vintage/2007/root_receipt.json';assert read(initial)['finite_horizon_market_fiscal_converged']
    oldterminal=root/'batches/night_surprises_recovery_20260912/results/arm_1/trial_00_2007/terminal/root_receipt.json'
    newterminal=matched/'terminal/root_receipt.json';ratio=read(newterminal)['final']['prices'][0]/read(oldterminal)['final']['prices'][0]
    base=read(root/'batches/night_warm_followup_20260912/plan_0.json')
    for arm,count in enumerate((28,56)):
        p=copy.deepcopy(base);p.update(output_root=str(batch/'results'),seed_steps=[-.01414,-.01414],
            total_seconds=9000,policy_reserve_seconds=0,native_smoke_only=False,forecast_diagnostic_only=True,
            forecast_counts=[count],verified_native_smoke=str(native),initialization_receipt=str(initial),initial_price_multiplier=ratio,
            diagnostic_question='How does the first-window fit change with the forecast horizon? Fixed preferences; full initial targets unchanged; no accepted historical sequence or policy.')
        p['history_root_controls']['max_evaluations']=8 if count==28 else 6
        p['file_sha256'].update({str(f):sha(f) for f in [*batch.glob('*.py'),native,initial,newterminal,oldterminal]})
        (batch/f'plan_{arm}.json').write_text(json.dumps(p,indent=2)+'\n')
    env=dict(os.environ,PYTHONPATH=f'{batch}:{root}/code/model/tools:{root}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','test_e5f_successive_surprises','-q'],env=env,check=True,timeout=180)
    s=f'''#!/bin/bash
#SBATCH --job-name=e5f_matched_horizon
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:50:00
#SBATCH --output={batch}/slurm_%A_%a.out
#SBATCH --error={batch}/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export PYTHONPATH="{batch}:{root}/code/model/tools:{root}/code/model"
export NUMBA_CACHE_DIR="{root}/output/cache/numba"
python -B {batch}/run_e5f_successive_surprises_overnight.py --plan {batch}/plan_${{SLURM_ARRAY_TASK_ID}}.json --arm "$SLURM_ARRAY_TASK_ID"
'''
    script=batch/'submit.sh';script.write_text(s);subprocess.run(['bash','-n',str(script)],check=True)
    job=subprocess.check_output(['sbatch','--parsable',str(script)],text=True).strip().split(';')[0]
    receipt=dict(job=job,counts=[28,56],budget_per_arm_seconds=9000,native_smoke_sha256=sha(native),plan_sha256=[sha(batch/f'plan_{i}.json') for i in range(2)],initialization='Nearby converged 28-date path, endpoint price ratio, exact matched short prefix; old Jacobian only for matching dimension; numerical guesses only')
    (batch/'submission.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
