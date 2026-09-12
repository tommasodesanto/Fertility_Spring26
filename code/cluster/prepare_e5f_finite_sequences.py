"""Launch bounded parallel diagnostic surprise sequences from pinned model sources."""
from pathlib import Path
import copy,hashlib,json,os,subprocess,sys

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    batch=root/'batches/finite_sequences_20260912'
    prior=root/'batches/matched_continuation_20260912'
    base=read(prior/'plan_0.json')
    native=Path(base['verified_native_smoke'])
    long=prior/'results/arm_0/trial_00_2007/forecast_28_1/vintage/2007/root_receipt.json'
    patch=root/'batches/stationary_2019_pf_test_20260912'
    patch_receipt=patch/'results/arm_0/trial_00_2019/forecast_6_2/vintage/2019/root_receipt.json'
    for receipt in (native,long,patch_receipt):
        assert read(receipt)['finite_horizon_market_fiscal_converged'],str(receipt)
    controls=[]
    for i,(n,hours) in enumerate(((6,3),(12,5),(28,8),(6,3))):
        p=copy.deepcopy(base if i<3 else read(patch/'plan.json'))
        receipt=(native if i==0 else long) if i<3 else patch_receipt
        p.update(output_root=str(batch/'results'),total_seconds=hours*3600,
            finite_sequence_diagnostic=True,skip_policies=True,policy_reserve_seconds=0,
            native_smoke_only=False,forecast_diagnostic_only=False,forecast_counts=[n],
            maximum_trials_per_window=6,seed_steps=[-.01414]*4,proposal_step=-.005,
            seed_steps_by_year={'2011':-.010,'2015':-.012,'2019':-.014},
            initialization_receipt=str(receipt),initial_price_multiplier=1.,initialization_native_prefix=False,
            maximum_terminal_solves_per_arm=24 if i<3 else 6,
            maximum_forecast_root_mappings_per_trial=24,
            diagnostic_question='Fit successive unanticipated preference shocks with finite forecast market/PAYGO/replay gates; retain failed terminal-distance status. Only first-period states carry. No policy or production admission.')
        p['history_root_controls']['max_evaluations']=8
        if i==3:
            p['initial_psi_by_year']={'2019':.10239514522037684}
            p.pop('verified_native_smoke',None)
            p['diagnostic_question']='Fit last shock from conditional stationary2019households;2023is on the transition. Finite-horizon patch only.'
        pinned=[receipt,*batch.glob('*.py')]
        p['file_sha256'].update({str(q):sha(q) for q in pinned})
        plan=batch/f'plan_{i}.json';plan.write_text(json.dumps(p,indent=2)+'\n')
        controls.append(dict(arm=i,dates=n,mode='successive_sequence' if i<3 else 'stationary2019_patch',
            maximum_trials=24 if i<3 else 6,root_evaluations_per_trial_maximum=24,
            wall_budget_seconds=hours*3600,observed_seconds_per_mapping={6:170,12:350,28:814}[n],
            estimate='Planning range 1-3 trials per window and 2-6 mappings per trial; time cap may prevent completion.',
            plan_sha256=sha(plan)))
    env=dict(os.environ,PYTHONPATH=f'{batch}:{root}/code/model/tools:{root}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','test_e5f_successive_surprises','test_e5f_finite_sequence','-q'],env=env,check=True,timeout=180)
    script=f'''#!/bin/bash
#SBATCH --job-name=e5f_sequence_fit
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:10:00
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
    result=dict(job=job,arms=controls,production_eligible=False,policies_launched=False,
        smoke='25 pure tests including two-vintage loop, finite admission and carry; pinned numerical native loop. Each forecast requires exact reproduction and each carried period requires original-expectations replay.',
        stop='Time cap, six trials/window, three root continuations/trial, numerical failure or inability to fit within0.005. No inaccurate window carried forward.',
        progress='heartbeat.json every60s; latest_completed/best_so_far per root mapping; realized_fit and exact-reloaded next_state checkpoints per fitted window;17 standard graphs and2023ageallocation per valid trial.')
    (batch/'submission.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
