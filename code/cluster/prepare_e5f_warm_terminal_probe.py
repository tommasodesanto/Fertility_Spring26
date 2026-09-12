"""Generate a source-pinned numerical warm-start entrypoint; economic sources stay fixed."""
from pathlib import Path
import argparse,copy,difflib,hashlib,json,subprocess
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--no-shock',action='store_true');parser.add_argument('--price-multiplier',type=float,choices=(1.,1.05),default=1.);args=parser.parse_args()
    if args.price_multiplier!=1. and not args.no_shock:raise ValueError('Price perturbation is restricted to the no-shock feasibility probe')
    root=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
    name='irf_no_shock_terminal_20260912' if args.no_shock else 'night_warm_terminal_20260912'
    if args.price_multiplier!=1.:name+='_price105'
    out=root/'batches'/name;out.mkdir(exist_ok=False)
    plan=read(root/'batches/night_surprises_recovery_20260912/plan.json')
    warm=root/'batches/night_surprise_frontier_20260912/results/arm_0/trial_00_2007/terminal/root_receipt.json'
    if args.no_shock:
        warm=root/'batches/night_warm_followup_20260912/results/arm_0/trial_00_2007/terminal/root_receipt.json'
    receipt=read(warm);assert receipt['converged'];final=receipt['final']
    start=dict(asset_price=final['prices'][0]*args.price_multiplier,asset_price_multiplier=args.price_multiplier,pension_period=final['fiscal_values'][0],source_receipt=dict(path=str(warm),sha256=sha(warm)))
    original=(root/'code/model/tools/run_e5f_candidate_terminal.py').read_text();source=original
    changes={
        "ROOT = Path(__file__).resolve().parents[3]":f"ROOT = Path({str(root)!r})",
        "SCHEMA = 'e5f_candidate_terminal_v1'":"SCHEMA = 'e5f_candidate_terminal_warm_start_v1'",
        "'endpoint_controls', 'audit_controls', 'root_controls', 'standard_graph_count'}":"'endpoint_controls', 'audit_controls', 'root_controls', 'standard_graph_count', 'diagnostic_root_start'}",
        "        start_pension = float(P.pension)":"""        start_pension = float(P.pension)
        numerical_start = c['diagnostic_root_start']
        warm_spec = numerical_start['source_receipt']
        verify(warm_spec['path'], warm_spec['sha256'])
        warm_receipt = json.loads(Path(warm_spec['path']).read_text())
        if (not warm_receipt['converged']
                or numerical_start.get('asset_price_multiplier',1.) not in (1.,1.05)
                or numerical_start['asset_price'] != warm_receipt['final']['prices'][0]*numerical_start.get('asset_price_multiplier',1.)
                or numerical_start['pension_period'] != warm_receipt['final']['fiscal_values'][0]):
            raise ValueError('Numerical warm start must reproduce its pinned receipt and explicit price multiplier')
        start_price = float(numerical_start['asset_price'])
        start_pension = float(numerical_start['pension_period'])""",
        "status='verified initial equilibrium values; not a solved terminal state')":"status='explicit numerical warm start from a pinned nearby terminal; not a solution at current preferences', provenance=c['diagnostic_root_start'])",
    }
    for before,after in changes.items():
        assert source.count(before)==1,before
        source=source.replace(before,after)
    compile(source,'warm_terminal_driver.py','exec')
    driver=out/'warm_terminal_driver.py';driver.write_text(source)
    (out/'entrypoint.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),source.splitlines(True),fromfile='native',tofile='warm_entrypoint')))
    deltas=(0.,) if args.no_shock else (-.01414,-.0175)
    for i,delta in enumerate(deltas):
        c=copy.deepcopy(plan['terminal_template']);c.update(schema='e5f_candidate_terminal_warm_start_v1',psi_change_from_initial=delta,diagnostic_root_start=start)
        path=out/f'contract_{i}.json';path.write_text(json.dumps(c,indent=2)+'\n');(out/f'contract_{i}.sha256').write_text(sha(path)+'\n')
    driver_pin=sha(driver)
    script=f'''#!/bin/bash
#SBATCH --job-name=e5f_warm_terminal
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-{len(deltas)-1}
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output={out}/slurm_%A_%a.out
#SBATCH --error={out}/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export PYTHONPATH="{root}/code/model/tools:{root}/code/model"
export NUMBA_CACHE_DIR="{root}/output/cache/numba"
echo '{driver_pin}  {driver}' | sha256sum --check --status
python -B {driver} --contract {out}/contract_${{SLURM_ARRAY_TASK_ID}}.json --contract-sha256 "$(cat {out}/contract_${{SLURM_ARRAY_TASK_ID}}.sha256)" --output {out}/arm_${{SLURM_ARRAY_TASK_ID}}
'''
    launch=out/'submit.sh';launch.write_text(script);subprocess.run(['bash','-n',str(launch)],check=True)
    job=subprocess.check_output(['sbatch','--parsable',str(launch)],text=True).strip().split(';')[0]
    result=dict(job=job,deltas=list(deltas),price_multiplier=args.price_multiplier,initial_prices_and_pension=start,driver_sha256=driver_pin,warm_source_sha256=sha(warm),maximum_mappings_per_case=8,maximum_seconds_per_case=1800,changes='Only numerical root starting coordinates, explicit schema and entrypoint provenance; original scientific source manifest retained')
    (out/'submission.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
