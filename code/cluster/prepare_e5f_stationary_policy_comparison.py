"""Pinned baseline replay then two parallel stationary rebate comparisons."""
from pathlib import Path
from types import SimpleNamespace as NS
import gzip,hashlib,json,os,pickle,subprocess,sys
R=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
B=R/'batches/stationary_policy_comparison_fit_20260912'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=read(R/'batches/finite_sequences_20260912/plan_0.json');c=p['terminal_template']
    sys.path[:0]=[str(B),str(R/'code/model/tools'),str(R/'code/model')]
    import run_e5f_candidate_terminal as td
    td.validate_contract(c);td.verify_sources(c)
    for k in ('initial_checkpoint','initial_summary','initial_contract'):td.verify(c[k]['path'],c[k]['sha256'])
    for k,v in p['file_sha256'].items():td.verify(k,v)
    with gzip.open(c['initial_checkpoint']['path'],'rb') as f:seed=pickle.load(f)
    selected=read(R/'batches/finite_sequences_20260912/patch_retry/results/arm_3/summary.json')['realized'][0]
    assert selected['error_abs']<=.005
    prior=Path(selected['folder']).parent/'terminal'
    with gzip.open(prior/'terminal_state.pkl.gz','rb') as f:t=pickle.load(f)
    tr=read(prior/'root_receipt.json');assert tr['converged']
    q=t['parameters'];e=t['evaluation'];g=e.g_current
    assert q.psi_child==selected['psi'] and q.tau_H==.04 and q.property_tax_lump_sum_transfer==0
    old=NS(parameters=seed['parameters'],b_grid=seed['b_grid'],policy=seed['evaluation'].policy,supply_rule=seed['supply_rule'])
    packet=dict(plan=p,old=old,psi=q.psi_child,demographics=t['demographic_seed'],terminal_only=True,
        terminal_guess=[float(t['policy'].price[0]),float(q.pension)],
        rebate_guess=float(t['endpoint'].residuals['tax_revenue']/g.sum()),
        baseline_reference=dict(asset_price=float(t['policy'].price[0]),pension_period=float(q.pension),
            household_heads=float(g.sum()),resident_persons=float(t['fixed_point'].persons.persons.sum()),ownership=float(g[:,1:].sum()/g.sum())))
    with gzip.open(B/'inputs.pkl.gz','wb',compresslevel=1) as f:pickle.dump(packet,f,protocol=5)
    contract=dict(source_root=str(R),initial_checkpoint=c['initial_checkpoint'],terminal_checkpoint=str(prior/'terminal_state.pkl.gz'),terminal_checkpoint_sha256=sha(prior/'terminal_state.pkl.gz'),
        target_fingerprint=p['target_fingerprint'],psi=q.psi_child,annual_taxes=[.01,.01,.02],equal_rebate=[False,True,True],
        source_files={str(f):sha(f) for f in B.glob('*.py')},maximum_outer_mappings_per_case=24,seconds_per_case=1200,
        production_eligible=False,estimated='Structural parameters inherited; no re-estimation',external='Supply elasticity.63, payroll.179, same demographics and fiscal units',outstanding='Preference fitted to finalpatchwindow, but horizon uncertified; entry/migration normalization diagnostic; no policytransition',
        run_size='Three cases, at most24mappings each. Baseline warm replay expected2mappings; rebates2-24mappings each. Prior stationary mappings about40seconds.20min cap percase.')
    (B/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    env=dict(os.environ,PYTHONPATH=f'{B}:{R}/code/model/tools:{R}/code/model')
    subprocess.run([sys.executable,'-B','-m','unittest','test_e5f_surprise_overnight','-q'],env=env,check=True,timeout=120)
    record=dict(jobs=[],contract_sha256=sha(B/'contract.json'))
    for i,case in enumerate(['baseline','equal-rebate-1pct','equal-rebate-2pct']):
        path=B/f'run_{i}.sh';path.write_text(f'''#!/bin/bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export PYTHONPATH="{B}:{R}/code/model/tools:{R}/code/model"
export NUMBA_CACHE_DIR="{R}/output/cache/numba"
python -B {B}/run_e5f_successive_surprise_policy.py --inputs {B}/inputs.pkl.gz --case {case} --output {B}/results/{case} --seconds 1200
''')
        subprocess.run(['bash','-n',str(path)],check=True)
        args=['sbatch','--parsable',f'--job-name=e5f_ss_policy_{i}','--account=torch_pr_570_general','--cpus-per-task=1','--mem=8G','--time=00:25:00',f'--output={B}/slurm_%j.out',f'--error={B}/slurm_%j.err']
        if i:args.append(f"--dependency=afterok:{record['jobs'][0]['job']}")
        job=subprocess.check_output([*args,str(path)],text=True).strip().split(';')[0]
        record['jobs'].append(dict(case=case,job=job));(B/'submission.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))
if __name__=='__main__':main()
