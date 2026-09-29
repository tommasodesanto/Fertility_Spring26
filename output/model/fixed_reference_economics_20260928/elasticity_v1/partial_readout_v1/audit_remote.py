"""Read-only compact audit of six passed cases and failed +1 attempt on Torch."""
import csv, hashlib, json
from pathlib import Path
root=Path('/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929/results_v2')
main=root/'solve_v4'; old=root/'solve_v2'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def read(p): return json.loads(p.read_text())
def count(p): return len(list(csv.DictReader(p.open()))) if p.exists() else 0
records=read(main/'latest_completed.json')['completed']; rows=[]
for item in records:
    case=(old if item['case'] in ('grid_control','grid_control_repeat','credit','credit_repeat') else main)/item['case']
    receipt=case/'receipt.json'; r=read(receipt); cp=case/'conditional_cohort_state.pkl.gz'
    plots=sorted((case/'standard_diagnostics').glob('*.png'))
    rows.append(dict(case=item['case'],regime=item['regime'],factor=item['factor'],wall_seconds=item['wall_seconds'],
        receipt_sha256=item['receipt_sha256'],receipt_pin_ok=sha(receipt)==item['receipt_sha256'],
        receipt_status=r['status'],checkpoint_exists=cp.exists(),
        checkpoint_pin_ok=cp.exists() and sha(cp)==r['checkpoint']['sha256'],
        checkpoint_bytes=cp.stat().st_size if cp.exists() else 0,
        fit_rows=count(case/'target_fit.csv'),parameter_rows=count(case/'parameters.csv'),
        plot_count=len(plots),plot_names=[p.name for p in plots],plot_receipt_count=r['standard_plot_count'],
        source_manifest=r['reference_manifest_sha256'],plan_sha=r['plan_sha256'],
        cohort_gates_exist=(case/'gates.json').exists(),
        impact_gates_exist=(case/'baseline_state_impact/gates.json').exists(),
        exact_q0_repeat=r.get('exact_q0_repeat'),numerical_refinement=r.get('numerical_refinement'),
        completed_fertility=r['completed_fertility'],cohort_summary=r['cohort_summary'],
        impact_summary=r['impact_summary'],price=r['price'],rent=r['rent']))
failed=main/'reference_1010'; cp=failed/'conditional_cohort_state.pkl.gz'
info=dict(progress=read(failed/'progress.json') if (failed/'progress.json').exists() else None,
    checkpoint_exists=cp.exists(),checkpoint_bytes=cp.stat().st_size if cp.exists() else 0,
    checkpoint_gzip_header=cp.open('rb').read(2).hex() if cp.exists() else None,
    receipt_exists=(failed/'receipt.json').exists(),fit_rows=count(failed/'target_fit.csv'),
    parameter_rows=count(failed/'parameters.csv'),gates_exists=(failed/'gates.json').exists(),
    impact_gates_exists=(failed/'baseline_state_impact/gates.json').exists(),
    plot_count=len(list((failed/'standard_diagnostics').glob('*.png'))),
    case_failure=read(failed/'failure.json') if (failed/'failure.json').exists() else None,
    controller_failure=read(main/'failure.json') if (main/'failure.json').exists() else None,
    log_tail=(main/'reference_1010.log').read_text(errors='replace').splitlines()[-18:])
print(json.dumps(dict(reference_label='2007 stationary reference — block0506, September 28 verified export',
    original_job=18801318,continuation_job=18803216,original_deadline_epoch=1790700222.6872504,
    plan_sha=sha(root/'plan_v2.json'),latest_completed_count=len(rows),cases=rows,failed_plus_one=info),
    indent=2,sort_keys=True,default=str))
