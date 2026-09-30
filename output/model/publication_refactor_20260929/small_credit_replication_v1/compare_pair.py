"""Compare complete pair locally on Torch; no model imports or solves."""
import csv,hashlib,json,sys,subprocess
from pathlib import Path
import numpy as np
a,b,out=map(Path,sys.argv[1:4])
historical_mode="--historical" in sys.argv[4:]
def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
ra,rb=load(a/'completed.json'),load(b/'completed.json')
assert ra['status']==rb['status']=='full_passed'
assert ra['lifecycle_solves']==rb['lifecycle_solves']==6
for key in ['selected_d_bar','input_identity']: assert ra[key]==rb[key],key
# Compare economic closure at every price, excluding elapsed/source path metadata.
closures=[p.relative_to(a) for p in (a/'phase_b_ge').rglob('closure.json')]
assert set(closures)=={p.relative_to(b) for p in (b/'phase_b_ge').rglob('closure.json')}
for rel in closures: assert load(a/rel)==load(b/rel),str(rel)
arrays=[p.relative_to(a) for p in a.rglob('solution_arrays.npz')]
assert set(arrays)=={p.relative_to(b) for p in b.rglob('solution_arrays.npz')}
counts={}
for rel in arrays:
    with np.load(a/rel) as x,np.load(b/rel) as y:
        assert set(x.files)==set(y.files)
        for key in x.files:
            assert np.isfinite(x[key]).all() and np.isfinite(y[key]).all()
            assert np.array_equal(x[key],y[key]),(str(rel),key)
        counts[str(rel)]=len(x.files)
for final in ['selected_root','selected_repeat_final']:
    for name,n in [('target_fit.csv',14),('parameters.csv',31)]:
        rel=Path('phase_b_ge')/final/name
        with (a/rel).open() as x,(b/rel).open() as y:
            xx,yy=list(csv.DictReader(x)),list(csv.DictReader(y))
        assert len(xx)==len(yy)==n and xx==yy,str(rel)
    pa=a/'phase_b_ge'/final;pb=b/'phase_b_ge'/final
    png=[p.relative_to(pa) for p in pa.rglob('*.png')]
    assert len(png)==17 and set(png)=={p.relative_to(pb) for p in pb.rglob('*.png')}
    for rel in png: assert sha(pa/rel)==sha(pb/rel),str(rel)
r=dict(status='passed',closure_paths=len(closures),array_paths=counts,target_rows=14,parameter_rows=31,standard_pngs_per_final=17,lifecycle_evaluations=12,scope='corrected-credit matched scalar/indexed; not pre-correction economics')
timings={}
for arm,root in ([] if historical_mode else [('scalar',a),('indexed',b)]):
    cases={}
    for path in root.rglob('summary.json'):
        item=load(path)
        if 'elapsed_seconds' in item and 'label' in item:
            cases[item['label']]=float(item['elapsed_seconds'])
    assert len(cases)==6,(arm,cases)
    total=float(load(root/'workflow_timing.json')['elapsed_seconds'])
    timing_sum=sum(cases.values())
    timings[arm]=dict(workflow_seconds=total,case_solve_seconds=timing_sum,
                     residual_overhead_seconds=total-timing_sum,cases=cases)
if not historical_mode:
    r['timing']=timings
    r['workflow_speedup']=timings['scalar']['workflow_seconds']/timings['indexed']['workflow_seconds']
    r['solve_speedup']=timings['scalar']['case_solve_seconds']/timings['indexed']['case_solve_seconds']
    historical=Path(sys.argv[4]) if len(sys.argv)>4 else Path('/scratch/td2248/projects/small_credit_v1/results/full')
    assert historical.is_dir(), 'Historical source missing'
    h_out=out.with_name('historical_indexed_comparison.json')
    subprocess.run([sys.executable,__file__,str(b),str(historical),str(h_out),'--historical'],check=True)
    r['historical_indexed_comparison']=load(h_out)
else:
    r['scope']='historical18869900/indexed identity; no historical timing comparison'
out.write_text(json.dumps(r,indent=2)+'\n'); print(json.dumps(r))
