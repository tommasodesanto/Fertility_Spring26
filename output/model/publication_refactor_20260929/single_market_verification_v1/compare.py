"""Read-only matched control comparison; no imports of model or solves."""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
control,cleanup,out=map(Path,sys.argv[1:4])
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a,b=load(control/'completed.json'),load(cleanup/'completed.json')
assert a['result']['status']==b['result']['status']=='passed'
assert a['lifecycle_solves']==b['lifecycle_solves']==6
assert a['result']['selected_d_bar']==b['result']['selected_d_bar']==.14
reference=load(control/'effective_input_contract.json')['reference']
assert reference==b['input_identity'],'Input identities differ'
closures=sorted(p.relative_to(control) for p in (control/'phase_b_ge').rglob('closure.json'))
assert set(closures)=={p.relative_to(cleanup) for p in (cleanup/'phase_b_ge').rglob('closure.json')}
for rel in closures:assert load(control/rel)==load(cleanup/rel),str(rel)
arrays=sorted(p.relative_to(control) for p in control.rglob('solution_arrays.npz'))
assert set(arrays)=={p.relative_to(cleanup) for p in cleanup.rglob('solution_arrays.npz')}
counts={}
for rel in arrays:
    with np.load(control/rel) as x,np.load(cleanup/rel) as y:
        assert set(x.files)==set(y.files),(str(rel),'array fields')
        for key in x.files:
            assert np.isfinite(x[key]).all() and np.isfinite(y[key]).all(),(str(rel),key,'nonfinite')
            assert np.array_equal(x[key],y[key]),(str(rel),key,'not exact')
        counts[str(rel)]=len(x.files)
for name in ['selected_root','selected_repeat_final']:
    rel=Path('phase_b_ge')/name
    for filename,n in [('target_fit.csv',14),('parameters.csv',31)]:
        with (control/rel/filename).open() as x,(cleanup/rel/filename).open() as y:
            xx,yy=list(csv.DictReader(x)),list(csv.DictReader(y))
        assert len(xx)==len(yy)==n and xx==yy,str(rel/filename)
    png=sorted(p.relative_to(control/rel) for p in (control/rel).rglob('*.png'))
    assert len(png)==17 and set(png)=={p.relative_to(cleanup/rel) for p in (cleanup/rel).rglob('*.png')}
    for p in png:assert sha(control/rel/p)==sha(cleanup/rel/p),str(rel/p)
result=dict(status='passed',baseline_job=18879780,cleanup_lifecycle_solves=6,control_lifecycle_solves=6,closure_paths=len(closures),array_paths=counts,target_rows=14,parameter_rows=31,standard_pngs_per_final=17,scope='single-market phase1 cleanup vs matched160x15 D14 grid control; excludes common_support_policies.npz; no claim of full location-axis removal or measured speedup')
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
