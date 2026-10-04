"""Zero-solve verification of an unapplied stationary diagnostic-routing patch."""
import ast
import copy
import difflib
import hashlib
import json
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[5]
WORK = Path(__file__).resolve().parent
SOURCE = ROOT/'code/model/experiments/transition_readiness/floor_runtime.py'
CALENDAR = ROOT/'code/model/tools/run_dynamic_population_transition.py'
NATIVE_PRICE = ROOT/'code/model/experiments/birth_count_choice/model/native_price.py'
paths = [SOURCE, CALENDAR, NATIVE_PRICE]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
old = SOURCE.read_text()
needle = '        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);P=copy.deepcopy(self.P);P.psi_child=float(psi)\n'
addition = '''        saved_output=getattr(P,'native_inherited_distribution_evidence_dir',None)
        P.native_inherited_distribution_evidence_dir=str(folder/'inherited_state_evidence')
        write(folder/'output_override.json',dict(field='native_inherited_distribution_evidence_dir',
            saved=saved_output,effective=P.native_inherited_distribution_evidence_dir,economic_change=False))
'''
assert old.count(needle) == 1
new = old.replace(needle, needle+addition)
relative = str(SOURCE.relative_to(ROOT))
(WORK/'failure_output_fix_v1.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+relative,tofile='b/'+relative)))
# Execute the exact patched stationary prefix, stopping before native_bindings
# and all numerical callbacks. This verifies routing on its copied parameters.
mod = ast.parse(new)
cls = next(n for n in mod.body if isinstance(n,ast.ClassDef) and n.name=='FloorRuntime')
method = copy.deepcopy(next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='stationary'))
method.body = method.body[:next(i for i,n in enumerate(method.body) if isinstance(n,ast.With))]
method.body.append(ast.Return(value=ast.Name(id='P',ctx=ast.Load())))
namespace = dict(Path=Path,copy=copy,write=lambda p,d:p.write_text(json.dumps(d)),float=float)
exec(compile(ast.fix_missing_locations(ast.Module(body=[method],type_ignores=[])),'patched_stationary_prefix','exec'),namespace)
# Execute unchanged real feasibility helper and exception, with model tolerance
# stub only. Neither import nor execute the full calendar/model runtime.
calendar_tree = ast.parse(CALENDAR.read_text())
nodes = [copy.deepcopy(n) for n in calendar_tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in ('InheritedDistributionInfeasible','_require_exact_inherited_distribution')]
namespace.update(np=np,hashlib=hashlib,json=json,time=time,model=SimpleNamespace(DEAD_MASS_TOL=1e-12))
exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),'unchanged_calendar_guard','exec'),namespace)
guard = namespace['_require_exact_inherited_distribution']
exc_type = namespace['InheritedDistributionInfeasible']
results=[]
with tempfile.TemporaryDirectory() as temp:
    base=Path(temp); stale=base/'read_only_saved_case/root_06/stage/inherited_state_failures'
    saved=SimpleNamespace(psi_child=.17, native_inherited_distribution_evidence_dir=str(stale),age_start=20,da=4)
    obj=SimpleNamespace(P=saved,ctx={},pf=SimpleNamespace(calendar=SimpleNamespace(SolveCounter=lambda:None)))
    original_mkdir=Path.mkdir
    def readonly_mkdir(path,*args,**kwargs):
        if path==stale or stale in path.parents:
            raise OSError(30,'Read-only file system',str(path))
        return original_mkdir(path,*args,**kwargs)
    with patch.object(Path,'mkdir',readonly_mkdir):
        for label,mass,value,should_reject in [('below_tolerance',5e-13,-1e10,False),('above_tolerance',2e-12,-1e10,True),('nonfinite',5e-13,float('nan'),True),('no_bad_values',2e-12,0.,False)]:
            folder=base/label/'endpoint/point_002'
            P=namespace['stationary'](obj,.3,.77,folder)
            assert P is not saved and saved.psi_child==.17 and P.psi_child==.3
            assert saved.native_inherited_distribution_evidence_dir==str(stale)
            assert P.native_inherited_distribution_evidence_dir==str(folder/'inherited_state_evidence')
            route=json.loads((folder/'output_override.json').read_text())
            assert route['saved']==str(stale) and route['economic_change'] is False
            g=np.full((1,1,1,1,1,1,1),mass); unchanged=g.copy()
            policy=SimpleNamespace(V=np.full(g.shape,value),price=np.array([.77]))
            rejected=False
            try: guard(g,policy,P,np.array([0.]))
            except exc_type as error:
                rejected=True; assert error.dead_mass==mass and Path(error.evidence_path).is_file()
            assert rejected==should_reject
            assert np.array_equal(g,unchanged)
            files=list((folder/'inherited_state_evidence').glob('*.json'))
            if value==0.: assert not files
            else:
                assert len(files)==1
                evidence=json.loads(files[0].read_text())
                assert evidence['projection_mass']==0. and evidence['distribution_modified'] is False
                assert evidence['feasibility_mass_tolerance']==1e-12 and evidence['value_cutoff']==-1e9
                assert evidence['status']==('rejected' if should_reject else 'retained_below_existing_feasibility_tolerance')
            results.append(dict(case=label,status='passed',rejected=rejected,evidence_files=len(files)))
        # An unchanged stale parameter reproduces errno30, including a tail
        # that otherwise falls below the existing tolerance.
        try: guard(np.full((1,1,1,1,1,1,1),5e-13),SimpleNamespace(V=np.full((1,1,1,1,1,1,1),-1e10),price=np.array([.77])),saved,np.array([0.]))
        except OSError as error: assert error.errno==30
        else: raise AssertionError('Old stale route did not reproduce errno30')
        results.append(dict(case='old_route_below_tolerance_reproduces_errno30',status='passed'))
after = {str(p.relative_to(ROOT)): sha(p) for p in paths}
assert before==after
receipt=dict(status='passed',new_native_calls=0,active_sources_unchanged=True,source_sha256=before,cases=results,patch_sha256=sha(WORK/'failure_output_fix_v1.patch'))
(WORK/'failure_output_fix_v1_mock_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(status='passed',mock_cases=len(results),new_native_calls=0,active_sources_unchanged=True)))
