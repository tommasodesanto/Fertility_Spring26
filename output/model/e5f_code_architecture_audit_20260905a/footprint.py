from pathlib import Path
import ast,collections,csv,hashlib,json,sys,platform
root=Path.cwd(); model=root/'code/model'; frozen=root/'tmp/e5f_overnight_20260905a/frozen_production/code/model'; out=root/'output/model/e5f_code_architecture_audit_20260905a'
sys.path[:0]=[str(model/'tools'),str(model)]
import run_e5f_transition_calibration as calibration
chain,solver=calibration.transition.configure_sequential_model()
loaded=sorted({Path(m.__file__).resolve() for m in sys.modules.values() if getattr(m,'__file__',None) and str(Path(m.__file__).resolve()).startswith(str(model)) and str(m.__file__).endswith('.py')})
rows=[]; funcs=[]
for path in sorted((model/'intergen_eqscale_seq_optimized').glob('*.py'))+sorted((model/'intergen_eqscale_seq_optimized/tests').glob('*.py'))+loaded:
 if any(r['path']==str(path.relative_to(root)) for r in rows): continue
 s=path.read_text(); t=ast.parse(s); fs=[n for n in ast.walk(t) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))]
 rel=str(path.relative_to(root)); f=frozen/path.relative_to(model)
 rows.append(dict(path=rel,lines=len(s.splitlines()),nonblank=sum(bool(l.strip()) for l in s.splitlines()),functions=len(fs),loaded_during_initialization=path in loaded,sha256=hashlib.sha256(s.encode()).hexdigest(),frozen_exists=f.exists(),same_as_frozen=f.exists() and f.read_bytes()==path.read_bytes()))
 for n in fs: funcs.append(dict(path=rel,name=n.name,line=n.lineno,end=n.end_lineno,lines=n.end_lineno-n.lineno+1))
for name,data in [('footprint.csv',rows),('functions.csv',funcs)]:
 with (out/name).open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(data[0])); w.writeheader();w.writerows(data)
current=calibration.code_fingerprint_contract(solver)
summary=json.loads((root/'output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/summary.json').read_text())
prod=summary['code_fingerprints']
print('production fingerprints keys',list(prod))
payload={'python':sys.version,'numpy':calibration.np.__version__,'numba_available':solver.NUMBA_AVAILABLE,'loaded_modules':[str(p.relative_to(root)) for p in loaded],'loaded_lines':sum(r['lines'] for r in rows if r['loaded_during_initialization']),'package_top_level_lines':sum(r['lines'] for r in rows if Path(r['path']).parent.name=='intergen_eqscale_seq_optimized'),'package_top_level_files':sum(Path(r['path']).parent.name=='intergen_eqscale_seq_optimized' for r in rows),'current_bundle':current,'production_bundle':prod,'largest_functions':sorted(funcs,key=lambda x:-x['lines'])[:20]}
(out/'architecture_summary.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({k:v for k,v in payload.items() if k not in ['current_bundle','production_bundle','largest_functions']},indent=2));print(json.dumps(payload['largest_functions'][:10],indent=2))
