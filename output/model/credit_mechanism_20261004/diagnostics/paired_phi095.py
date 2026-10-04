"""Prepare one phi=.95 solve against an authenticated retained phi=.8 control.

Default is zero-solve preparation. --run requires successful exact input/grid/
price and the lead-reviewed reached-runtime checks; no economic/native gate
changes or fallback allowed. The failed strict-inventory receipt is retained.
"""
from pathlib import Path
import argparse, hashlib, importlib.metadata, json, os, signal, sys, time, resource
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'code/model'))
from production.inputs import load_inputs, DEFAULT_PRICE
from production.credit import bind_engine_credit
REF=ROOT/'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation'
# Lead approved these exact exclusions after checking the fixed-price callgraph.
# Calibration/workflow frontends and their tests are outside this runtime.
EXCLUSIONS={
 'code/model/production/calibration.py':'Calibration adapter is neither imported nor called by solve_at_price.',
 'code/model/production/workflow.py':'GE publication frontend is neither imported nor called by solve_at_price.',
 'code/model/production/test_storage.py':'Test file is not imported or executed by this diagnostic.',
 'code/model/production/parameter_files.py':'Parameter-file frontend is not used; inputs.load_inputs is called directly.',
 'code/model/production/test_parameter_files.py':'Test file is not imported or executed by this diagnostic.',
 'code/model/production/test_parameter_frontend.py':'Test file is not imported or executed by this diagnostic.',
}

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for v in iter(lambda:f.read(1<<20),b''): h.update(v)
    return h.hexdigest()

def encode(v):
    if isinstance(v,np.ndarray): return v.tolist()
    if isinstance(v,np.generic): return v.item()
    return str(v)

def write(p,v): p.write_text(json.dumps(v,indent=2,default=encode)+'\n')

def imported_runtime(expected):
    modules={}
    base=(ROOT/'code/model/production').resolve()
    for name,mod in list(sys.modules.items()):
        filename=getattr(mod,'__file__',None)
        if not filename: continue
        path=Path(filename).resolve()
        if not path.is_relative_to(base): continue
        rel=str(path.relative_to(ROOT))
        if rel in EXCLUSIONS: raise RuntimeError('Excluded module was loaded: '+rel)
        if rel not in expected or sha(path)!=expected[rel]: raise RuntimeError('Unpinned loaded runtime: '+rel)
        modules[name]=dict(path=rel,sha256=expected[rel])
    return modules

def prepare(kind='credit095'):
    old=json.loads((REF/'prepared.json').read_text()); mismatch=[]
    production={k:v for k,v in old['source_sha256'].items() if k.startswith('code/model/production/')}
    expected={k:v for k,v in old['source_sha256'].items() if k not in EXCLUSIONS}
    inventory={str(p.relative_to(ROOT)) for p in (ROOT/'code/model/production').rglob('*.py')}
    inventory.update(('code/model/production/reference_inputs/bundle.json','code/model/production/reference_inputs/arrays.npz'))
    if set(production)-set(EXCLUSIONS)!=inventory-set(EXCLUSIONS): mismatch.append('unexcluded production source inventory differs')
    for k,v in expected.items():
        if sha(ROOT/k)!=v: mismatch.append(k)
    P,grid=load_inputs(); bind_engine_credit(P,'corrected',float(P.unsecured_credit_limit))
    saved=json.loads((REF/'phi_080/executed_P.json').read_text())
    live=json.loads(json.dumps(vars(P),default=encode))
    public0={k:v for k,v in saved.items() if not k.startswith('_')}
    public1={k:v for k,v in live.items() if not k.startswith('_')}
    fields=sorted(k for k in set(public0)|set(public1) if public0.get(k)!=public1.get(k))
    with np.load(REF/'phi_080/solution_arrays.npz',allow_pickle=False) as z:
        grid_equal=bool(np.array_equal(z['b_grid'],grid)); price_equal=bool(z['p_eq'][0]==DEFAULT_PRICE)
    manifest=json.loads((ROOT/'output/model/experiments/credit_distribution_report_v1/manifest.json').read_text())['cases']['reference_phi08']
    if sha(REF/'phi_080/solution_arrays.npz')!=manifest['source_sha256'] or sha(REF/'phi_080/executed_P.json')!=manifest['parameter_sha256']:
        mismatch.append('retained baseline arrays/input hashes differ from original manifest')
    Q,qgrid=load_inputs(external_inputs={'phi':[.95]*4} if kind=='credit095' else {}); bind_engine_credit(Q,'corrected',float(Q.unsecured_credit_limit))
    qlive=json.loads(json.dumps(vars(Q),default=encode))
    changes=sorted(k for k in set(live)|set(qlive) if live.get(k)!=qlive.get(k))
    expected_changes=['phi'] if kind=='credit095' else []
    target_price=DEFAULT_PRICE if kind=='credit095' else 1.1*DEFAULT_PRICE
    prefix='phi095' if kind=='credit095' else 'price110'
    passed=not mismatch and not fields and grid_equal and price_equal and changes==expected_changes and np.array_equal(grid,qgrid)
    from production.equilibrium import solve_at_price
    from production.engine.diagnostics import write_diagnostics
    runtime=imported_runtime(expected)
    versions={k:importlib.metadata.version(k) for k in ('numpy','matplotlib','numba','llvmlite','scipy')}
    versions['python']=sys.version; versions['executable']=sys.executable
    receipt=dict(status='prepared_passed' if passed else 'blocked_identity',kind=kind,contract='lead-approved exact reached-runtime identity; economic and native gates unchanged',strict_failure_receipt=str(HERE/'phi095_strict_inventory_failure.json'),exclusions=EXCLUSIONS,production_mismatches=mismatch,cached_public_P_mismatches=fields,grid_identical=grid_equal,baseline_price_identical=price_equal,changed_fields=changes,price_change=dict(baseline=DEFAULT_PRICE,diagnostic=target_price,multiplier=target_price/DEFAULT_PRICE),economic_change=('uniform financed share .8 to .95' if kind=='credit095' else 'price +10% jointly changes rents, purchase cost, collateral and inherited owner housing wealth; all P, grid, entry and taxes unchanged'),source_sha256=expected,imported_modules_before=runtime,versions=versions,baseline_arrays_sha256=sha(REF/'phi_080/solution_arrays.npz'),baseline_parameters_sha256=sha(REF/'phi_080/executed_P.json'),price=target_price,case_budget_seconds=600,memory_budget_gib=8,case_count=1,phase_budget_seconds=1500,threads={k:os.environ.get(k) for k in ['NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},standard_graph_command='production.engine.diagnostics.write_diagnostics(solution,P,destination/standard_diagnostics)',expected_standard_graphs=17)
    write(HERE/(prefix+'_preparation.json'),receipt)
    if not passed: raise RuntimeError('Paired source/input identity blocked; see phi095_preparation.json')
    return Q,qgrid,receipt

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run',action='store_true'); ap.add_argument('--supervised-rss',action='store_true'); ap.add_argument('--kind',choices=['credit095','price110'],default='credit095'); args=ap.parse_args()
    Q,grid,receipt=prepare(args.kind)
    if not args.run: print(json.dumps(receipt,indent=2)); return
    if not args.supervised_rss or os.environ.get('CREDIT_RSS_SUPERVISED')!='1':
        raise RuntimeError('Numerical launch requires approved external RSS supervisor')
    if any(v!='1' for v in receipt['threads'].values()): raise RuntimeError('All thread caps must be 1')
    dest=HERE/('phi_095_run1' if args.kind=='credit095' else 'price_110_run1')
    if dest.exists(): raise RuntimeError('Refusing an existing case; no automatic retry')
    dest.mkdir(); started=time.monotonic(); write(dest/'status.json',dict(status='started',time_epoch=time.time()))
    def alarm(*_): raise TimeoutError('600-second exact paired case budget')
    signal.signal(signal.SIGALRM,alarm); signal.setitimer(signal.ITIMER_REAL,600)
    try:
        from production.equilibrium import solve_at_price
        from production.engine.diagnostics import write_diagnostics
        result=solve_at_price(Q,grid,receipt['price']); sol=result['solution']; P=result['P']
        arrays={k:v for k,v in vars(sol).items() if isinstance(v,np.ndarray) and v.dtype!=object}
        arrays.update({'shared.'+k:v for k,v in vars(result['shared']).items() if isinstance(v,np.ndarray) and v.dtype!=object})
        np.savez_compressed(dest/'solution_arrays.npz',**arrays); write(dest/'executed_P.json',vars(P))
        write_diagnostics(sol,P,dest/'standard_diagnostics')
        plots=list((dest/'standard_diagnostics').glob('*.png'))
        if len(plots)!=17: raise RuntimeError('Expected all 17 standard diagnostics')
        for k,v in receipt['source_sha256'].items():
            if sha(ROOT/k)!=v: raise RuntimeError('Source drift during paired solve: '+k)
        write(dest/'runtime_after.json',dict(imported_modules_after=imported_runtime(receipt['source_sha256']),versions=receipt['versions']))
        write(dest/'status.json',dict(status='completed_fixed_price_diagnostic',elapsed_seconds=time.monotonic()-started,standard_figures=len(plots),arrays_sha256=sha(dest/'solution_arrays.npz'),parameters_sha256=sha(dest/'executed_P.json'),g_mass=float(sol.g.sum()),array_count=len(arrays),resource_peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    except Exception as e:
        write(dest/'failure.json',dict(error=repr(e),elapsed_seconds=time.monotonic()-started)); raise
    finally: signal.setitimer(signal.ITIMER_REAL,0)
    print((dest/'status.json').read_text())

if __name__=='__main__': main()
