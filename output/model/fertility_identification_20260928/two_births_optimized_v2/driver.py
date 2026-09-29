"""Isolated Torch experiment: at most two optimized births per four-year cell.

Original files remain immutable. Authenticated full-source copies describe the
effective overlay; changed function bodies are installed only in this process.
The reference is never promoted or overwritten. No transition is run.
"""
from __future__ import annotations
import argparse
import ast
import copy
import csv
import difflib
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import time
import types
import __future__

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE = ROOT/'output/model/fertility_identification_20260928'
HERE = BASE/'two_births_optimized_v2'
REF = BASE/'resume_v1/selected_export/primary'
LABEL = '2007 stationary reference — block0506, September 28 verified export'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def read(path): return json.loads(Path(path).read_text())


def graft(module, original, changed, generated):
    """Preserve function identity so previously imported aliases use new code."""
    old_tree=ast.parse(original); new_tree=ast.parse(changed)
    classes=(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)
    other=lambda tree:[ast.dump(n,include_attributes=False) for n in tree.body if not isinstance(n,classes)]
    if other(old_tree)!=other(new_tree):
        raise RuntimeError('Overlay changes top-level imports/assignments: '+str(generated))
    old={n.name:n for n in old_tree.body if isinstance(n,classes)}
    altered=[n for n in new_tree.body if isinstance(n,classes) and
             (n.name not in old or ast.dump(n)!=ast.dump(old[n.name]))]
    installed=[]
    for node in altered:
        previous=module.__dict__.get(node.name)
        unit=ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[]))
        exec(compile(unit,str(generated),'exec',flags=__future__.annotations.compiler_flag),module.__dict__)
        current=module.__dict__[node.name]
        if isinstance(previous,types.FunctionType) and isinstance(current,types.FunctionType):
            if previous.__closure__ or current.__closure__:
                raise RuntimeError('Unexpected closure in overlay '+node.name)
            previous.__code__=current.__code__
            previous.__defaults__=current.__defaults__
            previous.__kwdefaults__=current.__kwdefaults__
            previous.__annotations__=current.__annotations__
            previous.__doc__=current.__doc__
            module.__dict__[node.name]=previous
        elif previous is not None:
            if not isinstance(node,ast.ClassDef):
                raise RuntimeError('Non-function existing overlay '+node.name)
            # Existing defaults do not refer to this class; refresh direct imports.
            for imported in list(sys.modules.values()):
                mapping=getattr(imported,'__dict__',None)
                if mapping is None: continue
                for key,value in list(mapping.items()):
                    if value is previous: mapping[key]=current
        installed.append(node.name)
    return installed


def load_runtime(out):
    assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit()
    sys.path.insert(0,str(ROOT/'code/model/tools'))
    import run_e5f_fertility_identification as identification
    import e5f_evening_calibration_runtime as runtime
    c,objectives=identification.verify(BASE/'contract_v1/contract.json')
    for name,digest in read(REF/'artifact_hashes.json').items():
        assert sha(REF/name)==digest,name
    receipt=read(REF/'receipt.json')
    assert sha(REF/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256']
    evaluator=runtime.setup(dict(c,objective=c['lanes']['primary']['objective']),
                            objectives['primary'],out/'preparation')
    with gzip.open(REF/'initial_state.pkl.gz','rb') as stream:
        packet=pickle.load(stream)
    point={r['parameter']:float(r['estimate']) for r in csv.DictReader((REF/'parameters.csv').open())
           if r['parameter'] in runtime.FREE}
    assert len(point)==10
    return c,objectives['primary'],evaluator,packet,point,receipt


def install_overlay(out,evaluator):
    import patch_solver
    import patch_reporting
    # Both patchers own a narrow set of authenticated source files.
    paths=sorted(set(patch_solver.SOURCE_PATHS)|set(patch_reporting.SOURCE_PATHS))
    originals={path:(ROOT/path).read_text() for path in paths}
    changed=patch_solver.patch_sources(dict(originals))
    changed=patch_reporting.patch_sources(changed)
    # Runtime sometimes imports a source under a private module name. Load the
    # public names before grafting so later ordinary imports cannot bypass it.
    for path in paths:
        if path.startswith('code/model/tools/'):
            canonical=importlib.import_module(Path(path).stem)
            assert Path(canonical.__file__).resolve()==(ROOT/path).resolve()
    manifest={}
    for path in paths:
        before=originals[path];after=changed[path]
        generated=out/'effective_sources'/path
        generated.parent.mkdir(parents=True,exist_ok=True);generated.write_text(after)
        generated.with_suffix('.diff').write_text(''.join(difflib.unified_diff(
            before.splitlines(True),after.splitlines(True),fromfile=path,tofile=str(generated))))
        modules=[m for m in list(sys.modules.values()) if getattr(m,'__file__',None)
                 and Path(m.__file__).resolve()==(ROOT/path).resolve()]
        if Path(path).name=='e5f_initial_fertility_observer.py':
            # The native runtime deliberately pins this observer to a frozen
            # source copy under a private name. Authenticate equality before
            # applying the same transformation to its live function object.
            frozen=sys.modules[evaluator.rt['observe_initial_fertility'].__module__]
            assert Path(frozen.__file__).read_text()==before,'Frozen fertility definition differs'
            if frozen not in modules: modules.append(frozen)
        if not modules:
            if path.startswith('code/model/tools/'):
                modules=[importlib.import_module(Path(path).stem)]
            else: raise RuntimeError('No loaded model module for '+path)
        names=[]
        for module in modules:
            names.extend(graft(module,before,after,generated))
        manifest[path]=dict(original_sha256=sha(ROOT/path),effective_sha256=sha(generated),
                            generated_source=str(generated),functions_or_classes=names,
                            module_names=[module.__name__ for module in modules],
                            loaded_source_paths={module.__name__:str(module.__file__) for module in modules})
    fertility_function=evaluator.rt['observe_initial_fertility']
    assert Path(fertility_function.__code__.co_filename).is_relative_to(out/'effective_sources')
    write(out/'effective_source_manifest.json',manifest)
    return manifest


def audit_extra_cache(packet, model):
    import numpy as np
    P=packet['parameters'];sol=packet['solution'];evaluation=packet['evaluation']
    arrays=[np.asarray(sol.fert_extra_probs),np.asarray(evaluation.policy.fert_extra_probs),
            np.asarray(P._fert_extra_probs)]
    expected=np.shape(evaluation.policy.fert_probs)[:-1]+(2,2,int(P.n_child_states))
    for a in arrays:
        assert a.shape==expected
        assert np.isfinite(a).all() and np.min(a)>=0. and np.max(a)<=1.
    for i in range(3):
        for k in range(i):
            assert np.array_equal(arrays[i],arrays[k]),'extra cache differs across owners'
            assert not np.shares_memory(arrays[i],arrays[k]),'extra cache ownership alias'
    sums=arrays[0].sum(axis=5)
    maximum_sum_error=0.;occupied_dead_mass=0.;zero_menu_mass=0.;bad_nonzero_mass=0.
    for j in range(int(P.A_f_start)-1,int(P.A_f_end)):
        for n in range(2):
            for m in range(n+1):
                total=sums[:,:,:,j,:,n,m]
                maximum_sum_error=max(maximum_sum_error,float(np.max(np.minimum(abs(total),abs(total-1.)))))
                outer=(evaluation.policy.fert_probs[:,:,:,j,:,1] if n==0 else
                       evaluation.policy.fert2_probs[:,:,:,j,:,1,n-1,m])
                reached=evaluation.g_pre[:,:,:,j,:,n,m]*outer*model.get_fecundity_by_age(P)[j]
                zero_menu_mass+=float(np.sum(reached[total==0.]))
                bad_nonzero_mass+=float(np.sum(reached[(total!=0.) & (abs(total-1.)>1e-12)]))
                occupied_dead_mass=zero_menu_mass+bad_nonzero_mass
    assert maximum_sum_error<1e-12 and occupied_dead_mass<2e-10
    return dict(status='passed',shape=list(expected),minimum=float(arrays[0].min()),
        maximum=float(arrays[0].max()),maximum_menu_sum_error=maximum_sum_error,
        reached_invalid_menu_mass=occupied_dead_mass,reached_zero_menu_mass=zero_menu_mass,
        reached_nonzero_bad_sum_mass=bad_nonzero_mass,three_equal_independent_owners=True)


def worker(mode,out):
    import numpy as np
    out.mkdir(parents=True,exist_ok=False)
    c,objective,evaluator,packet,point,reference_receipt=load_runtime(out)
    manifest=install_overlay(out,evaluator)
    model=evaluator.rt['model']
    if mode=='smoke':
        import unittest
        suites=[]
        for name in ('tests_solver','tests_reporting'):
            testmodule=importlib.import_module(name)
            suites.append(unittest.defaultTestLoader.loadTestsFromModule(testmodule))
        result=unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suites))
        assert result.testsRun>0 and result.wasSuccessful(),'synthetic source/kernel tests failed'
        # Full flag-off Bellman/KFE at the authenticated reference price and psi.
        P=copy.deepcopy(packet['parameters']);P.two_births_per_period=False
        sol=model.solve_markov_income_at_prices(packet['solution'].p_eq,P,packet['b_grid'],verbose=False)
        comparisons={}
        for name in ('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs',
                     'fert_probs','fert2_probs','g','bp_pol_stay','c_pol_stay'):
            a=np.asarray(getattr(sol,name));b=np.asarray(getattr(packet['solution'],name))
            assert a.shape==b.shape,name
            difference=float(np.max(np.abs(a.astype(float)-b.astype(float))))
            assert difference<1e-9,(name,difference)
            comparisons[name]=difference
        write(out/'receipt.json',dict(status='passed',scope='synthetic kernels plus flag-off full reference replay',
            model_solves=1,unit_tests=result.testsRun,reference_label=LABEL,
            reference_checkpoint_sha256=reference_receipt['case_checkpoint_sha256'],comparisons=comparisons,
            effective_source_manifest_sha256=sha(out/'effective_source_manifest.json'),
            source_pins={name:sha(HERE/name) for name in ('driver.py','patch_solver.py','patch_reporting.py','tests_solver.py','tests_reporting.py')}))
        return
    smoke=read(Path(os.environ['TWO_BIRTH_SMOKE_RECEIPT']))
    assert smoke['status']=='passed'
    assert smoke['reference_checkpoint_sha256']==reference_receipt['case_checkpoint_sha256']
    for name,digest in smoke['source_pins'].items(): assert sha(HERE/name)==digest,('changed after smoke',name)
    original_bind=evaluator.bind
    def experimental_bind(self,trial_point):
        P=original_bind(trial_point)
        P.two_births_per_period=True
        return P
    evaluator.bind=types.MethodType(experimental_bind,evaluator)
    changes=[
        'Experimental: one extra optional birth attempt after a successful birth; maximum two new births per four-year cell.',
        'Experimental: extra conditional binary choice has existing later-birth Gumbel scale and its inclusive value; no new fitted parameter.',
        'Experimental: an independent second conception draw uses the existing age-specific success probability.',
        'Measurement convention retained: pre/post linear age interpolation, without separately dated within-period births.',
        'Derived: child-benefit level renormalized to completed fertility2.1 with original demographic renewal check.',
        'All ten reference calibration coordinates, income, entry distribution, housing/fiscal primitives, targets, weights and bounds held fixed.'
    ]
    evaluator.c=copy.deepcopy(evaluator.c);evaluator.c['economic_changes']=changes
    # Numerical start only: use the previously normalized saved diagnostic.
    warm=BASE/'two_births_optimized_v1/run_v1/worker/case'
    assert sha(warm/'receipt.json')=='bf88d21049591d56a3168863253c72844f0855f3a1a87bb3c37d055dc6f361f3'
    assert sha(warm/'parameters.csv')=='0e8a18804eacf827c3f649f9293a85004fa666f1cbda70e232c1662560419a0f'
    warm_receipt=read(warm/'receipt.json')
    assert warm_receipt['case_checkpoint_sha256']=='7eb14f50ae468f1e7b6fe9c46fe5df31b767af8b0996a4cdb14a1299c2050337'
    warm_rows={row['parameter']:row for row in csv.DictReader((warm/'parameters.csv').open())}
    for name,value in point.items(): assert float(warm_rows[name]['estimate'])==value,name
    initial_psi=float(warm_rows['psi_child']['estimate'])
    evaluator.c['normalization']['initial_psi']=initial_psi
    evaluator.c['normalization']['maximum_stationary_solves']=8
    write(out/'numerical_initialization.json',dict(initial_psi=initial_psi,
        source_case=str(warm),source_receipt_sha256=sha(warm/'receipt.json'),
        method='previous normalized diagnostic as numerical starting guess',
        probability_repair='shifted exponential inner probabilities; inclusive value and gates unchanged'))
    receipt=evaluator.evaluate(point,out/'case',deadline_epoch=time.time()+1050,graphs=True)
    with gzip.open(out/'case/initial_state.pkl.gz','rb') as stream:
        experimental_packet=pickle.load(stream)
    cache_audit=audit_extra_cache(experimental_packet,model)
    write(out/'case/extra_cache_audit.json',cache_audit)
    receipt.update(status='verified_experimental_two_birth_equilibrium_not_adopted',reference_label=LABEL,
        reference_checkpoint_sha256=reference_receipt['case_checkpoint_sha256'],
        effective_source_manifest_sha256=sha(out/'effective_source_manifest.json'),
        source_identity_note='Inherited source hash identifies original ancestry; effective overlay separately pinned.',
        held_reference_coordinate_count=10,estimated_coordinates_this_experiment=0,normalized_coordinates=1,
        free_count=0,reference_calibration_coordinate_count=10,extra_cache_audit=cache_audit,
        fertility_plot_scope='Standard probability panels show the outer attempt; extra conditional probabilities are separately audited, not plotted.',
        economic_changes=changes,scientific_promotion=False)
    write(out/'case/receipt.json',receipt)
    rows=list(csv.DictReader((out/'case/parameters.csv').open()))
    for row in rows:
        if row['parameter'] in point: row['status']='held at reference estimate during two-birth diagnostic'
    with (out/'case/parameters.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    write(out/'case/artifact_hashes.json',{str(path.relative_to(out/'case')):sha(path)
        for path in sorted((out/'case').rglob('*')) if path.is_file() and path.name!='artifact_hashes.json'})
    write(out/'receipt.json',dict(status='passed',case=str(out/'case'),loss=receipt['loss'],
        objective_stationary_solves=receipt['objective_stationary_solves'],reference_label=LABEL,
        source_pins=smoke['source_pins'],scientific_promotion=False))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('smoke','evaluate'))
    parser.add_argument('--worker',action='store_true');parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.worker:
        worker(args.mode,args.output);return
    budget=1200
    args.output.mkdir(parents=True,exist_ok=False)
    for name in ('latest_completed.json','best_so_far.json'):
        write(args.output/name,dict(status='none_completed',mode=args.mode))
    log=args.output/'worker.log'
    started=time.time()
    with log.open('w') as stream:
        proc=subprocess.Popen([sys.executable,__file__,args.mode,'--worker','--output',str(args.output/'worker')],
            stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            elapsed=time.time()-started
            write(args.output/'heartbeat.json',dict(status='running',elapsed_seconds=elapsed,pid=proc.pid,
                stage=args.mode,budget_seconds=budget,progress_ledger=str(args.output/'worker/case/stationary_solves.json')))
            if elapsed>=budget:
                os.killpg(proc.pid,signal.SIGTERM)
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                write(args.output/'completion.json',dict(status='timeout',elapsed_seconds=time.time()-started))
                raise SystemExit(2)
            time.sleep(15)
    status='completed' if proc.returncode==0 else 'failed'
    if proc.returncode==0:
        for name in ('latest_completed.json','best_so_far.json'):
            write(args.output/name,read(args.output/'worker/receipt.json'))
    write(args.output/'completion.json',dict(status=status,exit_code=proc.returncode,elapsed_seconds=time.time()-started))
    raise SystemExit(proc.returncode)


if __name__=='__main__':main()
