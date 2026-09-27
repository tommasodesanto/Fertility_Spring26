"""Authenticated current-source baseline/PF runtime, with no generated adapters.

Preparation imports and validates only. ``replay`` is an explicit one-steady-
state baseline regression, not a calibration or a policy/transition approval.
Run in a fresh interpreter; module-cache contamination is rejected.
"""
from __future__ import annotations
import argparse
import copy
import csv
import gzip
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import shutil
import signal
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PORTABLE = ROOT / 'tmp/e5f_overnight_local_20260927/portable'
CONTRACT = PORTABLE / 'night_launch_v4/primary_continuation/production_contract.json'
REFERENCE = PORTABLE / 'night_launch_v4/primary_continuation/search/de_0093/case'
CONTRACT_SHA = '3b770d8c8c22d2b0449b34a575d6353b063bc015d74ce11016dad7e22ed7ca5e'
CHECKPOINT_SHA = '090c9ebda662bf7837c4f4cf1d816159bc9d203a9babe7c70d00d0c4be1e575e'


class EstateAuditContract:
    """Reject origin-blind pinned ledgers when a DUE policy is requested."""
    def __init__(self, module):
        self.module = module
        self.EstateFundingShortfall = module.EstateFundingShortfall

    def audit(self, evaluation, P, b_grid, **kwargs):
        if (bool(getattr(P, 'native_due_stayer_credit', False))
                and not bool(getattr(self.module, 'SUPPORTS_NATIVE_DUE_STAYER_CREDIT', False))):
            raise RuntimeError('DUE requires a newly authenticated origin-specific estate audit contract')
        return self.module.audit(evaluation, P, b_grid, **kwargs)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1<<20), b''): h.update(block)
    return h.hexdigest()


def read(path): return json.loads(Path(path).read_text())


def write(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def require_current_model(model):
    expected=ROOT/'code/model/intergen_eqscale_seq_optimized/solver.py'
    if Path(model.__file__).resolve()!=expected:
        raise RuntimeError('Current-source model identity failed')
    for name,module in list(sys.modules.items()):
        if name.startswith('intergen_eqscale_seq_optimized') and getattr(module,'__file__',None):
            if not Path(module.__file__).resolve().is_relative_to(expected.parent):
                raise RuntimeError('Mixed model package: '+name)


def require_current_runtime(rt):
    require_current_model(rt['model'])
    modules = {'pf':rt['primitive'].pf, 'calendar':rt['primitive'].pf.calendar,
               'transition':rt['primitive'].pf.transition, 'chain':rt['chain'],
               'primitive':rt['primitive']}
    for name,module in modules.items():
        path=Path(module.__file__).resolve()
        if not path.is_relative_to(ROOT/'code/model'):
            raise RuntimeError('Noncurrent runtime module '+name+': '+str(path))


def verify_current_sources(preparation):
    for relative,expected in preparation['current_source_files'].items():
        if sha(ROOT/relative)!=expected:
            raise RuntimeError('Source changed during replay: '+relative)


def compare_arrays(reference, candidate):
    """Full numeric-array census; disagreements remain visible for lead review."""
    def census(packet):
        result={}
        def visit(value,prefix,depth):
            if isinstance(value,np.ndarray):
                result[prefix]=value
            elif depth>0 and (isinstance(value,dict) or hasattr(value,'__dict__')):
                mapping=value if isinstance(value,dict) else vars(value)
                for name,item in mapping.items():
                    if not str(name).startswith('_'): visit(item,prefix+'.'+str(name),depth-1)
        for key in ('b_grid','shared','solution','evaluation','stationary_g_pre'):
            visit(packet[key],key,3)
        return result
    old,new=census(reference),census(candidate)
    records={}
    mass=np.asarray(reference['evaluation'].g_current,dtype=float)
    for name in sorted(set(old)|set(new)):
        if name not in old or name not in new:
            records[name]={'status':'missing_array'};continue
        a,b=old[name],new[name]
        if a.shape!=b.shape:
            records[name]={'status':'shape_mismatch','reference_shape':a.shape,'candidate_shape':b.shape};continue
        if a.dtype.kind not in 'biufc' or b.dtype.kind not in 'biufc':
            records[name]={'status':'nonnumeric','exact':bool(np.array_equal(a,b))};continue
        finite=bool(np.isfinite(a).all() and np.isfinite(b).all())
        delta=np.abs(a.astype(float)-b.astype(float))
        row={'status':'compared','shape':a.shape,'finite':finite,
             'exact':bool(np.array_equal(a,b)),
             'max_abs':float(delta.max(initial=0)) if finite else None,
             'l1':float(delta.sum()) if finite else None}
        if a.shape==mass.shape and finite:
            row['mass_weighted_abs']=float((mass*delta).sum()/mass.sum())
            row['occupied_max_abs']=float(delta[mass>1e-12].max(initial=0))
        records[name]=row
    return {'arrays':records,'array_count':len(records),
            'status':'review_required',
            'scope':'full saved numeric shared, solution, evaluation, grid and distribution arrays; no silent inactive-frontier exclusion'}


def verify_current_reporter(current_path, frozen_path):
    """Authenticate the sole reviewed reporting extension, without executing patches."""
    source=Path(current_path).read_text()
    additions=[(
        "    def evaluate_point(self,*,tax,objective,selected,runtime,point,output,deadline_epoch,graphs=False,\n"
        "                       report_population_scale=1.0):\n"
        "        # Stationary packets have unit household mass. For a closed endpoint,\n"
        "        # express unchanged absolute supply per household only in this reporter.\n"
        "        if not math.isfinite(float(report_population_scale)) or report_population_scale <= 0:\n"
        "            raise ValueError('Report population scale must be finite and positive')\n",
        "    def evaluate_point(self,*,tax,objective,selected,runtime,point,output,deadline_epoch,graphs=False):\n"),(
        ")**P.xi_supply[0])/report_population_scale,float(P.xi_supply[0]))",
        ")**P.xi_supply[0]),float(P.xi_supply[0]))"),(
        "        if report_population_scale != 1.0:\n"
        "            receipt['stationary_report_units'] = dict(\n"
        "                population_scale=float(report_population_scale),\n"
        "                supply='absolute supply divided by population; economic H0 unchanged',\n"
        "                distribution='unit household mass; multiply by population for terminal-distance checks')\n", "")]
    for new,old in additions:
        if source.count(new)!=1: raise RuntimeError('Current reporter differs from reviewed extension')
        source=source.replace(new,old,1)
    if source!=Path(frozen_path).read_text():
        raise RuntimeError('Unreviewed reporting/scoring source change')


def verify_current_recent_observer(current_path, frozen_path):
    """Verify the passive audit extension; every measurement byte stays frozen."""
    source=Path(current_path).read_text()
    start=source.index('# BEGIN REVIEWED DEAD-TAIL LOCATION AUDIT\n')
    end=source.index('# END REVIEWED DEAD-TAIL LOCATION AUDIT\n\n',start)+len('# END REVIEWED DEAD-TAIL LOCATION AUDIT\n\n')
    extension=source[start:end]
    if hashlib.sha256(extension.encode()).hexdigest()!='069a28a7baeb9fd8ec1fce92725e7340d40d67e4c64a93b99ac95e8eecb8ff6c':
        raise RuntimeError('Unreviewed retained-dead-tail audit extension')
    source=source[:start]+source[end:]
    replacements=[('    diagnostic_allow_retained_dead_tail: bool = False,\n',''),('    location_sum_error, location_audit = _audit_location_lottery(\n        location_probs, post, getattr(policy, "V", None),\n        allow_dead_tail=diagnostic_allow_retained_dead_tail,\n        dead_mass_tolerance=model.DEAD_MASS_TOL,\n        dead_value_cutoff=model.DEAD_VALUE_CUTOFF)','    location_sum_error = float(np.max(np.abs(location_probs.sum(axis=3)[post > 0] - 1.),\n                                      initial=0.))\n    _check(location_sum_error, "occupied location probability sum")'),("    if diagnostic_allow_retained_dead_tail:\n        accounting['retained_dead_tail_location_audit'] = location_audit\n",'')]
    for new,old in replacements:
        if source.count(new)!=1:raise RuntimeError('Recent observer audit seam differs')
        source=source.replace(new,old,1)
    if source!=Path(frozen_path).read_text():
        raise RuntimeError('Recent-parent measurement source changed')


def economic_contract(P):
    """Classifications identify retained assumptions without approving new ones."""
    if float(getattr(P, 'property_tax_lump_sum_transfer', 0.0)) != 0.0:
        raise ValueError('Current reference requires zero property-tax rebates')
    return {
        'preferences': {'classification':'estimated_reference_fixed', 'source':'de_0093 checkpoint',
                        'psi_child':float(P.psi_child), 'renormalize_fertility':False},
        'entry_wealth_income_joint': {'classification':'empirically_normalized_reference_fixed',
                                    'source':'authenticated checkpoint conditional entry arrays',
                                    'joint_rank_coupling':'retained approximation; not newly estimated'},
        'paygo': {'classification':'empirically_normalized', 'payroll_tax':float(P.tau_pay),
                 'rule':'fixed tax; dated benefit must balance actual current population'},
        'property_tax_rebate': {'classification':'externally_fixed', 'transfer':0.0},
        'estate_funding': {'classification':'provisional',
                          'rule':'net liquidation estates fund actual next entrants; residual sink; shortage fails'},
        'baseline_credit': {'classification':'reference_fixed',
                            'rule':'purchase income admitted; final owner collateral floor; renter debt taper'},
        'natural_credit': {'classification':'experimental_not_enabled',
                           'rule':'caller must explicitly enable native solvency mode and record approximation'},
        'population_and_geography': {'classification':'outstanding_for_transition',
                                    'rule':'caller must supply closed endpoint and inherited initial population; no default migration'},
    }


def setup(output, *, contract=CONTRACT, reference=REFERENCE, fixed_reference_price=False):
    """Return current PF runtime plus authenticated reference and observer writer.

    No equilibrium solve, generated source, or household-function replacement.
    Dated closure, native natural credit and acceptance remain caller-owned.
    """
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    if any(name.startswith('intergen_eqscale_seq_optimized') for name in sys.modules):
        raise RuntimeError('Fresh interpreter required: model already imported')
    if sha(contract)!=CONTRACT_SHA: raise RuntimeError('Reference contract differs')
    c=read(contract)
    # Authenticate all declared reference tools/data before executing any of them.
    for item in list(c['files'].values())+[c['base_contract'],c['objective'],c['source_manifest']]:
        if sha(item['path'])!=item['sha256']: raise RuntimeError('Reference pin differs: '+item['path'])
    source=Path(c['source_root'])
    inventory=read(c['source_manifest']['path'])
    for relative, expected in inventory['files'].items():
        path=(source/relative).resolve()
        if not path.is_relative_to(source.resolve()) or sha(path)!=expected:
            raise RuntimeError('Reference inventory differs: '+relative)
    checkpoint=Path(reference)/'initial_state.pkl.gz'; receipt=read(Path(reference)/'receipt.json')
    if sha(checkpoint)!=CHECKPOINT_SHA or receipt['case_checkpoint_sha256']!=CHECKPOINT_SHA:
        raise RuntimeError('Selected reference checkpoint differs')
    os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=CONTRACT_SHA
    sys.path[:0]=[c['runtime_tools']]
    driver=load('current_reference_driver',c['files']['driver']['path'])
    verified,obj=driver.verify(contract)
    runner=load('current_reference_recovery',c['files']['recovery_runner']['path'])
    pair,ancestor,*_=runner.pair_runtime(read(c['base_contract']['path'])['reference_root'])
    tax=ancestor.load_tax_driver()
    current_reporter=ROOT/'code/model/tools/e5f_calibration_runtime.py'
    verify_current_reporter(current_reporter,c['files']['calibration_runtime']['path'])
    observer_runtime=load('current_reference_observer_runtime',current_reporter)
    # Current main imports are authoritative before any checkpoint is unpickled.
    sys.path[:0]=[str(ROOT/'code/model/tools'),str(ROOT/'code/model')]
    primitive=importlib.import_module('run_e5f_matched_pf_smoke')
    chain,model=primitive.pf.transition.configure_sequential_model()
    require_current_model(model)
    primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution
    primitive.pf.transition.calendar.model=model
    with gzip.open(checkpoint,'rb') as stream: selected=pickle.load(stream)
    P=copy.deepcopy(selected['parameters'])
    model.configure_current_household_contract(P)
    P.native_exact_inherited_distribution = True
    P.native_inherited_distribution_evidence_dir = str(output/'inherited_state_failures')
    np.testing.assert_array_equal(model.make_grid(P),selected['b_grid'])
    np.testing.assert_array_equal(P.fixed_reference_entry_grid,selected['b_grid'])
    # Pure observer definitions stay pinned to the actual calibration contract.
    frozen_tools=source/'code/model/tools'
    fert=load('current_pinned_fertility_observer',frozen_tools/'e5f_initial_fertility_observer.py')
    housing=load('current_pinned_housing_observer',frozen_tools/'e5f_initial_housing_observer.py')
    current_recent=ROOT/'code/model/tools/e5f_recent_parent_flow_observer.py'
    verify_current_recent_observer(current_recent,frozen_tools/'e5f_recent_parent_flow_observer.py')
    recent=load('current_pinned_recent_observer',current_recent)
    def observe_recent_with_retained_tail(*args, **kwargs):
        kwargs['diagnostic_allow_retained_dead_tail']=True
        return recent.observe_recent_parent_flow(*args, **kwargs)
    accounting=load('current_pinned_purchase_audit',Path(c['runtime_tools'])/'e5f_earnings_wealth_contract.py')
    estate=EstateAuditContract(load('current_pinned_stationary_estate',c['files']['estate_audit']['path']))
    paygo=importlib.import_module('e5f_stationary_paygo')
    audit=importlib.import_module('run_e5f_independent_numerical_audit')
    require_current_model(model)
    rt=dict(model=model,primitive=primitive,chain=chain,accounting=accounting,audit=audit,
        solve_balanced_initial_equilibrium=paygo.solve_balanced_initial_equilibrium,
        certify_initial_pension=paygo.certify_initial_pension,
        observe_initial_fertility=fert.observe_initial_fertility,
        observe_initial_housing_wealth=housing.observe_initial_housing_wealth,
        observe_recent_parent_flow=observe_recent_with_retained_tail,
        SNAPSHOT=recent.SNAPSHOT,AGE_PROJECTION=recent.AGE_PROJECTION)
    require_current_runtime(rt)
    objective=observer_runtime.StationaryObjective(c,obj,selected,tax,ancestor,rt,runner.adapter,estate)
    # Fixed parameters from de_0093, not constructors and no fertility normalization.
    def fixed_solve(point, case_output, deadline):
        if time.time()>=deadline: raise TimeoutError('Baseline replay deadline reached')
        parameters=copy.deepcopy(P)
        payroll,fiscal_rule=runner.adapter.pension_tax_from_demographics(parameters)
        if abs(payroll-float(parameters.tau_pay))>1e-14:
            raise RuntimeError('Inherited payroll tax does not match adopted stationary rule')
        start=time.monotonic()
        if fixed_reference_price:
            price=np.asarray(selected['solution'].p_eq).copy()
            sol=model.solve_markov_income_at_prices(price,parameters,selected['b_grid'],
                                                   verbose=False,fast_stats=False)
            fiscal=paygo.certify_initial_pension(sol.g,parameters,
                marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        else:
            sol,parameters,price,fiscal=paygo.solve_balanced_initial_equilibrium(
                model=model,parameters=parameters,b_grid=selected['b_grid'],
                initial_prices=selected['solution'].p_eq,payroll_tax=payroll,
                marginal_tolerance=1e-9,fiscal_tolerance=1e-6,warm_price_state={})
        if time.time()>=deadline: raise TimeoutError('Baseline replay exceeded deadline')
        seconds=time.monotonic()-start
        runner.adapter.verify_pension_ratio(fiscal)
        completed=float(chain.extract_moments(sol,parameters)['tfr'])
        normalization=dict(status='fixed_reference_psi_no_normalization',psi_child=float(P.psi_child),
            target=2.1,completed_fertility=completed,absolute_gap=abs(completed-2.1),
            stationary_solves=1,stationary_solve_seconds=seconds)
        renewal=dict(status='reference_conditional_stationary_replay',
            relative_replacement_gap=float(sol.adult_entry_stationary_relative_gap))
        return (sol,parameters,price,seconds,normalization),fiscal_rule,renewal
    objective.normalize=fixed_solve
    source_receipt={}
    for name,module in list(sys.modules.items()):
        path=getattr(module,'__file__',None)
        if path and Path(path).suffix=='.py' and Path(path).resolve().is_relative_to(ROOT/'code'):
            source_receipt[str(Path(path).resolve().relative_to(ROOT))]=sha(path)
    evidence=dict(status='prepared_no_solve',reference_contract_sha256=CONTRACT_SHA,
        reference_checkpoint_sha256=CHECKPOINT_SHA,reference_source_manifest_sha256=c['source_manifest']['sha256'],
        current_source_files=source_receipt,generated_accounting_installers_used=False,
        native_flags={k:v for k,v in vars(P).items() if k.startswith('native_')},
        fixed_reference_price=bool(fixed_reference_price),
        economic_object_classification=economic_contract(P),
        target_rows=len(obj['target_rows']),reference_parameter_rows=len(list(csv.DictReader((Path(reference)/'parameters.csv').open()))),
        limitations=['Native/reference numerical replay required','No dated fiscal or terminal closure certification','No natural-credit adapter installed'])
    write(output/'preparation.json',evidence)
    for relative, expected in source_receipt.items():
        destination=output/'source_snapshot'/relative
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/relative,destination)
        if sha(destination)!=expected:
            raise RuntimeError('Source changed during snapshot: '+relative)
    return dict(parameters=P,selected=selected,reference_receipt=receipt,reference_path=Path(reference),
        runtime=rt,objective=objective,objective_definition=obj,tax=tax,preparation=evidence)


def replay(prepared,output,*,deadline_epoch,graphs=True):
    """One fixed-parameter stationary replay with unchanged full numerical gates."""
    rt=prepared['runtime']; reference=prepared['reference_receipt']
    require_current_runtime(rt)
    verify_current_sources(prepared['preparation'])
    result=prepared['objective'].evaluate_point(tax=prepared['tax'],
        objective=prepared['objective_definition'],selected=prepared['selected'],runtime=rt,
        point=reference['point'],output=Path(output),deadline_epoch=deadline_epoch,graphs=graphs)
    differences={}
    for filename,key,value in [('target_fit.csv','moment','model'),('parameters.csv','parameter','estimate')]:
        old={r[key]:float(r[value]) for r in csv.DictReader((prepared['reference_path']/filename).open())}
        new={r[key]:float(r[value]) for r in csv.DictReader((Path(output)/filename).open())}
        if set(old)!=set(new): raise RuntimeError('Complete reference table rows differ: '+filename)
        differences[filename]={k:new[k]-old[k] for k in old}
    tolerances={'target_fit.csv':1e-5,'parameters.csv':1e-12}
    comparison=dict(differences=differences,tolerances=tolerances,
        passed=all(abs(v)<=tolerances[name] for name,table in differences.items() for v in table.values()))
    with gzip.open(Path(output)/'initial_state.pkl.gz','rb') as stream:
        candidate=pickle.load(stream)
    array_comparison=compare_arrays(prepared['selected'],candidate)
    write(Path(output)/'native_reference_arrays.json',array_comparison)
    verify_current_sources(prepared['preparation'])
    write(Path(output)/'native_reference_comparison.json',comparison)
    # Frozen reporting labels referred to estimation; this replay changes no parameters.
    parameter_path=Path(output)/'parameters.csv'
    rows=list(csv.DictReader(parameter_path.open()))
    for row in rows:
        if row['parameter'] in reference['point'] or row['parameter']=='psi_child':
            row['status']='fixed reference value; no recalibration or fertility normalization'
    with parameter_path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    result.update(status='native_reference_tables_pass_array_review_required' if comparison['passed'] else 'native_reference_replay_failed',
        native_source_identity=prepared['preparation'],fixed_reference_parameters=True,
        economic_changes_relative_to_de0093=[])
    write(Path(output)/'receipt.json',result)
    if not comparison['passed']: raise RuntimeError('Native reference numeric replay failed')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--fixed-reference-price',action='store_true')
    parser.add_argument('--seconds',type=int,default=600)
    args=parser.parse_args()
    if not 1<=args.seconds<=1800: raise ValueError('Explicit bounded 1–1800 second budget required')
    prepared=setup(args.output/'preparation',fixed_reference_price=args.fixed_reference_price)
    if args.replay:
        def timed_out(signum,frame):
            raise TimeoutError('Native replay hard wall-clock budget exhausted')
        prior=signal.signal(signal.SIGALRM,timed_out)
        signal.setitimer(signal.ITIMER_REAL,args.seconds)
        try:
            replay(prepared,args.output/'case',deadline_epoch=time.time()+args.seconds)
        finally:
            signal.setitimer(signal.ITIMER_REAL,0)
            signal.signal(signal.SIGALRM,prior)


if __name__=='__main__': main()
