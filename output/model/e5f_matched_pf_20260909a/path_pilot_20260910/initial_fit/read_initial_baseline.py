#!/usr/bin/env python3
"""Read-only initial-state diagnostic, never a calibration or model solve.

Validates immutable source and checkpoint hashes before imports/unpickling.
Uses existing observers only. Unimplemented/mismatched rows remain explicit,
with proxy values in a separate column and no invented objective weights.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
HERE=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def read_rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def write_rows(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def validate_files(source_root,contract_path):
    manifest=json.loads((HERE/'manifest.json').read_text())
    if digest(contract_path)!=manifest['contract_sha256']:
        raise ValueError('Historical source/input contract hash changed')
    contract=json.loads(contract_path.read_text())
    rows=read_rows(HERE/'initial_restriction_observer_map.csv')
    if (digest(HERE/'initial_restriction_observer_map.csv')!=manifest['table_sha256']
            or len(rows)!=13 or len({r['restriction'] for r in rows})!=13):
        raise ValueError('Diagnostic observation ledger changed or incomplete')
    if digest(HERE/'inherited_parameters.csv')!=manifest['parameter_source_sha256']:
        raise ValueError('Inherited parameter receipt changed')
    if any(r['weight'] for r in rows):raise ValueError('This readout cannot set weights')
    failures=[]
    for relative,expected in contract['source_sha256'].items():
        p=(source_root/relative).resolve()
        if not p.is_relative_to(source_root.resolve()):raise ValueError('Unsafe source path')
        if not p.is_file() or digest(p)!=expected:failures.append(relative)
    if failures:raise ValueError('Source fingerprint mismatch: '+', '.join(failures[:8]))
    if (contract['normalized_checkpoint']!=manifest['normalized_checkpoint'] or
            contract['normalized_checkpoint_sha256']!=manifest['normalized_checkpoint_sha256']):
        raise ValueError('Wrong normalized old state')
    return manifest,contract,rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',required=True,type=Path)
    ap.add_argument('--contract',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--validate-only',action='store_true')
    args=ap.parse_args();started=time.monotonic()
    manifest,c,rows=validate_files(args.source_root,args.contract)
    # The supplied output directory must be fresh: never overwrite a readout.
    args.output.mkdir(parents=True,exist_ok=False)
    validation={'source_files_verified':len(c['source_sha256']),
        'contract_sha256':manifest['contract_sha256'],'target_rows':13,
        'weights_set':False,'checkpoint_loaded':False,'model_solve_count':0}
    if args.validate_only:
        validation['status']='PASS_pure_input_validation_no_model_import'
        (args.output/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
        print(json.dumps(validation));return

    sys.path[:0]=[str(args.source_root/'code/model/tools'),str(args.source_root/'code/model')]
    import numpy as np
    import run_e5f_matched_pf_baseline as baseline
    import run_e5f_transition_calibration as measurement
    from intergen_eqscale_seq_optimized.calibration import extract_moments
    old,_=baseline.load_normalized(c,'sequential')
    validation['checkpoint_loaded']=True
    P=old.parameters
    if float(old.diagnostics['normalization']['target'])!=2.1:
        raise ValueError('Author-selected initial benchmark changed')
    if float(old.parameters.period_years)!=4. or bool(P.joint_nested_choice):
        raise ValueError('Expected retained sequential four-year old state')
    params=read_rows(HERE/'inherited_parameters.csv')
    for row in params:
        name=row['parameter'];expected=float(row['value'])
        if name=='beta_annual':actual=float(P.beta)**.25
        elif name=='psi_child_2007':actual=float(P.psi_child)
        elif name=='psi_child_2023':actual=float(old.psi_path[-1])
        elif name=='psi_child_change_2023':actual=float(old.psi_path[-1]-old.psi_path[0])
        elif name=='housing_supply_elasticity':actual=float(old.supply_rule.elasticity)
        else:actual=float(getattr(P,name))
        if not math.isclose(actual,expected,rel_tol=0,abs_tol=1e-12):
            raise ValueError('Parameter receipt does not match checkpoint: '+name)
        row['readout_role']=('historical shock; excluded from initial search' if name=='psi_child_change_2023'
            else 'inherited input; no estimation performed')
    write_rows(args.output/'parameters.csv',params)

    pf=baseline.pf
    pf.transition.configure_sequential_model()
    pf.calendar.apply_fertility=pf.transition.apply_sequential_fertility
    pf.calendar.advance_calendar_distribution=pf.transition.advance_sequential_calendar_distribution
    counter=pf.calendar.SolveCounter()
    # Use the stationary distribution and stationary supply primitive. The
    # ACS-age-reweighted initial_state has a separately normalized dated supply.
    evaluation=pf.calendar.evaluate_period(old.policy.price,old.stationary_g_pre,
        P,old.b_grid,old.shared,counter,supply_rule=None,supplied_policy=old.policy)
    if counter.total!=0:raise RuntimeError('A readout unexpectedly solved a model')
    cached=extract_moments(old.solution,P)
    period=measurement.period_fertility_diagnostics(evaluation,P)
    housing=measurement.first_birth_housing_response(evaluation,P,old.b_grid,old.shared)
    names=tuple(dict.fromkeys(r['proxy_key'] for r in rows
        if r['proxy_key'] and not r['proxy_key'].startswith('period_') and r['proxy_key']!='tfr'))
    cross=measurement.transition_cross_section_moments(evaluation,P,old.b_grid,old.shared,
        names,housing_increment_override=housing)
    values={**cached,**cross,**period}
    values['tfr']=float(old.diagnostics['normalization']['completed_fertility'])
    outrows=[]
    for row in rows:
        q=dict(row);key=row['proxy_key'];value=values.get(key)
        finite=bool(key and value is not None and np.ndim(value)==0 and math.isfinite(float(value)))
        comparable=row['comparison_status'].startswith('available_')
        q['model_value']=float(value) if finite and comparable else ''
        q['available_proxy_value']=float(value) if finite and not comparable else ''
        q['gap']=float(value)-float(row['target']) if finite and comparable else ''
        q['loss_contribution']=''
        if not finite:q['execution_note']='No existing named observer; left uncomputed'
        elif comparable:q['execution_note']='Existing observer evaluated; stated approximation remains'
        else:q['execution_note']='Proxy only; no data-model gap or objective assigned'
        outrows.append(q)
    write_rows(args.output/'all_13_initial_rows.csv',outrows)
    supplemental={
        'old_completed_fertility_normalization':old.diagnostics['normalization'],
        'cached_terminal_pooled_parity_dist':np.asarray(old.solution.parity_dist).tolist(),
        'cached_terminal_pooled_parity_warning':'Not CPS women40–44 and not a female exposure rate',
        'old_cps_completed_fertility_capped5_diagnostic':1.8566081212581498,
        'old_cps_completed_fertility_uncapped_diagnostic':1.8783838839209661,
        'old_age_wealth_median_existing_observer':float(cached['old_total_wealth_to_annual_income_median_7684']),
        'young_ownership_existing_observer':float(cached['own_rate_2534']),
        'period_fertility_diagnostics':{k:np.asarray(v).tolist() for k,v in period.items()},
        'stationary_mass':float(old.stationary_g_pre.sum()),
        'age_reweighted_2007_mass':float(old.initial_state.g_pre.sum()),
        'stationary_supply_elasticity':np.asarray(P.xi_supply).tolist(),
        'dated_supply_elasticity':float(old.supply_rule.elasticity),
        'stationary_evaluation_market_residual':float(evaluation.relative_market_residual),
        'choice_distribution':'old stationary, before ACS2007 age reweight and announcement',
        'no_new_female_or_capped_room_or_family_group_observer':True}
    (args.output/'supplemental_existing_observers.json').write_text(json.dumps(supplemental,indent=2)+'\n')
    validation.update(status='complete_existing_observer_readout_not_calibration',
        elapsed_seconds=time.monotonic()-started,model_solve_count=counter.total,
        objective_computed=False,production_promoted=False,
        unsupported_rows=[r['restriction'] for r in outrows if r['model_value']==''])
    (args.output/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation))

if __name__=='__main__':main()
