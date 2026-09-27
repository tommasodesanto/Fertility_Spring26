#!/usr/bin/env python3
"""One authenticated old-floor objective replay through the optional-off new source.

Torch-only regression, not the new utility calibration. Original adapters,
objective, fiscal closure, normalizer start/step and exact-loop evaluator remain.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import sys
import time

BASE = Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
OLD = BASE/'nightpair_20260925_v1'
NEW = BASE/'calibration_code_integration_20260927_v2'
REFERENCE = OLD/'results/run_001/oasi_087510/worker_05/point_04'
PAIR_PIN = '6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1'
NEW_PIN = '399abb6e9eab0d447dca627f94e3de4e6a8006d8920d02241fadac21f6c6ebae'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda:stream.read(1<<20),b''): h.update(part)
    return h.hexdigest()


def write(path, data):
    Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); sys.modules[name]=module
    spec.loader.exec_module(module); return module


def array_comparison(old,new):
    """Exact direct ndarray fields, including nonfinite positions; no timings."""
    import numpy as np
    result=[]
    for group in ('solution','shared','evaluation'):
        left,right=vars(old[group]),vars(new[group])
        names=sorted(k for k,v in left.items() if isinstance(v,np.ndarray))
        for name in names:
            a,b=left[name],right.get(name)
            equal=isinstance(b,np.ndarray) and a.dtype==b.dtype and a.shape==b.shape and np.array_equal(a,b,equal_nan=True)
            row=dict(group=group,array=name,shape=list(a.shape),dtype=str(a.dtype),exact_equal=bool(equal))
            if isinstance(b,np.ndarray) and a.shape==b.shape and np.issubdtype(a.dtype,np.number):
                mask=np.isfinite(a)&np.isfinite(b)
                row['maximum_finite_absolute_gap']=float(np.max(np.abs(a[mask]-b[mask]))) if mask.any() else 0.
            result.append(row)
    for name in ('b_grid','stationary_g_pre'):
        a,b=old[name],new[name]
        result.append(dict(group='packet',array=name,shape=list(a.shape),dtype=str(a.dtype),exact_equal=bool(a.dtype==b.dtype and np.array_equal(a,b,equal_nan=True))))
    return result


def table_comparison(old,new,key):
    with Path(old).open() as stream: left=list(csv.DictReader(stream))
    with Path(new).open() as stream: right=list(csv.DictReader(stream))
    assert [r[key] for r in left]==[r[key] for r in right], 'Changed table row identities'
    records=[]
    for a,b in zip(left,right):
        assert set(a)==set(b), 'Changed table columns'
        for column in a:
            records.append(dict(row=a[key],column=column,reference=a[column],replay=b[column],exact_equal=a[column]==b[column]))
    return records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=('preflight','run'))
    parser.add_argument('output',type=Path)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    os.environ['EXPECTED_PAIR_LOCK_SHA256']=PAIR_PIN
    contract_path=NEW/'launch_v1/contract.json'
    assert sha(contract_path)==NEW_PIN
    contract=json.loads(contract_path.read_text())
    manifest=contract['source_manifest']; assert sha(manifest['path'])==manifest['sha256']
    source=Path(contract['source_root']).resolve()
    inventory=json.loads(Path(manifest['path']).read_text())
    for relative,pin in inventory['files'].items():
        path=(source/relative).resolve()
        assert path.is_relative_to(source) and sha(path)==pin,relative
    pair=load('default_replay_original_pair',OLD/'run_pair.py')
    old=pair.load_ancestor(); lock,objective,_=pair.read_contract(old)
    # Source substitution occurs only after complete original-package validation.
    pair.SOURCE=source
    pair.configure(old,'oasi_087510',lock)
    tax,plan,selected,objective=pair.prepare(old,lock)
    tax.SOURCE_SHA=manifest['sha256']
    runtime=old.setup_runtime(tax,plan,selected,args.output)
    assert Path(runtime['model'].__file__).resolve()==source/'code/model/intergen_eqscale_seq_optimized/solver.py'
    reference=json.loads((REFERENCE/'receipt.json').read_text())
    point=reference['point']; P=old.apply_point(selected,point)
    assert getattr(P,'child_benefit_curvature',0.)==0.
    assert not getattr(P,'compensated_child_housing_shares',False)
    assert P.hbar_first_child_jump==point['h_P'] and P.child_room_floor
    assert old.TAX==reference['payroll_tax']==0.08751017424959717
    with gzip.open(REFERENCE/'initial_state.pkl.gz','rb') as stream: retained=pickle.load(stream)
    assert sha(REFERENCE/'initial_state.pkl.gz')==reference['case_checkpoint_sha256']
    repeated_shared=runtime['model'].precompute_shared(retained['parameters'],retained['b_grid'])
    import numpy as np
    preflight_arrays=[]
    for key,value in vars(retained['shared']).items():
        if isinstance(value,np.ndarray):
            fresh=getattr(repeated_shared,key)
            preflight_arrays.append(dict(array=key,exact_equal=bool(value.dtype==fresh.dtype and np.array_equal(value,fresh,equal_nan=True))))
    assert preflight_arrays and all(row['exact_equal'] for row in preflight_arrays)
    provenance=dict(old_pair_lock_sha256=PAIR_PIN,new_contract_sha256=NEW_PIN,
                    new_source_inventory_sha256=manifest['sha256'],
                    reference_receipt_sha256=sha(REFERENCE/'receipt.json'),
                    reference_checkpoint_sha256=sha(REFERENCE/'initial_state.pkl.gz'),
                    source_root=str(source),reference=str(REFERENCE),point=point,
                    child_benefit_curvature=0.,compensated_child_housing_shares=False,
                    warm_price=False,normalization_start=old.START_PSI,
                    normalization_step=.25,payroll_tax=old.TAX,
                    floor=point['h_P'],preflight_arrays=preflight_arrays)
    write(args.output/'preflight.json',provenance)
    print(json.dumps(dict(status='default_path_preflight_passed',arrays=len(preflight_arrays))),flush=True)
    if args.stage=='preflight': return
    started=time.monotonic()
    receipt=old.evaluate_point(tax=tax,objective=objective,selected=selected,runtime=runtime,
                point=point,output=args.output/'case',deadline_epoch=time.time()+1650,graphs=True)
    elapsed=time.monotonic()-started
    with gzip.open(args.output/'case/initial_state.pkl.gz','rb') as stream: replay=pickle.load(stream)
    arrays=array_comparison(retained,replay)
    fits=table_comparison(REFERENCE/'target_fit.csv',args.output/'case/target_fit.csv','moment')
    params=table_comparison(REFERENCE/'parameters.csv',args.output/'case/parameters.csv','parameter')
    scalars={key:dict(reference=reference[key],replay=receipt[key],exact_equal=reference[key]==receipt[key])
             for key in ('loss','price','market_residual','objective_stationary_solves')}
    for key in ('completed_fertility','psi_child','absolute_gap','stationary_solves'):
        a,b=reference['normalization'][key],receipt['normalization'][key]
        scalars['normalization.'+key]=dict(reference=a,replay=b,exact_equal=a==b)
    passed=all(r['exact_equal'] for r in arrays+fits+params) and all(r['exact_equal'] for r in scalars.values())
    result=dict(status='exact_default_path_replay_passed' if passed else 'default_path_replay_difference',
                elapsed_seconds=elapsed,original_stationary_seconds=reference['objective_stationary_solve_seconds'],
                replay_stationary_seconds=receipt['objective_stationary_solve_seconds'],
                original_solve_count=reference['objective_stationary_solves'],
                replay_solve_count=receipt['objective_stationary_solves'],
                target_rows=13,parameter_rows=len({r['row'] for r in params}),
                diagnostics=len(list((args.output/'case/standard_diagnostics').glob('*.png'))),
                array_comparison=arrays,scalar_comparison=scalars,
                provenance=provenance)
    old.table(args.output/'target_comparison.csv',fits)
    old.table(args.output/'parameter_comparison.csv',params)
    write(args.output/'comparison.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('array_comparison','scalar_comparison','provenance')}),flush=True)
    assert passed,'Default-path numerical regression: retained differences in comparison.json'


if __name__=='__main__': main()
