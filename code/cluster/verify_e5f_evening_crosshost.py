#!/usr/bin/env python3
"""One frozen-source fixed-price replay on Torch; diagnose cross-host equality.

Never changes calibration tolerances or launches/restarts any search. The plan
must pin this wrapper, both contracts, and both original/current checkpoints.
An external parent must enforce the 900-second cap inside the global deadline.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);sys.modules[name]=result;spec.loader.exec_module(result)
    return result


def verify(plan):
    if plan.get('execution_authorized') is not True:raise ValueError('Explicit plan authorization required')
    if plan['runner_sha256']!=sha(__file__):raise ValueError('Runner changed')
    for key in ('original_contract','evening_contract','original_checkpoint','native_checkpoint'):
        item=plan[key]
        if sha(item['path'])!=item['sha256']:raise ValueError('Pinned input changed: '+key)
    if plan['maximum_model_solves']!=1 or not 0<plan['seconds_per_case']<=900:raise ValueError('Only one bounded model solve allowed')
    parent=plan['parent_supervision']
    if (parent['deadline_owner']!='parent' or int(parent['pid'])!=os.getppid()
            or parent['case_deadline_epoch']>plan['global_end_epoch']
            or parent['case_deadline_epoch']-parent['start_epoch']>plan['seconds_per_case']
            or not parent['start_epoch']<=time.time()<parent['case_deadline_epoch']):
        raise ValueError('Parent ownership/deadline mismatch')
    os.kill(parent['pid'],0)


def solution_arrays(sol):
    import numpy as np
    return {key:value for key,value in vars(sol).items()
            if not key.startswith('_') and isinstance(value,np.ndarray)}


def difference(left,right,mass):
    import numpy as np
    if left.shape!=right.shape:return dict(status='shape_mismatch',exact=False)
    finite=bool(np.isfinite(left).all() and np.isfinite(right).all())
    if not finite:return dict(status='nonfinite',exact=False)
    delta=np.abs(left.astype(float)-right.astype(float))
    row=dict(status='compared',exact=bool(np.array_equal(left,right)),max_abs=float(delta.max(initial=0)),l1=float(delta.sum()))
    if left.shape==mass.shape:
        row['occupied_max_abs']=float(delta[mass>1e-12].max(initial=0))
        row['mass_weighted_abs']=float((delta*mass).sum()/mass.sum())
    return row


def compare_three(mac,old_torch,native_torch):
    import numpy as np
    objects={key:solution_arrays(value) for key,value in
             [('mac_original',mac),('torch_frozen',old_torch),('torch_native',native_torch)]}
    names=set.union(*(set(x) for x in objects.values()));mass=np.asarray(mac.g)
    rows=[];all_same_host=True;discrete_same=True;platform_different=False
    for name in sorted(names):
        for first,second in [('mac_original','torch_frozen'),('mac_original','torch_native'),('torch_frozen','torch_native')]:
            if name not in objects[first] or name not in objects[second]:
                comparison=dict(status='missing_array',exact=False)
            else:comparison=difference(objects[first][name],objects[second][name],mass)
            rows.append(dict(array=name,first=first,second=second,**comparison))
            if first=='torch_frozen':
                all_same_host &= comparison['exact']
                if name=='tenure_choice' or (name in objects[first] and objects[first][name].dtype.kind in 'biu'):
                    discrete_same &= comparison['exact']
            elif second=='torch_frozen':platform_different |= not comparison['exact']
    return rows,dict(array_count=len(names),same_host_all_solution_arrays_exact=bool(all_same_host),
        same_host_discrete_choices_exact=bool(discrete_same),mac_to_frozen_torch_difference=bool(platform_different),
        platform_explanation_demonstrated=bool(all_same_host and platform_different),
        scientific_gates_or_tolerances_changed=False)


def run(plan_path,output):
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
    plan=read(plan_path);verify(plan)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    stage='authentication';solve_count=0
    def progress(name):
        nonlocal stage
        stage=name;write(out/'heartbeat.json',dict(stage=stage,epoch=time.time(),model_solves=solve_count))
    try:
        import numpy as np
        import gzip
        import pickle
        cpath=Path(plan['original_contract']['path']);c=read(cpath)
        if any(name.startswith('intergen_eqscale_seq_optimized') for name in sys.modules):raise RuntimeError('Fresh interpreter required')
        sys.path.insert(0,c['runtime_tools'])
        os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=plan['original_contract']['sha256']
        driver_pin=c['files']['driver']
        if sha(driver_pin['path'])!=driver_pin['sha256']:raise RuntimeError('Original driver changed')
        driver=module('samehost_original_driver',driver_pin['path']);verified,obj=driver.verify(cpath)
        runtime_pin=c['files']['calibration_runtime']
        if sha(runtime_pin['path'])!=runtime_pin['sha256']:raise RuntimeError('Original runtime changed')
        runtime=module('samehost_original_runtime',runtime_pin['path'])
        with (Path(plan['original_checkpoint']['path']).parent/'parameters.csv').open(newline='') as stream:
            reference_rows=list(csv.DictReader(stream))
        if len(reference_rows)!=31:raise RuntimeError('Expected complete31parameter reference')
        reference={r['parameter']:float(r['estimate']) for r in reference_rows}
        point={key:reference[key] for key in runtime.FREE}
        progress('frozen_source_setup_zero_solve')
        (out/'preparation').mkdir()
        objective,tax,seed,rt,runner,native_type,evidence,binding=runtime.setup(c,obj,point,out/'preparation')
        model=rt['model'];expected=Path(c['source_root'])/'code/model/intergen_eqscale_seq_optimized/solver.py'
        if Path(model.__file__).resolve()!=expected.resolve():raise RuntimeError('Wrong frozen solver')
        # The original runtime's load_seed is deliberately restricted to its
        # ancestral seed. These two distinct replay packets were authenticated
        # by the outer plan before setup; load them directly after frozen model
        # imports, as the reviewed native Torch runtime does for de_0093.
        with gzip.open(plan['original_checkpoint']['path'],'rb') as stream:
            original=pickle.load(stream)
        with gzip.open(plan['native_checkpoint']['path'],'rb') as stream:
            native=pickle.load(stream)
        P=copy.deepcopy(original['parameters']);before=copy.deepcopy(P)
        grid=original['b_grid'];price=np.asarray(original['solution'].p_eq).copy()
        np.testing.assert_array_equal(seed['b_grid'],grid)
        progress('one_frozen_fixed_price_fixed_psi_solve');solve_count+=1;start=time.monotonic()
        sol=model.solve_markov_income_at_prices(price,P,grid,verbose=False,fast_stats=False)
        seconds=time.monotonic()-start
        with gzip.open(out/'frozen_torch_solution.pkl.gz','wb',compresslevel=1) as stream:
            pickle.dump(dict(parameters=P,b_grid=grid,solution=sol,source_contract_sha256=plan['original_contract']['sha256']),stream,protocol=5)
        for key,value in vars(before).items():
            if key.startswith('_') or key=='eq_iter':continue
            try:equal=np.array_equal(value,getattr(P,key),equal_nan=True)
            except (TypeError,ValueError):equal=repr(value)==repr(getattr(P,key))
            if not equal:raise RuntimeError('Fixed parameter changed: '+key)
        actual=dict(tax.actual_parameters(P),delta_alpha_jump=P.delta_alpha_jump,child_benefit_curvature=P.child_benefit_curvature,
            psi_child=P.psi_child,child_benefit_CRRA_coefficient=(1-P.child_benefit_curvature)*P.psi_child,
            theta1=P.theta1,sigma=P.sigma,alpha_cons=P.alpha_cons,delta_alpha=P.delta_alpha,h_P=P.hbar_first_child_jump,
            utility_reference_rent=P.utility_reference_rent,tenure_choice_kappa=P.tenure_choice_kappa,
            q_annual=(1+P.q)**(1/P.period_years)-1,financed_share=P.phi[0],housing_supply_elasticity=P.xi_supply[0],
            payroll_tax=P.tau_pay,pension_period=P.pension,annual_depreciation=objective.ancestor.ANNUAL_DEP,
            period_depreciation=P.delta,annual_property_tax=objective.ancestor.ANNUAL_PROPERTY_TAX,
            period_property_tax=P.tau_H,selling_cost=P.psi,rental_cap=P.hR_max,wealth_grid_nodes=len(grid),income_states=len(P.z_grid))
        if not set(reference)<=set(actual):raise RuntimeError('Incomplete31parameter comparison')
        # Compare all31 reporter parameters in addition to the full primitive object.
        param_comparison={name:dict(original=value,current=float(actual[name]),exact=value==float(actual[name]))
                          for name,value in reference.items()}
        write(out/'parameter_comparison.json',param_comparison)
        if not all(row['exact'] for row in param_comparison.values()):raise RuntimeError('Reported fixed parameter differs')
        progress('three_way_saved_solution_comparison')
        rows,summary=compare_three(original['solution'],sol,native['solution'])
        fields=list(dict.fromkeys(key for row in rows for key in row))
        with (out/'array_comparison.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        summary.update(status='diagnosis_complete_not_calibration_approval',model_solves=solve_count,
            solve_seconds=seconds,price=price.tolist(),psi_child=float(P.psi_child),
            frozen_model_source_path=str(expected),frozen_model_source_sha256=sha(expected),
            original_checkpoint_sha256=plan['original_checkpoint']['sha256'],native_checkpoint_sha256=plan['native_checkpoint']['sha256'],
            plan_sha256=sha(plan_path),all_parameter_object_primitives_unchanged=True,
            scope='One same-host original-source solution comparison; no normalization, target rescore, or acceptance gate changes')
        verify(plan);write(out/'summary.json',summary)
        return summary
    except BaseException as exc:
        write(out/'failure.json',dict(stage=stage,model_solves=solve_count,error_type=type(exc).__name__,error=str(exc),audit=getattr(exc,'audit',None)))
        raise


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    run(args.plan,args.output)
if __name__=='__main__':main()
