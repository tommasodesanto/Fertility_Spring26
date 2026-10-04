"""Bound-wide Estate-A exploration; no optimizer or economic-code changes."""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys, time
from pathlib import Path
import numpy as np

DEPLOY=Path('/work/deployment')
PARENT=Path('/work/parent_v2')
MODEL=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/experiments/birth_count_choice')
sys.path.insert(0,str(MODEL))
from model import estate_contract as estate
from model.inputs import load_inputs
from model.engine.shared import InfeasibleThetaError
from cluster_calibrate import ANCHOR, check_native

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
    p=Path(p);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n');tmp.replace(p)
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    a=argparse.ArgumentParser();a.add_argument('--mode',choices=['smoke','preflight','production'],required=True)
    a.add_argument('--task',type=int,required=True);a.add_argument('--out',type=Path,required=True)
    a.add_argument('--deadline-epoch',type=float,required=True);a.add_argument('--plan-sha256',required=True)
    z=a.parse_args();require(0<=z.task<16,'Task index outside 0..15')
    require(not z.out.exists(),'Refusing existing run output');z.out.mkdir(parents=True);out=z.out.resolve()
    require(sha(DEPLOY/'control/plan.json')==z.plan_sha256,'Plan SHA drift')
    plan=json.loads((DEPLOY/'control/plan.json').read_text())
    require(plan['stage']=='global_exploration_v1' and len(plan['points'])==64 and plan['task_count']==16 and plan['cases_per_task']==4,'Plan cardinality drift')
    require(sha(PARENT/'inventory.json')==plan['parent_inventory_sha256'] and sha(PARENT/'control/starts.json')==plan['parent_starts_sha256'],'Parent source/starts drift')
    targets,targetpin,weightpin,bounds=estate.contract()
    require(targets==plan['target_contract'] and targetpin==plan['target_fingerprint'] and weightpin==plan['weight_fingerprint'],'Target/weight contract drift')
    require({k:list(v) for k,v in bounds.items()}==plan['bounds'],'Bounds drift')
    anchor=json.loads(ANCHOR.read_text());require(anchor['status']=='selected_numerically_verified' and anchor['chain']==13 and anchor['arm']=='alternative','Anchor drift')
    start=time.time();deadline=min(z.deadline_epoch,start+plan['task_wall_seconds'])
    require(deadline-start>150,'Insufficient run budget')
    write(out/'search_contract.json',dict(stage=plan['stage'],mode=z.mode,task=z.task,plan_sha256=z.plan_sha256,
        target_fingerprint=targetpin,weight_fingerprint=weightpin,task_wall_seconds=plan['task_wall_seconds'],case_budget_seconds=plan['case_budget_seconds'],
        first_case_index=4*z.task,case_count=2 if z.mode=='smoke' else 1 if z.mode=='preflight' else 4,
        search_evaluator='mock_execution_only' if z.mode=='smoke' else 'full_native',max_lifecycle_per_full_GE=32))
    cases=[];best=None
    def beat(status,**kw):write(out/'heartbeat.json',dict(epoch=time.time(),status=status,task=z.task,completed_cases=len(cases),**kw))
    beat('initialized')
    if z.mode!='smoke':
        P,grid=load_inputs(parameters=plan['incumbent_control']['parameters'])
        evaluate=estate.make_evaluator(out,'floor_s0',P,grid,deadline,anchor['selected']['price'],birth_cap=1,
            target_fingerprint=targetpin,weight_fingerprint=weightpin)
    indexes=[0,1] if z.mode=='smoke' else [-1] if z.mode=='preflight' else list(range(4*z.task,4*z.task+4))
    for ix in indexes:
        if time.time()>=deadline-120:break
        point=plan['incumbent_control']['parameters'] if ix<0 else plan['points'][ix]
        require(set(point)==set(bounds) and all(math.isfinite(v) and bounds[k][0]<=v<=bounds[k][1] for k,v in point.items()),'Point outside pinned bounds')
        label=f'case_{ix:03d}' if ix>=0 else 'incumbent_control'
        beat('running_case',case_index=ix,label=label)
        if z.mode=='smoke':
            result=dict(status='passed',residual=([0.1+0.01*ix]*10),lifecycle_solves=0,mock=True)
        else:
            case_end=min(deadline-120,time.time()+plan['case_budget_seconds'])
            try:result=evaluate(label,point,case_end)
            except InfeasibleThetaError as exc:
                result=dict(status='inadmissible_numerical',reason=str(exc),error_type=type(exc).__name__,
                            rejection_kind='native_infeasible_theta_gate',lifecycle_solves=None)
            except RuntimeError as exc:
                if str(exc)=='native GE acceptance failed: uncomputed_bounded_budget':
                    result=dict(status='budget_exhausted',reason=str(exc),error_type=type(exc).__name__,
                                rejection_kind='native_solve_cap',lifecycle_solves=None)
                else:raise
        status=result.get('status')
        row=dict(index=ix,label=label,parameters=point,status=status,case_result=result)
        if status=='passed':
            rr=np.asarray(result['residual'],dtype=float)
            require(rr.shape==(10,) and np.isfinite(rr).all(),'Residual shape/finite drift')
            loss=float(rr@rr)
            if z.mode!='smoke':
                fits,params,plots,native_loss=check_native(result,point,targets,plan['bounds'])
                require(abs(loss-native_loss)<1e-8 and result['target_fingerprint']==targetpin and result['weight_fingerprint']==weightpin,'Native fit drift')
            row['loss']=loss;row['valid_loss']=True
            if best is None or loss<best['loss']:best=row
        elif status=='inadmissible_numerical':
            row.update(valid_loss=False,rejection_kind=result.get('rejection_kind','native_numerical'),reason=result.get('reason',''),error_type=result.get('error_type','inadmissible_numerical'))
        elif status=='budget_exhausted':
            row.update(valid_loss=False,rejection_kind=result.get('rejection_kind','case_time_budget'),reason=result.get('reason',''),error_type=result.get('error_type','budget_exhausted'))
        else:raise RuntimeError('Unexpected native evaluator status: '+str(status))
        cases.append(row)
        write(out/'cases.json',cases);write(out/'latest_completed.json',dict(case=row,completed_cases=len(cases),epoch=time.time()))
        write(out/'best_so_far.json',dict(status='provisional_until_fresh_native_selected_repeat',best=best,completed_cases=len(cases)))
        beat('case_completed',case_index=ix,best_loss=best['loss'] if best else None)
    completed=len(cases)==len(indexes)
    summary=dict(status='smoke_loop_passed_zero_solves' if z.mode=='smoke' and completed else
                 'native_preflight_passed' if z.mode=='preflight' and completed and best is not None else
                 'exploration_completed_provisional' if z.mode=='production' and completed else 'stage_budget_exhausted',
                 mode=z.mode,task=z.task,completed_cases=len(cases),expected_cases=len(indexes),best=best,
                 plan_sha256=z.plan_sha256,target_fingerprint=targetpin,weight_fingerprint=weightpin,
                 no_auto_extension=True,elapsed_seconds=time.time()-start)
    write(out/'completed.json',summary);beat('completed',best_loss=best['loss'] if best else None)
    if z.mode in ('smoke','preflight'): require(completed and best is not None,'Smoke/preflight incomplete or failed')
if __name__=='__main__': main()
