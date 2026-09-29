"""Bounded independent local-search stream; Torch only at the CLI boundary.

No model imports. Residuals are the worker's ten weighted moment gaps. A failed
probe never becomes a zero derivative, a penalty, or an imputed observation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback
import numpy as np

PARAMETERS=('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert',
 'kappa_fert_continuation','theta0','delta_alpha_jump','child_benefit_curvature','tenure_choice_kappa')
RIDGES=(1e-4,1e-2,1.)


def require(ok,message):
    if not ok: raise RuntimeError(message)


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):result.update(block)
    return result.hexdigest()


def key(point):return json.dumps(point,sort_keys=True,separators=(',',':'),allow_nan=False)


def check_point(point,names,bounds):
    require(set(point)==set(names) and len(point)==10,'Exactly ten coordinates required')
    for name in names:
        value=float(point[name]);lo,hi=bounds[name]
        require(math.isfinite(value) and lo<=value<=hi,'Coordinate outside bounds: '+name)


def same_point(left,right,names):
    """The native annual-to-period conversion can move beta by one ULP only."""
    require(set(left)==set(right)==set(names),'Point coordinate set changed')
    for name in names:
        a,b=float(left[name]),float(right[name])
        require(math.isfinite(a) and math.isfinite(b),'Nonfinite point coordinate')
        if name=='beta_annual':
            require(abs(a-b)<=2e-12,'Annual beta changed beyond inherited tolerance')
        else:
            require(a==b,'Proposal coordinate changed: '+name)


def scales(point,names):
    recipe=dict(H0=max(abs(point['H0']),1.),beta_annual=.01,chi=max(abs(point['chi']),.5),
      first_birth_fixed_cost=max(abs(point['first_birth_fixed_cost']),.25),
      kappa_fert=max(abs(point['kappa_fert']),.1),kappa_fert_continuation=max(abs(point['kappa_fert_continuation']),.1),
      theta0=max(abs(point['theta0']),.1),delta_alpha_jump=.05,child_benefit_curvature=.1,
      tenure_choice_kappa=max(abs(point['tenure_choice_kappa']),.005))
    return np.array([recipe[name] for name in names],dtype=float)


def probe_design(center,names,bounds):
    check_point(center,names,bounds);scale=scales(center,names);result=[]
    for i,name in enumerate(names):
        lo,hi=bounds[name];x=center[name];wanted=.05*scale[i]
        up,down=hi-x,x-lo
        delta=min(wanted,up) if up>=wanted or up>=down else -min(wanted,down)
        require(delta!=0 and math.isfinite(delta),'No feasible finite-difference displacement')
        point=dict(center);point[name]=x+delta;check_point(point,names,bounds)
        actual=(point[name]-x)/scale[i]
        require(actual!=0,'Finite-difference displacement rounded to zero')
        result.append(dict(parameter=name,point=point,scaled_step=float(actual)))
    return scale,result


def jacobian(center,probes,design):
    require(len(probes)==len(design)==10,'Incomplete Jacobian must not be used')
    r=np.asarray(center['residuals'],float)
    require(r.shape==(10,) and np.isfinite(r).all(),'Invalid center residuals')
    J=np.column_stack([(np.asarray(probe['residuals'],float)-r)/row['scaled_step'] for probe,row in zip(probes,design)])
    dpsi=np.array([(probe['psi']-center['psi'])/row['scaled_step'] for probe,row in zip(probes,design)])
    require(J.shape==(10,10) and np.isfinite(J).all() and np.isfinite(dpsi).all(),'Nonfinite Jacobian')
    return J,dpsi


def gn_proposals(center,names,bounds,scale,J,dpsi):
    r=np.asarray(center['residuals'],float);u,s,vt=np.linalg.svd(J,full_matrices=False)
    largest=float(s[0]);rank=int(np.sum(s>largest*1e-6)) if largest>0 else 0
    diagnostics=dict(singular_values=s.tolist(),numerical_rank=rank,rank_relative_threshold=1e-6,
       condition=float(s[0]/s[-1]) if s[-1]>0 else None,
       interpretation='Local numerical sensitivity, not statistical or global identification.',
       finite_difference_scaled=.05,finite_difference_step_noise='Not estimated by one-sided probes; final repeat differences reported separately')
    result=[];x=np.array([center['point'][name] for name in names])
    if largest==0:return result,diagnostics
    for ridge in RIDGES:
        step=-(vt.T@((s/(s*s+ridge*largest*largest))*(u.T@r)))
        step/=max(1.,float(np.max(np.abs(step))))
        point={name:float(np.clip(x[i]+scale[i]*step[i],*bounds[name])) for i,name in enumerate(names)}
        actual=np.array([(point[name]-center['point'][name])/scale[i] for i,name in enumerate(names)])
        predicted=r+J@actual
        psi=float(np.clip(center['psi']+dpsi@actual,.5*center['psi'],1.5*center['psi']))
        require(psi>0 and np.isfinite(predicted).all(),'Invalid warm prediction')
        result.append(dict(point=point,initial_psi=psi,ridge=ridge,scaled_step=actual.tolist(),
           predicted_residuals=predicted.tolist(),predicted_loss=float(predicted@predicted)))
    return result,diagnostics


def exploratory(best,names,bounds,rng,index):
    radius=(.25,.5,1.)[index%3]
    direction=rng.normal(size=10);direction/=max(np.max(np.abs(direction)),1e-300)
    scale=scales(best['point'],names)
    point={name:float(np.clip(best['point'][name]+radius*scale[i]*direction[i],*bounds[name])) for i,name in enumerate(names)}
    check_point(point,names,bounds)
    return point


def nearest_psi(point,successes,names):
    scale=scales(point,names)
    closest=min(successes,key=lambda r:sum(((point[n]-r['point'][n])/scale[i])**2 for i,n in enumerate(names)))
    return closest['psi']


def run_process(command,log,deadline,heartbeat,poll=.25):
    """One owned process group, no retry; parent owns hard timeout classification."""
    require(time.time()<deadline,'Case deadline already exhausted')
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[k]='1'
    start=time.monotonic();remaining=deadline-time.time();timed_out=False
    with Path(log).open('xb') as stream:
        process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        def kill():
            try:os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            process.wait(timeout=10)
        try:
            while True:
                code=process.poll();expired=time.time()>=deadline or time.monotonic()-start>=remaining
                if code is None and expired:
                    timed_out=True;kill();code=process.returncode
                heartbeat()
                if code is not None:
                    return dict(returncode=code,timed_out=timed_out,late_completion=expired and not timed_out,
                                elapsed_seconds=time.monotonic()-start,pid=process.pid)
                time.sleep(min(poll,max(.001,remaining-(time.monotonic()-start))))
        finally:
            if process.poll() is None:kill()


def validate_success(result,request,folder,config,synthetic=False):
    require(result['status']=='passed','Worker success is not passed')
    for name in ('candidate_id','lane','config_sha256','source_fingerprint','target_fingerprint'):
        require(result[name]==request[name],'Worker identity mismatch: '+name)
    same_point(result['point'],request['point'],config['parameters'])
    check_point(result['point'],config['parameters'],config['bounds'])
    residuals=np.asarray(result['residuals'],float)
    require(residuals.shape==(10,) and np.isfinite(residuals).all(),'Ten finite scored residuals required')
    loss=float(result['loss']);require(math.isfinite(loss) and math.isclose(loss,float(residuals@residuals),rel_tol=1e-10,abs_tol=1e-10),'Loss/residual inconsistency')
    require(math.isfinite(result['psi']) and result['psi']>0,'Invalid normalized psi')
    require(isinstance(result['model_evaluations'],int) and 1<=result['model_evaluations']<=request['maximum_stationary_solves'],'Stationary solve count outside cap')
    require(math.isfinite(result['elapsed_seconds']) and result['elapsed_seconds']>=0,'Invalid worker duration')
    case=Path(result['case_path']).resolve();require(case.is_relative_to(folder.resolve()),'Case escaped owned output')
    artifacts=result['artifacts'];names=[];seen=set()
    for artifact in artifacts:
        path=Path(artifact['path']).resolve()
        require(path.is_relative_to(case) and str(path) not in seen,'Invalid or duplicate artifact')
        require(sha(path)==artifact['sha256'],'Artifact hash changed: '+str(path))
        seen.add(str(path));names.append(path.name)
    if not synthetic:
        require(names.count('target_fit.csv')==names.count('parameters.csv')==1,'Complete table artifacts required')
        require(sum(name.endswith('.png') for name in names)==17,'All17 standard plot hashes required')
        require(result.get('complete_target_rows')==14 and result.get('complete_parameter_rows')==31 and result.get('standard_plots')==17,
                'Worker did not attest complete14/31/17 export')
        case_receipt=case/'receipt.json';identity=case/'scientific_identity.json';checkpoint=case/'initial_state.pkl.gz'
        require(case_receipt.exists() and identity.exists() and checkpoint.exists(),'Required native receipt/identity/checkpoint missing')
        native=read(case_receipt);scientific=read(identity)
        for item in (native,scientific):
            for name in ('candidate_id','lane','config_sha256','source_fingerprint','target_fingerprint'):
                require(item[name]==request[name],'Case scientific identity mismatch: '+name)
        require(native['case_checkpoint_sha256']==result['checkpoint_sha256']==scientific['checkpoint_sha256'],
                'SUCCESS/native/identity checkpoint disagreement')
        require(math.isclose(float(native['normalization']['psi_child']),float(result['psi']),rel_tol=0,abs_tol=2e-12),
                'SUCCESS normalized psi differs from native receipt')
        require(math.isclose(float(scientific['normalized_psi']),float(result['psi']),rel_tol=0,abs_tol=2e-12),
                'SUCCESS normalized psi differs from scientific identity')
    return result


def validate_failure(failure,request,budget):
    """Authenticate a native failure before controller timeout classification."""
    for name in ('candidate_id','lane','config_sha256','source_fingerprint','target_fingerprint'):
        require(failure[name]==request[name],'Failure identity mismatch: '+name)
    require(failure['status'] in ('inadmissible','fatal','censored'),'Unknown failure classification')
    require(failure.get('authenticated') is True or failure['status']=='fatal',
            'Unauthenticated scientific/numerical rejection')
    if failure['status']=='censored':
        require(failure.get('model_evaluations')==budget['maximum_stationary_solves'],
                'Solve-budget censor requires complete owned ledger')
    return failure


class Stream:
    def __init__(self,config,lane,output,config_sha):
        self.c=config;self.lane=lane;self.spec=config['lanes'][lane];self.output=Path(output)
        self.pin=config_sha;self.names=config['parameters'];self.bounds=config['bounds'];self.b=config['budget']
        require(len(self.names)==10 and set(self.names)==set(PARAMETERS),'Wrong parameter ordering/set')
        require(len(config['scored_moments'])==10 and len(set(config['scored_moments']))==10,'Ten distinct scored moments required')
        require(self.b['max_evaluations']==36,'Evaluation design changed')
        if not config.get('synthetic',False):
            require(self.b==dict(total_seconds=25200,final_reserve_seconds=4200,case_seconds=2100,maximum_stationary_solves=8,max_evaluations=36),'Production budget differs from reviewed design')
        check_point(self.spec['initial_point'],self.names,self.bounds)
        self.output.mkdir(parents=True,exist_ok=False)
        self.start=time.time();self.end=min(self.start+self.b['total_seconds'],config['hard_end_epoch'])
        self.cutoff=self.end-self.b['final_reserve_seconds']
        require(self.start+self.b['case_seconds']<=self.cutoff,'Insufficient time for warm replay and final reserve')
        self.records=[];self.successes=[];self.best=None;self.seen=set();self.consecutive=0;self.stop_reason=None
        self.last_heartbeat=0.;self.active=None;self.selected_checkpoint_audited=False
        write(self.output/'clock.json',dict(start=self.start,end=self.end,search_cutoff=self.cutoff,config_sha256=self.pin))
        write(self.output/'latest_completed.json',dict(status='none_completed'))
        write(self.output/'best_so_far.json',dict(status='none_completed',promotion='disabled'))
        self.heartbeat(force=True)

    def heartbeat(self,force=False):
        now=time.time()
        if force or now-self.last_heartbeat>=15:
            write(self.output/'heartbeat.json',dict(epoch=now,elapsed=now-self.start,end=self.end,search_cutoff=self.cutoff,active=self.active,
               attempts=len(self.records),best_loss=None if self.best is None else self.best['loss'],stop_reason=self.stop_reason))
            self.last_heartbeat=now

    def can_start(self,final=False):
        limit=self.b['max_evaluations'] if final else self.b['max_evaluations']-2
        end=self.end if final else self.cutoff
        allowed=self.stop_reason is None or (final and self.stop_reason=='two_consecutive_censored_or_inadmissible')
        return allowed and len(self.records)<limit and time.time()+self.b['case_seconds']<=end

    def save(self):
        write(self.output/'records.json',self.records)
        if self.records:write(self.output/'latest_completed.json',self.records[-1])
        if self.best:write(self.output/'best_so_far.json',dict(self.best,promotion='disabled'))
        self.heartbeat(force=True)

    def evaluate(self,point,psi,role,final=False,proposal=None):
        if not self.can_start(final):return None
        check_point(point,self.names,self.bounds)
        require(math.isfinite(psi) and psi>0,'Invalid starting psi')
        if not final and key(point) in self.seen:return None
        case_id=f'{self.lane}_{len(self.records):03d}_{role}'
        folder=self.output/case_id;deadline=min(time.time()+self.b['case_seconds'],self.end if final else self.cutoff)
        request=dict(candidate_id=case_id,lane=self.lane,config_sha256=self.pin,point=point,initial_psi=float(psi),role=role,
            deadline_epoch=deadline,case_budget_seconds=self.b['case_seconds'],maximum_stationary_solves=self.b['maximum_stationary_solves'],
            source_fingerprint=self.spec['source_fingerprint'],target_fingerprint=self.spec['target_fingerprint'],scored_moments=self.c['scored_moments'])
        path=self.output/(case_id+'.request.json');write(path,request)
        record=dict(candidate_id=case_id,role=role,point=point,initial_psi=float(psi),status='started',request_path=str(path),proposal=proposal)
        self.records.append(record);self.seen.add(key(point));self.active=case_id;self.save()
        try:
            process=run_process(self.c['worker_command']+['--request',str(path),'--output',str(folder)],self.output/(case_id+'.log'),deadline,self.heartbeat)
            record.update(process)
            failure_path=folder/'FAILURE.json';success_path=folder/'SUCCESS.json'
            if failure_path.exists():
                failure=validate_failure(read(failure_path),request,self.b)
                require(not success_path.exists(),'Both success and failure exist')
                record.update(status=failure['status'],failure=failure)
            elif process['timed_out'] or process['late_completion']:
                record.update(status='censored',failure='owned_process_timeout' if process['timed_out'] else 'late_completion')
            elif process['returncode']==0:
                receipt=success_path
                result=validate_success(read(receipt),request,folder,self.c,self.c.get('synthetic',False))
                result=dict(result,success_path=str(receipt),success_sha256=sha(receipt))
                record.update(status='success',result=result)
                if not final:
                    self.successes.append(result)
                    if self.best is None or result['loss']<self.best['loss']:self.best=result
                self.consecutive=0;self.active=None;self.save();return result
            else:
                failure=validate_failure(read(failure_path),request,self.b)
                record.update(status=failure['status'],failure=failure)
        except Exception as exc:
            record.update(status='fatal',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
        if record['status']=='fatal':self.stop_reason='fatal_integrity_or_scientific_failure'
        else:
            self.consecutive+=1
            if self.consecutive>=2:self.stop_reason='two_consecutive_censored_or_inadmissible'
        self.active=None;self.save();return None

    def authenticate_best(self):
        require(self.best is not None,'No valid best point')
        require(sha(self.best['success_path'])==self.best['success_sha256'],'Selected receipt changed')
        if not self.selected_checkpoint_audited:
            case=Path(self.best['case_path']);native=read(case/'receipt.json');identity=read(case/'scientific_identity.json')
            require(sha(case/'initial_state.pkl.gz')==self.best['checkpoint_sha256']==native['case_checkpoint_sha256']==identity['checkpoint_sha256'],
                    'Selected checkpoint authentication failed')
            self.selected_checkpoint_audited=True

    def run(self):
        self.evaluate(self.spec['initial_point'],self.spec['initial_psi'],'initial_replay')
        if self.best is None:
            self.stop_reason=self.stop_reason or 'initial_replay_unavailable';self.save();return self.finish([])
        for round_index in range(2):
            if not self.can_start():break
            center=dict(self.best);scale,design=probe_design(center['point'],self.names,self.bounds)
            probes=[]
            for i,row in enumerate(design):
                result=self.evaluate(row['point'],center['psi'],f'jac{round_index}_{i}')
                if result is None:break
                probes.append(result)
            info=dict(center=center,coordinate_scales=scale.tolist(),design=design,complete=len(probes)==10,valid_probe_count=len(probes))
            if len(probes)!=10:
                info['reason']='Incomplete derivative round; no imputation or zero columns'
                write(self.output/f'jacobian_{round_index}.json',info);break
            J,dpsi=jacobian(center,probes,design)
            proposals,diagnostics=gn_proposals(center,self.names,self.bounds,scale,J,dpsi)
            info.update(matrix=J.tolist(),psi_derivative=dpsi.tolist(),diagnostics=diagnostics,proposals=proposals)
            write(self.output/f'jacobian_{round_index}.json',info)
            for i,proposal in enumerate(proposals):
                self.evaluate(proposal['point'],proposal['initial_psi'],f'gn{round_index}_{i}',proposal=proposal)
                if self.stop_reason:break
            if self.stop_reason:break
        rng=np.random.default_rng(self.spec['seed'])
        for index in range(7):
            if not self.can_start():break
            point=exploratory(self.best,self.names,self.bounds,rng,index)
            self.evaluate(point,nearest_psi(point,self.successes,self.names),f'explore{index}')
        repeats=[]
        if self.stop_reason in (None,'two_consecutive_censored_or_inadmissible'):
            self.authenticate_best();selected=dict(self.best)
            write(self.output/'selected_before_repeats.json',selected)
            for i in range(2):
                self.authenticate_best()
                result=self.evaluate(selected['point'],selected['psi'],f'final_repeat{i}',final=True)
                if result is None:break
                repeats.append(result)
        return self.finish(repeats)

    def finish(self,repeats):
        screens=[]
        if self.best:
            for result in repeats:
                loss=abs(result['loss']-self.best['loss'])
                residual=float(np.max(np.abs(np.asarray(result['residuals'])-self.best['residuals'])))
                screens.append(dict(candidate_id=result['candidate_id'],loss_difference=loss,max_residual_difference=residual,passed=loss<=.05 and residual<=.01))
        if len(repeats)==2:
            loss=abs(repeats[0]['loss']-repeats[1]['loss'])
            residual=float(np.max(np.abs(np.asarray(repeats[0]['residuals'])-repeats[1]['residuals'])))
            screens.append(dict(candidate_id='repeat_pair',loss_difference=loss,max_residual_difference=residual,passed=loss<=.05 and residual<=.01))
        certified=len(repeats)==2 and all(r['passed'] for r in screens)
        completed=[];incomplete=[]
        for path in sorted(self.output.glob('jacobian_*.json')):
            record=read(path)
            (completed if record['complete'] else incomplete).append((path.name,record))
        if completed:
            name,last=completed[-1];rank=last['diagnostics']['numerical_rank']
            identification=dict(status='local_numerical_rank_deficiency_identification_unverified' if rank<10 else 'local_numerical_full_rank_only',
                numerical_rank=rank,latest_complete_jacobian=name,incomplete_rounds=[x[0] for x in incomplete],
                statistical_identification_established=False,calibrated_smm_certified=False,
                caveat='Ten scored moments satisfy a count requirement; local numerical rank does not certify identification.')
        else:
            identification=dict(status='no_complete_local_jacobian',numerical_rank=None,
                incomplete_rounds=[x[0] for x in incomplete],statistical_identification_established=False,calibrated_smm_certified=False)
        result=dict(status='complete_with_repeat_screens' if certified else 'incomplete_or_unrepeated',lane=self.lane,
             config_sha256=self.pin,attempts=len(self.records),successful=sum(r['status']=='success' for r in self.records),
             elapsed_seconds=time.time()-self.start,stop_reason=self.stop_reason,selected=self.best,repeat_screens=screens,
             repeat_count=len(repeats),numerical_repeat_screens_passed=certified,promotion='disabled',
             remaining_evaluation_slots=self.b['max_evaluations']-len(self.records),identification_status=identification)
        write(self.output/'FINAL.json',result);self.save();return result


def main():
    require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm execution required')
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--lane',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();pin=sha(args.config)
    require(pin==os.environ.get('EXPECTED_TWO_STREAM_CONFIG_SHA256'),'Explicit config pin differs')
    config=read(args.config);require(config['schema']=='two_stream_search_v1','Wrong search schema')
    require(args.lane in config['lanes'],'Unknown stream')
    verify_launch(config,args.lane,pin)
    stream=Stream(config,args.lane,args.output,pin)
    try:result=stream.run()
    except Exception as exc:
        write(args.output/'CONTROLLER_FAILURE.json',dict(error=str(exc),error_type=type(exc).__name__,traceback=traceback.format_exc()))
        raise
    return 1 if stream.stop_reason=='fatal_integrity_or_scientific_failure' else 0


def verify_launch(config,lane,config_sha):
    pins=config['pins'];pins=list(pins.values()) if isinstance(pins,dict) else pins
    own=False
    for pin in pins:
        require(sha(pin['path'])==pin['sha256'],'Pinned source/input changed: '+pin['path'])
        own=own or Path(pin['path']).resolve()==Path(__file__).resolve()
    require(own,'Executing search source must be pinned')
    if config.get('synthetic',False):return
    approval_path=Path(os.environ['TWO_STREAM_LAUNCH_APPROVAL'])
    require(sha(approval_path)==os.environ.get('EXPECTED_TWO_STREAM_LAUNCH_APPROVAL_SHA256'),'Launch approval hash differs')
    approval=read(approval_path)
    require(approval.get('schema')=='two_stream_launch_approval_v1','Wrong launch approval schema')
    require(approval['status']=='approved_two_stream_overnight' and approval['config_sha256']==config_sha,'Wrong launch approval')
    require(lane in approval['allowed_lanes'],'Stream not authorized by approval')
    require(approval['not_before_epoch']<=time.time()<=approval['expires_epoch'],'Launch outside approved time window')
    require(approval['source_fingerprints']=={name:config['lanes'][name]['source_fingerprint'] for name in ('one_birth','two_birth')},
            'Approval source fingerprints differ')
    required=[('synthetic',approval['synthetic_receipt'],None),
              ('integration_one_birth',approval['integration_receipts']['one_birth'],'one_birth'),
              ('integration_two_birth',approval['integration_receipts']['two_birth'],'two_birth')]
    for label,pin,expected_lane in required:
        require(sha(pin['path'])==pin['sha256'],'Smoke receipt changed')
        receipt=read(pin['path'])
        require(receipt['status']==pin.get('expected_status','passed'),'Smoke did not pass')
        require(receipt['config_sha256']==config_sha,'Smoke used another configuration')
        if expected_lane is None:
            require(receipt.get('synthetic') is True,'Typed synthetic receipt required')
        else:
            require(receipt.get('lane')==expected_lane and receipt.get('real_model_evaluations')==2 and
                    receipt.get('source_fingerprint')==config['lanes'][expected_lane]['source_fingerprint'],
                    'Typed lane-specific integration receipt required: '+label)


if __name__=='__main__':sys.exit(main())
