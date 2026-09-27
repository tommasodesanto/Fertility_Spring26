#!/usr/bin/env python3
"""Finite single-arm controller with pinned science and owned subprocess deadlines.

The calibration runtime owns mathematics and all scientific gates. This module
owns dispatch, immutable normalization inputs, artifact verification and export.
No source rewriting, implicit retries or model imports occur in the controller.
"""
from __future__ import annotations
import argparse,csv,dataclasses,hashlib,importlib.util,json,math,os,random,shutil,signal,sys,time,traceback
from pathlib import Path

COMMON=('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation','theta0')
FREE=COMMON+('delta_alpha_jump','child_benefit_curvature')
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')
SCHEMA='e5f_single_share_overnight_v1'
FIT_NUMERIC=('target','model','gap','weight','loss_contribution')
PARAM_NUMERIC=('estimate','lower','upper')


def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()
def canon(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+f'.{os.getpid()}.tmp')
    q.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n');q.replace(p)
def module(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def science_id(c):
    # Approval metadata may change after smoke; numerical/source inputs may not.
    return canon({k:v for k,v in c.items() if k not in ('status','approval','verified_smoke','production_blockers')})

def norm_inputs(c):
    n=c['normalization']
    assert math.isfinite(n['initial_psi']) and math.isfinite(n['initial_step']) and n['initial_step']>0
    assert isinstance(n['warm_price'],bool) and n['maximum_stationary_solves']==23
    return dict(n)

def candidate_id(c,point):return canon(dict(science=science_id(c),point=point,normalization_inputs=norm_inputs(c)))

def verify(path):
    c=read(path);assert c['schema']==SCHEMA and not sys.flags.optimize
    assert os.environ.get('EXPECTED_UTILITY_OVERNIGHT_SHA256')==sha(path),'Unreviewed contract fingerprint'
    assert c['status'] in ('reviewed_smoke','approved_production')
    for item in list(c['files'].values())+[c['base_contract'],c['objective']]:assert sha(item['path'])==item['sha256'],item['path']
    assert c['files']['driver']['sha256']==sha(__file__)
    obj=read(c['objective']['path']);assert canon(obj)==c['objective_canonical_sha256']
    assert canon(obj['target_rows'])==c['target_weight_fingerprint']
    restrictions=obj['parameter_restrictions'];assert len(restrictions)==len(FREE)
    assert set(r['parameter'] for r in restrictions)==set(FREE)
    for r in restrictions:
        assert math.isfinite(r['lower']) and math.isfinite(r['upper']) and r['lower']<r['upper']
        assert r['lower']<=c['initial_point'][r['parameter']]<=r['upper']
    assert c['fixed']['delta_alpha']==0 and c['fixed']['sigma']==2
    budget=c['budget']
    assert budget['total_seconds']==28800 and 1<=budget['workers']<=10
    assert budget['search_seconds']==23400 and budget['repeat_seconds']==3600 and budget['export_seconds']==1800
    assert 1<=budget['rounds']<=30 and 1<=budget['points_per_round']<=10
    assert budget['rounds']*budget['points_per_round']<=300
    assert 0<budget['objective_cap_seconds']<=budget['repeat_seconds']
    assert 0<budget['smoke_seconds']<=budget['total_seconds']
    assert set(c['proposal_widths'])==set(FREE) and all(math.isfinite(v) and v>0 for v in c['proposal_widths'].values())
    assert c['pending_observer_mismatches'] and c['economic_changes'];norm_inputs(c)
    return c,obj

def setup(c,obj,point,out):
    runtime=module('overnight_calibration_runtime',c['files']['calibration_runtime']['path'])
    return runtime.setup(c,obj,point,out)

def evaluate(a,c,obj):
    request=read(a.request)
    assert request['contract_sha256']==sha(a.contract) and request['controller_pid']==os.getppid()
    assert request['point_sha256']==canon(request['point']) and request['normalization_inputs']==norm_inputs(c)
    assert request['scientific_candidate_id']==candidate_id(c,request['point'])
    a.output.mkdir(parents=True,exist_ok=False);context=request['context']
    assert context['candidate_id']==a.output.name and context['contract_sha256']==sha(a.contract)
    assert context['source_sha256']==c['files']['driver']['sha256'] and context['target_sha256']==c['objective']['sha256']
    assert context['point_sha256']==canon(request['point']) and context['stage']==request['stage']
    write(a.output/'startup.json',dict(context=context,pid=os.getpid(),parent_pid=os.getppid(),normalization_inputs=norm_inputs(c)))
    runtime=None;start=time.monotonic()
    try:
        old,tax,selected,rt,runner,native,evidence,binding=setup(c,obj,request['point'],a.output)
        runtime=(runner,native,evidence)
        write(a.output/'binding.json',dict(context=context,native=evidence,binding=binding))
        result=old.evaluate_point(tax=tax,objective=obj,selected=selected,runtime=rt,point=request['point'],output=a.output/'case',deadline_epoch=request['deadline_epoch'],graphs=request['graphs'])
        assert result['normalization']['psi_child']>0
        assert result['normalization_inputs']==norm_inputs(c),'Runtime normalization inputs differ'
        assert 1<=result['objective_stationary_solves']<=23
        assert result['target_system_sha256']==c['objective']['sha256']
        result.update(overnight_contract_sha256=sha(a.contract),scientific_identity=science_id(c),scientific_candidate_id=candidate_id(c,request['point']),target_weight_fingerprint=c['target_weight_fingerprint'],objective_canonical_sha256=c['objective_canonical_sha256'],wall_seconds=time.monotonic()-start)
        old.write(a.output/'case/receipt.json',result)
        write(a.output/'success.json',dict(context=context,loss=result['loss'],receipt_sha256=sha(a.output/'case/receipt.json'),checkpoint_sha256=result['case_checkpoint_sha256']))
    except Exception as exc:
        captured=None;status='fatal';classification='unknown_or_integrity_failure'
        if runtime:
            runner,native,evidence=runtime
            captured=runner.recovery.capture_native_failure(exc,expected_native_type=native,native_gate_tolerance=evidence['native_gate_tolerance'],context=runner.recovery.Context(**context))
            if runner.failure_status(exc)=='inadmissible_parameter_proposal':status='inadmissible';classification='existing_native_prefix'
            elif captured.narrow_infeasibility_verified:status='inadmissible';classification='authenticated_structured_native_gate'
        ledger=getattr(exc,'audit',getattr(exc,'ledger',None))
        if ledger is not None:write(a.output/'error_ledger.json',dict(context=context,ledger=ledger))
        write(a.output/'failure.json',dict(context=context,status=status,classification=classification,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),native=dataclasses.asdict(captured) if captured else None,error_ledger=ledger))
        raise

def keyed_csv(path,key):
    with Path(path).open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    mapped={r[key]:r for r in rows}
    assert len(mapped)==len(rows),'Duplicate table row';return mapped

def compare_tables(original,repeat):
    """Compare every numeric fit and parameter cell, preserving missing cells."""
    evidence={}
    for name,key,numeric in (('target_fit.csv','moment',FIT_NUMERIC),('parameters.csv','parameter',PARAM_NUMERIC)):
        left,right=keyed_csv(Path(original)/name,key),keyed_csv(Path(repeat)/name,key)
        assert left.keys()==right.keys(),name
        for row in left:
            for field in numeric:
                x,y=left[row][field],right[row][field]
                if x=='' or y=='':assert x==y,(name,row,field)
                else:assert math.isfinite(float(x)) and float(x)==float(y),(name,row,field,x,y)
            if name=='parameters.csv':assert left[row].get('near_bound')==right[row].get('near_bound'),row
        evidence[name]=dict(rows=len(left),numeric_fields=list(numeric),exact=True)
    return evidence

def validate_success(folder,c,request):
    folder=Path(folder);result=read(folder/'success.json')
    assert result['context']==request['context']
    receipt=read(folder/'case/receipt.json');assert result['receipt_sha256']==sha(folder/'case/receipt.json')
    assert receipt['overnight_contract_sha256']==request['contract_sha256']
    assert receipt['scientific_identity']==science_id(c) and receipt['scientific_candidate_id']==candidate_id(c,request['point'])
    assert receipt['target_system_sha256']==c['objective']['sha256'] and receipt['target_weight_fingerprint']==c['target_weight_fingerprint']
    assert receipt['objective_canonical_sha256']==c['objective_canonical_sha256']
    assert receipt['point']==request['point'] and receipt['normalization_inputs']==norm_inputs(c)
    assert receipt['normalization']['psi_child']>0 and 1<=receipt['objective_stationary_solves']<=23
    assert receipt['objective_stationary_solves']==receipt['normalization']['stationary_solves']
    assert sha(folder/'case/initial_state.pkl.gz')==result['checkpoint_sha256']==receipt['case_checkpoint_sha256']
    ledger=read(folder/'case/stationary_solves.json')
    assert len(ledger)==receipt['objective_stationary_solves'] and all(x['status']=='completed' for x in ledger)
    assert ledger[0]['psi_child']==norm_inputs(c)['initial_psi']
    if len(ledger)>1:assert math.isclose(abs(ledger[1]['psi_child']-ledger[0]['psi_child']),norm_inputs(c)['initial_step'],rel_tol=0,abs_tol=1e-14)
    fits=keyed_csv(folder/'case/target_fit.csv','moment');target_rows=read(c['objective']['path'])['target_rows']
    targets={r['restriction_id']:r for r in target_rows};assert fits.keys()==targets.keys()
    for name,row in fits.items():
        assert float(row['target'])==float(targets[name]['target'])
        assert math.isfinite(float(row['model'])) and math.isfinite(float(row['gap']))
        assert math.isclose(float(row['gap']),float(row['model'])-float(row['target']),rel_tol=1e-12,abs_tol=1e-12)
        w=targets[name]['actual_weight']
        if w is None:assert row['weight']==row['loss_contribution']==''
        else:
            assert float(row['weight'])==float(w)>0
            assert math.isclose(float(row['loss_contribution']),float(w)*float(row['gap'])**2,rel_tol=1e-12,abs_tol=1e-12)
    loss=sum(float(r['loss_contribution']) for r in fits.values() if r['loss_contribution']!='')
    assert math.isfinite(loss) and math.isclose(loss,float(result['loss']),rel_tol=1e-12) and float(result['loss'])==float(receipt['loss'])
    params=keyed_csv(folder/'case/parameters.csv','parameter');assert set(FREE)<=params.keys()
    restrictions={r['parameter']:r for r in read(c['objective']['path'])['parameter_restrictions']}
    for name in FREE:
        actual=float(params[name]['estimate']);requested=request['point'][name]
        if name=='beta_annual':assert math.isclose(actual,requested,rel_tol=2e-12,abs_tol=2e-12)
        else:assert actual==requested
        assert float(params[name]['lower'])==restrictions[name]['lower'] and float(params[name]['upper'])==restrictions[name]['upper']
    if request['graphs']:assert len(list((folder/'case/standard_diagnostics').glob('*.png')))==17
    return dict(loss=float(result['loss']),case_path=str((folder/'case').resolve()),checkpoint_sha256=result['checkpoint_sha256'],receipt_sha256=result['receipt_sha256'],scientific_candidate_id=candidate_id(c,request['point']))

def proposals(c,obj,center,count,round_id,seen):
    bounds={r['parameter']:(r['lower'],r['upper']) for r in obj['parameter_restrictions']}
    rng=random.Random(c['seed']+round_id);rows=[];attempts=0
    while len(rows)<count and attempts<count*100:
        attempts+=1;p=dict(center)
        names=list(FREE) if len(rows)%3==0 else [FREE[(len(rows)+round_id)%len(FREE)]]
        for k in names:
            lo,hi=bounds[k];width=c['proposal_widths'][k]*(1 if round_id<2 else .5)
            p[k]=min(hi,max(lo,p[k]+rng.gauss(0,width)))
        key=candidate_id(c,p)
        if key in seen:continue
        seen.add(key);rows.append(p)
    assert len(rows)==count,'Finite proposal construction exhausted';return rows

def classify(folder,c,req,process,code):
    """Integrity failures dominate censoring; timeout needs positive owner evidence."""
    data={};error=None
    try:
        startup=read(folder/'startup.json') if (folder/'startup.json').exists() else None
        if startup is not None:
            assert startup['context']==req['context'] and startup['parent_pid']==os.getpid()
            assert startup['pid']==process.process.pid and startup['normalization_inputs']==norm_inputs(c)
        failure=read(folder/'failure.json') if (folder/'failure.json').exists() else None
        if failure is not None:assert failure['context']==req['context']
        if failure and failure['status']=='fatal':return 'fatal',{},failure
        if process.observed_running_at_expiry and process.deadline_kill_reaped and code==-signal.SIGKILL:return 'censored_timeout',{},None
        assert startup is not None,'Missing authenticated startup'
        if process.deadline_expired:
            if code==0 and failure is None:return 'censored_late_completion',{},None
            return 'fatal',{},'Late nonzero exit without owned timeout'
        if code==0 and failure is None:return 'success',validate_success(folder,c,req['request']),None
        if code!=0 and failure and failure['status']=='inadmissible':
            assert failure['classification'] in ('existing_native_prefix','authenticated_structured_native_gate')
            if failure['classification']=='authenticated_structured_native_gate':
                assert failure['native']['narrow_infeasibility_verified'] and failure['native']['context']==req['context']
            return 'inadmissible',{},failure
        error=failure or 'Unclassified child exit'
    except Exception as exc:error=f'{type(exc).__name__}: {exc}'
    return 'fatal',data,error

def load_smoke(path,pin,c):
    assert sha(path)==pin;smoke=read(path)
    assert smoke['status']=='exact_loop_smoke_passed' and smoke['scientific_identity']==science_id(c)
    assert sha(smoke['contract_path'])==smoke['contract_sha256']
    assert len(smoke['records'])==2 and all(r['status']=='success' for r in smoke['records'])
    for r in smoke['records']:
        request=read(r['request_path']);assert request['point']==c['initial_point'] and request['graphs']
        verified=validate_success(Path(r['case_path']).parent,c,request)
        assert all(r[key]==value for key,value in verified.items())
    compare_tables(smoke['records'][0]['case_path'],smoke['records'][1]['case_path'])
    return smoke

def export_selected(selected,repeats,out,c,deadline):
    assert time.time()<deadline and len(repeats)==2
    comparisons=[]
    for r in repeats:
        assert r['status']=='success' and r['point']==selected['point'] and r['loss']==selected['loss']
        comparisons.append(compare_tables(selected['case_path'],r['case_path']))
    out.mkdir(exist_ok=False)
    source=Path(selected['case_path'])
    for path in source.iterdir():
        if time.time()>=deadline:raise TimeoutError('Export reserve exhausted')
        if path.is_file() and path.suffix in ('.csv','.json'):shutil.copy2(path,out/path.name)
    graphs=Path(repeats[0]['case_path'])/'standard_diagnostics'
    assert len(list(graphs.glob('*.png')))==17
    shutil.copytree(graphs,out/'standard_diagnostics')
    (out/'initial_state.pkl.gz').symlink_to((source/'initial_state.pkl.gz').resolve())
    if time.time()>=deadline:raise TimeoutError('Export finished beyond the original total deadline')
    write(out/'export_receipt.json',dict(status='verified_selected_export',selected=selected,repeats=repeats,comparisons=comparisons,graph_source=str(graphs),scientific_identity=science_id(c),target_weight_fingerprint=c['target_weight_fingerprint']))

def controller(a,c,obj):
    sys.path.insert(0,c['runtime_tools']);sys.path.insert(0,str(Path(c['files']['recovery_policy']['path']).parent))
    supervision=module('overnight_reviewed_supervision',c['files']['recovery_search']['path'])
    smoke=load_smoke(a.smoke_receipt,a.smoke_sha256,c) if a.stage=='search' else None
    a.output.mkdir(parents=True,exist_ok=False);start=time.time();budget=c['budget']
    cutoff=start+(budget['smoke_seconds'] if a.stage=='smoke' else budget['search_seconds'])
    records=[];best=dict(smoke['records'][0]) if smoke else None;seen={candidate_id(c,c['initial_point'])};fatal=False;batches=[]
    clock=dict(start=start,search_cutoff=cutoff,repeat_cutoff=start+budget['search_seconds']+budget['repeat_seconds'],end=start+budget['total_seconds'],contract_sha256=sha(a.contract))
    write(a.output/'clock.json',clock)
    if best:write(a.output/'best_so_far.json',best)
    last_heartbeat=0.
    def heartbeat(**progress):
        nonlocal last_heartbeat
        if time.time()-last_heartbeat>=30 or progress.get('force'):
            write(a.output/'heartbeat.json',dict(epoch=time.time(),stage=a.stage,completed=len(records),best_loss=best['loss'] if best else None,**progress));last_heartbeat=time.time()
    def batch(points,stage,deadline,graphs=False):
        nonlocal best,fatal
        requests=[]
        for p in points:
            name=f'{stage}_{len(records)+len(requests):04d}'
            context=dict(candidate_id=name,stage=stage,contract_sha256=sha(a.contract),source_sha256=c['files']['driver']['sha256'],target_sha256=c['objective']['sha256'],point_sha256=canon(p))
            requests.append(dict(id=name,point=p,context=context))
        def launch(req,batch_deadline):
            verify(a.contract);end=min(batch_deadline,time.time()+budget['objective_cap_seconds']);folder=a.output/req['id']
            request=dict(req,stage=stage,point_sha256=canon(req['point']),scientific_candidate_id=candidate_id(c,req['point']),normalization_inputs=norm_inputs(c),contract_sha256=sha(a.contract),controller_pid=os.getpid(),deadline_epoch=end,graphs=graphs)
            path=a.output/(req['id']+'.request.json');write(path,request);req['request']=request;req['request_path']=str(path.resolve())
            env=os.environ.copy();env.update({k:'1' for k in THREADS})
            cmd=[sys.executable,__file__,'--stage','evaluate','--contract',str(a.contract),'--output',str(folder),'--request',str(path)]
            return supervision.ManagedProcess(cmd,a.output/(req['id']+'.log'),end,env)
        def finish(req,process,code):
            nonlocal best,fatal
            try:verify(a.contract);status,data,error=classify(a.output/req['id'],c,req,process,code)
            except Exception as exc:status,data,error='fatal',{},str(exc)
            record=dict(case=req['id'],status=status,point=req['point'],request_path=req['request_path'],execution=dict(returncode=code,deadline=process.deadline,observed=process.observed_epoch,owned_timeout=process.observed_running_at_expiry and process.deadline_kill_reaped),error=error)
            record.update(data);records.append(record);write(a.output/'latest_completed.json',record)
            if status=='success' and stage!='repeat' and (best is None or record['loss']<best['loss']):best=record;write(a.output/'best_so_far.json',best)
            if status=='fatal' or stage in ('smoke','repeat') and status!='success':fatal=True
            record['halt_new_dispatch']=fatal
            write(a.output/'checkpoint.json',dict(records=records,best=best,seen=sorted(seen),clock=clock,scientific_identity=science_id(c),contract_sha256=sha(a.contract)))
            heartbeat(force=True);return record
        result=supervision.run_batch(requests,workers=min(budget['workers'],len(requests)),deadline=deadline,launch=launch,finish=finish,heartbeat=heartbeat,allowed_statuses={'success'} if stage in ('smoke','repeat') else {'success','inadmissible','censored_timeout','censored_late_completion'},guard=lambda:'fatal_stop' if fatal else None)
        batches.append(dict(stage=stage,requested=len(points),**result));write(a.output/'batches.json',batches)
        if any(r['status']=='failed' for r in result['results']):fatal=True
        return result
    try:
        if a.stage=='smoke':
            result=batch([c['initial_point'],c['initial_point']],'smoke',cutoff,graphs=True)
            assert not fatal and result['complete'] and len(records)==2
            comparison=compare_tables(records[0]['case_path'],records[1]['case_path'])
            assert records[0]['loss']==records[1]['loss']
            write(a.output/'complete.json',dict(status='exact_loop_smoke_passed',records=records,comparison=comparison,scientific_identity=science_id(c),contract_path=str(a.contract.resolve()),contract_sha256=sha(a.contract)));return
        for round_id in range(budget['rounds']):
            if fatal or time.time()>=cutoff:break
            points=proposals(c,obj,best['point'],budget['points_per_round'],round_id,seen)
            result=batch(points,'initial' if round_id==0 else 'de',cutoff)
            if not result['complete']:break
        selected=dict(best);write(a.output/'selected.json',dict(selected=selected,scientific_identity=science_id(c),frozen_before_repeats=True))
        if fatal:raise RuntimeError('Fatal candidate failure; no further dispatch')
        result=batch([selected['point'],selected['point']],'repeat',clock['repeat_cutoff'],graphs=True)
        assert not fatal and result['complete'] and len(result['results'])==2,'Selected repetitions incomplete'
        export_selected(selected,result['results'],a.output/'selected_export',c,clock['end'])
        write(a.output/'complete.json',dict(status='bounded_search_complete',records=records,selected=selected,elapsed_seconds=time.time()-start,clock=clock,scientific_identity=science_id(c),contract_sha256=sha(a.contract)))
    except Exception as exc:
        write(a.output/'complete.json',dict(status='incomplete_or_fatal_stop',error_type=type(exc).__name__,error=str(exc),records=records,best=best,clock=clock,scientific_identity=science_id(c)));raise

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True,choices=['preflight','smoke','search','evaluate']);ap.add_argument('--contract',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--request',type=Path);ap.add_argument('--smoke-receipt',type=Path);ap.add_argument('--smoke-sha256');a=ap.parse_args()
    assert os.environ.get('SLURM_JOB_ID'),'Torch Slurm execution only'
    c,obj=verify(a.contract)
    if a.stage=='search':
        assert (c['status']=='approved_production' and a.smoke_receipt and a.smoke_sha256
                and c['approval']['production_authorized'] is True
                and not c.get('production_blockers')), 'Author approval and verified preparation are required'
    if a.stage=='preflight':
        a.output.mkdir(parents=True,exist_ok=False)
        *_,binding=setup(c,obj,c['initial_point'],a.output)
        write(a.output/'preflight.json',dict(status='zero_solve_passed',contract_sha256=sha(a.contract),binding=binding));return
    if a.stage=='evaluate':evaluate(a,c,obj)
    else:controller(a,c,obj)
if __name__=='__main__':main()
