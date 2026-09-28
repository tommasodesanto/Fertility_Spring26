#!/usr/bin/env python3
"""Versioned overnight controller; unchanged evening household runtime/gates."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,importlib.util,json,math,os,pickle,random,shutil,signal,sys,time,traceback
from pathlib import Path
FREE=('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation','theta0','delta_alpha_jump','child_benefit_curvature','tenure_choice_kappa')
LANES=('primary','identity','block')
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')
SCHEMA='e5f_night_v1'

def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()
def canon(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+f'.{os.getpid()}.tmp');q.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n');q.replace(p)
def module(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def table(p,key):
    with Path(p).open(newline='') as f:rows=list(csv.DictReader(f))
    result={r[key]:r for r in rows};assert len(result)==len(rows),'Duplicate table row';return result
def science(c):return canon(c)
def identity(c,lane,point):return canon(dict(science=science(c),lane=lane,point=point,normalization=c['normalization']))
def verify_anchor_pins(c):
    anchor=c['anchor'];assert len(anchor['repeats'])==2
    for case in [anchor]+anchor['repeats']:
        for name in ('receipt','target_fit','parameters'):
            item=case[name];assert sha(item['path'])==item['sha256'],'Anchor artifact changed: '+item['path']
def clock(c,now):
    b=c['budget'];start=b['absolute_start_epoch'];end=b['absolute_end_epoch']
    assert end==1790596800 and start<=now<end,'Outside authorized overnight window'
    return dict(start=start,end=end,search_cutoff=end-3600,repeat_cutoff=end-600)
def verify(path):
    c=read(path);assert c['schema']==SCHEMA and not sys.flags.optimize
    assert os.environ.get('EXPECTED_E5F_NIGHT_SHA256')==sha(path),'Contract changed'
    pins=list(c['files'].values())+[c['source_manifest']]
    for item in pins:assert sha(item['path'])==item['sha256'],item['path']
    assert c['files']['driver']['sha256']==sha(__file__)
    assert set(c['lanes'])==set(LANES) and set(c['initial_point'])==set(FREE)
    objectives={};common=None
    assert len(c['validation_rows'])==3 and len(set(c['validation_rows']))==3
    for lane in LANES:
        spec=c['lanes'][lane];assert sha(spec['objective']['path'])==spec['objective']['sha256']
        obj=read(spec['objective']['path']);objectives[lane]=obj
        assert canon(obj)==spec['canonical_sha256'] and canon(obj['target_rows'])==spec['target_weight_fingerprint']
        rows={r['restriction_id']:r for r in obj['target_rows']};assert len(rows)==len(obj['target_rows'])==14
        assert all(rows[k]['actual_weight']==0 for k in c['validation_rows'])
        assert sum(r['actual_weight'] is None for r in rows.values())==1
        assert sum(r['actual_weight'] is not None and r['actual_weight']>0 for r in rows.values())==10
        restrictions={r['parameter']:r for r in obj['parameter_restrictions']};assert set(restrictions)==set(FREE)
        for k,r in restrictions.items():assert math.isfinite(r['lower']) and math.isfinite(r['upper']) and r['lower']<r['upper'] and r['lower']<=c['initial_point'][k]<=r['upper']
        invariant=([(r['restriction_id'],r['target']) for r in obj['target_rows']],obj['parameter_restrictions'])
        if common is None:common=invariant
        else:assert invariant==common,'Lanes differ beyond weights'
    b=c['budget'];assert b['total_seconds']==b['absolute_end_epoch']-b['absolute_start_epoch'] and 1<=b['workers']<=24
    assert b['absolute_end_epoch']==1790596800
    assert b['search_reserve_seconds']==3600 and b['export_reserve_seconds']==600
    assert b['max_objective_cases']==746 and b['maximum_diagnostic_objectives']<=12
    assert 3<=b['max_search_cases']<=720 and b['max_search_cases']%3==0
    assert b['objective_cap_seconds']==1800
    assert c['log_parameters']==['tenure_choice_kappa']
    tenure={r['parameter']:r for r in objectives['primary']['parameter_restrictions']}['tenure_choice_kappa'];assert tenure['lower']>0
    assert set(c['proposal_widths'])==set(FREE) and all(math.isfinite(x) and x>0 for x in c['proposal_widths'].values())
    n=c['normalization'];assert n['initial_psi']>0 and n['initial_step']>0 and n['maximum_stationary_solves']==23 and isinstance(n['warm_price'],bool)
    assert c['fixed']['due'] is True and c['fixed']['delta_alpha']==0 and c['fixed']['sigma']==2
    assert len(c['standard_diagnostic_names'])==len(set(c['standard_diagnostic_names']))==17
    parent=read(c['evening_contract']['path']);assert sha(c['evening_contract']['path'])==c['evening_contract']['sha256']
    assert c['normalization']==parent['normalization'] and c['fixed']==parent['fixed']
    assert all(read(c['lanes'][lane]['objective']['path'])==read(parent['lanes'][lane]['objective']['path']) for lane in LANES)
    verify_anchor_pins(c)
    assert read(c['anchor']['receipt']['path'])['point']==c['initial_point']
    history=c['timeout_replay_history'];assert sha(history['path'])==history['sha256']
    old_records={r['case']:r for r in read(history['path'])['records']}
    for design,name in zip(c['initial_design'],('initial_0044_block','initial_0186_primary')):
        assert old_records[name]['status']=='censored_timeout' and design['point']==old_records[name]['point']
    assert len(c['initial_design'])==2
    return c,objectives

def compare_tables(left,right):
    for name,key,fields in [('target_fit.csv','moment',('target','model','gap','weight','loss_contribution')),('parameters.csv','parameter',('estimate','lower','upper'))]:
        a,b=table(Path(left)/name,key),table(Path(right)/name,key);assert a.keys()==b.keys()
        for k in a:
            for field in fields:
                x,y=a[k][field],b[k][field]
                if x=='' or y=='':assert x==y,(name,k,field)
                else:assert math.isfinite(float(x)) and float(x)==float(y),(name,k,field)
            if name=='parameters.csv':assert a[k].get('near_bound')==b[k].get('near_bound')
    return {'target_rows':14,'parameter_rows':31,'exact_numeric_cells':True}

def compare_anchor(anchor,case):
    """Weights may differ by lane; physical moments and all31 parameters may not."""
    for filename,key,fields,size in [('target_fit.csv','moment',('target','model','gap'),14),('parameters.csv','parameter',('estimate','lower','upper'),31)]:
        left,right=table(Path(anchor)/filename,key),table(Path(case)/filename,key)
        assert len(left)==len(right)==size and left.keys()==right.keys()
        for row in left:
            for field in fields:
                x,y=left[row][field],right[row][field]
                assert x==y if x=='' or y=='' else math.isfinite(float(x)) and float(x)==float(y),(filename,row,field)
    return dict(target_rows=14,parameter_rows=31,physical_numeric_cells_exact=True)

def render_saved(a,c):
    req=read(a.request);ctx=req['context'];source=Path(req['saved_case']['case_path'])
    assert req['contract_sha256']==sha(a.contract) and req['controller_pid']==os.getppid()
    assert sha(source/'receipt.json')==req['saved_case']['receipt_sha256']
    assert sha(source/'initial_state.pkl.gz')==req['saved_case']['checkpoint_sha256']
    assert read(source/'receipt.json')['point']==req['point']
    a.output.mkdir(parents=True,exist_ok=False)
    write(a.output/'startup.json',dict(context=ctx,pid=os.getpid(),parent_pid=os.getppid(),normalization_inputs=c['normalization']))
    runtime=module('night_saved_runtime',c['files']['runtime']['path'])
    lane=req['lane'];obj=read(c['lanes'][lane]['objective']['path'])
    evaluator=runtime.setup(dict(c,objective=c['lanes'][lane]['objective']),obj,a.output)
    with gzip.open(source/'initial_state.pkl.gz','rb') as stream:packet=pickle.load(stream)
    case=a.output/'case';case.mkdir()
    for p in source.iterdir():
        if p.is_file() and p.suffix in ('.json','.csv'):shutil.copy2(p,case/p.name)
    evaluator.rt['audit'].standard_diagnostics(packet,case,validate_production_young=False)
    assert sorted(p.name for p in (case/'standard_diagnostics').glob('*.png'))==sorted(c['standard_diagnostic_names'])
    compare_tables(source,case)
    assert time.time()<req['deadline_epoch']
    write(a.output/'render_receipt.json',dict(status='saved_checkpoint_rendered',context=ctx,source=req['saved_case'],objective_evaluations=0,case_path=str(case.resolve())))

def evaluate(a,c,objs):
    req=read(a.request);lane=req['lane'];ctx=req['context']
    assert req['contract_sha256']==sha(a.contract) and req['controller_pid']==os.getppid()
    assert req['scientific_candidate_id']==identity(c,lane,req['point']) and req['normalization_inputs']==c['normalization']
    assert ctx['candidate_id']==a.output.name and ctx['point_sha256']==canon(req['point'])
    assert ctx['contract_sha256']==sha(a.contract) and ctx['source_sha256']==c['files']['driver']['sha256'] and ctx['target_sha256']==c['lanes'][lane]['objective']['sha256']
    a.output.mkdir(parents=True,exist_ok=False);write(a.output/'startup.json',dict(context=ctx,pid=os.getpid(),parent_pid=os.getppid(),normalization_inputs=c['normalization']))
    evaluator=None;phase='setup'
    try:
        runtime=module('evening_runtime',c['files']['runtime']['path'])
        runtime_contract=dict(c,objective=c['lanes'][lane]['objective'])
        evaluator=runtime.setup(runtime_contract,objs[lane],a.output)
        phase='objective';receipt=evaluator.evaluate(req['point'],a.output/'case',req['deadline_epoch'],graphs=req['graphs'],due=True,fixed_reference=False)
        receipt.update(contract_sha256=sha(a.contract),scientific_identity=science(c),scientific_candidate_id=identity(c,lane,req['point']),lane=lane,target_weight_fingerprint=c['lanes'][lane]['target_weight_fingerprint'])
        write(a.output/'case/receipt.json',receipt)
        write(a.output/'success.json',dict(context=ctx,receipt_sha256=sha(a.output/'case/receipt.json'),loss=receipt['loss'],checkpoint_sha256=receipt['case_checkpoint_sha256']))
    except Exception as exc:
        failure=dict(status='fatal',classification='unknown_or_integrity_failure')
        if evaluator is not None:
            try:failure.update(evaluator.classify_failure(exc,ctx))
            except Exception as classification_error:failure['classifier_error']=str(classification_error);failure['status']='fatal'
        assert failure['status'] in ('fatal','inadmissible')
        failure.update(context=ctx,phase=phase,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
        audit=getattr(exc,'audit',getattr(exc,'ledger',None))
        if audit is not None:write(a.output/'error_ledger.json',dict(context=ctx,audit=audit))
        write(a.output/'failure.json',failure);raise

def validate(folder,c,objs,req):
    folder=Path(folder);case=folder/'case';lane=req['lane'];spec=c['lanes'][lane]
    side=read(folder/'success.json');r=read(case/'receipt.json');assert side['context']==req['context'] and side['receipt_sha256']==sha(case/'receipt.json')
    assert r['contract_sha256']==req['contract_sha256'] and r['scientific_identity']==science(c) and r['scientific_candidate_id']==identity(c,lane,req['point'])
    assert r['lane']==lane and r['point']==req['point'] and r['normalization_inputs']==c['normalization']
    assert r['target_system_sha256']==spec['objective']['sha256'] and r['target_weight_fingerprint']==spec['target_weight_fingerprint']
    assert r['source_manifest_sha256']==c['source_manifest']['sha256']
    assert sha(case/'initial_state.pkl.gz')==side['checkpoint_sha256']==r['case_checkpoint_sha256']
    assert r['normalization']['psi_child']>0 and 1<=r['objective_stationary_solves']<=23
    ledger=read(case/'stationary_solves.json');assert len(ledger)==r['objective_stationary_solves']==r['normalization']['stationary_solves']
    assert all(x['status']=='completed' for x in ledger) and ledger[0]['psi_child']==c['normalization']['initial_psi']
    if len(ledger)>1:assert math.isclose(abs(ledger[1]['psi_child']-ledger[0]['psi_child']),c['normalization']['initial_step'],rel_tol=0,abs_tol=1e-14)
    fits=table(case/'target_fit.csv','moment');targets={x['restriction_id']:x for x in objs[lane]['target_rows']};assert len(fits)==14 and fits.keys()==targets.keys()
    total=0.;primary=0.;primary_rows={x['restriction_id']:x for x in objs['primary']['target_rows']}
    for k,row in fits.items():
        target=float(row['target']);model=float(row['model']);gap=float(row['gap']);assert target==targets[k]['target'] and math.isfinite(model) and math.isclose(gap,model-target,rel_tol=1e-12,abs_tol=1e-12)
        w=targets[k]['actual_weight']
        if w is None:assert row['weight']==row['loss_contribution']==''
        else:
            loss=float(row['loss_contribution']);assert float(row['weight'])==w>=0 and math.isclose(loss,w*gap*gap,rel_tol=1e-12,abs_tol=1e-12);total+=loss
            if k in c['validation_rows']:assert row['role']=='validation' and loss==0
        wp=primary_rows[k]['actual_weight'];primary+=(wp or 0)*gap*gap
    assert math.isfinite(total) and math.isclose(total,r['loss'],rel_tol=1e-12,abs_tol=1e-12) and r['loss']==side['loss']
    params=table(case/'parameters.csv','parameter');assert len(params)==31 and set(FREE)<=params.keys()
    assert all(math.isfinite(float(row['estimate'])) for row in params.values())
    bounds={x['parameter']:x for x in objs[lane]['parameter_restrictions']}
    for k in FREE:
        actual=float(params[k]['estimate']);wanted=req['point'][k]
        assert math.isclose(actual,wanted,rel_tol=2e-12,abs_tol=2e-12) if k=='beta_annual' else actual==wanted
        assert float(params[k]['lower'])==bounds[k]['lower'] and float(params[k]['upper'])==bounds[k]['upper']
    if req['graphs']:assert sorted(p.name for p in (case/'standard_diagnostics').glob('*.png'))==sorted(c['standard_diagnostic_names'])
    return dict(loss=total,primary_rescore=primary,case_path=str(case.resolve()),receipt_sha256=side['receipt_sha256'],checkpoint_sha256=side['checkpoint_sha256'])

def classify(folder,c,objs,req,proc,code):
    try:
        startup=read(folder/'startup.json') if (folder/'startup.json').exists() else None
        if startup:assert startup['context']==req['context'] and startup['pid']==proc.process.pid and startup['parent_pid']==os.getpid() and startup['normalization_inputs']==c['normalization']
        failure=read(folder/'failure.json') if (folder/'failure.json').exists() else None
        if failure:assert failure['context']==req['context']
        if failure and failure['status']=='fatal':return 'fatal',{},failure
        if proc.observed_running_at_expiry and proc.deadline_kill_reaped and code==-signal.SIGKILL:return 'censored_timeout',{},None
        assert startup,'Missing startup'
        if proc.deadline_expired:return ('censored_late_completion',{},None) if code==0 and not failure else ('fatal',{},'Late nonzero exit without owned timeout')
        if code==0 and not failure:return 'success',validate(folder,c,objs,req),None
        if code!=0 and failure and failure['status']=='inadmissible':
            assert failure['classification'] in ('inherited_native_prefix','authenticated_structured_native_gate','explicit_economic_gate')
            assert failure['authenticated'] is True,'Missing runtime authentication'
            return 'inadmissible',{},failure
        return 'fatal',{},failure or 'Unknown child exit'
    except Exception as exc:return 'fatal',{},f'{type(exc).__name__}: {exc}'

def propose(c,obj,lane,center,index,seen):
    design=c.get('initial_design',[])
    if index<len(design):
        p=dict(design[index]['point']);key=identity(c,lane,p)
        assert key not in seen,'Declared design duplicates an already evaluated point'
        seen.add(key);return p
    index-=len(design)
    rng=random.Random(c['seed']+100000*LANES.index(lane)+index);bounds={x['parameter']:x for x in obj['parameter_restrictions']}
    for attempt in range(100):
        p=dict(center);mode=index%8
        names=FREE if mode in (0,1,2,3) else ('chi','delta_alpha_jump','tenure_choice_kappa') if mode==4 else ('child_benefit_curvature','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation') if mode==5 else (FREE[(index//8+mode)%len(FREE)],)
        for k in names:
            lo,hi=bounds[k]['lower'],bounds[k]['upper'];old=p[k]
            if k in c['log_parameters']:lo,hi,old=math.log(lo),math.log(hi),math.log(old)
            # No full-box draws: evening broad trials had no successful solves.
            width=c['proposal_widths'][k]*(2 if mode==0 else 1)
            value=min(hi,max(lo,old+rng.gauss(0,width)))
            p[k]=math.exp(value) if k in c['log_parameters'] else value
            p[k]=min(bounds[k]['upper'],max(bounds[k]['lower'],p[k]))
        key=identity(c,lane,p)
        if key not in seen:seen.add(key);return p
    raise RuntimeError('Finite distinct-proposal budget exhausted')

def proposal_label(c,index):
    design=c.get('initial_design',[])
    if index<len(design):return design[index]['label']
    index-=len(design)
    mode=index%8
    return 'moderate_joint' if mode==0 else 'joint_local' if mode in (1,2,3) else 'housing_subspace' if mode==4 else 'fertility_subspace' if mode==5 else 'coordinate_'+FREE[(index//8+mode)%len(FREE)]

def final_selection(best,common_best):
    selected={lane:dict(best[lane]) for lane in LANES}
    key=next((name for name,row in selected.items() if row['point']==common_best['point']),None)
    if key is None:selected['common_primary']=dict(common_best);key='common_primary'
    return selected,key

def verified_smoke(c,objs,path,pin,contract_sha256):
    assert sha(path)==pin;s=read(path);assert s['status']=='exact_loop_smoke_passed' and s['contract_sha256']==contract_sha256 and s['scientific_identity']==science(c)
    assert len(s['records'])==6
    for lane in LANES:
        rows=[r for r in s['records'] if r['lane']==lane];assert len(rows)==2
        for r in rows:
            req=read(r['request_path']);assert req['point']==c['initial_point'] and req['graphs'] and r['status']=='success'
            actual=validate(Path(r['case_path']).parent,c,objs,req);assert all(r[k]==v for k,v in actual.items())
            compare_anchor(c['anchor']['case_path'],r['case_path'])
        compare_tables(rows[0]['case_path'],rows[1]['case_path'])
    return s

def controller(a,c,objs):
    _CONTRACT_SHA=sha(a.contract);timing=clock(c,time.time());budget=c['budget']
    sys.path.insert(0,str(Path(c['files']['recovery_search']['path']).parent))
    supervisor=module('evening_supervisor',c['files']['recovery_search']['path'])
    smoke=None
    if a.stage=='search':
        assert sha(a.approval)==a.approval_sha256;approval=read(a.approval)
        assert approval['status']=='approved_search' and approval['contract_sha256']==_CONTRACT_SHA
        smoke=verified_smoke(c,objs,approval['smoke_receipt']['path'],approval['smoke_receipt']['sha256'],_CONTRACT_SHA)
    a.output.mkdir(parents=True,exist_ok=False);write(a.output/'clock.json',timing)
    records=[];best={lane:next((dict(r) for r in smoke['records'] if r['lane']==lane),None) if smoke else None for lane in LANES}
    common_best=min(smoke['records'],key=lambda row:row['primary_rescore']) if smoke else None
    seen={identity(c,lane,c['initial_point']) for lane in LANES};fatal=False;last=0.
    reports=[];next_report=time.time()+3600
    if smoke:write(a.output/'best_so_far.json',best)
    def heartbeat(**extra):
        nonlocal last
        if time.time()-last>=5 or extra.get('force'):
            write(a.output/'heartbeat.json',dict(epoch=time.time(),stage=a.stage,completed=len(records),best={k:v['loss'] if v else None for k,v in best.items()},clock=timing,**extra));last=time.time()
    def batch(points,stage,deadline,graphs):
        nonlocal fatal
        reqs=[]
        for item in points:
            lane,p=item[:2];design=item[2] if len(item)>2 else stage
            name=f'{stage}_{len(records)+len(reqs):04d}_{lane}'
            ctx=dict(candidate_id=name,stage=stage,contract_sha256=_CONTRACT_SHA,source_sha256=c['files']['driver']['sha256'],target_sha256=c['lanes'][lane]['objective']['sha256'],point_sha256=canon(p))
            reqs.append(dict(id=name,lane=lane,point=p,context=ctx,design=design,**({'saved_case':item[3]} if stage=='render' else {})))
        def launch(req,batch_deadline):
            end=min(batch_deadline,time.time()+budget['objective_cap_seconds']);folder=a.output/req['id']
            payload=dict(req,stage=stage,contract_sha256=_CONTRACT_SHA,scientific_candidate_id=identity(c,req['lane'],req['point']),normalization_inputs=c['normalization'],controller_pid=os.getpid(),deadline_epoch=end,graphs=graphs)
            path=a.output/(req['id']+'.request.json');write(path,payload);req['payload']=payload;req['request_path']=str(path.resolve())
            env=os.environ.copy();env.update({k:'1' for k in THREADS})
            return supervisor.ManagedProcess([sys.executable,__file__,'--stage','render' if stage=='render' else 'evaluate','--contract',str(a.contract),'--output',str(folder),'--request',str(path)],a.output/(req['id']+'.log'),end,env)
        def finish(req,proc,code):
            nonlocal fatal,common_best
            try:
                verify(a.contract)
                if stage=='render':
                    folder=a.output/req['id'];r=read(folder/'render_receipt.json');startup=read(folder/'startup.json')
                    assert code==0 and not proc.deadline_expired and r['context']==req['context']
                    assert startup['context']==req['context'] and startup['pid']==proc.process.pid and startup['parent_pid']==os.getpid()
                    assert r['source']==req['saved_case'] and r['objective_evaluations']==0
                    compare_tables(req['saved_case']['case_path'],r['case_path'])
                    assert sorted(p.name for p in (Path(r['case_path'])/'standard_diagnostics').glob('*.png'))==sorted(c['standard_diagnostic_names'])
                    status,data,error='success',dict(case_path=r['case_path']),None
                else:status,data,error=classify(a.output/req['id'],c,objs,req['payload'],proc,code)
            except Exception as exc:status,data,error='fatal',{},str(exc)
            row=dict(case=req['id'],lane=req['lane'],point=req['point'],design=req['design'],request_path=req['request_path'],status=status,error=error,returncode=code,deadline=proc.deadline,**data)
            if stage=='render':
                reports.append(row);write(a.output/'hourly_reports.json',reports)
                if status!='success':fatal=True
                row['halt_new_dispatch']=fatal;heartbeat(force=True);return row
            records.append(row);write(a.output/'latest_completed.json',row);write(a.output/f"latest_{req['lane']}.json",row)
            if status=='success' and stage!='repeat' and (best[req['lane']] is None or row['loss']<best[req['lane']]['loss']):best[req['lane']]=row;write(a.output/'best_so_far.json',best);write(a.output/f"best_{req['lane']}.json",row)
            if status=='success' and stage!='repeat' and (common_best is None or row['primary_rescore']<common_best['primary_rescore']):common_best=row;write(a.output/'best_common_primary.json',row)
            if status=='fatal' or stage in ('smoke','repeat') and status!='success':fatal=True
            row['halt_new_dispatch']=fatal
            write(a.output/'checkpoint.json',dict(records=records,best=best,common_primary_best=common_best,reports=reports,seen=sorted(seen),clock=timing,contract_sha256=_CONTRACT_SHA));heartbeat(force=True);return row
        result=supervisor.run_batch(reqs,workers=min(budget['workers'],len(reqs)),deadline=deadline,launch=launch,finish=finish,heartbeat=heartbeat,poll_seconds=1,allowed_statuses={'success'} if stage in ('smoke','repeat') else {'success','inadmissible','censored_timeout','censored_late_completion'},guard=lambda:'fatal_stop' if fatal else None)
        if any(r['status']=='failed' for r in result['results']):fatal=True
        write(a.output/f'batch_{stage}_{len(records):04d}.json',result);return result
    try:
        if a.stage=='smoke':
            result=batch([(lane,c['initial_point']) for lane in LANES for _ in range(2)],'smoke',timing['search_cutoff'],True)
            assert not fatal and result['complete'] and len(records)==6
            for lane in LANES:
                pair=[r for r in records if r['lane']==lane];assert pair[0]['loss']==pair[1]['loss'];compare_tables(pair[0]['case_path'],pair[1]['case_path'])
                for row in pair:compare_anchor(c['anchor']['case_path'],row['case_path'])
            write(a.output/'complete.json',dict(status='exact_loop_smoke_passed',records=records,contract_sha256=_CONTRACT_SHA,scientific_identity=science(c)));return
        counts=dict.fromkeys(LANES,0);limit=budget['max_search_cases']//3
        while not fatal and time.time()<timing['search_cutoff'] and any(x<limit for x in counts.values()):
            points=[]
            for _ in range(max(1,budget['workers']//3)):
                for lane in LANES:
                    if counts[lane]>=limit:continue
                    points.append((lane,propose(c,objs[lane],lane,best[lane]['point'],counts[lane],seen),proposal_label(c,counts[lane])));counts[lane]+=1
            result=batch(points,'initial',timing['search_cutoff'],False)
            if not result['complete']:break
            if not fatal and time.time()>=next_report and time.time()<timing['search_cutoff']:
                snapshot={lane:dict(best[lane]) for lane in LANES}
                write(a.output/f'hourly_snapshot_{len(reports):04d}.json',dict(best=snapshot,common_primary_best=common_best,epoch=time.time()))
                batch([(row['lane'],row['point'],'hourly_'+lane,row) for lane,row in snapshot.items()],'render',timing['search_cutoff'],True)
                next_report=time.time()+3600
        selected,common_key=final_selection(best,common_best)
        write(a.output/'selected.json',dict(selected=selected,common_primary_key=common_key,frozen_before_repeats=True,contract_sha256=_CONTRACT_SHA))
        if fatal:raise RuntimeError('Fatal case; dispatch stopped and siblings drained')
        result=batch([(row['lane'],row['point'],key) for key,row in selected.items() for _ in range(2)],'repeat',timing['repeat_cutoff'],True)
        assert not fatal and result['complete'] and len(result['results'])==2*len(selected)
        for lane,original in selected.items():
            repeats=[r for r in result['results'] if r['design']==lane]
            for repeat in repeats:
                assert repeat['point']==original['point'] and repeat['loss']==original['loss'];compare_tables(original['case_path'],repeat['case_path'])
            out=a.output/'selected_export'/lane;out.mkdir(parents=True,exist_ok=False);source=Path(original['case_path'])
            for p in source.iterdir():
                if time.time()>=timing['end']:raise TimeoutError('Original global deadline reached during export')
                if p.is_file() and p.suffix in ('.json','.csv'):shutil.copy2(p,out/p.name)
            shutil.copytree(Path(repeats[0]['case_path'])/'standard_diagnostics',out/'standard_diagnostics');(out/'initial_state.pkl.gz').symlink_to((source/'initial_state.pkl.gz').resolve())
            write(out/'export_receipt.json',dict(selected=original,repeats=repeats,graph_source=repeats[0]['case_path'],contract_sha256=_CONTRACT_SHA))
        assert time.time()<timing['end']
        write(a.output/'complete.json',dict(status='bounded_search_complete',records=records,selected=selected,common_primary_key=common_key,hourly_reports=reports,clock=timing,contract_sha256=_CONTRACT_SHA,scientific_identity=science(c)))
    except Exception as exc:write(a.output/'complete.json',dict(status='incomplete_or_fatal_stop',error_type=type(exc).__name__,error=str(exc),records=records,best=best,clock=timing,contract_sha256=_CONTRACT_SHA));raise

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('prepare','smoke','search','evaluate','render'));p.add_argument('--contract',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--request',type=Path);p.add_argument('--approval',type=Path);p.add_argument('--approval-sha256');a=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only'
    c,objs=verify(a.contract)
    if a.stage=='prepare':a.output.mkdir(parents=True,exist_ok=False);write(a.output/'prepared.json',dict(status='source_contract_verified_no_solves',contract_sha256=sha(a.contract),clock=clock(c,time.time()),free=list(FREE),lanes=list(LANES)));return
    if a.stage=='render':render_saved(a,c)
    elif a.stage=='evaluate':evaluate(a,c,objs)
    else:controller(a,c,objs)
if __name__=='__main__':main()
