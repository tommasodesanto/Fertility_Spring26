#!/usr/bin/env python3
"""Ten-case exploratory early-fertility line probe, not calibrated SMM.

No model, weighting, normalization or scientific-gate changes. This external
controller dispatches only the byte-pinned night driver's existing evaluator.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,os,sys,time,traceback
from pathlib import Path
SCHEMA='e5f_early_frontier_v1'
CONTRACT_SHA='c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d'
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')

def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1<<20),b''):h.update(part)
    return h.hexdigest()
def write(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.tmp');q.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n');q.replace(p)
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def csvwrite(path,rows):
    assert rows
    with Path(path).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def verify_plan(path,pin):
    assert sha(path)==pin,'Plan fingerprint changed';p=read(path);assert p['schema']==SCHEMA
    for name in ('frontier_driver','night_driver','contract'):
        item=p[name];assert sha(item['path'])==item['sha256'],name
    assert p['frontier_driver']['sha256']==sha(__file__)
    assert p['contract']['sha256']==CONTRACT_SHA
    assert p['workers']==10 and p['max_cases']==10 and p['objective_cap_seconds']==1800 and p['total_seconds']==2400
    assert (p['not_before'],p['objective_end'],p['hard_end'])==(1790593200,1790595600,1790595900)
    return p

def timing(p,now):
    assert now>=p['not_before'],'Do not start before main search cutoff'
    end=min(now+p['total_seconds'],p['hard_end']);objective_end=min(p['objective_end'],end-300)
    assert now+1800<objective_end,'Insufficient time for one full parallel batch and reporting'
    return dict(start=now,end=end,objective_end=objective_end,latest_dispatch=objective_end-1800)

def completion_proof(p,path,pin,now):
    if path is None:return None
    assert pin and sha(path)==pin
    proof=read(path);complete=read(p['main_complete'])
    assert proof['job_id']==p['main_job_id'] and proof['state']=='COMPLETED' and proof['exit_code']=='0:0'
    assert 0<=now-proof['verified_epoch']<=60
    assert proof['complete_sha256']==sha(p['main_complete'])
    assert complete['status']=='bounded_search_complete' and complete['contract_sha256']==CONTRACT_SHA
    return dict(method='wrapper_verified_ended_main',proof=proof)

def wait_resources(p,clock,out,proof_path=None,proof_pin=None):
    while time.time()<clock['latest_dispatch']:
        now=time.time();assert now>=p['not_before']
        heartbeat=Path(p['main_heartbeat'])
        if heartbeat.exists():
            h=read(heartbeat)
            if 0<=now-h['epoch']<=30 and isinstance(h.get('active'),int) and 0<=h['active']<=8:
                return dict(method='fresh_main_heartbeat',active=h['active'],epoch=h['epoch'],maximum_simultaneous_household_workers=19)
        if proof_path:return completion_proof(p,proof_path,proof_pin,now)
        write(out/'heartbeat.json',dict(stage='waiting_for_main_at_most8',epoch=now,latest_dispatch=clock['latest_dispatch']))
        time.sleep(5)
    raise TimeoutError('No verified capacity before full-case dispatch cutoff')

def compact_case(row,c,objs,night):
    source=Path(row['case_path']);req=read(row['request_path']);receipt=read(source/'receipt.json');side=read(source.parent/'success.json');lane=row['lane']
    assert row['point']==req['point']==receipt['point'] and lane==req['lane']==receipt['lane']
    assert sha(source/'receipt.json')==row['receipt_sha256']==side['receipt_sha256']
    assert side['context']==req['context'] and receipt['contract_sha256']==req['contract_sha256']==CONTRACT_SHA
    assert receipt['scientific_identity']==night.science(c) and receipt['scientific_candidate_id']==night.identity(c,lane,row['point'])
    assert receipt['normalization_inputs']==c['normalization'] and receipt['source_manifest_sha256']==c['source_manifest']['sha256']
    assert receipt['target_system_sha256']==c['lanes'][lane]['objective']['sha256'] and receipt['target_weight_fingerprint']==c['lanes'][lane]['target_weight_fingerprint']
    fits=night.table(source/'target_fit.csv','moment');targets={x['restriction_id']:x for x in objs[lane]['target_rows']};primary={x['restriction_id']:x for x in objs['primary']['target_rows']}
    assert len(fits)==14 and fits.keys()==targets.keys();loss=0.;primary_loss=0.
    for k,x in fits.items():
        target=float(x['target']);model=float(x['model']);gap=float(x['gap']);weight=targets[k]['actual_weight']
        assert target==targets[k]['target'] and math.isfinite(model) and math.isclose(gap,model-target,rel_tol=1e-12,abs_tol=1e-12)
        if weight is None:assert x['weight']==x['loss_contribution']==''
        else:
            assert float(x['weight'])==weight and math.isclose(float(x['loss_contribution']),weight*gap*gap,rel_tol=1e-12,abs_tol=1e-12);loss+=float(x['loss_contribution'])
        primary_loss+=(primary[k]['actual_weight'] or 0)*gap*gap
    assert side['loss']==receipt['loss'] and math.isclose(row['loss'],side['loss'],rel_tol=1e-12,abs_tol=1e-12) and math.isclose(loss,row['loss'],rel_tol=1e-12,abs_tol=1e-12)
    assert math.isclose(primary_loss,row['primary_rescore'],rel_tol=1e-12,abs_tol=1e-12)
    early=fits['early_fertility']
    return dict(record=row,early_model=float(early['model']),early_target=float(early['target']),early_gap=float(early['gap']),primary_loss=primary_loss)

def snapshot_anchors(p,c,objs,night,out,deadline):
    snapshot=Path(p['main_checkpoint']).read_bytes();target=out/'main_snapshot.json';target.write_bytes(snapshot);data=json.loads(snapshot)
    assert data['contract_sha256']==CONTRACT_SHA
    candidates=list(data['records'])+list(data.get('best',{}).values())+[data.get('common_primary_best')]
    unique={};audited=[]
    for row in candidates:
        if not row or row.get('status')!='success':continue
        unique.setdefault(row['case_path'],row)
    for row in unique.values():
        if time.time()>=deadline:raise TimeoutError('Snapshot audit consumed dispatch reserve')
        audited.append(compact_case(row,c,objs,night))
    assert audited,'No authenticated successful cases'
    primary=min(audited,key=lambda x:(x['primary_loss'],x['record']['case_path']))
    early=min(audited,key=lambda x:(abs(x['early_gap']),x['primary_loss'],x['record']['case_path']))
    anchors={'primary_best':primary,'early_best':early}
    for item in anchors.values():
        row=item['record'];source=Path(row['case_path']);actual=night.validate(source.parent,c,objs,read(row['request_path']))
        assert all(row[key]==value for key,value in actual.items()),'Chosen anchor full audit differs'
        item['table_pins']={name:sha(source/name) for name in ('target_fit.csv','parameters.csv','receipt.json')}
    write(out/'snapshot_audit.json',dict(snapshot_sha256=sha(target),successful_cases=len(audited),anchors=anchors,pool=audited))
    return anchors

def requests_for(anchors,c,obj):
    base=anchors['primary_best']['record']['point'];rows=[]
    for name in ('primary_best','early_best'):rows.append(dict(id='anchor_'+name,kind='anchor',anchor=name,point=dict(anchors[name]['record']['point']),label='Exact physical replay under primary weights'))
    recipes=[('continuation_x2',{'kappa_fert_continuation':base['kappa_fert_continuation']*2}),('continuation_x4',{'kappa_fert_continuation':base['kappa_fert_continuation']*4}),('initial_x_half',{'kappa_fert':base['kappa_fert']*.5}),('initial_x_quarter',{'kappa_fert':base['kappa_fert']*.25}),('paired_half_double',{'kappa_fert':base['kappa_fert']*.5,'kappa_fert_continuation':base['kappa_fert_continuation']*2}),('paired_quarter_quadruple',{'kappa_fert':base['kappa_fert']*.25,'kappa_fert_continuation':base['kappa_fert_continuation']*4})]
    strong=dict(recipes[-1][1],first_birth_fixed_cost=0.)
    recipes.extend([('boundary_zero_first_birth_cost',strong),('boundary_zero_first_birth_cost_and_loading',dict(strong,delta_alpha_jump=0.))])
    bounds={r['parameter']:(r['lower'],r['upper']) for r in obj['parameter_restrictions']}
    seen={json.dumps(r['point'],sort_keys=True):r['id'] for r in rows}
    for name,overrides in recipes:
        point=dict(base);changes={}
        for k,requested in overrides.items():
            lo,hi=bounds[k];applied=min(hi,max(lo,requested));point[k]=applied;changes[k]=dict(requested=requested,applied=applied,clipped=applied!=requested)
        assert all(math.isfinite(point[k]) and lo<=point[k]<=hi for k,(lo,hi) in bounds.items())
        key=json.dumps(point,sort_keys=True);duplicate=seen.get(key);seen.setdefault(key,name)
        rows.append(dict(id=name,kind='proposal',anchor='primary_best',point=point,label='Parameter-bound mechanism shutdown DIAGNOSTIC' if name.startswith('boundary_') else 'Exploratory timing-taste line probe',overrides=overrides,parameter_changes=changes,duplicate_of=duplicate))
    assert len(rows)==10
    return rows

def nondominated(rows):
    return [r['id'] for r in rows if not any(abs(s['early_gap'])<=abs(r['early_gap']) and s['other_moments_primary_loss']<=r['other_moments_primary_loss'] and (abs(s['early_gap'])<abs(r['early_gap']) or s['other_moments_primary_loss']<r['other_moments_primary_loss']) for s in rows if s['id']!=r['id'])]

def reports(out,records,anchors,c,objs,night,verified):
    pool=[];fits_all=[];params_all=[];initial_loss=anchors['primary_best']['primary_loss'];initial_fits=night.table(Path(anchors['primary_best']['record']['case_path'])/'target_fit.csv','moment')
    items=[('initial_primary_best','reference',anchors['primary_best']['record']),('initial_early_best','reference',anchors['early_best']['record'])]+[(r['id'],r['kind'],r) for r in records if r['status']=='success']
    weights={r['restriction_id']:r['actual_weight'] for r in objs['primary']['target_rows']}
    for ident,kind,record in items:
        fits=night.table(Path(record['case_path'])/'target_fit.csv','moment');early=fits['early_fertility'];loss=sum((weights[k] or 0)*float(r['gap'])**2 for k,r in fits.items())
        pool.append(dict(id=ident,kind=kind,source_lane=record['lane'],primary_loss=loss,other_moments_primary_loss=loss-weights['early_fertility']*float(early['gap'])**2,loss_minus_initial=loss-initial_loss,early_target=float(early['target']),early_model=float(early['model']),early_gap=float(early['gap']),absolute_early_gap=abs(float(early['gap'])),accepted_for_diagnostic=verified,case_path=record['case_path']))
        for name,row in fits.items():
            weight=weights[name];contribution='' if weight is None else weight*float(row['gap'])**2
            fits_all.append(dict(case_id=ident,kind=kind,source_lane=record['lane'],moment=name,target=row['target'],model=row['model'],gap=row['gap'],model_minus_main_anchor=float(row['model'])-float(initial_fits[name]['model']),primary_weight='' if weight is None else weight,primary_loss_contribution=contribution,share_of_primary_loss='' if weight is None or loss==0 else contribution/loss,role=row.get('role','')))
        for name,row in night.table(Path(record['case_path'])/'parameters.csv','parameter').items():params_all.append(dict(case_id=ident,**row))
    csvwrite(out/'full_target_comparison.csv',fits_all);csvwrite(out/'full_parameter_comparison.csv',params_all)
    eligible=pool if verified else []
    frontier=nondominated(eligible);ordered=sorted(pool,key=lambda r:(r['absolute_early_gap'],r['primary_loss']))
    csvwrite(out/'early_gap_ranking.csv',[dict(r,nondominated_among_observed_points=r['id'] in frontier) for r in ordered])
    write(out/'tradeoff.json',dict(accepted=verified,nondominated_ids=frontier,criterion='Absolute early-fertility gap and unchanged primary-weighted loss from all OTHER moments among these observed points only',not_calibrated_SMM=True,identification_or_global_optimum_claim=False))
    return dict(ranking=ordered,nondominated_ids=frontier)

def run(p,c,objs,night,supervisor,out,clock,anchors):
    proposals=requests_for(anchors,c,objs['primary']);write(out/'declared_cases.json',proposals)
    records=[dict(spec,status='duplicate_unrun',lane='primary') for spec in proposals if spec.get('duplicate_of')];dispatch=[x for x in proposals if not x.get('duplicate_of')];fatal=False;started=time.monotonic()
    def heartbeat(**extra):write(out/'heartbeat.json',dict(stage='ten_case_batch',epoch=time.time(),elapsed_seconds=time.monotonic()-started,completed=len(records),**extra))
    def launch(spec,batch_deadline):
        assert time.time()<=clock['latest_dispatch'],'No short late objective attempts'
        context=dict(candidate_id=spec['id'],stage='initial',contract_sha256=CONTRACT_SHA,source_sha256=p['night_driver']['sha256'],target_sha256=c['lanes']['primary']['objective']['sha256'],point_sha256=night.canon(spec['point']))
        req=dict(id=spec['id'],lane='primary',point=spec['point'],context=context,stage='initial',contract_sha256=CONTRACT_SHA,scientific_candidate_id=night.identity(c,'primary',spec['point']),normalization_inputs=c['normalization'],controller_pid=os.getpid(),deadline_epoch=min(time.time()+1800,batch_deadline),graphs=True)
        path=out/(spec['id']+'.request.json');write(path,req);spec.update(request=req,request_path=str(path))
        env=os.environ.copy();env.update({key:'1' for key in THREADS});env['EXPECTED_E5F_NIGHT_SHA256']=CONTRACT_SHA
        return supervisor.ManagedProcess([sys.executable,p['night_driver']['path'],'--stage','evaluate','--contract',p['contract']['path'],'--output',str(out/spec['id']),'--request',str(path)],out/(spec['id']+'.log'),req['deadline_epoch'],env)
    def finish(spec,proc,code):
        nonlocal fatal
        try:
            night.verify(p['contract']['path']);status,data,error=night.classify(out/spec['id'],c,objs,spec['request'],proc,code)
        except Exception as exc:status,data,error='fatal',{},str(exc)
        row=dict(id=spec['id'],kind=spec['kind'],anchor=spec['anchor'],label=spec['label'],point=spec['point'],lane='primary',status=status,error=error,request_path=spec['request_path'],**data)
        records.append(row);fatal|=status=='fatal';row['halt_new_dispatch']=fatal
        write(out/'latest_completed.json',row);write(out/'checkpoint.json',dict(records=records,clock=clock,contract_sha256=CONTRACT_SHA,max_cases=10,attempted=len([x for x in proposals if 'request' in x])));return row
    result=supervisor.run_batch(dispatch,workers=10,deadline=clock['objective_end'],launch=launch,finish=finish,heartbeat=heartbeat,poll_seconds=1,allowed_statuses={'success','inadmissible','censored_timeout','censored_late_completion'},guard=lambda:'fatal_stop' if fatal else None)
    fatal|=any(row['status'] in ('failed','fatal') for row in result['results'])
    anchor_checks={}
    for name in ('primary_best','early_best'):
        row=next((r for r in records if r['id']=='anchor_'+name),None)
        try:
            assert row and row['status']=='success'
            anchor=anchors[name];source=Path(anchor['record']['case_path'])
            assert all(sha(source/k)==v for k,v in anchor['table_pins'].items())
            anchor_checks[name]=night.compare_anchor(source,row['case_path'])
        except Exception as exc:anchor_checks[name]=dict(failed=True,error=str(exc))
    verified=not fatal and result['complete'] and all(not x.get('failed',False) for x in anchor_checks.values())
    summary=reports(out,records,anchors,c,objs,night,verified)
    assert time.time()<clock['end'],'Report exceeded fixed total budget'
    complete=dict(status='verified_exploratory_line_probe' if verified else 'unverified_exploratory_line_probe',accepted=verified,anchor_checks=anchor_checks,records=records,batch=result,summary=summary,clock=clock,contract_sha256=CONTRACT_SHA,weight_changes=False,not_calibrated_SMM=True,identification_or_reachability_claim=False)
    write(out/'complete.json',complete);return complete

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--plan-sha256',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--completion-proof',type=Path);ap.add_argument('--completion-proof-sha256');a=ap.parse_args()
    assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit() and not sys.flags.optimize
    p=verify_plan(a.plan,a.plan_sha256);assert int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=10
    clock=timing(p,time.time());a.output.mkdir(parents=True,exist_ok=False);write(a.output/'clock.json',clock)
    try:
        evidence=wait_resources(p,clock,a.output,a.completion_proof,a.completion_proof_sha256);write(a.output/'resource_gate.json',evidence)
        os.environ['EXPECTED_E5F_NIGHT_SHA256']=CONTRACT_SHA
        night=load('frontier_verified_night_driver',p['night_driver']['path']);c,objs=night.verify(p['contract']['path'])
        assert c['files']['driver']==p['night_driver']
        anchors=snapshot_anchors(p,c,objs,night,a.output,clock['latest_dispatch']-10)
        sys.path.insert(0,str(Path(c['files']['recovery_search']['path']).parent));supervisor=load('frontier_owned_supervision',c['files']['recovery_search']['path'])
        result=run(p,c,objs,night,supervisor,a.output,clock,anchors)
        print(json.dumps(dict(status=result['status'],cases=len(result['records']),accepted=result['accepted'],not_calibrated_SMM=True)))
        if not result['accepted']:raise RuntimeError('Both exact anchor checks are required; diagnostic is unverified')
    except BaseException as exc:
        if not (a.output/'complete.json').exists():write(a.output/'complete.json',dict(status='unverified_setup_or_execution_failure',accepted=False,error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),clock=clock))
        raise
if __name__=='__main__':main()
