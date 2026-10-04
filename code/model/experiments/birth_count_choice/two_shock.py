#!/usr/bin/env python3
"""Two unanticipated permanent preference levels on the frozen Estate-A baseline.

Empirical scalar fits use local fertility index one; smoke runs fixed baseline levels. Only the accepted 2007 vintage's
first two dates are implemented before the 2015 surprise. Import/preflight solve
nothing; native execution requires authenticated staged sources and a passed
same-source two-stage smoke.
"""
from __future__ import annotations
import argparse
import contextlib
import copy
import csv
import gzip
import hashlib
import importlib
import json
import math
from pathlib import Path
import pickle
import sys
import time
import threading
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASELINE_PSI=.17892072066041628
SCHEMA='estate_a_two_unanticipated_psi_v1'
TARGETS=[1.974875,1.861,1.755375,1.64575]
STAGES=[dict(start_year=2007,accepted_periods=2,target_index=1),dict(start_year=2015,accepted_periods=0,target_index=1)]
SMOKE_PROTOCOL='native_two_stage_execution_only_v1'
GATES=dict(market_tolerance=2e-4,fiscal_tolerance=2e-5,final_reproduction_tolerance=1e-10,
 stationary_renewal_tolerance=1e-6,terminal_tolerance=1e-3,raw_queue_relative_tolerance=1e-3,horizon_relative_tolerance=1e-3)
def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    def plain(x):
        if isinstance(x, np.ndarray): return plain(x.tolist())
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, dict): return {k: plain(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [plain(v) for v in x]
        if isinstance(x, float) and not math.isfinite(x): return None
        if isinstance(x, (str, int, float, bool)) or x is None: return x
        # Fit receipts retain scalar/root evidence; native arrays and dated
        # states are retained by their own exported checkpoints, never JSON
        # reconstructed here.
        raise TypeError('Nonserializable native object in receipt: '+type(x).__name__)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(plain(value), indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=sha(path))


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Exact file pin required')
    path = Path(item['path']); require(path.is_file() and sha(path) == item['sha256'], 'Missing or changed pin: '+str(path))
    return path




def digest(value): return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def state_hash(state,queue_values=lambda x:x):
    values=[]
    for key in ('g_pre','scheduled_entries','scheduled_raw_entries'):
        x=np.asarray(getattr(state,key) if key=='g_pre' else queue_values(getattr(state,key)))
        require(x.dtype.kind in 'biuf' and np.isfinite(x).all() and (x>=0).all(),'Finite nonnegative state required: '+key)
        values.append(dict(name=key,shape=list(x.shape),dtype=str(x.dtype),sha256=hashlib.sha256(x.tobytes()).hexdigest()))
    return digest(values)


def lightweight_original(base):
    directory=pinned(base['source_files']['original_estimator']).parent
    sys.path.insert(0,str(directory))
    modules=[importlib.import_module(n) for n in ('run_e5f_preference_estimation','e5f_preference_shock_fit')]
    for key,module in zip(('original_estimator','original_fitter'),modules):
        require(Path(module.__file__).resolve()==pinned(base['source_files'][key]).resolve(),'Cached original source path differs: '+key)
    return modules


def empirical_controls(base):
    controls=dict(gates=copy.deepcopy(GATES),seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),
        fit=copy.deepcopy(base['fit']),path=copy.deepcopy(base['path']),endpoint=copy.deepcopy(base['endpoint']),
        budget=copy.deepcopy(base['budget']),horizons=[24,32],smoke_seed_endpoint_padding=False)
    controls['fit']['max_evaluations']=12;controls['endpoint']['max_evaluations']=48
    controls['path']['max_evaluations']=12
    controls['budget'].update(total_seconds=21480,maximum_policy_calls=20000,endpoint_seconds=1800)
    return controls


def smoke_controls(empirical):
    controls=copy.deepcopy(empirical)
    controls['horizons']=[6];controls['smoke_seed_endpoint_padding']=True
    controls['gates'].update(market_tolerance=.05,fiscal_tolerance=.005)
    controls['path']['max_evaluations']=3;controls['endpoint']['max_evaluations']=8
    controls['budget'].update(total_seconds=1680,maximum_policy_calls=400,endpoint_seconds=240,
        candidate_seconds=360,path_seconds=240,seed_seconds=360,mapping_seconds=180,render_seconds=120)
    return controls


def numerical_controls(plan):
    return copy.deepcopy({k:plan[k] for k in ('gates','seed','fit','endpoint','path','budget','horizons','smoke_seed_endpoint_padding')})


def prepare_manifest(base_plan_pin,source_pins=None,*,smoke=False,starts=(.14736308634876963,.12),legacy_source_overlay=None):
    base=json.loads(pinned(base_plan_pin).read_text())
    sources=source_pins or dict(two_shock_driver=pin(__file__),two_shock_runtime=pin(HERE/'two_shock_runtime.py'))
    empirical=empirical_controls(base)
    plan=dict(schema=SCHEMA,kind='two_unanticipated_permanent',mode='diagnostic',smoke=bool(smoke),
        baseline_psi=BASELINE_PSI,psi_bound_ratios=[.01,2.],stages=copy.deepcopy(STAGES),
        rows=copy.deepcopy(base['target_contract']['rows']),weights=[0,1,0,1],target_contract=copy.deepcopy(base['target_contract']),
        base_plan=base_plan_pin,source_pins=sources,identity=copy.deepcopy(base['identity']),
        standard_plot_names=base['standard_plot_names'],execution_enabled=True,reviewed=True,
        smoke_targets=None,smoke_protocol=SMOKE_PROTOCOL,empirical_controls=empirical,
        disclosure='Experimental two surprises; unchanged one-birth Estate-A economics. Provisional estate closure; no production certification.')
    plan.update(smoke_controls(empirical) if smoke else copy.deepcopy(empirical))
    if legacy_source_overlay is not None: plan['legacy_source_overlay']=copy.deepcopy(legacy_source_overlay)
    plan['stage_starts']=[dict(initial=BASELINE_PSI if smoke else float(v),bounds=[BASELINE_PSI*.01,BASELINE_PSI*2]) for v in starts]
    return plan


def fingerprints(plan):
    controls=numerical_controls(plan)
    empirical=plan.get('empirical_controls',controls)
    return dict(source=digest(plan['source_pins']),contract=digest(dict({k:plan[k] for k in
        ('schema','kind','identity','baseline_psi','psi_bound_ratios','stages','rows','weights','target_contract')},legacy_source_overlay=plan.get('legacy_source_overlay'))),
        controls=digest(controls),empirical_controls=digest(empirical),
        smoke_controls=digest(dict(protocol=SMOKE_PROTOCOL,numerical_controls=smoke_controls(empirical),
            fixed_levels=[BASELINE_PSI]*2,scalar_optimizer=False,derivative_probes=False)))


def validate_mode_controls(plan):
    # This check runs before any source import or native work. Smoke never sets
    # empirical controls: its receipt authenticates this original contract.
    empirical=plan.get('empirical_controls')
    require(isinstance(empirical,dict),'Explicit original empirical-control fingerprint required')
    require(plan.get('smoke_protocol')==SMOKE_PROTOCOL,'Explicit execution-only smoke protocol required')
    require(empirical['gates']==GATES and empirical['seed']==dict(horizon=12,perturbed_date=5,log_step=1e-5), 'Original empirical gates/seed required')
    require(empirical['horizons']==[24,32] and empirical['smoke_seed_endpoint_padding'] is False,'Original empirical horizons/padding required')
    b=empirical['budget']
    require(b['total_seconds']==21480 and b['maximum_policy_calls']==20000 and b['endpoint_seconds']==1800,'Original empirical budget required')
    require(empirical['fit']['max_evaluations']==12 and empirical['fit']['fertility_tolerance']==.005 and
        empirical['endpoint']['max_evaluations']==48 and empirical['path']['max_evaluations']==12,'Original empirical iteration controls required')
    for name in ('path','endpoint'):
        require(empirical[name]['price_bound_ratios']==[.05,20.] and empirical[name]['max_log_step']==.15,'Original empirical numerical domains required')
    require(empirical['path']['pension_bound_ratios']==[.05,20.] and empirical['endpoint']['slope']==1.,'Original empirical fiscal domains/endpoint slope required')
    require(numerical_controls(plan)==(smoke_controls(empirical) if plan['smoke'] else empirical),
        'Exact smoke execution controls or original empirical controls required')
    require(plan.get('smoke_targets') is None,'Synthetic fit targets are forbidden in execution-only smoke')
    if plan['smoke']:require(all(x['initial']==BASELINE_PSI for x in plan['stage_starts']),'Fixed baseline smoke levels required')


def preflight(plan):
    require(plan.get('schema')==SCHEMA and plan.get('kind')=='two_unanticipated_permanent','Two-surprise schema required')
    require(plan.get('mode')=='diagnostic' and type(plan.get('smoke')) is bool,'Explicit diagnostic smoke/run mode required')
    require(plan.get('baseline_psi')==BASELINE_PSI and plan.get('psi_bound_ratios')==[.01,2.],'Original absolute preference bounds required')
    require(plan.get('weights')==[0,1,0,1] and [r['target'] for r in plan['rows']]==TARGETS,'Full four targets and weights required')
    require(plan.get('stages')==STAGES,'Both local target indices must be one; only stage1 advances')
    validate_mode_controls(plan)
    b=plan['budget']
    require(all(s['bounds']==[BASELINE_PSI*.01,BASELINE_PSI*2] and s['bounds'][0]<s['initial']<s['bounds'][1] for s in plan['stage_starts']),'Same absolute bounds and interior starts required')
    require(plan.get('reviewed') is True and plan.get('execution_enabled') is True,'Reviewed execution-enabled package required')
    if plan.get('legacy_source_overlay') is not None:
        require(set(plan['legacy_source_overlay'])=={'helper','manifest'},'Exact pinned legacy overlay fields required')
        for item in plan['legacy_source_overlay'].values(): pinned(item)
    base=json.loads(pinned(plan['base_plan']).read_text())
    require(plan['empirical_controls']==empirical_controls(base),'Original empirical-control contract changed')
    require(base['initial_psi']==BASELINE_PSI and base['identity']==plan['identity'] and base['target_contract']==plan['target_contract'], 'Original baseline/target identity required')
    for item in base['source_files'].values(): pinned(item)
    frozen_root=pinned(base['source_files']['runtime']).parents[4]
    for name,value in base['identity']['source_pins'].items():
        path=(frozen_root/name).resolve();require(path.is_relative_to(frozen_root) and sha(path)==value,'Frozen source changed: '+name)
    for item in plan['source_pins'].values(): pinned(item)
    require(plan['source_pins']['two_shock_driver']==pin(__file__),'Executing two-shock driver path/hash differs')
    require(plan['source_pins']['two_shock_runtime']==pin(HERE/'two_shock_runtime.py'),'Concrete runtime path/hash differs')
    require(ROOT==frozen_root,'Driver must execute inside the authenticated frozen source root')
    handoff=json.loads(pinned(base['handoff']).read_text());require(handoff['source_pins']==base['identity']['source_pins'],'Handoff source identity differs')
    case=(frozen_root/handoff['saved_case']).resolve();require(case.is_relative_to(frozen_root),'Saved case escapes frozen root')
    for name,value in {**handoff['saved_files'],**handoff.get('standard_plot_pins',{})}.items():
        path=(case/name).resolve();require(path.is_relative_to(case) and sha(path)==value,'Saved input/plot differs: '+name)
    require(plan['standard_plot_names']==base['standard_plot_names'],'Original standard plot names required')
    require(len(plan['standard_plot_names'])==len(set(plan['standard_plot_names']))==17,'Exact 17 standard names required')
    require(all(type(v) in (int,float) and math.isfinite(v) and v>0 for v in b.values()),'Every budget positive and finite')
    original,_=lightweight_original(base)
    contract=original.target_contract(pinned(base['target_contract']['blocks']),pinned(base['target_contract']['annual']))
    require(contract==plan['target_contract'] and contract['rows']==plan['rows'],'Pinned annual-builder provenance differs')
    return dict(status='PASS',schema=SCHEMA,native_calls=0,scientific_validation=False,production_ready=False,
        fingerprints=fingerprints(plan),horizons=plan['horizons'],total_seconds=b['total_seconds'],policy_call_stop_cap=b['maximum_policy_calls'])


def queue_array(queue):
    if hasattr(queue,'due_in_16') and hasattr(queue,'due_in_20'):
        return np.asarray(list(queue.due_in_16)+list(queue.due_in_20),float)
    return np.asarray(queue)


def validate_state_export(exported,identity,stage):
    require(exported.get('calendar_year')==2023 and exported.get('period_index')==(4 if stage==0 else 2) and
        exported.get('exact_native_state') is True and exported.get('reconstructed_or_rescaled') is False and
        exported.get('queue_lags')==[16,20] and exported.get('forecast_and_continuation_saved') is True,
        'Actual smoke 2023 export evidence required')
    with gzip.open(pinned(dict(path=exported['path'],sha256=exported['sha256'])),'rb') as stream:packet=pickle.load(stream)
    require(packet.get('schema')=='current_floor_actual_2023_v1' and packet.get('calendar_year')==2023 and
        packet.get('period')==(4 if stage==0 else 2) and packet.get('continuation_calendar_year')==2027 and
        packet.get('reference_identity')==identity,'Actual exported checkpoint clock/identity differs')
    state=packet['initial_state'];actual_hash=state_hash(state,queue_array)
    require(actual_hash==exported.get('state_sha256'),'Actual exported distribution/queue hash differs')
    for name in ('scheduled_entries','scheduled_raw_entries'):
        require(np.array_equal(queue_array(packet[name]),queue_array(getattr(state,name))), 'Actual exported queue copy differs')
    require(float(np.asarray(state.g_pre).sum())==packet['initial_population']>0,'Actual exported physical population differs')
    for key in ('forecast_prices','forecast_pensions','forecast_psi','b_grid'):
        values=np.asarray(packet[key]);require(values.size>0 and np.isfinite(values).all(),'Actual exported finite continuation required: '+key)
    require(np.asarray(packet['current_2023_V']).shape==np.asarray(state.g_pre).shape and
        np.asarray(packet['continuation_V']).shape==np.asarray(state.g_pre).shape,'Actual exported value shape differs')
    occupied=np.asarray(state.g_pre)>0
    for key in ('current_2023_V','continuation_V'):
        require(np.isfinite(np.asarray(packet[key])[occupied]).all(),'Actual exported occupied value integrity differs')
    require(np.asarray(packet['terminal_V']).shape==np.asarray(state.g_pre).shape,'Actual terminal value shape differs')
    lengths=[len(packet[k]) for k in ('forecast_prices','forecast_pensions','forecast_psi')]
    require(len(set(lengths))==1 and lengths[0]==6-packet['period'] and
        np.all(np.asarray(packet['forecast_psi'])==BASELINE_PSI),'Actual reduced forecast/baseline level differs')
    return packet


def validate_plot_receipt(item,names):
    path=pinned(item);plots=json.loads(path.read_text())
    require(plots.get('sampled_dates'),'Actual native dated plots required')
    for dated in plots['sampled_dates']:
        require(set(dated['plots'])==set(names),'Exact 17 smoke plot names required')
        directory=path.parent/f"date_{dated['period']:03d}"/'standard_diagnostics'
        for name,value in dated['plots'].items():require((directory/name).is_file() and sha(directory/name)==value,'Actual smoke plot hash differs')


def validate_smoke_pin(plan,item):
    validate_mode_controls(plan)
    receipt=json.loads(pinned(item).read_text())
    require(receipt.get('status')=='execution_passed' and receipt.get('schema')==SCHEMA and receipt.get('smoke') is True and
        all(receipt.get(k) is False for k in ('scientific_validation','production_ready','empirical_fitted')),
        'Passed actual execution-only two-stage smoke required')
    smoke=json.loads(pinned(receipt['smoke_manifest_pin']).read_text());preflight(smoke)
    require(smoke['smoke'] is True,'Pinned actual smoke manifest required')
    expected=fingerprints(plan);observed=receipt.get('fingerprints',{})
    require(observed==fingerprints(smoke) and all(observed.get(k)==expected[k] for k in
        ('source','contract','empirical_controls','smoke_controls')),'Smoke source/contract/empirical fingerprints differ')
    require(receipt.get('semantic_gates')==dict(fixed_baseline_candidates=True,scalar_optimizer=False,derivative_probes=False,
        nonanticipating_prefix=True,both_queues=True,fresh_replays=True,exact_2023=True,standard_diagnostics=True,
        empirical_horizon_gates_tested=False,empirical_forecast_slice_tested=False),'Execution-only smoke semantic gates required')
    reference=receipt['reference_evidence'];proof=json.loads(pinned(reference['reconstruction_receipt']).read_text())
    checkpoint=pinned(reference['reference_checkpoint'])
    require(reference.get('accounting_valid') is True and type(reference.get('policy_calls')) is int and reference['policy_calls']>0 and
        proof.get('status')=='passed' and proof.get('policy_calls')==reference['policy_calls'] and
        proof.get('checkpoint_sha256')==sha(checkpoint),'Actual fresh reference reconstruction/repeat required')
    calls=receipt.get('actual_policy_calls');require(type(calls) is int and 0<calls<=400 and
        receipt.get('native_actual_policy_calls')==calls,'Actual smoke native call cap/accounting differs')
    for evidence in receipt['evidence_pins']:pinned(evidence)
    require(len(receipt.get('stage_seed_pins',[]))==len(receipt.get('selected_candidate_pins',[]))==2,'Both concrete smoke seeds/candidates required')
    for stage,(seed_pin,candidate_pin) in enumerate(zip(receipt['stage_seed_pins'],receipt['selected_candidate_pins'])):
        seed=json.loads(pinned(seed_pin).read_text());candidate=json.loads(pinned(candidate_pin).read_text())
        require(seed['mapping_count']==5 and seed['horizon']==12 and
            (stage==0 or seed.get('nonstationary_inherited_state') is True),'Actual stage-bound five-map seeds required')
        require(len(seed.get('source_evidence',[]))==5,'Five actual native seed map pins required')
        for proof in seed['source_evidence']:
            record=json.loads(pinned(proof).read_text())
            require(record.get('accounting_valid') is True and bool(record.get('gates')) and all(record['gates'].values()) and
                len(record['rows'])==12 and all(row['calendar_year']==STAGES[stage]['start_year']+4*i for i,row in enumerate(record['rows'])) and
                record['two_shock_provenance']['inherited_state_sha256']==seed['identity']['inherited_state_sha256'],
                'Actual smoke seed native map clock/state/accounting differs')
        require(candidate.get('status')=='execution_passed' and candidate.get('certified') is False and
            candidate['stage']==stage and candidate['psi']==BASELINE_PSI and candidate['horizon']==6 and
            candidate['fresh_native_replay'] is True and candidate['accounting_valid'] is True and
            candidate['inherited_state_sha256']==seed['identity']['inherited_state_sha256'],'Actual smoke fixed candidate/replay evidence differs')
        for evidence in candidate['source_evidence']:pinned(evidence)
        terminal=json.loads(pinned(candidate['terminal_diagnostics_pin']).read_text())
        require(type(candidate.get('terminal_pass')) is bool and candidate.get('terminal_diagnostics_gating') is False and
            terminal.get('all_checks_pass') is candidate['terminal_pass'] and receipt.get('terminal_diagnostics_gating') is False and
            receipt.get('terminal_passes_by_stage',{}).get(str(stage+1)) is candidate['terminal_pass'],
            'Actual nongating smoke terminal diagnostics differ')
        root=json.loads(pinned(candidate['root_pin']).read_text())
        require(root.get('converged') is True and bool(root.get('gates')) and all(root['gates'].values()) and
            root.get('final_reproduction_max_abs',float('inf'))<=GATES['final_reproduction_tolerance'], 'Actual smoke root fresh replay failed')
        record=json.loads(pinned(candidate['final_mapping_pin']).read_text())
        require(2<=candidate['path_evaluations']<=3 and record.get('accounting_valid') is True and
            bool(record.get('gates')) and all(record['gates'].values()) and len(record['rows'])==6 and
            all(row['calendar_year']==STAGES[stage]['start_year']+4*i for i,row in enumerate(record['rows'])) and
            record['two_shock_provenance']['inherited_state_sha256']==candidate['inherited_state_sha256'] and
            record['two_shock_provenance']['start_year']==STAGES[stage]['start_year'], 'Actual smoke final native mapping differs')
        for key,tolerance in (('market_residual',.05),('fiscal_residual',.005)):
            require(len(record[key])==6 and all(math.isfinite(x) and abs(x)<=tolerance for x in record[key]), 'Actual smoke final physical residual differs')
        final=root['final'];require(final['mapping_valid'] is True and len(final['prices'])==len(final['fiscal_values'])==6 and
            all(abs(row['asset_price']-q)<=1e-10 and abs(row['pension_period_units']-b)<=1e-10
                for row,q,b in zip(record['rows'],final['prices'],final['fiscal_values'])), 'Actual smoke fresh replay point differs')
        validate_state_export(candidate['exact_2023_state'],plan['identity'],stage)
    handoff=json.loads(pinned(receipt['handoff_pin']).read_text())
    require(handoff.get('year')==2015 and handoff.get('start_year')==2007 and handoff.get('both_queues_preserved') is True and
        handoff.get('no_rescaling') is True and handoff.get('nonanticipating_prefix') is True and
        handoff.get('source_pins')==plan['identity']['source_pins'],'Actual smoke nonanticipating handoff required')
    with gzip.open(pinned(handoff['checkpoint']),'rb') as stream:packet=pickle.load(stream)
    require(packet['year']==2015 and packet['reference_identity']==plan['identity'] and
        state_hash(packet['households'],queue_array)==handoff['state_sha256'] and
        packet['boundary_price']==handoff['boundary_q'] and packet['boundary_pension']==handoff['boundary_b'] and
        hashlib.sha256(np.asarray(packet['boundary_value']).tobytes()).hexdigest()==handoff['boundary_V_sha256'],
        'Actual smoke handoff checkpoint differs')
    stage2=json.loads(pinned(receipt['stage_seed_pins'][1]).read_text())
    require(stage2['identity']['inherited_state_sha256']==handoff['state_sha256'],'Smoke stage2 seed not inherited from accepted prefix')
    require(len(receipt.get('diagnostic_pins',[]))==2,'Both actual smoke diagnostic packets required')
    for item in receipt['diagnostic_pins']:validate_plot_receipt(item,plan['standard_plot_names'])
    validate_state_export(receipt['exact_2023_state'],plan['identity'],1)
    return receipt


def diagnostic_horizon_comparison(short,long):
    require(long['horizon']>short['horizon']>=5,'Two original diagnostic horizons required')
    for key in ('identity','source_pins','housing','shock_contract'):require(short[key]==long[key],'Horizon identity changed '+key)
    gaps={}
    for key in ('asset_price','renter_price','adult_population','birth_children','housing_demand','pension_period_units'):
        a,b=([r[key] for r in x['rows'][:5]] for x in (short,long));require(len(a)==len(b)==5,'Five original macro rows required')
        gaps[key]=max(abs(x-y)/max(abs(x),abs(y),1e-12) for x,y in zip(a,b))
    a,b=(x['fertility'][:4] for x in (short,long));require(len(a)==len(b)==4,'Four original fertility windows required')
    fertility=[abs(x['period_tfr_topcode_adjusted']-y['period_tfr_topcode_adjusted']) for x,y in zip(a,b)]
    require(all(math.isfinite(v) for v in [*gaps.values(),*fertility]),'Finite horizon gaps required')
    year=short['shock_contract']['start_year']
    return dict(passed=all(v<=.001 for v in gaps.values()) and max(fertility)<=.001,relative_gaps=gaps,
        fertility_absolute_gaps=fertility,early_macro_periods=5,early_macro_years=[year+4*i for i in range(5)],
        fertility_windows=4,macro_tolerance=.001,fertility_tolerance=.001,
        terminal_passes=[short['terminal_pass'],long['terminal_pass']],classification='Strict finite diagnostic horizon comparison; terminal gates remain separate')


def compare_dated_state(short,long,index,queue_values=lambda x:x):
    # Preserve the original physical, value and forecast diagnostics at the actual local index.
    a,b=(x['native_reply'] for x in (short,long));states=[x.dated_states[index]['state'] for x in (a,b)]
    g=[np.asarray(s.g_pre,float) for s in states]
    require(g[0].shape==g[1].shape and all(np.isfinite(x).all() and (x>=0).all() for x in g),'Finite nonnegative common distributions required')
    pops=[float(x.sum()) for x in g];require(min(pops)>0,'Positive physical populations required')
    gaps=dict(population_absolute_gap=abs(pops[0]-pops[1]),population_relative_gap=abs(pops[0]-pops[1])/max(pops),normalized_distribution_l1=float(np.abs(g[0]/pops[0]-g[1]/pops[1]).sum()))
    for key in ('scheduled_entries','scheduled_raw_entries'):
        x,y=(np.asarray(queue_values(getattr(s,key)),float) for s in states)
        require(x.shape==y.shape and all(np.isfinite(v).all() and (v>=0).all() for v in (x,y)),'Both finite nonnegative common queues required')
        gaps[key+'_relative_l1']=float(np.abs(x-y).sum()/max(np.abs(x).sum(),np.abs(y).sum(),1e-12))
    def info(x):
        x=np.asarray(x);return dict(shape=list(x.shape),dtype=str(x.dtype),finite=bool(np.isfinite(x).all()),sha256=hashlib.sha256(x.tobytes()).hexdigest())
    weights=.5*(g[0]/pops[0]+g[1]/pops[1]);values={};integrity=True
    for i,label in ((index,'current_2023_V'),(index+1,'continuation_2027_V')):
        x,y=(np.asarray(r.values[i]) for r in (a,b));record=dict(short=info(x),long=info(y),gating=False)
        if x.shape==y.shape==weights.shape:
            occupied=weights>0;integrity=integrity and bool(np.isfinite(x[occupied]).all() and np.isfinite(y[occupied]).all())
            w=weights[occupied];dx=x[occupied];dy=y[occupied]
            record.update(occupied_finite_weight=float(w.sum()),occupied_weighted_absolute_gap=float(np.sum(w*np.abs(dx-dy))),
                occupied_weighted_relative_gap=float(np.sum(w*np.abs(dx-dy)/np.maximum(1.,np.maximum(np.abs(dx),np.abs(dy))))))
        else: integrity=False;record['weighted_gap_unavailable']='Native value/distribution shape differs'
        values[label]=record
    forecasts={}
    for key in ('prices','pensions','psi_path'):
        x,y=(np.asarray(r.floor_runtime_paths[key],float)[index:] for r in (a,b));n=min(len(x),len(y));require(n>0,'Actual continuation required')
        forecasts[key]=dict(short=info(x),long=info(y),overlap_periods=n,maximum_overlap_relative_gap=float(np.max(np.abs(x[:n]-y[:n])/np.maximum(1e-12,np.maximum(np.abs(x[:n]),np.abs(y[:n]))))),gating=False)
    return dict(passed=all(v<=.001 for k,v in gaps.items() if k!='population_absolute_gap'),available=True,index=index,
        state_value_integrity=integrity,tolerance=.001,state_gaps=gaps,values=values,forecasts=forecasts,
        value_forecast_disclosure='Reported without new acceptance tolerance; exact actual arrays exported')


def table(path,rows):
    require(bool(rows),'Nonempty table required');Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


class Controller:
    def __init__(self,plan,runtime,output):
        self.plan,self.runtime,self.out=plan,runtime,Path(output);self.calls=0;self.count=0;self.last=None
        self.deadline=time.monotonic()+plan['budget']['total_seconds'];self.prefix=[];self.prefix_fertility=[];self.current_state=None
        self.best={};self.fit_receipts=[];self.evidence=[];self.selected_receipts={};self.selected_candidate_pins=[];self.progress_lock=threading.Lock()
    def remaining(self):
        value=self.deadline-time.monotonic();require(value>0,'Two-stage deadline exhausted');return value
    def account(self,reply):
        require(reply.get('accounting_valid') is True,'Native accounting failed')
        n=reply.get('policy_calls');require(type(n) is int and n>=0,'Actual policy-call count required')
        self.calls+=n;require(self.calls<=self.plan['budget']['maximum_policy_calls'],'Shared native-call cap exhausted')
    def progress(self,phase,**extra):
        with self.progress_lock:
            write(self.out/'heartbeat.json',dict(phase=phase,epoch=time.time(),policy_calls=self.calls,candidate=self.count,**extra))
    def evaluate(self,stage,psi):
        self.remaining();self.last=None;self.count+=1;candidate=self.count;folder=self.out/f'stage{stage+1}'/f'candidate_{candidate:04d}'
        deadline=min(self.deadline,time.monotonic()+self.plan['budget']['candidate_seconds']);replies=[];year=STAGES[stage]['start_year']
        for horizon in self.plan['horizons']:
            self.progress('candidate',stage=stage+1,psi=psi,horizon=horizon)
            reply=self.runtime.evaluate_stage(stage=stage,psi=float(psi),horizon=horizon,start_year=year,inherited_state=self.current_state,deadline=deadline,folder=folder/f'horizon_{horizon:03d}')
            self.account(reply)
            require(reply['shock_contract']==dict(start_year=year,psi=float(psi),expectations='permanent_until_next_surprise'),'Permanent stage clock/expectations changed')
            require(reply['source_pins']==self.plan['identity']['source_pins'] and reply['housing']=='static-elastic','Original source/housing identity differs')
            require(reply['identity']==dict(self.plan['identity'],stage_start_year=year,inherited_state_sha256=state_hash(self.current_state,self.runtime.queue_values)), 'Native state identity differs')
            require(all(r['calendar_year']==year+4*i for i,r in enumerate(reply['rows'])),'Actual native stage calendar differs')
            valid=all(reply.get(k) is True for k in ('root_pass','replay_pass','stationary_pass','accounting_valid')) and all(reply[k]<=GATES[g] for k,g in
                (('market_maximum_residual','market_tolerance'),('fiscal_maximum_residual','fiscal_tolerance'),('replay_maximum_gap','final_reproduction_tolerance'),('stationary_renewal_gap','stationary_renewal_tolerance')))
            if not valid:
                fail=dict(certified=False,model=None,payload=dict(candidate=candidate,stage=stage,psi=psi,path=str(folder),error='Original native gates failed'))
                write(folder/'failure.json',fail);return fail
            replies.append(reply)
        comparison=diagnostic_horizon_comparison(*replies);state=compare_dated_state(*replies,4 if stage==0 else 2,self.runtime.queue_values)
        index=4 if stage==0 else 2
        exported=self.runtime.export_state(replies[-1],index=index,year=2023,folder=folder/'state_2023_checkpoint')
        require(exported.get('calendar_year')==2023 and exported.get('period_index')==index and
            exported.get('exact_native_state') is True and exported.get('reconstructed_or_rescaled') is False and
            exported.get('queue_lags')==[16,20] and exported.get('forecast_and_continuation_saved') is True,
            'Exact candidate 2023 state/queues/forecast required')
        export_pin=dict(path=exported['path'],sha256=exported['sha256']);pinned(export_pin)
        checkpoint_receipt=dict(schema=SCHEMA,status='complete',candidate=candidate,stage=stage,
            state_experiment_ready=False,production_ready=False,full_path_certified=False,
            state_physical_horizon_stable=state['passed'],state_value_integrity=state['state_value_integrity'],
            horizon_comparison=comparison,state_comparison=state,export=exported)
        write(folder/'state_2023_checkpoint'/'checkpoint_receipt.json',checkpoint_receipt)
        candidate_evidence=[export_pin,pin(folder/'state_2023_checkpoint'/'checkpoint_receipt.json')]
        for reply in replies:
            require(bool(reply.get('source_evidence')),'Actual root/mapping source evidence required')
            for item in reply['source_evidence']:pinned(item)
            candidate_evidence.extend(reply['source_evidence'])
        write(folder/'horizon_comparison.json',comparison);write(folder/'state_comparison.json',state)
        if not comparison['passed']:
            fail=dict(certified=False,model=None,payload=dict(candidate=candidate,path=str(folder),error='Original strict horizon gates failed'))
            write(folder/'failure.json',fail);return fail
        model=float(replies[-1]['fertility'][1]['period_tfr_topcode_adjusted']);require(math.isfinite(model),'Finite local index1 fertility required')
        target=self.plan['smoke_targets'][stage] if self.plan['smoke'] else TARGETS[1 if stage==0 else 3]
        payload=dict(candidate=candidate,stage=stage,psi=float(psi),path=str(folder),models=[r['period_tfr_topcode_adjusted'] for r in replies[-1]['fertility'][:4]])
        receipt=dict(certified=True,model=model,payload=payload,loss_contribution=(model-target)**2,horizon_comparison=comparison,
            state_comparison=state,terminal_passes=comparison['terminal_passes'],production_ready=False,
            checkpoint_2023=pin(folder/'state_2023_checkpoint'/'checkpoint_receipt.json'),exact_2023_state=exported,evidence_pins=candidate_evidence)
        write(folder/'complete.json',receipt);write(self.out/f'stage{stage+1}'/'latest_completed.json',receipt);write(self.out/'latest_completed.json',receipt)
        if stage not in self.best or receipt['loss_contribution']<self.best[stage]['loss_contribution']:
            self.best[stage]=receipt;write(self.out/f'stage{stage+1}'/'best_so_far.json',receipt)
        self.last=dict(candidate=candidate,stage=stage,psi=float(psi),reply=replies[-1],receipt=receipt)
        return dict(certified=True,model=model,payload=payload)
    def selected(self,stage,result):
        final=result['root']['final'];require(result['converged'] and final is not None and final['mapping_valid'] and self.last is not None and
            self.last['stage']==stage and self.last['psi']==result['parameter']['estimate'] and
            final['payload']['candidate']==self.last['candidate'] and final['payload']['psi']==self.last['psi'], 'Fresh selected replay differs from export candidate')
        return self.last['reply']
    def _handoff(self,accepted,candidate):
        n=2;folder=self.out/'accepted_2007';state=accepted['dated_states'][n]['state'];q=accepted['prices'][n];V=accepted['values'][n]
        self.remaining();replay=self.runtime.replay_prefix(start_year=2007,inherited_state=self.current_state,
            prices=accepted['prices'][:n],pensions=accepted['pensions'][:n],values=accepted['values'][:n],psi=[accepted['psi']]*n,
            boundary_price=q,boundary_value=V,folder=folder,deadline=min(self.deadline,time.monotonic()+self.plan['budget']['mapping_seconds']))
        self.account(replay);require(all(replay['gates'].values()) and max(map(abs,replay['market_residual']))<=self.plan['gates']['market_tolerance'] and
            max(map(abs,replay['fiscal_residual']))<=self.plan['gates']['fiscal_tolerance'],'Original prefix accounting/housing/fiscal gates failed')
        require(len(replay['rows'])==len(replay['fertility'])==n,'Exact two-period prefix required')
        gaps=[]
        for actual,wanted in zip(replay['rows'],accepted['rows'][:n]):
            require(set(actual)==set(wanted),'Prefix macro row fields differ')
            for key,value in wanted.items():
                if isinstance(value,(int,float,np.number)):
                    delta=abs(float(actual[key])-float(value))/max(1.,abs(float(value)));require(math.isfinite(delta),'Finite prefix macro gaps required');gaps.append(delta)
                else:require(actual[key]==value,'Prefix macro metadata differs')
        require(gaps and max(gaps)<=2e-10,'Accepted prefix macro rows changed')
        for i in range(n):
            require(np.asarray(replay['values'][i]).shape==np.asarray(accepted['values'][i]).shape and
                np.allclose(replay['values'][i],accepted['values'][i],rtol=0,atol=2e-10),'Accepted prefix values changed')
            require(abs(replay['fertility'][i]['period_tfr_topcode_adjusted']-accepted['fertility'][i]['period_tfr_topcode_adjusted'])<=self.plan['fit']['reproduction_tolerance'],'Accepted prefix fertility changed')
        for key in ('g_pre','scheduled_entries','scheduled_raw_entries'):
            a,b=(getattr(x,key) for x in (replay['terminal_state'],state))
            if key!='g_pre':a,b=self.runtime.queue_values(a),self.runtime.queue_values(b)
            a,b=np.asarray(a),np.asarray(b);require(a.shape==b.shape and np.allclose(a,b,rtol=0,atol=2e-10),'Exact inherited handoff failed: '+key)
        folder.mkdir(parents=True,exist_ok=True);checkpoint=folder/'inherited_2015.pkl.gz'
        packet=dict(year=2015,households=copy.deepcopy(replay['terminal_state']),reference_identity=self.plan['identity'],
            candidate=candidate,boundary_price=float(q),boundary_pension=float(accepted['pensions'][n]),boundary_value=np.asarray(V).copy())
        with gzip.open(checkpoint,'wb') as stream:pickle.dump(packet,stream,protocol=pickle.HIGHEST_PROTOCOL)
        with gzip.open(checkpoint,'rb') as stream:restored=pickle.load(stream)
        for key in ('g_pre','scheduled_entries','scheduled_raw_entries'):
            a,b=(getattr(x,key) for x in (restored['households'],replay['terminal_state']))
            if key!='g_pre':a,b=self.runtime.queue_values(a),self.runtime.queue_values(b)
            require(np.array_equal(a,b),'Checkpoint array roundtrip changed '+key)
        require(np.array_equal(restored['boundary_value'],V),'Checkpoint forecast boundary value changed')
        receipt=dict(checkpoint=pin(checkpoint),year=2015,start_year=2007,candidate=candidate,state_sha256=state_hash(restored['households'],self.runtime.queue_values),
            boundary_q=float(q),boundary_b=float(accepted['pensions'][n]),boundary_V_sha256=hashlib.sha256(np.asarray(V).tobytes()).hexdigest(),
            source_pins=self.plan['identity']['source_pins'],both_queues_preserved=True,no_rescaling=True,nonanticipating_prefix=True)
        write(folder/'handoff.json',receipt);self.evidence.extend([pin(checkpoint),pin(folder/'handoff.json')])
        self.prefix=copy.deepcopy(replay['rows']);self.prefix_fertility=copy.deepcopy(replay['fertility']);self.current_state=restored['households']
        self.runtime.set_stage2_initialization(accepted,candidate)
        prefix_reply=dict(accepted,native_reply=replay['native_reply']);self.render(prefix_reply,folder/'diagnostics')
        return receipt
    def render(self,reply,folder):
        plots=self.runtime.render_standard(reply,folder,deadline=min(self.deadline,time.monotonic()+self.plan['budget']['render_seconds']))
        names=set(self.plan['standard_plot_names']);require(isinstance(plots,dict) and plots.get('sampled_dates'),'Actual native dated diagnostic packet required')
        for dated in plots['sampled_dates']:
            require(set(dated['plots'])==names,'Exact 17 standard plot names required')
            directory=Path(folder)/f"date_{dated['period']:03d}"/'standard_diagnostics'
            for name,value in dated['plots'].items():require((directory/name).is_file() and sha(directory/name)==value,'Diagnostic plot hash differs')
        write(Path(folder)/'diagnostics.json',plots);self.evidence.append(pin(Path(folder)/'diagnostics.json'));return plots
    def run(self):
        preflight(self.plan);self.runtime.bind_budget(self.deadline,self.progress)
        with self.runtime.budget_context():return self._smoke_run() if self.plan['smoke'] else self._run()
    def _smoke_run(self):
        # Execution integration only: one fixed candidate per stage. Native
        # path roots perform their own fresh replay; no scalar fitter is called.
        write(self.out/'executed_manifest.json',self.plan)
        prepared=self.runtime.prepare_reference(min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']),self.out/'reference');self.account(prepared)
        require(all(k in prepared for k in ('reference_checkpoint','reconstruction_receipt')),'Fresh actual reference evidence required')
        self.evidence.extend([prepared[k] for k in ('reference_checkpoint','reconstruction_receipt')])
        self.current_state=self.runtime.initial_state();seed_pins=[];candidate_pins=[];diagnostic_pins=[];terminal_passes={}
        for stage,settings in enumerate(STAGES):
            self.remaining();folder=self.out/f'stage{stage+1}'
            seed=self.runtime.measure_seed(stage=stage,start_year=settings['start_year'],inherited_state=self.current_state,
                folder=folder/'seed',deadline=min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']))
            self.account(seed);require(seed['mapping_count']==5 and seed['horizon']==12 and
                (stage==0 or seed['nonstationary_inherited_state'] is True),'Original fresh stage-bound seed required')
            require(seed['identity']['inherited_state_sha256']==state_hash(self.current_state,self.runtime.queue_values),'Fresh seed state identity differs')
            require(len(seed.get('source_evidence',[]))==5,'Five actual native seed map pins required')
            self.evidence.extend(seed['source_evidence'])
            write(folder/'seed_receipt.json',seed);seed_pins.append(pin(folder/'seed_receipt.json'))
            self.count+=1;candidate=self.count;path=folder/f'candidate_{candidate:04d}'/'horizon_006'
            self.progress('fixed_baseline_smoke_candidate',stage=stage+1,psi=BASELINE_PSI,horizon=6)
            reply=self.runtime.evaluate_stage(stage=stage,psi=BASELINE_PSI,horizon=6,start_year=settings['start_year'],
                inherited_state=self.current_state,deadline=min(self.deadline,time.monotonic()+self.plan['budget']['candidate_seconds']),folder=path)
            self.account(reply)
            require(reply['shock_contract']==dict(start_year=settings['start_year'],psi=BASELINE_PSI,expectations='permanent_until_next_surprise') and
                reply['identity']==dict(self.plan['identity'],stage_start_year=settings['start_year'],inherited_state_sha256=state_hash(self.current_state,self.runtime.queue_values)) and
                reply['source_pins']==self.plan['identity']['source_pins'] and reply['housing']=='static-elastic','Actual fixed smoke source/state/clock differs')
            require(reply['horizon']==6 and len(reply['rows'])==len(reply['fertility'])==6 and
                all(r['calendar_year']==settings['start_year']+4*i for i,r in enumerate(reply['rows'])),'Exact reduced smoke horizon/calendar required')
            require(all(reply.get(k) is True for k in ('root_pass','replay_pass','stationary_pass','accounting_valid')) and
                all(math.isfinite(reply[k]) and reply[k]<=self.plan['gates'][g] for k,g in
                    (('market_maximum_residual','market_tolerance'),('fiscal_maximum_residual','fiscal_tolerance'),
                     ('replay_maximum_gap','final_reproduction_tolerance'),('stationary_renewal_gap','stationary_renewal_tolerance'))),
                'Smoke native root/replay/accounting gates failed')
            require(reply.get('source_evidence'),'Actual smoke root/mapping evidence required')
            for item in reply['source_evidence']:pinned(item)
            exported=self.runtime.export_state(reply,index=4 if stage==0 else 2,year=2023,folder=folder/'state_2023_checkpoint')
            exported['state_sha256']=state_hash(reply['dated_states'][4 if stage==0 else 2]['state'],self.runtime.queue_values)
            validate_state_export(exported,self.plan['identity'],stage)
            require(type(reply.get('terminal_pass')) is bool and isinstance(reply.get('terminal'),dict),'Actual terminal diagnostics required')
            write(folder/'terminal_diagnostics.json',reply['terminal']);terminal_pin=pin(folder/'terminal_diagnostics.json')
            self.evidence.append(terminal_pin);terminal_passes[str(stage+1)]=reply['terminal_pass']
            receipt=dict(status='execution_passed',certified=False,scientific_validation=False,production_ready=False,
                stage=stage,candidate=candidate,psi=BASELINE_PSI,horizon=6,fresh_native_replay=True,accounting_valid=True,
                inherited_state_sha256=state_hash(self.current_state,self.runtime.queue_values),root_pin=pin(path/'root.json'),
                source_evidence=reply['source_evidence'],final_mapping_pin=reply['final_mapping_pin'],path_evaluations=reply['path_evaluations'],
                exact_2023_state=exported,empirical_horizon_gates_tested=False,terminal_pass=reply['terminal_pass'],
                terminal_diagnostics_pin=terminal_pin,terminal_diagnostics_gating=False)
            write(folder/'candidate_execution.json',receipt);candidate_pins.append(pin(folder/'candidate_execution.json'))
            self.evidence.extend(reply['source_evidence']+[dict(path=exported['path'],sha256=exported['sha256'])])
            write(folder/'latest_completed.json',receipt);write(self.out/'latest_completed.json',receipt)
            if stage==0:
                self._handoff(reply,candidate);diagnostic_pins.append(pin(self.out/'accepted_2007'/'diagnostics'/'diagnostics.json'))
                self.runtime.release_stage1();reply=None
        self.render(reply,self.out/'diagnostics');diagnostic_pins.append(pin(self.out/'diagnostics'/'diagnostics.json'))
        actual=getattr(getattr(self.runtime,'rt',None),'total_native_calls',self.calls)
        require(actual==self.calls,'Controller/native actual smoke call accounting differs')
        receipt=dict(status='execution_passed',schema=SCHEMA,smoke=True,scientific_validation=False,production_ready=False,empirical_fitted=False,
            fingerprints=fingerprints(self.plan),smoke_manifest_pin=pin(self.out/'executed_manifest.json'),
            semantic_gates=dict(fixed_baseline_candidates=True,scalar_optimizer=False,derivative_probes=False,nonanticipating_prefix=True,
                both_queues=True,fresh_replays=True,exact_2023=True,standard_diagnostics=True,
                empirical_horizon_gates_tested=False,empirical_forecast_slice_tested=False),
            reference_evidence=prepared,stage_seed_pins=seed_pins,selected_candidate_pins=candidate_pins,handoff_pin=pin(self.out/'accepted_2007'/'handoff.json'),
            diagnostic_pins=diagnostic_pins,exact_2023_state=exported,evidence_pins=self.evidence,
            actual_policy_calls=self.calls,native_actual_policy_calls=actual,terminal_passes_by_stage=terminal_passes,
            terminal_diagnostics_gating=False,
            classification='Execution-only native integration; no scalar fit or empirical convergence certification',
            untested=['empirical 24/32 horizon stability','full scalar fits and derivatives','actual stage1 forecast slice [2:14]'],
            initializer='Smoke-only accepted short forecast plus stationary endpoint numerical padding',disclosure=self.plan['disclosure'])
        write(self.out/'complete.json',receipt);return receipt

    def _run(self):
        prepared=self.runtime.prepare_reference(min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']),self.out/'reference');self.account(prepared)
        self.evidence.extend([prepared[k] for k in ('reference_checkpoint','reconstruction_receipt') if k in prepared])
        self.current_state=self.runtime.initial_state();results=[];parameters=[];accepted=None
        for stage,settings in enumerate(STAGES):
            deadline=min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds'])
            seed=self.runtime.measure_seed(stage=stage,start_year=settings['start_year'],inherited_state=self.current_state,folder=self.out/f'stage{stage+1}'/'seed',deadline=deadline)
            self.account(seed);require(seed['mapping_count']==5 and seed['horizon']==12 and (stage==0 or seed['nonstationary_inherited_state'] is True),'Original fresh stage-bound seed required')
            write(self.out/f'stage{stage+1}'/'seed_receipt.json',seed)
            self.evidence.append(pin(self.out/f'stage{stage+1}'/'seed_receipt.json'))
            target=self.plan['smoke_targets'][stage] if self.plan['smoke'] else TARGETS[1 if stage==0 else 3]
            result=self.runtime.fitter.fit_one(evaluate=lambda x:self.evaluate(stage,x),target=target,
                initial_level=self.plan['stage_starts'][stage]['initial'],bounds=self.plan['stage_starts'][stage]['bounds'],
                controls=dict(self.plan['fit'],total_seconds=self.remaining()),callback=lambda row:self.progress('scalar_fit',stage=stage+1,root_row=row))
            fit_path=self.out/f'stage{stage+1}'/'fit.json';write(fit_path,result);self.evidence.append(pin(fit_path));accepted=self.selected(stage,result)
            self.selected_receipts[stage]=copy.deepcopy(self.last['receipt'])
            self.selected_candidate_pins.append(pin(Path(self.last['receipt']['payload']['path'])/'complete.json'))
            self.evidence.extend(self.last['receipt']['evidence_pins'])
            parameter=dict(parameter='psi_'+str(settings['start_year']),**result['parameter']);parameters.append(parameter);results.append(result)
            if stage==0:
                self._handoff(accepted,self.last['candidate'])
                accepted=None
                self.last['reply']=None
                self.runtime.release_stage1()
            write(self.out/f'stage{stage+1}'/'complete.json',dict(status='matched',parameter=parameter,candidate=self.last['candidate'],selected_replay_verified=True))
        # No further state advance or reevaluation: fit_one's final fresh replay is the export.
        models=self.prefix_fertility+accepted['fertility'][:2];require(len(models)==4,'Four calendar-lineage fertility rows required')
        rows=[]
        for i,(target,model) in enumerate(zip(self.plan['rows'],models)):
            value=float(model['period_tfr_topcode_adjusted']);gap=value-target['target'];weight=self.plan['weights'][i]
            rows.append(dict(target,model=value,gap=gap,weight=weight,loss_contribution=weight*gap*gap,
                source_stage=1 if i<2 else 2,local_index=i if i<2 else i-2,role='fitted' if weight else 'validation'))
        table(self.out/'fertility_fit.csv',rows);table(self.out/'parameters.csv',parameters)
        projection=[dict(row,source_stage=1 if i<2 else 2,local_index=i if i<2 else i-2,global_period=i) for i,row in enumerate(self.prefix+accepted['rows'])]
        table(self.out/'final_projection.csv',projection)
        exported=self.last['receipt']['exact_2023_state']
        pinned(self.last['receipt']['checkpoint_2023'])
        require(exported.get('exact_native_state') is True and exported.get('reconstructed_or_rescaled') is False and
            exported.get('calendar_year')==2023 and exported.get('period_index')==2 and exported.get('queue_lags')==[16,20] and
            exported.get('forecast_and_continuation_saved') is True,'Exact native 2023 state/queues/continuation required')
        self.evidence.append(dict(path=exported['path'],sha256=exported['sha256']));pinned(self.evidence[-1])
        plots=self.render(accepted,self.out/'diagnostics')
        actual=getattr(getattr(self.runtime,'rt',None),'total_native_calls',self.calls)
        if hasattr(self.runtime,'adapters'):require(actual==self.calls,'Controller/native actual-call accounting differs')
        receipt=dict(status='matched',schema=SCHEMA,smoke=self.plan['smoke'],scientific_validation=False,production_ready=False,
            fingerprints=fingerprints(self.plan),semantic_gates=dict(two_fits=True,nonanticipating_prefix=True,both_queues=True,
                selected_replays=True,original_horizon_gates=True,exact_2023=True,standard_diagnostics=True),
            selected_candidate_pins=self.selected_candidate_pins,parameters=parameters,fit_table=rows,total_loss=sum(r['loss_contribution'] for r in rows),actual_policy_calls=self.calls,
            native_actual_policy_calls=getattr(getattr(self.runtime,'rt',None),'total_native_calls',self.calls),
                        stage_fit_pins=[pin(self.out/f'stage{i+1}'/'fit.json') for i in range(2)],exact_2023_state=exported,evidence_pins=self.evidence,
            final_horizon_comparison=self.last['receipt']['horizon_comparison'],final_state_comparison=self.last['receipt']['state_comparison'],
            terminal_passes_by_stage={str(i+1):self.selected_receipts[i]['terminal_passes'] for i in range(2)},
            continuation_validation='Outstanding; physical comparison and values/forecast diagnostics do not certify production convergence',
            classification='Synthetic zero-surprise two-stage smoke; not an empirical fit' if self.plan['smoke'] else 'Two-surprise 24/32 finite diagnostic fit; terminal failures explicitly retained',
            disclosure=self.plan['disclosure'])
        write(self.out/'complete.json',receipt);return receipt


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',required=True);parser.add_argument('--output',required=True)
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--preflight',action='store_true');mode.add_argument('--smoke',action='store_true');mode.add_argument('--run',action='store_true')
    parser.add_argument('--smoke-receipt-pin',help='JSON object with exact path and sha256');args=parser.parse_args(argv)
    plan=json.loads(Path(args.manifest).read_text());out=Path(args.output);controller=None;runtime=None
    try:
        receipt=preflight(plan)
        if args.preflight:write(out/'preflight.json',receipt);return receipt
        require(plan['smoke'] is args.smoke,'Manifest smoke/run mode differs from CLI')
        if args.run:
            require(args.smoke_receipt_pin is not None,'Run requires passed two-stage smoke pin')
        module=importlib.import_module('two_shock_runtime');require(Path(module.__file__).resolve()==pinned(plan['source_pins']['two_shock_runtime']).resolve(),'Concrete runtime import differs')
        runtime=module.build_runtime(plan=plan,output=out/'runtime',smoke=args.smoke)
        if args.run:
            require(runtime.rt.total_native_calls==0,'Smoke release check must precede every native call')
            validate_smoke_pin(plan,json.loads(args.smoke_receipt_pin))
        controller=Controller(plan,runtime,out);result=controller.run();write(out/('smoke_receipt.json' if args.smoke else 'run_receipt.json'),result);return result
    except BaseException as exc:
        write(out/'failure.json',dict(status='failed',schema=SCHEMA,exception=type(exc).__name__,reason=str(exc),
            policy_calls=None if controller is None else controller.calls,native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',None),
            completed_cases_retained=True,production_ready=False));raise


if __name__=='__main__':main()
