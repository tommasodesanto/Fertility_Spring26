#!/usr/bin/env python3
"""Isolated, finite-H16 rough two-surprise estimator; never a production certificate."""
from __future__ import annotations
import argparse
import csv
import contextlib
import copy
import importlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
import numpy as np
import two_shock as original

SCHEMA='estate_a_two_unanticipated_psi_rough_h16_v1'
ROOT=original.ROOT
HERE=Path(__file__).resolve().parent
BOUNDS=[0.001789207206604163,0.35784144132083257]
BASELINE_EXPECTED=0.176948250201189
ROUGH_GATES=dict(original.GATES,market_tolerance=.005,fiscal_tolerance=.001)
ROUGH_BUDGET=dict(total_seconds=14280,maximum_policy_calls=19897,candidate_seconds=7200,path_seconds=6000,
                  seed_seconds=1800,mapping_seconds=1800,render_seconds=120,endpoint_seconds=1800)
SMOKE_BUDGET=dict(total_seconds=900,maximum_policy_calls=100,candidate_seconds=360,path_seconds=240,
                  seed_seconds=360,mapping_seconds=180,render_seconds=120,endpoint_seconds=240)
STAGES=copy.deepcopy(original.STAGES)
TARGETS=original.TARGETS


def check_new_wealth_target(case):
    path=Path(case)/'target_fit_new_contract.csv'
    with path.open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    wealth=[row for row in rows if row.get('moment')=='wealth_earnings']
    original.require(len(rows)==14 and len(wealth)==1 and
                     abs(float(wealth[0]['target'])-4.45838713455674)<1e-10,
                     'New Estate-A 4.458387 wealth target required')
    contract=json.loads((Path(case)/'input_contract.json').read_text())
    original.require(contract.get('target_fingerprint')=='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70' and
                     contract.get('weight_fingerprint')=='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4',
                     'Bridge input contract target/weight fingerprints differ')
    return original.pin(path)


def prepare_manifest(base_plan_pin,*,smoke=False,deadline_epoch=None):
    base=json.loads(original.pinned(base_plan_pin).read_text())
    psi=float(base['initial_psi'])
    original.require(math.isclose(psi,BASELINE_EXPECTED,rel_tol=0,abs_tol=1e-12),'Selected new-baseline psi differs')
    controls=dict(gates=copy.deepcopy(ROUGH_GATES),seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),
                  fit=copy.deepcopy(base['fit']),path=copy.deepcopy(base['path']),endpoint=copy.deepcopy(base['endpoint']),
                  budget=dict(base['budget'],**ROUGH_BUDGET),horizons=[16],smoke_seed_endpoint_padding=False)
    controls['fit'].update(max_evaluations=12,fertility_tolerance=.005)
    controls['path']['max_evaluations']=12
    controls['endpoint']['max_evaluations']=48
    if smoke:
        controls['seed']=dict(horizon=4,perturbed_date=1,log_step=1e-5)
        controls['horizons']=[6]
        controls['smoke_seed_endpoint_padding']=True
        controls['gates']=dict(original.GATES,market_tolerance=.05,fiscal_tolerance=.005)
        controls['path']['max_evaluations']=3
        controls['endpoint']['max_evaluations']=8
        controls['budget']=dict(controls['budget'],**SMOKE_BUDGET)
    plan=dict(schema=SCHEMA,kind='two_unanticipated_permanent_rough',mode='rough_diagnostic',smoke=bool(smoke),
              baseline_psi=psi,absolute_psi_bounds=BOUNDS,stages=copy.deepcopy(STAGES),
              rows=copy.deepcopy(base['target_contract']['rows']),weights=[0,1,0,1],
              target_contract=copy.deepcopy(base['target_contract']),base_plan=base_plan_pin,
              new_wealth_target_pin=copy.deepcopy(base['new_wealth_target_pin']),
              source_pins=dict(rough_driver=original.pin(__file__),rough_runtime=original.pin(HERE/'rough_two_shock_runtime.py')),
              identity=copy.deepcopy(base['identity']),standard_plot_names=base['standard_plot_names'],
              stage_starts=[dict(initial=x,bounds=BOUNDS) for x in ((psi,psi) if smoke else (.14736308634876963,.12))],
              deadline_epoch=deadline_epoch,execution_enabled=True,reviewed=True,
              disclosure='Experimental new-baseline two-surprise rough H16 diagnostic; original H24/H32 horizon gates untested; estate closure provisional.')
    plan.update(controls)
    if base.get('legacy_source_overlay') is not None:
        plan['legacy_source_overlay']=copy.deepcopy(base['legacy_source_overlay'])
    return plan


def prepare_base_plan(bridge_case,output,source_template_manifest):
    """Pin a new selected Estate-A case using the retained transition constructor."""
    import transition
    original.require(transition.ROOT==ROOT,'Transition builder source root differs')
    template=json.loads(Path(source_template_manifest).read_text())
    original.require('identity' in template and bool(template['identity']['source_pins']),
                     'Authenticated source-template manifest required')
    source_paths=[ROOT/name for name in template['identity']['source_pins']]
    overlay=template.get('legacy_source_overlay')
    if overlay is not None:
        original.require(set(overlay)=={'helper','manifest'},'Exact legacy overlay pins required')
        helper=original.pinned(overlay['helper']);metadata=original.pinned(overlay['manifest'])
        spec=importlib.util.spec_from_file_location('rough_two_shock_authenticated_overlay',helper)
        module=importlib.util.module_from_spec(spec)
        sys.modules[spec.name]=module
        spec.loader.exec_module(module)
        module.install_overlay(metadata)
    wealth_pin=check_new_wealth_target(bridge_case)
    output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
    handoff_path=output.with_name(output.stem+'_handoff.json')
    handoff=transition.build_handoff(bridge_case,source_paths)
    original.write(handoff_path,handoff)
    handoff_pin=original.pin(handoff_path)
    runtime=transition.runtime_module().CurrentEstateARuntime.from_handoff(
        handoff_pin,output.parent/'base_constructor')
    original.require(runtime.total_native_calls==0,'Base constructor made a native call')
    plan=transition.build_plan(handoff_pin,runtime,mode='diagnostic',horizons=[24,32],
        budget=dict(ROUGH_BUDGET),fit_max_evaluations=12,path_max_evaluations=12,
        endpoint_max_evaluations=48,fiscal_relaxation_authorized=True,execution_enabled=True)
    original.require(math.isclose(plan['initial_psi'],BASELINE_EXPECTED,rel_tol=0,abs_tol=1e-12),
                     'Bridge case has wrong selected new-baseline preference')
    original.require([x['target'] for x in plan['target_contract']['rows']]==TARGETS,
                     'Bridge target-builder four-window contract differs')
    plan['new_wealth_target_pin']=wealth_pin
    if overlay is not None:plan['legacy_source_overlay']=copy.deepcopy(overlay)
    original.write(output,plan)
    return dict(base_plan=original.pin(output),handoff=handoff_pin,native_calls=0,
                source_identity=plan['identity'],baseline_psi=plan['initial_psi'])


def preflight(plan):
    original.require(plan.get('schema')==SCHEMA and plan.get('kind')=='two_unanticipated_permanent_rough' and
                     plan.get('mode')=='rough_diagnostic' and type(plan.get('smoke')) is bool,'Isolated rough schema required')
    original.require(plan.get('execution_enabled') is True and plan.get('reviewed') is True,'Reviewed rough plan required')
    original.require(plan['stages']==STAGES and plan['weights']==[0,1,0,1] and
                     [x['target'] for x in plan['rows']]==TARGETS,'Exact four-row target contract required')
    original.require(plan['absolute_psi_bounds']==BOUNDS and len(plan['stage_starts'])==2 and
                     all(x['bounds']==BOUNDS and BOUNDS[0]<x['initial']<BOUNDS[1] for x in plan['stage_starts']),
                     'Original absolute psi bounds required')
    base=json.loads(original.pinned(plan['base_plan']).read_text())
    original.require(base['initial_psi']==plan['baseline_psi'] and
                     math.isclose(plan['baseline_psi'],BASELINE_EXPECTED,rel_tol=0,abs_tol=1e-12) and
                     base['identity']==plan['identity'] and base['target_contract']==plan['target_contract'],
                     'New Estate-A baseline identity/psi/target differs')
    expected=prepare_manifest(plan['base_plan'],smoke=plan['smoke'],deadline_epoch=plan.get('deadline_epoch'))
    original.require(plan==expected,'Rough manifest differs from exact prepared controls/source pins')
    if plan.get('legacy_source_overlay') is not None:
        original.require(set(plan['legacy_source_overlay'])=={'helper','manifest'},'Legacy overlay fields differ')
        for item in plan['legacy_source_overlay'].values():original.pinned(item)
    for item in base['source_files'].values():original.pinned(item)
    frozen_root=original.pinned(base['source_files']['runtime']).parents[4]
    original.require(ROOT==frozen_root,'Driver must run inside authenticated staged source root')
    for name,hash_value in base['identity']['source_pins'].items():
        path=(frozen_root/name).resolve()
        original.require(path.is_relative_to(frozen_root) and original.sha(path)==hash_value,'Frozen native source differs: '+name)
    handoff=json.loads(original.pinned(base['handoff']).read_text())
    original.require(handoff['source_pins']==base['identity']['source_pins'],'Handoff source pins differ')
    case=(frozen_root/handoff['saved_case']).resolve()
    original.require(case.is_relative_to(frozen_root),'Saved case escapes frozen root')
    original.require(plan['new_wealth_target_pin']==base['new_wealth_target_pin']==check_new_wealth_target(case),
                     'Pinned new Estate-A wealth target differs')
    for name,value in {**handoff['saved_files'],**handoff.get('standard_plot_pins',{})}.items():
        path=(case/name).resolve()
        original.require(path.is_relative_to(case) and original.sha(path)==value,'Saved input differs: '+name)
    original.require(len(plan['standard_plot_names'])==len(set(plan['standard_plot_names']))==17,'Standard plot contract differs')
    if plan['smoke']:
        original.require(plan['horizons']==[6] and plan['seed']==dict(horizon=4,perturbed_date=1,log_step=1e-5) and
                         plan['smoke_seed_endpoint_padding'] is True,'Tiny smoke controls required')
    else:
        original.require(plan['horizons']==[16] and plan['seed']==dict(horizon=12,perturbed_date=5,log_step=1e-5) and
                         plan['smoke_seed_endpoint_padding'] is False,'Actual H16 forecast controls required')
    original.require(plan['gates']['final_reproduction_tolerance']==1e-10 and
                     plan['gates']['stationary_renewal_tolerance']==1e-6 and
                     plan['fit']['max_evaluations']==12 and plan['fit']['fertility_tolerance']==.005 and
                     plan['budget']['maximum_policy_calls']<=20000,'Replay, scalar or shared-call guard changed')
    if plan.get('deadline_epoch') is not None:
        original.require(type(plan['deadline_epoch']) in (int,float) and math.isfinite(plan['deadline_epoch']),
                         'Finite shared external deadline required')
    return dict(status='PASS',schema=SCHEMA,native_calls=0,source_identity=plan['identity'],
                rough_diagnostic=True,production_ready=False,original_horizon_gates_tested=False)


class Controller(original.Controller):
    def __init__(self,plan,runtime,output):
        super().__init__(plan,runtime,output)
        if plan.get('deadline_epoch') is not None:
            self.deadline=min(self.deadline,time.monotonic()+max(0.,float(plan['deadline_epoch'])-time.time()))

    def run(self):
        preflight(self.plan)
        self.runtime.bind_budget(self.deadline,self.progress)
        with self.runtime.budget_context():
            return self._smoke() if self.plan['smoke'] else self._run_rough()

    def render(self,reply,folder):
        if self.plan['smoke']:
            receipt=dict(status='not_tested',reason='Tiny execution smoke excludes standard diagnostic rendering',
                         standard_diagnostics_tested=False,production_ready=False)
            original.write(Path(folder)/'diagnostics.json',receipt)
            return receipt
        return super().render(reply,folder)

    def _prune_superseded(self,stage,current_candidate):
        keep={current_candidate}
        if stage in self.best:keep.add(self.best[stage]['candidate'])
        for item in self.selected_candidate_pins:
            path=Path(item['path'])
            if path.parent.parent.name==f'stage{stage+1}':
                keep.add(int(path.parent.name.split('_')[-1]))
        removed=[]
        for folder in sorted((self.out/f'stage{stage+1}').glob('candidate_*')):
            if not folder.is_dir() or not (folder/'summary.json').is_file():continue
            number=int(folder.name.split('_')[-1])
            if number in keep:continue
            count=self.runtime.prune_candidate_packets(folder)
            if count:removed.append(dict(candidate=number,diagnostic_packets_removed=count))
        if removed:
            original.write(self.out/f'stage{stage+1}'/'storage_pruning.json',
                           dict(removed=removed,kept_candidate_numbers=sorted(keep),
                                numeric_records_and_summaries_retained=True))

    def evaluate(self,stage,psi):
        self.remaining();self.last=None;self.count+=1;candidate=self.count
        folder=self.out/f'stage{stage+1}'/f'candidate_{candidate:04d}'
        year=STAGES[stage]['start_year'];horizon=self.plan['horizons'][0]
        deadline=min(self.deadline,time.monotonic()+self.plan['budget']['candidate_seconds'])
        self.progress('rough_candidate',stage=stage+1,psi=psi,horizon=horizon)
        reply=self.runtime.evaluate_stage(stage=stage,psi=float(psi),horizon=horizon,start_year=year,
            inherited_state=self.current_state,deadline=deadline,folder=folder)
        self.account(reply)
        original.require(reply['shock_contract']==dict(start_year=year,psi=float(psi),expectations='permanent_until_next_surprise') and
                         reply['source_pins']==self.plan['identity']['source_pins'] and reply['housing']=='static-elastic' and
                         reply['identity']==dict(self.plan['identity'],stage_start_year=year,
                                                inherited_state_sha256=original.state_hash(self.current_state,self.runtime.queue_values)),
                         'Rough native source/state/shock identity differs')
        original.require(len(reply['rows'])==len(reply['fertility'])==horizon and
                         all(r['calendar_year']==year+4*i for i,r in enumerate(reply['rows'])),
                         'Dated rough path clock differs')
        valid=all(reply.get(k) is True for k in ('root_pass','replay_pass','stationary_pass','accounting_valid')) and all(
            math.isfinite(reply[k]) and reply[k]<=self.plan['gates'][g] for k,g in
            (('market_maximum_residual','market_tolerance'),('fiscal_maximum_residual','fiscal_tolerance'),
             ('replay_maximum_gap','final_reproduction_tolerance'),('stationary_renewal_gap','stationary_renewal_tolerance')))
        payload=dict(candidate=candidate,stage=stage,psi=float(psi),path=str(folder),horizon=horizon,
                     original_horizon_gates_tested=False)
        compact=dict(candidate=candidate,stage=stage,psi=float(psi),rough_objective_eligible=bool(valid),
                     housing_residual=reply['market_maximum_residual'],fiscal_residual=reply['fiscal_maximum_residual'],
                     replay_gap=reply['replay_maximum_gap'],stationary_renewal_gap=reply['stationary_renewal_gap'],
                     actual_policy_calls=self.calls,terminal_diagnostic_pass=reply.get('terminal_pass'),
                     original_horizon_gates_tested=False,production_ready=False)
        if valid:
            target=TARGETS[1 if stage==0 else 3]
            model=float(reply['fertility'][1]['period_tfr_topcode_adjusted'])
            original.require(math.isfinite(model),'Finite local fertility required')
            compact.update(model=model,target=target,gap=model-target,loss_contribution=(model-target)**2)
            receipt=dict(certified=True,model=model,payload=payload,loss_contribution=compact['loss_contribution'],
                         terminal_passes=[reply.get('terminal_pass')],production_ready=False,
                         original_horizon_gates_tested=False,horizon_comparison=None,state_comparison=None,
                         evidence_pins=reply['source_evidence'])
            self.last=dict(candidate=candidate,stage=stage,psi=float(psi),reply=reply,receipt=receipt)
            if stage not in self.best or compact['loss_contribution']<self.best[stage]['loss_contribution']:
                self.best[stage]=dict(compact,source_evidence=reply['source_evidence'])
                original.write(self.out/f'stage{stage+1}'/'best_so_far.json',self.best[stage])
        original.write(folder/'summary.json',compact)
        original.write(self.out/f'stage{stage+1}'/'latest_completed.json',compact)
        original.write(self.out/'latest_completed.json',compact)
        if hasattr(self.runtime,'prune_candidate_packets'):
            self._prune_superseded(stage,candidate)
        return dict(certified=bool(valid),model=compact.get('model'),payload=payload)

    def _prepare(self):
        prepared=self.runtime.prepare_reference(min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']),self.out/'reference')
        self.account(prepared)
        self.evidence.extend([prepared[k] for k in ('reference_checkpoint','reconstruction_receipt') if k in prepared])
        self.current_state=self.runtime.initial_state()
        return prepared

    def _seed(self,stage):
        settings=STAGES[stage];folder=self.out/f'stage{stage+1}'/'seed'
        seed=self.runtime.measure_seed(stage=stage,start_year=settings['start_year'],inherited_state=self.current_state,
             folder=folder,deadline=min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']))
        self.account(seed)
        original.require(seed['mapping_count']==5 and seed['horizon']==self.plan['seed']['horizon'] and
                         seed.get('perturbed_date')==self.plan['seed']['perturbed_date'] and
                         seed.get('perturbation_log_step')==self.plan['seed']['log_step'] and
                         len(seed.get('source_evidence',[]))==5 and
                         (stage==0 or seed['nonstationary_inherited_state'] is True),'Actual five-map state-bound seed required')
        for evidence in seed['source_evidence']:original.pinned(evidence)
        original.require(seed['identity']==dict(self.plan['identity'],stage_start_year=settings['start_year'],
                         inherited_state_sha256=original.state_hash(self.current_state,self.runtime.queue_values)),
                         'Seed source/state identity differs')
        original.write(self.out/f'stage{stage+1}'/'seed_receipt.json',seed)
        self.evidence.append(original.pin(self.out/f'stage{stage+1}'/'seed_receipt.json'))
        return seed

    def _selected(self,stage,result):
        accepted=self.selected(stage,result)
        candidate=self.last['candidate'];folder=self.out/f'stage{stage+1}'/f'candidate_{candidate:04d}'
        exported=self.runtime.export_state(accepted,index=4 if stage==0 else 2,year=2023,
                                           folder=folder/'state_2023_checkpoint')
        original.require(exported.get('exact_native_state') is True and exported.get('reconstructed_or_rescaled') is False and
                         exported.get('calendar_year')==2023 and exported.get('period_index')==(4 if stage==0 else 2) and
                         exported.get('queue_lags')==[16,20] and exported.get('forecast_and_continuation_saved') is True,
                         'Exact selected 2023 state/queues/forecast required')
        original.pinned({k:exported[k] for k in ('path','sha256')})
        self.evidence.append({k:exported[k] for k in ('path','sha256')})
        self.last['receipt']['exact_2023_state']=exported
        selected=dict(self.last['receipt'],source_evidence=accepted['source_evidence'],
                      original_horizon_gates_tested=False,rough_h16_only=True)
        original.write(folder/'selected.json',selected)
        self.selected_candidate_pins.append(original.pin(folder/'selected.json'))
        return accepted,exported

    def _run_rough(self):
        self._prepare();results=[];parameters=[];selected_receipts={};accepted=None;exported=None
        for stage,settings in enumerate(STAGES):
            self._seed(stage)
            result=self.runtime.fitter.fit_one(evaluate=lambda x:self.evaluate(stage,x),
                target=TARGETS[1 if stage==0 else 3],initial_level=self.plan['stage_starts'][stage]['initial'],
                bounds=BOUNDS,controls=dict(self.plan['fit'],total_seconds=self.remaining()),
                callback=lambda row:self.progress('rough_scalar_fit',stage=stage+1,root_row=row))
            original.write(self.out/f'stage{stage+1}'/'fit.json',result)
            original.require(result['converged'],'Rough scalar fit did not meet 0.005 target and fresh replay')
            accepted,exported=self._selected(stage,result)
            selected_receipts[stage]=copy.deepcopy(self.last['receipt'])
            parameters.append(dict(parameter='psi_'+str(settings['start_year']),**result['parameter']))
            results.append(result)
            if stage==0:
                self._handoff(accepted,self.last['candidate'])
                self.runtime.release_stage1();accepted=None
            original.write(self.out/f'stage{stage+1}'/'complete.json',dict(status='rough_matched',parameter=parameters[-1],
                selected_candidate=self.last['candidate'],original_horizon_gates_tested=False))
        models=self.prefix_fertility+accepted['fertility'][:2]
        original.require(len(models)==4,'Four calendar-lineage fertility rows required')
        rows=[]
        for i,(target,model) in enumerate(zip(self.plan['rows'],models)):
            value=float(model['period_tfr_topcode_adjusted']);gap=value-target['target'];weight=self.plan['weights'][i]
            rows.append(dict(target,model=value,gap=gap,weight=weight,loss_contribution=weight*gap*gap,
                source_stage=1 if i<2 else 2,local_index=i if i<2 else i-2,role='fitted' if weight else 'validation'))
        original.table(self.out/'fertility_fit.csv',rows);original.table(self.out/'parameters.csv',parameters)
        projection=[dict(row,source_stage=1 if i<2 else 2,local_index=i if i<2 else i-2,global_period=i)
                    for i,row in enumerate(self.prefix+accepted['rows'])]
        original.table(self.out/'final_projection.csv',projection)
        self.render(accepted,self.out/'diagnostics')
        actual=getattr(getattr(self.runtime,'rt',None),'total_native_calls',self.calls)
        original.require(actual==self.calls,'Controller/native actual-call accounting differs')
        receipt=dict(status='rough_matched',schema=SCHEMA,scientific_validation=False,production_ready=False,
            original_horizon_gates_tested=False,empirical_h16_seed_tested=True,terminal_diagnostics_gating=False,
            source_identity=self.plan['identity'],parameters=parameters,fit_table=rows,
            total_loss=sum(r['loss_contribution'] for r in rows),actual_policy_calls=self.calls,
            native_actual_policy_calls=actual,selected_candidate_pins=self.selected_candidate_pins,
            exact_2023_state=exported,terminal_passes_by_stage={str(i+1):selected_receipts[i]['terminal_passes'] for i in range(2)},
            classification='Sequential two-shock rough H16 fit only; original H24/H32 horizon and terminal certification untested',
            disclosure=self.plan['disclosure'])
        original.write(self.out/'complete.json',receipt)
        return receipt

    def _smoke(self):
        self._prepare();selected=[]
        for stage,settings in enumerate(STAGES):
            self._seed(stage)
            result=self.evaluate(stage,self.plan['baseline_psi'])
            original.require(result['certified'],'Reduced fixed-baseline exact-loop smoke candidate failed')
            accepted=self.last['reply']
            if stage==0:
                self._handoff(accepted,self.last['candidate']);self.runtime.release_stage1()
            selected.append(dict(stage=stage,psi=self.plan['baseline_psi'],model=result['model'],candidate=self.last['candidate']))
        actual=getattr(getattr(self.runtime,'rt',None),'total_native_calls',self.calls)
        original.require(actual==self.calls,'Smoke native call ledger differs')
        receipt=dict(status='execution_passed',schema=SCHEMA,smoke=True,scientific_validation=False,
            production_ready=False,empirical_fitted=False,empirical_h16_seed_tested=False,
            original_horizon_gates_tested=False,standard_diagnostics_tested=False,
            selected=selected,actual_policy_calls=self.calls,
            source_identity=self.plan['identity'],classification='Tiny fixed-baseline two-stage execution smoke only')
        original.write(self.out/'smoke_receipt.json',receipt)
        return receipt


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest');parser.add_argument('--base-plan');parser.add_argument('--output',required=True)
    parser.add_argument('--bridge-case');parser.add_argument('--source-template-manifest')
    parser.add_argument('--deadline-epoch',type=float)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare-base-plan',action='store_true')
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--prepare-smoke',action='store_true')
    mode.add_argument('--preflight',action='store_true')
    mode.add_argument('--smoke',action='store_true');mode.add_argument('--run',action='store_true')
    args=parser.parse_args(argv);out=Path(args.output)
    if args.prepare_base_plan:
        original.require(args.bridge_case is not None and args.source_template_manifest is not None,
                         '--bridge-case and --source-template-manifest required')
        return prepare_base_plan(args.bridge_case,out,args.source_template_manifest)
    if args.prepare or args.prepare_smoke:
        original.require(args.base_plan is not None,'--base-plan required for manifest preparation')
        plan=prepare_manifest(original.pin(args.base_plan),smoke=args.prepare_smoke,deadline_epoch=args.deadline_epoch)
        original.write(out,plan);return plan
    original.require(args.manifest is not None,'--manifest required')
    plan=json.loads(Path(args.manifest).read_text())
    preflight(plan)
    if args.preflight:
        receipt=preflight(plan);original.write(out/'preflight.json',receipt);return receipt
    original.require(plan['smoke'] is args.smoke,'Manifest smoke/run mode differs')
    original.require(not out.exists(),'Execution requires a fresh output directory')
    module=importlib.import_module('rough_two_shock_runtime')
    original.require(Path(module.__file__).resolve()==original.pinned(plan['source_pins']['rough_runtime']).resolve(),
                     'Rough runtime import path differs')
    runtime=module.build_runtime(plan=plan,output=out/'runtime',smoke=args.smoke)
    original.require(runtime.rt.total_native_calls==0,'Native constructor made calls')
    controller=Controller(plan,runtime,out)
    try:return controller.run()
    except BaseException as exc:
        original.write(out/'failure.json',dict(status='failed',schema=SCHEMA,exception=type(exc).__name__,reason=str(exc),
            policy_calls=controller.calls,native_calls=getattr(runtime.rt,'total_native_calls',None),
            completed_cases_retained=True,production_ready=False))
        raise

if __name__=='__main__':main()
