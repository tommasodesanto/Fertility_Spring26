#!/usr/bin/env python3
"""Evaluate one permanent shock guess per worker; collect authenticated panels.

The maintained Controller prepares a fresh reference/seed and runs its original
scalar fitter from a distinct numerical start on each worker.
Each worker repeats setup deliberately; no historical state or seed is reused.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.birth_count_choice import transition as driver

SCHEMA = 'current_estate_a_shock_candidate_v1'
COLLECTION_SCHEMA = 'current_estate_a_shock_panel_v1'
control = driver.controller


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def table(path, rows):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def contract_identity(plan):
    return dict(identity=plan['identity'], target_contract=plan['target_contract'], gates=plan['gates'],
        horizons=plan['horizons'], mode=plan['mode'], source_files=plan['source_files'],
        initial_psi=plan['initial_psi'], psi_bound_ratios=plan['psi_bound_ratios'],
        seed=plan['seed'], endpoint=plan['endpoint'], path=plan['path'], fit=plan['fit'], budget=plan['budget'])


def validate_candidate(plan, psi, index, panel_config):
    driver.preflight(plan)
    control.require(plan.get('execution_enabled') is True,'Pinned plan must enable candidate execution')
    control.require(plan['mode']=='diagnostic' and plan['horizons']==[24,32],
        'Current shock panel requires diagnostic 24/32 horizons')
    control.require(plan['budget']['total_seconds']<=21480,'Candidate budget exceeds six-hour worker allowance')
    control.require(type(index) is int and 0<=index<48,'Candidate index must lie in 0..47')
    control.require(type(psi) in (int,float) and math.isfinite(psi) and
        plan['initial_psi']*.01<=psi<=plan['initial_psi']*2.,'Shock guess outside retained bounds')
    config=json.loads(control.pinned(panel_config).read_text())
    control.require(config.get('schema')=='estate_a_transition_panel_v1','Pinned panel configuration required')
    control.require(config['identity']==plan['identity'],'Panel native identity differs')
    control.require(json.loads(control.pinned(config['plan']).read_text())==plan,'Panel plan differs')
    control.require(config['panel_source']==driver.pin(__file__),'Panel orchestration source differs')
    guesses=config['guesses']
    control.require(len(guesses)==12 and [r['index'] for r in guesses]==list(range(12)),
        'Exactly twelve deterministic panel guesses required')
    control.require(index<len(guesses) and guesses[index]['psi']==float(psi),'Panel index/psi differs')


def evaluate_candidate(plan, psi, index, output, *, panel_config):
    validate_candidate(plan,psi,index,panel_config)
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    control.require(not (out/'candidate_receipt.json').exists(),'Refusing to overwrite candidate receipt')
    started=time.monotonic()
    rt=driver.runtime_module().CurrentEstateARuntime.from_handoff(plan['handoff'],out/'runtime')
    effective_plan=copy.deepcopy(plan);effective_plan['fit_start_psi']=float(psi)
    control.write(out/'effective_plan.json',effective_plan)
    adapter=control.NativeAdapter(rt,effective_plan)
    runner=control.Controller(effective_plan,adapter,out)
    common=dict(schema=SCHEMA,index=index,psi=float(psi),
        contract=contract_identity(plan),contract_sha256=fingerprint(contract_identity(plan)),
        panel_source_sha256=control.sha(__file__),panel_config_sha256=panel_config['sha256'],production_ready=False,scientific_validation=False,
        setup_strategy='fresh current reference and five-map seed independently on every worker',
        setup_redundant_across_workers=True,scalar_optimizer_used=True,
        search_strategy='independent scalar fit from selected numerical start; baseline and bounds unchanged',
        effective_plan_sha256=control.sha(out/'effective_plan.json'))
    try:
        fitted=runner.run()
        reply=runner.last
        final_psi=float(fitted['fitted_parameter']['estimate'])
        models=[row['period_tfr_topcode_adjusted'] for row in reply['fertility'][:4]]
        target=plan['target_contract']['rows'][3]['target']
        summary=dict(certified=fitted['status']=='matched',psi=final_psi,model=models[3],gap=models[3]-target,
            loss_contribution=(models[3]-target)**2,payload=dict(models=models))
        seed_calls=runner.seed.get('policy_calls',0)
        reference_receipt=out/'native_reference/reference_reconstruction.json'
        reference_calls=json.loads(reference_receipt.read_text())['policy_calls'] if reference_receipt.is_file() else None
        setup_calls=None if reference_calls is None else seed_calls+reference_calls
        _,fitter=control.original_modules()
        models=summary.get('payload',{}).get('models',[None]*4)
        rows=fitter.fit_rows(plan['target_contract']['rows'],models,'one_permanent')
        table(out/'fertility_fit.csv',rows)
        lower,upper=[plan['initial_psi']*v for v in plan['psi_bound_ratios']]
        parameter=dict(parameter='psi_child',estimate=final_psi,lower=lower,upper=upper,
            near_bound=min(final_psi-lower,upper-final_psi)<=.01*(upper-lower),
            role='fitted permanent 2007 level; finite-horizon diagnostic')
        table(out/'shock_parameters.csv',[parameter])
        _,case=driver.validate_handoff(plan['handoff'])
        (out/'baseline_parameters.csv').write_bytes((case/'parameters.csv').read_bytes())
        files={name:control.sha(out/name) for name in ('fertility_fit.csv','shock_parameters.csv','baseline_parameters.csv','effective_plan.json')}
        # Authenticate scalar path/root/accounting receipts without hashing policy caches.
        path_receipts={str(path.relative_to(out)):control.sha(path) for path in
            sorted(out.glob('candidate_*/*.json'))+sorted(out.glob('candidate_*/**/*.json')) if path.name in
            ('root.json','native_record.json','complete.json','horizon_comparison.json',
             'state_2023_horizon_comparison.json','checkpoint_receipt.json','failure.json')}
        reply=runner.last
        gates=None if reply is None else {key:reply[key] for key in
            ('root_pass','replay_pass','accounting_valid','stationary_pass','market_maximum_residual',
             'fiscal_maximum_residual','replay_maximum_gap','stationary_renewal_gap','terminal_passes',
             'horizon_comparison','state_horizon_comparison')}
        receipt=dict(common,status='accepted_diagnostic' if summary.get('certified') is True else 'rejected',
            accepted=summary.get('certified') is True,summary=summary,fit_rows=rows,parameter=parameter,
            gates=gates,files=files,path_receipts=path_receipts,setup_policy_calls=setup_calls,fit_receipt=fitted,
            actual_policy_calls=runner.policy_calls,candidate_policy_calls=None if setup_calls is None else runner.policy_calls-setup_calls,
            elapsed_seconds=time.monotonic()-started,shock_fit_complete=True,
            classification='Independent scalar-fit start; matched finite-horizon diagnostic, production withheld')
    except Exception as exc:
        receipt=dict(common,status='failed',accepted=False,error=control.failure_receipt(exc,runner),
            elapsed_seconds=time.monotonic()-started,shock_fit_complete=False)
        control.write(out/'candidate_receipt.json',receipt)
        raise
    control.write(out/'candidate_receipt.json',receipt)
    return receipt


def collect_candidates(inputs, output):
    receipts=[];expected=None;seen=set();rejected=[]
    for item in inputs:
        path=Path(item);path=path/'candidate_receipt.json' if path.is_dir() else path
        receipt=json.loads(path.read_text())
        control.require(receipt.get('schema')==SCHEMA,'Wrong candidate receipt schema')
        digest=fingerprint(receipt['contract'])
        control.require(digest==receipt['contract_sha256'],'Candidate contract fingerprint differs')
        identity=(digest,receipt['panel_source_sha256'],receipt['panel_config_sha256'])
        if expected is None:expected=identity
        control.require(identity==expected,'Mixed target/source/numerical fingerprints in panel')
        key=receipt['index'];control.require(key not in seen,'Duplicate panel candidate index');seen.add(key)
        control.require(0<=key<48 and math.isfinite(receipt['psi']),'Invalid candidate identity')
        if receipt.get('accepted') is not True:
            rejected.append(dict(index=key,psi=receipt['psi'],status=receipt['status']));continue
        control.require(receipt['status']=='accepted_diagnostic' and receipt['summary']['certified'] is True,
            'Accepted status conflicts with retained controller certificate')
        for relative,digest in {**receipt['files'],**receipt['path_receipts']}.items():
            file=(path.parent/relative).resolve();file.relative_to(path.parent.resolve())
            control.require(control.sha(file)==digest,'Candidate output receipt changed: '+relative)
        effective=json.loads((path.parent/'effective_plan.json').read_text())
        control.require(effective.get('fit_start_psi')==receipt['psi'] and
            contract_identity(effective)==receipt['contract'],'Effective numerical start changed baseline contract')
        control.require(receipt['fit_receipt']['status']=='matched' and
            receipt['fit_receipt']['fitted_parameter']['estimate']==receipt['parameter']['estimate']==receipt['summary']['psi'],
            'Selected parameter differs from fresh reproduced scalar fit')
        rows=receipt['fit_rows'];control.require(len(rows)==4,'Complete four-window fit required')
        _,fitter=control.original_modules()
        rebuilt=fitter.fit_rows(receipt['contract']['target_contract']['rows'],
            [row['model'] for row in rows],'one_permanent')
        control.require(rows==rebuilt and rows[3]['gap']==receipt['summary']['gap'] and
            rows[3]['loss_contribution']==receipt['summary']['loss_contribution'], 'Candidate fit arithmetic differs')
        gates=receipt['gates'];limits=receipt['contract']['gates']
        control.require(all(gates[key] is True for key in ('root_pass','replay_pass','accounting_valid','stationary_pass')) and
            gates['market_maximum_residual']<=limits['market_tolerance'] and
            gates['fiscal_maximum_residual']<=limits['fiscal_tolerance'] and
            gates['replay_maximum_gap']<=limits['final_reproduction_tolerance'] and
            gates['stationary_renewal_gap']<=limits['stationary_renewal_tolerance'] and
            gates['horizon_comparison']['passed'] is True,'Candidate retained numerical gates failed')
        receipts.append(receipt)
    control.require(bool(expected),'At least one candidate receipt required')
    receipts.sort(key=lambda r:(abs(r['fit_rows'][3]['gap']),r['index']))
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    ranking=[dict(rank=i+1,index=r['index'],initial_psi=r['psi'],psi=r['parameter']['estimate'],model=r['fit_rows'][3]['model'],
        target=r['fit_rows'][3]['target'],gap=r['fit_rows'][3]['gap'],
        loss_contribution=r['fit_rows'][3]['loss_contribution']) for i,r in enumerate(receipts)]
    all_rows=[dict(index=r['index'],initial_psi=r['psi'],psi=r['parameter']['estimate'],**row) for r in receipts for row in r['fit_rows']]
    if ranking:table(out/'ranking.csv',ranking);table(out/'fertility_fit.csv',all_rows)
    result=dict(schema=COLLECTION_SCHEMA,contract_sha256=expected[0],panel_source_sha256=expected[1],
        panel_config_sha256=expected[2],accepted_candidates=receipts,rejected_candidates=rejected,ranking=ranking,
        scientific_validation=False,production_ready=False,shock_fit_complete=False,
        classification='Authenticated parallel scalar-fit panel; matched chains retain original fresh replay, production certification outstanding')
    control.write(out/'collection.json',result)
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collect',action='store_true')
    parser.add_argument('--inputs',nargs='+')
    parser.add_argument('--plan');parser.add_argument('--psi',type=float);parser.add_argument('--index',type=int)
    parser.add_argument('--panel-config');parser.add_argument('--panel-config-sha256')
    parser.add_argument('--output',required=True)
    args=parser.parse_args(argv)
    if args.collect:
        if not args.inputs:parser.error('--collect requires --inputs receipt paths')
        receipt=collect_candidates(args.inputs,args.output)
        print(json.dumps(dict(accepted=len(receipt['ranking']),rejected=len(receipt['rejected_candidates']))))
    else:
        if any(v is None for v in (args.plan,args.psi,args.index,args.panel_config,args.panel_config_sha256)):
            parser.error('--plan, --psi, --index, --panel-config and --panel-config-sha256 required')
        receipt=evaluate_candidate(json.loads(Path(args.plan).read_text()),args.psi,args.index,args.output,
            panel_config=dict(path=args.panel_config,sha256=args.panel_config_sha256))
        print(json.dumps(dict(status=receipt['status'],index=receipt['index'],psi=receipt['psi'])))


if __name__=='__main__':main()
