#!/usr/bin/env python3
"""One isolated continuation using authenticated roots from a failed local fit.

The retained panel/controller remain unchanged. Scalar fitting starts fresh at
an explicit numeric guess; only converged dated-root initializations are reused.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/model'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def pinned(pin):
    require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'A path/hash pin is required')
    require(sha(pin['path']) == pin['sha256'], 'Pinned file changed: ' + str(pin['path']))
    return Path(pin['path'])


def _finite_tree(x):
    if isinstance(x, list):
        return all(_finite_tree(v) for v in x)
    return type(x) in (int, float) and math.isfinite(float(x))


def _shape(x):
    if not isinstance(x, list):
        return ()
    if not x:
        return (0,)
    child = _shape(x[0])
    require(all(_shape(v) == child for v in x), 'Ragged saved root array')
    return (len(x),) + child


def _contract(plan):
    return dict(identity=plan['identity'], target_contract=plan['target_contract'],
                gates=plan['gates'], horizons=plan['horizons'], mode=plan['mode'],
                source_files=plan['source_files'], initial_psi=plan['initial_psi'],
                psi_bound_ratios=plan['psi_bound_ratios'], seed=plan['seed'],
                endpoint=plan['endpoint'], path=plan['path'], fit=plan['fit'], budget=plan['budget'])


def _root_path(root, name):
    p = Path(root) / name / 'root.json'
    p.resolve().relative_to(Path(root).resolve())
    return p


def authenticate(plan_pin, panel_pin, failed_pin, candidate, new_psi, index):
    """Verify source bundle without importing native model code."""
    plan_path = pinned(plan_pin); plan = read(plan_path)
    panel_path = pinned(panel_pin); panel_config = read(panel_path)
    require(plan.get('mode') == 'diagnostic' and plan.get('horizons') == [24, 32], 'Expected unchanged diagnostic 24/32 plan')
    require(plan.get('budget') == dict(total_seconds=21480, seed_seconds=1800, candidate_seconds=7200,
            endpoint_seconds=1800, mapping_seconds=1800, path_seconds=6000, render_seconds=900,
            maximum_policy_calls=20000) and plan['fit']['max_evaluations'] == 12 and
            plan['path']['max_evaluations'] == 12 and plan['endpoint']['max_evaluations'] == 48 and
            plan['fit']['fertility_tolerance'] == .005 and plan['target_contract']['rows'][3]['target'] == 1.64575,
            'Continuation contract/budget differs from authorized v7 limits')
    lo, hi = (plan['initial_psi'] * x for x in plan['psi_bound_ratios'])
    require(type(new_psi) in (int, float) and math.isfinite(new_psi) and lo <= new_psi <= hi,
            'New scalar start is outside unchanged bounds')
    require(type(index) is int and 0 <= index < 12, 'Continuation panel index must be 0..11')
    require(panel_config.get('schema') == 'estate_a_transition_panel_v1' and
            panel_config.get('identity') == plan['identity'], 'Panel source identity differs')
    require(read(pinned(panel_config['plan'])) == plan, 'Panel plan differs')
    failed_path = pinned(failed_pin); failed = read(failed_path)
    require(failed.get('schema') == 'current_estate_a_shock_candidate_v1' and
            failed.get('status') == 'failed' and failed.get('accepted') is False and
            failed.get('contract') == _contract(plan) and
            failed.get('contract_sha256') == fingerprint(_contract(plan)),
            'Prior failed candidate receipt is not an exact-contract failure')

    candidate = Path(candidate).resolve()
    require(candidate.is_dir() and candidate.name == 'candidate_0005',
            'Only the specified completed candidate_0005 is eligible; partial candidates are forbidden')
    complete_path = candidate / 'complete.json'
    complete = read(complete_path)
    require(complete.get('certified') is True and complete.get('mode') == 'diagnostic' and
            complete.get('production_ready') is False, 'Source candidate is not a completed diagnostic')
    source_psi = float(complete['psi'])
    require(math.isfinite(source_psi), 'Source candidate psi is nonfinite')
    require(complete.get('horizon_comparison', {}).get('passed') is True,
            'Source candidate horizon comparison failed')
    expected_model = float(complete['model'])
    expected_gap = expected_model - float(plan['target_contract']['rows'][3]['target'])
    require(math.isfinite(expected_model) and math.isfinite(float(complete['gap'])) and
            abs(float(complete['gap']) - expected_gap) <= 1e-12,
            'Completed candidate does not match the pinned fertility target')
    cp_path = candidate / 'state_2023_checkpoint' / 'checkpoint_receipt.json'
    checkpoint = read(cp_path)
    require(checkpoint.get('status') == 'complete' and checkpoint.get('identity') == plan['identity'] and
            checkpoint.get('psi') == source_psi and checkpoint.get('horizon') == 32 and
            checkpoint.get('candidate') == 5 and
            checkpoint.get('shock_fit_complete') is False,
            'Source candidate checkpoint identity/status differs')
    export_pin = checkpoint.get('export')
    require(isinstance(export_pin, dict) and set(export_pin) >= {'path', 'sha256'} and
            sha(export_pin['path']) == export_pin['sha256'], 'Source native checkpoint export pin differs')
    roots = {}
    for horizon in plan['horizons']:
        name = f'horizon_{horizon:03d}'
        path = _root_path(candidate, name)
        root = read(path)
        expected_gates = {'fiscal_replay', 'housing', 'mapping', 'market_replay', 'social_security'}
        require(root.get('converged') is True and set(root.get('gates', {})) == expected_gates and
                all(root['gates'].get(key) is True for key in expected_gates) and
                root.get('final', {}).get('mapping_valid') is True and
                root.get('final', {}).get('market_gate') is True and
                root.get('final', {}).get('fiscal_gate') is True,
                f'{name}: root did not converge with every gate')
        require(float(root.get('final_reproduction_max_abs', math.inf)) <= plan['gates']['final_reproduction_tolerance'],
                f'{name}: fresh replay residual exceeds unchanged tolerance')
        require(float(root.get('market_reproduction_max_abs', math.inf)) <= plan['gates']['market_tolerance'] and
                float(root.get('fiscal_reproduction_max_abs', math.inf)) <= plan['gates']['fiscal_tolerance'],
                f'{name}: market/fiscal replay residual fails')
        prices = root['final']['prices']; fiscal = root['final']['fiscal_values']; jac = root['final_jacobian']
        require(_finite_tree(prices) and _shape(prices) == (horizon,), f'{name}: price path shape/values invalid')
        require(_finite_tree(fiscal) and _shape(fiscal) == (horizon,), f'{name}: fiscal path shape/values invalid')
        require(_finite_tree(jac) and _shape(jac) == (2*horizon, 2*horizon), f'{name}: Jacobian shape/values invalid')
        require(all(float(v) > 0 for v in prices) and all(float(v) > 0 for v in fiscal),
                f'{name}: price/fiscal starts must be positive')
        roots[horizon] = dict(path=str(path), sha256=sha(path), data=root)
    # checkpoint is a candidate-level receipt and both horizon roots are pinned
    # independently; never glob for latest/partial roots.
    return dict(plan=plan, plan_pin=plan_pin, panel_config=panel_config,
                panel_pin=panel_pin, failed=failed, failed_pin=failed_pin,
                candidate=str(candidate), complete_pin=dict(path=str(complete_path), sha256=sha(complete_path)),
                checkpoint_pin=dict(path=str(cp_path), sha256=sha(cp_path)), source_psi=source_psi,
                roots=roots, new_psi=float(new_psi), index=index)


def create_manifest(args):
    wrapper = Path(__file__).resolve()
    evidence = authenticate(args.plan_pin, args.panel_config_pin, args.failed_receipt_pin,
                            args.source_candidate, args.start_psi, args.index)
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    wrapper_pin = dict(path=str(wrapper), sha256=sha(wrapper))
    continuation_id = fingerprint(dict(wrapper=wrapper_pin, plan=evidence['plan_pin'],
        panel_config=evidence['panel_pin'], failed=evidence['failed_pin'],
        candidate=evidence['candidate'], complete=evidence['complete_pin'],
        checkpoint=evidence['checkpoint_pin'], roots={h:{'path':r['path'],'sha256':r['sha256']}
                                                       for h,r in evidence['roots'].items()},
        source_psi=evidence['source_psi'], start_psi=evidence['new_psi'], index=evidence['index']))
    config = copy.deepcopy(evidence['panel_config'])
    config['guesses'][args.index]['psi'] = evidence['new_psi']
    config['guesses'][args.index]['ratio'] = evidence['new_psi'] / evidence['plan']['initial_psi']
    config['continuation_identity'] = continuation_id
    config['local_indices'] = [args.index]
    config['cluster_indices'] = []
    config['search_strategy'] = 'One isolated local scalar fit from an explicit numeric start with authenticated warm roots; optimizer state restarts'
    config['warm_manifest'] = dict(path=str(out/'manifest.json'), wrapper_sha256=wrapper_pin['sha256'])
    config_path = out/'panel_config.json'; write(config_path, config)
    evidence.update(schema='estate_birth_transition_local_continuation_v1', wrapper=wrapper_pin,
                    continuation_identity=continuation_id, output=str(out),
                    panel_config_path=str(config_path),
                    panel_config_pin=dict(path=str(config_path), sha256=sha(config_path)))
    # Replace loaded roots with compact pins; load them again only after authentication.
    evidence['roots'] = {h:{'path':r['path'],'sha256':r['sha256']} for h,r in evidence['roots'].items()}
    write(out/'manifest.json', evidence)
    return evidence


def load_manifest(path):
    manifest_path = Path(path).resolve(); manifest = read(manifest_path)
    require(manifest.get('schema') == 'estate_birth_transition_local_continuation_v1', 'Continuation manifest schema differs')
    require(sha(manifest['wrapper']['path']) == manifest['wrapper']['sha256'] == sha(__file__), 'Continuation wrapper source changed')
    pinned(manifest['plan_pin']); pinned(manifest['panel_pin']); pinned(manifest['failed_pin'])
    pinned(manifest['complete_pin']); pinned(manifest['checkpoint_pin']); pinned(manifest['panel_config_pin'])
    for pin in manifest['roots'].values(): pinned(pin)
    plan = read(manifest['plan_pin']['path'])
    require(manifest['continuation_identity'] == fingerprint(dict(wrapper=manifest['wrapper'], plan=manifest['plan_pin'],
        panel_config=manifest['panel_pin'], failed=manifest['failed_pin'], candidate=manifest['candidate'],
        complete=manifest['complete_pin'], checkpoint=manifest['checkpoint_pin'], roots=manifest['roots'],
        source_psi=manifest['source_psi'], start_psi=manifest['new_psi'], index=manifest['index'])),
        'Continuation identity fingerprint differs')
    require(plan['identity'] == read(manifest['checkpoint_pin']['path'])['identity'], 'Checkpoint/runtime identity differs')
    # Re-authenticate files and semantics before each native entrypoint.
    check = authenticate(manifest['plan_pin'], manifest['panel_pin'], manifest['failed_pin'],
                         manifest['candidate'], manifest['new_psi'], manifest['index'])
    require(manifest['plan'] == check['plan'] and
            check['complete_pin'] == manifest['complete_pin'] and check['checkpoint_pin'] == manifest['checkpoint_pin'] and
            {str(h):{'path':x['path'],'sha256':x['sha256']} for h,x in check['roots'].items()} == manifest['roots'],
            'Warm source evidence changed since manifest creation')
    require(manifest['output'] == str(manifest_path.parent), 'Manifest moved from its fresh output directory')
    return manifest, check


def _adapter_factory(control, roots, source_psi):
    base = control.NativeAdapter
    class WarmAdapter(base):
        def __init__(self, runtime, plan):
            super().__init__(runtime, plan)
            require(self.identity() == plan['identity'], 'Fresh runtime identity differs from authenticated roots')
            for horizon, record in roots.items():
                path = pinned({k:record[k] for k in ('path','sha256')})
                root = read(path)
                self.retain_warm(int(horizon), source_psi, root, path.parent)
    return base, WarmAdapter


def verify_smoke_outputs(smoke_dir, plan, source_psi, source_models, reply, result):
    require(result.get('certified') is True and reply is not None, 'Warm-start smoke candidate did not complete')
    require(reply.get('identity') == plan['identity'] and reply.get('psi') == source_psi and
            reply.get('horizon') == 32 and reply.get('root_pass') is True and reply.get('replay_pass') is True and
            reply.get('accounting_valid') is True and reply.get('stationary_pass') is True and
            reply.get('horizon_comparison', {}).get('passed') is True,
            'Warm-start smoke final horizon identity/replay/gates failed')
    comparison = result.get('horizon_comparison', {})
    require(comparison.get('passed') is True, 'Smoke 24/32 horizon comparison failed')
    models = result.get('payload', {}).get('models')
    require(isinstance(models, list) and len(models) >= 4 and isinstance(source_models, list) and len(source_models) >= 4,
            'Fresh smoke and retained source need four fertility values')
    absolute_gaps = [abs(float(a)-float(b)) for a,b in zip(models[:4],source_models[:4])]
    require(max(absolute_gaps) <= plan['fit']['fertility_tolerance']/5,
            'Fresh same-psi smoke fertility differs from retained source beyond the absolute fertility gate')
    for horizon in (24, 32):
        folder = Path(smoke_dir) / 'candidate_0001' / f'horizon_{horizon:03d}'
        root = read(folder/'root.json'); mapped = read(folder/'latest_completed.json')
        expected_root_gates = {'fiscal_replay','housing','mapping','market_replay','social_security'}
        require(root.get('converged') is True and set(root.get('gates',{})) == expected_root_gates and
                all(root['gates'].get(key) is True for key in expected_root_gates) and
                root.get('final', {}).get('mapping_valid') is True and
                root.get('final', {}).get('market_gate') is True and
                root.get('final', {}).get('fiscal_gate') is True and
                float(root.get('final_reproduction_max_abs', math.inf)) <=
                plan['gates']['final_reproduction_tolerance'],
                f'Warm-start smoke horizon {horizon} did not pass all root gates/fresh replay')
        required_mapping_gates = {'mass','policy_reproduction','projection','dated_audits'}
        require(mapped.get('accounting_valid') is True and set(mapped.get('gates',{})) == required_mapping_gates and
                all(mapped['gates'].get(key) is True for key in required_mapping_gates),
                f'Warm-start smoke horizon {horizon} accounting/mapping gates failed')


def smoke(manifest_path, output_path):
    manifest, evidence = load_manifest(manifest_path)
    from experiments.birth_count_choice import transition as driver
    driver.preflight(manifest['plan'])
    original, wrapped = _adapter_factory(driver.controller, evidence['roots'], evidence['source_psi'])
    driver.controller.NativeAdapter = wrapped
    smoke_dir = Path(output_path).resolve()
    smoke_dir.mkdir(exist_ok=False)
    try:
        runtime = driver.runtime_module().CurrentEstateARuntime.from_handoff(manifest['plan']['handoff'], smoke_dir/'runtime')
        effective = copy.deepcopy(manifest['plan']); effective['fit_start_psi'] = manifest['source_psi']
        effective['budget'] = dict(effective['budget'], total_seconds=min(3600, effective['budget']['total_seconds']))
        adapter = driver.controller.NativeAdapter(runtime, effective)
        runner = driver.controller.Controller(effective, adapter, smoke_dir)
        runner.prepare()
        result = runner.evaluate(manifest['source_psi'])
        reply = runner.last
        source_models = read(manifest['complete_pin']['path'])['payload']['models']
        verify_smoke_outputs(smoke_dir, manifest['plan'], manifest['source_psi'], source_models, reply, result)
        root_pins = {}
        mapped_pins = {}
        for h in (24, 32):
            folder = smoke_dir/'candidate_0001'/f'horizon_{h:03d}'
            root_pins[str(h)] = dict(path=str(folder/'root.json'), sha256=sha(folder/'root.json'))
            mapped_pins[str(h)] = dict(path=str(folder/'latest_completed.json'), sha256=sha(folder/'latest_completed.json'))
        complete_pin = dict(path=str(smoke_dir/'candidate_0001/complete.json'),
                            sha256=sha(smoke_dir/'candidate_0001/complete.json'))
        checkpoint_path = smoke_dir/'candidate_0001/state_2023_checkpoint/checkpoint_receipt.json'
        checkpoint_pin = dict(path=str(checkpoint_path), sha256=sha(checkpoint_path))
        receipt = dict(schema='estate_birth_transition_warm_smoke_v1', passed=True,
            continuation_identity=manifest['continuation_identity'], wrapper_sha256=manifest['wrapper']['sha256'],
            source_psi=manifest['source_psi'], identity=reply['identity'], horizons=[24,32],
            root_pins=root_pins, mapped_pins=mapped_pins, complete_pin=complete_pin,
            checkpoint_pin=checkpoint_pin, fertility_absolute_gaps=[
                abs(float(a)-float(b)) for a,b in zip(result['payload']['models'][:4],source_models[:4])],
            policy_calls=runner.policy_calls,
            latest_candidate=result, output=str(smoke_dir))
        receipt_path = smoke_dir/'smoke_receipt.json'; write(receipt_path, receipt)
        return receipt
    finally:
        driver.controller.NativeAdapter = original


def run(manifest_path, output_path, smoke_pin):
    manifest, evidence = load_manifest(manifest_path)
    smoke_path = pinned(smoke_pin); smoke_receipt = read(smoke_path)
    require(smoke_receipt.get('schema') == 'estate_birth_transition_warm_smoke_v1' and
            smoke_receipt.get('passed') is True and smoke_receipt.get('continuation_identity') == manifest['continuation_identity'] and
            smoke_receipt.get('wrapper_sha256') == manifest['wrapper']['sha256'] and
            smoke_receipt.get('source_psi') == manifest['source_psi'] and
            smoke_receipt.get('identity') == manifest['plan']['identity'] and
            smoke_receipt.get('horizons') == [24,32], 'Passed exact-source warm smoke receipt required')
    require(smoke_receipt.get('complete_pin') and smoke_receipt.get('checkpoint_pin'),
            'Smoke receipt must pin its completed candidate and checkpoint')
    require(smoke_receipt.get('fertility_absolute_gaps') and
            max(smoke_receipt['fertility_absolute_gaps']) <= manifest['plan']['fit']['fertility_tolerance']/5,
            'Smoke same-psi fertility comparison fails retained horizon tolerance')
    smoke_complete = read(pinned(smoke_receipt['complete_pin']))
    require(smoke_complete.get('certified') is True and smoke_complete.get('psi') == manifest['source_psi'] and
            smoke_complete.get('horizon_comparison',{}).get('passed') is True,
            'Pinned smoke completion record no longer passes')
    pinned(smoke_receipt['checkpoint_pin'])
    smoke_checkpoint = read(smoke_receipt['checkpoint_pin']['path'])
    require(smoke_checkpoint.get('identity') == manifest['plan']['identity'] and
            smoke_checkpoint.get('psi') == manifest['source_psi'] and smoke_checkpoint.get('status') == 'complete' and
            smoke_checkpoint.get('candidate') == 1 and smoke_checkpoint.get('horizon') == 32,
            'Smoke checkpoint identity differs')
    export_pin = smoke_checkpoint.get('export')
    require(isinstance(export_pin, dict) and sha(export_pin['path']) == export_pin['sha256'],
            'Smoke native export checkpoint changed')
    for h in ('24','32'):
        root_path = pinned(smoke_receipt['root_pins'][h]); root = read(root_path)
        expected_gates = {'fiscal_replay', 'housing', 'mapping', 'market_replay', 'social_security'}
        require(root.get('converged') is True and set(root.get('gates',{})) == expected_gates and
                all(root['gates'].get(key) is True for key in expected_gates) and
                root.get('final',{}).get('mapping_valid') is True and
                root.get('final',{}).get('market_gate') is True and
                root.get('final',{}).get('fiscal_gate') is True and
                all(float(v)>0 for v in root.get('final',{}).get('prices',[])) and
                all(float(v)>0 for v in root.get('final',{}).get('fiscal_values',[])) and
                float(root.get('final_reproduction_max_abs', math.inf)) <= manifest['plan']['gates']['final_reproduction_tolerance'],
                'Smoke root gates/replay changed after receipt')
        mapped = read(pinned(smoke_receipt['mapped_pins'][h]))
        require(mapped.get('accounting_valid') is True and
                set(mapped.get('gates',{})) == {'mass','policy_reproduction','projection','dated_audits'} and
                all(mapped['gates'].get(k) is True for k in ('mass','policy_reproduction','projection','dated_audits')),
                'Smoke mapping/accounting gates changed after receipt')
    run_dir = Path(output_path).resolve(); run_dir.mkdir(parents=True, exist_ok=False)
    from experiments.birth_count_choice import transition_panel as panel
    base, wrapped = _adapter_factory(panel.control, evidence['roots'], evidence['source_psi'])
    panel.control.NativeAdapter = wrapped
    try:
        # The panel still owns plan validation, Controller.run, output writing,
        # and standard collector-compatible receipts. Optimizer starts fresh.
        result = panel.evaluate_candidate(manifest['plan'], manifest['new_psi'], manifest['index'],
                                          run_dir, panel_config=manifest['panel_config_pin'])
        write(run_dir/'continuation_manifest_link.json', dict(continuation_identity=manifest['continuation_identity'],
            manifest_sha256=sha(manifest_path), smoke_receipt=smoke_pin, numeric_start=manifest['new_psi'],
            scalar_optimizer_restarted=True, warm_roots=manifest['roots']))
        return result
    finally:
        panel.control.NativeAdapter = base


def main():
    ap = argparse.ArgumentParser()
    modes = ap.add_mutually_exclusive_group(required=True)
    modes.add_argument('--prepare', action='store_true')
    modes.add_argument('--smoke', action='store_true')
    modes.add_argument('--run', action='store_true')
    ap.add_argument('--manifest')
    ap.add_argument('--plan-pin', type=json.loads)
    ap.add_argument('--panel-config-pin', type=json.loads)
    ap.add_argument('--failed-receipt-pin', type=json.loads)
    ap.add_argument('--source-candidate')
    ap.add_argument('--start-psi', type=float)
    ap.add_argument('--index', type=int, default=0)
    ap.add_argument('--output')
    ap.add_argument('--smoke-receipt-pin', type=json.loads)
    args = ap.parse_args()
    if args.prepare:
        for key in ('plan_pin','panel_config_pin','failed_receipt_pin','source_candidate','start_psi','output'):
            require(getattr(args,key) is not None, '--'+key.replace('_','-')+' required')
        result=create_manifest(args)
        result = dict(status='prepared', manifest=str(Path(args.output).resolve()/'manifest.json'),
                      continuation_identity=result['continuation_identity'], numeric_start=result['new_psi'])
    elif args.smoke:
        require(args.manifest and args.output, '--manifest and --output required')
        result=smoke(args.manifest,args.output)
    else:
        require(args.manifest and args.smoke_receipt_pin and args.output,
                '--manifest, --smoke-receipt-pin and --output required')
        require(args.smoke_receipt_pin is not None, '--smoke-receipt-pin required for --run')
        result=run(args.manifest,args.output,args.smoke_receipt_pin)
    print(json.dumps({k:result[k] for k in ('status','passed','accepted','continuation_identity','numeric_start','source_psi','output') if k in result},
                     sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
