#!/usr/bin/env python3
"""Current Estate-A one-permanent-shock plan and unchanged controller entrypoint.

Preparation loads a saved case without solving. Execution delegates every seed,
root, fitting, checkpoint and acceptance operation to one_shock_floor unchanged.
Finite horizons are diagnostics; no policy-closure or production claim is implied.
"""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SHARED = HERE.parent / 'transition_readiness'
sys.path.insert(0, str(SHARED))
import one_shock_floor as controller

HANDOFF_SCHEMA = 'current_estate_a_transition_handoff_v1'
SAVED_FILES = ('metadata.json', 'native_result.npz', 'input_contract.json',
               'target_fit.csv', 'parameters.csv', 'native/phase_b_ge/selected_root/closure.json')
BUDGET_KEYS = ('total_seconds', 'seed_seconds', 'candidate_seconds', 'endpoint_seconds',
               'mapping_seconds', 'path_seconds', 'render_seconds', 'maximum_policy_calls')


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=controller.sha(path))


def build_handoff(saved_case, source_paths):
    """Pin explicitly selected source files and the complete saved-case contract."""
    case = Path(saved_case).resolve()
    # Relative identities make a staged, root-preserving package portable.
    case.relative_to(ROOT)
    sources = {}
    for item in source_paths:
        path = Path(item).resolve()
        sources[str(path.relative_to(ROOT))] = controller.sha(path)
    controller.require(bool(sources), 'Explicit current-runtime source pins required')
    plots = sorted((case / 'standard_diagnostics').glob('*.png'))
    controller.require(len(plots) == 17, 'Saved case must have the retained 17 standard plots')
    return dict(schema=HANDOFF_SCHEMA, saved_case=str(case.relative_to(ROOT)),
        saved_files={name: controller.sha(case/name) for name in SAVED_FILES},
        source_pins=sources,
        standard_plot_pins={str(p.relative_to(case)): controller.sha(p) for p in plots})


def runtime_module():
    # Keep the current sibling runtime separate from the older floor constructor.
    sys.path.insert(0, str(ROOT/'code/model'))
    return importlib.import_module('experiments.birth_count_choice.transition_runtime')


def validate_handoff(item):
    path = controller.pinned(item)
    data = json.loads(path.read_text())
    controller.require(data.get('schema') == HANDOFF_SCHEMA, 'Current Estate-A handoff required')
    case = (ROOT/data['saved_case']).resolve()
    case.relative_to(ROOT)
    controller.require(set(SAVED_FILES).issubset(data['saved_files']), 'Complete saved-case pins required')
    for name, digest in data['saved_files'].items():
        target = (case/name).resolve()
        target.relative_to(case)
        controller.require(controller.sha(target) == digest, 'Saved case changed: '+name)
    for name, digest in data.get('standard_plot_pins', {}).items():
        target = (case/name).resolve()
        target.relative_to(case)
        controller.require(controller.sha(target) == digest, 'Standard plot changed: '+name)
    controller.require(bool(data['source_pins']), 'Current runtime source pins required')
    for name, digest in data['source_pins'].items():
        target = (ROOT/name).resolve()
        target.relative_to(ROOT)
        controller.require(controller.sha(target) == digest, 'Current source changed: '+name)
    return data, case


def build_plan(handoff, runtime, *, mode, horizons, budget, fit_max_evaluations,
               path_max_evaluations, endpoint_max_evaluations,
               fiscal_relaxation_authorized, execution_enabled=False):
    """Build the retained controller contract with explicit numerical budgets."""
    data, case = validate_handoff(handoff)
    controller.require(set(budget) == set(BUDGET_KEYS), 'Every stage budget must be explicit')
    original, _ = controller.original_modules()
    blocks = ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/inputs/empirical_blocks.csv'
    annual = ROOT/'output/model/e5f_matched_pf_20260909a/path_pilot_20260910/fertility_data/annual_fertility_2007_2023.csv'
    contract = original.target_contract(blocks, annual)
    controller.require(abs(contract['rows'][3]['target'] - 1.64575) < 1e-12,
                       'Retained final-window target changed')
    plots = sorted(p.name for p in (case/'standard_diagnostics').glob('*.png'))
    identity = runtime.identity()
    controller.require(identity['source_pins'] == data['source_pins'], 'Runtime source identity differs from handoff')
    plan = dict(schema='current_floor_one_permanent_v1',
        current_context_schema='current_estate_a_one_permanent_v1', kind='one_permanent', start_year=2007,
        mode=mode, horizons=list(horizons), execution_enabled=execution_enabled,
        fiscal_relaxation_authorized=fiscal_relaxation_authorized, gates=dict(controller.GATES),
        identity=identity, handoff=handoff,
        source_files=dict(controller=pin(SHARED/'one_shock_floor.py'),
            original_estimator=pin(controller.PINNED/'run_e5f_preference_estimation.py'),
            original_fitter=pin(controller.PINNED/'e5f_preference_shock_fit.py'),
            runtime=pin(HERE/'transition_runtime.py'), driver=pin(__file__)),
        target_contract=contract, budget=dict(budget), seed=dict(horizon=12, perturbed_date=5, log_step=1e-5),
        initial_psi=float(runtime.P.psi_child), psi_bound_ratios=[.01, 2.], standard_plot_names=plots,
        fit=dict(max_evaluations=fit_max_evaluations, log_difference_step=.01, fertility_tolerance=.005,
            max_log_step=.15, damping=.7, max_condition_number=1e8, worsening_factor=1.5, reproduction_tolerance=1e-8),
        path=dict(max_evaluations=path_max_evaluations, price_bound_ratios=[.05,20.],
            pension_bound_ratios=[.05,20.], max_log_step=.15, damping=.7),
        endpoint=dict(max_evaluations=endpoint_max_evaluations, price_bound_ratios=[.05,20.],
            max_log_step=.15, damping=.7, slope=1.),
        policy_contract_closed=False,
        context_disclosure='Current chain-13 one-birth Estate-A, post-interest transactions, soft credit; provisional estate/recipient accounting retained.',
        fitted_window='2020–2023', validation_windows=['2008–2011','2012–2015','2016–2019'],
        production_disclosure='Finite-horizon diagnostics do not close policy or numerical production gates.')
    if list(horizons) == [6,8]:
        plan['native_smoke_psi'] = plan['initial_psi']*1.001
    preflight(plan)
    return plan


def preflight(plan):
    controller.require(plan.get('current_context_schema') == 'current_estate_a_one_permanent_v1',
                       'Current-context controller plan required')
    controller.require(not any(plan.get(key) is not None for key in
        ('prepared_native_inputs','prepared_consumer_compatibility','diagnostic_measurement_reuse')),
        'Current adapter requires fresh reference and measured seed; historical reuse is forbidden')
    data, _ = validate_handoff(plan['handoff'])
    controller.require(plan['identity']['source_pins'] == data['source_pins'], 'Plan source pins differ')
    controller.require(plan['source_files']['driver'] == pin(__file__), 'Transition driver changed')
    controller.require(plan['source_files']['runtime'] == pin(HERE/'transition_runtime.py'), 'Current runtime changed')
    return controller.preflight(plan)


def execute(plan, output, *, native_smoke=False):
    receipt = preflight(plan)
    controller.require(plan.get('execution_enabled') is True, 'Pinned plan must explicitly enable execution')
    rt = runtime_module().CurrentEstateARuntime.from_handoff(plan['handoff'], Path(output)/'runtime')
    adapter = controller.NativeAdapter(rt, plan)
    runner = controller.Controller(plan, adapter, output)
    try:
        if native_smoke:
            controller.require(plan['mode'] == 'diagnostic' and plan['horizons'] == [6,8] and
                plan['path']['max_evaluations'] == 6, 'Smoke requires explicit diagnostic 6/8 horizons and six path maps')
            controller.require(plan.get('native_smoke_psi') == plan['initial_psi']*1.001,
                'Smoke must disclose its 0.1-percent preference perturbation')
            runner.prepare()
            started = time.monotonic()
            reply = adapter.evaluate(psi=plan['native_smoke_psi'], start_year=2007, horizon=6, seed=runner.seed,
                gates=plan['gates'], budget=plan['budget'], endpoint_controls=plan['endpoint'], path_controls=plan['path'],
                deadline=min(runner.deadline,time.monotonic()+plan['budget']['candidate_seconds']), folder=Path(output)/'six_date_smoke')
            runner.account(reply)
            receipt = dict(controller.smoke_readiness(reply), actual_policy_calls=runner.policy_calls,
                six_date_elapsed_seconds=time.monotonic()-started, identity=adapter.identity(),
                fresh_twelve_date_seed=True, maximum_maps=6)
            controller.write(Path(output)/'smoke_receipt.json', receipt)
        else:
            receipt = runner.run()
    except Exception as exc:
        controller.write(Path(output)/'failure.json', controller.failure_receipt(exc,runner))
        raise
    controller.write(Path(output)/'run_receipt.json', receipt)
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    for action in ('prepare-plan','preflight','native-smoke','execute'):
        actions.add_argument('--'+action, action='store_true')
    parser.add_argument('--plan', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--handoff')
    parser.add_argument('--mode', choices=('diagnostic','production'))
    parser.add_argument('--horizons', type=int, nargs=2)
    parser.add_argument('--execution-enabled', action='store_true')
    parser.add_argument('--fiscal-relaxation-authorized', action='store_true')
    for name in BUDGET_KEYS:
        parser.add_argument('--'+name.replace('_','-'), type=int if name == 'maximum_policy_calls' else float)
    for name in ('fit','path','endpoint'):
        parser.add_argument('--'+name+'-max-evaluations',type=int)
    args = parser.parse_args(argv)
    if args.prepare_plan:
        required = ('handoff','mode','horizons',*BUDGET_KEYS,'fit_max_evaluations','path_max_evaluations','endpoint_max_evaluations')
        missing = [key for key in required if getattr(args,key) is None]
        if missing: parser.error('Preparation requires explicit '+', '.join(missing))
        handoff = pin(args.handoff)
        validate_handoff(handoff)
        rt = runtime_module().CurrentEstateARuntime.from_handoff(handoff,Path(args.output)/'runtime')
        plan = build_plan(handoff,rt,mode=args.mode,horizons=args.horizons,
            budget={key:getattr(args,key) for key in BUDGET_KEYS}, fit_max_evaluations=args.fit_max_evaluations,
            path_max_evaluations=args.path_max_evaluations, endpoint_max_evaluations=args.endpoint_max_evaluations,
            fiscal_relaxation_authorized=args.fiscal_relaxation_authorized, execution_enabled=args.execution_enabled)
        controller.write(args.plan,plan)
        receipt = preflight(plan)
        controller.write(Path(args.output)/'preflight.json',receipt)
    else:
        plan = json.loads(Path(args.plan).read_text())
        if args.preflight:
            receipt = preflight(plan)
            controller.write(Path(args.output)/'preflight.json',receipt)
        else:
            receipt = execute(plan,args.output,native_smoke=args.native_smoke)
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
