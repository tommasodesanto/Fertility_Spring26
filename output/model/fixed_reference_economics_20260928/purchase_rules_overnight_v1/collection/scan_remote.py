"""Read only completed, fresh-postchecked Torch calibration chains.

This standard-library script is sent to Torch on stdin by collect.py. It does
not solve a model or alter remote outputs.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path

REMOTE = Path(os.environ.get('PURCHASE_RESULTS_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/results'))
SOURCE_KIND = os.environ.get('PURCHASE_SOURCE_KIND', 'original')
if SOURCE_KIND not in ('original', 'restart', 'local_restart', 'regions'):
    raise ValueError('Unknown collection source kind')
PARENT_RESULTS = Path(os.environ.get('PURCHASE_PARENT_RESULTS_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/results'))
PACKET = Path(os.environ.get('PURCHASE_PACKET_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1'))
PLAN = json.loads((PACKET / 'plan.json').read_text())
CONTRACT = PLAN['base_target_contract']
REGION_ROOT = Path(os.environ.get('PURCHASE_REGION_ROOT', '/scratch/td2248/projects/purchase_broader_regions_v1'))


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def chain_root(index):
    """Local original chains may use chain48; Torch and restarts use chain_48."""
    if SOURCE_KIND == 'regions':
        return REMOTE / f'slot_{index}'
    names = ((f'chain{index}',) if SOURCE_KIND == 'local_restart' else
             (f'chain_{index}',) if SOURCE_KIND == 'restart' else
             (f'chain_{index}', f'chain{index}'))
    existing = [REMOTE / name for name in names if (REMOTE / name).is_dir()]
    if len(existing) > 1:
        raise ValueError('Ambiguous duplicate chain directory')
    return existing[0] if existing else REMOTE / names[0]


def parent_chain_root(index):
    names = (f'chain_{index}', f'chain{index}')
    existing = [PARENT_RESULTS / name for name in names if (PARENT_RESULTS / name).is_dir()]
    if len(existing) != 1:
        raise ValueError('Missing or ambiguous original parent chain')
    return existing[0]


def restart_provenance(root, index):
    """Bind a continuation to its unchanged parent and the original budgets."""
    if SOURCE_KIND == 'local_restart':
        return local_restart_provenance(root, index)
    contract_path = root / 'restart_contract.json'
    summary_path = root / 'restart_summary.json'
    if not contract_path.exists() or not summary_path.exists():
        if (root / 'postcheck/completed.json').exists():
            raise ValueError('Unbound restart postcheck without contract and summary')
        return None
    parent = parent_chain_root(index)
    contract, summary = json.loads(contract_path.read_text()), json.loads(summary_path.read_text())
    launch = json.loads((parent / 'launcher_start.json').read_text())
    original = json.loads((parent / 'search/search_completed.json').read_text())
    continuation = json.loads((root / 'search/search_completed.json').read_text())
    if (contract['chain'] != index or summary['chain'] != index or
            contract['parent'] not in (str(parent), '/work/parent') or
            contract['parent_search_receipt_sha256'] != sha(parent / 'search/search_completed.json') or
            contract['parent_launcher_receipt_sha256'] != sha(parent / 'launcher_start.json') or
            contract['parent_postcheck_receipt_sha256'] != (
                sha(parent / 'postcheck/completed.json') if original.get('selected') else None) or
            contract['original_start_epoch'] != launch['start_epoch'] or
            contract['original_deadline_epoch'] != launch['deadline_epoch'] or
            launch['deadline_epoch'] != launch['start_epoch'] + 14400 or
            contract['original_objective_calls'] != original['objective_calls'] or
            contract['remaining_objective_calls'] != 250 - original['objective_calls'] or
            contract.get('no_clock_reset') is not True or contract.get('no_call_count_reset') is not True or
            re.fullmatch(r'[0-9a-f]{64}', contract.get('optimizer_source_sha256', '')) is None or
            contract['original_target_contract_sha256'] != canonical(CONTRACT) or
            original.get('search_stop_reason') != 'native_evaluation_budget_exhausted' or
            summary['original_deadline_epoch'] != launch['deadline_epoch'] or
            summary['cumulative_objective_calls'] != original['objective_calls'] + continuation['objective_calls'] or
            summary['cumulative_objective_calls'] > 250 or
            summary['parent_selected_loss'] != (original['selected']['loss'] if original.get('selected') else None) or
            summary['restart_selected_loss'] != (continuation['selected']['loss'] if continuation.get('selected') else None) or
            summary['selected_requires_fresh_postcheck'] != (summary['winner'] == 'restart') or
            summary['winner'] not in ('restart', 'parent')):
        raise ValueError('Restart/parent provenance or objective budget drift')
    if summary['winner'] == 'restart':
        if (continuation.get('selected') is None or
                original.get('selected') is not None and
                continuation['selected']['loss'] >= original['selected']['loss']):
            raise ValueError('Restart is not a strict improvement over parent')
    return dict(parent_remote_root=str(parent),restart_contract_sha256=sha(contract_path),
                restart_summary_sha256=sha(summary_path),restart_winner=summary['winner'],
                optimizer_source_sha256=contract['optimizer_source_sha256'],
                cumulative_objective_calls=summary['cumulative_objective_calls'],
                original_deadline_epoch=summary['original_deadline_epoch'])


def local_restart_provenance(root, index):
    """Local continuations have worker-terminal and PID-registry receipts."""
    contract_path, summary_path = root / 'restart_contract.json', root / 'restart_summary.json'
    if not contract_path.exists() or not summary_path.exists():
        if (root / 'postcheck/completed.json').exists():
            raise ValueError('Unbound local restart postcheck')
        return None
    parent = parent_chain_root(index)
    contract, summary = json.loads(contract_path.read_text()), json.loads(summary_path.read_text())
    terminal = json.loads((parent / 'worker_terminal.json').read_text())
    original = json.loads((parent / 'search/completed.json').read_text())
    continuation = json.loads((root / 'search/completed.json').read_text())
    registry_name = 'pids.json' if index in (48, 49, 50, 53, 54, 55) else 'resume_pids.json'
    registry = json.loads((PARENT_RESULTS / registry_name).read_text())
    launch = registry[str(index)]
    winner = summary.get('winner')
    if (contract['chain'] != index or summary['chain'] != index or
            contract['parent'] != str(parent) or terminal['chain'] != index or
            contract['parent_search_receipt_sha256'] != sha(parent / 'search/completed.json') or
            contract['parent_worker_receipt_sha256'] != sha(parent / 'worker_terminal.json') or
            contract['parent_postcheck_receipt_sha256'] != (
                sha(parent / 'postcheck/completed.json') if original.get('selected') else None) or
            contract['original_start_epoch'] != launch['started_epoch'] or
            contract['original_deadline_epoch'] != terminal['deadline_epoch'] or
            terminal['deadline_epoch'] != launch['deadline_epoch'] or
            abs(terminal['deadline_epoch'] - launch['started_epoch'] - 14400) >= 2 or
            contract['original_objective_calls'] != original['objective_calls'] or
            contract['remaining_objective_calls'] != 250 - original['objective_calls'] or
            contract.get('no_clock_reset') is not True or contract.get('no_call_count_reset') is not True or
            contract.get('no_native_cap_change') is not True or
            contract.get('reviewed_optimizer_function_sha256') !=
                '9f72add8cced622aae2532ca231231a27386d22f601cdcf3f2e6f2fa3ac96c1b' or
            re.fullmatch(r'[0-9a-f]{64}', contract.get('optimizer_source_sha256', '')) is None or
            contract['original_target_contract_sha256'] != canonical(CONTRACT) or
            original.get('search_stop_reason') != 'native_evaluation_budget_exhausted' or
            summary['original_deadline_epoch'] != terminal['deadline_epoch'] or
            summary['cumulative_objective_calls'] != original['objective_calls'] + continuation['objective_calls'] or
            summary['cumulative_objective_calls'] > 250 or
            summary['parent_selected_loss'] != (original['selected']['loss'] if original.get('selected') else None) or
            summary['restart_selected_loss'] != (continuation['selected']['loss'] if continuation.get('selected') else None) or
            summary['selected_requires_fresh_postcheck'] != (winner == 'restart') or
            winner not in ('restart', 'parent')):
        raise ValueError('Local restart/parent provenance or budget drift')
    if winner == 'restart' and (continuation.get('selected') is None or
            original.get('selected') is not None and
            continuation['selected']['loss'] >= original['selected']['loss']):
        raise ValueError('Local restart is not a strict improvement')
    status = summary['status']
    if status == 'restart_search_finished':
        state = 'postcheck_running' if winner == 'restart' else 'final_receipt_pending'
    elif winner == 'restart' and status == 'restart_selected_numerically_verified':
        if summary.get('postcheck_exit_code') != 0 or not (root / 'postcheck/completed.json').is_file():
            raise ValueError('Local restart claims verification without successful native postcheck')
        state = 'verified'
    elif winner == 'restart' and status in ('restart_selected_unverified_deadline_elapsed',
                                             'restart_selected_postcheck_failed'):
        state = status
    elif winner == 'parent' and status == ('parent_selected_retained' if original.get('selected')
                                           else 'no_verified_selection'):
        state = 'parent_retained'
    else:
        raise ValueError('Local restart status inconsistent with selected winner')
    return dict(parent_remote_root=str(parent),restart_contract_sha256=sha(contract_path),
                restart_summary_sha256=sha(summary_path),restart_winner=winner,
                restart_state=state,
                optimizer_source_sha256=contract['optimizer_source_sha256'],
                cumulative_objective_calls=summary['cumulative_objective_calls'],
                original_deadline_epoch=summary['original_deadline_epoch'])


def region_design():
    """Pin the approved initial-vector-only array to its launched source and design."""
    design_path = REGION_ROOT / 'source/design.json'
    pins_path = REGION_ROOT / 'source_sha256.json'
    expected_design = os.environ.get('PURCHASE_REGION_DESIGN_SHA256')
    expected_pins = os.environ.get('PURCHASE_REGION_PINS_SHA256')
    if not expected_design or not expected_pins or sha(design_path) != expected_design or sha(pins_path) != expected_pins:
        raise ValueError('Approved broad-region design or source-pin identity drift')
    pins = json.loads(pins_path.read_text())
    for name, digest in pins.items():
        if sha(REGION_ROOT / 'source' / name) != digest:
            raise ValueError('Broad-region source hash drift: ' + name)
    for name in ('launch_torch.sh', 'preflight_torch.sh'):
        if sha(REGION_ROOT / name) != pins[name]:
            raise ValueError('Launched broad-region script hash drift: ' + name)
    design = json.loads(design_path.read_text())
    submit = json.loads((REGION_ROOT / 'submission_receipt.json').read_text())
    if (design['status'] != 'experimental_initialization_only' or
            design['economic_changes'] != ['initial optimizer vectors only; no model, target, weight or bound change'] or
            design['original_plan_sha256'] != sha(PACKET / 'plan.json') or
            design['original_target_contract_sha256'] != PLAN['target_fingerprint'] or
            design['original_weight_fingerprint'] != PLAN['weight_fingerprint'] or
            design['full_approved_bounds'] != PLAN['bounds'] or
            design['free_coordinates'] != PLAN['free_coordinates'] or
            design['scored_targets'] != 10 or design['total_targets'] != 14 or
            design['maximum_objective_calls_per_chain'] != 80 or
            design['wall_seconds_per_chain'] != 4500 or
            design['final_reserve_seconds'] != 900 or
            design['absolute_deadline_epoch'] != 1790932500 or
            len(design['starts']) != 16 or
            submit['status'] != 'submitted_reviewed_broad_starts' or
            str(submit['job_id']) != '19019947' or
            submit['slot_ids'] != list(range(16)) or
            submit['source_sha256'] != pins or
            submit['design_sha256'] != expected_design or
            submit['original_plan_sha256'] != design['original_plan_sha256'] or
            submit['target_contract_sha256'] != PLAN['target_fingerprint'] or
            submit['weight_fingerprint'] != PLAN['weight_fingerprint'] or
            submit['maximum_objective_calls_per_chain'] != 80 or
            submit['wall_seconds_per_chain'] != 4500 or
            submit['final_reserve_seconds'] != 900 or
            submit['absolute_deadline_epoch'] != 1790932500):
        raise ValueError('Broad-region design, target, budget or submission drift')
    for slot, start in enumerate(design['starts']):
        chain = slot if slot < 8 else slot + 16
        arm = 'hard' if slot < 8 else 'quarter'
        if (start['slot'] != slot or start['original_chain_index'] != chain or
                start['arm'] != arm or start['weight_profile'] != 'base_control' or
                set(start['parameters']) != set(PLAN['free_coordinates']) or
                submit['mapping'][str(slot)] != dict(arm=arm, original_chain_index=chain)):
            raise ValueError('Broad-region slot/original-chain mapping drift')
        for name, value in start['parameters'].items():
            lo, hi = PLAN['bounds'][name]
            if not lo <= value <= hi:
                raise ValueError('Broad-region start outside approved bounds')
    return design, pins


def region_provenance(root, slot, design):
    start = design['starts'][slot]
    chain = start['original_chain_index']
    search_path = root / 'search_region_contract.json'
    if not search_path.is_file():
        if (root / 'postcheck/completed.json').exists():
            raise ValueError('Unbound broad-region postcheck')
        return None
    search = json.loads(search_path.read_text())
    for stage, receipt in [('search', search)]:
        if (receipt['stage'] != stage or receipt['slot'] != slot or
                receipt['original_chain_index'] != chain or receipt['arm'] != start['arm'] or
                receipt['parameters'] != start['parameters'] or
                receipt['design_sha256'] != sha(REGION_ROOT / 'source/design.json') or
                receipt['runner_sha256'] != sha(PACKET / 'run_psi.py') or
                receipt['reviewed_censor_sha256'] != design['reviewed_censor_source_sha256'] or
                receipt['original_plan_sha256'] != design['original_plan_sha256'] or
                receipt['target_contract_sha256'] != PLAN['target_fingerprint'] or
                receipt['weight_fingerprint'] != PLAN['weight_fingerprint'] or
                receipt['max_objective_calls'] != 80 or receipt['final_reserve_seconds'] != 900 or
                receipt['absolute_deadline_epoch'] != 1790932500 or
                receipt.get('experimental_initialization_only') is not True or
                receipt['deadline_epoch'] > 1790932500 or
                re.fullmatch(r'[0-9a-f]{64}', receipt.get('optimizer_source_sha256', '')) is None):
            raise ValueError('Broad-region search contract drift')
    if sha(Path(os.environ.get('PURCHASE_REVIEWED_CENSOR', '/scratch/td2248/projects/purchase_restart_controller_v2/source/controller.py'))) != design['reviewed_censor_source_sha256']:
        raise ValueError('Reviewed candidate-censor source drift')
    native_search = root / 'search/search_contract.json'
    if native_search.is_file():
        native = json.loads(native_search.read_text())
        if (native['bounds'] != PLAN['bounds'] or native['free_coordinates'] != PLAN['free_coordinates'] or
                native['seed'] != start['parameters'] or
                native['maximum_objective_calls'] != 80 or native['final_reserve_seconds'] != 900 or
                native['targets_scored'] != 10 or native['targets_total'] != 14 or
                native['exact_vector_cache'] is not True):
            raise ValueError('Broad-region native search contract drift')
    input_path = root / 'search/input_contract.json'
    if input_path.is_file():
        inp = json.loads(input_path.read_text())
        if (inp['bounds'] != PLAN['bounds'] or inp['purchase_rule'] != start['arm'] or
                inp['owner_financed_share'] != .8 or inp['normalized_population'] != 1. or
                inp['weight_contract_sha256'] != PLAN['weight_fingerprint'] or
                inp['base_target_contract_sha256'] != PLAN['target_fingerprint']):
            raise ValueError('Broad-region native search target or purchase contract drift')
    post_path = root / 'postcheck_region_contract.json'
    if post_path.is_file():
        post = json.loads(post_path.read_text())
        for key in ('slot','arm','original_chain_index','design_sha256','runner_sha256',
                    'reviewed_censor_sha256','deadline_epoch','absolute_deadline_epoch',
                    'max_objective_calls','final_reserve_seconds','parameters',
                    'original_plan_sha256','target_contract_sha256','weight_fingerprint'):
            if post.get(key) != search.get(key):
                raise ValueError('Broad-region search/postcheck contract drift: ' + key)
        if post.get('stage') != 'postcheck' or post.get('optimizer_source_sha256') is not None:
            raise ValueError('Broad-region postcheck stage drift')
    elif (root / 'postcheck/completed.json').exists():
        raise ValueError('Unbound broad-region selected postcheck')
    terminal_path = root / 'launcher_terminal.json'
    if terminal_path.is_file():
        terminal = json.loads(terminal_path.read_text())
        if (terminal['slot'] != slot or terminal['deadline_epoch'] != search['deadline_epoch'] or
                terminal['absolute_deadline_epoch'] != 1790932500 or
                terminal['maximum_objective_calls'] != 80 or terminal['final_reserve_seconds'] != 900 or
                terminal['deadline_epoch'] > terminal['start_epoch'] + 4500 or
                terminal['deadline_epoch'] > 1790932500):
            raise ValueError('Broad-region launcher clock or call budget drift')
    else:
        terminal = None
    receipt_path = root / 'search/search_completed.json'
    if receipt_path.is_file():
        result = json.loads(receipt_path.read_text())
        if result['objective_calls'] > 80:
            raise ValueError('Broad-region search exceeded approved objective budget')
    claimed = (root / 'postcheck/completed.json').is_file()
    if claimed and (not post_path.is_file() or not receipt_path.is_file() or
                    not native_search.is_file() or not input_path.is_file() or
                    terminal is not None and terminal['exit_code'] != 0):
        raise ValueError('Claimed broad-region postcheck lacks successful source/launcher/search receipts')
    return dict(region_slot=slot, original_chain_index=chain,
                region_state='awaiting_terminal' if claimed and terminal is None else 'ready',
                region_design_sha256=search['design_sha256'],
                region_search_contract_sha256=sha(search_path),
                region_postcheck_contract_sha256=sha(post_path) if post_path.is_file() else None,
                region_launcher_terminal_sha256=sha(terminal_path) if terminal_path.is_file() else None)


def check_chain(index):
    arm = ('hard' if index < 8 else 'quarter') if SOURCE_KIND == 'regions' else ('hard' if index < 24 or 48 <= index <= 52 else 'quarter')
    root = chain_root(index)
    if SOURCE_KIND == 'regions':
        provenance = region_provenance(root, index, REGION_DESIGN)
        chain = REGION_DESIGN['starts'][index]['original_chain_index']
        if provenance is None:
            return dict(chain=chain, arm=arm, status='region_pending', region_slot=index, remote_root=str(root))
        if provenance['region_state'] == 'awaiting_terminal':
            return dict(chain=chain, arm=arm, status='postcheck_terminal_pending',
                        remote_root=str(root), **provenance)
    else:
        provenance = restart_provenance(root, index) if SOURCE_KIND in ('restart', 'local_restart') else {}
        chain = index
    if SOURCE_KIND in ('restart', 'local_restart') and provenance is None:
        return dict(chain=index, arm=arm, status='restart_pending', remote_root=str(root))
    if SOURCE_KIND == 'local_restart' and provenance['restart_state'] != 'verified':
        return dict(chain=index, arm=arm, status=provenance['restart_state'],
                    remote_root=str(root), **provenance)
    if SOURCE_KIND in ('restart', 'local_restart') and provenance['restart_winner'] == 'parent':
        return dict(chain=index, arm=arm, status='restart_parent_retained', remote_root=str(root), **provenance)
    completed = root / 'postcheck/completed.json'
    if not completed.exists():
        search_receipt = root / ('search/completed.json' if SOURCE_KIND == 'local_restart'
                                 else 'search/search_completed.json')
        return dict(chain=chain, arm=arm, status='postcheck_pending',
                    search_completed=search_receipt.exists(), remote_root=str(root), **provenance)
    report = json.loads(completed.read_text())
    if report.get('status') != 'selected_numerically_verified':
        return dict(chain=chain, arm=arm, status=report.get('status', 'invalid_postcheck'), error='Selected postcheck did not pass')
    selected = report['selected']
    checked = report['selected_postcheck']
    search = json.loads((root / ('search/completed.json' if SOURCE_KIND == 'local_restart'
                                 else 'search/search_completed.json')).read_text())
    if (search.get('selected') is None or selected['parameters'] != search['selected']['parameters'] or
            selected['loss'] != search['selected']['loss'] or
            selected.get('weight_fingerprint') != PLAN['weight_fingerprint']):
        raise ValueError('Fresh selected postcheck/search coordinates or weights drift')
    if checked['status'] != 'passed' or checked['weight_fingerprint'] != PLAN['weight_fingerprint']:
        raise ValueError('postcheck status or weight fingerprint drift')
    input_contract = json.loads((root / 'postcheck/input_contract.json').read_text())
    if (input_contract['weight_contract_sha256'] != PLAN['weight_fingerprint'] or
            input_contract['base_target_contract_sha256'] != PLAN['target_fingerprint'] or
            input_contract['purchase_rule'] != arm or
            input_contract['owner_financed_share'] != .8 or
            input_contract['normalized_population'] != 1.):
        raise ValueError('Selected postcheck target, weight, or purchase contract drift')
    if canonical(CONTRACT) != PLAN['target_fingerprint']:
        raise ValueError('target contract fingerprint drift')
    base = root / 'postcheck/selected_postcheck/phase_b_ge'
    full = base / 'selected_root'
    target = rows(full / 'target_fit.csv')
    parameter = rows(full / 'parameters.csv')
    if len(target) != 14 or len(parameter) != 31:
        raise ValueError('target or parameter row count drift')
    identity = [{k: x[k] for k in ('moment', 'target', 'weight', 'role')} for x in target]
    if identity != CONTRACT or target != checked['base_target_fit']:
        raise ValueError('target table identity or postcheck mismatch')
    if len([x for x in target if x['role'] == 'scored']) != 10:
        raise ValueError('ten scored targets required')
    loss = sum(float(x['loss_contribution']) for x in target if x['loss_contribution'])
    if not math.isclose(loss, float(checked['weighted_loss']), rel_tol=0, abs_tol=1e-8):
        raise ValueError('loss contribution mismatch')
    if not math.isclose(loss, float(selected['loss']), rel_tol=0, abs_tol=1e-8):
        raise ValueError('search/postcheck loss mismatch')
    bound = PLAN['bounds']
    for row in parameter:
        key = row['parameter']
        if key in bound:
            if not (math.isclose(float(row['lower']), float(bound[key][0]), abs_tol=1e-12) and
                    math.isclose(float(row['upper']), float(bound[key][1]), abs_tol=1e-12)):
                raise ValueError('free parameter bound drift: ' + key)
    if len([x for x in parameter if x['parameter'] in bound]) != 10:
        raise ValueError('ten free coordinates required')
    closure = json.loads((full / 'closure.json').read_text())
    if closure['population_scale'] != 1. or closure['normalized_population'] != 1.:
        raise ValueError('population normalization drift')
    if closure['standard_plot_count'] != 17 or len(list((full / 'standard_diagnostics').glob('*.png'))) != 17:
        raise ValueError('standard diagnostic set incomplete')
    if abs(float(closure['renewal_residual'])) > 1e-6 or abs(float(closure['absolute_housing_residual'])) > 1e-6:
        raise ValueError('renewal or housing closure failed')
    arrays = base / 'selected_repeat/stage/solution_arrays.npz'
    if not arrays.is_file():
        raise ValueError('selected native checkpoint missing')
    # NPZ members are checked by the buyer and policy workers when loaded.
    return dict(chain=chain, arm=arm, status='postchecked', loss=loss,
                source_run=SOURCE_KIND, **provenance,
                price=float(checked['price']), H0=float(checked['H0_derived']),
                target_fingerprint=PLAN['target_fingerprint'], weight_fingerprint=PLAN['weight_fingerprint'],
                selected_parameters=selected['parameters'], remote_root=str(root),
                remote_report=str(full), remote_arrays=str(arrays), target_fit=target,
                parameters=parameter, closure=closure,
                report_sha256={name: sha(full / name) for name in ('target_fit.csv', 'parameters.csv', 'closure.json')},
                native_arrays_bytes=arrays.stat().st_size)


if __name__ == '__main__':
    if SOURCE_KIND == 'regions':
        REGION_DESIGN, _ = region_design()
    result = {'target_fingerprint': PLAN['target_fingerprint'],
              'weight_fingerprint': PLAN['weight_fingerprint'], 'chains': [], 'errors': []}
    first = int(os.environ.get('PURCHASE_CHAIN_FIRST', '0'))
    last = int(os.environ.get('PURCHASE_CHAIN_LAST', '48'))
    for index in range(first, last):
        try:
            result['chains'].append(check_chain(index))
        except Exception as error:
            result['errors'].append(dict(chain=index, error=f'{type(error).__name__}: {error}'))
    if result['errors']:
        result['status'] = 'rejected_invalid_postcheck'
    else:
        result['status'] = 'valid_snapshot'
    print(json.dumps(result, sort_keys=True, allow_nan=False))
