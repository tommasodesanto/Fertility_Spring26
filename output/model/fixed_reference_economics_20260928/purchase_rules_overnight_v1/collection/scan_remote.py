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
if SOURCE_KIND not in ('original', 'restart', 'local_restart'):
    raise ValueError('Unknown collection source kind')
PARENT_RESULTS = Path(os.environ.get('PURCHASE_PARENT_RESULTS_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/results'))
PACKET = Path(os.environ.get('PURCHASE_PACKET_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1'))
PLAN = json.loads((PACKET / 'plan.json').read_text())
CONTRACT = PLAN['base_target_contract']


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
    expected_status = ('restart_selected_numerically_verified' if winner == 'restart' else
                       'parent_selected_retained' if original.get('selected') else 'no_verified_selection')
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
            winner not in ('restart', 'parent') or summary['status'] != expected_status):
        raise ValueError('Local restart/parent provenance, final status, or budget drift')
    if winner == 'restart' and (continuation.get('selected') is None or
            original.get('selected') is not None and
            continuation['selected']['loss'] >= original['selected']['loss']):
        raise ValueError('Local restart is not a strict improvement')
    return dict(parent_remote_root=str(parent),restart_contract_sha256=sha(contract_path),
                restart_summary_sha256=sha(summary_path),restart_winner=winner,
                optimizer_source_sha256=contract['optimizer_source_sha256'],
                cumulative_objective_calls=summary['cumulative_objective_calls'],
                original_deadline_epoch=summary['original_deadline_epoch'])


def check_chain(index):
    arm = 'hard' if index < 24 or 48 <= index <= 52 else 'quarter'
    root = chain_root(index)
    provenance = restart_provenance(root, index) if SOURCE_KIND in ('restart', 'local_restart') else {}
    if SOURCE_KIND in ('restart', 'local_restart') and provenance is None:
        return dict(chain=index, arm=arm, status='restart_pending', remote_root=str(root))
    if SOURCE_KIND in ('restart', 'local_restart') and provenance['restart_winner'] == 'parent':
        return dict(chain=index, arm=arm, status='restart_parent_retained', remote_root=str(root), **provenance)
    completed = root / 'postcheck/completed.json'
    if not completed.exists():
        search_receipt = root / ('search/completed.json' if SOURCE_KIND == 'local_restart'
                                 else 'search/search_completed.json')
        return dict(chain=index, arm=arm, status='postcheck_pending',
                    search_completed=search_receipt.exists(), remote_root=str(root), **provenance)
    report = json.loads(completed.read_text())
    if report.get('status') != 'selected_numerically_verified':
        return dict(chain=index, arm=arm, status=report.get('status', 'invalid_postcheck'), error='Selected postcheck did not pass')
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
    return dict(chain=index, arm=arm, status='postchecked', loss=loss,
                source_run=SOURCE_KIND, **provenance,
                price=float(checked['price']), H0=float(checked['H0_derived']),
                target_fingerprint=PLAN['target_fingerprint'], weight_fingerprint=PLAN['weight_fingerprint'],
                selected_parameters=selected['parameters'], remote_root=str(root),
                remote_report=str(full), remote_arrays=str(arrays), target_fit=target,
                parameters=parameter, closure=closure,
                report_sha256={name: sha(full / name) for name in ('target_fit.csv', 'parameters.csv', 'closure.json')},
                native_arrays_bytes=arrays.stat().st_size)


if __name__ == '__main__':
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
