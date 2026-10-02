"""Continue a stopped local chain within its original clock and call budget."""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import time
from pathlib import Path

LOCAL = Path(__file__).resolve().parents[1]
PACKET = LOCAL.parent
RUNS = LOCAL/'runs/local10_v1'
HERE = Path(__file__).resolve().parent
REVIEWED_OPTIMIZER_SHA256 = '9f72add8cced622aae2532ca231231a27386d22f601cdcf3f2e6f2fa3ac96c1b'
sys.path.insert(0, str(PACKET))
sys.path.insert(0, str(LOCAL))
from restart_controller_v2.controller import install_censored_optimizer  # noqa: E402


def read(path: Path):
    return json.loads(path.read_text())


def write(path: Path, payload: dict):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
    temp.replace(path)


def original_contract(chain: int):
    assert hashlib.sha256(inspect.getsource(install_censored_optimizer).encode()).hexdigest() == REVIEWED_OPTIMIZER_SHA256
    parent = RUNS/f'chain{chain}'
    terminal = read(parent/'worker_terminal.json')
    search = read(parent/'search/completed.json')
    cases = read(parent/'search/cases.json')
    registry = read(RUNS/('pids.json' if chain in (48, 49, 50, 53, 54, 55) else 'resume_pids.json'))
    launch = registry[str(chain)]
    assert terminal['chain'] == chain
    assert terminal['deadline_epoch'] == launch['deadline_epoch'] == search['deadline_epoch']
    assert abs(terminal['deadline_epoch']-launch['started_epoch']-14400) < 2
    assert terminal['search_exit_code'] == 0
    assert search['search_stop_reason'] == 'native_evaluation_budget_exhausted'
    assert 0 < search['objective_calls'] < 250
    assert len(cases) == search['completed_full_ge'] <= search['objective_calls']
    assert cases[-1]['status'] == 'budget_exhausted'
    assert cases[-1]['reason'] == 'uncomputed_bounded_budget'
    assert cases[-1]['lifecycle_solves'] >= 31
    assert time.time()+900 < terminal['deadline_epoch']
    selected = search.get('selected')
    if selected is None:
        assert terminal['status'] == 'search_no_selected_candidate'
    else:
        assert terminal['status'] == 'selected_numerically_verified'
        assert selected['status'] == 'passed' and selected['label'] in {case['label'] for case in cases}
        postcheck = read(parent/'postcheck/completed.json')
        assert postcheck['status'] == 'selected_numerically_verified'
        assert postcheck['selected']['label'] == selected['label']
        assert postcheck['selected']['loss'] == selected['loss']
    target = read(parent/'search/input_contract.json')['parameter_contract']['target_contract_sha256']
    assert target == read(LOCAL/'local_plan.json')['target_fingerprint'] or target == read(PACKET/'plan.json')['target_fingerprint']
    return parent, launch, terminal, search, selected, target


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chain', type=int, choices=range(48, 58), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    parent, launch, terminal, old, selected, target = original_contract(args.chain)
    out = args.out.resolve()
    assert out == (HERE/'runs'/f'chain{args.chain}').resolve()
    remaining = 250-old['objective_calls']
    contract = dict(status='prepared', chain=args.chain, parent=str(parent),
                    parent_search_receipt_sha256=hashlib.sha256((parent/'search/completed.json').read_bytes()).hexdigest(),
                    parent_worker_receipt_sha256=hashlib.sha256((parent/'worker_terminal.json').read_bytes()).hexdigest(),
                    parent_postcheck_receipt_sha256=hashlib.sha256((parent/'postcheck/completed.json').read_bytes()).hexdigest() if selected else None,
                    original_start_epoch=launch['started_epoch'], original_deadline_epoch=terminal['deadline_epoch'],
                    original_objective_calls=old['objective_calls'], remaining_objective_calls=remaining,
                    original_selected_loss=selected['loss'] if selected else None,
                    original_target_contract_sha256=target,
                    reviewed_optimizer_function_sha256=REVIEWED_OPTIMIZER_SHA256,
                    no_clock_reset=True, no_call_count_reset=True, no_native_cap_change=True)
    if args.dry_run:
        print(json.dumps(contract, indent=2)); return
    assert os.environ.get('ALLOW_LOCAL_CALIBRATION') == '1'
    assert not out.exists(), 'Refusing existing restart results'
    out.mkdir(parents=True)
    write(out/'restart_contract.json', contract)
    sys.argv = ['run_local_psi.py', '--chain', str(args.chain), '--out', str(out/'search'),
                '--deadline-epoch', str(terminal['deadline_epoch']), '--fast-objective']
    import run_local_psi as original
    if selected is not None:
        original.CONFIG['nearby_starts'][args.chain]['parameters'] = selected['parameters']
    contract['optimizer_source_sha256'] = install_censored_optimizer(original, remaining)
    write(out/'restart_contract.json', contract)
    original.main()
    new = read(out/'search/completed.json')
    assert old['objective_calls']+new['objective_calls'] <= 250
    assert new['selected'] is None or new['selected']['status'] == 'passed'
    winner = 'restart' if new['selected'] is not None and (selected is None or new['selected']['loss'] < selected['loss']) else 'parent'
    receipt = dict(status='restart_search_finished', chain=args.chain, winner=winner,
                   cumulative_objective_calls=old['objective_calls']+new['objective_calls'],
                   original_deadline_epoch=terminal['deadline_epoch'],
                   parent_selected_loss=selected['loss'] if selected else None,
                   restart_selected_loss=new['selected']['loss'] if new['selected'] else None,
                   original_search_stop=old['search_stop_reason'],
                   restart_search_stop=new['search_stop_reason'],
                   selected_requires_fresh_postcheck=winner == 'restart')
    write(out/'restart_summary.json', receipt)
    if winner == 'restart':
        if time.time() >= terminal['deadline_epoch']:
            receipt['status'] = 'restart_selected_unverified_deadline_elapsed'
        else:
            command = [sys.executable, str(LOCAL/'bootstrap.py'), '--chain', str(args.chain),
                       '--out', str(out/'postcheck'), '--deadline-epoch', str(terminal['deadline_epoch']),
                       '--verify-only', str(out/'search/completed.json')]
            with (out/'postcheck.log').open('x') as log:
                code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
            check = out/'postcheck/completed.json'
            receipt['postcheck_exit_code'] = code
            receipt['status'] = ('restart_selected_numerically_verified'
                                 if code == 0 and check.is_file() and read(check)['status'] == 'selected_numerically_verified'
                                 else 'restart_selected_postcheck_failed')
    else:
        receipt['status'] = 'parent_selected_retained' if selected else 'no_verified_selection'
    write(out/'restart_summary.json', receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
