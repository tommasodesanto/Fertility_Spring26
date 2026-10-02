"""Launch only terminal candidate-capped local chains, with live RAM gates."""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import psutil

from controller import HERE, RUNS, original_contract, read, write

ONE_GIB = 1024**3
MIN_FREE = 3*ONE_GIB
WORKER_BUDGET = int(3.8*ONE_GIB)
SWAP_LIMIT = ONE_GIB
HERE_RUNS = HERE/'runs'
FINAL_STATUSES = {'restart_selected_numerically_verified',
                  'restart_selected_postcheck_failed',
                  'restart_selected_unverified_deadline_elapsed',
                  'parent_selected_retained', 'no_verified_selection'}


def eligible():
    result = []
    for chain in range(48, 58):
        terminal = RUNS/f'chain{chain}/worker_terminal.json'
        if not terminal.is_file():
            continue  # Active original search is never interrupted.
        if (HERE_RUNS/f'chain{chain}').exists() or (HERE_RUNS/f'chain{chain}.log').exists():
            continue  # No restart may be retried.
        try:
            parent, launch, worker, search, selected, target = original_contract(chain)
        except (AssertionError, KeyError, FileNotFoundError):
            continue
        result.append(dict(chain=chain, original_deadline_epoch=worker['deadline_epoch'],
                           original_objective_calls=search['objective_calls'],
                           valid_selected=selected is not None))
    # Prefer chains with no valid point, then early original deadline.
    return sorted(result, key=lambda row: (row['valid_selected'], row['original_deadline_epoch']))


def snapshot(launched, swapout_base, swap_used_base):
    mem = psutil.virtual_memory()
    swap = psutil.swap_memory()
    workers = []
    for chain, record in launched.items():
        try:
            proc = psutil.Process(record['pid'])
            rss = proc.memory_info().rss + sum(child.memory_info().rss
                                               for child in proc.children(recursive=True)
                                               if child.is_running())
            state = proc.status()
        except psutil.Error:
            rss, state = 0, 'exited'
        summary = HERE_RUNS/f'chain{chain}/restart_summary.json'
        workers.append(dict(chain=chain, pid=record['pid'], state=state,
                            rss_bytes=rss, summary_status=read(summary)['status'] if summary.is_file() else None))
    return dict(time_epoch=time.time(), available_bytes=mem.available,
                swap_used_since_start_bytes=max(0, swap.used-swap_used_base),
                swapout_since_start_bytes=max(0, swap.sout-swapout_base),
                workers=workers)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    assert os.environ.get('ALLOW_LOCAL_CALIBRATION') == '1' or args.dry_run
    planned = eligible()
    assert planned, 'No terminal candidate-capped local chains with time remaining'
    assert not (HERE_RUNS/'launch_contract.json').exists(), 'Restart launcher already used'
    if args.dry_run:
        print(json.dumps(dict(status='filesystem_preflight_passed_zero_solves', planned=planned,
                              existing_results_refused=True, active_chains_untouched=True), indent=2))
        return
    HERE_RUNS.mkdir(parents=True, exist_ok=True)
    with (HERE_RUNS/'launch_contract.json').open('x') as file:
        json.dump(dict(status='active', planned=planned, started_epoch=time.time(),
                       no_retries=True, original_deadlines_preserved=True,
                       cumulative_250_call_cap_preserved=True, one_thread_each=True,
                       minimum_memory_headroom_bytes=MIN_FREE,
                       per_worker_budget_bytes=WORKER_BUDGET,
                       active_original_workers_untouched=True), file, indent=2)
        file.write('\n')
    original_swap = psutil.swap_memory()
    launched = {}
    stop_reason = None
    for row in planned:
        chain = row['chain']
        if launched:
            time.sleep(60)  # Let preceding worker reach its normal resident size.
        state = snapshot(launched, original_swap.sout, original_swap.used)
        if time.time()+900 >= row['original_deadline_epoch']:
            stop_reason = f'original_deadline_reserve_chain_{chain}'
            break
        if state['available_bytes'] < WORKER_BUDGET+MIN_FREE:
            stop_reason = f'memory_headroom_chain_{chain}'
            break
        if (state['swapout_since_start_bytes'] > SWAP_LIMIT or
                state['swap_used_since_start_bytes'] > SWAP_LIMIT):
            stop_reason = f'active_swap_growth_chain_{chain}'
            break
        out = HERE_RUNS/f'chain{chain}'
        assert not out.exists(), 'Existing scientific results are never overwritten'
        log = HERE_RUNS/f'chain{chain}.log'
        command = [sys.executable, str(HERE/'bootstrap.py'), '--chain', str(chain),
                   '--out', str(out)]
        with log.open('x') as handle:
            process = subprocess.Popen(command, cwd=HERE.parents[5],
                                       env=dict(os.environ, ALLOW_LOCAL_CALIBRATION='1'),
                                       stdin=subprocess.DEVNULL, stdout=handle,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        launched[chain] = dict(pid=process.pid, command=command,
                               original_deadline_epoch=row['original_deadline_epoch'],
                               started_epoch=time.time(), out=str(out), log=str(log))
        write(HERE_RUNS/'pids.json', launched)
        write(HERE_RUNS/'latest_progress.json', dict(status='launching',
              planned=len(planned), launched=len(launched), stop_reason=None,
              state=snapshot(launched, original_swap.sout, original_swap.used)))
    while launched:
        state = snapshot(launched, original_swap.sout, original_swap.used)
        write(HERE_RUNS/'latest_progress.json', dict(status='monitor', planned=len(planned),
              launched=len(launched), stop_reason=stop_reason, state=state))
        if all(worker['summary_status'] in FINAL_STATUSES or worker['state'] in ('exited', psutil.STATUS_ZOMBIE)
               for worker in state['workers']):
            break
        if time.time() >= read(RUNS/'launch_contract.json')['hard_deadline_epoch']:
            stop_reason = 'original_eight_hour_batch_cap'
            for record in launched.values():
                try:
                    os.killpg(record['pid'], signal.SIGTERM)
                except ProcessLookupError:
                    pass
            break
        time.sleep(60)
    write(HERE_RUNS/'terminal.json', dict(status='all_planned_launched' if len(launched)==len(planned)
          else 'partial_launch', stop_reason=stop_reason, planned=len(planned),
          launched=len(launched), state=snapshot(launched, original_swap.sout, original_swap.used),
          finished_epoch=time.time(), no_retries=True))


if __name__ == '__main__':
    main()
