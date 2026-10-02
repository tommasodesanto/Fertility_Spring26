"""Append the four untouched local starts after a numerical rejection in chain 50.

This never restarts a chain and never changes the pinned objective or its budget.
Run once in a persistent terminal after lead review.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import time
from pathlib import Path

import psutil

from supervisor import (CHAIN_SECONDS, HERE, INITIAL_NEW_WORKER_BUDGET,
                        MIN_FREE, PACKET, SWAPOUT_LIMIT, write)

ROOT = HERE/'runs/local10_v1'
PAIRS = ((51, 56), (52, 57))


def read(path: Path):
    return json.loads(path.read_text())


def snapshot(new_workers, swapout_base, swap_used_base):
    mem = psutil.virtual_memory()
    swap = psutil.swap_memory()
    rows = []
    for chain, info in new_workers.items():
        try:
            process = psutil.Process(info['pid'])
            rss = process.memory_info().rss + sum(
                child.memory_info().rss for child in process.children(recursive=True)
                if child.is_running())
            state = process.status()
        except psutil.Error:
            rss, state = 0, 'exited'
        terminal = ROOT/f'chain{chain}/worker_terminal.json'
        rows.append(dict(chain=chain, pid=info['pid'], rss_bytes=rss, state=state,
                         terminal_status=read(terminal)['status'] if terminal.is_file() else None))
    return dict(time_epoch=time.time(), available_bytes=mem.available,
                swap_used_since_resume_bytes=max(0, swap.used-swap_used_base),
                swapout_since_resume_bytes=max(0, swap.sout-swapout_base),
                workers=rows)


def launch_pair(pair, new_workers, hard_deadline):
    for chain in pair:
        out = ROOT/f'chain{chain}'
        assert not out.exists(), f'Existing chain is never retried: {out}'
        log = ROOT/f'chain{chain}.log'
        deadline = min(time.time()+CHAIN_SECONDS, hard_deadline)
        command = [str(HERE/'launch_one.sh'), '--chain', str(chain), '--out',
                   str(out), '--deadline-epoch', str(deadline)]
        with log.open('x') as handle:
            process = subprocess.Popen(command, cwd=PACKET.parents[3],
                                       env=dict(os.environ, ALLOW_LOCAL_CALIBRATION='1'),
                                       stdin=subprocess.DEVNULL, stdout=handle,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        new_workers[chain] = dict(pid=process.pid, out=str(out), log=str(log),
                                  command=command, started_epoch=time.time(),
                                  deadline_epoch=deadline)
        write(ROOT/'resume_pids.json', new_workers)


def stop_own_workers(new_workers):
    for info in new_workers.values():
        try:
            os.killpg(info['pid'], signal.SIGTERM)
        except ProcessLookupError:
            pass


def main():
    assert os.environ.get('ALLOW_LOCAL_CALIBRATION') == '1'
    assert ROOT.is_dir() and not (ROOT/'resume_contract.json').exists()
    original = read(ROOT/'pids.json')
    assert {int(key) for key in original} == {48, 49, 50, 53, 54, 55}
    assert read(ROOT/'chain50/worker_terminal.json')['status'] == 'search_no_selected_candidate'
    assert read(ROOT/'chain48/worker_terminal.json')['status'] == 'selected_numerically_verified'
    packet_manifest = read(HERE/'preflight_chain48.json')['packet_manifest_sha256']
    import hashlib
    assert hashlib.sha256((PACKET/'manifest.json').read_bytes()).hexdigest() == packet_manifest
    contract = read(ROOT/'launch_contract.json')
    hard_deadline = contract['hard_deadline_epoch']
    assert time.time() < hard_deadline
    start_swap = psutil.swap_memory()
    swapout_base, swap_used_base = start_swap.sout, start_swap.used
    with (ROOT/'resume_contract.json').open('x') as handle:
        json.dump(dict(status='active', requested_chain_ids=[51, 56, 52, 57],
                       original_chain_ids=sorted(map(int, original)),
                       original_chain50_status='search_no_selected_candidate',
                       no_retries=True, same_batch_hard_deadline_epoch=hard_deadline,
                       started_epoch=time.time(), one_thread_each=True,
                       minimum_memory_headroom_bytes=MIN_FREE,
                       initial_worker_rss_budget_bytes=INITIAL_NEW_WORKER_BUDGET,
                       swap_delta_stop_bytes=SWAPOUT_LIMIT,
                       packet_manifest_sha256=packet_manifest), handle, indent=2)
        handle.write('\n')
    new_workers = {}
    stop_reason = None
    for stage, pair in enumerate(PAIRS, 1):
        if stage > 1:
            # Let the previous pair reach its normal RSS before the next headroom check.
            time.sleep(60)
        state = snapshot(new_workers, swapout_base, swap_used_base)
        budget = max(INITIAL_NEW_WORKER_BUDGET,
                     int(max((x['rss_bytes'] for x in state['workers']), default=0)*1.25))
        if time.time() >= hard_deadline:
            stop_reason = 'original_eight_hour_batch_cap'
        elif state['available_bytes'] < 2*budget+MIN_FREE:
            stop_reason = f'predicted_headroom_below_3GiB_before_resume_pair_{stage}'
        elif (state['swapout_since_resume_bytes'] > SWAPOUT_LIMIT or
              state['swap_used_since_resume_bytes'] > SWAPOUT_LIMIT):
            stop_reason = 'active_swap_growth_exceeded_1GiB'
        if stop_reason:
            break
        launch_pair(pair, new_workers, hard_deadline)
        write(ROOT/'resume_latest_progress.json', dict(status='active',
              launched=len(new_workers), planned=4, stop_reason=None,
              state=snapshot(new_workers, swapout_base, swap_used_base)))
    while new_workers and time.time() < hard_deadline:
        state = snapshot(new_workers, swapout_base, swap_used_base)
        write(ROOT/'resume_latest_progress.json', dict(status='monitor',
              launched=len(new_workers), planned=4, stop_reason=stop_reason,
              state=state))
        if all(row['terminal_status'] is not None or row['state'] in ('exited', psutil.STATUS_ZOMBIE)
               for row in state['workers']):
            break
        time.sleep(60)
    if time.time() >= hard_deadline:
        stop_reason = 'original_eight_hour_batch_cap'
        stop_own_workers(new_workers)
    write(ROOT/'resume_terminal.json', dict(status='all_four_additional_starts_launched'
          if len(new_workers) == 4 else 'partial_resume', launched=len(new_workers),
          planned=4, stop_reason=stop_reason, finished_epoch=time.time(),
          state=snapshot(new_workers, swapout_base, swap_used_base), no_retries=True))


if __name__ == '__main__':
    main()
