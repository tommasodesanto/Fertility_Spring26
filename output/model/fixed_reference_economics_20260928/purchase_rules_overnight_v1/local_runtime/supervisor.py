"""Supervise ten distinct local chains with a measured two-worker ramp."""
from __future__ import annotations

import argparse
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import psutil

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
PAIRS = [(48, 53), (49, 54), (50, 55), (51, 56), (52, 57)]
ONE_GIB = 1024 ** 3
BATCH_SECONDS = 8 * 3600
CHAIN_SECONDS = 4 * 3600
MIN_FREE = 3 * ONE_GIB
INITIAL_NEW_WORKER_BUDGET = int(3.8 * ONE_GIB)
SWAPOUT_LIMIT = ONE_GIB
STAGE_WAIT = 900


def write(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def first_passed_case(folder):
    candidates = []
    best = folder/'best_so_far.json'
    if best.is_file():
        candidates.append(read(best).get('best'))
    cases = folder/'cases.json'
    if cases.is_file():
        candidates.extend(read(cases))
    latest = folder/'latest_completed.json'
    if latest.is_file():
        candidates.append(read(latest).get('latest'))
    for row in candidates:
        if not isinstance(row, dict) or row.get('status') != 'passed':
            continue
        if not math.isfinite(float(row.get('loss', float('nan')))):
            continue
        report = Path(row.get('report', ''))
        if not report.is_dir():
            continue
        # The first case is a useful native SMM evaluation, not a separate smoke.
        if all((report/name).is_file() and len((report/name).read_text().splitlines()) == count + 1
               for name, count in [('target_fit.csv', 14), ('parameters.csv', 31)]):
            return dict(label=row['label'], loss=float(row['loss']), report=str(report))
    return None


def current_state(workers, swapout_base, swap_used_base):
    memory = psutil.virtual_memory()
    swap = psutil.swap_memory()
    rows = []
    for chain, info in workers.items():
        pid = info['pid']
        try:
            proc = psutil.Process(pid)
            state = proc.status()
            rss = proc.memory_info().rss + sum(
                child.memory_info().rss for child in proc.children(recursive=True)
                if child.is_running())
        except psutil.Error:
            state, rss = 'exited', 0
        folder = Path(info['out'])
        receipt = first_passed_case(folder/'search')
        terminal = read(folder/'worker_terminal.json') if (folder/'worker_terminal.json').is_file() else None
        rows.append(dict(chain=chain, pid=pid, state=state, rss_bytes=rss,
                         first_passed_case=receipt,
                         completed=terminal is not None and terminal['status'] == 'selected_numerically_verified',
                         failed=terminal is not None and terminal['status'] != 'selected_numerically_verified'))
    return dict(time_epoch=time.time(), available_bytes=memory.available,
                memory_percent=memory.percent, swap_used_bytes=swap.used,
                swap_used_since_launch_bytes=max(0, swap.used - swap_used_base),
                swapout_since_launch_bytes=max(0, swap.sout - swapout_base), workers=rows)


def launch_pair(pair, root, workers, batch_deadline):
    for chain in pair:
        out = root/f'chain{chain}'
        assert not out.exists(), out
        log = root/f'chain{chain}.log'
        handle = log.open('w')
        env = dict(os.environ, ALLOW_LOCAL_CALIBRATION='1')
        deadline = min(time.time() + CHAIN_SECONDS, batch_deadline)
        command = [str(HERE/'launch_one.sh'), '--chain', str(chain), '--out', str(out),
                   '--deadline-epoch', str(deadline)]
        proc = subprocess.Popen(command, cwd=PACKET.parents[3], env=env,
                                stdin=subprocess.DEVNULL, stdout=handle,
                                stderr=subprocess.STDOUT, start_new_session=True)
        handle.close()
        workers[chain] = dict(pid=proc.pid, out=str(out), log=str(log),
                              started_epoch=time.time(), deadline_epoch=deadline,
                              command=command)
    write(root/'pids.json', workers)


def stop_own_workers(workers):
    for info in workers.values():
        try:
            os.killpg(info['pid'], signal.SIGTERM)
        except ProcessLookupError:
            pass
    time.sleep(10)
    for info in workers.values():
        try:
            os.killpg(info['pid'], signal.SIGKILL)
        except ProcessLookupError:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-root', type=Path, required=True)
    args = parser.parse_args()
    root = args.out_root.resolve()
    assert root.is_relative_to(HERE.resolve()) and not root.exists()
    assert os.getenv('ALLOW_LOCAL_CALIBRATION') == '1'
    plan = read(HERE/'local_plan.json')
    assert plan['local_chain_ids'] == list(range(48, 58))
    assert len({tuple(plan['nearby_starts'][j]['parameters'][k] for k in plan['free_coordinates'])
                for j in range(48, 58)}) == 10
    parent = read(PACKET/'plan.json')
    for key in ('base_target_contract', 'target_fingerprint', 'weight_fingerprint', 'profiles'):
        assert plan[key] == parent[key], key
    manifest_bytes = (PACKET/'manifest.json').read_bytes()
    import hashlib
    manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()
    for chain in (48, 53):
        receipt = read(HERE/f'preflight_chain{chain}.json')
        assert receipt['status'] == 'authenticated_zero_solve_context'
        assert receipt['packet_manifest_sha256'] == manifest_sha
        assert receipt['lifecycle_solves'] == 0
    root.mkdir(parents=True)
    started = time.time()
    batch_deadline = started + BATCH_SECONDS
    initial_swap = psutil.swap_memory()
    swapout_base = initial_swap.sout
    swap_used_base = initial_swap.used
    workers = {}
    launch_history = []
    write(root/'launch_contract.json', dict(status='active', started_epoch=started,
          hard_deadline_epoch=batch_deadline, chain_seconds=CHAIN_SECONDS,
          reserve_seconds=900, ordered_pairs=PAIRS, one_thread_each=True,
          initial_worker_rss_budget_bytes=INITIAL_NEW_WORKER_BUDGET,
          minimum_memory_headroom_bytes=MIN_FREE,
          swapout_delta_stop_bytes=SWAPOUT_LIMIT, no_automatic_retry=True,
          packet_manifest_sha256=manifest_sha))
    stop_reason = None
    for stage, pair in enumerate(PAIRS, 1):
        gate_wait_started = time.time()
        while True:
            state = current_state(workers, swapout_base, swap_used_base)
            write(root/'latest_progress.json', dict(stage=stage, launched=len(workers),
                  planned=10, next_pair=pair, state=state, launch_history=launch_history))
            if time.time() >= batch_deadline:
                stop_reason = 'eight_hour_batch_hard_cap'
                break
            if state['available_bytes'] < MIN_FREE:
                stop_reason = 'actual_free_headroom_below_3GiB'
                break
            if state['swapout_since_launch_bytes'] > SWAPOUT_LIMIT or state['swap_used_since_launch_bytes'] > SWAPOUT_LIMIT:
                stop_reason = 'swap_exceeded_1GiB'
                break
            if stage > 1:
                prior = PAIRS[stage - 2]
                if any(next(r for r in state['workers'] if r['chain'] == j)['failed'] for j in prior):
                    stop_reason = 'prior_pair_failed'
                    break
                if not all(next(r for r in state['workers'] if r['chain'] == j)['first_passed_case'] for j in prior):
                    if time.time() - gate_wait_started > STAGE_WAIT:
                        stop_reason = 'prior_pair_no_useful_case_within_15m'
                        break
                    time.sleep(30)
                    continue
            # Before each new pair, budget their expected RSS as well as 3 GiB OS/app headroom.
            observed = [r['rss_bytes'] for r in state['workers'] if r['rss_bytes'] > 0]
            budget = max(INITIAL_NEW_WORKER_BUDGET,
                         int(max(observed, default=0) * 1.25))
            if state['available_bytes'] < 2 * budget + MIN_FREE:
                stop_reason = f'predicted_headroom_below_3GiB_before_pair_{stage}'
                break
            launch_pair(pair, root, workers, batch_deadline)
            launch_history.append(dict(stage=stage, pair=pair, launched_epoch=time.time(),
                                       available_before_bytes=state['available_bytes'],
                                       per_worker_rss_budget_bytes=budget))
            write(root/'latest_progress.json', dict(stage=stage, launched=len(workers),
                  planned=10, next_pair=PAIRS[stage] if stage < len(PAIRS) else None,
                  state=current_state(workers, swapout_base, swap_used_base), launch_history=launch_history))
            break
        if stop_reason:
            break
    # Once ramping ends, keep a durable heartbeat until own workers finish.
    while workers and time.time() < batch_deadline:
        state = current_state(workers, swapout_base, swap_used_base)
        write(root/'latest_progress.json', dict(stage='monitor', launched=len(workers),
              planned=10, stop_reason=stop_reason, state=state,
              launch_history=launch_history))
        if all(row['completed'] or row['failed'] or row['state'] in ('exited', psutil.STATUS_ZOMBIE)
               for row in state['workers']):
            break
        time.sleep(60)
    if time.time() >= batch_deadline:
        stop_reason = 'eight_hour_batch_hard_cap'
        stop_own_workers(workers)
    state = current_state(workers, swapout_base, swap_used_base)
    write(root/'terminal.json', dict(status='all_ten_launched' if len(workers) == 10 else 'partial_ramp',
          stop_reason=stop_reason, launched=len(workers), planned=10,
          finished_epoch=time.time(), state=state, launch_history=launch_history,
          no_automatic_retry=True))


if __name__ == '__main__':
    main()
