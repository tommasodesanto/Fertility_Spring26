"""Read-only collection and health check for the 48 matched timing chains.

``--fetch`` copies small receipts and CSVs from Torch, omitting every NPZ and
PNG. ``--winner-artifacts`` subsequently fetches only each verified arm
winner's native root/repeat packets, including arrays and the standard plots.
No command submits, cancels, resumes, or changes a remote run.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import subprocess
import time
from pathlib import Path

from collect_torch import read, sha

PROJECT = Path(__file__).resolve().parents[3]
PACKET = PROJECT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1'
PLAN = PACKET / 'collection_plan.json'
ARMS = ('original', 'alternative')


def slots(plan: dict) -> list[dict]:
    rows = []
    for arm in ARMS:
        for chain in range(24):
            stage = 'v2' if chain < 4 else 'v3'
            task = (chain if arm == 'original' else chain + 4) if stage == 'v2' else (
                chain - 4 if arm == 'original' else chain + 16)
            rows.append(dict(arm=arm, chain=chain, stage=stage, task=task,
                name=f'production_{arm}_chain_{chain}',
                remote=f"{plan['stages'][stage]['remote_root']}/results/production_{arm}_chain_{chain}"))
    assert len(rows) == len({(r['stage'], r['task']) for r in rows}) == 48
    return rows


def report_path(folder: Path, remote_report: str, expected_remote: str) -> Path:
    prefixes = (expected_remote + '/run/', '/work/results/run/')
    prefix = next((p for p in prefixes if remote_report.startswith(p)), None)
    if prefix is None or not remote_report.endswith('/selected_root'):
        raise ValueError('unexpected native report path')
    relative = Path(remote_report[len(prefix):])
    if not relative.parts or '..' in relative.parts or relative.is_absolute():
        raise ValueError('unsafe native report path')
    return folder / 'run' / relative


def table(path: Path) -> list[dict]:
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def classify(folder: Path, slot: dict, plan: dict) -> tuple[dict, list[dict] | None, list[dict] | None]:
    row = dict(arm=slot['arm'], chain=slot['chain'], stage=slot['stage'], task=slot['task'],
               name=slot['name'], status='unfinished', reason='no local terminal receipt')
    if not folder.is_dir():
        return row, None, None
    start_path = folder / 'launcher_start.json'
    terminal_path = folder / 'launcher_terminal.json'
    contract_path = folder / 'run/start_contract.json'
    completed_path = folder / 'run/completed.json'
    heartbeat_path = folder / 'run/heartbeat.json'
    if heartbeat_path.exists():
        heartbeat = read(heartbeat_path)
        row['heartbeat_status'] = heartbeat.get('status')
        row['heartbeat_age_seconds'] = max(0, round(time.time() - float(heartbeat['epoch'])))
    try:
        if start_path.exists():
            start = read(start_path)
            if (start.get('mode'), start.get('arm'), int(start.get('chain', -1))) != (
                'production', slot['arm'], slot['chain']):
                raise ValueError('launcher identity mismatch')
            # Slurm may assign a distinct SLURM_JOB_ID to each array element;
            # the submission's base array ID is not an element identity check.
            row['slurm_job_id'] = start.get('slurm_job_id')
        if contract_path.exists():
            contract = read(contract_path)
            if (contract.get('arm'), int(contract.get('chain', -1))) != (slot['arm'], slot['chain']):
                raise ValueError('contract chain identity mismatch')
            for key in ('target_fingerprint', 'weight_fingerprint', 'selected_source_sha256'):
                row[key] = contract.get(key)
                if row[key] != plan[key]:
                    raise ValueError(f'mixed {key}')
            expected_start = plan['stages'][slot['stage']]['starts_file_sha256']
            expected_count = 4 if slot['stage'] == 'v2' else 24
            if slot['stage'] == 'v3' and (contract.get('starts_file_sha256') != expected_start or
                                          contract.get('starts_count') != expected_count):
                raise ValueError('expanded start-table identity mismatch')
            if slot['stage'] == 'v2' and ('starts_file_sha256' in contract or 'starts_count' in contract):
                raise ValueError('legacy start-table schema drift')
            if len(contract.get('all_starts', [])) != expected_count or len(contract.get('free_coordinates', [])) != 10:
                raise ValueError('start count or coordinate count mismatch')
            if contract.get('objective_calls_max') != 250 or contract.get('reserve_seconds') != 1800:
                raise ValueError('search budget mismatch')
        if terminal_path.exists():
            terminal = read(terminal_path)
            if (terminal.get('mode'), terminal.get('arm'), int(terminal.get('chain', -1))) != (
                'production', slot['arm'], slot['chain']):
                raise ValueError('terminal identity mismatch')
            row['exit_code'] = terminal.get('exit_code')
            if terminal.get('exit_code') != 0:
                row.update(status='terminal_failure', reason=f"launcher exit {terminal.get('exit_code')}")
                return row, None, None
        elif (folder / 'run/failure.json').exists():
            row.update(status='failure_receipt_without_terminal', reason=read(folder / 'run/failure.json').get('message'))
            return row, None, None
        else:
            row['reason'] = 'running or terminal receipt not fetched'
            return row, None, None
        if not all(p.exists() for p in (start_path, contract_path, completed_path)):
            row.update(status='incomplete_terminal', reason='missing start, contract, or completion receipt')
            return row, None, None
        completed = read(completed_path)
        row['objective_calls'] = completed.get('objective_calls')
        for key in ('target_fingerprint', 'weight_fingerprint'):
            if completed.get(key) != plan[key]:
                raise ValueError(f'mixed {key} in completion receipt')
        if completed.get('status') == 'no_admissible_candidate':
            row.update(status='no_admissible_candidate', reason='search found no admissible selected point')
            return row, None, None
        if completed.get('status') == 'selected_native_passed_search_loss_differs':
            row.update(status='review_needed', reason='search and native losses differ',
                       native_loss=completed.get('native_loss'))
            return row, None, None
        if completed.get('status') != 'selected_numerically_verified':
            row.update(status='unverified_terminal', reason=f"completion status {completed.get('status')}")
            return row, None, None
        # Reuse the original collector's receipt/report contract, with compact
        # copies of the embedded repeat/plot evidence instead of 816 PNG files.
        latest = read(folder / 'run/latest_completed.json')
        best = read(folder / 'run/best_so_far.json')
        search = read(folder / 'run/search_completed.json')
        inputs = read(folder / 'run/input_contract.json')
        if latest.get('completed_full_ge', 0) < 1 or best.get('best') is None or search.get('selected') is None:
            raise ValueError('missing admissible search checkpoint')
        for source in (inputs, completed, search):
            for key in ('target_fingerprint', 'weight_fingerprint'):
                if source.get(key) != plan[key]:
                    raise ValueError(f'mixed {key}')
        if slot['stage'] == 'v3' and search.get('starts_file_sha256') != plan['stages'][slot['stage']]['starts_file_sha256']:
            raise ValueError('search start-table identity mismatch')
        if slot['stage'] == 'v2' and 'starts_file_sha256' in search:
            raise ValueError('legacy search start-table schema drift')
        post = completed.get('selected_postcheck', {})
        repeat = completed.get('repeat', {})
        if post.get('status') != 'passed' or repeat.get('status') != 'exact_full_ge_repeat_passed':
            raise ValueError('native postcheck or exact repeat failed')
        if (repeat.get('target_rows') != 14 or repeat.get('parameter_rows') != 31 or
                len(repeat.get('standard_plot_hashes', {})) != 17):
            raise ValueError('exact repeat report counts differ')
        child_path = folder / 'run/native_postcheck/completed.json'
        child = read(child_path)
        if (child.get('status') != 'full_native_postcheck_passed' or
                child.get('search_receipt_sha256') != sha(folder / 'run/search_completed.json') or
                child.get('target_fingerprint') != plan['target_fingerprint'] or
                child.get('weight_fingerprint') != plan['weight_fingerprint']):
            raise ValueError('fresh native child receipt mismatch')
        if slot['stage'] == 'v3' and child.get('starts_file_sha256') != plan['stages']['v3']['starts_file_sha256']:
            raise ValueError('expanded child start-table identity mismatch')
        if slot['stage'] == 'v2' and 'starts_file_sha256' in child:
            raise ValueError('legacy child start-table schema drift')
        report = report_path(folder, post['report'], slot['remote'])
        fits, params = table(report / 'target_fit.csv'), table(report / 'parameters.csv')
        if len(fits) != 14 or len(params) != 31:
            raise ValueError('native report row count differs')
        if not {'parameter', 'estimate', 'lower', 'upper', 'near_bound', 'status'} <= set(params[0]):
            raise ValueError('parameter bounds/status columns missing')
        if [r['moment'] for r in fits] != plan['target_rows']:
            raise ValueError('native target rows differ from fixed contract')
        native_loss = sum(float(r['loss_contribution'] or 0) for r in fits)
        if not math.isfinite(native_loss) or abs(native_loss - float(completed['native_loss'])) > 1e-8:
            raise ValueError('native loss/table discrepancy')
        search_loss = float(best['best']['loss'])
        row.update(native_loss=native_loss, search_loss=search_loss,
            completed_sha256=sha(completed_path), terminal_sha256=sha(terminal_path),
            report=str(report), repeat_status=repeat['status'])
        if abs(native_loss - search_loss) > 1e-8:
            row.update(status='review_needed', reason='search/native loss discrepancy')
            return row, None, None
        row.update(status='verified', reason='fresh native point and exact repeat passed')
        return row, fits, params
    except (OSError, KeyError, TypeError, ValueError, AssertionError) as exc:
        row.update(status='contract_failure', reason=str(exc))
        return row, None, None


def fetch_one(slot: dict, root: Path, winner: bool = False) -> None:
    folder = root / slot['name']
    folder.mkdir(parents=True, exist_ok=True)
    if winner:
        # The parent chain may contain hundreds of exploratory candidate arrays.
        # Fetch only the two native reports for a selected winner.
        completed = read(folder / 'run/completed.json')
        remote_report = completed['selected_postcheck']['report']
        local_report = report_path(folder, remote_report, slot['remote'])
        relative_parent = local_report.relative_to(folder / 'run').parent
        remote_parent = f"{slot['remote']}/run/{relative_parent}"
        for name in ('selected_root', 'selected_repeat_final'):
            remote = f'{remote_parent}/{name}'
            local = local_report.parent / name
            local.mkdir(parents=True, exist_ok=True)
            subprocess.run(['rsync', '-a', '--', f'torch:{remote}/', str(local) + '/'], check=True)
        return
    cmd = [
        'rsync', '-a', '--exclude=*.npz', '--exclude=*.png', '--exclude=*.npy',
        '--exclude=*.pkl', '--exclude=*.pickle', '--exclude=numba_cache/', '--']
    subprocess.run(cmd + [f"torch:{slot['remote']}/", str(folder) + '/'], check=True)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict], columns: tuple[str, ...]) -> list[str]:
    def cell(value: object) -> str:
        return str(value).replace('|', '\\|').replace('\n', ' ')
    return ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join('---' for _ in columns) + ' |',
            *['| ' + ' | '.join(cell(row.get(column, '')) for column in columns) + ' |'
              for row in rows]]


def collect(root: Path, plan: dict, *, fetch: bool = False, winner_artifacts: bool = False) -> dict:
    all_slots = slots(plan)
    if fetch:
        for slot in all_slots:
            fetch_one(slot, root)
    assessments = [classify(root / s['name'], s, plan) for s in all_slots]
    rows = [r for r, _, _ in assessments]
    if any(r['status'] == 'contract_failure' for r in rows):
        state = 'failed_contract'
    elif all(r['status'] in ('verified', 'no_admissible_candidate') for r in rows):
        state = 'complete'
    elif any(r['status'] in ('terminal_failure', 'failure_receipt_without_terminal', 'incomplete_terminal',
                            'unverified_terminal', 'review_needed') for r in rows):
        state = 'review_needed'
    else:
        state = 'unfinished'
    winners = {}
    winner_tables = {}
    for arm in ARMS:
        candidates = [item for item in assessments if item[0]['arm'] == arm and item[0]['status'] == 'verified']
        if candidates:
            best = min(candidates, key=lambda item: (item[0]['native_loss'], item[0]['chain']))
            winners[arm] = {k: best[0][k] for k in ('arm', 'chain', 'stage', 'native_loss', 'search_loss', 'report')}
            winner_tables[arm] = (best[1], best[2])
            write_csv(root / f'{arm}_target_fit.csv', best[1])
            write_csv(root / f'{arm}_parameters.csv', best[2])
            if winner_artifacts:
                slot = next(s for s in all_slots if s['arm'] == arm and s['chain'] == best[0]['chain'])
                fetch_one(slot, root, winner=True)
                report = Path(best[0]['report'])
                if len(list((report / 'standard_diagnostics').glob('*.png'))) != 17:
                    raise ValueError(f'{arm} winner plot packet incomplete')
    receipt = dict(status=state, expected_chains=48, counts={s: sum(r['status'] == s for r in rows)
                   for s in sorted({r['status'] for r in rows})}, winners=winners, chains=rows,
                   target_fingerprint=plan['target_fingerprint'], weight_fingerprint=plan['weight_fingerprint'],
                   selected_source_sha256=plan['selected_source_sha256'])
    (root / 'collection.json').write_text(json.dumps(receipt, indent=2) + '\n')
    lines = ['# Matched soft timing calibration readout', '',
             f"Collection status: **{state}**; {receipt['counts']}. No calibration is adopted by this collector.", '']
    for arm in ARMS:
        winner = winners.get(arm)
        lines.append(f"## {arm}")
        if winner:
            lines += [f"Lowest verified native loss: **{winner['native_loss']:.12g}** "
                      f"(chain {winner['chain']}, {winner['stage']}).",
                      f"Full CSVs: `{arm}_target_fit.csv` and `{arm}_parameters.csv`.", '']
            fits, params = winner_tables[arm]
            lines += ['### Target fit', ''] + markdown_table(fits, (
                'moment', 'role', 'target', 'model', 'gap', 'weight', 'loss_contribution')) + ['',
                '### Parameters and restrictions', ''] + markdown_table(params, (
                'parameter', 'estimate', 'lower', 'upper', 'near_bound', 'status')) + ['']
        else:
            lines += ['No numerically verified candidate yet.', '']
    lines += ['Search optimization convergence is not certified by a selected-point native check.', '']
    (root / 'RESULTS.md').write_text('\n'.join(lines))
    return receipt


def status(plan: dict, root: Path) -> None:
    jobs = [str(plan['stages'][s].get('job_id')) for s in ('v2', 'v3') if plan['stages'][s].get('job_id')]
    if jobs:
        result = subprocess.run(['ssh', 'torch', 'squeue', '-h', '-j', ','.join(jobs),
                                 '-o', shlex.quote('%i|%T|%M')], capture_output=True, text=True, check=True)
        print(result.stdout.rstrip() or 'No matching active Slurm tasks')
    roots = {stage: plan['stages'][stage]['remote_root'] for stage in ('v2', 'v3')}
    script = '''import json,sys,time
from pathlib import Path
roots=json.loads(sys.argv[1])
now=time.time()
for stage,root in roots.items():
 for arm in ('original','alternative'):
  chains=range(4) if stage=='v2' else range(4,24)
  for chain in chains:
   folder=Path(root)/'results'/f'production_{arm}_chain_{chain}'
   path=folder/'run/heartbeat.json'
   terminal=folder/'launcher_terminal.json'
   if path.exists():
    try:
     h=json.loads(path.read_text())
     print(f"{stage} {arm} {chain:02d}: {h.get('status')} age={max(0,round(now-float(h['epoch'])))}s terminal={terminal.exists()}")
    except (OSError,ValueError,KeyError):
     print(f"{stage} {arm} {chain:02d}: malformed heartbeat")
   elif terminal.exists():
    print(f"{stage} {arm} {chain:02d}: terminal, heartbeat unavailable")
   else:
    print(f"{stage} {arm} {chain:02d}: no heartbeat or terminal receipt")
'''
    result = subprocess.run(['ssh', 'torch', 'python3', '-c', shlex.quote(script),
                             shlex.quote(json.dumps(roots))],
                            capture_output=True, text=True, check=True)
    print(result.stdout.rstrip() or 'No remote heartbeat receipts yet')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=PACKET / 'collection')
    parser.add_argument('--fetch', action='store_true')
    parser.add_argument('--status', action='store_true')
    parser.add_argument('--winner-artifacts', action='store_true')
    args = parser.parse_args()
    plan = read(PLAN)
    if args.status:
        status(plan, args.root)
        return
    args.root.mkdir(parents=True, exist_ok=True)
    receipt = collect(args.root, plan, fetch=args.fetch, winner_artifacts=args.winner_artifacts)
    print(json.dumps({k: receipt[k] for k in ('status', 'counts', 'winners')}))
    if receipt['status'] == 'failed_contract':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
