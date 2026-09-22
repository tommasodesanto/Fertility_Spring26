#!/usr/bin/env python3
"""Collect pinned earnings/entry battery evidence without solving or changing it."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text())


def optional(path):
    return read(path) if Path(path).is_file() else {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def solve_counts(ev):
    records = []
    for path in ev.glob('raw/repetition_*/stationary_solves.json'):
        records.extend(read(path))
    hb = optional(ev / 'raw/heartbeat.json')
    started = max(len(records), int(hb.get('stationary_solves', 0)))
    completed = sum(x.get('status') == 'completed' for x in records)
    return dict(started_solves=started, completed_solves=completed,
                unfinished_solves=max(0, started - completed),
                count_evidence='native heartbeat and per-solve records' if hb or records else 'no native solve receipt')


def verify_case(ev, plan, pin, manifest):
    score_dir = ev / 'scored_repetition_01'
    summary, score = read(ev / 'summary.json'), read(score_dir / 'score.json')
    run_path, initial_path = ev.parent / 'run_contract.json', ev.parent / 'initial_contract.json'
    run, initial = read(run_path), read(initial_path)
    receipt_path = score_dir / 'verified_evaluation_receipt.json'
    receipt = read(receipt_path)
    assert summary['status'] == 'verified_scored_candidate'
    assert summary['repetitions'] == 1
    assert summary['objective_canonical_sha256'] == manifest['objective_canonical_sha256']
    assert score['contract_sha256'] == manifest['objective_canonical_sha256']
    assert sha(run_path) == summary['run_contract_sha256']
    assert sha(initial_path) == summary['initial_solve_contract_sha256']
    assert initial['earnings_wealth_plan_sha256'] == pin['plan_sha256']
    assert initial['structural_candidate'] == plan['structural_parameters']
    assert canonical(initial['source_sha256']) == manifest['source_inventory_sha256']
    assert len(initial['source_sha256']) == 641
    objective = read(run['working_objective']['path'])
    assert sha(run['working_objective']['path']) == run['working_objective']['sha256']
    assert canonical(objective) == manifest['objective_canonical_sha256']
    assert canonical({k: v for k, v in objective.items() if k != 'parameter_restrictions'}) == manifest['target_system_sha256']
    assert plan['target_system_sha256'] == manifest['target_system_sha256']
    expected_sources = {}
    for name, item in run['objective_source_files'].items():
        expected_sources[name] = canonical(read(item['path'])) if item['hash_kind'] == 'canonical_json' else sha(item['path'])
    assert score['source_fingerprints'] == expected_sources
    assert receipt['source_fingerprints'] == expected_sources
    assert receipt['status'] == 'verified' and receipt['numerical_gates_verified'] is True
    assert canonical(receipt) == summary['evaluation_receipt_sha256'][0] == score['evaluation_receipt_sha256']
    assert receipt['checkpoint_sha256'] == score['checkpoint_sha256']
    assert score['schema'] == 'e5f_initial_minimum_distance_result_v1'
    targets, parameters = score['target_fit'], score['parameters']
    assert len(targets) == 13 and len(parameters) == 17
    assert score['scored_moment_count'] == 12 and sum(bool(x['scored']) for x in targets) == 12
    assert len({x['restriction_id'] for x in targets}) == 13
    assert len({x['parameter'] for x in parameters}) == 17
    with (score_dir / 'target_fit.csv').open() as stream:
        target_csv = {x['restriction_id']: x for x in csv.DictReader(stream)}
    with (score_dir / 'parameters.csv').open() as stream:
        parameter_csv = {x['parameter']: x for x in csv.DictReader(stream)}
    assert len(target_csv) == 13 and len(parameter_csv) == 17
    for row in targets:
        for field in ['target', 'model', 'gap', 'actual_weight', 'loss_contribution']:
            if row.get(field) is not None:
                assert float(target_csv[row['restriction_id']][field]) == float(row[field])
    for row in parameters:
        assert float(parameter_csv[row['parameter']]['estimate']) == float(row['estimate'])
    loss = float(score['loss'])
    assert math.isfinite(loss) and loss == float(summary['loss'])
    assert math.isclose(sum(float(x['loss_contribution']) for x in targets if x['scored']), loss, rel_tol=1e-12, abs_tol=1e-9)
    for row in targets:
        assert math.isclose(float(row['model']) - float(row['target']), float(row['gap']), rel_tol=1e-12, abs_tol=1e-12)
    for row in parameters:
        name = row['parameter']
        if name in plan['structural_parameters']:
            assert math.isclose(float(row['estimate']), float(plan['structural_parameters'][name]), rel_tol=1e-12, abs_tol=1e-12)
    graph_rows = summary['original_graphs']
    expected_names = set(run['expected_graph_filenames'])
    assert len(graph_rows) == len(expected_names) == 17
    assert {Path(x['path']).name for x in graph_rows} == expected_names
    for row in graph_rows:
        assert sha(ev / 'raw' / row['path']) == row['sha256']
    return score, summary


def collect(manifest_path, results_root, output):
    m = read(manifest_path)
    assert m['schema'] in {'e5f_earnings_entry_battery_manifest_v1', 'e5f_earnings_entry_resolution_manifest_v1'}
    output.mkdir(parents=True, exist_ok=False)
    inventory, valid = [], {}
    for stage in ('smoke', 'production'):
        for worker in m[stage]['cases']:
            cell, task = worker['cell_id'], worker['task_id']
            for index, pin in enumerate(worker['proposals'], 1):
                case = results_root / stage / cell / f'worker_{task:02d}' / f'proposal_{index:02d}_{pin["proposal_id"]}'
                ev = case / 'result/evaluation'
                row = dict(cell_id=cell, stage=stage, task_id=task, proposal_index=index,
                           case_id=pin['case_id'], case_directory=str(case), status='not_started', loss=None)
                if case.exists():
                    state = optional(case / 'status.json')
                    failure = optional(ev / 'failure.json')
                    row.update(solve_counts(ev))
                    row['status'] = 'failed_unscored' if failure or state.get('status') == 'incomplete' else 'running_or_incomplete'
                    row['failure_type'] = failure.get('error_type', state.get('error_type', ''))
                    row['failure_message'] = failure.get('error', state.get('error', ''))
                    if (ev / 'summary.json').is_file():
                        try:
                            assert sha(pin['plan_path']) == pin['plan_sha256']
                            plan = read(pin['plan_path'])
                            score, summary = verify_case(ev, plan, pin, m)
                            row.update(status='verified_scored', loss=float(score['loss']),
                                       normalized_objective_seconds=float(summary['elapsed_seconds']))
                            valid[pin['case_id']] = (ev, plan, pin, score, summary, row)
                        except (AssertionError, KeyError, ValueError, OSError, TypeError) as error:
                            row.update(status='collection_rejected', collection_error=f'{type(error).__name__}: {error}')
                inventory.append(row)
    selected, figure_manifest, common_targets = [], [], []
    for cell in m['cells']:
        choices = [v for v in valid.values() if v[-1]['cell_id'] == cell]
        if not choices:
            continue
        best = min(choices, key=lambda x: x[-1]['loss'])
        for tag, chosen in [('selected', best), *[('smoke', x) for x in choices if x[-1]['stage'] == 'smoke']]:
            ev, plan, pin, score, summary, row = chosen
            dest = output / cell / tag
            dest.mkdir(parents=True)
            score_dir = ev / 'scored_repetition_01'
            for name in ['score.json', 'target_fit.csv', 'parameters.csv', 'verified_evaluation_receipt.json']:
                shutil.copy2(score_dir / name, dest / name)
            shutil.copy2(ev / 'summary.json', dest / 'wrapper_summary.json')
            shutil.copy2(pin['plan_path'], dest / 'plan.json')
            for name in ['initial_contract.json', 'run_contract.json', 'runtime_contract.json']:
                if (ev.parent / name).is_file():
                    shutil.copy2(ev.parent / name, dest / name)
            for name in ['summary.json', 'lifecycle_2023.csv', 'stationary_solves.json']:
                path = ev / 'raw/repetition_01' / name
                if path.is_file():
                    shutil.copy2(path, dest / ('native_' + name))
            for name in ['entry_wealth.json', 'income_process.json', 'wealth_grid.json']:
                path = ev / name
                if path.is_file():
                    shutil.copy2(path, dest / name)
            annotated = []
            for item in score['parameters']:
                item = dict(item)
                name = item['parameter']
                if name in plan['parameter_bounds']:
                    lo, hi = plan['parameter_bounds'][name]
                    item.update(actual_lower=lo, actual_upper=hi,
                                near_actual_bound=min(float(item['estimate'])-lo, hi-float(item['estimate'])) <= .01*(hi-lo))
                else:
                    item.update(actual_lower=None, actual_upper=None, near_actual_bound=None)
                annotated.append(item)
            write_csv(dest / 'parameters_actual_bounds.csv', annotated)
            if tag == 'smoke':
                common_targets.extend(dict(cell_id=cell, **x) for x in score['target_fit'])
            if tag == 'selected':
                checkpoint = ev / 'raw/repetition_01/initial_state.pkl.gz'
                assert sha(checkpoint) == score['checkpoint_sha256']
                graphs = dest / 'standard_diagnostics'
                graphs.mkdir()
                for item in summary['original_graphs']:
                    source = ev / 'raw' / item['path']
                    target = graphs / source.name
                    shutil.copy2(source, target)
                    figure_manifest.append(dict(cell_id=cell, path=str(target.relative_to(output)), sha256=sha(target)))
                selected.append(dict(row, output=str(dest), checkpoint_sha256=score['checkpoint_sha256'],
                                     exact_repetitions_completed=0, grid_robustness_completed=False))
    attempted = [x for x in inventory if x['status'] != 'not_started']
    totals = {cell: {status: sum(x['cell_id']==cell and x['status']==status for x in inventory)
                    for status in sorted({x['status'] for x in inventory})} for cell in m['cells']}
    counts = {k: sum(x.get(k, 0) for x in attempted) for k in ['started_solves','completed_solves','unfinished_solves']}
    write(output / 'inventory.json', inventory)
    write_csv(output / 'inventory.csv', inventory)
    write(output / 'selection.json', selected)
    write(output / 'figure_manifest.json', figure_manifest)
    write_csv(output / 'common_smoke_targets.csv', common_targets)
    write(output / 'collection_summary.json', dict(cells=totals, stationary_solves=counts,
          manifest_sha256=sha(manifest_path), source_inventory_sha256=m['source_inventory_sha256'],
          objective_canonical_sha256=m['objective_canonical_sha256'], selected=selected,
          note='Finite paired design; selection includes smoke. Unfinished solves include currently running solves. No exact repeats or grid convergence asserted.'))
    # This permits matched-vector comparisons without treating unequal search pools as paired.
    matched = []
    for task in range(1,11):
        for index in range(1,7):
            rows = [x for x in attempted if x['stage']=='production' and x['task_id']==task and x['proposal_index']==index and x['status']=='verified_scored']
            if len(rows)>=2:
                matched.append(dict(task_id=task, proposal_index=index, cells={x['cell_id']:x['case_id'] for x in rows}, losses={x['cell_id']:x['loss'] for x in rows}))
    write(output / 'paired_completed_points.json', matched)
    (output / 'README.md').write_text('# Earnings and entry-wealth collection\n\nSelected points include smoke and production. Each selected cell has its complete 13-row target table, 17-parameter table with actual plan bounds, and the original 17 figures. The common smoke target table holds parameters fixed across cells. Missing cells have no verified scored result. See inventory and collection summary for failures, unstarted cases and unfinished solves. No model adoption, convergence, exact repetitions or grid robustness is implied.\n')
    return dict(attempted=len(attempted), verified=len(valid), selected_cells=len(selected), counts=counts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--results-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(collect(a.manifest,a.results_root,a.output),sort_keys=True))

if __name__=='__main__':
    main()
