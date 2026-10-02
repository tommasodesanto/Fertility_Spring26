"""Forward-only reconstruction of saved dated terminal states; no model solves."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import pickle
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REMOTE = Path('/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a')
FROZEN = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
PACKET = Path('output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1')
RESULTS = REMOTE / 'results'
CASES = {
    0: ('hard', 'control', 12), 1: ('hard', 'temporary', 12),
    2: ('hard', 'permanent', 12), 3: ('hard', 'control', 16),
    4: ('hard', 'temporary', 16), 5: ('hard', 'permanent', 16),
    6: ('quarter', 'control', 12), 7: ('quarter', 'temporary', 12),
    8: ('quarter', 'permanent', 12), 9: ('quarter', 'control', 16),
    10: ('quarter', 'temporary', 16), 11: ('quarter', 'permanent', 16),
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path, digest=None):
    if digest is not None and sha(path) != digest:
        raise ValueError(f'Packet SHA mismatch: {path}')
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def verify_forward_sources(selected_arm):
    """Match production's saved native preparation receipt for every case."""
    source_checks = {}
    for rel in (
        'code/model/tools/run_e5f_perfect_foresight_transition.py',
        'code/model/tools/run_e5f_open_population_transition.py',
        'code/model/tools/run_dynamic_population_transition.py',
    ):
        actual_path = FROZEN / rel
        actual = sha(actual_path)
        # Production binds the frozen tree first, then these two overlays.
        for overlay in ('/scratch/td2248/projects/grid_resolution_credit053_v2',
                        '/scratch/td2248/projects/normalized_floor_calibration_v1'):
            if (Path(overlay) / 'source' / rel).exists():
                raise ValueError(f'Production overlay replaces frozen helper: {rel}')
        if (REMOTE / 'source' / rel).exists():
            raise ValueError(f'Reviewed overlay replaces frozen helper: {rel}')
        receipts = {}
        for index, (arm, kind, horizon) in CASES.items():
            if arm != selected_arm:
                continue
            name = f'case_{index:02d}_{arm}_{kind}_h{horizon}'
            receipt = RESULTS / name / 'run/runtime/runtime_auth/runtime_preparation/native_preparation/preparation.json'
            source = json.loads(receipt.read_text())['current_source_files']
            if rel not in source or source[rel] != actual:
                raise ValueError(f'Production native source identity differs: {name}: {rel}')
            receipts[name] = source[rel]
        source_checks[rel] = {'loaded_path': str(actual_path), 'actual': actual,
                              'pinned': actual, 'matches': True,
                              'production_receipts': receipts}
    return source_checks


def main(out, selected_arm):
    out.mkdir(parents=True, exist_ok=True)
    # These are the pinned native forward functions used by the reviewed run.
    assert selected_arm in ('hard', 'quarter')
    engine = REMOTE / 'source' / PACKET / 'engines' / selected_arm
    sys.path[:0] = [str(engine), str(FROZEN / 'code/model/tools'), str(FROZEN / 'code/model')]
    from small_credit_lab.engine import solver, parameters, utils
    import run_e5f_perfect_foresight_transition as pf
    from run_e5f_perfect_foresight_transition import calendar, transition
    calendar.model = solver
    for name in ('independent_child_maturation_active', 'get_fecundity_by_age',
                 'readiness_settled_state', 'parent_age_maturation_active',
                 'readiness_childless_states', 'readiness_gate_active'):
        setattr(solver, name, getattr(parameters, name))
    solver.interp_indices = utils.interp_indices
    source_checks = verify_forward_sources(selected_arm)
    results = {}
    for index, (arm, kind, horizon) in CASES.items():
        if arm != selected_arm:
            continue
        name = f'case_{index:02d}_{arm}_{kind}_h{horizon}'
        run = RESULTS / name / 'run'
        if kind == 'control':
            mapping_path = run / 'control/mapping.json'
        else:
            root = json.loads((run / 'dated_path/root.json').read_text())
            assert root['converged'] and all(root['gates'].values()), name
            mapping_paths = sorted((run / 'dated_path').glob('mapping_*/mapping.json'))
            mapping_path = mapping_paths[-1]
        mapping = json.loads(mapping_path.read_text())
        assert len(mapping['rows']) == horizon and len(mapping['diagnostic_packets']) == 3
        assert mapping['diagnostic_packets'][-1]['period'] == horizon - 1
        last_pin = mapping['diagnostic_packets'][-1]
        last_path = run / Path(last_pin['path']).relative_to('/work/results/run')
        last = load(last_path, last_pin['sha256'])
        P, grid, ev, shared = (last[k] for k in ('parameters', 'b_grid', 'evaluation', 'shared'))
        assert last['period'] == horizon - 1
        # Starting queues come from the authenticated stationary 80% reference.
        reference_json = json.loads((run / 'reference/reference_reconstruction.json').read_text())
        initial_path = run / 'reference/selected_native_packet.pkl.gz'
        initial = load(initial_path, reference_json['checkpoint_sha256'])
        initial_state = pf.stationary_initial_state(
            initial['stationary_g_pre'],
            float(initial['stationary_g_pre'][:, :, :, 0].sum()),
            float(initial['evaluation'].births), initial['parameters'], 1 / 2.1)
        adjusted = pf.copy_birth_queue(initial_state.scheduled_entries)
        raw = pf.copy_birth_queue(initial_state.scheduled_raw_entries)
        queue_rows = []
        for row in mapping['rows']:
            due, adjusted = transition.advance_adult_entry_clock(
                adjusted, float(row['birth_children_topcode_adjusted']), 1 / 2.1,
                pf.entry_clock_timing(P))
            due_raw, raw = transition.advance_adult_entry_clock(
                raw, float(row['birth_children']), 1 / 2.1,
                pf.entry_clock_timing(P))
            queue_rows.append({'period': row['period'],
                'due_gap': due - float(row['effective_mature_entrant_flow_B']),
                'raw_due_gap': due_raw - float(row['raw_state_scheduled_mature_entrant_flow_B'])})
        if max(abs(q['due_gap']) + abs(q['raw_due_gap']) for q in queue_rows) > 1e-10:
            raise ValueError(f'Saved queue-flow replay differs: {name}')
        next_g, _, deaths, _ = transition.advance_sequential_calendar_distribution(
            ev, np.zeros(int(P.I)), P, grid, shared)
        shares = np.asarray(P.entry_shares, dtype=float).reshape(-1)
        shares /= float(shares.sum())
        entrants = float(mapping['rows'][-1]['entrant_flow_next']) * shares
        next_g[:, :, :, 0, :, :, :] = calendar.entrant_cohort(entrants, P, grid)
        mass_gap = float(next_g.sum() - (ev.g_post_fertility.sum() - deaths + entrants.sum()))
        if abs(mass_gap) > 1e-10:
            raise ValueError(f'Forward mass conservation differs: {name}: {mass_gap}')
        state = pf.PFInitialState(next_g, adjusted, raw)
        if kind == 'permanent':
            accepted = json.loads((run / 'permanent_terminal/accepted.json').read_text())
            endpoint = accepted['endpoint']
            stationary_packet = load(run / 'permanent_terminal/one_step/date_000/diagnostic_packet.pkl.gz')
            ref_g = stationary_packet['evaluation'].g_pre
            ref_P = stationary_packet['parameters']
            ref_births = float(stationary_packet['evaluation'].births)
            ref_price = float(endpoint['price'])
        else:
            ref_g = initial['stationary_g_pre']
            ref_P = initial['parameters']
            ref_births = float(initial['evaluation'].births)
            ref_price = float(json.loads((run /
                'reference/repeat_0/phase_b_ge/selected_root/closure.json').read_text())['price'])
        ref_state = pf.stationary_initial_state(ref_g,
            float(ref_g[:, :, :, 0].sum()), ref_births, ref_P, 1 / 2.1)
        native = SimpleNamespace(terminal_state=state,
            prices=np.asarray([row['asset_price'] for row in mapping['rows']]),
            rents=np.asarray([row['renter_price'] for row in mapping['rows']]))
        tolerance = 1e-6 if kind == 'control' else 1e-3
        check = pf.terminal_convergence_diagnostics(
            evaluation=native, psi_path=np.asarray([row['psi_child'] for row in mapping['rows']]),
            reference_state=ref_state, reference_entry_flow=float(ref_g[:, :, :, 0].sum()),
            reference_price=ref_price, reference_psi=float(ref_P.psi_child),
            base_parameters=ref_P,
            tolerances={k: tolerance for k in pf.DEFAULT_TERMINAL_TOLERANCES})
        entry = float(ref_g[:, :, :, 0].sum())
        raw_gap = float(np.max(np.abs(pf.birth_queue_values(raw) -
            pf.birth_queue_values(ref_state.scheduled_raw_entries)))) / max(entry, 1e-15)
        check['raw_queue_maximum_relative_gap'] = raw_gap
        check['raw_queue_pass'] = raw_gap <= tolerance
        check['all_checks_pass'] = bool(check['all_checks_pass'] and check['raw_queue_pass'])
        if kind == 'control':
            official = json.loads((run / 'completed.json').read_text())['terminal']
            for key in official['metrics']:
                if abs(float(check['metrics'][key]) - float(official['metrics'][key])) > 1e-10:
                    raise ValueError(f'Control terminal metric differs: {name}: {key}')
            if abs(raw_gap - official['raw_queue_maximum_relative_gap']) > 1e-10:
                raise ValueError(f'Control raw queue metric differs: {name}')
        results[name] = {'mapping_path': str(mapping_path), 'last_packet_sha256': last_pin['sha256'],
            'forward_only': True, 'model_solves': 0,
            'maximum_queue_replay_gap': max(abs(q['due_gap']) + abs(q['raw_due_gap']) for q in queue_rows),
            'mass_gap': mass_gap, 'terminal': check,
            'control_exact_metric_reproduction': kind == 'control'}
        print(name, check['all_checks_pass'], check['metrics'], 'raw_queue', raw_gap, flush=True)
    (out / 'recovered_terminal.json').write_text(json.dumps({
        'source_checks': source_checks, 'results': results}, indent=2, default=lambda x: x.tolist()
        if isinstance(x, np.ndarray) else x) + '\n')


if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[3] == '--provenance-only':
        path = Path(sys.argv[1]) / 'recovered_terminal.json'
        result = json.loads(path.read_text())
        result['source_checks'] = verify_forward_sources(sys.argv[2])
        path.write_text(json.dumps(result, indent=2) + '\n')
        print(sys.argv[2], 'all production native helper hashes match')
    else:
        main(Path(sys.argv[1]), sys.argv[2])
