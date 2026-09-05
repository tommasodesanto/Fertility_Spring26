"""Bounded compiled replay of hash-pinned saved solutions after policy repair.

Use the isolated repaired snapshot via PYTHONPATH and pass a pinned manifest.
Old incomplete policies are never reused: a fresh Bellman solve must reproduce
all their arrays before its self-contained policy is subjected to interference.
"""
import argparse
import gzip
import hashlib
import json
import pickle
import time
from pathlib import Path

import numpy as np
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as calibration
import run_e5f_independent_numerical_audit as audit
from intergen_eqscale_seq_optimized import solver


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    assert sha(args.manifest) == args.manifest_sha256
    manifest = json.loads(args.manifest.read_text())
    assert sha(__file__) == manifest['replay_sha256']
    actual = calibration.code_fingerprint_contract(solver)
    assert actual == manifest['scientific_contract']
    for relative, digest in manifest['helper_sha256'].items():
        assert sha(Path(calendar.__file__).parent / relative) == digest
    assert len(manifest['cases']) == 5
    args.outdir.mkdir(parents=True, exist_ok=False)
    calendar.model = solver
    calendar.apply_fertility = transition.apply_sequential_fertility
    rows = []
    for case in manifest['cases']:
        started = time.perf_counter()
        path = Path(case['checkpoint'])
        assert sha(path) == case['sha256']
        print('CASE_START', case['name'], flush=True)
        with (gzip.open if path.suffix == '.gz' else open)(path, 'rb') as stream:
            packet = pickle.load(stream)
        P, grid = packet['parameters'], packet['b_grid']
        reference = packet['evaluation']
        shared = solver.precompute_shared(P, grid)
        continuation_reference = P._fert2_probs.copy()
        counter = calendar.SolveCounter()
        policy = calendar.solve_policy(reference.policy.price, P, grid, shared, counter)
        compared = []
        for name in ('V', 'c_pol', 'hR_pol', 'bp_pol', 'tenure_choice',
                     'tenure_probs', 'loc_probs', 'fert_probs', 'fert_value', 'price'):
            left, right = getattr(policy, name), getattr(reference.policy, name)
            np.testing.assert_array_equal(left, right, err_msg=case['name'] + '/' + name)
            compared.append(name)
        np.testing.assert_array_equal(policy.fert2_probs, continuation_reference)
        # Use the exact inherited input if retained; otherwise the saved prepared
        # input. Report that distinction, and require exact reproduction below.
        inherited = packet['state'].g_pre if 'state' in packet else reference.g_pre
        def evaluate(pol):
            return calendar.evaluate_period(pol.price, inherited, P, grid, shared,
                counter, packet['supply_rule'], supplied_policy=pol)
        first = evaluate(policy)
        for name in ('g_pre', 'g_post_fertility', 'g_current', 'births_by_loc',
                     'demand_by_loc', 'supply_by_loc'):
            np.testing.assert_array_equal(getattr(first, name), getattr(reference, name),
                                          err_msg=case['name'] + '/' + name)
        assert first.births == reference.births
        assert first.relative_market_residual == reference.relative_market_residual
        # Deliberately corrupt only the last-solve scratch probabilities.
        P._fert2_probs[...] = .99
        for pol in (policy, pickle.loads(pickle.dumps(policy, protocol=5))):
            again = evaluate(pol)
            for name in ('g_pre', 'g_post_fertility', 'g_current', 'births_by_loc',
                         'demand_by_loc', 'supply_by_loc'):
                np.testing.assert_array_equal(getattr(first, name), getattr(again, name))
            assert first.births == again.births
            assert first.relative_market_residual == again.relative_market_residual
        assert counter.bellman == 1
        case_out = args.outdir / case['name']
        case_out.mkdir()
        packet.update(evaluation=first, shared=shared)
        audit.standard_diagnostics(packet, case_out, validate_production_young=False)
        graphs = sorted((case_out / 'standard_diagnostics').glob('*.png'))
        assert len(graphs) == 17
        row = dict(name=case['name'], status='pass', exact_policy_arrays=compared,
                   exact_continuation=True, exact_saved_distributions=True,
                   exact_saved_births=True, exact_interference_and_pickle_replays=True,
                   bellman_solves=counter.bellman, checkpoint_sha256=case['sha256'],
                   used_original_preparation_input='state' in packet,
                   graphs={p.name: sha(p) for p in graphs},
                   elapsed_seconds=time.perf_counter()-started)
        rows.append(row)
        (args.outdir / 'latest_completed_case.json').write_text(json.dumps(row, indent=2)+'\n')
        (args.outdir / 'summary.json').write_text(json.dumps(dict(
            status='complete' if len(rows)==5 else 'running', cases=rows,
            manifest_sha256=args.manifest_sha256, no_calibration_search=True,
            production_promoted=False), indent=2)+'\n')
        print('CASE_PASS', case['name'], row['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
