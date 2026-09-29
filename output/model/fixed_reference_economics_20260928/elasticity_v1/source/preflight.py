#!/usr/bin/env python3
"""Zero-solve common-grid support check for the five prescribed prices."""
import argparse
import copy
import json
import os
from pathlib import Path
import sys
import time
import traceback

import run_credit as credit
import run_fixed_price as base
import natural_credit as adapter

FACTORS = (.98, .99, 1., 1.01, 1.02)
GRID_SHA = 'a3b81c459290554f634d0c28409ab846364c0ecf7270bd164f975f092901cbe2'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    base.require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    base.require(not out.exists(), 'Preflight output exists')
    out.mkdir(parents=True)
    started = time.monotonic()
    source = Path(__file__).resolve().parent
    grid_path = source / 'credit_grid.json'
    base.require(base.sha(grid_path) == GRID_SHA, 'Grid differs from authenticated credit-v1 262-node input')
    try:
        import numpy as np
        manifest, contract, objective, runtime, prepared, reference = base.authenticate(out)
        base.require(base.sha(base.MANIFEST) == base.MANIFEST_SHA, 'Frozen manifest mismatch')
        original = np.asarray(reference['b_grid'])
        P = copy.deepcopy(reference['parameters'])
        grid, inherited, numerical, embedding = credit.embed_credit_grid(base, P, reference, grid_path)
        base.require(len(grid) == 262 and len(original) == 160, 'Unexpected common-grid size')
        base.require(np.array_equal(inherited[np.asarray(embedding['old_indices'])], reference['stationary_g_pre']),
                     'Inherited occupied atoms changed')
        q0 = float(np.asarray(reference['solution'].p_eq)[0])
        model = prepared.rt['model']
        rows = []
        for factor in FACTORS:
            price = q0 * factor
            q = adapter.construct(model, P, grid, price)
            # Independent resource recurrence over all tenure alternatives.
            M = np.zeros((P.J, len(P.z_grid), len(q['cost'])))
            for j in range(P.J - 1, -1, -1):
                L = np.empty(len(q['cost']))
                for h in range(len(L)):
                    candidates = [-q['sale'][h]] if q['survival'][j] < 1 else []
                    if q['survival'][j] > 0:
                        candidates.append(M[j + 1, :, h].max())
                    L[h] = max(candidates)
                np.testing.assert_allclose(L, q['human'][j] - q['sale'], atol=2e-13, rtol=0)
                for z in range(len(P.z_grid)):
                    for old in range(len(L)):
                        M[j, z, old] = min(
                            (L[new] + q['oc'][new] - q['income'][j, z]) / P.R_gross
                            - (0 if new == old else q['sale'][old] - q['cost'][new])
                            for new in range(len(L)))
                np.testing.assert_allclose(M[j], q['minimum'][j, :, None] - q['sale'][None, :],
                                           atol=2e-13, rtol=0)
            # Original point masses are embedded exactly once, independent of price.
            support = q['pre'].transpose(3, 2, 0, 1)[:, :, None, :, :, None, None]
            bad = np.broadcast_to(~support, inherited.shape)
            inherited_bad = float(inherited[bad].sum())
            economic_ok = grid[:, None, None, None] > (
                q['minimum'].T[None, None, :, :] - q['sale'][None, :, None, None])
            economic_bad = np.broadcast_to(
                ~economic_ok.transpose(0, 1, 3, 2)[:, :, None, :, :, None, None], inherited.shape)
            economic_badmass = float(inherited[economic_bad].sum())
            entry = np.asarray(P.fixed_reference_entry_conditional)
            entry_bad = [float(entry[:, z][~q['pre'][0, z, 0]].sum()) for z in range(len(P.z_grid))]
            tightening = q['floors'] - (q['human'][:, None] - q['sale'][None, :])
            base.require(np.isfinite(tightening).all() and np.min(tightening) >= -1e-9,
                         'Numerical floor is below economic solvency floor')
            # Candidate-specific nodes are diagnostic only; they never enter a solve.
            candidate_grid = adapter.build_grid(model, reference['parameters'], original, price)
            if factor == 1.:
                base.require(np.array_equal(candidate_grid, grid), 'q0 fixed grid differs from authenticated natural grid')
            candidate = adapter.construct(model, P, candidate_grid, price)
            candidate_tightening = candidate['floors'] - (candidate['human'][:, None] - candidate['sale'][None, :])
            candidate_indices = np.searchsorted(candidate_grid, original)
            candidate_inherited = np.zeros((len(candidate_grid),) + inherited.shape[1:], dtype=inherited.dtype)
            candidate_inherited[candidate_indices] = reference['stationary_g_pre']
            candidate_support = candidate['pre'].transpose(3, 2, 0, 1)[:, :, None, :, :, None, None]
            candidate_bad = float(candidate_inherited[np.broadcast_to(~candidate_support, candidate_inherited.shape)].sum())
            extra_tightening = tightening - candidate_tightening
            worst = np.unravel_index(np.argmax(extra_tightening), extra_tightening.shape)
            inherited_wealth_equal_floor = 0.
            for j in range(P.J):
                for h in range(len(q['sale'])):
                    inherited_wealth_equal_floor += float(inherited[:, h, :, j][
                        np.abs(grid - q['floors'][j, h]) < 1e-8].sum())
            rows.append(dict(price_factor=factor, price=price, inherited_unsupported_mass=inherited_bad,
                             inherited_economically_insolvent_mass=economic_badmass,
                             entrant_unsupported_mass_by_income=entry_bad,
                             maximum_numerical_floor_tightening=float(np.max(tightening)),
                             minimum_numerical_floor_tightening=float(np.min(tightening)),
                             candidate_specific_grid_nodes=len(candidate_grid),
                             candidate_specific_maximum_floor_tightening=float(np.max(candidate_tightening)),
                             maximum_elementwise_extra_tightening_from_fixed_grid=float(extra_tightening[worst]),
                             worst_extra_tightening_age_index=int(worst[0]),
                             worst_extra_tightening_tenure_index=int(worst[1]),
                             elementwise_extra_tightening_by_age_tenure=extra_tightening.tolist(),
                             candidate_specific_inherited_unsupported_mass=candidate_bad,
                             fixed_minus_candidate_inherited_unsupported_mass=inherited_bad-candidate_bad,
                             maximum_economic_renter_debt=float(-np.min(q['human'])),
                             inherited_liquid_wealth_equal_saving_floor_mass=inherited_wealth_equal_floor,
                             support_ok=inherited_bad <= 2e-10 and economic_badmass <= 2e-10 and
                                        max(entry_bad) == 0.))
        ready = all(row['support_ok'] for row in rows)
        receipt = dict(status='ready' if ready else 'blocked_before_solve', model_solves=0,
            reference_label=base.LABEL, reference_manifest_sha256=base.MANIFEST_SHA,
            grid_sha256=base.sha(grid_path), grid_nodes=len(grid), original_grid_nodes=len(original),
            fixed_common_grid_all_prices=True, original_entry_and_inherited_atoms_exact=True,
            independent_general_recursion_all_prices_pass=True,
            source_sha256={p.name: base.sha(p) for p in sorted(source.glob('*.py'))},
            factors=list(FACTORS), rows=rows, elapsed_seconds=time.monotonic()-started,
            interpretation='Support screen only. No utility/policy convergence or market clearing certified.')
        base.write(out / 'preflight.json', receipt)
        base.require(ready, 'Common 262-node grid fails price support; stop before any solve')
        print(json.dumps(dict(status=receipt['status'], rows=rows)))
    except BaseException as exc:
        base.write(out / 'failure.json', dict(status='failed', error=str(exc), traceback=traceback.format_exc()))
        raise


if __name__ == '__main__':
    main()
