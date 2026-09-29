"""Zero-solve Torch preflight for the stationary-price natural-credit adapter.

This only authenticates the frozen export and evaluates algebra/Boolean support.
It deliberately never calls a lifecycle or equilibrium solver.
"""
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import sys
import time


LABEL = '2007 stationary reference — block0506, September 28 verified export'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
PRICE_FACTORS = (1., 1.05, 1.35)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def check_case(adapter, P, original_grid, inherited, entry, model, price):
    import numpy as np
    grid = adapter.build_grid(model, P, original_grid, price)
    indices = np.searchsorted(grid, original_grid)
    np.testing.assert_array_equal(grid[indices], original_grid)
    q = adapter.construct(model, P, grid, price)

    # Independent general minimum-over-tenures recurrence, not the adapter's
    # Boolean construction.  It permits every tenure choice and hence checks
    # immediate sale plus renting as the resource-minimum fallback.
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

    shape = q['pre'].transpose(3, 2, 0, 1)[:, :, None, :, :, None, None]
    g = np.zeros((len(grid),) + inherited.shape[1:])
    g[indices] = inherited
    badmass = float(g[np.broadcast_to(~shape, g.shape)].sum())
    expanded_entry = np.zeros((len(grid), entry.shape[1]))
    expanded_entry[indices] = entry
    conditional_bad = [float(expanded_entry[:, z][~q['pre'][0, z, 0]].sum())
                       for z in range(len(P.z_grid))]
    economic_ok = grid[:, None, None, None] > (
        q['minimum'].T[None, None, :, :] - q['sale'][None, :, None, None])
    economic_bad = np.broadcast_to(~economic_ok.transpose(0, 1, 3, 2)[:, :, None, :, :, None, None],
                                   g.shape)
    economic_badmass = float(g[economic_bad].sum())
    assert badmass <= 2e-10 and max(conditional_bad) == 0 and economic_badmass <= 2e-10
    return dict(price=float(price), grid=grid, old_indices=indices,
                nodes=len(grid), minimum=float(grid[0]), maximum=float(grid[-1]),
                independent_general_recursion_pass=True,
                inherited_mass_outside_discrete_solvency=badmass,
                inherited_mass_outside_economic_solvency=economic_badmass,
                infeasible_entry_conditional_mass_by_income=conditional_bad)


def main():
    assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
    import numpy as np
    start = time.monotonic()
    source = Path(__file__).parent
    out = Path(sys.argv[1]); out.mkdir(exist_ok=False)
    base = load(source / 'run_fixed_price.py', 'credit_ge_preflight_base')
    adapter = load(source / 'natural_credit.py', 'credit_ge_preflight_adapter')
    assert base.LABEL == LABEL and base.MANIFEST_SHA == MANIFEST_SHA
    assert adapter.API_VERSION == 'block0506_credit_adapter_v1'
    assert adapter.MODE == 'natural_solvency_boolean_grid_stationary_ge_v1'
    manifest, contract, objective, runtime, prepared, ref = base.authenticate(out)
    P = ref['parameters']
    original_grid = np.asarray(ref['b_grid'])
    inherited = np.asarray(ref['stationary_g_pre'])
    entry = np.asarray(P.fixed_reference_entry_conditional)
    model = prepared.rt['model']
    q0 = float(ref['solution'].p_eq[0])
    cases = [check_case(adapter, P, original_grid, inherited, entry, model, q0 * factor)
             for factor in PRICE_FACTORS]
    fixtures = [(-.1, False), (0., False), (.25, False), (.5, True), (.75, True),
                (1., True), (1.1, False), (.5, True), (.6, True), (.4, False),
                (.001, False), (.999, True)]
    for x, want in fixtures:
        got = adapter.reachable(np.array([0., .5, 1.]), np.array([False, True, True]),
                                np.array([x]))[0]
        assert bool(got) == want
    pinned_grid = json.loads((source / 'credit_grid.json').read_text())
    assert cases[0]['nodes'] == 262
    np.testing.assert_array_equal(cases[0]['grid'], np.asarray(pinned_grid['grid']))
    assert cases[0]['old_indices'].tolist() == pinned_grid['old_indices']
    for case in cases:
        del case['grid']; del case['old_indices']
    result = dict(reference_label=LABEL, status='ready', model_solves=0,
        slurm_job=os.environ['SLURM_JOB_ID'], reference_manifest_sha256=MANIFEST_SHA,
        price_factors=list(PRICE_FACTORS), boolean_support_fixtures=len(fixtures),
        source_sha256={path.name: sha(path) for path in sorted(source.glob('*')) if path.is_file()},
        q0_grid_reproduces_prior_262_nodes_exactly=True, cases=cases,
        elapsed_seconds=time.monotonic() - start)
    (out / 'preflight_ge.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
