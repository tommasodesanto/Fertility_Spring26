#!/usr/bin/env python3
"""Zero-solve union-grid support and independent solvency check on Torch."""
import argparse, copy, json, os, sys, time, traceback
from pathlib import Path
import run_credit as credit
import run_fixed_price as base
import natural_credit as adapter

FACTORS = (.98, .99, 1., 1.01, 1.02)

def independent(q, P):
    import numpy as np
    M = np.zeros((P.J, len(P.z_grid), len(q['cost'])))
    for j in range(P.J-1, -1, -1):
        L = np.empty(len(q['cost']))
        for h in range(len(L)):
            candidates = [-q['sale'][h]] if q['survival'][j] < 1 else []
            if q['survival'][j] > 0: candidates.append(M[j+1, :, h].max())
            L[h] = max(candidates)
        np.testing.assert_allclose(L, q['human'][j]-q['sale'], atol=2e-13, rtol=0)
        for z in range(len(P.z_grid)):
            for old in range(len(L)):
                M[j,z,old] = min((L[new]+q['oc'][new]-q['income'][j,z])/P.R_gross
                    -(0 if new == old else q['sale'][old]-q['cost'][new]) for new in range(len(L)))
        np.testing.assert_allclose(M[j], q['minimum'][j,:,None]-q['sale'][None,:], atol=2e-13, rtol=0)

def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path, required=True)
    out = parser.parse_args().output.resolve()
    base.require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    base.require(not out.exists(), 'Preflight output exists'); out.mkdir(parents=True)
    start = time.monotonic(); source = Path(__file__).resolve().parent
    try:
        import numpy as np
        manifest, contract, objective, runtime, prepared, ref = base.authenticate(out)
        base.require(base.sha(base.MANIFEST) == base.MANIFEST_SHA, 'Frozen manifest differs')
        original = np.asarray(ref['b_grid']); original_P = ref['parameters']
        model = prepared.rt['model']; q0 = float(np.asarray(ref['solution'].p_eq)[0])
        candidate = {f: adapter.build_grid(model, original_P, original, q0*f) for f in FACTORS}
        grid = np.unique(np.concatenate([original]+[candidate[f] for f in FACTORS]))
        indices = np.searchsorted(grid, original)
        base.require(np.array_equal(grid[indices], original) and len(grid)>len(original),
                     'Union grid loses original atoms')
        grid_path = out/'common_grid.json'
        base.write(grid_path, dict(grid=grid.tolist(),old_indices=indices.tolist(),
            method='Exact union of original and all five price-specific natural-credit grids',
            factors=list(FACTORS),reference_label=base.LABEL,original_grid_preserved_exactly=True))
        P=copy.deepcopy(original_P)
        common, inherited, numerical, embedding=credit.embed_credit_grid(base,P,ref,grid_path)
        added=np.ones(len(grid),dtype=bool); added[indices]=False
        base.require(np.array_equal(common,grid) and np.array_equal(inherited[indices],ref['stationary_g_pre'])
                     and not np.any(inherited[added]),'Inherited atom embedding changed')
        rows=[]
        for f in FACTORS:
            price=q0*f; q=adapter.construct(model,P,grid,price); independent(q,P)
            support=q['pre'].transpose(3,2,0,1)[:,:,None,:,:,None,None]
            bad=float(inherited[np.broadcast_to(~support,inherited.shape)].sum())
            econ=grid[:,None,None,None]>(q['minimum'].T[None,None,:,:]-q['sale'][None,:,None,None])
            ebad=float(inherited[np.broadcast_to(
                ~econ.transpose(0,1,3,2)[:,:,None,:,:,None,None],inherited.shape)].sum())
            entry=np.asarray(P.fixed_reference_entry_conditional)
            entry_bad=[float(entry[:,z][~q['pre'][0,z,0]].sum()) for z in range(len(P.z_grid))]
            tight=q['floors']-(q['human'][:,None]-q['sale'][None,:])
            qc=adapter.construct(model,original_P,candidate[f],price); independent(qc,original_P)
            ctight=qc['floors']-(qc['human'][:,None]-qc['sale'][None,:])
            extra=tight-ctight; worst=np.unravel_index(np.argmax(extra),extra.shape)
            base.require(np.isfinite(tight).all() and np.min(tight)>=-1e-9,'Numerical floor below economic floor')
            fm=float(np.max(tight)); cm=float(np.max(ctight))
            rows.append(dict(price_factor=f,price=price,candidate_specific_grid_nodes=len(candidate[f]),
                union_grid_nodes=len(grid),inherited_unsupported_mass=bad,
                inherited_economically_insolvent_mass=ebad,entrant_unsupported_mass_by_income=entry_bad,
                maximum_union_numerical_floor_tightening=fm,
                maximum_candidate_specific_floor_tightening=cm,
                maximum_elementwise_extra_tightening=float(extra[worst]),
                worst_extra_tightening_age_index=int(worst[0]),
                worst_extra_tightening_tenure_index=int(worst[1]),
                floor_tightening_gate=fm<=cm+1e-8,
                support_gate=bad<=2e-10 and ebad<=2e-10 and max(entry_bad)==0.))
        ready=all(r['floor_tightening_gate'] and r['support_gate'] for r in rows)
        receipt=dict(status='ready' if ready else 'blocked_before_solve',model_solves=0,
            reference_label=base.LABEL,reference_manifest_sha256=base.MANIFEST_SHA,
            common_grid_sha256=base.sha(grid_path),common_grid_nodes=len(grid),
            original_grid_nodes=len(original),original_entry_and_inherited_atoms_exact=True,
            independent_general_recursion_all_prices_pass=True,fixed_common_union_grid_all_prices=True,
            source_sha256={x.name:base.sha(x) for x in sorted(source.glob('*.py'))},
            factors=list(FACTORS),rows=rows,elapsed_seconds=time.monotonic()-start,
            interpretation='No solve; support/floor check only, not policy convergence or market clearing.')
        base.write(out/'preflight.json',receipt)
        base.require(ready,'Common union grid fails support or maximum floor parity; no solve')
        print(json.dumps(dict(status=receipt['status'],grid_nodes=len(grid),rows=rows)))
    except BaseException as exc:
        base.write(out/'failure.json',dict(status='failed',error=str(exc),traceback=traceback.format_exc()))
        raise

if __name__=='__main__': main()
