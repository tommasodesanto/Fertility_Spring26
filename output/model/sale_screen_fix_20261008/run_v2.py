"""COPY for the Oct 8 sale-screen fix (engine root tmp/sale_screen_fix_20261008; SALE_SCREEN_LEGACY=1 restores the old screen;
CASES_DIR names the output folder). Original: Full chain-11-style set at the Mac round-3 best 14.402 (chain 0), rebate on in every run. Writes cases_v2/<case>/.

Per case: solve (fixed price: solve_balanced_at_price; ge_*: solve_stationary_ge, closure fixed_h0), saved arrays,
slim executed_P.json (chain-11 format), 17 standard diagnostics, native observation (fixed price: renewal tolerance
relaxed, diagnostic) and the new-contract 14-row fit. Headline cases also get the frontend packet in packet/:
native_result.npz + metadata.json (save_case), 8 policy and 7 aggregate plots, explorer arrays and config.
"""
import argparse, hashlib, json, os, shutil, time, traceback
from pathlib import Path
import numpy as np
from common import *

PACKET = set()   # no frontend packets in this verification copy


def write(path, obj):
    def enc(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.floating, np.integer)): return o.item()
        if isinstance(o, (set, tuple)): return list(o)
        return str(o)
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, default=enc) + '\n')


def slim_P(Q):
    out = {}
    for k, v in vars(Q).items():
        if isinstance(v, np.ndarray) and v.size > 5000 and not v.dtype.hasobject:
            out[k] = dict(shape=list(v.shape), dtype=str(v.dtype), sha256=hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest())
        else:
            out[k] = v
    return out


def packet(dest, case, sol, Q, grid, price, params):
    from model.storage import StoredResult, save_case
    from model import workflow
    pk = dest / 'packet'; pk.mkdir()
    res = StoredResult(sol, Q, grid, price, parameters=params, label=f'Mac round-3 best 14.402 (chain 0), rebate on, {case}')
    save_case(res, pk, metadata=dict(case=case, source='Mac round-3 best 14.402 (chain 0)', rebate='balanced, on',
                                     mode='GE fixed_h0' if case.startswith('ge_') else 'fixed price'))
    shutil.copytree(dest / 'standard_diagnostics', pk / 'standard_diagnostics')
    workflow._cached_plots(pk)
    workflow._write_explorer_assets(pk, res)
    cfg = json.loads((pk / 'explorer_cases.json').read_text())
    cfg['cases'][0].update(id=case, label=res.label)
    (pk / 'explorer_cases.json').write_text(json.dumps(cfg, indent=2) + '\n')


def run(case):
    from model.equilibrium import solve_balanced_at_price, solve_stationary_ge
    from model.reporting import build_context
    from model.engine.diagnostics import write_diagnostics
    from model.estate_contract import rescore_report
    from model import native_phase_b
    dest = HERE / os.environ.get('CASES_DIR', 'cases') / case; dest.mkdir(parents=True, exist_ok=False)
    P, grid, price, spec, P0, grid0 = build_inputs(case)
    spec['changed_fields'] = changed_fields(P0, P); spec['entry_mean_ratio'] = entry_mean(P, grid) / entry_mean(P0, grid0)
    write(dest / 'spec.json', spec)
    t0 = time.monotonic()
    if case.startswith('ge_'):
        r = solve_stationary_ge(P, grid, out=dest / 'ge', price_start=price, budget_seconds=3000, max_lifecycle=32, closure='fixed_h0')
        sol, Q, sd, price = r['solution'], r['P'], r['shared'], float(r['price'])
        rep = Path(r['report_directory']); rescore_report(rep)
        shutil.copytree(rep / 'standard_diagnostics', dest / 'standard_diagnostics')
        write(dest / 'closure.json', dict(r['closure'], rebate=r.get('property_tax_rebate'), price=price,
                                          lifecycle_solves=r['lifecycle_solves'], price_trials=r['price_trials']))
    else:
        (dest / 'reporting').mkdir()
        ctx = build_context(P, grid, dest / 'reporting', price_start=price, deadline=time.time() + 1500, max_lifecycle=32, closure='fixed_h0')
        out = solve_balanced_at_price(P, grid, price, start=spec['rebate_start'], seconds_per_solve=600)
        sol, Q, sd, rec = out['solution'], out['P'], out['shared'], out['property_tax_rebate']
        write(dest / 'fiscal_record.json', rec)
        write_diagnostics(sol, Q, dest / 'standard_diagnostics')
        native_phase_b.RENEWAL_TOL = 1e9   # diagnostic: a fixed price is not a renewal root
        try:
            obs = native_phase_b.observe_price(ctx, dict(P=Q, b_grid=np.asarray(grid), sd=sd, sol=sol, price=np.asarray([price]),
                                                         case_deadline_epoch=time.time() + 600, fiscal=rec), case, final=True)
            write(dest / 'native_observation.json', obs)
            rescore_report(dest / 'reporting/phase_b_ge' / case)
        except Exception as exc:
            write(dest / 'observation_failure.json', dict(error=repr(exc), traceback=traceback.format_exc()))
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and not v.dtype.hasobject}
    arrays.update({'shared.' + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and not v.dtype.hasobject})
    np.savez_compressed(dest / 'solution_arrays.npz', **arrays)
    write(dest / 'executed_P.json', slim_P(Q))
    write(dest / 'solve_completed.json', dict(case=case, price=price, seconds=time.monotonic() - t0,
                                              transfer=float(Q.property_tax_lump_sum_transfer), phi=np.asarray(Q.phi).tolist()))
    print(case, 'solved', round(time.monotonic() - t0, 1), 's', flush=True)
    if case in PACKET:
        try:
            packet(dest, case, sol, Q, np.asarray(grid), price, spec.get('params', point()[0]))
        except Exception as exc:
            write(dest / 'packet_failure.json', dict(error=repr(exc), traceback=traceback.format_exc()))
            print(case, 'packet failed', repr(exc)[:300], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--cases', required=True); a = ap.parse_args()
    for c in a.cases.split(','):
        assert c in CASES, c
        try: run(c)
        except Exception as exc:
            d = HERE / os.environ.get('CASES_DIR', 'cases') / c; d.mkdir(parents=True, exist_ok=True)
            write(d / 'FAILED.json', dict(error=repr(exc), traceback=traceback.format_exc()))
            print('CASE FAILED', c, repr(exc)[:300], flush=True)
