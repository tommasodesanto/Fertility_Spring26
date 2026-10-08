"""Sale-screen fix verification at the saved 14.402 parameters (fixed price, balanced rebate re-solved, as the saved case).
Compares legacy (old screen) and fix against the saved packet: solution arrays bitwise, the 14-row target table, residuals."""
import json, csv
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
SAVED = HERE.parent / 'jmp_draft_deck_20261005/refit_best_14p40_mac_r3_chain0/cases_v2'
def arrays_equal(a, b):
    A, B = np.load(a), np.load(b); keys = sorted(set(A.files) & set(B.files))
    diff = [k for k in keys if not np.array_equal(A[k], B[k], equal_nan=True)]
    return len(keys), diff
def fit(d, case):
    rows = list(csv.DictReader(open(d / 'reporting/phase_b_ge' / case / 'target_fit_new_contract.csv')))
    return rows
def obs(d):
    o = json.loads((d / 'native_observation.json').read_text()); r = json.loads((d / 'fiscal_record.json').read_text())
    return dict(renewal_residual=o.get('renewal_residual'), absolute_housing_residual=o.get('absolute_housing_residual'),
                rebate_T=r.get('transfer'), rebate_residual=r.get('relative_residual', r.get('residual')))
out = []
for case in ['baseline', 'phi095']:
    l, f = HERE / 'cases_legacy' / case, HERE / 'cases_fix' / case
    s = SAVED / case if (SAVED / case).exists() else l   # phi095 has no saved packet: reference = legacy rerun
    out.append(f'reference for {case}: ' + ('saved packet' if s != l else 'legacy rerun (no saved packet)'))
    n, d = arrays_equal(l / 'solution_arrays.npz', s / 'solution_arrays.npz')
    out.append(f'## {case}\n\nLegacy screen vs saved packet: {n} solution arrays compared, {len(d)} differ {d[:8]}.')
    n2, d2 = arrays_equal(f / 'solution_arrays.npz', s / 'solution_arrays.npz')
    out.append(f'Fix vs saved packet: {len(d2)} of {n2} arrays differ.\n')
    for lab, dd in [('saved', s), ('legacy', l), ('fix', f)]:
        out.append(f'- {lab} residuals: {json.dumps(obs(dd), default=float)}')
    S, L, F = fit(s, case), fit(l, case), fit(f, case)
    keys = list(S[0].keys())
    out.append('\nkeys: ' + ', '.join(keys) + '\n')
    hdr = ['moment', 'target', 'model saved', 'model legacy', 'model fix', 'fix - saved', 'weight', 'loss saved', 'loss fix']
    out.append('| ' + ' | '.join(hdr) + ' |\n|' + '---|' * len(hdr))
    def g(r, *names):
        for nm in names:
            if nm in r and r[nm] not in ('', None): return r[nm]
        return ''
    tot_s = tot_f = 0.0
    for a, b_, c in zip(S, L, F):
        mid = g(a, 'moment', 'name', 'row', 'key')
        ms, ml, mf = (float(g(r, 'model', 'model_value') or 'nan') for r in (a, b_, c))
        ls_, lf = (float(g(r, 'loss_contribution', 'loss', 'contribution') or 0) for r in (a, c))
        tot_s += ls_; tot_f += lf
        out.append(f"| {mid} | {g(a, 'target', 'target_value')} | {ms:.6g} | {ml:.6g} | {mf:.6g} | {mf - ms:+.3g} | {g(a, 'weight')} | {ls_:.4g} | {lf:.4g} |")
    out.append(f'\nTotal loss saved {tot_s:.6f}, fix {tot_f:.6f}.\n')
txt = '\n'.join(out); (HERE / 'baseline_verification.md').write_text(txt); print(txt)
