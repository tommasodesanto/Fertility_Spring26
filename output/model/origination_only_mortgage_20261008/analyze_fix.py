"""Tables for the origination-only vs revolving-collateral comparison (fixed price, 14.402). Writes tables.md.

Definitions (all from each cell's own stationary distribution unless 'impact'):
  births            KFE births per period, sum(pre * realized birth count); first births: childless states only.
  impact            the 95% cell's birth policies on the SAME regime's 80% start-of-period distribution, % of that 80% cell's births.
  completed CEB     mean children ever born at 46-49 (3+ top-coded); childless share at 46-49.
  ownership         current (post-decision) tenure; rooms = mean occupied rooms (all households).
  debt              start-of-period owners (origin tenure owned): share with b < 0, mean debt -b among them, mean LTV -b/(pH).
  margin            childless households 22-33 (start of period), weighted by pre-mass x p(1-p), p = first-birth attempt probability.
  liquid resources  L = b - stayer floor: revolving b - min(b, -phi p H) for owners, origination-only max(b, 0); renters max(b, 0).
                    (what the household can spend next period without selling; in model units, per four-year period)
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); import run_cells as rc
K = np.arange(4)
CELLS4 = ['R80F', 'R95F', 'O80F', 'O95F']
P, grid, price, T, _ = rc.build('R80')
from model.engine.shared import income_at_state, income_transition_values
z, _, _ = income_transition_values(P); z = np.asarray(z, float); J = int(P.J)
Y = np.array([[income_at_state(P, 0, j, float(zz)) for zz in z] for j in range(J)])     # (J, Nz)
R = float(P.R_gross); H = np.asarray(P.H_own, float); b = np.asarray(grid, float).reshape(-1)
ages = 18 + 4 * np.arange(J)
MED22 = float(Y[1, 4])


def load(c):
    S = dict(np.load(HERE / 'cells' / c / 'solution_arrays.npz'))
    f = (S['birth_count_realized_probs'] * K).sum(-1); f[..., 3, :] = 0
    for n in range(4): f[..., n, n + 1:] = 0
    S['f'] = f; S['st'] = json.loads((HERE / 'cells' / c / 'scalar_stats.json').read_text())
    S['phi'] = json.loads((HERE / 'cells' / c / 'spec.json').read_text())['phi']; S['orig'] = c[0] == 'O'
    return S


def outcomes(S):
    pre, f, g = S['birth_count_pre_distribution'], S['f'], S['g']
    B = float((pre * f).sum()); FB = float((pre[..., 0, 0] * f[..., 0, 0]).sum())
    nm = pre[:, :, :, 7].sum(axis=(0, 1, 2, 3, 5)); ceb = float((nm * K).sum() / nm.sum()); cl = float(nm[0] / nm.sum())
    first = (pre[..., 0, 0] * f[..., 0, 0]).sum(axis=(0, 1, 2, 4)); afb = float((first * (ages + 2)).sum() / first.sum())
    t = g.sum(axis=(0, 2, 4, 5, 6))
    own = lambda js: float(t[1:, js].sum() / t[:, js].sum())
    o = dict(B=B, FB=FB, ceb=ceb, childless=cl, afb=afb, own=own(slice(None)), own1829=own(slice(0, 3)), own3055=own(slice(3, 10)),
             rooms=S['st']['aggregate_housing_demand'] / S['st']['total_mass'], fb_rooms=S['st']['housing_increment_0to1_eventstudy_t3'])
    # debt among start-of-period owners
    w = pre.sum(axis=(2, 4, 5, 6))                          # (b, ten, j)
    own_w = w[:, 1:, :]; value = price * H[None, :, None]
    d = np.maximum(-b, 0)[:, None, None]
    def debt(js):
        m = own_w[:, :, js]; dd = np.broadcast_to(d, own_w.shape)[:, :, js]; v = np.broadcast_to(value, own_w.shape)[:, :, js]
        has = m * (dd > 0)
        return float(has.sum() / m.sum()), float((has * dd).sum() / max(has.sum(), 1e-300)), float((has * dd / v).sum() / max(has.sum(), 1e-300))
    o['debt'] = {int(ages[j]): debt([j]) for j in [1, 2, 3, 4, 5, 7, 9, 11, 13]}
    o['debt_2545'] = debt([2, 3, 4, 5, 6])
    # fertility margin
    js = [1, 2, 3]
    p = S['birth_count_action_probs'][:, :, 0, 1:4, :, 0, 0, 1].astype(float)   # (b, ten, j, z)
    m = pre[:, :, 0, 1:4, :, 0, 0].astype(float)
    wt = m * p * (1 - p)
    B4 = b[:, None, None, None]; TEN = np.arange(1 + len(H))[None, :, None, None]
    fl = np.concatenate([[0.0], -S['phi'] * price * H])[None, :, None, None]
    if S['orig']: L = np.where(TEN > 0, np.maximum(B4, 0), np.maximum(B4, 0)) + 0 * fl
    else: L = np.where(TEN > 0, B4 - np.minimum(B4, fl), np.maximum(B4, 0))
    x = R * B4 + Y[js, :][None, None, :, :]
    avg = lambda a, ww=wt: float((ww * np.broadcast_to(a, ww.shape)).sum() / ww.sum())
    o['margin'] = dict(mass=float(wt.sum() / m.sum()), p=avg(p, m), own=avg((TEN > 0).astype(float)), debt=avg((B4 < -1e-9) & (TEN > 0)),
                       b=avg(B4), x=avg(x), L=avg(L), L_all=avg(L, m))
    return o


C = {c: load(c) for c in CELLS4}; O = {c: outcomes(S) for c, S in C.items()}
def impact(c80, c95):
    pre = C[c80]['birth_count_pre_distribution']
    return [100 * float((pre * (C[c95]['f'] - C[c80]['f'])).sum()) / O[c80]['B'],
            100 * float((pre[..., 0, 0] * (C[c95]['f'] - C[c80]['f'])[..., 0, 0]).sum()) / O[c80]['FB']]
IR, IO = impact('R80F', 'R95F'), impact('O80F', 'O95F')

def pct(a, b_): return f'{100 * (b_ / a - 1):+.2f}%'
def pp(a, b_): return f'{100 * (b_ - a):+.1f} pp'
rows = []
def row(label, key, fmt, diff):
    v = [O[c][key] for c in CELLS4]
    rows.append([label, fmt(v[0]), fmt(v[1]), diff(v[0], v[1]), fmt(v[2]), fmt(v[3]), diff(v[2], v[3])])
f3 = lambda v: f'{v:.3f}'; f4 = lambda v: f'{v:.4f}'; pc = lambda v: f'{100 * v:.1f}%'
rows.append(['Births on impact (95% policies, 80% distribution)', '', '', f'{IR[0]:+.2f}%', '', '', f'{IO[0]:+.2f}%'])
rows.append(['First births on impact', '', '', f'{IR[1]:+.2f}%', '', '', f'{IO[1]:+.2f}%'])
row('Births per period, stationary', 'B', f4, pct)
row('First births per period, stationary', 'FB', f4, pct)
row('Completed fertility (CEB 46-49)', 'ceb', f3, pct)
row('Childless at 46-49', 'childless', pc, pp)
row('Mean age at first birth', 'afb', lambda v: f'{v:.2f}', lambda a, b_: f'{b_ - a:+.2f} y')
row('Ownership, all', 'own', pc, pp)
row('Ownership 18-29', 'own1829', pc, pp)
row('Ownership 30-55', 'own3055', pc, pp)
row('Mean rooms, all households', 'rooms', f3, pct)
row('First-birth room response (-1 to +3)', 'fb_rooms', f3, lambda a, b_: f'{b_ - a:+.3f}')
for lab, k, fmt, diff in [('Owners 26-45 with debt', 0, pc, pp), ('Mean debt of indebted owners 26-45', 1, f3, pct), ('Mean LTV of indebted owners 26-45', 2, pc, pp)]:
    v = [O[c]['debt_2545'][k] for c in CELLS4]; rows.append([lab, fmt(v[0]), fmt(v[1]), diff(v[0], v[1]), fmt(v[2]), fmt(v[3]), diff(v[2], v[3])])
for lab, k, fmt, diff in [('Margin: owner share', 'own', pc, pp), ('Margin: share owing a mortgage', 'debt', pc, pp), ('Margin: mean b', 'b', f3, lambda a, b_: f'{b_ - a:+.3f}'),
                          ('Margin: cash on hand R b + y', 'x', f3, pct), ('Margin: liquid resources L', 'L', f3, pct), ('All childless 22-33: liquid resources L', 'L_all', f3, pct),
                          ('Margin: mean attempt probability (childless 22-33)', 'p', f4, pct)]:
    v = [O[c]['margin'][k] for c in CELLS4]; rows.append([lab, fmt(v[0]), fmt(v[1]), diff(v[0], v[1]), fmt(v[2]), fmt(v[3]), diff(v[2], v[3])])

def md(rows, h): return '| ' + ' | '.join(h) + ' |\n|' + '---|' * len(h) + '\n' + ''.join('| ' + ' | '.join(str(x) for x in r) + ' |\n' for r in rows)
out = [f'Fixed price {price:.6f}, rebate T {T:.5f} held, base 14.402. Median age-22 four-year earnings = {MED22:.3f}. R = {R:.4f} per period on both signs of b.\n',
       md(rows, ['', 'Revolving 80%', 'Revolving 95%', 'LTV effect', 'Orig.-only 80%', 'Orig.-only 95%', 'LTV effect']),
       '\n### Debt by age, start-of-period owners: share with debt / mean debt / mean LTV\n',
       md([[a] + [f'{100 * O[c]["debt"][a][0]:.0f}% / {O[c]["debt"][a][1]:.2f} / {100 * O[c]["debt"][a][2]:.0f}%' for c in CELLS4] for a in O['R80F']['debt']],
          ['age', 'Revolving 80%', 'Revolving 95%', 'Orig.-only 80%', 'Orig.-only 95%'])]
(HERE / 'tables_fix.md').write_text('\n'.join(out)); print('\n'.join(out))
json.dump({c: O[c] for c in O}, open(HERE / 'outcomes_fix.json', 'w'), indent=1, default=float)

# regime effect at 80% and old-age ownership
pre = C['R80F']['birth_count_pre_distribution']
reg = [100 * float((pre * (C['O80F']['f'] - C['R80F']['f'])).sum()) / O['R80F']['B'],
       100 * float((pre[..., 0, 0] * (C['O80F']['f'] - C['R80F']['f'])[..., 0, 0]).sum()) / O['R80F']['FB']]
extra = ['\n### Regime effect at 80% (origination-only vs revolving)\n',
         f'Births on impact (O80 policies on the R80 distribution): {reg[0]:+.2f}%; first births {reg[1]:+.2f}%. '
         f'Stationary births {pct(O["R80F"]["B"], O["O80F"]["B"])}, completed fertility {pct(O["R80F"]["ceb"], O["O80F"]["ceb"])}.',
         '\n| | ' + ' | '.join(CELLS4) + ' |\n|---|---|---|---|---|\n| Ownership 65-75 | ' + ' | '.join(f'{100 * C[c]["st"]["old_age_own_rate_6575"]:.1f}%' for c in CELLS4) + ' |\n'
         '| Wealth / annual earnings | ' + ' | '.join(f'{C[c]["st"]["aggregate_wealth_to_annual_gross_labor_earnings"]:.2f}' for c in CELLS4) + ' |\n'
         '| Young owners 25-34 with b < 0 | ' + ' | '.join(f'{100 * C[c]["st"]["owner_neg_liquid_share_2534"]:.1f}%' for c in CELLS4) + ' |\n']
(HERE / 'tables_fix.md').write_text((HERE / 'tables_fix.md').read_text() + '\n'.join(extra)); print('\n'.join(extra))
