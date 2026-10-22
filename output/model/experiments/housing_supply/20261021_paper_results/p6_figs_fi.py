"""Section 6 numbers and figures from the FREE-ENTRY eta 0.63 policy arms (Oct 22 2026; main text since his decision Oct 21 22:36Z).
Saved Torch results only, no solve: 2026-10-18_accounting_recal/overnight/pol_063fi/out/runs/P063fi_<ARM>/result.json, copied to
20261017_accounting/transitions/data_fi_torch/P063fi_<ARM>.json. Only converged arms are used; an arm with no file is skipped and
stays a placeholder in the text. The +10% supply arm is the converged run SUP110p. Matched supply: SUPMC when it exists.
Same definitions as p6_figs.py / p6_figs_oldsetup.py: % against P063fi_BASE from the same 2023 state; completed fertility by cohort =
sum of age-specific rates along the diagonal from the cohort aged 18-21 in the year shown (7 cells); long run = last date (2403).
Run: python p6_figs_fi.py <data dir> <out dir>"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = Path(sys.argv[1]); OUT = Path(sys.argv[2]); ETA = 0.63
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True, 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
NAVY, RUST, TEAL, GOLD, GREY = '#1f4e79', '#a84432', '#4f9a9a', '#c08a1e', '#6e6e6e'
FILE = {'BASE': 'BASE', 'PHI100': 'PHI100', 'SUP110': 'SUP110p', 'SUPMC': 'SUPMC', 'PPHI100': 'PPHI100', 'FLAT': 'FLAT', 'INS': 'INS',
        'PTAX': 'PTAX'}
def load(arm):
    f = D / f'P063fi_{FILE[arm]}.json'
    if not f.exists(): return None
    d = json.load(open(f))
    if not d.get('converged'): return None
    L = d['last']; R = L['rows']; F = L['fertility']
    s = dict(y=np.array([r['calendar_year'] for r in R]))
    for k, n in (('b', 'birth_children'), ('p', 'asset_price'), ('r', 'renter_price'), ('N', 'adult_population'), ('H', 'housing_supply')):
        s[k] = np.array([r[n] for r in R])
    s['tfr'] = np.array([f['period_tfr_topcode_adjusted'] for f in F]); s['asfr'] = np.array([f['age_specific_birth_rate_topcode_adjusted'] for f in F])
    s['fb'] = np.array([f['first_birth_rate_all_households'] for f in F]); s['mid'] = np.array(F[0]['age_cell_midpoint'], float)
    s['o1829'] = np.array([st['own_18_29'] for st in L['stats']])
    return s
def ccf(s):
    A = s['asfr']; return np.array([sum(A[k + j, j] for j in range(7)) for k in range(A.shape[0] - 7)])
def cohort(s):
    a = s['asfr'][np.arange(7), np.arange(7)]; fb = s['fb'][np.arange(7), np.arange(7)]; mid = s['mid'][:7]
    return dict(cf=float(a.sum()), ever_first=float(fb.sum()), mean_age_first=float((fb * mid).sum() / fb.sum()))
YRS = (2023, 2063, 2103, 2403)
B = load('BASE'); assert B is not None
ix = [int(np.where(B['y'] == y)[0][0]) for y in YRS]; cb = ccf(B)
vals = dict(base=dict(tfr=[float(B['tfr'][i]) for i in ix], cf2023=float(cb[0]), price=[float(B['p'][i]) for i in ix]),
            cohort2023=dict(BASE=cohort(B)), arms_used={})
X_ = {}
for a in ('PHI100', 'SUP110', 'SUPMC', 'PPHI100', 'FLAT', 'INS', 'PTAX'):
    X = load(a)
    if X is None: continue
    X_[a] = X; c = ccf(X); i0 = ix[0]; vals['arms_used'][a] = f'P063fi_{FILE[a]}.json'
    vals[a] = dict(births=[100 * (X['b'][i] / B['b'][i] - 1) for i in ix], cf2023=100 * (c[0] / cb[0] - 1),
                   price=[100 * (X['p'][i] / B['p'][i] - 1) for i in ix], rent=[100 * (X['r'][i] / B['r'][i] - 1) for i in ix],
                   households=[100 * (X['N'][i] / B['N'][i] - 1) for i in ix], own1829_pp=[100 * (X['o1829'][i] - B['o1829'][i]) for i in ix],
                   supply_multiplier=float((X['H'][i0] / X['p'][i0] ** ETA) / (B['H'][i0] / B['p'][i0] ** ETA)))
    if a in ('PHI100', 'SUPMC', 'PPHI100'): vals['cohort2023'][a] = cohort(X)
json.dump(vals, open(OUT / 'p6_values_fi.json', 'w'), indent=1)
m = B['y'] <= 2103; yc = B['y'][:len(cb)]; mc = yc <= 2063
def panel(arms, keys, fname, figsize, ncols):
    fig, axs = plt.subplots(len(keys) // ncols, ncols, figsize=figsize); axs = np.atleast_1d(axs).ravel()
    for (k, ttl), ax in zip(keys, axs):
        for a, lab, col, mk in arms:
            if a not in X_: continue
            X = X_[a]
            if k == 'cf': ax.plot(yc[mc], 100 * (ccf(X)[mc] / cb[mc] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
            else: ax.plot(X['y'][m], 100 * (X[k][m] / B[k][m] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
        ax.axhline(0, color='0.6', lw=0.6); ax.set_title(ttl, fontsize=10, loc='left'); ax.set_ylabel('% from no policy')
        ax.set_xticks([2023, 2043, 2063, 2083, 2103] if k != 'cf' else [2023, 2033, 2043, 2053, 2063]); ax.grid(axis='y', alpha=.25, lw=.5)
        ax.set_xlabel('Year the cohort is aged 18 to 21' if k == 'cf' else 'Year', fontsize=9 if k == 'cf' else 10)
    axs[0].legend(frameon=False, fontsize=8, loc='upper left')
    fig.tight_layout(); fig.savefig(OUT / fname); plt.close(fig)
panel((('PHI100', 'Mortgage credit', RUST, 'o'), ('SUPMC', 'Matched supply', NAVY, 's'), ('PPHI100', 'Credit for parents', TEAL, '^')),
      (('b', '(a) Births'), ('cf', '(b) Completed fertility by cohort'), ('N', '(c) Households'), ('p', '(d) House price')),
      'fig_p6_paths_fi_eta063.pdf', (7.4, 5.0), 2)
panel((('FLAT', 'Cash transfer', GOLD, 'o'), ('INS', 'Earnings insurance', TEAL, '^'), ('PTAX', 'Property tax, rebated', GREY, 'd')),
      (('b', '(a) Births'), ('p', '(b) House price')), 'fig_p6_transfers_fi_eta063.pdf', (7.4, 2.8), 2)
print('arms', vals['arms_used']); print('base', {k: np.round(x, 4) for k, x in vals['base'].items()})
print('cohort2023', {a: {k: round(x, 3) for k, x in v.items()} for a, v in vals['cohort2023'].items()})
for a in X_:
    v = vals[a]
    print(f'{a:8s}', 'births', np.round(v['births'], 2), 'cf', round(v['cf2023'], 2), 'price', np.round(v['price'], 2), 'rent', np.round(v['rent'], 1),
          'hh', np.round(v['households'], 1), 'own1829', np.round(v['own1829_pp'], 1), 'mult', round(v['supply_multiplier'], 4))
