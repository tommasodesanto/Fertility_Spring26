"""Appendix D numbers and figure (Oct 24 2026): free entry eta 1.5 (data_fi_torch/P15fi_*, 200-period arms) and slow stock (KMV)
eta 0.63 at the 200-period fit (data_k200_torch/P063k_*). Saved converged runs only. % vs the arm's own BASE from the same 2023 state.
Long run = the stationary endpoint ('endpoint' of the result: population_scale and price), where births per household are at
replacement, so long-run births and households move by the same %. Completed fertility as in p6_figs.py.
Run: python p6_appendix.py <transitions dir> <out dir>"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
T = Path(sys.argv[1]); OUT = Path(sys.argv[2])
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True, 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
NAVY, RUST, TEAL = '#1f4e79', '#a84432', '#4f9a9a'
def load(f):
    if not f.exists(): return None
    d = json.load(open(f))
    if not d.get('converged'): return None
    L = d['last']; R = L['rows']; F = L['fertility']
    s = dict(ep=d['endpoint'], y=np.array([r['calendar_year'] for r in R]))
    for k, n in (('b', 'birth_children'), ('p', 'asset_price'), ('N', 'adult_population')): s[k] = np.array([r[n] for r in R])
    s['asfr'] = np.array([f['age_specific_birth_rate_topcode_adjusted'] for f in F]); s['tfr'] = np.array([f['period_tfr_topcode_adjusted'] for f in F])
    return s
def ccf(s):
    A = s['asfr']; return np.array([sum(A[k + j, j] for j in range(7)) for k in range(A.shape[0] - 7)])
SETS = {'fi15': (T / 'data_fi_torch', 'P15fi'), 'k063': (T / 'data_k200_torch', 'P063k')}
ARMS = ('PHI100', 'SUP110', 'SUPMC', 'PPHI100', 'FLAT', 'INS', 'PTAX', 'LR_BASE_RS')
vals = {}
for key, (D, pre) in SETS.items():
    B = load(D / f'{pre}_BASE.json'); cb = ccf(B); i63 = int(np.where(B['y'] == 2063)[0][0])
    v = dict(base=dict(tfr2023=float(B['tfr'][0]), tfr2063=float(B['tfr'][i63]), last_year=int(B['y'][-1])))
    X_ = {}
    for a in ARMS:
        X = load(D / f'{pre}_{a}.json')
        if X is None: continue
        X_[a] = X; c = ccf(X)
        v[a] = dict(file=f'{pre}_{a}.json', births2023=100 * (X['b'][0] / B['b'][0] - 1), births2063=100 * (X['b'][i63] / B['b'][i63] - 1),
                    births_lr=100 * (X['ep']['population_scale'] / B['ep']['population_scale'] - 1), cf2023=100 * (c[0] / cb[0] - 1),
                    price2023=100 * (X['p'][0] / B['p'][0] - 1), price2063=100 * (X['p'][i63] / B['p'][i63] - 1),
                    price_lr=100 * (X['ep']['price'] / B['ep']['price'] - 1))
    vals[key] = v
    if key == 'k063':   # paths figure, slow stock: credit, matched (or +10%) supply, parents, to 2103
        m = B['y'] <= 2103; yc = B['y'][:len(cb)]; mc = yc <= 2063
        sup = 'SUPMC' if 'SUPMC' in X_ else 'SUP110'
        arms = (('PHI100', 'Mortgage credit', RUST, 'o'), (sup, 'Housing supply' + ('' if sup == 'SUPMC' else ', +10%'), NAVY, 's'),
                ('PPHI100', 'Credit for parents', TEAL, '^'))
        fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.0))
        for (k, ttl), ax in zip((('b', '(a) Births'), ('cf', '(b) Completed fertility by cohort'), ('N', '(c) Households'), ('p', '(d) House price')), axs.ravel()):
            for a, lab, col, mk in arms:
                if a not in X_: continue
                X = X_[a]
                if k == 'cf': ax.plot(yc[mc], 100 * (ccf(X)[mc] / cb[mc] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
                else: ax.plot(X['y'][m], 100 * (X[k][m] / B[k][m] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
            ax.axhline(0, color='0.6', lw=0.6); ax.set_title(ttl, fontsize=10, loc='left'); ax.set_ylabel('% from no policy')
            ax.set_xticks([2023, 2043, 2063, 2083, 2103] if k != 'cf' else [2023, 2033, 2043, 2053, 2063]); ax.grid(axis='y', alpha=.25, lw=.5)
            ax.set_xlabel('Year the cohort is aged 18 to 21' if k == 'cf' else 'Year', fontsize=9 if k == 'cf' else 10)
        axs[0, 0].legend(frameon=False, fontsize=8, loc='upper left')
        fig.tight_layout(); fig.savefig(OUT / 'fig_appD_paths_k200_eta063.pdf'); plt.close(fig)
json.dump(vals, open(OUT / 'appD_values.json', 'w'), indent=1)
for key, v in vals.items():
    print('==', key, {k: round(x, 4) if isinstance(x, float) else x for k, x in v['base'].items()})
    for a, x in v.items():
        if a == 'base': continue
        print(f'  {a:10s}', ' '.join(f'{k} {x[k]:+.3f}' for k in ('births2023', 'births2063', 'births_lr', 'cf2023', 'price2023', 'price2063', 'price_lr')))
