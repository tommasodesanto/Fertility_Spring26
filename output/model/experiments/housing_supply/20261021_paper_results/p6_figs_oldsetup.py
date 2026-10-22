"""Section 6 placeholder numbers and figure, Oct 22 2026 (his word 03:52Z: fill the free-entry cells with the slides' runs until
the free-entry runs land). OLD SETUP: instant supply at eta 0.63, old accounting, 7.9 calibration, decay history
2026-10-12_psi_decay_7p9; saved runs 20261017_accounting/transitions/data/prod_<ARM>.json (the JMP slides' policy runs). No solve.
Same definitions as p6_figs.py: % against prod_BASE from the same 2023 state; completed fertility by cohort = sum of age-specific
rates along the diagonal from the cohort aged 18-21 in the year shown (7 cells). Supply multiplier of an arm = (H/P^eta) vs BASE in 2023.
Run: python p6_figs_oldsetup.py <data dir> <out dir>"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = Path(sys.argv[1]); OUT = Path(sys.argv[2]); ETA = 0.63
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True, 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
NAVY, RUST, TEAL = '#1f4e79', '#a84432', '#4f9a9a'
def load(arm):
    d = json.load(open(D / f'prod_{arm}.json')); assert d['converged']; L = d['last']; R = L['rows']; F = L['fertility']
    s = dict(y=np.array([r['calendar_year'] for r in R]))
    for k, n in (('b', 'birth_children'), ('p', 'asset_price'), ('r', 'renter_price'), ('N', 'adult_population'), ('own', 'owner_rate'),
                 ('H', 'housing_supply')):
        s[k] = np.array([r[n] for r in R])
    s['tfr'] = np.array([f['period_tfr_topcode_adjusted'] for f in F]); s['asfr'] = np.array([f['age_specific_birth_rate_topcode_adjusted'] for f in F])
    s['o1829'] = np.array([st['own_18_29'] for st in L['stats']])
    return s
def ccf(s):
    A = s['asfr']; return np.array([sum(A[k + j, j] for j in range(7)) for k in range(A.shape[0] - 7)])
YRS = (2023, 2063, 2103, 2403)   # 2403 = last date (long run)
B = load('BASE'); ix = [int(np.where(B['y'] == y)[0][0]) for y in YRS]; cb = ccf(B)
vals = dict(base=dict(tfr=[float(B['tfr'][i]) for i in ix], cf2023=float(cb[0]), price=[float(B['p'][i]) for i in ix]))
for a in ('PHI100', 'SUP110', 'SUPM', 'PPHI100', 'INS', 'FLAT', 'PTAX2'):
    X = load(a); c = ccf(X); i0 = ix[0]
    vals[a] = dict(births=[100 * (X['b'][i] / B['b'][i] - 1) for i in ix], cf2023=100 * (c[0] / cb[0] - 1),
                   price=[100 * (X['p'][i] / B['p'][i] - 1) for i in ix], rent=[100 * (X['r'][i] / B['r'][i] - 1) for i in ix],
                   households=[100 * (X['N'][i] / B['N'][i] - 1) for i in ix], own1829_pp=[100 * (X['o1829'][i] - B['o1829'][i]) for i in ix],
                   supply_multiplier=float((X['H'][i0] / X['p'][i0] ** ETA) / (B['H'][i0] / B['p'][i0] ** ETA)))
# 2023 cohort, quantum and tempo: completed fertility, share with a first birth, mean age at first birth and at any birth
# (diagonal sums of per-household rates from the cohort aged 18-21 in 2023; mean ages are rate-weighted age-cell midpoints)
def cohort(arm):
    d = json.load(open(D / f'prod_{arm}.json')); F = d['last']['fertility']
    A = np.array([f['age_specific_birth_rate_topcode_adjusted'] for f in F]); FB = np.array([f['first_birth_rate_all_households'] for f in F])
    mid = np.array(F[0]['age_cell_midpoint'], float); a = A[np.arange(7), np.arange(7)]; fb = FB[np.arange(7), np.arange(7)]
    return dict(cf=float(a.sum()), ever_first=float(fb.sum()), mean_age_first=float((fb * mid[:7]).sum() / fb.sum()),
                mean_age_birth=float((a * mid[:7]).sum() / a.sum()))
vals['cohort2023'] = {a: cohort(a) for a in ('BASE', 'PHI100', 'SUPM', 'PPHI100')}
json.dump(vals, open(OUT / 'p6_values_oldsetup.json', 'w'), indent=1)
# quantum and tempo figure: births and households, credit vs matched supply, 2023 to the last date
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.8))
for ax, k, ttl in zip(axs, ('b', 'N'), ('(a) Births', '(b) Households')):
    for a, lab, col, mk in (('PHI100', 'Mortgage credit', RUST, 'o'), ('SUPM', 'Matched supply', NAVY, 's')):
        X = load(a); ax.plot(X['y'], 100 * (X[k] / B[k] - 1), color=col, marker=mk, ms=2, lw=1.4, label=lab)
    ax.axhline(0, color='0.6', lw=0.6); ax.set_title(ttl, fontsize=10, loc='left'); ax.set_ylabel('% from no policy'); ax.set_xlabel('Year')
    ax.set_xticks([2023, 2103, 2203, 2303, 2403]); ax.grid(axis='y', alpha=.25, lw=.5)
axs[0].legend(frameon=False, fontsize=8, loc='lower right')
fig.tight_layout(); fig.savefig(OUT / 'fig_p6_qt_oldsetup_eta063.pdf'); plt.close(fig)
# cash, insurance and property tax: births and house price, 2023-2103
GOLD, GREY = '#c08a1e', '#6e6e6e'; m = B['y'] <= 2103
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.8))
for ax, k, ttl in zip(axs, ('b', 'p'), ('(a) Births', '(b) House price')):
    for a, lab, col, mk in (('FLAT', 'Cash transfer', GOLD, 'o'), ('INS', 'Earnings insurance', TEAL, '^'), ('PTAX2', 'Property tax, rebated', GREY, 'd')):
        X = load(a); ax.plot(X['y'][m], 100 * (X[k][m] / B[k][m] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
    ax.axhline(0, color='0.6', lw=0.6); ax.set_title(ttl, fontsize=10, loc='left'); ax.set_ylabel('% from no policy'); ax.set_xlabel('Year')
    ax.set_xticks([2023, 2043, 2063, 2083, 2103]); ax.grid(axis='y', alpha=.25, lw=.5)
axs[0].legend(frameon=False, fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig(OUT / 'fig_p6_transfers_oldsetup_eta063.pdf'); plt.close(fig)
m = B['y'] <= 2103; yc = B['y'][:len(cb)]; mc = yc <= 2063
ARMS = (('PHI100', 'Mortgage credit', RUST, 'o'), ('SUPM', 'Matched supply', NAVY, 's'), ('PPHI100', 'Credit for parents', TEAL, '^'))
fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.0))
for (k, ttl), ax in zip((('b', '(a) Births'), ('cf', '(b) Completed fertility by cohort'), ('N', '(c) Households'), ('p', '(d) House price')), axs.ravel()):
    for a, lab, col, mk in ARMS:
        X = load(a)
        if k == 'cf':
            ax.plot(yc[mc], 100 * (ccf(X)[mc] / cb[mc] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
        else:
            ax.plot(X['y'][m], 100 * (X[k][m] / B[k][m] - 1), color=col, marker=mk, ms=3, lw=1.4, label=lab)
    ax.axhline(0, color='0.6', lw=0.6); ax.set_title(ttl, fontsize=10, loc='left'); ax.set_ylabel('% from no policy')
    ax.set_xticks([2023, 2043, 2063, 2083, 2103] if k != 'cf' else [2023, 2033, 2043, 2053, 2063]); ax.grid(axis='y', alpha=.25, lw=.5)
axs[1, 0].set_xlabel('Year'); axs[1, 1].set_xlabel('Year'); axs[0, 1].set_xlabel('Year the cohort is aged 18 to 21', fontsize=9)
axs[0, 0].legend(frameon=False, fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig(OUT / 'fig_p6_paths_oldsetup_eta063.pdf'); plt.close(fig)
print('base', {k: np.round(x, 3) for k, x in vals['base'].items()})
print('cohort2023', {a: {k: round(x, 3) for k, x in v.items()} for a, v in vals['cohort2023'].items()})
for a, v in vals.items():
    if a in ('base', 'cohort2023'): continue
    print(f'{a:8s}', 'births', np.round(v['births'], 1), 'cf', round(v['cf2023'], 1), 'price', np.round(v['price'], 1), 'rent', np.round(v['rent'], 1),
          'own1829', np.round(v['own1829_pp'], 1), 'mult', round(v['supply_multiplier'], 4))
