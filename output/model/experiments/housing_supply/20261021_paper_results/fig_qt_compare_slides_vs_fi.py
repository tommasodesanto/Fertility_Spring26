"""Oct 23 2026 (his ask 15:02Z): the slides' quantum-tempo figure (credit vs matched supply, births and households, % vs no policy,
to 2200; decay_ship_20261014/figures/test_quantum_tempo/propB_births_households_two_panels.pdf) redrawn side by side for the
slides' runs (legacy: 20261017_accounting/transitions/data/prod_{BASE,PHI100,SUPM}.json, H0 x1.063) and free entry eta 0.63
(data_fi_torch/P063fi_{BASE,PHI100,SUPMC}.json, H0 x1.0435). Same style as the slide, shared y axes. Saved results only."""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; T = HERE.parent / '20261017_accounting/transitions'
def load(f):
    R = json.load(open(f))['last']['rows']
    return np.array([r['calendar_year'] for r in R]) + 2, {k: np.array([r[n] for r in R]) for k, n in (('b', 'birth_children'), ('N', 'adult_population'))}
SETS = (("Slides' runs (supply +6.3%)", T / 'data/prod_BASE.json', T / 'data/prod_PHI100.json', T / 'data/prod_SUPM.json'),
        ('Free entry (supply +4.35%)', T / 'data_fi_torch/P063fi_BASE.json', T / 'data_fi_torch/P063fi_PHI100.json', T / 'data_fi_torch/P063fi_SUPMC.json'))
RED, BLUE = '#c0392b', '#2c5f99'
fig, axs = plt.subplots(2, 2, figsize=(11, 7.4), sharey='row')
for c, (ttl, fb, fc, fs) in enumerate(SETS):
    y, B = load(fb); _, C = load(fc); _, S = load(fs); m = y <= 2200
    for r, (k, lab) in enumerate((('b', 'Births in each period'), ('N', 'Number of households'))):
        ax = axs[r, c]
        ax.plot(y[m], 100 * (S[k][m] / B[k][m] - 1), color=RED, lw=2.6, label='Housing supply')
        ax.plot(y[m], 100 * (C[k][m] / B[k][m] - 1), color=BLUE, lw=2.6, label='Credit, $\\phi=1$')
        if k == 'b':
            ax.plot(y[0], 100 * (S[k][0] / B[k][0] - 1), 'o', color=RED, ms=7); ax.plot(y[0], 100 * (C[k][0] / B[k][0] - 1), 'o', color=BLUE, ms=7)
        ax.axhline(0, color='0.75', lw=1); ax.set_xlim(2023, 2200); ax.set_title(lab if c == 0 else lab, loc='left', fontsize=12)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        if c == 0: ax.set_ylabel('% vs no policy')
        if r == 1: ax.set_xlabel('Year')
    axs[0, c].annotate(ttl, xy=(0.5, 1.16), xycoords='axes fraction', ha='center', fontsize=13, fontweight='bold')
axs[0, 0].legend(frameon=False, loc='upper left')
fig.tight_layout(); fig.savefig(HERE / 'fig_qt_compare_slides_vs_fi.pdf'); plt.close(fig)
