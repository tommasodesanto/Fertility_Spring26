#!/usr/bin/env python3
"""Two main-text figures for the first-birth event study (A2h sample).

Figure 1: rooms (common support, the calibration target) and ownership.
Figure 2: moved for more space and moved for neighborhood (all households).
Reads saved coefficients only; no estimation. Style follows
build_first_birth_paper_figures.py. Writes one PDF per figure into
latex/JMP_DS_mock/figures/.
"""
import csv
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FIG = ROOT / 'latex/JMP_DS_mock/figures'
ROOMS_CS = ROOT / 'output/model/first_birth_rooms_event_20261010/measurement_bridge/original/rooms_coefficients.csv'
V3 = HERE / 'output/sa_first_birth_outcomes_v3'
WINDOWS = ['Wleft', 'Wm2', 'Wm1', 'REF', 'Wp1', 'Wp2', 'Wp3', 'Wp4', 'Wp5', 'Wp6', 'Wright']
LABELS = ['$\\leq-8$', '$-7/-6$', '$-5/-4$', '$-3/-2$', '$-1/0$', '$+1/+2$', '$+3/+4$',
          '$+5/+6$', '$+7/+8$', '$+9/+10$', '$\\geq+11$']
COLOR = '#1866b4'


def coefs(path, scale=1.0):
    d = {r['coefficient']: (float(r['estimate']) * scale, float(r['standard_error']) * scale)
         for r in csv.DictReader(open(path))}
    d['REF'] = (0.0, 0.0)
    m = np.array([d[w][0] for w in WINDOWS])
    s = np.array([d[w][1] for w in WINDOWS])
    return m, s


def panel(ax, m, s, title, ylabel):
    """Old-figure style: equally spaced windows, light horizontal grid, dotted line at the reference."""
    ref = WINDOWS.index('REF')
    x = np.arange(len(WINDOWS))
    keep = x != ref
    ax.grid(axis='y', color='#dddddd', linewidth=0.6)
    ax.set_axisbelow(True)
    ax.axvline(ref, color='#8a8a8a', linewidth=0.8, linestyle=':')
    ax.plot(x, m, color=COLOR, linewidth=1.3)
    ax.errorbar(x[keep], m[keep], yerr=1.96 * s[keep], color=COLOR, marker='o', markersize=4,
                linewidth=0, elinewidth=0.8, capsize=2)
    ax.plot([ref], [0], marker='D', markersize=7, markerfacecolor='white', markeredgecolor=COLOR,
            markeredgewidth=1.4, linestyle='none', zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, rotation=35, fontsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_title(title, fontsize=9.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_xlabel('Years relative to first birth', fontsize=8.5)


def main():
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/first_birth_main_figures_mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    FIG.mkdir(exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    panel(axes[0], *coefs(ROOMS_CS), '(a) Rooms', 'Change in rooms')
    panel(axes[1], *coefs(V3 / 'A2h_own/coefficients.csv', 100.0), '(b) Ownership',
          'Change in percentage points')
    fig.tight_layout()
    fig.savefig(FIG / 'first_birth_rooms_ownership_2x1.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    panel(axes[0], *coefs(V3 / 'A2h_moved_space/coefficients.csv', 100.0), '(a) Moved for space',
          'Change in percentage points')
    panel(axes[1], *coefs(V3 / 'A2h_moved_nbhd/coefficients.csv', 100.0), '(b) Moved for neighborhood',
          'Change in percentage points')
    fig.tight_layout()
    fig.savefig(FIG / 'first_birth_move_reasons_2x1.pdf')
    plt.close(fig)
    print(FIG / 'first_birth_rooms_ownership_2x1.pdf')
    print(FIG / 'first_birth_move_reasons_2x1.pdf')


if __name__ == '__main__':
    main()
