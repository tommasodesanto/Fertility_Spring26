"""Figure: ACS 2023, all US, owner-occupied homes with 6+ rooms by age of householder and child under 18.
Same style as code/data/ahs_supply_snapshot/plot_appendix_A_figures_v1.py (Figure A.2). Usage: python plot_large_homes_national.py ROOT OUTDIR"""
import csv, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT, OUT = Path(sys.argv[1]), Path(sys.argv[2])
BLUE = '#1866b4'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': 0.8, 'pdf.fonttype': 42,
                     'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
rows = list(csv.DictReader(open(ROOT / 'output/empirical_acs_large_homes_national_20261112/large_owner_national_2023.csv')))
order = [('22-39', 'False'), ('22-39', 'True'), ('40-59', 'False'), ('40-59', 'True'), ('60-85', 'False'), ('60-85', 'True')]
sh = {(r['age_group'], r['children_under_18']): 100 * float(r['share']) for r in rows}
v = [sh[k] for k in order]
xl = ['22-39\nno child', '22-39\nchild', '40-59\nno child', '40-59\nchild', '60-85\nno child', '60-85\nchild']
fig, ax = plt.subplots(figsize=(5.6, 2.7))
ax.bar(np.arange(6), v, 0.6, color=[BLUE if k[1] == 'True' else '#7fa6d2' for k in order])
for i, y in enumerate(v):
    ax.text(i, y + 0.8, f'{y:.1f}', ha='center', va='bottom', fontsize=8.5)
ax.set_xticks(np.arange(6)); ax.set_xticklabels(xl, fontsize=8.5)
ax.set_xlabel('Age of householder and child under 18 at home')
ax.set_ylabel('Share of large owned homes (%)'); ax.set_ylim(0, 48)
ax.grid(axis='y', color='#dddddd', linewidth=0.6); ax.set_axisbelow(True)
fig.savefig(OUT / 'large_homes_by_age_children_acs2023_national.pdf')
print('sum', sum(v), 'no child', v[0] + v[2] + v[4])
