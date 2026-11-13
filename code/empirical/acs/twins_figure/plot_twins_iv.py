"""Twins IV figure for Appendix B.4: first stage, IV rooms, IV ownership by child age, two samples.
Reads output/national_acs_comparison/headship_iv_20261112/headship_iv_table.csv (estimates and SEs as rounded there).
Usage: python plot_twins_iv.py ROOT OUTDIR"""
import csv, re, sys
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT, OUT = Path(sys.argv[1]), Path(sys.argv[2])
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True,
                     'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8,
                     'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
rows = list(csv.DictReader(open(ROOT / 'output/national_acs_comparison/headship_iv_20261112/headship_iv_table.csv')))
def pe(s):
    m = re.match(r'\s*(-?[\d.]+)\s*\(([\d.]+)\)', s); return float(m.group(1)), float(m.group(2))
cases = [('Twins at first birth, child ages 0-5 pooled', 'Ages 0-5'), ('Twins, child age 3', 'Age 3'), ('Twins, child age 5', 'Age 5')]
samples = [('All mothers (saved run 18274624)', '#1866b4', -0.1, 'All mothers'), ('Householder/spouse/partner', '0.55', 0.1, 'Householder, spouse or partner')]
cols = [('first stage, pp more likely extra child (SE)', 'First stage (pp)'), ('IV rooms (SE)', 'Rooms'), ('IV ownership, pp (SE)', 'Owner-occupied (pp)')]
fig, axes = plt.subplots(1, 3, figsize=(9, 2.8))
for ax, (c, title) in zip(axes, cols):
    for ver, color, off, lab in samples:
        for i, (d, _) in enumerate(cases):
            r = next(x for x in rows if x['design'] == d and x['version'] == ver)
            b, se = pe(r[c])
            ax.errorbar(i + off, b, yerr=1.96 * se, fmt='o', ms=4, color=color, capsize=0, lw=1.2, label=lab if i == 0 else None)
    ax.axhline(0, color='k', lw=0.5); ax.set_xticks(range(3)); ax.set_xticklabels([x[1] for x in cases])
    ax.set_xlabel('Age of the twins'); ax.set_title(title); ax.set_xlim(-0.5, 2.5)
axes[0].set_ylim(0, 75); axes[0].legend(frameon=False, loc='lower left', fontsize=8)
fig.tight_layout(); fig.savefig(OUT / 'twins_iv_by_child_age.pdf')
