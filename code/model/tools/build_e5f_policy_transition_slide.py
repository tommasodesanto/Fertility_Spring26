"""Plot frozen matched tax paths in model units, without a new model solve."""
from pathlib import Path
import csv
import hashlib
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / 'output/model/e5f_original_queue_20260913a'
OUT = PACKET / 'inherited_2023_tax/transition_readout'


def main():
    source = OUT / 'frozen_policy_paths.json'
    frozen = json.loads(source.read_text())
    paths, diagnostics = {}, {}
    for arm in ['tax1', 'tax2']:
        item = frozen[arm]
        rows = item['files']['rows.json']
        fertility = item['files']['fertility.json']
        assert len(rows) == len(fertility) == 100
        coordinates = [r[k] for k in ['asset_price', 'pension_period_units',
                                     'equal_transfer_period_units'] for r in rows]
        assert coordinates == item['best_receipt']['prices']
        assert item['best_receipt']['mapping_valid']
        gates = item['files']['native_gates.json']
        assert gates['mass'] < 1e-12 and gates['policy'] == gates['feasibility'] == 0
        paths[arm] = []
        for r, f in zip(rows, fertility):
            assert r['calendar_year'] == f['calendar_year']
            assert r['psi_child'] == .09221854783921073
            assert r['payroll_tax_rate'] == .179
            assert r['net_migrant_heads_over_period'] == 0
            paths[arm].append(dict(year=r['calendar_year'],
                households=r['adult_population'],
                tfr=f['period_tfr_topcode_adjusted'],
                births=r['birth_children_topcode_adjusted']))
        diagnostics[arm] = dict(native_gates=gates,
            max_housing_relative_gap=max(abs(r['relative_market_residual']) for r in rows),
            max_pension_relative_gap=max(abs(r['scaled_pension_budget_residual']) for r in rows),
            max_rebate_relative_gap=max(abs(r['scaled_government_budget_residual']) for r in rows),
            terminal_distance=item['files']['terminal_distance.json'])
    assert math.isclose(paths['tax1'][0]['households'], paths['tax2'][0]['households'], abs_tol=1e-12)
    rows = []
    for b, p in zip(paths['tax1'], paths['tax2']):
        assert b['year'] == p['year']
        rows.append(dict(year=b['year'], baseline_tfr=b['tfr'], policy_tfr=p['tfr'],
            baseline_households=b['households'], policy_households=p['households'],
            baseline_births=b['births'], policy_births=p['births'],
            household_difference_percent=100*(p['households']/b['households']-1),
            birth_difference_percent=100*(p['births']/b['births']-1)))
    with (OUT/'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    shown = [r for r in rows if r['year'] <= 2063]
    years = [r['year'] for r in shown]
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':12,
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'pdf.fonttype':42})
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.25))
    plotted = []
    for ax, key, title, ylabel in zip(axes, ['tfr','households'],
            ['Fertility', 'Household population'],
            ['Total fertility rate', 'Model units']):
        for prefix, label, color, style in [('baseline','1% tax','#246395','-'),
                                           ('policy','2% tax','#ce6c26','--')]:
            values = [r[f'{prefix}_{key}'] for r in shown]
            line, = ax.plot(years, values, label=label, color=color, ls=style, lw=2.5)
            np.testing.assert_array_equal(line.get_ydata(), values)
            plotted.append({'series':f'{prefix}_{key}', 'values':values})
        ax.set(title=title, ylabel=ylabel, xlim=(2023,2063))
        ax.set_xticks([2023,2031,2039,2047,2055,2063])
        ax.tick_params(axis='x', labelsize=10)
        ax.grid(axis='y', alpha=.18)
    axes[0].set_ylim(1.68, 1.80)
    axes[1].set_ylim(.85, 1.01)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.52,1.025), ncol=2, frameon=False)
    fig.subplots_adjust(left=.07, right=.965, bottom=.15, top=.84, wspace=.31)
    for suffix in ['pdf','png']:
        fig.savefig(OUT/f'policy_transition.{suffix}', dpi=180)
    plt.close(fig)
    tex = r'''\documentclass[11pt,aspectratio=169]{beamer}
\usepackage[T1]{fontenc}
\usepackage{lmodern,graphicx}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\begin{document}
\begin{frame}{Policy Results}
\small Annual property tax from 1\% to 2\% in 2023; equal rebates in both paths.
\par\smallskip
\includegraphics[width=\textwidth]{../output/model/e5f_original_queue_20260913a/inherited_2023_tax/transition_readout/policy_transition.pdf}
\par\smallskip
\footnotesize Preliminary transitions: market convergence is incomplete.\\
Common inherited 2023 economy after one permanent preference shock; no immigration.\\
Households are in model units (2023 $\simeq 1$), not resident-person counts.
\end{frame}
\end{document}
'''
    (ROOT/'latex/appendix_property_tax_transition.tex').write_text(tex)
    verification = dict(artifact_checks='PASS', equilibrium_status='NOT CONVERGED',
        model_solves=0, population_transformation_in_graphs=False,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        plotted_arrays=plotted, diagnostics=diagnostics,
        effects_2063=shown[-1])
    (OUT/'verification.json').write_text(json.dumps(verification, indent=2)+'\n')
    print(json.dumps(shown[-1], indent=2))


if __name__ == '__main__':
    main()
