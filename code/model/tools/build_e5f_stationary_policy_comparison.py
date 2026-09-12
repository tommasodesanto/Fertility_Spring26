"""Collect the three verified, conditional stationary property-tax equilibria."""
from pathlib import Path
import argparse
import csv
import json

import numpy as np

CASES = ('baseline', 'equal-rebate-1pct', 'equal-rebate-2pct')
LABELS = ('1% tax, no rebate', '1% tax, equal rebate', '2% tax, equal rebate')


def build(root, pdf=None):
    root = Path(root)
    contract = json.loads((root / 'contract.json').read_text())
    summaries = [json.loads((root / 'results' / case / 'summary.json').read_text()) for case in CASES]
    for case, summary in zip(CASES, summaries):
        assert summary['case'] == case
        assert summary['status'] == 'stationary_equilibrium_passed'
        assert summary['psi'] == contract['psi']
        assert np.max(np.abs(summary['root_residual'])) <= 2e-4
        assert not summary['production_eligible']
        receipts = [json.loads(p.read_text()) for p in sorted((root / 'results' / case).glob('terminal_*/root_receipt.json'))]
        receipt = receipts[-1]
        assert receipt['converged']
        assert np.array_equal(receipt['final']['prices'], summary['root_coordinates'])
        payload = receipt['final']['payload']
        assert all(payload['household_gates'].values())
        assert len(list((root / 'results' / case / 'graphs/standard_diagnostics').glob('*.png'))) == 17
        summary['property_tax_revenue'] = payload['rebate_revenue']
        summary['rebate_outlays'] = payload['rebate_outlays']
        if summary['equal_rebate']:
            assert abs(summary['property_tax_revenue'] - summary['rebate_outlays']) / max(abs(summary['property_tax_revenue']), abs(summary['rebate_outlays']), 1e-12) <= 1e-6
        summary['rooms_per_household'] = summary['housing_demand'] / summary['household_heads']
        summary['births_per_household'] = summary['births'] / summary['household_heads']
        summary['ownership_percent'] = 100 * summary['ownership']
        summary['with_children_percent'] = 100 * summary['with_children']
    rows = []
    specifications = (
        ('period_fertility', 'Period fertility', 4),
        ('ownership_percent', 'Homeownership (%)', 2),
        ('rooms_per_household', 'Rooms per household', 4),
        ('asset_price', 'House price per room', 4),
        ('rent_per_room', 'Period rent per room', 4),
        ('rebate_per_head_period', 'Period rebate per household', 4),
        ('pension_period', 'Period pension benefit', 4),
        ('property_tax_revenue', 'Property-tax revenue (model scale)', 6),
        ('rebate_outlays', 'Equal-rebate outlays (model scale)', 6),
        ('with_children_percent', 'Households with dependents (%)', 2),
        ('births_per_household', 'Births per household per period', 4),
        ('births', 'Total period births (model scale)', 6),
        ('household_heads', 'Households (model scale)', 6),
        ('resident_persons', 'Persons (model scale)', 6),
        ('housing_demand', 'Total rooms (model scale)', 6),
    )
    lines = ['# Conditional stationary property-tax comparison', '',
             'All three cases use the same structural parameters and fitted preference. '
             'Payroll taxes finance pensions in each equilibrium. Equal rebates return '
             'property-tax revenue per household. These are long-run stationary comparisons; '
             'they are not policy effects from the inherited 2023 economy.', '',
             '| Outcome | ' + ' | '.join(LABELS) + ' |', '|---|---:|---:|---:|']
    for key, label, digits in specifications:
        values = [float(s[key]) for s in summaries]
        rows.append(dict(outcome=key, label=label, **dict(zip(CASES, values))))
        lines.append('| ' + label + ' | ' + ' | '.join(f'{v:.{digits}f}' for v in values) + ' |')
    lines += ['', 'The demographic and entry closure is retained from the diagnostic '
              'patch endpoint. The historical patch still fails terminal-distance certification. '
              'Population levels use that retained normalization.', '',
              f'Fertility preference: {contract["psi"]:.15g}. '
              'The initial calibration and its target fingerprint are pinned in contract.json.', '']
    (root / 'comparison.md').write_text('\n'.join(lines))
    with (root / 'comparison.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=('outcome', 'label', *CASES))
        writer.writeheader()
        writer.writerows(rows)
    changes = []
    for summary in summaries:
        changes.append(dict(case=summary['case'], **{
            key + '_change_percent': 100 * (summary[key] / summaries[0][key] - 1)
            for key in ('period_fertility', 'births', 'asset_price', 'rent_per_room',
                        'rooms_per_household', 'household_heads', 'resident_persons')},
            ownership_change_pp=100 * (summary['ownership'] - summaries[0]['ownership'])))
    (root / 'comparison.json').write_text(json.dumps(dict(summaries=summaries, changes=changes), indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), constrained_layout=True)
    short = ('1% tax\nno rebate', '1% tax\nequal rebate', '2% tax\nequal rebate')
    for ax, key, title in zip(axes.flat,
                             ('period_fertility', 'ownership_percent', 'rooms_per_household', 'asset_price'),
                             ('Period fertility', 'Homeownership (%)', 'Rooms per household', 'House price per room')):
        values = [s[key] for s in summaries]
        ax.bar(short, values, color=('#555555', '#437DA0', '#C78648'), width=.6)
        ax.set_title(title)
        ax.set_ylim(0, max(values) * 1.18)
        ax.spines[['top', 'right']].set_visible(False)
        for i, value in enumerate(values):
            ax.text(i, value + max(values) * .025, f'{value:.3f}', ha='center', fontsize=10)
    fig.suptitle('Stationary policy comparison — conditional on the fitted patch preference', fontsize=13)
    fig.savefig(root / 'comparison.png', dpi=160)
    plt.close(fig)
    if pdf is not None:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
        pdf = Path(pdf)
        pdf.parent.mkdir(parents=True, exist_ok=True)
        styles = getSampleStyleSheet()
        styles['BodyText'].fontSize = 9
        styles['BodyText'].leading = 12
        story = [Paragraph('Stationary property-tax comparison', styles['Title']),
                 Paragraph('September 12, 2026 | Same structural calibration and fertility preference in every case', styles['BodyText']),
                 Spacer(1, 10)]
        table = [['Outcome', '1% tax\nNo rebate', '1% tax\nEqual rebate', '2% tax\nEqual rebate']]
        for row, (_, _, digits) in zip(rows, specifications):
            table.append([row['label'], *[f'{row[case]:.{digits}f}' for case in CASES]])
        t = Table(table, colWidths=[247, 85, 85, 85], repeatRows=1)
        t.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (1, 0), (-1, -1), 'RIGHT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LINEBELOW', (0, 0), (-1, 0), .7, colors.black),
            ('LINEBELOW', (0, -1), (-1, -1), .7, colors.black),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f3f5f7')]),
        ]))
        incremental = 100 * (summaries[2]['period_fertility'] / summaries[1]['period_fertility'] - 1)
        story += [t, Spacer(1, 12),
                  Paragraph(f'<b>Fertility:</b> rebating the existing tax raises period fertility by {changes[1]["period_fertility_change_percent"]:.2f}%. The higher tax with rebate raises it by {changes[2]["period_fertility_change_percent"]:.2f}% relative to the unrebated baseline, or {incremental:.2f}% relative to the already rebated 1% case.', styles['BodyText']),
                  Paragraph('<b>Housing:</b> total rooms increase, but households increase faster, so rooms per household fall. Ownership changes are small. This comparison does not establish an improvement in housing allocation.', styles['BodyText']),
                  Paragraph('<b>Checks:</b> all three stationary equilibria pass housing clearance, pension balance, household feasibility and exact reproduction. Both rebate budgets balance separately. The standard 17 diagnostic figures are retained for every case.', styles['BodyText']),
                  Paragraph('<b>Scope:</b> long-run stationary comparisons under the retained diagnostic demographic and entry closure. These are not impacts from the inherited 2023 economy. The historical patch still fails terminal-distance certification. Population and total-birth levels use the retained model normalization; period fertility is distinct from total births.', styles['BodyText']),
                  Paragraph(f'<b>Specification:</b> annual tax rates 1%, 1%, 2%; equal per-household rebate where shown; payroll tax 17.9%; housing supply elasticity 0.63; fertility preference {contract["psi"]:.10f}. Monetary flows cover the four-year model period. Source jobs: 17499515, 17499516, 17499517. Full receipts and inputs: stationary_policy_comparison_fit/contract.json and results/.', styles['BodyText'])]
        SimpleDocTemplate(str(pdf), pagesize=A4, leftMargin=42, rightMargin=42, topMargin=35, bottomMargin=35).build(story)
    print(json.dumps(changes, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--pdf', type=Path)
    args = parser.parse_args()
    build(args.root, args.pdf)
