#!/usr/bin/env python3
"""Render the reviewed target catalogue; never alter a calibration contract."""
import argparse
import csv
import json
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, PageBreak, Image, Spacer


def build(packet, output):
    packet = Path(packet)
    review = json.loads((packet / 'review.json').read_text())
    rows = review['rows']
    assert len(rows) == 13 and len({r['id'] for r in rows}) == 13
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('SmallReview', fontName='Helvetica', fontSize=10, leading=13, spaceAfter=6))
    styles.add(ParagraphStyle('ReviewTitle', fontName='Helvetica-Bold', fontSize=19, leading=23, textColor=colors.HexColor('#17324d'), spaceAfter=12))
    styles.add(ParagraphStyle('ReviewHeading', fontName='Helvetica-Bold', fontSize=11, leading=14, textColor=colors.HexColor('#17324d'), spaceBefore=9, spaceAfter=5))
    def p(text, style='SmallReview'):
        return Paragraph(escape(str(text)), styles[style])
    def table(data, widths):
        t = Table([[p(x) for x in row] for row in data], colWidths=widths, repeatRows=1, hAlign='LEFT')
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#dce8f0')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f3f6f8')]),
            ('GRID', (0, 0), (-1, -1), .3, colors.HexColor('#bcc9d3')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 5), ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        return t
    story = [p('Calibration targets: what is verified?', 'ReviewTitle'), p(review['date'])]
    for line in review['summary']:
        story.append(p(line))
    story += [p('Complete current target list', 'ReviewHeading'), table(
        [['Moment', 'Current target', 'Remaining decision']] +
        [[r['label'], r['display_target'], r['short_issue']] for r in rows], [205, 83, 236])]
    for group, heading in [('housing', 'Housing: source and model measurement'), ('fertility', 'Fertility: source and model measurement'), ('wealth', 'Wealth: source and model measurement')]:
        story += [PageBreak(), p(heading, 'ReviewTitle')]
        for r in [x for x in rows if x['group'] == group]:
            story += [p(r['label'] + ' | ' + r['display_target'], 'ReviewHeading'),
                      p('Data: ' + r['data']), p('Model: ' + r['model']),
                      p('Assessment: ' + r['assessment'])]
    story += [PageBreak(), p('Weights and decisions before launch', 'ReviewTitle'),
              p('Current objective: sum of twelve squared gaps divided by fixed working scales squared. The separate 2.1 normalization is unscored. These scales are not all sampling standard errors; cross-moment covariance is not used.'),
              table([['Moment', 'Working scale', 'Basis']] + [[r['label'], '-' if r['working_scale'] is None else f"{r['working_scale']:.6g}", r['weight_basis']] for r in rows], [207, 87, 230])]
    story.append(p('Recommended closure order', 'ReviewHeading'))
    for i, line in enumerate(review['next_steps'], 1):
        story.append(p(f'{i}. {line}'))
    story += [PageBreak(), p('Evidence and limits of this review', 'ReviewTitle')]
    for label, source in review['sources']:
        story += [p(label, 'ReviewHeading'), p(source)]
    story.append(p('This review reads saved empirical receipts and the actual scored observer/source. It does not rerun microdata, change estimates or weights, certify structural identification, or submit a model job. Detailed per-row receipts and the complete target catalogue accompany this PDF.'))
    def footer(canvas, doc):
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#52616b'))
        canvas.drawString(44, 25, 'Target review | no target changes or calibration submission')
        canvas.drawRightString(568, 25, str(doc.page))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), pagesize=letter, leftMargin=44, rightMargin=44, topMargin=40, bottomMargin=42,
                      title='Calibration targets: verification and decisions').build(story, onFirstPage=footer, onLaterPages=footer)


def branch_diagnostic_plot(packet, diagnostic):
    """Supplemental two-date graph; never replace the standard 17 plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.2, 3.8))
    for index, item in enumerate(diagnostic['branch_results']):
        result = json.loads((packet / item['receipt']).read_text())
        if result.get('d1_gate') != 'PASS' or result.get('status') != 'complete':
            raise ValueError('Branch plot requires verified complete measurement')
        ys = [result['d0_origin_postbirth_minus_control'], result['d1_destination_birth_minus_control']]
        color = ['#245d80', '#bc6c24'][index]
        axis.plot([0, 4], ys, 'o--', color=color, label=item['label'], linewidth=1.4)
        for x, y in zip([0, 4], ys):
            axis.annotate(f'{y:.3f}', (x, y), xytext=(0, -17 if index == 0 else 9),
                          textcoords='offset points', ha='center', color=color, fontsize=10)
    axis.set_xticks([0, 4], ['Origin decision (D0)', 'Next decision: 4 years later (D1)'])
    axis.set_xlim(-.6, 4.6)
    axis.set_ylim(0, 1.3)
    axis.set_ylabel('Birth minus childless control (rooms)')
    axis.set_title('Housing response at two model dates', fontsize=13)
    axis.grid(axis='y', color='#dddddd', linewidth=.6)
    axis.spines[['top', 'right']].set_visible(False)
    axis.legend(loc='upper left', frameon=False)
    figure.tight_layout()
    path = packet / 'lunch_first_birth_branch_comparison.png'
    figure.savefig(path, dpi=180, bbox_inches='tight')
    figure.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(figure)
    return path


def build_decision_review(packet, output):
    """Render a lead-authored decision JSON without altering the legacy report mode."""
    packet = Path(packet)
    review = json.loads((packet / 'lunch_decision_review.json').read_text())
    contract_path = packet / 'overnight/report_tables/target_contract.csv'
    national_path = packet.parent / 'housing_profiles_v1/full/target_recomputed.json'
    with contract_path.open(newline='') as handle:
        contract = {row['internal_key']: float(row['current_value']) for row in csv.DictReader(handle)}
    national = json.loads(national_path.read_text())['recomputed']['national']
    national_override_keys = {
        'mean_rooms', 'ownership_30_55', 'family_rooms', 'recent_parent_ownership'
    }
    if set(national) != national_override_keys:
        raise ValueError('National source must contain exactly the four authorized housing overrides')
    checks = []
    for row in review['rows']:
        key, value = row['key'], float(row['target'])
        if key not in contract:
            raise ValueError(f'Missing target in contract: {key}')
        expected = float(national[key]) if key in national_override_keys else contract[key]
        gap = value - expected
        checks.append({'key': key, 'review_value': value, 'expected_value': expected,
                       'source': 'national_override' if key in national_override_keys else 'frozen_contract',
                       'absolute_gap': abs(gap), 'matches': abs(gap) <= 1e-12})
    if len(checks) != 13 or {r['key'] for r in review['rows']} != set(contract) or not all(x['matches'] for x in checks):
        raise ValueError('Lunch review rows do not match the authorized national/frozen target sources')
    (packet / 'lunch_decision_table_checks.json').write_text(json.dumps({
        'schema': 'lunch_decision_table_checks_v2_national_housing', 'contract_path': str(contract_path),
        'national_override_path': str(national_path), 'authorized_national_override_keys': sorted(national_override_keys),
        'rows_checked': len(checks), 'all_match_within_1e-12': True,
        'frozen_rows_unchanged': all(x['matches'] for x in checks if x['source'] == 'frozen_contract'),
        'national_overrides_match_source': all(x['matches'] for x in checks if x['source'] == 'national_override'),
        'checks': checks
    }, indent=2) + '\n')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('LunchTitle', fontName='Helvetica-Bold', fontSize=18, leading=21,
                              textColor=colors.HexColor('#17324d'), spaceAfter=7))
    styles.add(ParagraphStyle('LunchHead', fontName='Helvetica-Bold', fontSize=10.2, leading=12,
                              textColor=colors.HexColor('#17324d'), spaceBefore=5, spaceAfter=2))
    styles.add(ParagraphStyle('Lunch', fontName='Helvetica', fontSize=8.2, leading=10.1, spaceAfter=3))
    styles.add(ParagraphStyle('LunchTiny', fontName='Helvetica', fontSize=7.25, leading=8.55, spaceAfter=2))
    def p(text, style='Lunch'):
        return Paragraph(escape(str(text)), styles[style])
    def make_table(data, widths, small='LunchTiny'):
        result = Table([[p(cell, small) for cell in row] for row in data], colWidths=widths,
                       repeatRows=1, hAlign='LEFT')
        result.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#dce8f0')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f3f6f8')]),
            ('GRID', (0, 0), (-1, -1), .25, colors.HexColor('#bcc9d3')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
        ]))
        return result
    story = [p(review['title'], 'LunchTitle'), p(review['date'] + ' | ' + review['status'])]
    for line in review['summary']:
        story.append(p(line))
    story.append(p('Six author decisions', 'LunchHead'))
    for item in review['decisions']:
        story += [p(item['decision'], 'LunchHead'), p('Recommendation: ' + item['recommendation']),
                  p('Next action: ' + item['next_operation'], 'LunchTiny')]
    story += [PageBreak(), p('Current target table', 'LunchTitle'),
              p('Parameter connections describe economic margins, not one-to-one identification. Confidence concerns measurement comparability, not model fit.')]
    target_data = [['Target', 'Value', 'Parameter connection', 'Confidence / decision']]
    for row in review['rows']:
        target_data.append([row['label'], row['display'], row['parameter_connection'],
                            row['confidence'] + ': ' + row['recommendation']])
    story.append(make_table(target_data, [132, 68, 139, 173]))
    story += [PageBreak(), p('Decision evidence and implementation boundary', 'LunchTitle')]
    for item in review['decisions']:
        story += [p(item['decision'], 'LunchHead'), p('Evidence: ' + item['evidence']),
                  p('Alternative: ' + item['alternative'], 'LunchTiny'), p('Consequence: ' + item['consequence'], 'LunchTiny')]
    story += [PageBreak(), p('Other documented choices', 'LunchTitle')]
    for line in review['other_decisions']:
        story.append(p(line))
    story.append(p('Scope and evidence', 'LunchHead'))
    for item in review['source_evidence']:
        story.append(p(item, 'LunchTiny'))
    story.append(p('Not done or claimed', 'LunchHead'))
    for item in review['not_done_or_claimed']:
        story.append(p(item))
    diagnostic = review.get('first_birth_diagnostic', {})
    comparisons = diagnostic.get('comparison_table', [])
    plot = diagnostic.get('plot_path')
    if diagnostic.get('branch_results'):
        plot = branch_diagnostic_plot(packet, diagnostic)
    if comparisons or plot:
        story += [PageBreak(), p('Supplemental first-birth diagnostic', 'LunchTitle'),
                  p(diagnostic.get('estimator', ''))]
        for paragraph in diagnostic.get('paragraphs', []):
            story.append(p(paragraph))
        if comparisons:
            story.append(make_table(comparisons, diagnostic.get('table_widths', [150, 115, 247])))
        if plot and Path(plot).exists():
            story.append(Spacer(1, 10))
            story.append(Image(str(plot), width=470, height=265, kind='proportional'))
            story.append(p(diagnostic.get('plot_caption', ''), 'LunchTiny'))
    def footer(canvas, doc):
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#52616b'))
        canvas.drawString(44, 25, 'Calibration decisions | September 23, 2026')
        canvas.drawRightString(568, 25, str(doc.page))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), pagesize=letter, leftMargin=44, rightMargin=44, topMargin=46, bottomMargin=42,
                      title=review['title']).build(story, onFirstPage=footer, onLaterPages=footer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--packet', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--decision-review', action='store_true', help='render lunch_decision_review.json')
    args = parser.parse_args()
    if args.decision_review:
        build_decision_review(args.packet, args.output)
    else:
        build(args.packet, args.output)
