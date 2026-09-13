#!/usr/bin/env python3
"""Build the verified initial readout and unchanged native figure appendix; no solve."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--as-of', required=True)
    args = ap.parse_args()
    initial = args.packet / 'corrected_initial'
    candidate = json.loads((initial / 'candidate_result.json').read_text())
    if candidate['status'] != 'verified' or candidate['second_signature_equal'] is not True:
        raise ValueError('Verified double replay required')
    with (initial / 'target_fit.csv').open() as stream:
        targets = list(csv.DictReader(stream))
    with (initial / 'parameters.csv').open() as stream:
        parameters = list(csv.DictReader(stream))
    if len(targets) != 13 or len(parameters) != 17:
        raise ValueError('Complete retained target and parameter tables required')
    loss = sum(float(r['loss_contribution']) for r in targets if r['loss_contribution'])
    if not math.isclose(loss, candidate['loss'], abs_tol=1e-7):
        raise ValueError('Full table does not reproduce the verified loss')

    navy = colors.HexColor('#19374b')
    teal = colors.HexColor('#16786f')
    gray = colors.HexColor('#eaf0f3')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('TitleCustom', fontName='Helvetica-Bold', fontSize=24,
                              leading=28, textColor=navy, spaceAfter=14))
    styles.add(ParagraphStyle('SectionCustom', fontName='Helvetica-Bold', fontSize=15,
                              leading=18, textColor=navy, spaceAfter=10))
    styles.add(ParagraphStyle('BodyCustom', fontSize=10, leading=14, spaceAfter=10))
    styles.add(ParagraphStyle('SmallCustom', fontSize=8, leading=11, spaceAfter=7))
    styles.add(ParagraphStyle('CellCustom', fontSize=8, leading=10))
    styles.add(ParagraphStyle('NumberCustom', fontSize=8, leading=10, alignment=TA_RIGHT))
    width = A4[0] - 80
    story = []

    def p(text, style='BodyCustom'):
        return Paragraph(text, styles[style])

    def add(text, style='BodyCustom'):
        story.append(p(text, style))

    def table(rows, widths, numeric_start=None):
        cells = [[p(escape(str(value)), 'NumberCustom' if numeric_start is not None
                    and i > 0 and j >= numeric_start else 'CellCustom')
                  for j, value in enumerate(row)] for i, row in enumerate(rows)]
        tab = Table(cells, colWidths=widths, repeatRows=1, hAlign='LEFT')
        tab.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), gray),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8fafb')]),
            ('LINEBELOW', (0, 0), (-1, 0), .7, teal),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('RIGHTPADDING', (0, 0), (-1, -1), 6),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ]))
        story.append(tab)
        story.append(Spacer(1, 12))

    def number(value):
        return '-' if value in ('', None) else f'{float(value):.6g}'

    add('Quantitative model<br/>Verified initial readout', 'TitleCustom')
    add(escape(args.as_of) + ' | Historical and policy work remains in progress.', 'SmallCustom')
    add('<b>The fiscal and distribution repairs pass. The economic fit remains weak.</b> '
        'The corrected initial equilibrium reproduces exactly in two independent numerical '
        'repetitions. Its objective is %.6f. This is a verified candidate from an incomplete '
        'search, not a converged calibration optimizer.' % candidate['loss'])
    table([
        ['Question', 'Verified evidence'],
        ['Does the pension budget balance?', 'Yes in the initial equilibrium and accepted first historical forecasts. Every later accepted date must pass the same PAYGO check.'],
        ['Is property tax rebated?', 'Yes. Revenue returns equally per current household head; the rebate budget is checked separately from pensions.'],
        ['Was the probability error fixed?', 'Yes. Float64 normalization in both Markov distribution operators passes the captured mass-error replay and full initial reconstruction.'],
        ['Is the entire shock path fitted?', 'Not yet. Both short-horizon variants have accepted the first historical window; subsequent surprises are being fitted.'],
        ['Are policy results ready?', 'Not yet. Baseline and higher-tax rebated policies follow each completed history from its inherited 2023 households.'],
    ], [158, width - 158])
    add('What deserves attention', 'SectionCustom')
    add('Mean rooms and ownership at ages 30-55 are the largest remaining fit problems. '
        'Annual beta reaches its 0.99 cap. The numerical repair changes the objective by '
        'less than one millionth; it is not evidence of improved economic fit.')
    add('The 6-, 24- and 100-period forecasts are separate horizon checks. Passing a finite '
        'market-clearing root does not establish that a shorter horizon is adequate. The '
        'alternative demographic closure remains unresolved and has not been promoted.')
    add('How the historical exercise works', 'SectionCustom')
    add('At each surprise, households expect the new fertility preference to remain constant. '
        'For guessed paths of house prices, pensions and rebates, the solver works backward '
        'over household choices and forward over inherited households. It adjusts those paths '
        'until markets and both fiscal budgets clear, then adjusts the preference to match '
        'that observation window. Only the first realized period is carried to the next surprise.')

    story.append(PageBreak())
    add('Complete initial target fit', 'TitleCustom')
    add('All 13 retained rows: 12 scored moments and the separate completed-fertility '
        'normalization. Objective %.6f. The target system and weights are unchanged.' % loss)
    rows = [['Moment', 'Target', 'Model', 'Gap', 'Weight', 'Loss']]
    rows += [[r['label'].replace('\u2212', '-'), *[number(r[k]) for k in
              ('target', 'model', 'gap', 'actual_weight', 'loss_contribution')]] for r in targets]
    table(rows, [211, 61, 61, 61, 65, width - 459], numeric_start=1)
    add('The normalization row does not enter the loss. Units follow each empirical moment: '
        'ownership is a fraction, room responses are in rooms, and first-birth age is in years. '
        'The CSV preserves every row\'s builder, sample, vintage, observation definition, '
        'uncertainty and weighting provenance.', 'SmallCustom')
    add('Interpretation', 'SectionCustom')
    add('Initial fertility levels are close by construction and the fertility-composition '
        'moments fit relatively well. Mean rooms are too high and ownership is too low. '
        'The first-birth room response overshoots its target. A numerically valid equilibrium '
        'therefore does not imply an adequate fit across the housing mechanisms.')
    add('Historical observation clock', 'SectionCustom')
    table([['Decision vintage', 'Observed TFR window', 'Target'],
           ['2007', '2008-2011', '1.974875'], ['2011', '2012-2015', '1.861000'],
           ['2015', '2016-2019', '1.755375'], ['2019', '2020-2023', '1.645750']],
          [150, 220, width - 370])
    add('A 2023 decision-vintage fertility flow belongs to the 2024-2027 forecast. It is '
        'not an additional fitted historical observation. Other 2023 model/data comparisons '
        'must use the correctly dated inherited economy.', 'SmallCustom')

    story.append(PageBreak())
    add('Parameters and restrictions', 'TitleCustom')
    add('Nine structural search coordinates. The other rows are fixed restrictions, '
        'normalizations or endogenous fiscal outcomes; they are not additional structural estimates.')
    rows = [['Parameter', 'Value', 'Lower', 'Upper', 'Near bound', 'Search']]
    rows += [[r['parameter'], number(r['estimate']), number(r['lower']), number(r['upper']),
              'Yes' if r['near_bound'] == 'True' else 'No',
              'Yes' if r['structural_coordinate'] == 'True' else 'No'] for r in parameters]
    table(rows, [171, 76, 64, 64, 77, width - 452], numeric_start=1)
    add('The annual-beta upper bound is the enforced 0.99. The scorer\'s older metadata '
        'is preserved in parameters_raw.csv. Near-bound flags use the retained reporting '
        'rule and should not be read as proof of binding first-order conditions.', 'SmallCustom')
    add('Numerical acceptance evidence', 'SectionCustom')
    table([['Initial check', 'Observed value'],
           ['Stationary nesting, L1', '2.2747454e-14'],
           ['One-period nesting, L1', '2.6512812e-13'],
           ['PAYGO relative budget gap', '3.0689413e-12'],
           ['Equal-rebate relative budget gap', '3.8384145e-7'],
           ['Second numerical signature', 'Exactly equal']], [285, width - 285])
    add('Numerical source: corrected_initial_source_v2. Job: 17655042. Source and '
        'target fingerprints are saved with the packet. Different serialized checkpoint '
        'hashes do not contradict identical numerical arrays.', 'SmallCustom')

    graphs = ['ownership_by_age', 'fertility_by_age', 'housing_prices', 'housing_market',
              'market_clearing_residuals', 'market_clearing_by_market',
              'tenure_services', 'owner_rungs', 'ownership_by_age_income_state',
              'housing_by_age_income_state', 'liquid_wealth_by_age_income_state',
              'fertility_policy_by_age_income_state', 'income_state_outcomes',
              'policy_childless_renter_age30', 'wealth_dist_childless_renter_age30',
              'policy_childless_renter_age42', 'wealth_dist_childless_renter_age42']
    for i, stem in enumerate(graphs):
        if i % 2 == 0:
            story.append(PageBreak())
            add('Native initial diagnostics', 'SectionCustom')
            add('The unchanged 17-figure packet from the verified initial solution. '
                'These are initial model diagnostics, not fitted 2023 comparisons. '
                'Native graph statistics can use different caps or age weights from '
                'the calibration observers.', 'SmallCustom')
        path = initial / 'standard_diagnostics' / (stem + '.png')
        iw, ih = ImageReader(str(path)).getSize()
        scale = min(width / iw, 274 / ih)
        story.append(Image(str(path), width=iw * scale, height=ih * scale, hAlign='CENTER'))
        add(escape(stem.replace('_', ' ')), 'SmallCustom')
        story.append(Spacer(1, 8))

    def footer(canvas, doc):
        canvas.setStrokeColor(gray)
        canvas.line(40, 35, A4[0] - 40, 35)
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(navy)
        canvas.drawString(40, 23, 'Fertility and housing | Verified initial readout')
        canvas.drawRightString(A4[0] - 40, 23, str(doc.page))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(args.output), pagesize=A4, leftMargin=40, rightMargin=40,
                            topMargin=38, bottomMargin=46,
                            title='Fertility and housing: verified initial readout',
                            author='Quantitative model working report')
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
