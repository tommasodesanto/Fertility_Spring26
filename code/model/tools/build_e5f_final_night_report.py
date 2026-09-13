#!/usr/bin/env python3
"""Build the verified initial readout and unchanged native figure appendix; no solve."""
from __future__ import annotations

import argparse
import csv
import hashlib
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
from pypdf import PdfReader, PdfWriter


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--as-of', required=True)
    ap.add_argument('--history-case', type=Path, default=None,
                    help='Optional verified four-window history packet')
    ap.add_argument('--policy-dir', type=Path, default=None,
                    help='Optional verified policy_comparison directory')
    ap.add_argument('--policy-state-verification', type=Path, default=None,
                    help='Optional passed policy state verification JSON')
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

    history_rows = None
    policy_rows = None
    policy_verification = None
    policy_state = None
    if args.policy_dir or args.policy_state_verification:
        if not (args.policy_dir and args.policy_state_verification and args.history_case):
            raise ValueError('--policy-dir, --policy-state-verification and --history-case are required together')
        policy_dir = args.policy_dir
        with (policy_dir / 'comparison.csv').open() as stream:
            policy_rows = list(csv.DictReader(stream))
        policy_verification = json.loads((policy_dir / 'verification.json').read_text())
        policy_state = json.loads(args.policy_state_verification.read_text())
        if policy_verification.get('status') != 'passed' or policy_verification.get('horizon_verified') is not False:
            raise ValueError('Policy comparison must be passed with horizon_verified=false')
        if policy_state.get('status') != 'passed' or policy_state.get('horizon_verified') is not False:
            raise ValueError('Policy state verification must be passed with horizon_verified=false')
        if len(policy_verification.get('years', [])) != 6 or policy_verification['years'] != [2023, 2027, 2031, 2035, 2039, 2043]:
            raise ValueError('Expected six policy dates')
        if len(policy_rows) != 66 or {int(r['decision_year']) for r in policy_rows} != {2023, 2027, 2031, 2035, 2039, 2043}:
            raise ValueError('Expected 6 x 11 policy comparison rows')
        if policy_state.get('policy_root_sha256', {}).keys() != {'baseline_rebate', 'tax2_rebate'}:
            raise ValueError('Both policy root provenance hashes are required')
        roots = {}
        for name in ('baseline_rebate', 'tax2_rebate'):
            root = args.history_case / 'policies' / name / 'root_receipt.json'
            if not root.exists() or hashlib.sha256(root.read_bytes()).hexdigest() != policy_state['policy_root_sha256'][name]:
                raise ValueError('Policy root provenance hash mismatch: ' + name)
            roots[name] = root
        source_sha = policy_verification.get('source_sha256', {})
        for path, digest in source_sha.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed policy comparison input: '+path)
        if source_sha.get(str(args.policy_state_verification.resolve())) != hashlib.sha256(args.policy_state_verification.read_bytes()).hexdigest():
            raise ValueError('Policy comparison must use this native state proof')
        for name in ('comparison.csv','policy_comparison.pdf','policy_comparison.png'):
            if policy_verification.get('output_sha256', {}).get(name) != hashlib.sha256((policy_dir/name).read_bytes()).hexdigest():
                raise ValueError('Policy comparison output hash mismatch: '+name)
        for root in roots.values():
            if str(root.resolve()) not in source_sha or source_sha[str(root.resolve())] != hashlib.sha256(root.read_bytes()).hexdigest():
                raise ValueError('Policy comparison source hash does not pin ' + root.name)
    if args.history_case:
        case = args.history_case
        realized = json.loads((case / 'realized_fit.json').read_text())
        finite = json.loads((case / 'finite_history_complete.json').read_text())
        verification = json.loads((case / 'readout' / 'verification.json').read_text())
        contract = json.loads((case / 'contract_receipt.json').read_text())
        if contract.get('case') != 'A0' or contract.get('count') != 6:
            raise ValueError('This report interpretation requires the A0 six-period case')
        if len(realized) != 4 or finite.get('realized') != realized:
            raise ValueError('History packet must contain the exact four-window fit list')
        if any(abs(float(row['gap'])) > .005 for row in realized):
            raise ValueError('History fertility fit gap exceeds 0.005')
        if verification.get('status') != 'PASS' or verification.get('native_2023_snapshot_sha256') is None:
            raise ValueError('Native history readout verification did not pass')
        if verification.get('historical_fit_status', {}).get('realized') != realized:
            raise ValueError('Readout must verify this complete realized history')
        provenance = json.loads((case / 'validation' / 'validation_2023_manifest.json').read_text())
        for key, path in [('model_source_sha256', case/'readout/model_2023.json'),
                          ('readout_verification_sha256', case/'readout/verification.json')]:
            if provenance[key] != hashlib.sha256(path.read_bytes()).hexdigest():
                raise ValueError('Validation source pin failed: '+key)
        if provenance.get('deterministic_fixture') is not False:
            raise ValueError('An actual native readout is required')
        validation = case / 'validation' / 'validation_2023.csv'
        with validation.open() as stream:
            history_rows = list(csv.DictReader(stream))
        if len(history_rows) != 13 or any(str(r.get('model_year')) != '2023' for r in history_rows):
            raise ValueError('Expected 13 rows for model year 2023')
        figure_pdf = case / 'figures' / 'e5f_final_history_figures.pdf'
        if not figure_pdf.exists():
            raise ValueError('Existing history figure PDF is required')

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

    add('Quantitative model<br/>' + ('Verified history and model fit' if args.history_case else 'Verified initial readout'), 'TitleCustom')
    add(escape(args.as_of) + (' | Six-period history and paired tax policies verified; horizon adequacy pending.' if policy_rows is not None else ' | Six-period history verified; policies and longer horizons continue.'
        if args.history_case else ' | Historical and policy work remains in progress.'), 'SmallCustom')
    add('<b>The fiscal and distribution repairs pass. The economic fit remains weak.</b> '
        'The corrected initial equilibrium reproduces exactly in two independent numerical '
        'repetitions. Its objective is %.6f. This is a verified candidate from an incomplete '
        'search, not a converged calibration optimizer.' % candidate['loss'])
    table([
        ['Question', 'Verified evidence'],
        ['Does the pension budget balance?', 'Yes in the initial equilibrium and every accepted historical forecast date.' if args.history_case else 'Yes in the initial equilibrium and accepted first historical forecasts. Every later accepted date must pass the same PAYGO check.'],
        ['Is property tax rebated?', 'Yes. Revenue returns equally per current household head; the rebate budget is checked separately from pensions.'],
        ['Was the probability error fixed?', 'Yes. Float64 normalization in both Markov distribution operators passes the captured mass-error replay and full initial reconstruction.'],
        ['Is the entire shock path fitted?', 'Yes for the provisional six-period A0 history. All four observed fertility windows pass the 0.005 fit requirement. Horizon adequacy remains unverified.' if args.history_case else 'Not yet. Both short-horizon variants have accepted the first historical window; subsequent surprises are being fitted.'],
        ['Are policy results ready?', 'Yes for the six-period A0 case: baseline and doubled property tax both return revenue equally, balance pensions, and start from identical 2023 households.' if policy_rows is not None else 'Not yet. Baseline and higher-tax rebated policies follow each completed history from its inherited 2023 households.'],
    ], [158, width - 158])
    add('What deserves attention', 'SectionCustom')
    add('Mean rooms and ownership at ages 30-55 are the largest remaining fit problems. '
        'Annual beta reaches its 0.99 cap. The numerical repair changes the objective by '
        'less than one millionth; it is not evidence of improved economic fit.')
    add('The 6-, 24- and 100-period forecasts are separate horizon checks. Passing a finite '
        'market-clearing root does not establish that a shorter horizon is adequate. The '
        'alternative demographic closure remains unresolved and has not been promoted.')
    if policy_rows is not None:
        policy_map={(int(r['decision_year']),r['metric']):r for r in policy_rows}
        add('The paired tax exercise changes the first birth flow by %.2f%% and the final simulated '
            'birth flow by %+.2f%%. These are small, provisional effects. A0 assumes zero net migration '
            'after 2023; six periods cover 24 years, followed by a constant-condition boundary '
            'for remaining household lifetimes.' % (float(policy_map[(2023,'birth_children_topcode_adjusted')]['percent_change']),
            float(policy_map[(2043,'birth_children_topcode_adjusted')]['percent_change'])))
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
    table(rows, [201, 61, 61, 71, 65, width - 459], numeric_start=1)
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

    if history_rows is not None:
        story.append(PageBreak())
        add('Verified four-window history', 'TitleCustom')
        add('All four fertility windows are matched in a six-period forecast. Each accepted date '
            'passes housing, PAYGO and equal-rebate gates. Horizon adequacy remains unverified; '
            'longer forecasts continue. Matching the fertility decline does '
            'not resolve weak cross-sectional fit.')
        table([['Window end', 'Data', 'Model', 'Gap', 'psi']] +
              [[int(r['year'])+4, number(r['target']), number(r['model']), number(r['gap']), number(r['psi'])]
               for r in realized], [105, 95, 95, 95, width - 390], numeric_start=1)
        story.append(PageBreak())
        add('Full 2023 validation: Data / Model (untargeted)', 'TitleCustom')
        add('These comparisons are a validation readout; they are not additional fitted targets.')
        table([['Moment', 'Data', 'Model', 'Gap', 'Vintage']] +
              [[r['moment'], number(r['data']), number(r['model']), number(r['gap']), r['data_vintage']]
               for r in history_rows], [185, 63, 63, 63, width - 374], numeric_start=1)
        caveats = [r for r in history_rows if r.get('measurement_note')]
        if caveats:
            story.append(PageBreak())
            add('2023 validation measurement caveats', 'TitleCustom')
            for r in caveats:
                add('<b>%s.</b> %s' % (escape(r['moment']), escape(r['measurement_note'])))

    if policy_rows is not None:
        story.append(PageBreak())
        add('Policy comparison: 2023–2043', 'TitleCustom')
        add('A0 means no net migration after 2023. The six periods span 24 years, with a remaining-lifetime constant-condition boundary. The paired reforms apply a 1% annual property tax versus a 2% annual property tax, with equal rebates. Both start from the same inherited 2023 households and balance PAYGO separately. Births are per-period flows, not cumulative; the last decision in 2043 corresponds to the 2044–47 flow. These are provisional short-horizon results because horizon adequacy is unverified.')
        metrics = ['asset_price','renter_price','housing_demand','resident_persons','household_heads','birth_children_topcode_adjusted','owner_rate','pension_period_units','equal_transfer_period_units','period_tfr_topcode_adjusted','rooms_per_head']
        labels = {'asset_price':'Asset price','renter_price':'Renter price','housing_demand':'Housing demand','resident_persons':'Resident persons','household_heads':'Household heads','birth_children_topcode_adjusted':'Births (flow)','owner_rate':'Owner rate','pension_period_units':'Pension (period units)','equal_transfer_period_units':'Equal transfer (period units)','period_tfr_topcode_adjusted':'Period TFR (flow)','rooms_per_head':'Rooms per head'}
        by = {(int(r['decision_year']), r['metric']): r for r in policy_rows}
        for year in (2023, 2043):
            if year == 2043:
                story.append(PageBreak())
            add('Decision year %d' % year, 'SectionCustom')
            rows = [['Metric','Baseline','2% tax','Percent change','Change (pp)']]
            for metric in metrics:
                r = by[(year, metric)]
                pp = r.get('percentage_point_change','')
                rows.append([labels[metric], number(r['baseline']), number(r['reform']), number(r['percent_change']) + '%', (number(pp) + ' pp') if metric == 'owner_rate' else ''])
            table(rows, [150, 78, 78, 92, width-398], numeric_start=1)
        add('Housing supply uses the same price-dependent curve in both policies, with elasticity '
            '0.63. The lower house price reduces supply, so the tax experiment includes a contraction '
            'in total housing as well as changes in its allocation. Reported rooms per head are '
            'physical rooms; the empirical validation table uses its separate top-coded room measure.')

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

    policy_figure_paths = []
    if policy_rows is not None:
        for policy_name, label in (('baseline_rebate', 'Baseline policy'), ('tax2_rebate', '2% tax policy')):
            diag = args.history_case / 'policies' / policy_name / 'graphs' / 'standard_diagnostics'
            for stem in graphs:
                path = diag / (stem + '.png')
                if not path.exists():
                    raise ValueError('Missing policy diagnostic: ' + str(path))
                policy_figure_paths.append((label, stem, path))
        for i, (label, stem, path) in enumerate(policy_figure_paths):
            if (i % len(graphs)) % 2 == 0:
                story.append(PageBreak())
                add('Policy native diagnostics: ' + label, 'SectionCustom')
                add('The unchanged 17-figure standard diagnostic packet for this policy.', 'SmallCustom')
            iw, ih = ImageReader(str(path)).getSize()
            scale = min(width / iw, 274 / ih)
            story.append(Image(str(path), width=iw * scale, height=ih * scale, hAlign='CENTER'))
            add(escape(label + ' — ' + stem.replace('_', ' ')), 'SmallCustom')
            story.append(Spacer(1, 8))

    def footer(canvas, doc):
        canvas.setStrokeColor(gray)
        canvas.line(40, 35, A4[0] - 40, 35)
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(navy)
        canvas.drawString(40, 23, 'Fertility and housing | ' + ('Verified history and model fit' if args.history_case else 'Verified initial readout'))
        canvas.drawRightString(A4[0] - 40, 23, str(doc.page))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(args.output), pagesize=A4, leftMargin=40, rightMargin=40,
                            topMargin=38, bottomMargin=46,
                            title='Fertility and housing: '+('verified history and model fit' if args.history_case else 'verified initial readout'),
                            author='Quantitative model working report')
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    if args.history_case:
        writer = PdfWriter()
        for page in PdfReader(str(args.output)).pages:
            writer.add_page(page)
        for page in PdfReader(str(args.history_case / 'figures' / 'e5f_final_history_figures.pdf')).pages:
            writer.add_page(page)
        if policy_rows is not None:
            for page in PdfReader(str(args.policy_dir / 'policy_comparison.pdf')).pages:
                writer.add_page(page)
        writer.add_metadata({'/Title': 'Fertility and housing: verified history and model fit',
                             '/Author': 'Quantitative model working report'})
        with args.output.open('wb') as stream:
            writer.write(stream)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
