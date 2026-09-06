#!/usr/bin/env python3
"""Read-only discussion PDF from collected simultaneous-choice verification."""
import argparse, hashlib, json, math
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Spacer, PageBreak, SimpleDocTemplate
import build_e5f_overnight_numerical_report as layout


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fmt(x):
    value = float(x)
    return f'{value:.6g}' if math.isfinite(value) else '-'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--history-results', type=Path, required=True)
    ap.add_argument('--policy-results', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    hist, policies = a.history_results.resolve(), a.policy_results.resolve()
    anchor = hist/'smoke/smoke_anchor/task_001'
    fit = layout.read_csv(anchor/'target_fit_long.csv')
    pars = layout.read_csv(anchor/'parameter_table.csv')
    summary = layout.read_json(anchor/'summary.json')
    cross = layout.read_json(hist/'smoke/cross_snapshot_verification.json')
    smoke = layout.read_json(hist/'smoke/smoke_verification.json')
    verification = layout.read_json(policies/'equilibrium_receipt.json')
    handoff = layout.read_json(policies/'inherited_state_verification.json')
    reporting = layout.read_json(hist/'reporting_after/reporting_check.json')
    reporting_before = layout.read_json(hist/'reporting_before/reporting_check.json')
    if (reporting['status'] != 'pass' or reporting['exact_unchanged_arrays'] != 15 or
        reporting['exact_standard_graphs'] != 17 or
        reporting['budget']['budget_excess_mass'] > 2e-10 or
        reporting['reference_sha256'] != digest(hist/'reporting_before/reporting_check.json') or
        reporting['array_hashes'] != reporting_before['array_hashes'] or
        reporting['graph_hashes'] != reporting_before['graph_hashes'] or
        reporting['quantities'] != reporting_before['quantities']):
        raise RuntimeError('Owner-consumption reporting repair is not independently verified')
    if cross['status'] != 'pass' or smoke['status'] != 'pass':
        raise RuntimeError('Verified historical repetitions are required')
    expected_policies = {'baseline','supply-plus-20','dependent-child-ltv95','property-tax-2pct-no-rebate'}
    if (not verification['smoke'] or verification['status'] != 'complete' or verification['failures'] or
        set(verification['cases']) != expected_policies or handoff['status'] != 'exact_feasibility_replay'):
        raise RuntimeError('This report requires the verified two-date policy experiment')
    for result in verification['cases'].values():
        if (result['status'] != 'complete' or result['dates'] != 2 or
            result['gates']['maximum_market_residual'] > 2e-4 or
            result['gates']['maximum_mass_residual'] > 2e-10):
            raise RuntimeError('Incomplete or invalid policy-path gates')
    if verification['selected_summary_sha256'] != digest(anchor/'summary.json'):
        raise RuntimeError('Policy paths belong to a different anchor')
    if digest(policies/'inherited_state_verification.json') != verification['inherited_state_verification_sha256']:
        raise RuntimeError('Policy handoff receipt changed')
    loss = float(summary['best_candidate']['transition_loss'])
    if len(fit) != 12 or sum(p['is_free_parameter'].lower() == 'true' for p in pars) != 11:
        raise RuntimeError('Changed moment or free-parameter count')
    for row in fit:
        gap = float(row['model'])-float(row['target'])
        if abs(gap-float(row['gap'])) > 1e-12 or not math.isclose(gap*gap*float(row['weight']),float(row['loss_contribution']),abs_tol=1e-10):
            raise RuntimeError('Inconsistent target-fit arithmetic')
    if abs(sum(float(x['loss_contribution']) for x in fit)-loss) > 1e-9:
        raise RuntimeError('Objective does not sum to full fit table')
    receipt = layout.read_json(anchor/'case_receipt.json')
    checked = 0
    for rel, sha in receipt['artifact_sha256'].items():
        path = anchor/rel
        if not path.exists() or digest(path) != sha:
            raise RuntimeError(f'Missing or changed selected artifact: {path}')
        checked += 1

    fonts = Path('/System/Library/Fonts/Supplemental')
    pdfmetrics.registerFont(TTFont('Review', str(fonts/'Arial.ttf')))
    pdfmetrics.registerFont(TTFont('ReviewBold', str(fonts/'Arial Bold.ttf')))
    pdfmetrics.registerFontFamily('Review',normal='Review',bold='ReviewBold')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('ReviewTitle',fontName='ReviewBold',fontSize=23,leading=28,textColor=layout.BLUE,spaceAfter=12))
    styles.add(ParagraphStyle('ReviewHead',fontName='ReviewBold',fontSize=15,leading=19,textColor=layout.BLUE,spaceAfter=10))
    styles.add(ParagraphStyle('ReviewBody',fontName='Review',fontSize=10,leading=14,spaceAfter=9))
    styles.add(ParagraphStyle('ReviewSmall',fontName='Review',fontSize=8,leading=10.5,spaceAfter=6))
    styles.add(ParagraphStyle('ReviewTableHead',fontName='ReviewBold',fontSize=8,leading=10.5,textColor=layout.colors.white))
    story = []
    def p(text, kind='ReviewBody'):
        return Paragraph(text,styles[kind])
    def add(text, kind='ReviewBody'):
        story.append(p(text,kind))
    def heading(text):
        story.append(PageBreak());add(text,'ReviewHead')
    def table(rows,widths):
        body = [[p(escape(str(x)).replace('\n','<br/>'),'ReviewTableHead' if i==0 else 'ReviewSmall') for x in row] for i,row in enumerate(rows)]
        story.append(layout.table(body,widths,font_size=8))
    moments = {r['moment']:r for r in fit}
    own = moments['own_family_gap']; birth = moments['housing_increment_0to1']
    status = verification['status']
    add('Simultaneous choice:<br/>review before calibration','ReviewTitle')
    add('6 September 2026 | Experimental specification | Full overnight search stopped','ReviewSmall')
    add('<b>The lifecycle calculation is reproducible. A successful recalibration has not yet been established.</b> '
        'Two full histories reproduce all twelve target rows, all parameters, 253 historical entries and seventeen standard graphs exactly. '
        'The checkpoint and consumption-reporting repairs preserve the economic choices and fit.')
    add(f'<b>The starting point fits poorly.</b> Its objective is {loss:.3f}. The ownership gap between parents and nonparents is '
        f'{100*float(own["model"]):.3f} percentage points, against {100*float(own["target"]):.3f} in the data. '
        f'The first-birth housing response is {float(birth["model"]):.3f} rooms, against {float(birth["target"]):.3f}. '
        'These are values at starting parameters, not the result of a calibration search.')
    add('<b>The economic issue to discuss is the nesting restriction.</b> The experiment groups waiting and attempting a birth '
        'within each tenure. The outer shock scale is κ; the inner scale is λκ, with 0 &lt; λ ≤ 1. '
        'For example, κ = 2 and λ = 0.8 imply an inner scale of 1.6. '
        'This links tenure and fertility dispersion. Both must be re-estimated jointly with the other nine parameters.')
    add('<b>Simultaneous tastes still require a conception contract.</b> The household chooses tenure and whether to attempt a birth '
        'together, after observing the joint tastes. Conception succeeds afterward. Housing size, consumption and saving may adjust '
        'within the chosen tenure. The same nesting parameter applies across birth orders; this is an experimental restriction.')
    add(f'<b>Policy-stage verification:</b> {escape(status.replace("_"," "))}. '
        'The experiment checks 2023 and 2027 for the baseline, additional supply, relaxed parental LTV and higher property tax. '
        'These short paths test the workflow; they are not calibrated policy estimates or a 2063 forecast.')
    add('<b>Proposed next step after our discussion:</b> search all eleven parameters against the unchanged twelve targets and weights; '
        'reserve time for a local sensitivity matrix, two exact final repetitions and full policy paths. '
        'Having twelve moments for eleven parameters is not an identification test. Retain the current expectation method for this run and discuss '
        'perfect foresight separately. No production replacement has been made.')

    heading('Complete target fit at the starting point')
    add('Every original target is retained. Shares and share gaps are fractions: 0.01 equals one percentage point. '
        'The objective is the sum of weight × gap²; the weights do not all represent estimated sampling precision.','ReviewSmall')
    labels = dict(layout.MOMENT_LABELS)
    labels.update(own_family_gap='Parent/nonparent ownership gap',own_rate='Prime-age ownership',
                  prime30_55_parent_3plus_minus_1to2_mean_rooms='Rooms: 3+ minus 1-2 children')
    rows = [['Moment','Target','Model','Gap','Weight','Loss']]
    rows += [[labels[r['moment']],*[fmt(r[k]) for k in ('target','model','gap','weight','loss_contribution')]] for r in fit]
    table(rows,[181,62,62,66,75,77])
    story.append(Spacer(1,10));add(f'<b>Total objective: {loss:.9f}.</b> The starting ownership-gap miss contributes '
        f'{float(own["loss_contribution"]):.3f}, or {100*float(own["loss_contribution"])/loss:.1f}% of the total. '
        'A small fit improvement in the verification probes is not a searched optimum.','ReviewSmall')

    heading('Parameters, search bounds and fixed restrictions')
    add('Starting values of the eleven free parameters. The two shock parameters are experimental starting choices; '
        'the remaining starting values come from the inherited candidate. None is presented as a newly estimated optimum.','ReviewSmall')
    labels = {'beta_annual':'Annual discount factor','tenure_choice_kappa':'Outer tenure scale κ','joint_nest_lambda':'Dissimilarity λ',
              'chi':'Owner housing-service premium χ','H0':'Housing-supply normalization H0','theta0':'Bequest utility weight θ0',
              'theta1':'Bequest wealth shift θ1','hbar_child_rooms':'Per-child room floor','first_birth_fixed_cost':'First-birth utility cost',
              'hbar_first_child_jump':'First-child room-floor jump','psi_child_change_2023':'2007-2023 child-value change'}
    free = [r for r in pars if r['is_free_parameter'].lower()=='true']
    table([['Free parameter','Starting value','Lower','Upper','Near bound?']]+[
        [labels[r['parameter']],fmt(r['value']),fmt(r['lower_bound']),fmt(r['upper_bound']),r['near_bound']] for r in free],
        [221,84,69,69,80])
    story.append(Spacer(1,12))
    fixed = [r for r in pars if r not in free]
    table([['Other object','Value','Restriction']]+[[
        {'psi_child_2007':'Child value in 2007','psi_child_2023':'Child value in 2023','housing_supply_elasticity':'Supply elasticity'}[r['parameter']],
        fmt(r['value']),{'psi_child_2007':'Normalized to old completed fertility = 2.1','psi_child_2023':'Derived from 2007 level plus fitted change','housing_supply_elasticity':'Externally fixed'}[r['parameter']]] for r in fixed], [190,85,248])
    story.append(Spacer(1,10));add('Near-bound flags use the existing convention: within 2% of the full physical parameter range. '
        'Thus a flag for a log-transformed parameter need not imply a small distance in its search coordinate. '
        'The lower restriction λ ≥ 0.02 is numerical; the upper restriction λ ≤ 1 comes from the chosen nested-GEV law.','ReviewSmall')

    heading('Standard policy diagnostics: starting-point weakness')
    add('The ownership gradient is weak at these starting scales. This is economically important despite passing probability, '
        'market-clearing and occupied-value checks. It does not show that the target is unreachable after recalibration.','ReviewSmall')
    for name in ('ownership_by_age_income_state.png','housing_by_age_income_state.png'):
        story.append(layout.image_fit(anchor/'standard_diagnostics'/name,523,282));story.append(Spacer(1,10))

    heading('Policy-loop check and limits of the current evidence')
    effects = layout.read_csv(policies/'policy_effects.csv')
    if len(effects) != 6 or {(r['policy'],int(r['year'])) for r in effects} != {
        (name,year) for name in expected_policies-{'baseline'} for year in (2023,2027)}:
        raise RuntimeError('Missing or duplicated policy-effect rows')
    names = {'supply-plus-20':'Housing supply +20%','dependent-child-ltv95':'Dependent-child LTV 95%',
             'property-tax-2pct-no-rebate':'Property tax doubled'}
    add('The table below is a workflow check at uncalibrated starting parameters. Each row compares the policy against its '
        'same-date baseline, starting from the same original 2023 population. Do not put these values into the presentation as policy estimates.','ReviewSmall')
    table([['Policy','Year','Births (%)','Ownership (pp)','Rooms (%)']]+[
        [names[r['policy']],r['year'],fmt(r['births_percent']),fmt(r['ownership_pp']),fmt(r['rooms_percent'])] for r in effects],
        [195,48,88,102,90])
    story.append(Spacer(1,12))
    add('<b>Maintained closures.</b> Markets clear separately at each date; households treat each current price as permanent. '
        'This is temporary equilibrium, not perfect foresight. After 2023 the national population is closed: no outside entry, '
        'retention one, the inherited four-slot birth queue and births converted into entrant households by 1/2.1. '
        'The property-tax experiment discards revenue and has no rebate or grant. These are finite-horizon paths, not stationary endpoints.','ReviewSmall')
    add('<b>Verified population handoff.</b> The original inherited population is saved separately from its price-specific feasibility '
        f'adjustment. The fitted gate replays exactly; adjusted mass is {handoff["fitted_projection_mass"]:.4g}. '
        'The queue entering 2023 comes from the completed 2019 row, not the already advanced 2023 row.','ReviewSmall')
    add('<b>Verified reporting repair.</b> A legacy owner-consumption output floor reported 0.040 where one budget supported 0.015366. '
        'The experimental output now reports consumption implied by the optimizer. At the exact failed policy price, all fifteen '
        'value, choice and distribution arrays, market quantities and seventeen graphs are unchanged. Budget-violating mass falls '
        f'from {reporting_before["budget"]["budget_excess_mass"]:.3g} to {reporting["budget"]["budget_excess_mass"]:.3g}, '
        'below the unchanged 2e-10 gate.','ReviewSmall')
    add('<b>Still outstanding.</b> Full calibration; local identification and parameter trade-offs; exact final repetitions of a '
        'selected optimum; calibrated 2023-2063 policy paths; and author adoption of the nesting/commitment restrictions. '
        'No completed smoke settles those questions.','ReviewSmall')
    add('Evidence: the adjacent experiment README indexes both smoke revisions, immutable contracts, numerical checks, full tables '
        'and all seventeen unchanged standard graphs. Rebuild this PDF with build_e5f_joint_nested_review.py and the same collected '
        'history/policy directories.','ReviewSmall')

    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('Review',7);canvas.setFillColor(layout.MID_GREY)
        canvas.drawString(36,20,'Simultaneous-choice experiment | Discussion copy | Production unchanged')
        canvas.drawRightString(A4[0]-36,20,str(doc.page));canvas.restoreState()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    SimpleDocTemplate(str(a.output),pagesize=A4,rightMargin=36,leftMargin=36,topMargin=36,bottomMargin=36).build(story,onFirstPage=footer,onLaterPages=footer)
    checks = dict(status='numerical_source_checks_passed_visual_review_pending',pdf=str(a.output.resolve()),pdf_sha256=digest(a.output),
                  fit_rows=len(fit),free_parameters=len(free),verified_selected_artifacts=checked,loss=loss,
                  history_results=str(hist),policy_results=str(policies),builder_sha256=digest(__file__),production_promoted=False,
                  reporting_check_sha256=digest(hist/'reporting_after/reporting_check.json'))
    (a.output.parent/(a.output.stem+'_verification.json')).write_text(json.dumps(checks,indent=2)+'\n')
    print(a.output)


if __name__=='__main__':
    main()
