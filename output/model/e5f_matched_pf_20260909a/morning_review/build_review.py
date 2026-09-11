#!/usr/bin/env python3
"""Refreshable four-page initial-equilibrium working review; no model solves."""
from pathlib import Path
import csv, hashlib, json, math
from datetime import datetime, timezone
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.enums import TA_RIGHT
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
INPUT=HERE.parent/'initial_fit_readout'
STATUS=json.loads((HERE/'report_status.json').read_text())
summary=json.loads((INPUT/'summary.json').read_text())
smoke=json.loads((INPUT/'source_initial_smoke_summary.json').read_text())
obs=json.loads((INPUT/'source_initial_measurement_summary.json').read_text())
def readcsv(name): return list(csv.DictReader((INPUT/name).open()))
fit=readcsv('target_fit.csv'); valid=readcsv('validation_fit.csv')
params=readcsv('parameters.csv'); cps=readcsv('cps_projection_sensitivity.csv')
assert len(fit)==13 and len(valid)==2 and len(params)==17 and len(cps)==4
assert summary['actual_weight_profile'] is None and summary['summed_loss'] is None
assert not summary['calibrated_smm'] and not summary['parameter_search_performed']
for row in fit+valid+cps:
    assert not row['actual_weight'] and not row['loss_contribution']
    if row['model']:
        assert math.isclose(float(row['model'])-float(row['target']),float(row['gap']),abs_tol=1e-12)
assert [r['restriction_id'] for r in fit if not r['model']]==['recent_parent_ownership']
for entry in summary['source_snapshots'].values():
    p=Path(entry['durable_path']); assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
assert obs['checkpoint_sha256']==summary['checkpoint_sha256']==smoke['final']['checkpoint_sha256']
assert smoke['status']=='passed_initial_candidate_loop' and obs['fiscal']['fiscal_gate']
assert not smoke['perfect_foresight_solved']
assert smoke['final']['household_budget']['budget_excess_mass']==0
assert smoke['final']['policy_array_gates']['occupied_negative_steps']==0

FONT=Path('/System/Library/Fonts/Supplemental')
for name,file in [('Body','Arial.ttf'),('Bold','Arial Bold.ttf'),('Italic','Arial Italic.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(FONT/file)))
pdfmetrics.registerFontFamily('Body',normal='Body',bold='Bold',italic='Italic',boldItalic='Bold')
ink=colors.HexColor('#20354B'); gray=colors.HexColor('#5C6269'); pale=colors.HexColor('#EDF2F6')
styles={
 'title':ParagraphStyle('title',fontName='Bold',fontSize=21,leading=24,textColor=ink,spaceAfter=12),
 'sub':ParagraphStyle('sub',fontName='Body',fontSize=10,leading=13,textColor=gray,spaceAfter=13),
 'head':ParagraphStyle('head',fontName='Bold',fontSize=12,leading=15,textColor=ink,spaceBefore=11,spaceAfter=6),
 'body':ParagraphStyle('body',fontName='Body',fontSize=10.7,leading=14.5,spaceAfter=8),
 'small':ParagraphStyle('small',fontName='Body',fontSize=9,leading=12,spaceAfter=6),
 'cell':ParagraphStyle('cell',fontName='Body',fontSize=8.8,leading=11),
 'num':ParagraphStyle('num',fontName='Body',fontSize=8.6,leading=11,alignment=TA_RIGHT),
}
W=518
story=[]
def P(t,style='body'): return Paragraph(t,styles[style])
def add(t,style='body'): story.append(P(t,style))
def head(t): add(t,'head')
def title(t,subtitle): add(t,'title');add(subtitle,'sub')
def fmt(x):
    if x is None or x=='':return '—'
    v=float(x)
    return f'{v:.6f}' if abs(v)>=.0001 or v==0 else f'{v:.4g}'
def table(headers,rows,widths,numeric=()):
    data=[[P(escape(str(v)),'cell') for v in headers]]
    for row in rows:data.append([P(escape(str(v)),'num' if i in numeric else 'cell') for i,v in enumerate(row)])
    t=Table(data,colWidths=widths,hAlign='LEFT',repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),pale),('LINEBELOW',(0,0),(-1,0),.6,ink),
       ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),
       ('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),
       ('LINEBELOW',(0,-1),(-1,-1),.5,gray)]))
    story.append(t);story.append(Spacer(1,8))
F={r['restriction_id']:r for r in fit};V={r['restriction_id']:r for r in valid}
final=smoke['final']; fiscal=obs['fiscal']['actual_accounts']; generated=datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')

# PAGE 1
title('Initial equilibrium: working review',f"Evidence as of {STATUS['as_of']} • Generated {generated}<br/>New parenthood-only housing requirement and balanced initial pensions")
add('<b>The new initial economy is numerically verified. It is not yet calibrated.</b> The nine structural coordinates are inherited or mapped starting values. Only the initial fertility-preference level has been adjusted to its separate normalization. No active weight profile or total objective exists for the proposed early target system.')
head('What changed economically')
add('Housing needs now depend on whether any dependent child is at home: a common parenthood requirement applies while dependents are present, with no per-child slope. The external equivalence scale, consumption share, curvature and linear child reward remain. The approved choice structure is sequential.')
add('Initial pensions now balance against actual payroll earnings and retiree exposure. The payroll-tax rate is externally fixed; the pension benefit is derived from the budget. This verifies the initial stationary fiscal calculation. A transition requires a separately balanced benefit at every date, entering both household expectations and realized budgets.')
head('What is established')
table(['Check','Verified evidence'],[
 ['Initial calculation',f"{smoke['repetitions']} fresh repeated normalizations; {smoke['stationary_solves']} stationary solves; same checkpoint read by the observer."],
 ['Markets and pensions',f"Market error {final['legacy_stationary_moments']['market_residual']:.3g}; absolute scaled pension imbalance {abs(fiscal['scaled_pension_budget_residual']):.3g} (gate 1e-6)."],
 ['Household and operator checks','Zero recorded household budget-violating mass and occupied negative value steps; probability, population and checkpoint checks pass.'],
 ['Numerical tests',f"{STATUS['compiled_tests_passed']} compiled utility/pension/anticipation tests pass; {STATUS['pure_history_tests_passed']} combined pure historical-adapter tests pass."],
 ['Sensitivity work',f"{STATUS['sensitivity_cases_verified']} of {STATUS['sensitivity_cases_planned']} initial cases verified at the stated checkpoint; no completed case failed. Full parameter Jacobian remains pending."],
], [132,386])
head('The economic fit that deserves attention')
add(f"At this starting point, mean occupied rooms are <b>{float(F['mean_rooms']['model']):.3f} versus {float(F['mean_rooms']['target']):.3f}</b>, while ownership at ages 30–55 is <b>{100*float(F['ownership_30_55']['model']):.2f}% versus {100*float(F['ownership_30_55']['target']):.2f}%</b>. The first-birth room response is <b>{float(F['first_birth_rooms']['model']):.3f} versus {float(F['first_birth_rooms']['target']):.3f}</b>. Excess average space coexists with insufficient ownership and a weak family-space response. This is a diagnosis of the mapped point, not a claim that these gaps are unavoidable.")
add('<b>The 2.1 restriction is model completed fertility.</b> It is not a measured female period total fertility rate, and its successful normalization does not establish the other fertility fits. No new perfect-foresight history or matched policy is certified by this packet.','small')
story.append(PageBreak())

# PAGE 2
title('Complete initial fit', '13 restrictions: one separate normalization and twelve proposed scored rows')
add('Gap = model − target. Shares are fractions, ages are years and housing is rooms. A dash denotes unavailable, not zero. Every actual weight and loss contribution is null; reference precisions in the source are not adopted weights.','small')
labels={
 'initial_completed_fertility':'Initial completed fertility (normalization)',
 'cps_childlessness':'Childless women, 40–44', 'cps_exactly_one':'One child among mothers, 40–44',
 'mean_age_first_birth':'Period mean first-birth age','share_first_births_age30plus':'First births at 30+',
 'wealth_earnings':'Wealth / annual gross labor earnings','bequest_wealth':'Annual bequests / wealth',
 'old_dispersion':'Old wealth/income p90 / median, 76–84','mean_rooms':'Mean rooms, capped at 9',
 'ownership_30_55':'Ownership, heads 30–55','first_birth_rooms':'First-birth rooms response, −1 to +3',
 'family_rooms':'Rooms: 3+ versus 1–2 resident children¹','recent_parent_ownership':'Recent-parent ownership gap²'}
rows=[[labels.get(r['restriction_id'],r['label']),fmt(r['target']),fmt(r['model']),fmt(r['gap']),fmt(r['actual_weight']),fmt(r['loss_contribution'])] for r in fit]
table(['Restriction','Target','Model','Gap','Weight','Loss'],rows,[217,63,63,63,56,56],(1,2,3,4,5))
add('¹ Model value uses dependent counts as an explicit proxy for resident-child groups. ² Exact ACS recent-parent groups cannot be recovered by a static mask of current counts. A passive realized-birth/empty-home observer is under development; this is not an intrinsic impossibility result.','small')
head('Validation only')
table(['Observation','Target','Model','Gap','Weight / loss'],[[r['label'],fmt(r['target']),fmt(r['model']),fmt(r['gap']),'— / —'] for r in valid],[260,66,66,66,60],(1,2,3,4))
head('CPS age-projection sensitivity')
rows=[]
for key,label in [('cps_childlessness','Childlessness'),('cps_exactly_one','One child among mothers')]:
    a=next(r for r in cps if r['restriction_id']==key and r['projection']=='uniform_birth_time')
    b=next(r for r in cps if r['restriction_id']==key and r['projection']=='constant_post_cell')
    rows.append([label,fmt(a['target']),fmt(a['model']),fmt(a['gap']),fmt(b['model']),fmt(b['gap'])])
table(['Observation','Target','Uniform birth','Gap','Post-cell fixed','Gap'],rows,[170,65,76,65,77,65],(1,2,3,4,5))
add('Both approximate ages 40–44 within four-year cells. The one-child gap changes sign; first-birth flow timing is unchanged. Neither stock projection is promoted to an empirical weight or a certified exact age mapping.','small')
story.append(PageBreak())

# PAGE 3
title('Every coordinate and restriction', 'Nine structural starting coordinates; eight fixed, normalized or derived rows')
add('The values below are copied from the same checkpoint receipt. They are not new structural estimates. The near-bound flags are the supplied flags, not a new criterion and not evidence by themselves of identification. Full machine precision remains in parameters.csv.')
name={'beta_annual':'Annual discount factor (beta)','kappa_fert':'First-attempt taste scale','kappa_fert_continuation':'Later-attempt taste scale','chi':'Owner-service premium (chi)','H0':'Initial supply scale (H0)','theta0':'Bequest parameter (theta0)','theta1':'Bequest wealth shift (theta1)','first_birth_fixed_cost':'First-birth utility cost','h_P':'Parenthood housing requirement'}
searched=params[:9]
table(['Structural coordinate','Value','Lower','Upper','Transform','Near bound'],[[name[r['parameter']],fmt(r['estimate']),fmt(r['lower']),fmt(r['upper']),r['transform'],'Yes' if r['near_bound']=='True' else 'No'] for r in searched],[209,75,57,57,60,60],(1,2,3))
head('Fixed, normalized and budget-derived quantities')
fixednames={'hbar_child_rooms':'Per-dependent-child housing slope','psi_child':'Initial child-preference intercept','payroll_tax':'Payroll-tax rate','pension_period':'Pension benefit, four-year units','housing_supply_elasticity':'Housing-supply elasticity','tenure_choice_kappa':'Tenure-choice taste scale','alpha_cons':'Consumption share','sigma':'Utility curvature'}
table(['Restriction','Value','Status'],[[fixednames[r['parameter']],fmt(r['estimate']),r['status']] for r in params[9:]],[248,89,181],(1,))
head('How to interpret the starting values')
add('The parenthood requirement combines the old first-child jump and per-child floor; the slope is then fixed to zero. Initial supply was mapped under the common 0.63 elasticity. These are starting-value mappings, not estimates of the revised specification.')
add('The preference intercept is normalized to 2.1 model completed fertility. Pension equals revenue divided by retiree benefit exposure. A changing population requires renewed fiscal balance checks.')
add('The bequest wealth shift retains its near-bound flag. The complete sensitivity panel is needed to assess rank and weak parameter combinations. Counting restrictions does not establish identification.','small')
story.append(PageBreak())

# PAGE 4
title('What must be resolved next', 'Measurement, identification and transition certification are distinct tasks')
head('Measurement qualifications that affect interpretation')
add('<b>Recent parenthood.</b> ACS requires every resident own child to be younger than four, with no-resident-child controls that include former parents. The static count distribution lacks those exact age/residence labels. A passive observer can instead tag actual births into empty-dependent homes, transport that mass through realized housing choices, and compare with currently empty homes. Its timing and residence mapping must be declared and tested; it must not silently replace the empirical target.')
add('<b>Age, residence and geography.</b> Four-year age overlap and CPS parity timing affect measured stocks. Dependent counts are a proxy for ACS resident own children, including adult children. Model households are compared with 42 housing cities; structure filters and the revised city footprint are not exact model states. Rooms are capped at nine for ACS comparisons, whereas the reviewed PSID birth response remains uncapped.')
add('<b>Income and uncertainty.</b> Modeled retirement income is a proxy for PSID family income; the empirical income cutoff and selection are not fully reproduced. ACS errors are metro-resampling errors, not official survey-design errors. Pooled CPS covariance, NCHS discrepancy scales, cross-block covariance and the synthetic bequest scale require explicit treatment before weighting.')
add('<b>Birth-housing response.</b> Preserve the reviewed Sun–Abraham estimate and its stated sample/contrast. Applying it to an initial stationary economy assumes stability; model risk-set weighting and four-year timing are approximations. The existing prepath and normalization qualifications remain disclosed. This packet does not rerun or replace that regression.')
head('Next gates, in order')
for text in [
 'Finish the nine-coordinate sensitivity panel and verify every direction. Diagnose which parameter combinations move ownership, average rooms and family responses together; do not infer the full Jacobian from an incomplete panel.',
 'Complete and test the missing observation mapping; settle the early target/uncertainty contract. Only then score a joint structural candidate with an explicit objective and identification assessment.',
 'Verify the actual-population terminal price/pension loop and the joined historical perfect-foresight loop. Every dated benefit must enter backward decisions, forward accounting and a pension-balance gate.',
 'Establish horizon stability, then compare jointly recalibrated histories. Run matched property-tax policies only from a certified shared historical state and fiscal contract. Initial test success alone does not establish policy readiness.'
]:add('• '+text,'small')
head('Sources and refresh boundary')
add('Numerical tables: <b>initial_fit_readout/target_fit.csv, validation_fit.csv, cps_projection_sensitivity.csv and parameters.csv</b>. Their summary.json and the two preserved initial observation/smoke summaries identify the checkpoint and verified source snapshots. The early empirical contracts and recent_parent_flow_design.md supply definitions and the passive-observer qualification.','small')
add(f"Operational checkpoint: <b>{escape(STATUS['status_source'])}</b>. Source {escape(STATUS['source_commit'])}; receipt backup {escape(STATUS['receipt_commit'])}. Tables refer to the mapped initial checkpoint, not to partial sensitivity candidates or the newest source commit. This report reads local receipts; it does not query the cluster or presume queued jobs have finished.",'small')
add('Refresh with the supplied builder after updating the receipt tables and report_status.json. The builder checks row counts, null weights/losses, arithmetic, common checkpoint identity and source-snapshot hashes. Source fingerprints and output QA are written beside the PDF.','small')

class NumberedCanvas(canvas.Canvas):
    def __init__(self,*args,**kw): super().__init__(*args,**kw);self.saved=[]
    def showPage(self): self.saved.append(dict(self.__dict__));self._startPage()
    def save(self):
        n=len(self.saved)
        for state in self.saved:
            self.__dict__.update(state);self.setFont('Body',8);self.setFillColor(gray)
            self.drawString(47,26,'Working review • Initial diagnostic, not calibrated SMM')
            self.drawRightString(565,26,f'{self._pageNumber} / {n}');super().showPage()
        super().save()
pdf=HERE/'review.pdf'
SimpleDocTemplate(str(pdf),pagesize=(612,792),leftMargin=47,rightMargin=47,topMargin=39,bottomMargin=42,title='Initial equilibrium working review',author='Quantitative research review').build(story,canvasmaker=NumberedCanvas)
r=PdfReader(pdf);assert len(r.pages)==4,f'Expected four pages, got {len(r.pages)}'
text='\n'.join(p.extract_text() for p in r.pages)
for row in fit+valid:
    for key in ['target','model','gap']:
        if row[key]:assert fmt(row[key]) in text,(row['restriction_id'],key)
for row in params:assert fmt(row['estimate']) in text,row['parameter']
(HERE/'extracted.txt').write_text(text)
files=['target_fit.csv','validation_fit.csv','cps_projection_sensitivity.csv','parameters.csv','summary.json','source_initial_measurement_summary.json','source_initial_smoke_summary.json']
qa={'generated_utc':generated,'evidence_as_of':STATUS['as_of'],'pages':len(r.pages),'restriction_rows':13,'validation_rows':2,'parameter_rows':17,'projection_rows':4,'all_actual_weights_null':True,'all_actual_losses_null':True,'arithmetic_checked':True,'displayed_table_values_extraction_checked':True,'source_snapshots_verified':True,'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'input_sha256':{f:hashlib.sha256((INPUT/f).read_bytes()).hexdigest() for f in files},'visual_review':'PENDING render inspection'}
(HERE/'qa.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps({k:qa[k] for k in ['generated_utc','pages','pdf_sha256']}))
