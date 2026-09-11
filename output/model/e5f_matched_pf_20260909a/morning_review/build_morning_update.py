"""Four-page saved-result morning update; no model solves or new scientific plots."""
from pathlib import Path
import json, hashlib, math
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
SCORES=HERE.parent/'initial_calibration_contract/saved_case_scores'
baseline=json.loads((SCORES/'baseline_score.json').read_text())
joint=json.loads((SCORES/'joint_smoke_score.json').read_text())
assert baseline['contract_sha256']==joint['contract_sha256']
assert all(len(x['target_fit'])==13 and len(x['parameters'])==17 for x in (baseline,joint))
assert math.isclose(baseline['loss'],1499.85182523238,abs_tol=1e-10)
assert math.isclose(joint['loss'],1428.1716060194092,abs_tol=1e-10)
style={
 'title':ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=21,leading=25,spaceAfter=15,textColor=colors.HexColor('#20354B')),
 'head':ParagraphStyle('head',fontName='Helvetica-Bold',fontSize=12,leading=16,spaceBefore=10,spaceAfter=6),
 'body':ParagraphStyle('body',fontName='Helvetica',fontSize=10,leading=14,spaceAfter=10),
 'small':ParagraphStyle('small',fontName='Helvetica',fontSize=8,leading=11,spaceAfter=8),
 'cell':ParagraphStyle('cell',fontName='Helvetica',fontSize=8,leading=10),
 'num':ParagraphStyle('num',fontName='Helvetica',fontSize=8,leading=10,alignment=TA_RIGHT),
}
def para(t,k='body'):return Paragraph(t,style[k])
def table(rows,widths):
 data=[[para(str(x),'cell' if j==0 else 'num') for j,x in enumerate(row)] for row in rows]
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E8EEF3')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.5,colors.gray),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F6F8FA')])]))
 return t
labels={
'initial_normalization':'Initial fertility normalization','cps_childlessness':'Childlessness, women 40-44','cps_exactly_one':'Exactly one, among mothers 40-44','nchs_mean_age':'Mean age at first birth','nchs_share30':'First births at ages 30+','wealth_earnings':'Wealth / annual gross earnings','bequest_wealth':'Annual bequests / wealth','old_dispersion':'Old wealth-income ratio: p90 / median','mean_rooms':'Mean occupied rooms, capped at 9','ownership_30_55':'Ownership, household heads 30-55','first_birth_rooms':'First-birth room response','family_rooms':'Rooms: 3+ versus 1-2 dependents','recent_parent_ownership':'Recent-parent ownership gap'}
checked=[]
def fmt(v):return '-' if v is None else f'{v:.6g}'
def fitpage(x,title):
 story=[para(title,'title'),para('All 13 restrictions: 12 scored moments and a separate normalization. Proportions and ownership gaps are in fraction units; 0.01 equals one percentage point.','small')]
 rows=[['Moment','Target','Model','Gap','Weight','Loss']]
 for r in x['target_fit']:
  cells=[fmt(r[k]) for k in ('target','model','gap','actual_weight','loss_contribution')]
  checked.extend(z for z in cells if z!='-')
  rows.append([labels[r['restriction_id']]]+cells)
 story+=[table(rows,[171,61,61,61,69,69]),Spacer(1,12),para('Total working minimum-distance loss: <b>'+fmt(x['loss'])+'</b>. The fertility normalization has no weight or loss contribution.','small'),para('Weights are fixed working choices, not an efficient GMM covariance matrix. CPS uses approximate pooled scales; NCHS uses annual 2003-2006 variation as a typical-year discrepancy scale. Housing and wealth use retained bootstrap scales; bequests retain a synthetic 5% scale.','small'),para('Recent-parent ownership uses actual current births into previously empty-dependent homes, compared with current empty homes including former parents. The synchronized snapshot and dependent-residence approximation are maintained explicitly; this is not an exact annual ACS interview reconstruction.','small')]
 return story
story=[para('Morning quantitative update','title'),para('September 11, 2026 | saved results checked after the 09:36 EDT restart','small'),para('<b>Useful repairs and some large fit improvements, but no finished calibration or certified policy package.</b>')]
story += [para('What is established','head'),para('The approved parenthood-only housing floor and equivalence scale are implemented in the isolated sequential model. The initial steady state reproduces exactly and balances Social Security at the fixed 17.9% payroll tax. A separate terminal steady-state test also balances; its preference change remains diagnostic.'),para('The 19-case sensitivity panel finished successfully. The 24-proposal joint round produced 23 valid candidates, including the twice-reproduced smoke point; one candidate failed the unchanged strict housing-equilibrium check. No failed candidate was admitted.')]
story += [para('The economic result that matters','head'),para('The tested joint candidate raises the first-birth room response from 0.439 to 0.692, against 0.720 in the data. Mean rooms improve from 6.246 to 5.609, against 5.561; ownership at ages 30-55 rises from 57.48% to 63.10%, against 64.83%.'),para('However, its recent-parent ownership gap worsens from -5.32 to -6.34 percentage points, against +16.29 in the data. This row is retained. Once all 12 moments are counted, loss falls only 4.8%: 1499.85 to 1428.17. Recent-parent ownership accounts for 84.3% of baseline loss and 97.0% of the joint candidate loss. This candidate is not a selected final estimate.')]
story += [para('Transition and interruption','head'),para('The last collected six-date transition iterate reduced maximum housing imbalance from 45.3% to 1.76%, still above the 0.02% gate. Pension residuals were below 4.4e-10. The final result has not been collected: cluster authentication was unavailable at 09:36. This does not establish a completed perfect-foresight path, an estimated historical shock, or policy effects.'),para('Agent work stopped when credits ran out around 03:27. Previously submitted cluster jobs could finish independently; repeated automatic wake-ups are not evidence of further research. Work has resumed with one lead and hourly monitoring. Next: collect the transition, finish recent-parent measurement for all joint candidates, then continue the full-objective search after its exact-loop check.','small'),PageBreak()]
story+=fitpage(baseline,'Starting point: complete fit')+[PageBreak()]
story+=fitpage(joint,'Tested joint candidate: complete fit')+[PageBreak()]
story += [para('Parameters and remaining claims','title'),para('All 17 parameter/restriction rows. Nine structural coordinates vary jointly. The initial child-utility intercept is separately normalized to 2.1; fixed and budget-derived objects are not additional searched coordinates.','small')]
bp={r['parameter']:r for r in baseline['parameters']}
paramlabels={'beta_annual':'Annual discount factor','kappa_fert':'First-birth taste scale','kappa_fert_continuation':'Later-birth taste scale','chi':'Tenure taste parameter','H0':'Housing-supply scale','theta0':'Bequest level','theta1':'Bequest curvature','first_birth_fixed_cost':'First-birth fixed cost','h_P':'Parenthood housing floor','hbar_child_rooms':'Additional per-child floor','psi_child':'Child-utility intercept','payroll_tax':'Payroll tax','pension_period':'Period pension','housing_supply_elasticity':'Supply elasticity','tenure_choice_kappa':'Tenure shock scale','alpha_cons':'Consumption share','sigma':'CRRA curvature'}
rows=[['Parameter','Starting','Joint','Bounds / restriction','Near bound?']]
for r in joint['parameters']:
 n=r['parameter']; structural=r['structural_coordinate']
 restriction=f"[{fmt(r['lower'])}, {fmt(r['upper'])}]" if structural else ('Normalized' if n=='psi_child' else 'Budget derived' if n=='pension_period' else 'Fixed')
 rows.append([paramlabels[n],fmt(bp[n]['estimate']),fmt(r['estimate']),restriction,'Yes' if structural and r['near_bound'] else 'No' if structural else '-'])
 checked.extend([fmt(bp[n]['estimate']),fmt(r['estimate'])])
story += [table(rows,[160,80,80,112,60]),Spacer(1,10),para('Annual beta enters a four-year period as beta to the fourth power. Near-bound flags use the inherited rule: within 1% of the stated bound span. The local 12-by-9 Jacobian has numerical rank 9 but condition number about 1856; this is not evidence of strong identification.','small'),para('Source boundary: initial_sensitivity_panel and initial_joint_round_01/smoke_readout preserve the original solved packets. initial_calibration_contract/saved_case_scores contains complete machine-readable fits, parameters, source hashes and scoring receipts. Original model/observer diagnostic flags remain unchanged; these working scores do not certify a production benchmark.','small')]
def footer(c,doc):
 c.setFont('Helvetica',8);c.setFillColor(colors.gray);c.drawString(36,23,'Fertility project | verified saved results; transition completion outstanding');c.drawRightString(576,23,str(doc.page))
pdf=HERE/'morning_update.pdf'
SimpleDocTemplate(str(pdf),pagesize=(612,792),leftMargin=36,rightMargin=36,topMargin=32,bottomMargin=37,title='Morning quantitative update - September 11').build(story,onFirstPage=footer,onLaterPages=footer)
r=PdfReader(pdf);assert len(r.pages)==4,len(r.pages)
text='\n'.join(p.extract_text() for p in r.pages)
assert all(s in text for s in checked),[s for s in checked if s not in text]
(HERE/'morning_update.txt').write_text(text)
(HERE/'morning_update_qa.json').write_text(json.dumps({'pages':4,'numeric_cells_verified':len(checked),'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'source_contract':joint['contract_sha256'],'visual_review':'pending','cluster_status':'authentication unavailable; final root uncollected'},indent=2)+'\n')
print(pdf)
