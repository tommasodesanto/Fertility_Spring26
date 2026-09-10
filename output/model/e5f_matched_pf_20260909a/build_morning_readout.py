"""Build the dated partial advisor PDF from independently verified receipts."""
from pathlib import Path
import csv,json,hashlib
from xml.sax.saxutils import escape
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from review_horizon100_root_progress import LABELS,main as verify
verify()
base=Path(__file__).resolve().parent
root=base/'meeting_receipts/historical_root_h100_01/sequential'
review=json.loads((base/'horizon100_root_progress_review.json').read_text())
trial=root/f"evaluation_{review['best_collected_evaluation']:03d}"
summary=json.loads((trial/'summary.json').read_text())
assert review['best_collected_evaluation']==5
with (trial/'target_fit.csv').open() as f: fits=list(csv.DictReader(f))
with (trial/'parameters.csv').open() as f: parameters=list(csv.DictReader(f))
out=base.parents[1]/'pdf/matched_pf_morning_readout_20260910.pdf'
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='BodySmall',fontName='Helvetica',fontSize=10,leading=14,spaceAfter=9))
styles.add(ParagraphStyle(name='CellSmall',fontName='Helvetica',fontSize=8,leading=10))
styles.add(ParagraphStyle(name='NumberCell',fontName='Helvetica',fontSize=7,leading=10))
styles.add(ParagraphStyle(name='NoteSmall',fontName='Helvetica',fontSize=8.5,leading=11,spaceAfter=8))
styles['Title'].fontSize=21;styles['Title'].leading=25;styles['Title'].alignment=TA_LEFT
styles['Heading2'].fontSize=13;styles['Heading2'].spaceBefore=12
story=[]
def p(s,style='BodySmall'):return Paragraph(s,styles[style])
def add(s,style='BodySmall'):story.append(p(s,style))
def table(rows,widths):
 t=Table([[p(escape(str(x)).replace('–','-').replace('—','-'),'CellSmall' if i==0 else 'NumberCell') for i,x in enumerate(row)] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5edf3')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#7890a0')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f9')])]))
 story.append(t)
def num(x):return '-' if x=='' else f'{float(x):.7g}'
add('Quantitative work: morning readout','Title')
add('10 September 2026 | Partial results | Status checked at 08:23 EDT','NoteSmall')
add('<b>We have made substantial progress in solving the historical price path, but we do not yet have a completed new calibration or a matched policy result.</b> The 100-date sequential path now reproduces exactly at its best prices. Its largest housing-market error is 0.04527%, still above our unchanged 0.02% requirement.')
add('What was solved overnight','Heading2')
add('At the eleven inherited parameter values, the solver adjusted all 100 dated prices jointly. For each candidate price path, households anticipate the entire future from 2007; their decisions are solved backward and the population is propagated forward. A verified 100-column price-response matrix guides the price updates. This matrix concerns market clearing; it is not evidence that the economic parameters are identified.')
table([['Path evaluation','Largest market error','Meaning']]+[[str(t['evaluation']),f"{100*t['maximum_market_residual']:.6f}%",'Fresh final replay' if t['evaluation']==6 else ('Initial prices' if t['evaluation']==1 else 'Updated prices')] for t in review['verified_trials']],[95,130,282])
story.append(Spacer(1,10))
add('All six evaluations passed their mapping, accounting and feasibility checks. The final replay has zero difference in prices and residuals; the target-fit, parameter, measurement and transition tables are byte-identical to the selected evaluation. This checks reproducibility, not convergence.')
add('The remaining work, in order','Heading2')
add('<b>1. Finish market clearing.</b> Continuation job 17319026 is running on Torch. At 08:23 it was in its first path, date 44 of 100, with a fresh heartbeat. The bounded round allows three paths including a final replay; there is no automatic further round.')
add('<b>2. Establish horizon stability.</b> The terminal unit-rent gap is 1.0867%, above its 1% gate. Other recorded endpoint checks pass. Historical prices and moments still need comparison across sufficiently long horizons; extending the endpoint alone is not a certificate.')
add('<b>3. Reconcile measurement and recalibrate.</b> Four pooled ACS targets are currently compared with a 2023 model cross-section, and two family-group definitions remain unresolved. Keep the reviewed childbirth regression unchanged. Then assess parameter sensitivities and conduct matched re-estimation before producing policies.')
story.append(PageBreak())
add('Fit at the reproduced provisional prices','Title')
add(f'Objective: {summary["loss"]:.8f}. This is an inherited-parameter diagnostic, not a calibrated-equilibrium loss. Every target and weight is unchanged. Shares and their gaps are in fraction units.','NoteSmall')
table([['Moment','Target','Model','Gap','Weight','Loss']]+[[LABELS[r['moment']]]+[num(r[k]) for k in ('target','model','gap','weight','loss_contribution')] for r in fits],[182,65,65,65,65,65])
story.append(Spacer(1,10))
add('What the fit tells us','Heading2')
add('Completed fertility is 1.777 versus 1.918, and childlessness is 23.36% versus 18.80%. These two rows contribute about 64.02 of the 94.48 objective. The first-birth housing response is 0.411 rooms against 0.720, while average occupied housing is 6.326 rooms against 5.780. Better market clearing alone has not repaired these economic discrepancies.')
add('Which years are being targeted','Heading2')
add('The initial economy approximates a pre-2007 steady state, normalized to completed fertility 2.1. Households learn the path immediately in 2007. The model is measured along that transition: 2023 is an observation date, not a steady state. Cohort and wealth moments retain their approved source-specific dates; four ACS rows pool 2012-2023, and the birth-response branch uses 2019-2023. This objective does not fit an annual historical time series.')
add('The historical age composition is imposed through 2023; the person-demographic model advances the population afterward. This remains a conditional historical fit. The numerical target fingerprint alone cannot resolve sample, family-group or timing differences.','NoteSmall')
add('Policy status','Heading2')
add('There is no new matched perfect-foresight property-tax comparison to report. The intended main comparison raises the annual tax from 1% to 2%, with equal household rebates in both regimes, starting from the same certified historical 2023 state. The unrebated historical baseline and the rebated policy baseline are distinct objects. Older policy effects must not be presented as results of this unfinished calibration.')
story.append(PageBreak())
add('Parameters and remaining checks','Title')
add('All eleven free coordinates are inherited inputs, not new estimates. Bounds and near-bound flags below reproduce the saved parameter contract. A target count of twelve against eleven coordinates does not itself establish identification.','NoteSmall')
rows=[['Parameter','Value','Lower','Upper','Role','Near bound']]
for r in parameters:
 role='Free input' if r['is_free_parameter']=='True' else ('Normalized' if r['parameter']=='psi_child_2007' else 'Derived' if r['parameter']=='psi_child_2023' else 'Fixed')
 rows.append([r['parameter'],num(r['value']),num(r['lower_bound']),num(r['upper_bound']),role,'Yes' if r['near_bound']=='True' else 'No'])
table(rows,[180,75,55,55,85,57])
add('Old fertility preferences are normalized to the initial completed-fertility target. The 2023 preference level equals that old intercept plus the transition coordinate. Tenure taste dispersion (0.005) and housing-supply elasticity (0.63) remain externally fixed. The near-bound flag concerns theta1; it is not a claim of parameter identification.','NoteSmall')
add('All endpoint-distance checks','Heading2')
tail=summary['terminal_distance']
table([['Distance','Value','Tolerance','Pass']]+[[k.replace('_',' '),num(v),num(tail['tolerances'][k]),'Yes' if tail['checks'][k] else 'No'] for k,v in tail['metrics'].items()],[285,80,80,62])
story.append(Spacer(1,8))
add('Evidence: HORIZON100_PROGRESS.md, horizon100_final_reproduction_review.json and overnight_target_mapping_review.md under output/model/e5f_matched_pf_20260909a/. The numerical snapshot is 96a41873; the original and continuation results remain separate. No source, target, weight or numerical gate changed for the continuation.','NoteSmall')
def footer(canvas,doc):
 canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#526675'));canvas.drawString(44,24,'Fertility project | Verified partial readout | 10 September 2026');canvas.drawRightString(A4[0]-44,24,str(doc.page))
SimpleDocTemplate(str(out),pagesize=A4,title='Quantitative work: morning readout',author='Fertility project',rightMargin=44,leftMargin=44,topMargin=38,bottomMargin=38).build(story,onFirstPage=footer,onLaterPages=footer)
print(out)
