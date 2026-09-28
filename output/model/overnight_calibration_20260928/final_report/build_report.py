"""Render the overnight memo on Torch from authenticated saved tables only."""
import argparse, csv, json, math, shutil, hashlib, os
from collections import Counter
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle

LABELS = {'initial_normalization':'Completed fertility (normalization)',
 'cps_childlessness':'Childlessness, ages 40-44','cps_exactly_one':'Exactly one child, among mothers',
 'nchs_mean_age':'Mean age at first birth','nchs_share30':'First births at age 30+ [V]',
 'wealth_earnings':'Aggregate wealth / annual earnings','bequest_wealth':'Annual bequest flow / wealth',
 'old_dispersion':'Older wealth/income: p90 / median [V]','mean_rooms':'Mean occupied rooms',
 'ownership_30_55':'Ownership, ages 30-55','first_birth_rooms':'Housing response to first birth',
 'family_rooms':'Rooms: 3+ vs 1-2 children at home [V]',
 'recent_parent_ownership':'Recent-parent ownership gap','early_fertility':'Children ever born at 25, capped at 3'}
NAMES={'H0':'Housing supply level','beta_annual':'Annual discount factor','chi':'Owner housing-service premium',
 'first_birth_fixed_cost':'First-birth fixed cost','kappa_fert':'First-birth taste scale',
 'kappa_fert_continuation':'Later-birth taste scale','theta0':'Bequest strength',
 'delta_alpha_jump':'First-child housing loading','child_benefit_curvature':'Child-benefit curvature',
 'tenure_choice_kappa':'Tenure-choice taste scale','psi_child':'Child-benefit level (normalized)'}
def read(p):return json.loads(Path(p).read_text())
def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def f(x):
 if x in ('',None):return '--'
 x=float(x)
 return f'{x:.2e}' if x and (abs(x)<.001 or abs(x)>=10000) else f'{x:.3f}'.rstrip('0').rstrip('.')
def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(ok, message):
 if not ok: raise ValueError(message)
def main(a):
 require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Render on Torch in a Slurm allocation only')
 root=Path(a.project_root);run=Path(a.run) if a.run else root/'output/model/overnight_calibration_20260928/gated_v1/search'
 out=root/'output/model/overnight_calibration_20260928/final_report';out.mkdir(parents=True,exist_ok=True)
 final=read(run/'complete.json');auth=read(a.authentication)
 require(final['status']=='bounded_search_complete','No completed search: use an honest failure memo, not this candidate renderer')
 require(auth['status']=='authenticated_final_export','Lead authentication required')
 require(sha(run/'complete.json')==auth['complete_sha256'],'Changed final completion receipt')
 require(auth['contract_sha256']==final['contract_sha256'],'Contract mismatch')
 selected=final['selected'];key=final['common_primary_key'];best=selected[key]
 require(set(selected)>={'primary','identity','block'},'Missing weight lane')
 require(auth['common_primary_best']['case']==best['case'],'Common best mismatch')
 for lane, original in selected.items():
  review=auth['lanes'][lane]
  require(review['all14_31_byte_exact'] and review['all17_repeat_export_plot_hashes_exact'] and review['plots_visually_reviewed'],'Incomplete final review: '+lane)
  packet=run/'selected_export'/lane;export=read(packet/'export_receipt.json')
  require(export['contract_sha256']==final['contract_sha256'] and export['selected']['case']==original['case'],'Export mismatch')
  require(len(export['repeats'])==2 and all(r['status']=='success' and r['point']==original['point'] for r in export['repeats']),'Two successful fresh repeats required')
  required=[packet/n for n in ('target_fit.csv','parameters.csv','receipt.json','export_receipt.json')]+list((packet/'standard_diagnostics').rglob('*.png'))
  require(len(required)==21,'Require precisely 17 standard PNGs')
  for path in required:
   rel=str(path.relative_to(run));require(auth['selected_files_sha256'].get(rel)==sha(path),'Missing/stale authenticated file: '+rel)
 fitpath=Path(a.fit);require(sha(fitpath)==auth['primary_fit_sha256'],'Primary fit mismatch')
 fit=rows(fitpath);params=rows(run/'selected_export'/key/'parameters.csv');pmap={r['parameter']:r for r in params}
 require(len(fit)==14 and set(LABELS)=={r['moment'] for r in fit} and len(params)==31,'Full 14/31 tables required')
 require(len([r for r in params if r['lower']!=''])==10 and set(NAMES)<=set(pmap),'Fitted parameter contract mismatch')
 require(sum(r['role']=='scored' for r in fit)==10 and sum(r['role']=='validation' for r in fit)==3,'Target roles mismatch')
 rawfit={r['moment']:r for r in rows(run/'selected_export'/key/'target_fit.csv')}
 total=0.
 for r in fit:
  require(all(r[k]==rawfit[r['moment']][k] for k in ('target','model','gap','role')),'Rescored table changed physical fit')
  require(math.isclose(float(r['model'])-float(r['target']),float(r['gap']),rel_tol=1e-10,abs_tol=1e-12),'Gap mismatch')
  if r['weight']!='':
   v=float(r['weight'])*float(r['gap'])**2;require(math.isclose(v,float(r['loss_contribution']),rel_tol=1e-10,abs_tol=1e-12),'Contribution mismatch');total+=v
 for r in params:
  if r['lower']!='':require(float(r['lower'])<=float(r['estimate'])<=float(r['upper']),'Parameter outside bounds')
 require(math.isclose(total,best['primary_rescore'],rel_tol=1e-11),'Common score mismatch')
 # Records contain search and repeat cases, but hourly render records live separately.
 search=[r for r in final['records'] if not r['case'].startswith('repeat_')]
 repeats=[r for r in final['records'] if r['case'].startswith('repeat_')]
 require(len(repeats)==2*len(selected) and all(r['status']=='success' for r in repeats),'Final repeat ledger mismatch')
 counts=Counter(r['status'] for r in search)
 findings=read(a.findings);require(findings['reviewed_by_lead'] is True,'Lead-reviewed interpretation required')
 for name in ('price_diagnostic','early_diagnostic'):
  d=findings[name];require(d['status'] in ('reviewed','pending','failed','unverified'),'Diagnostic status missing')
  if d['status']=='reviewed':
   require(d.get('evidence_sha256') and sha(d['evidence_path'])==d['evidence_sha256'],'Diagnostic evidence missing/stale')
 start=26.368192413904882;change=100*(1-total/start)
 for lane in selected:
  for filename in ('parameters.csv','target_fit.csv','export_receipt.json'):
   shutil.copy2(run/'selected_export'/lane/filename,out/(lane+'_'+filename))
 shutil.copy2(fitpath,out/'common_best_target_fit.csv')
 summary=dict(status='authenticated_bounded_exploration',common_best=best,common_primary_key=key,search_counts=dict(counts),search_attempts=len(search),fresh_successful_repeats=len(repeats),start_primary_loss=start,primary_loss=total,loss_reduction_percent=change,all_fit_rows=fit,all_parameter_rows=params,lane_winners=selected,findings=findings,source_contract_sha256=final['contract_sha256'],authentication_sha256=sha(a.authentication),model_solves_for_report=0)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 pdfpath=Path(a.pdf) if a.pdf else root/'output/pdf/overnight_calibration_20260928.pdf';pdfpath.parent.mkdir(parents=True,exist_ok=True)
 C=canvas.Canvas(str(pdfpath),pagesize=(612,792));C.setTitle('Overnight fertility calibration - September 28, 2026')
 navy=colors.HexColor('#172b43');teal=colors.HexColor('#167884')
 body=ParagraphStyle('body',fontName='Helvetica',fontSize=8.8,leading=11.4,textColor=navy)
 small=ParagraphStyle('small',parent=body,fontSize=8,leading=10.1)
 cell=ParagraphStyle('cell',parent=body,fontSize=7.7,leading=9.6)
 def paragraph(text,y,style=body):
  p=Paragraph(text,style);_,h=p.wrap(536,720);p.drawOn(C,38,y-h);return y-h-6
 def title(n,sub):
  C.setFillColor(navy);C.setFont('Helvetica-Bold',17);C.drawString(38,752,'Fertility and housing: overnight calibration')
  C.setFont('Helvetica',9);C.drawString(38,735,sub);C.setStrokeColor(teal);C.line(38,724,574,724)
  C.setFont('Helvetica',7);C.drawString(38,25,'September 28, 2026 | Torch job 18687184 | finite exploration, unchanged model')
  C.drawRightString(574,25,f'{n} / 2');return 710
 def heading(s,y):
  C.setFont('Helvetica-Bold',10);C.setFillColor(teal);C.drawString(38,y-11,s);return y-19
 def table(data,widths,y):
  t=Table([[Paragraph(escape(str(v)),cell) for v in row] for row in data],colWidths=widths)
  t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f0f3')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),2.8),('BOTTOMPADDING',(0,0),(-1,-1),2.8),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),('LINEBELOW',(0,0),(-1,0),.6,teal),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')])]))
  _,h=t.wrap(536,720);t.drawOn(C,38,y-h);return y-h-7
 y=title(1,'What improved and the complete target fit')
 y=paragraph(f'<b>Fit improved by {f(change)}%:</b> the same primary-weight score fell from {f(start)} to <b>{f(total)}</b>. The selected point ({escape(best["case"])}) passes the original gates and two fresh repeats; this does not establish an optimum.',y)
 counttext=', '.join(f'{n} {s.replace("_"," ")}' for s,n in sorted(counts.items()))
 y=paragraph(f'{len(search)} search attempts: {escape(counttext)}. {len(repeats)} final repeats cover {len(selected)} selected candidates. Up to 24 single-thread Torch workers; no local model computation. Model, targets, bounds and gates unchanged.',y,small)
 compare=[['Search / reference','Loss under the same primary weights'],['Starting point',f(start)]]
 compare.extend([[{'primary':'Primary weights','identity':'Relative-gap identity weights','block':'Equal-block weights','common_primary':'Additional common-primary winner'}[lane],f(row['primary_rescore'])] for lane,row in selected.items()])
 y=table(compare,[340,196],y)
 y=heading('All 14 rows: ten scored moments, three checks, one normalization',y)
 data=[['Moment','Target','Model','Gap','Weight','Loss']]
 for r in fit:data.append([LABELS[r['moment']],f(r['target']),f(r['model']),f(r['gap']),f(r['weight']),f(r['loss_contribution'])])
 y=table(data,[213,52,52,55,82,82],y)
 y=paragraph('Gap = model - target. Weights and loss contributions are the common primary system. [V] denotes an untargeted validation row. Completed fertility fixes the child-benefit level and is not scored. Raw identity and block objective values are not comparable.',y,small)
 y=paragraph('<b>Economic fit:</b> '+escape(findings['economic_fit']),y,small)
 require(y>43,f'Page 1 overflow at {y}: shorten findings, do not drop rows');C.showPage()
 y=title(2,'Fitted quantities, failure diagnosis and early-fertility tradeoffs')
 y=paragraph('Ten searched parameters plus one normalized child-benefit level. Ten scored moments plus the normalization give eleven restrictions by count; local identification rank is unproven.',y)
 data=[['Parameter','Estimate','Lower','Upper','Bound / restriction']]
 for name,label in NAMES.items():
  r=pmap[name];flag='normalized; positive' if name=='psi_child' else 'interior'
  if name!='psi_child' and r['near_bound']=='True':flag='near lower' if float(r['estimate'])<(float(r['lower'])+float(r['upper']))/2 else 'near upper'
  data.append([label,f(r['estimate']),f(r['lower']),f(r['upper']),flag])
 y=table(data,[222,65,63,63,123],y)
 y=paragraph('Near-bound = within 1% of the full interval. All 31 parameters remain in the supporting CSV. Fixed inputs include 2% annual real interest, B15 earnings and the adopted existing-owner credit rule. No experimental credit transition is adopted or solved here.',y,small)
 y=heading('Numerics and the earlier failures',y)
 y=paragraph(escape(findings['numerical_review']),y,small)
 y=paragraph('<b>Price-start diagnostic ('+escape(findings['price_diagnostic']['status'])+'):</b> '+escape(findings['price_diagnostic']['text']),y,small)
 y=heading('How much earlier fertility can these parameters buy?',y)
 y=paragraph('<b>Exploratory diagnostic ('+escape(findings['early_diagnostic']['status'])+'):</b> '+escape(findings['early_diagnostic']['text']),y,small)
 y=paragraph('Compare the early-fertility gap with primary-weighted loss from all other moments, excluding its own contribution. The small timing-preference probe is not a calibrated SMM estimate, global frontier, or proof of an unreachable target. Boundary proposals are diagnostics, not adoption.',y,small)
 y=heading('Next steps and limits',y)
 y=paragraph(escape(findings['next_steps']),y,small)
 y=paragraph('Buyer-conditional standard plots do not certify every occupied stayer policy. High-wealth housing/ownership downturns and retirement behavior remain caveats. First-birth housing uses a four-year model contrast; estate and empirical child-transfer measures also differ. Full tables and the stable 17 plots are separate supporting artifacts.',y,small)
 require(y>43,f'Page 2 overflow at {y}: shorten findings, do not drop parameters');C.showPage();C.save()
 import fitz
 doc=fitz.open(pdfpath);require(len(doc)==2,'Expected two pages')
 for i,page in enumerate(doc):page.get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(out/f'page_{i+1}.png')
 qa=dict(pages=len(doc),target_rows=14,parameter_rows=31,fitted_quantities=11,pdf=str(pdfpath),model_solves=0,fit_loss_recomputed=total,authentication_sha256=sha(a.authentication))
 (out/'report_qa.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',required=True);p.add_argument('--run');p.add_argument('--authentication',required=True);p.add_argument('--fit',required=True);p.add_argument('--findings',required=True);p.add_argument('--pdf');main(p.parse_args())
