"""Render the final evening memo on Torch from authenticated saved tables only."""
import argparse, csv, json, math, shutil
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
def read(p):return json.loads(p.read_text())
def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def f(x):
 if x in ('',None):return '--'
 x=float(x)
 return f'{x:.2e}' if x and (abs(x)<.001 or abs(x)>=10000) else f'{x:.3f}'.rstrip('0').rstrip('.')
def main(root):
 root=Path(root);run=root/'output/model/evening_calibration_20260927/gated_v2/search'
 out=root/'output/model/evening_calibration_20260927/final_report';out.mkdir(parents=True,exist_ok=True)
 pdfpath=root/'output/pdf/evening_calibration_20260927.pdf';pdfpath.parent.mkdir(parents=True,exist_ok=True)
 final=read(run/'complete.json');auth=read(run.parent/'final_review/authentication.json')
 assert final['status']=='bounded_search_complete' and auth['status']=='authenticated_final_export'
 search=[r for r in final['records'] if r['design']!='repeat'];repeats=[r for r in final['records'] if r['design']=='repeat']
 counts=Counter(r['status'] for r in search)
 assert len(search)==360 and len(repeats)==6 and all(r['status']=='success' for r in repeats)
 best=min((r for r in search if r['status']=='success'),key=lambda r:r['primary_rescore'])
 assert best['case']==auth['common_primary_best']['case']==final['selected']['block']['case']
 fit=rows(run.parent/'final_review/target_fit_primary_rescore.csv')
 params=rows(run/'selected_export/block/parameters.csv');pmap={r['parameter']:r for r in params}
 assert len(fit)==14 and len(params)==31 and set(LABELS)=={r['moment'] for r in fit}
 assert len([r for r in params if r['lower']!=''])==10
 assert sum(r['role']=='scored' for r in fit)==10 and sum(r['role']=='validation' for r in fit)==3
 total=0.
 for r in fit:
  assert math.isclose(float(r['model'])-float(r['target']),float(r['gap']),rel_tol=1e-10,abs_tol=1e-12)
  if r['weight']!='':
   v=float(r['weight'])*float(r['gap'])**2;assert math.isclose(v,float(r['loss_contribution']),rel_tol=1e-10,abs_tol=1e-12);total+=v
 assert math.isclose(total,best['primary_rescore'],rel_tol=1e-12)
 for lane in ('primary','identity','block'):
  assert auth['lanes'][lane]['all14_31_byte_exact'] and auth['lanes'][lane]['all17_repeat_export_plot_hashes_exact']
  for filename in ('parameters.csv','target_fit.csv','export_receipt.json'):
   shutil.copy2(run/'selected_export'/lane/filename,out/(lane+'_'+filename))
 shutil.copy2(run.parent/'final_review/target_fit_primary_rescore.csv',out/'common_best_target_fit.csv')
 start=read(run.parent/'smoke/smoke_0000_primary/case/receipt.json')['loss']
 change=100*(1-total/start)
 summary={'status':'completed_bounded_exploration','search_counts':dict(counts),'search_attempts':len(search),
  'fresh_successful_repeats':len(repeats),'common_best':best,'start_primary_loss':start,
  'common_primary_loss_reduction_percent':change,'all_fit_rows':fit,'all_parameter_rows':params,
  'source_contract_sha256':final['contract_sha256'],'scientific_identity':final['scientific_identity'],
  'authentication_source':str(run.parent/'final_review/authentication.json'),'model_solves_for_report':0,
  'lane_winners':final['selected']}
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 C=canvas.Canvas(str(pdfpath),pagesize=(612,792));C.setTitle('Evening calibration - September 27, 2026')
 navy=colors.HexColor('#172b43');teal=colors.HexColor('#167884')
 body=ParagraphStyle('body',fontName='Helvetica',fontSize=9.1,leading=12.1,textColor=navy)
 small=ParagraphStyle('small',parent=body,fontSize=8.3,leading=10.7)
 cell=ParagraphStyle('cell',parent=body,fontSize=8.1,leading=10.3)
 def paragraph(text,y,style=body):
  p=Paragraph(text,style);_,h=p.wrap(536,720);p.drawOn(C,38,y-h);return y-h-8
 def title(n,sub):
  C.setFillColor(navy);C.setFont('Helvetica-Bold',17);C.drawString(38,752,'Fertility and housing: evening calibration')
  C.setFont('Helvetica',9);C.drawString(38,735,sub);C.setStrokeColor(teal);C.line(38,724,574,724)
  C.setFont('Helvetica',7);C.drawString(38,25,'September 27, 2026 | job 18672459 completed at 21:23 EDT | finite exploration')
  C.drawRightString(574,25,f'{n} / 2');return 710
 def heading(s,y):
  C.setFont('Helvetica-Bold',10);C.setFillColor(teal);C.drawString(38,y-11,s);return y-20
 def table(data,widths,y):
  t=Table([[Paragraph(escape(str(v)),cell) for v in row] for row in data],colWidths=widths)
  t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f0f3')),('VALIGN',(0,0),(-1,-1),'TOP'),
    ('TOPPADDING',(0,0),(-1,-1),3.4),('BOTTOMPADDING',(0,0),(-1,-1),3.4),('LEFTPADDING',(0,0),(-1,-1),5),
    ('RIGHTPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.6,teal),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')])]))
  _,h=t.wrap(536,720);t.drawOn(C,38,y-h);return y-h-9
 y=title(1,'What improved, what remains missing, and the complete target fit')
 y=paragraph(f'<b>The run completed, and fit improved:</b> the best common-primary score fell from {f(start)} to '
  f'<b>{f(total)}</b> ({f(change)}%). The point comes from the block-weight search. It is a reproducible candidate, not an established optimum.',y)
 y=paragraph('Scope: 24 single-thread Torch workers; 360 search attempts: 221 accepted, 77 housing-equilibrium gate rejections, '
  '61 owned timeouts and one late completion excluded. No fatal error; six fresh final repeats passed. '
  'The scientific contract stayed fixed during the accepted search; no local model computation.',y,small)
 compare=[['Search / reference','Score under the same primary weights'],['Starting point',f(start)]]
 compare.extend([[{'primary':'Primary weights','identity':'Relative-gap identity weights','block':'Equal-block weights (common best)'}[lane],f(final['selected'][lane]['primary_rescore'])] for lane in ('primary','identity','block')])
 y=table(compare,[340,196],y)
 y=heading('All 14 rows: ten scored moments, three validations, one normalization',y)
 data=[['Moment','Target','Model','Gap','Weight','Loss']]
 for r in fit:data.append([LABELS[r['moment']],f(r['target']),f(r['model']),f(r['gap']),f(r['weight']),f(r['loss_contribution'])])
 y=table(data,[220,56,56,56,73,75],y)
 y=paragraph('Gap = model - target. All weights and loss contributions above are primary weights, including zero weights for '
  'the three untargeted validation rows [V]. Completed fertility fixes the benefit level and is not scored. '
  'Identity and block objective scalars are never directly compared.',y,small)
 y=paragraph('<b>The main economic misses persist:</b> early fertility is 0.534 versus 0.810; ownership is 63.341% versus 67.626%; '
  'wealth/earnings is 6.230 versus 6.927. Housing rises too much at first birth (1.673 versus 1.465 rooms). '
  'Childlessness, one-child frequency and mean first-birth age are close to their targets.',y,small)
 assert y>43,('page1 overflow',y);C.showPage()
 y=title(2,'Parameter estimates, numerical confidence, and the next decisions')
 y=paragraph('<b>Eleven fitted quantities:</b> ten searched parameters, including the tenure-choice scale, plus the child-benefit '
  'level normalized to completed fertility 2.1. Ten scored moments plus that normalization give eleven restrictions by count; identification is not established.',y)
 data=[['Parameter','Estimate','Lower','Upper','Bound / restriction']]
 for name,label in NAMES.items():
  r=pmap[name];flag='near lower' if r['near_bound']=='True' and float(r['estimate'])<(float(r['lower'])+float(r['upper']))/2 else 'near upper' if r['near_bound']=='True' else 'interior'
  if name=='psi_child':flag='positive; normalized'
  data.append([label,f(r['estimate']),f(r['lower']),f(r['upper']),flag])
 y=table(data,[222,65,63,63,123],y)
 y=paragraph('Near-bound flags use 1% of the full search interval; neither fertility taste scale hits its lower bound. '
  'The remaining fixed/derived parameters are preserved in the full 31-row CSV. Annual real interest stays 2%; the existing-owner '
  'borrowing rule and separate death-solvency condition are the author-adopted DUE specification.',y,small)
 y=heading('What the numerical evidence supports',y)
 receipt=read(run/'selected_export/block/receipt.json')
 y=paragraph('Each of the three selected candidates has two fresh successful repeats with identical full 14-row target and '
  '31-row parameter tables. All 17 figures per lane match repeat sources and passed final visual review. The common best has relative housing-market '
  f'residual {f(receipt["market_residual"])} and completed-fertility error {f(receipt["normalization"]["absolute_gap"])}. '
  'The cluster job ended normally; peak memory was about 102.1 GiB of 128 GiB.',y,small)
 y=heading('Limits that matter for interpretation',y)
 y=paragraph('The first-birth housing contrast is a four-year model comparison rather than the exact empirical event-study estimator. '
  'Bequest data measure child-directed transfers while the model counts positive estates. The three zero-weight rows remain visible validation comparisons. '
  'The housing/fertility tradeoff is therefore a mechanism question and a measurement question.',y,small)
 y=paragraph('The standard policy plots are buyer-conditional: they do not certify every occupied stayer policy. High-wealth policy downturns, '
  'late ownership and retirement decumulation remain diagnostic caveats. Timeouts and gate rejections limit the explored region; a finite search '
  'does not establish an optimum or an identification rank.',y,small)
 y=heading('Proposed next decisions - no automatic continuation',y)
 for s in ('Keep the same targets for a longer, better-initialized search; diagnose the 900-second timeouts without relaxing equilibrium gates.',
           'Measure local sensitivity around this candidate, including tenure scale, housing loading and fertility preferences, to assess parameter substitution.',
           'Review the early-fertility / first-birth-housing tension and the provisional data counterparts; retain the three validation rows throughout.'):
  y=paragraph('&#8226; '+s,y,small)
 y=paragraph('Supporting evidence: final_report/summary.json, common_best_target_fit.csv and per-lane full parameter/fit CSVs; '
  'the stable 17 figures remain separate. No transition result or policy conclusion is claimed.',y,small)
 assert y>43,('page2 overflow',y);C.showPage();C.save()
 import fitz
 doc=fitz.open(pdfpath);assert len(doc)==2
 for i,page in enumerate(doc):page.get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(out/f'page_{i+1}.png')
 qa={'pages':len(doc),'target_rows':len(fit),'all_parameter_rows':len(params),'searched_parameters':10,'normalized_parameters':1,
     'pdf':str(pdfpath),'model_solves':0,'source_contract_sha256':final['contract_sha256'],'fit_loss_recomputed':total}
 (out/'report_qa.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',required=True);main(p.parse_args().project_root)
