#!/usr/bin/env python3
"""Build the September 22 earnings-and-wealth diagnostic review PDF.

This is a reporting-only builder. It reads the frozen packet and does not run
the model or alter diagnostic images.
"""
from __future__ import annotations

import csv, json, os, re, hashlib
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph,
    Spacer, Table, TableStyle, PageBreak, Image, KeepTogether, HRFlowable)

ROOT = Path(__file__).resolve().parents[3]
PACK = ROOT / 'output/model/native_financing_diagnostic_20260919/specification_followup/earnings_wealth_v1'
OUT = ROOT / 'output/pdf/earnings_wealth_review.pdf'
VERIFY = ROOT / 'output/pdf/earnings_wealth_review_verification.json'
TMP = ROOT / 'tmp/pdfs/earnings_wealth_review'
PLOTS = PACK / 'search_readout/standard_diagnostics'

def read_csv(p):
    with open(p, newline='') as f: return list(csv.DictReader(f))
def read_json(p):
    with open(p) as f: return json.load(f)
def ascii_text(x):
    s = str(x or '')
    repl = {'\u2013':'-', '\u2014':'-', '\u2212':'-', '\u2011':'-', '\u00d7':'x', '\u2264':'<=', '\u2265':'>=', '\u03b2':'beta', '\u03c1':'rho', '\u03c3':'sigma', '\u03b8':'theta', '\u03ba':'kappa', '\u03c7':'chi', '\u2018':"'", '\u2019':"'", '\u201c':'"', '\u201d':'"'}
    for a,b in repl.items(): s=s.replace(a,b)
    return s.encode('ascii','ignore').decode('ascii')
def fmt(x, n=6):
    try:
        v=float(x)
        if abs(v) >= 1000: return f'{v:,.2f}'
        return f'{v:.{n}g}'
    except Exception: return ascii_text(x)
def p(txt, style='BodyText'):
    return Paragraph(ascii_text(txt), STY[style])

styles = getSampleStyleSheet()
STY = {
 'Title': ParagraphStyle('Title', parent=styles['Title'], fontName='Helvetica-Bold', fontSize=20, leading=24, textColor=colors.HexColor('#17324d'), spaceAfter=12),
 'H1': ParagraphStyle('H1', parent=styles['Heading1'], fontName='Helvetica-Bold', fontSize=14, leading=17, textColor=colors.HexColor('#17324d'), spaceBefore=5, spaceAfter=7),
 'H2': ParagraphStyle('H2', parent=styles['Heading2'], fontName='Helvetica-Bold', fontSize=10.5, leading=13, textColor=colors.HexColor('#2c526e'), spaceBefore=5, spaceAfter=4),
 'BodyText': ParagraphStyle('BodyText', parent=styles['BodyText'], fontName='Helvetica', fontSize=10.5, leading=14, textColor=colors.HexColor('#232323'), spaceAfter=6),
 'Small': ParagraphStyle('Small', parent=styles['BodyText'], fontName='Helvetica', fontSize=8.5, leading=11, textColor=colors.HexColor('#333333')),
 'Tiny': ParagraphStyle('Tiny', parent=styles['BodyText'], fontName='Helvetica', fontSize=8.5, leading=10.5, textColor=colors.HexColor('#333333')),
 'Caption': ParagraphStyle('Caption', parent=styles['BodyText'], fontName='Helvetica-Oblique', fontSize=8.5, leading=11, textColor=colors.HexColor('#555555'), spaceBefore=2, spaceAfter=5),
 'TableHead': ParagraphStyle('TableHead', parent=styles['BodyText'], fontName='Helvetica-Bold', fontSize=8.5, leading=10.5, textColor=colors.white),
 'TableCell': ParagraphStyle('TableCell', parent=styles['BodyText'], fontName='Helvetica', fontSize=8.5, leading=10.5),
 'TableCellR': ParagraphStyle('TableCellR', parent=styles['BodyText'], fontName='Helvetica', fontSize=8.5, leading=10.5, alignment=TA_RIGHT),
}

def table(data, widths, header=True, font='TableCell'):
    cells=[]
    for i,row in enumerate(data):
        rr=[]
        for val in row:
            sty='TableHead' if header and i==0 else (font if not isinstance(val, Paragraph) else None)
            rr.append(val if isinstance(val, Paragraph) else Paragraph(ascii_text(val), STY[sty]))
        cells.append(rr)
    t=Table(cells, colWidths=widths, repeatRows=1 if header else 0, hAlign='LEFT')
    t.setStyle(TableStyle([
      ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#2f607c')) if header else ('BACKGROUND',(0,0),(-1,-1),colors.white),
      ('TEXTCOLOR',(0,0),(-1,0),colors.white), ('GRID',(0,0),(-1,-1),0.25,colors.HexColor('#b8c5ce')),
      ('VALIGN',(0,0),(-1,-1),'TOP'), ('LEFTPADDING',(0,0),(-1,-1),4), ('RIGHTPADDING',(0,0),(-1,-1),4),
      ('TOPPADDING',(0,0),(-1,-1),3), ('BOTTOMPADDING',(0,0),(-1,-1),3),
      ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white, colors.HexColor('#f1f5f7')]),
    ]))
    return t

TARGET_LABELS={'initial_normalization':'Fertility normalization (unscored)','cps_childlessness':'Childlessness','cps_exactly_one':'Exactly one child','nchs_mean_age':'Mean first-birth age','nchs_share30':'First births at age 30+','wealth_earnings':'Wealth / annual labor earnings','bequest_wealth':'Annual bequests / wealth','old_dispersion':'Old-age wealth dispersion','mean_rooms':'Mean occupied rooms','ownership_30_55':'Ownership ages 30-55','first_birth_rooms':'First-birth rooms response','family_rooms':'Family rooms response','recent_parent_ownership':'Recent-parent ownership difference'}
def target_rows(path):
    rows=read_csv(path);out=[]
    for r in rows:
        name=TARGET_LABELS[r['restriction_id']]
        numeric=[fmt(r.get(k),7) for k in ['target','model','gap']]
        weight=fmt(r.get('actual_weight'),6) if r.get('actual_weight') else '-'
        loss=fmt(r.get('loss_contribution'),6) if r.get('loss_contribution') else '-'
        out.append([name,*numeric,weight,loss])
    assert len(out)==13
    return out

def param_rows(files):
    names=[]; allrows={}
    for label,path in files.items():
        for r in read_csv(path):
            name=r['parameter']; names.append(name) if name not in names else None
            allrows.setdefault(name,{})[label]=r
    return names,allrows

def parameter_table(files, names, allrows, subset):
    bounds=read_json(PACK/'staging/frozen_v5_plan.json')['parameter_bounds']
    out=[['Parameter','Anchor','Beta 0.985','Search 005','Actual bounds / restriction','Near: A / B / S']]
    labels={'beta_annual':'Annual discount factor','kappa_fert':'First-birth shock scale','kappa_fert_continuation':'Later-birth shock scale','chi':'Owner service multiplier','H0':'Supply intercept','theta0':'Bequest scale','theta1':'Bequest wealth shift','first_birth_fixed_cost':'First-birth utility cost','h_P':'Parenthood housing floor','hbar_child_rooms':'Extra floor per child','psi_child':'Fertility utility intercept','payroll_tax':'Payroll tax rate','pension_period':'Four-year pension','housing_supply_elasticity':'Supply elasticity','tenure_choice_kappa':'Tenure shock scale','alpha_cons':'Consumption exponent','sigma':'Utility curvature'}
    for name in names[subset[0]:subset[1]]:
        vals=[fmt(allrows[name][label]['estimate'],7) for label in files]
        if name in bounds:
            lo,hi=bounds[name]
            restriction=f'[{lo:g}, {hi:g}]'
            flags=['Y' if min(float(allrows[name][label]['estimate'])-lo,hi-float(allrows[name][label]['estimate'])) <= .01*(hi-lo) else 'N' for label in files]
            note=' / '.join(flags)
        else:
            restriction=allrows[name]['Anchor']['status'];note='-'
        out.append([labels[name],*vals,restriction,note])
    return out

class Doc(BaseDocTemplate):
    def __init__(self, fn):
        BaseDocTemplate.__init__(self,fn,pagesize=letter,leftMargin=.55*inch,rightMargin=.55*inch,topMargin=.58*inch,bottomMargin=.48*inch,title='Earnings and wealth: overnight review',author='Research review for Tommaso De Santo',subject='September 22, 2026: evidence, incomplete calibration and next decisions')
        frame=Frame(self.leftMargin,self.bottomMargin,self.width,self.height,id='normal')
        self.addPageTemplates([PageTemplate(id='main',frames=[frame],onPage=self.footer)])
    def footer(self, canv, doc):
        canv.saveState(); canv.setStrokeColor(colors.HexColor('#c8d2d8')); canv.setLineWidth(.4); canv.line(.55*inch,.38*inch,7.95*inch,.38*inch)
        canv.setFont('Helvetica',7); canv.setFillColor(colors.HexColor('#52616b')); canv.drawString(.55*inch,.23*inch,'Earnings and wealth diagnostic review - September 22, 2026 - diagnostic only')
        canv.drawRightString(7.95*inch,.23*inch,f'{doc.page}'); canv.restoreState()

def build():
    OUT.parent.mkdir(parents=True,exist_ok=True); TMP.mkdir(parents=True,exist_ok=True)
    files={'Anchor':PACK/'smoke_v5/anchor/parameters.csv','beta=.985 smoke':PACK/'smoke_v5/probe/parameters.csv','Search case005':PACK/'search_readout/parameters/parameters_scored_case005.csv'}
    target_files={'Anchor':PACK/'smoke_v5/anchor/target_fit.csv','beta=.985 smoke':PACK/'smoke_v5/probe/target_fit.csv','Search case005':PACK/'search_readout/target_fit/target_fit_case005.csv'}
    scores={k:read_json(v) for k,v in {'Anchor':PACK/'smoke_v5/anchor/score.json','beta=.985 smoke':PACK/'smoke_v5/probe/score.json','Search case005':PACK/'search_readout/receipts/score.json'}.items()}
    terminal=read_json(PACK/'staging/v5_terminal_lead_review.json')
    names,allrows=param_rows(files)
    comparison=[['Moment','Target','Anchor','Beta 0.985','Search 005']]
    by_point={label:{r['restriction_id']:r for r in read_csv(path)} for label,path in target_files.items()}
    for key,label,scale in [('nchs_mean_age','First-birth age',1),('cps_childlessness','Childlessness (%)',100),('mean_rooms','Mean rooms',1),('ownership_30_55','Ownership 30-55 (%)',100)]:
        comparison.append([label,f"{float(by_point['Anchor'][key]['target'])*scale:.2f}",*[f"{float(by_point[point][key]['model'])*scale:.2f}" for point in target_files]])
    story=[]
    story += [p('Earnings and wealth under four-year decisions','Title'),p('Operational diagnostic review of the native financing candidate'),p('September 22, 2026 | source packet frozen after the 07:50 EDT search stop','Small'),Spacer(1,10),HRFlowable(width='100%',color=colors.HexColor('#2f607c'),thickness=1),Spacer(1,10)]
    story += [p('Executive reading','H1'),p('The persistent-plus-iid earnings construction is operationally viable at the tested points: the full smoke passed, the direct four-year process uses 45 income states, and the expanded 160-node wealth grid reaches [-12, 3000] with zero assets at age 18. The evidence remains diagnostic. The longer search scored 7 of 8 proposals, stopped on a native timeout during proposal 006, and never completed final repeats or selection verification.'),p('The candidate has four-year persistence 0.77614467, persistent innovation SD 0.43743545, and iid SD 0.16499061. These are direct four-year estimates from the saved data construction. Bootstrap evidence is weak for the iid component: its 199-draw bootstrap interval for the transitory variance runs from approximately zero to 0.06398. The point fit of three covariance moments to three parameters is mechanical and is not a proof of identification.'),p('Across the three complete points, the beta=.985 smoke has the lowest loss (1040.447). The anchor is 1160.376 and search case005 is 1087.282. Case005 improves first-birth timing but worsens childlessness, occupied rooms, and ownership. It is a useful trade-off diagnostic, not an adopted calibration.'),p('The fit trade-off at a glance','H2'),table(comparison,[2.25*inch,1.1*inch,1.1*inch,1.1*inch,1.55*inch]),Spacer(1,8),p('The earnings candidate works at the tested points. A defensible final calibration is still unfinished. Full target and parameter tables follow; no specification has been adopted.','H2'),PageBreak()]
    story += [p('What is established, and what is not','H1'),p('Each of the three reported points completed a normalized solution requiring six stationary solves. The starting point reproduced exactly across two runs; the beta 0.985 smoke point and search proposal 005 were each evaluated once. The saved contract keeps 13 target rows, with 12 scored moments plus a separately imposed completed-fertility normalization. Nine structural coordinates are free; the remaining rows are fixed, derived, or restricted.'),p('The long search had 49 stationary solves started, 48 completed, and 1 incomplete. Including V5 smoke evidence, the honest count is 67 started, 66 completed, and 1 incomplete; the all-attempt count is 102 started, 98 completed, and 4 incomplete. The controller recorded a fatal contract error, but the observed cause is the native 3,100-second timeout; there is no evidence of a changed target or source fingerprint.'),p('The income process follows the architecture precedent in Sommer (2016 JME): an AR(1) persistent component plus an iid component without a fixed permanent type. De Nardi (2004 ReStud) and Bick (2016 JEEA) support direct estimation at multi-year model frequencies, but their five-year and three-year objects are not drop-in four-year parameters. The saved literature review therefore supports the architecture and timing choice while leaving external parameter comparisons qualified.'),p('What the earnings estimates measure','H2'),p('The PSID input measures gross reference-person plus spouse labor earnings. We sum four consecutive observed annual values, then remove block-start age and year effects from log block means. The selected annual-era sample (1984-1997 survey-year labels, ages 25-60) has 7,812 blocks for 3,932 persons, with 3,836 lag-one and 1,449 lag-two pairs. Complete positive-earnings histories are a selected sample; this is not a post-1997 estimate.'),p('There is no fixed permanent type. The persistent component follows an AR(1), and the second component is iid. Four-year income is known when choices are made. Current income enters purchase eligibility, with accounting checked to avoid spending it twice. Neither direct estimation nor that accounting check resolves within-period timing.'),p('Zero assets at entry age 18 and stationary initial earnings risk are external candidate assumptions: the estimation sample has no ages 18-24. The wealth grid is a numerical domain, not an economic saving cap. Its extension is not convergence evidence. Original figures retain crowded/clipped 45-state legends and wealth-axis compression; these inherited presentation limits are disclosed.'),PageBreak()]
    # Three target tables, one page each.
    for label,path in target_files.items():
        loss=scores[label]['loss']
        story += [p(f'Target fit: {label}', 'H1'), p(f'Loss = {loss:.6f}. Each row is read from the saved target-fit CSV; values are rounded for display directly from the saved records. The fertility normalization is separate from the scored objective.'), table([['Target moment','Target','Model','Gap','Actual weight','Loss contribution']]+target_rows(path),[2.2*inch,.9*inch,.9*inch,.9*inch,1.1*inch,1.1*inch]), Spacer(1,8), p('Interpretation: the table is descriptive evidence under an unchanged target contract. A lower scalar loss at one point does not establish a local optimum, identification, or adoption.', 'Caption'), PageBreak()]
    # parameters, two pages
    story += [p('All 17 parameters: estimates, bounds, and restrictions','H1'),p('The parameter comparison below distinguishes actual bounds from fixed or derived objects. Near-bound flags are recalculated from the actual frozen plan: within 1% of each raw bound width. A / B / S denote the anchor, beta 0.985 test and search proposal. Annual beta is bounded above by 0.99, and the parenthood floor by 2.3. These flags do not establish identification. At search case005, beta_annual is at its actual upper bound and theta0 is at zero, so theta1 is inactive in the implemented bequest utility.'), table(parameter_table(files,names,allrows,(0,9)),[1.65*inch,.85*inch,.85*inch,.85*inch,1.55*inch,.85*inch]), PageBreak()]
    story += [p('All 17 parameters continued','H1'),table(parameter_table(files,names,allrows,(9,18)),[1.65*inch,.85*inch,.85*inch,.85*inch,1.55*inch,.85*inch]),Spacer(1,7),p('Fixed and derived rows are included to make the full contract visible. The reported theta1 value is not identified when theta0=0; this is a structural inactivity statement, not a claim that the optimizer found a unique value.', 'Caption'),Spacer(1,10),p('Units and economic interpretation','H2'),p('At search case005, annual beta 0.99 implies a four-year discount factor of 0.96059601. The parenthood housing floor is 1.9146 rooms when children are resident; the additional floor per child is zero. H0 is the housing supply intercept, and chi multiplies owner housing services.'),p('The bequest scale theta0 is zero at this point, so the bequest wealth shift theta1 has no effect on utility here. Positive accidental bequests may still occur. The fertility and tenure kappa parameters are taste-shock scales: their numerical sizes depend on the utility normalization and have no directly portable literature benchmark.'),PageBreak()]
    story += [p('Numerical limits and next steps','H1'),p('The search result should be treated as a partial operational review. Seven valid proposals are insufficient to characterize a frontier, and the best valid proposal was evaluated once. The completed smoke point at beta=.985 is outside the frozen search candidate pool, so comparing its loss to case005 is informative but not a selection rule.'),p('Earlier attempts exposed three implementation problems: transactions clipped outside grid support, loss of a sale option at the old wealth-grid ceiling, and reported consumption/housing floors that differed from the optimizer allocation. Separately reviewed corrections passed targeted checks; older evidence is preserved. No economic target or numerical gate was relaxed.'),p('The principal numerical limits are the unverified continuation/interpolation behavior at tiny-mass states, wealth-grid density and domain sensitivity, and the absence of exact repetitions for case005. The all-17 graph packet remains the stable diagnostic set. No evidence supports dropping or demoting any target, and no unreachable-target conclusion is warranted.'),p('Recommended next sequence: preserve the direct four-year construction; repair and stub-test the timeout/failure path before another long run; declare the selection pool and any changed duration in a new frozen plan; resolve the continuation and grid diagnostics; then repeat a declared candidate exactly before a further focused calibration under the same targets. Inspect the full target and parameter tables before paper adoption.'),p('These are empirical estimates, with material uncertainty. The 199-draw person bootstrap gives persistence 0.7211-0.8348, persistent variance 0.4302-0.5249, and iid variance approximately 0-0.0640 (95% percentile intervals). The iid component is weakly estimated. Precision, the entrant distribution and within-period timing still need to be defended before paper adoption.', 'BodyText'),PageBreak()]
    story += [p('References and provenance','H1'),p('Sommer, Kamila (2016). “Fertility Choice in a Life Cycle Model with Idiosyncratic Uninsurable Earnings Risk.” Journal of Monetary Economics 83, 27-38. Author PDF: https://www.kamilasommer.net/Fertility.pdf'),p('De Nardi, Mariacristina (2004). “Wealth Inequality and Intergenerational Links.” Review of Economic Studies 71(3), 743-768. Primary PDF: https://users.nber.org/~denardim/research/denardi.pdf'),p('Bick, Alexander (2016). “The Quantitative Role of Child Care for Female Labor Force Participation and Fertility.” Journal of the European Economic Association 14(3), 639-668. Article record: https://onlinelibrary.wiley.com/doi/abs/10.1111/jeea.12143'),p('Saved project reviews: literature/annual_sources_review.md; literature/multiyear_review.md; literature/parameter_validation_rubric.md; terminal_review.md. These sources distinguish architecture precedent from directly estimated four-year parameters.'),p('Report provenance','H2'),p(f'Generated from {PACK}. Source data are read at build time. Graphs are copied by reference into the PDF from search_readout/standard_diagnostics and are not edited or regenerated.'),p('Frozen run fingerprints (SHA-256)','H2'),p('Objective: '+terminal['objective_canonical_sha256'],'Small'),p('Plan: '+terminal['plan_sha256'],'Small'),p('Economic source manifest: '+terminal['source_fingerprints']['economic_source_manifest_7e872053'],'Small'),p('Exact anchor repetition checks cover prices, value functions, population arrays and moments. No such completed final comparison exists for the improved smoke point or search case005. Full generated-code and target provenance are retained in the source packet.','Small'),PageBreak()]
    # 17 original plots, 2 per page (nine pages)
    plot_files=sorted(PLOTS.glob('*.png'))
    for i in range(0,len(plot_files),2):
        story += [p(f'Appendix: original standard diagnostic plots ({i+1}-{min(i+2,len(plot_files))} of {len(plot_files)})','H1')]
        for q in plot_files[i:i+2]:
            im=Image(str(q)); im._restrictSize(7.1*inch,3.85*inch); story += [p(q.stem.replace('_',' '),'H2'),im,p('Original case005 graph retained without redesign. Read wealth-axis compression and crowded/clipped 45-state legends with care.','Caption')]
        story.append(PageBreak())
    doc=Doc(str(OUT)); doc.build(story)
    # exact numeric sidecar
    from pypdf import PdfReader
    payload={'pdf':str(OUT),'source_packet':str(PACK),'page_count':len(PdfReader(str(OUT)).pages),'target_rows':{k:len(target_rows(v)) for k,v in target_files.items()},'parameter_rows':len(names),'plot_count':len(plot_files),'losses':{k:scores[k]['loss'] for k in scores},'source_sha256':{}}
    payload['actual_parameter_bounds']=read_json(PACK/'staging/frozen_v5_plan.json')['parameter_bounds']
    payload['pdf_sha256']=hashlib.sha256(OUT.read_bytes()).hexdigest()
    payload['visual_review']='pending lead page inspection'
    for q in list(target_files.values())+list(files.values())+[PACK/'staging/frozen_v5_plan.json']+plot_files: payload['source_sha256'][str(q.relative_to(ROOT))]=hashlib.sha256(q.read_bytes()).hexdigest()
    VERIFY.write_text(json.dumps(payload,indent=2)+'\n')
    return payload

if __name__=='__main__':
    print(json.dumps(build(),indent=2))
