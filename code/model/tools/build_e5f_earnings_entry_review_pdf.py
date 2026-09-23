#!/usr/bin/env python3
"""Build a readout-driven earnings-entry review PDF.

The collector owns the readout. This builder only formats its files and fails
with an explicit missing-cell report when the collection is incomplete.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, Table, TableStyle, PageBreak, Image, HRFlowable

def read_json(p):
    with open(p) as f: return json.load(f)
def read_csv(p):
    with open(p, newline='') as f: return list(csv.DictReader(f))
def ascii_text(x):
    repl={'\u2013':'-','\u2014':'-','\u2212':'-','\u2011':'-','\u2264':'<=','\u2265':'>=','\u03b2':'beta','\u03c1':'rho','\u03c3':'sigma','\u03b8':'theta','\u03ba':'kappa','\u03c7':'chi','\u2018':"'",'\u2019':"'",'\u201c':'"','\u201d':'"'}
    s=str(x or '')
    for a,b in repl.items(): s=s.replace(a,b)
    return s.encode('ascii','ignore').decode('ascii')
def fmt(x, n=5):
    try:
        v=float(x)
        return f'{v:,.2f}' if abs(v)>=1000 else f'{v:.{n}g}'
    except Exception: return ascii_text(x)

styles=getSampleStyleSheet()
STY={
 'Title':ParagraphStyle('Title',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=19,leading=23,textColor=colors.HexColor('#17324d'),spaceAfter=11),
 'H1':ParagraphStyle('H1',parent=styles['Heading1'],fontName='Helvetica-Bold',fontSize=14,leading=17,textColor=colors.HexColor('#17324d'),spaceBefore=5,spaceAfter=7),
 'H2':ParagraphStyle('H2',parent=styles['Heading2'],fontName='Helvetica-Bold',fontSize=10.5,leading=13,textColor=colors.HexColor('#2c526e'),spaceBefore=5,spaceAfter=4),
 'Body':ParagraphStyle('Body',parent=styles['BodyText'],fontName='Helvetica',fontSize=10.5,leading=14,textColor=colors.HexColor('#222222'),spaceAfter=6),
 'Small':ParagraphStyle('Small',parent=styles['BodyText'],fontName='Helvetica',fontSize=8.5,leading=10.5,textColor=colors.HexColor('#333333')),
 'Head':ParagraphStyle('Head',parent=styles['BodyText'],fontName='Helvetica-Bold',fontSize=8.5,leading=10,textColor=colors.white),
 'Cell':ParagraphStyle('Cell',parent=styles['BodyText'],fontName='Helvetica',fontSize=8.5,leading=10),
 'Cap':ParagraphStyle('Cap',parent=styles['BodyText'],fontName='Helvetica-Oblique',fontSize=8.5,leading=10.5,textColor=colors.HexColor('#555555'),spaceBefore=2,spaceAfter=5),
}
def para(x, sty='Body'): return Paragraph(escape(ascii_text(x)),STY[sty])
def make_table(rows,widths):
    cells=[]
    for i,row in enumerate(rows):
        cells.append([v if isinstance(v,Paragraph) else Paragraph(escape(ascii_text(v)),STY['Head' if i==0 else 'Cell']) for v in row])
    t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#2f607c')),('GRID',(0,0),(-1,-1),.25,colors.HexColor('#b8c5ce')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f1f5f7')])]))
    return t

class Doc(BaseDocTemplate):
    def __init__(self,fn,title):
        super().__init__(fn,pagesize=letter,leftMargin=.55*inch,rightMargin=.55*inch,topMargin=.58*inch,bottomMargin=.48*inch,title=title)
        self.addPageTemplates([PageTemplate(id='main',frames=[Frame(self.leftMargin,self.bottomMargin,self.width,self.height,id='normal')],onPage=self.footer)])
    def footer(self,canv,doc):
        canv.saveState(); canv.setStrokeColor(colors.HexColor('#c8d2d8')); canv.line(.55*inch,.38*inch,7.95*inch,.38*inch); canv.setFont('Helvetica',7); canv.setFillColor(colors.HexColor('#52616b')); canv.drawString(.55*inch,.23*inch,'Earnings entry review - diagnostic readout'); canv.drawRightString(7.95*inch,.23*inch,str(doc.page)); canv.restoreState()

def cell_path(readout,arm,cell): return Path(readout)/arm/cell
def moment_label(row):
    return 'Initial fertility normalization' if row.get('restriction_id')=='initial_normalization' else row.get('label',row.get('restriction_id',''))

def target_rows(path):
    rows=read_csv(path)
    if len(rows)!=13: raise ValueError(f'{path}: expected 13 target rows, found {len(rows)}')
    out=[['Target moment','Target','Model','Gap','Actual weight','Loss contribution']]
    for r in rows:
        out.append([moment_label(r),fmt(r.get('target')),fmt(r.get('model')),fmt(r.get('gap')),fmt(r.get('actual_weight')) if r.get('actual_weight') else '-',fmt(r.get('loss_contribution')) if r.get('loss_contribution') else '-'])
    return out
def parameter_rows(path):
    rows=read_csv(path)
    if len(rows)!=17: raise ValueError(f'{path}: expected 17 parameter rows, found {len(rows)}')
    out=[['Parameter','Estimate','Actual lower','Actual upper','Near bound','Restriction/status']]
    for r in rows:
        lo=r.get('actual_lower',r.get('lower_bound',r.get('lower',''))); hi=r.get('actual_upper',r.get('upper_bound',r.get('upper','')))
        status='search coordinate' if str(r.get('structural_coordinate','')).lower()=='true' else r.get('restriction_label',r.get('status',''))
        out.append([r.get('parameter',''),fmt(r.get('estimate')),fmt(lo) if lo else '-',fmt(hi) if hi else '-',ascii_text(r.get('near_actual_bound',r.get('near_bound_1pct_range',r.get('near_bound','-')))),ascii_text(status)])
    return out
def narrative(path):
    d=read_json(path)
    for k in ('title','date','summary','recommendations','limitations'):
        if k not in d: raise ValueError(f'narrative missing key {k}')
    return d

def selected_comparison_story(readout):
    grouped={}
    for arm in 'ABCD':
        path=Path(readout)/arm/'selected'/'target_fit.csv'
        if not path.exists(): continue
        for row in read_csv(path):
            item=grouped.setdefault(row['restriction_id'],{'label':moment_label(row),'target':row['target']})
            assert float(item['target'])==float(row['target'])
            item[arm]=row['model']
    rows=[['Moment','Target','A','B','C','D']]+[[r['label'],fmt(r['target'])]+[fmt(r[a]) if a in r else '-' for a in 'ABCD'] for r in grouped.values()]
    return [PageBreak(),para('Selected points: targets and model moments','H1'),para('A/B use one persistent earnings shock; C/D add iid earnings risk. A/C enter with zero assets; B/D use the inherited wealth marginal. Each column is its best verified observed point, including smoke. Search coverage can differ; this is not a comparison of converged calibrations.'),make_table(rows,[2.1*inch,.85*inch,.9*inch,.9*inch,.9*inch,.9*inch]),Spacer(1,8),para('Full gaps, weights, loss contributions, all17 parameter restrictions and original figures follow in the appendix. The next page holds the nine search parameters fixed across the four specifications.','Small')]

def resolution_story(readout, resolution_readout):
    fine=Path(resolution_readout)
    coarse=Path(readout)
    by_cell={}
    status=[]
    for arm in "ABCD":
        cp=coarse/arm/"smoke"
        fp=fine/arm/"smoke"
        if not (fp/"score.json").exists():
            status.append(f"{arm}: no verified finer-grid score")
            continue
        c=read_json(cp/"score.json"); f=read_json(fp/"score.json")
        assert c["contract_sha256"]==f["contract_sha256"], "Resolution comparison mixes objectives"
        assert read_json(cp/"plan.json")["structural_parameters"]==read_json(fp/"plan.json")["structural_parameters"], "Resolution comparison mixes parameters"
        cr={r["restriction_id"]:r for r in c["target_fit"]}
        fr={r["restriction_id"]:r for r in f["target_fit"]}
        assert cr.keys()==fr.keys()
        for key in cr:
            assert cr[key]["target"]==fr[key]["target"] and cr[key]["actual_weight"]==fr[key]["actual_weight"]
        by_cell[arm]={k:fr[k]["model"]-cr[k]["model"] for k in cr}
        status.append(f"{arm}: scored at finer resolution")
    labels={}
    for row in read_csv(coarse/"common_smoke_targets.csv"):
        labels.setdefault(row["restriction_id"],moment_label(row))
    rows=[["Moment","A change","B change","C change","D change"]]
    rows += [[label]+[fmt(by_cell[a][k]) if a in by_cell else "unavailable" for a in "ABCD"] for k,label in labels.items()]
    return [para("Income-grid sensitivity at common parameters","H1"),para("Entries are finer-grid minus coarse-grid model moments, in each moment's original units. A/B compare 15 versus 7 income states; C/D compare 45 versus 21. Preferences, targets and the wealth grid are held fixed; equilibrium and the existing fertility normalization are resolved."),make_table(rows,[2.55*inch]+[1.1*inch]*4),Spacer(1,9),para("; ".join(status)),para("Unavailable does not mean zero change. Failure or unfinished status is described in the brief review. These are single-point sensitivity checks, not grid-convergence certificates. Full finer-grid tables and original diagnostic figures are retained in the accompanying resolution readout.","Small"),PageBreak()]

def build(readout,narrative_path,output,resolution_readout=None):
    readout=Path(readout); output=Path(output); narr=narrative(narrative_path); output.parent.mkdir(parents=True,exist_ok=True)
    summary=read_json(readout/'collection_summary.json') if (readout/'collection_summary.json').exists() else {}
    arms=['A','B','C','D']; cells=['selected','smoke']; available=[]; missing=[]
    for arm in arms:
        for cell in cells:
            base=cell_path(readout,arm,cell)
            if (base/'score.json').exists() and (base/'target_fit.csv').exists() and (base/'parameters_actual_bounds.csv').exists(): available.append((arm,cell,base))
            else: missing.append(f'{arm}/{cell}')
    story=[para(narr['title'],'Title'),para(narr['date'],'Small'),HRFlowable(width='100%',color=colors.HexColor('#2f607c'),thickness=1),Spacer(1,10),para('Brief summary','H1')]
    for x in narr['summary']: story.append(para(x))
    story += [para('Collection status','H2'),para(f"Available cells: {', '.join(a+'/'+c for a,c,_ in available) or 'none'}. Explicitly unavailable or missing: {', '.join(missing) or 'none'}. No numeric value is inferred for an unavailable cell.")]
    if summary:
        status_rows=[['Cell','Scored','Failed/unscored','Running/incomplete','Unstarted']]
        for arm in arms:
            counts=summary.get('cells',{}).get(arm,{})
            status_rows.append([arm,str(counts.get('verified_scored',0)),str(counts.get('failed_unscored',0)+counts.get('collection_rejected',0)),str(counts.get('running_or_incomplete',0)),str(counts.get('not_started',0))])
        story.append(make_table(status_rows,[.6*inch,1.0*inch,1.6*inch,1.8*inch,1.3*inch]))
    story.append(para('Recommendations','H2'))
    for x in narr['recommendations']: story.append(para(x))
    story.append(para('Limitations','H2'))
    for x in narr['limitations']: story.append(para(x))
    story += selected_comparison_story(readout)
    if (readout/'common_smoke_targets.csv').exists():
        grouped={}
        for row in read_csv(readout/'common_smoke_targets.csv'):
            key=row['restriction_id']
            item=grouped.setdefault(key,{'label':moment_label(row),'target':row['target']})
            if float(item['target'])!=float(row['target']): raise ValueError('Mixed targets in common smoke comparison')
            item[row['cell_id']]=row['model']
        common=[['Moment','Target','A','B','C','D']]+[[r['label'],fmt(r['target'])]+[fmt(r[a]) if a in r else '-' for a in arms] for r in grouped.values()]
        story += [PageBreak(),para('Common parameters: smoke comparison','H1'),para('Each available smoke holds the nine search parameters fixed. Equilibrium and the fertility utility scale are resolved separately to match the same normalization target in each specification. Missing cells have no verified scored result.'),make_table(common,[2.1*inch,.85*inch,.9*inch,.9*inch,.9*inch,.9*inch]),PageBreak()]
    else: story.append(PageBreak())
    if resolution_readout is not None:
        story += resolution_story(readout,resolution_readout)
    for arm,cell,base in available:
        if cell!='selected': continue
        score=read_json(base/'score.json')
        story += [para(f'{arm} / {cell}: complete numeric readout','H1'),para(f"Loss: {fmt(score.get('loss'),9)}. Best verified observed point in this cell, including smoke; not a converged calibration. Exact values remain in the accompanying CSV files."),make_table(target_rows(base/'target_fit.csv'),[2.1*inch,.82*inch,.82*inch,.82*inch,1.05*inch,1.05*inch]),Spacer(1,8),para('All 17 parameters','H2'),para('Near bound means within 1% of the raw bound width; wide log-search bounds can flag small positive scales. This is not evidence of identification or convergence.','Small'),make_table(parameter_rows(base/'parameters_actual_bounds.csv'),[1.55*inch,.9*inch,.9*inch,.9*inch,.85*inch,1.25*inch]),PageBreak()]
        plots=sorted((base/'standard_diagnostics').glob('*.png'))
        if len(plots)!=17: raise ValueError(f'{base}: expected 17 standard diagnostics, found {len(plots)}')
        for i in range(0,17,2):
            story.append(para(f'Appendix: {arm} / {cell} original standard diagnostics ({i+1}-{min(i+2,17)} of 17)','H1'))
            for q in plots[i:i+2]:
                story += [para(q.stem.replace('_',' '),'H2')]
                im=Image(str(q)); im._restrictSize(7.1*inch,3.65*inch)
                caption='Original source legend is crowded or clipped; image retained unchanged.' if arm in 'CD' and q.stem in {'housing_by_age_income_state','ownership_by_age_income_state','wealth_dist_childless_renter_age30','wealth_dist_childless_renter_age42'} else 'Original diagnostic image retained without substitution or redesign.'
                story += [im,para(caption,'Cap')]
            if i + 2 < 17:
                story.append(PageBreak())
        story.append(PageBreak())
    if isinstance(story[-1],PageBreak): story.pop()
    Doc(str(output),ascii_text(narr['title'])).build(story)
    return {'output':str(output),'pages':len(__import__('pypdf').PdfReader(str(output)).pages),'available_cells':[f'{a}/{c}' for a,c,_ in available],'explicit_missing_cells':missing,'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--readout',required=True); ap.add_argument('--narrative',required=True); ap.add_argument('--output',required=True); ap.add_argument('--resolution-readout'); args=ap.parse_args(); print(json.dumps(build(args.readout,args.narrative,args.output,args.resolution_readout),indent=2))
if __name__=='__main__': main()
