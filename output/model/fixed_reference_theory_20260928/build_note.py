"""Torch-only, saved-table calculation and readable analytical PDF; zero solves."""
import os, sys, re, json, math, hashlib, csv
from pathlib import Path
from xml.sax.saxutils import escape
assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch allocation required'
import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties, findfont
from PIL import Image as PILImage
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
import fitz
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ECO=ROOT/'output/model/fixed_reference_economics_20260928/fixed_price_v1'
REF=ROOT/'output/model/fertility_identification_20260928'
SHA='147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(REF/'fixed_reference_manifest.json') == SHA
manifest=read(REF/'fixed_reference_manifest.json')
comp=read(ECO/'comparison.json'); complete=read(ECO/'completed.json')
inputs={}; receipts={}
for c in complete['completed']:
    name=c['case']; p=ECO/name/'receipt.json'
    assert sha(p)==c['receipt_sha256']; inputs[str(p.relative_to(ROOT))]=sha(p)
    r=read(p); receipts[name]=r
    assert r['status']=='passed' and r['normalization_performed'] is False
    assert r['reference_manifest_sha256']==SHA
    assert r['reference_checkpoint']['sha256']==manifest['checkpoint']['sha256']
    assert r['fixed_psi']==manifest['actual_serialized_parameters']['psi_child']
    for f in ['target_fit.csv','parameters.csv']:
        assert sha(ECO/name/f)==r['artifact_hashes'][f]
        with (ECO/name/f).open() as s: assert len(list(csv.DictReader(s))) == (14 if f=='target_fit.csv' else 31)
    assert r['standard_plot_count']==17
for case in ['control_1','control_2']:
    assert receipts[case]['completed_fertility']==receipts['control_1']['completed_fertility']
    c=receipts[case]['control']
    assert c['arrays']==113 and c['all_numeric_arrays_exact'] and c['all_14_fit_rows_exact'] and c['all_31_parameters_exact']
# Comparison JSON independently agrees with authenticated per-case summaries.
for phase,key in [('impact','baseline_state_impact_summary'),('cohort','cohort_summary')]:
    for item,vals in comp[phase].items():
        assert vals['reference']==receipts['control_1'][key][item]
        assert vals['price_110']==receipts['price_110'][key][item]
inputs[str((ECO/'comparison.json').relative_to(ROOT))]=sha(ECO/'comparison.json')
source_manifest=read(ROOT/'output/model/overnight_calibration_20260928/contract_v1/source_manifest.json')
assert sha(ROOT/'output/model/overnight_calibration_20260928/contract_v1/source_manifest.json')==manifest['source_manifest']['sha256']
source_pins={
'code/model/intergen_eqscale_seq_optimized/solver.py':'b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1',
'code/model/intergen_eqscale_seq_optimized/parameters.py':'66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464',
'code/model/intergen_eqscale_seq_optimized/child_preferences.py':'69544858b0137e09d37ac9e07ee6470ddc2d5849236b2c3e33d270737164f09d'}
for p,digest in source_pins.items(): assert sha(ROOT/p)==digest
for p in ['code/model/tools/run_e5f_perfect_foresight_transition.py','code/model/tools/run_e5f_open_population_transition.py']:
    source_pins[p]=sha(ROOT/p)
for p,digest in source_pins.items(): assert source_manifest['files'][p]==digest
P=manifest['actual_serialized_parameters']; qratio=receipts['price_110']['price'][0]/receipts['control_1']['price'][0]
assert abs(qratio-1.1)<1e-14
assert abs(receipts['price_110']['rent'][0]/receipts['control_1']['rent'][0]-qratio)<1e-14
metrics=[]
def add(label, phase, v0, v1):
    metrics.append(dict(outcome=label,distribution=phase,reference=v0,shock=v1,percent_change=100*(v1/v0-1),log_change_elasticity=math.log(v1/v0)/math.log(qratio),midpoint_arc_elasticity=((v1-v0)/((v1+v0)/2))/((qratio-1)/((qratio+1)/2))))
for phase,label in [('impact','Inherited states'),('cohort','Normalized cohort')]:
    for item,name in [('births_per_household','Raw birth flow'),('first_births','First-birth flow'),('second_births','Second-birth flow'),('third_bin_entries','Entry into 3+ flow'),('rooms_per_household','Rooms per household')]:
        v=comp[phase][item];add(name,label,v['reference'],v['price_110'])
add('Completed fertility','Normalized cohort',receipts['control_1']['completed_fertility'],receipts['price_110']['completed_fertility'])
# Adjustment for the unchanged top-bin child weight, explicitly distinct from event counts.
for phase,label in [('impact','Inherited states'),('cohort','Normalized cohort')]:
    vals=[comp[phase]['births']['reference']+(P['tfr_top_bin_weight']-3)*comp[phase]['third_bin_entries']['reference'],comp[phase]['births']['price_110']+(P['tfr_top_bin_weight']-3)*comp[phase]['third_bin_entries']['price_110']]
    add('Renewal-adjusted birth flow',label,*vals)
headers=['Outcome / distribution','Reference','Shock','Change %','Log elasticity']
table=['| '+' | '.join(headers)+' |','|---|---:|---:|---:|---:|']
# Compact seven-row table: detail for inherited states plus cohort birth/rooms/completed.
selected=[m for m in metrics if m['distribution']=='Inherited states' or m['outcome'] in ['Raw birth flow','Rooms per household','Completed fertility']]
for m in selected:
    table.append('| '+m['outcome']+' / '+('impact' if m['distribution']=='Inherited states' else 'cohort')+' | '+' | '.join(f'{m[k]:.3f}' for k in ['reference','shock','percent_change','log_change_elasticity'])+' |')
md=(HERE/'README.md').read_text()
md=re.sub(r'<!-- elasticities:start -->.*?<!-- elasticities:end -->','<!-- elasticities:start -->\n'+'\n'.join(table)+'\n<!-- elasticities:end -->',md,flags=re.S)
(HERE/'README.md').write_text(md)
qa=HERE/'qa';qa.mkdir(exist_ok=True); mathdir=qa/'math';mathdir.mkdir(exist_ok=True)
fontpath=findfont(FontProperties(family='DejaVu Sans'))
pdfmetrics.registerFont(TTFont('Body',fontpath))
pdfmetrics.registerFont(TTFont('BodyBold',findfont(FontProperties(family='DejaVu Sans',weight='bold'))))
styles={
'body':ParagraphStyle('body',fontName='Body',fontSize=9,leading=12,spaceAfter=7),
'small':ParagraphStyle('small',fontName='Body',fontSize=7.6,leading=10,spaceAfter=6),
'h1':ParagraphStyle('h1',fontName='BodyBold',fontSize=17,leading=21,spaceAfter=10),
'h2':ParagraphStyle('h2',fontName='BodyBold',fontSize=12,leading=15,spaceBefore=8,spaceAfter=7),
'h3':ParagraphStyle('h3',fontName='BodyBold',fontSize=10,leading=13,spaceBefore=7,spaceAfter=6),
'cell':ParagraphStyle('cell',fontName='Body',fontSize=7.3,leading=9.2),
}
cache={}
def mathimg(tex,display=False):
    key=(tex,display)
    if key not in cache:
        path=mathdir/f'eq_{len(cache):03d}.png'; size=13 if display else 9
        if display: tex=tex.replace('\\frac','\\dfrac')
        math_to_image('$'+tex+'$',str(path),prop=FontProperties(size=size),dpi=200,format='png',color='#122939')
        w,h=PILImage.open(path).size;w*=72/200;h*=72/200
        if display and w>495: h*=495/w;w=495
        cache[key]=(path,w,h)
    return cache[key]
def inline(text):
    # Protect equations before escaping prose.
    equations=[]
    def token(m):
        path,w,h=mathimg(m.group(1)); equations.append(f'<img src="{path}" width="{w:.2f}" height="{h:.2f}" valign="middle"/>');return f'ZZEQ{len(equations)-1}ZZ'
    text=re.sub(r'\\\((.*?)\\\)',token,text)
    links=[]
    def linktoken(m):
        target=m.group(2)
        if not target.startswith(('https://','http://')):
            relative=(HERE/target).resolve().relative_to(ROOT)
            target='file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/'+str(relative)
        links.append('<link href="'+escape(target,{'"':'&quot;'})+'" color="#245978">'+escape(m.group(1))+'</link>')
        return f'ZZLINK{len(links)-1}ZZ'
    text=re.sub(r'\[([^]]+)\]\(([^)]+)\)',linktoken,text)
    text=escape(text)
    text=re.sub(r'\*\*(.*?)\*\*',r'<font name="BodyBold">\1</font>',text)
    text=re.sub(r'`([^`]+)`',r'<font size="7">\1</font>',text)
    for i,x in enumerate(equations):text=text.replace(f'ZZEQ{i}ZZ',x)
    for i,x in enumerate(links):text=text.replace(f'ZZLINK{i}ZZ',x)
    return text
story=[];lines=md.splitlines();i=0;small=False
while i<len(lines):
    line=lines[i].strip();i+=1
    if line.startswith('<!-- handoff-index:'): break
    if not line:continue
    if line=='<!-- pagebreak -->':story.append(PageBreak());continue
    if line.startswith('<!--'):continue
    if line=='\\[':
        eq=[]
        while lines[i].strip()!='\\]':eq.append(lines[i].strip());i+=1
        i+=1;path,w,h=mathimg(' '.join(eq),True);story.append(Image(str(path),width=w,height=h,hAlign='CENTER'));story.append(Spacer(1,8));continue
    if line.startswith('|'):
        rows=[line]
        while i<len(lines) and lines[i].strip().startswith('|'):rows.append(lines[i].strip());i+=1
        values=[[Paragraph(inline(c.strip()),styles['cell']) for c in r.strip('|').split('|')] for r in rows if not r.startswith('|---')]
        t=Table(values,colWidths=[220,68,68,65,74],repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5edf3')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7f9')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.append(t);story.append(Spacer(1,8));continue
    if line.startswith('#'):
        n=len(line)-len(line.lstrip('#'));title=line[n:].strip()
        if title=='Evidence, limits and reproduction':small=True
        story.append(Paragraph(inline(title),styles['h'+str(min(n,3))]));continue
    para=[line]
    while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','<!--','\\[','|')):para.append(lines[i].strip());i+=1
    story.append(Paragraph(inline(' '.join(para)),styles['small' if small else 'body']))
def footer(c,doc):
    c.setFont('Body',7);c.setFillColor(colors.HexColor('#526674'));c.drawString(50,28,'Frozen block0506 | Analytical identities and saved finite changes');c.drawRightString(562,28,str(doc.page))
pdf=HERE/'theory_note.pdf'
SimpleDocTemplate(str(pdf),pagesize=(612,792),leftMargin=50,rightMargin=50,topMargin=38,bottomMargin=44,title='Fertility, housing costs and supply',author='Research note').build(story,onFirstPage=footer,onLaterPages=footer)
doc=fitz.open(pdf)
for k,page in enumerate(doc):page.get_pixmap(matrix=fitz.Matrix(1.25,1.25)).save(qa/f'page_{k+1}.png')
receipt=dict(reference_label=manifest['label'],slurm_job=os.environ['SLURM_JOB_ID'],model_solves=0,checkpoint_loaded=False,reference_manifest_sha256=SHA,checkpoint_sha256=manifest['checkpoint']['sha256'],inputs=inputs,verified_source_hashes=source_pins,source_locations={'attempts':'solver.py:3486-3550','continuation':'solver.py:3348-3453','forward_timing':'solver.py:5533','dated_rent':'run_e5f_perfect_foresight_transition.py:358-386','adjusted_births':'run_e5f_open_population_transition.py:767'},asset_price_ratio=qratio,metrics=metrics,first_birth_share_impact_decline=comp['impact']['first_births']['difference']/comp['impact']['births']['difference'],pdf_pages=len(doc),pdf_sha256=sha(pdf),readme_sha256=sha(HERE/'README.md'),status='calculation_and_render_pass_visual_review_pending')
(HERE/'calculation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(job=receipt['slurm_job'],pages=len(doc),metrics=selected,pdf_sha256=receipt['pdf_sha256']),indent=2))
