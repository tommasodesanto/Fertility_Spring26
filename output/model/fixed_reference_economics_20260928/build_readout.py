"""Torch-only PDF assembly from authenticated existing economic-analysis outputs."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image

LABEL = '2007 stationary reference — block0506, September 28 verified export'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE = ROOT / 'output/model/fixed_reference_economics_20260928'
REF = ROOT / 'output/model/fertility_identification_20260928'


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def fmt(value):
    if value in ('', None):
        return '-'
    try:
        x = float(value)
    except (ValueError, TypeError):
        return str(value)
    return f'{x:.3e}' if x and (abs(x) < .001 or abs(x) >= 10000) else f'{x:.3f}'


def main():
    assert os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Render on Torch only'
    parser = argparse.ArgumentParser()
    parser.add_argument('--findings', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path = REF / 'fixed_reference_manifest.json'
    assert sha(manifest_path) == MANIFEST_SHA
    manifest = read(manifest_path)
    findings = read(args.findings)
    assert findings['reference_label'] == LABEL and findings['lead_reviewed'] is True
    export = REF / 'resume_v1/selected_export/primary'
    assert sha(export/'target_fit.csv') == manifest['artifact_hashes']['target_fit.csv']
    assert sha(export/'parameters.csv') == manifest['artifact_hashes']['parameters.csv']
    args.output.parent.mkdir(parents=True, exist_ok=True)
    assert not args.output.exists(), 'Do not overwrite an existing readout'
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='ReportBody', fontName='Helvetica', fontSize=9, leading=12, spaceAfter=8))
    styles.add(ParagraphStyle(name='SmallCell', fontName='Helvetica', fontSize=6.8, leading=8.5))
    styles.add(ParagraphStyle(name='Note', fontName='Helvetica', fontSize=8, leading=10, spaceAfter=6))
    story = []

    def para(text, style='ReportBody'):
        story.append(Paragraph(escape(text), styles[style]))

    def title(text):
        story.append(Paragraph(escape(text), styles['Heading1']))

    def grid(headers, data, widths):
        values = [[Paragraph(escape(str(x)), styles['SmallCell']) for x in row] for row in [headers]+data]
        t = Table(values, colWidths=widths, repeatRows=1, hAlign='LEFT')
        t.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4edf1')),
            ('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),5),
            ('TOPPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#537c8a')),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8f9')])]))
        story.append(t)

    def pic(path, caption, height=500):
        width_px, height_px = ImageReader(str(path)).getSize()
        scale = min(516/width_px, height/height_px)
        story.append(Image(str(path), width=width_px*scale, height=height_px*scale))
        para(caption, 'Note')

    title('Household behavior at the frozen calibration')
    para(LABEL)
    para('Common-primary reference loss: 19.581. All preferences, including the saved child-benefit parameter, remain fixed. The replacement-fertility normalization belongs to calibration and is not reapplied to shocks.')
    for paragraph in findings['summary']:
        para(paragraph)
    title('Interpretation and limits')
    for paragraph in findings['limitations']:
        para(paragraph)
    story.append(PageBreak())
    if findings.get('comparison_rows'):
        title('Prescribed housing prices: impact and cohort composition')
        para('Only the asset price and implied rent rise by 10%. All preferences, earnings, taxes, pensions, entry endowments, survival, credit rules and supply primitives stay fixed. The impact column uses the exact reference pre-choice distribution and household policies for permanently higher prices. The cohort column separately recomputes composition with normalized entry. Neither is a cleared equilibrium or a transition.')
        grid(['Object','Reference','Impact','Cohort'],findings['comparison_rows'],[255,87,87,87])
        para('Birth flows refer to four-year model periods. Completed fertility is undefined for the one-date impact column. Housing supply follows the unchanged supply curve when reporting excess demand; physical stock is not fixed in this diagnostic. The separately requested fixed-physical-stock transition is still being prepared. The full 14-row target comparison and 31-row parameter table follow the supplemental figures. Counterfactual fit against old data is a comparison, not a recalibration objective.', 'Note')
        story.append(PageBreak())
    for panel in findings.get('supplemental_panels', []):
        title('Supplemental: '+panel['title'])
        path = Path(panel['path'])
        assert sha(path) == panel['sha256']
        pic(path, panel['caption'], 540)
        story.append(PageBreak())
    title('Complete reference target fit')
    para('Normalization, scored moments and zero-weight validation rows retain their original definitions. CSVs retain full precision; gaps equal model minus target. These are 2007-reference comparisons.', 'Note')
    fits = rows(export/'target_fit.csv')
    assert len(fits) == 14
    grid(['Moment / role','Target','Model','Gap','Weight','Loss'],
        [[r['moment']+' / '+r['role']]+[fmt(r[k]) for k in ('target','model','gap','weight','loss_contribution')] for r in fits],
        [171,63,63,65,80,74])
    story.append(PageBreak())
    title('Complete reference parameters and restrictions')
    params = rows(export/'parameters.csv')
    assert len(params) == 31
    grid(['Parameter','Estimate','Lower','Upper','Near bound','Role / restriction'],
         [[r['parameter']]+[fmt(r[k]) for k in ('estimate','lower','upper')]+[r['near_bound'],r['status']] for r in params],
         [128,57,48,48,43,192])
    story.append(PageBreak())
    for case in findings.get('cases', []):
        p = Path(case['path'])
        assert sha(p/'receipt.json') == case['receipt_sha256']
        receipt = read(p/'receipt.json')
        assert receipt['status'] == 'passed' and receipt['reference_label'] == LABEL
        assert receipt['normalization_performed'] is False
        for relative, digest in receipt['artifact_hashes'].items():
            assert sha(p/relative) == digest, relative
        title(case['title'])
        para(case['description'])
        fits = rows(p/'target_fit.csv')
        assert len(fits) == 14
        grid(['Moment / role','Target','Model','Gap','Weight','Loss'],
            [[r['moment']+' / '+r['role']]+[fmt(r[k]) for k in ('target','model','gap','weight','loss_contribution')] for r in fits],
            [171,63,63,65,80,74])
        story.append(PageBreak())
        title(case['title']+': complete parameters')
        params = rows(p/'parameters.csv')
        assert len(params) == 31
        grid(['Parameter','Estimate','Lower','Upper','Near bound','Role / restriction'],
             [[r['parameter']]+[fmt(r[k]) for k in ('estimate','lower','upper')]+[r['near_bound'],r['status']] for r in params],
             [128,57,48,48,43,192])
        story.append(PageBreak())
    packets = [('Reference',export)]+[(c['title'],Path(c['path'])) for c in findings.get('cases',[])]
    for label, packet in packets:
        names = manifest['standard_diagnostic_names']
        assert len(names) == 17
        for start in range(0,17,2):
            title(label+': standard diagnostics')
            for name in names[start:start+2]:
                path = packet/'standard_diagnostics'/name
                if packet == export:
                    assert sha(path) == manifest['artifact_hashes']['standard_diagnostics/'+name]
                pic(path, name.replace('.png','').replace('_',' '), 267)
            if start+2 < 17 or packet != packets[-1][1]:
                story.append(PageBreak())

    def footer(canvas, doc):
        canvas.setFont('Helvetica',6.5)
        canvas.drawString(42,22,LABEL)
        canvas.drawRightString(570,22,str(doc.page))

    doc = SimpleDocTemplate(str(args.output),pagesize=(612,792),leftMargin=42,rightMargin=54,topMargin=35,bottomMargin=38,
        title='Fixed-calibration household behavior: block0506',author='Fertility project economic analysis')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    import fitz
    pdf = fitz.open(args.output)
    qa = args.output.parent / (args.output.stem+'_qa')
    qa.mkdir(exist_ok=False)
    texts = []
    for i,page in enumerate(pdf):
        texts.append(page.get_text())
        page.get_pixmap(matrix=fitz.Matrix(1.25,1.25)).save(qa/f'page_{i+1:02d}.png')
    joined = '\n'.join(texts)
    for r in rows(export/'parameters.csv'):
        assert r['parameter'] in joined, r['parameter']
    for r in rows(export/'target_fit.csv'):
        assert r['moment'] in joined, r['moment']
    (qa/'receipt.json').write_text(json.dumps(dict(pdf=str(args.output),sha256=sha(args.output),pages=len(pdf),
        reference_manifest_sha256=MANIFEST_SHA,findings_sha256=sha(args.findings),model_solves=0,
        target_rows=14,parameter_rows=31,standard_plot_count=17,visual_review='pending'),indent=2)+'\n')
    print(json.dumps(dict(pdf=str(args.output),pages=len(pdf),qa=str(qa))),flush=True)


if __name__ == '__main__':
    main()
