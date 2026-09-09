#!/usr/bin/env python3
"""Render a supplied, verified tax-report manifest without importing model code.

The caller owns all results and economic interpretation. The builder requires
all 12 target rows and all 11 free parameters, preserves complete images, and
writes a machine-readable QA sidecar. Rendering inspection remains required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer,
    Table, TableStyle,
)

SCHEMA = 'e5f_rebated_tax_morning_report_v1'
BLUE = colors.HexColor('#204b6b')
PALE = colors.HexColor('#edf3f7')
GREY = colors.HexColor('#59636d')
PAGE_W, PAGE_H = A4
MARGIN = 40
WIDTH = PAGE_W - 2*MARGIN


def digest(path):
    hasher = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def display(value):
    """Numbers retain six significant digits; strings retain supplied precision."""
    if value is None:
        return 'Not available'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    if isinstance(value, (int, float)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite report number; supply an explicit missing-value label')
        return format(value, '.6g')
    return str(value)


def register_fonts():
    candidates = [
        ('/System/Library/Fonts/Supplemental/Arial.ttf',
         '/System/Library/Fonts/Supplemental/Arial Bold.ttf'),
        ('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
         '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'),
        ('/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf',
         '/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf'),
    ]
    for regular, bold in candidates:
        if Path(regular).is_file() and Path(bold).is_file():
            pdfmetrics.registerFont(TTFont('ReportBody', regular))
            pdfmetrics.registerFont(TTFont('ReportBold', bold))
            return 'ReportBody', 'ReportBold', [regular, bold]
    raise RuntimeError('No reviewed Unicode font pair found; install Arial, DejaVu Sans or Liberation Sans')


class Builder:
    def __init__(self, manifest, manifest_path):
        self.manifest = manifest
        self.base = manifest_path.parent
        self.inputs = {str(manifest_path): digest(manifest_path)}
        self.tables = []
        self.figures = []
        self.pages = 0
        normal, bold, font_paths = register_fonts()
        self.normal, self.bold = normal, bold
        self.font_paths = font_paths
        self.styles = {
            'title': ParagraphStyle('title', fontName=bold, fontSize=21, leading=25,
                                    textColor=BLUE, spaceAfter=12),
            'heading': ParagraphStyle('heading', fontName=bold, fontSize=15, leading=19,
                                      textColor=BLUE, spaceAfter=9),
            'body': ParagraphStyle('body', fontName=normal, fontSize=10, leading=14,
                                   spaceAfter=8),
            'small': ParagraphStyle('small', fontName=normal, fontSize=8, leading=10.5,
                                    textColor=GREY, spaceAfter=5),
            'cell': ParagraphStyle('cell', fontName=normal, fontSize=8, leading=10),
            'headcell': ParagraphStyle('headcell', fontName=bold, fontSize=8, leading=10,
                                       textColor=BLUE),
            'caption': ParagraphStyle('caption', fontName=normal, fontSize=8.5, leading=11,
                                      spaceAfter=7),
        }

    def paragraph(self, text, style='body'):
        return Paragraph(escape(display(text)).replace('\n', '<br/>'), self.styles[style])

    def input_path(self, spec):
        raw = Path(spec['path'])
        path = (raw if raw.is_absolute() else self.base/raw).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = digest(path)
        if spec.get('sha256') is not None and spec['sha256'] != actual:
            raise ValueError(f'Input hash mismatch: {path}')
        self.inputs[str(path)] = actual
        return path

    def table(self, spec):
        columns, rows = spec['columns'], spec['rows']
        if not columns or not rows or len(columns) > 8:
            raise ValueError('Tables need 1-8 columns and at least one row')
        keys = [column['key'] for column in columns]
        if len(keys) != len(set(keys)):
            raise ValueError('Duplicate table columns')
        for row in rows:
            missing = set(keys)-set(row)
            if missing:
                raise ValueError(f'Missing table cells: {sorted(missing)}')
        weights = [float(column.get('width', 1)) for column in columns]
        if any(not math.isfinite(w) or w <= 0 for w in weights):
            raise ValueError('Column widths must be positive finite relative weights')
        widths = [WIDTH*w/sum(weights) for w in weights]
        cells = [[self.paragraph(column['label'], 'headcell') for column in columns]]
        cells += [[self.paragraph(row[key], 'cell') for key in keys] for row in rows]
        table = Table(cells, colWidths=widths, repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), PALE),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LINEBELOW', (0, 0), (-1, 0), .6, BLUE),
            ('LINEBELOW', (0, 1), (-1, -1), .25, colors.HexColor('#dce2e6')),
            ('LEFTPADDING', (0, 0), (-1, -1), 5),
            ('RIGHTPADDING', (0, 0), (-1, -1), 5),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ]))
        self.tables.append(dict(title=spec.get('title', ''), rows=len(rows), columns=len(keys),
                                data_cells=len(rows)*len(keys), keys=keys,
                                numeric_cells=sum(isinstance(row[k], (int, float)) and not isinstance(row[k], bool)
                                                  for row in rows for k in keys)))
        result = []
        if spec.get('title'):
            result.append(self.paragraph(spec['title'], 'heading'))
        result.extend((table, Spacer(1, 8)))
        result.extend(self.paragraph(note, 'small') for note in spec.get('notes', []))
        return result

    def image(self, spec, width, height, *, supplemental, table_cell=False):
        path = self.input_path(spec)
        image = Image(str(path), lazy=2)
        factor = min(width/image.imageWidth, height/image.imageHeight)
        image.drawWidth = image.imageWidth*factor
        image.drawHeight = image.imageHeight*factor
        label = spec.get('label', '')
        caption = spec.get('caption', '')
        if supplemental:
            label = 'Supplemental' + (': '+label if label else '')
        text = ' - '.join(part for part in (label, caption) if part)
        self.figures.append(dict(path=str(path), caption=caption, label=label,
                                 supplemental=supplemental, width_points=image.drawWidth,
                                 height_points=image.drawHeight))
        content = [image]
        if text:
            content.append(self.paragraph(text, 'caption'))
        return content if table_cell else KeepTogether(content)

    def section(self, spec, *, first=False):
        story = [] if first else [PageBreak()]
        story.append(self.paragraph(spec['title'], 'heading'))
        story.extend(self.paragraph(text) for text in spec.get('paragraphs', []))
        for table in spec.get('tables', []):
            story.extend(self.table(table))
        for figure in spec.get('figures', []):
            story.extend((Spacer(1, 5), self.image(figure, WIDTH, 285, supplemental=True)))
        return story

    def target_table(self):
        spec = self.manifest['target_fit']
        rows = []
        for row in spec['rows']:
            name = row['moment'] + (f" [{row['units']}]" if row.get('units') else '')
            rows.append(dict(row, moment=name))
        return dict(title=spec.get('title', 'Complete calibration fit'), rows=rows,
                    columns=[dict(key='moment', label='Moment and units', width=3.6)] +
                            [dict(key=key, label=label, width=1.15) for key, label in
                             [('target', 'Target'), ('model', 'Model'), ('gap', 'Gap'),
                              ('weight', 'Weight'), ('loss', 'Loss contribution')]],
                    notes=spec.get('notes', []))

    def parameter_table(self):
        spec = self.manifest['parameters']
        rows = []
        for row in spec['rows']:
            bounds = row.get('restriction', '')
            if row['estimated']:
                interval = '['+display(row['lower_bound'])+', '+display(row['upper_bound'])+']'
                bounds = interval + ('; '+bounds if bounds else '')
            rows.append(dict(parameter=row['parameter'], estimate=row['estimate'],
                             restriction=bounds, status=('Estimated' if row['estimated'] else 'Fixed / normalized'),
                             near_bound=row['near_bound']))
        return dict(title=spec.get('title', 'Parameters and restrictions'), rows=rows,
                    columns=[dict(key='parameter', label='Parameter', width=3.4),
                             dict(key='estimate', label='Estimate', width=1.2),
                             dict(key='restriction', label='Bounds / restriction', width=2.7),
                             dict(key='status', label='Status', width=1.5),
                             dict(key='near_bound', label='Near bound', width=1.2)],
                    notes=spec.get('notes', []))

    def appendix(self, groups):
        story = []
        for group_index, group in enumerate(groups, 1):
            figures = group['figures']
            layout = group.get('layout', 'two_per_page')
            if layout not in ('two_per_page', 'four_per_page', 'one_per_page'):
                raise ValueError('Unknown diagnostic appendix layout')
            per_page = {'one_per_page':1, 'two_per_page':2, 'four_per_page':4}[layout]
            for offset in range(0, len(figures), per_page):
                story.append(PageBreak())
                title = f"Standard diagnostic {group_index}: {group['title']}"
                if offset:
                    title += ' (continued)'
                story.append(self.paragraph(title, 'heading'))
                if group.get('caption'):
                    story.append(self.paragraph(group['caption'], 'small'))
                subset = figures[offset:offset+per_page]
                if per_page == 4:
                    # Explicit compact opt-in: the caller must inspect readability.
                    cells = [self.image(f, (WIDTH-14)/2, 245, supplemental=False,
                                        table_cell=True) for f in subset]
                    if len(cells)%2:
                        cells.append('')
                    grid = Table([cells[j:j+2] for j in range(0, len(cells), 2)],
                                 colWidths=[WIDTH/2, WIDTH/2])
                    grid.setStyle(TableStyle([('VALIGN', (0,0), (-1,-1), 'TOP'),
                                              ('LEFTPADDING',(0,0),(-1,-1),3),
                                              ('RIGHTPADDING',(0,0),(-1,-1),3)]))
                    story.append(grid)
                else:
                    for figure in subset:
                        story.extend((self.image(figure, WIDTH, 560 if per_page==1 else 265,
                                                 supplemental=False), Spacer(1, 12)))
        return story

    def page(self, canvas, document):
        self.pages = max(self.pages, canvas.getPageNumber())
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor('#d7dee4'))
        canvas.line(MARGIN, 35, PAGE_W-MARGIN, 35)
        canvas.setFont(self.normal, 8)
        canvas.setFillColor(GREY)
        label = self.manifest.get('footer', self.manifest['title'])
        while pdfmetrics.stringWidth(label, self.normal, 8) > WIDTH-45:
            label = label[:-4]+'...' if len(label)>4 else ''
        canvas.drawString(MARGIN, 23, label)
        canvas.drawRightString(PAGE_W-MARGIN, 23, str(canvas.getPageNumber()))
        canvas.restoreState()

    def build(self, output):
        for source in self.manifest.get('source_files', []):
            self.input_path(source)
        story = [self.paragraph(self.manifest['title'], 'title'),
                 self.paragraph(self.manifest['date'], 'small')]
        if self.manifest.get('subtitle'):
            story.append(self.paragraph(self.manifest['subtitle']))
        sections = self.manifest['sections']
        if sections:
            story.extend(self.section(sections[0], first=True))
        story.append(PageBreak())
        story.extend(self.table(self.target_table()))
        story.append(PageBreak())
        story.extend(self.table(self.parameter_table()))
        for section in sections[1:]:
            story.extend(self.section(section))
        story.extend(self.appendix(self.manifest.get('diagnostic_groups', [])))
        document = SimpleDocTemplate(str(output), pagesize=A4, leftMargin=MARGIN,
                                     rightMargin=MARGIN, topMargin=37, bottomMargin=48,
                                     title=self.manifest['title'], author=self.manifest.get('author', ''),
                                     allowSplitting=True)
        document.build(story, onFirstPage=self.page, onLaterPages=self.page)
        return dict(schema=SCHEMA, status='built_requires_visual_review', output=str(output),
                    output_sha256=digest(output), page_count=self.pages,
                    input_sha256=self.inputs, font_files=self.font_paths, tables=self.tables,
                    table_data_cell_count=sum(t['data_cells'] for t in self.tables),
                    target_fit_row_count=len(self.manifest['target_fit']['rows']),
                    free_parameter_count=sum(row['estimated'] for row in self.manifest['parameters']['rows']),
                    figures=self.figures, diagnostic_group_count=len(self.manifest.get('diagnostic_groups', [])),
                    numeric_display='Six significant digits for JSON numbers; supplied strings unchanged',
                    visual_review_completed=False,
                    verification_limits='No scientific recomputation; manifest provenance and visual review remain caller responsibilities')


def validate(manifest):
    if manifest.get('schema') != SCHEMA:
        raise ValueError('Wrong report-manifest schema')
    for key in ('title', 'date', 'sections', 'target_fit', 'parameters'):
        if key not in manifest:
            raise ValueError(f'Missing manifest field: {key}')
    targets = manifest['target_fit']['rows']
    if len(targets) != 12 or len({row['moment'] for row in targets}) != 12:
        raise ValueError('Exactly 12 distinct target rows are required')
    for row in targets:
        for key in ('moment', 'target', 'model', 'gap', 'weight', 'loss'):
            if key not in row:
                raise ValueError(f'Missing target field: {key}')
            display(row[key])
    parameters = manifest['parameters']['rows']
    if any(type(row.get('estimated')) is not bool for row in parameters):
        raise ValueError('Every parameter needs an explicit estimated boolean')
    if sum(row['estimated'] for row in parameters) != 11:
        raise ValueError('Exactly 11 estimated parameter rows are required')
    if len({row['parameter'] for row in parameters}) != len(parameters):
        raise ValueError('Parameter names must be unique')
    for row in parameters:
        for key in ('parameter', 'estimate', 'near_bound'):
            if key not in row:
                raise ValueError(f'Missing parameter field: {key}')
        if row['estimated'] and any(key not in row for key in ('lower_bound', 'upper_bound')):
            raise ValueError('Estimated parameters need both search bounds')
        if not row['estimated'] and not row.get('restriction'):
            raise ValueError('Fixed/normalized parameters need a stated restriction')
    groups = manifest.get('diagnostic_groups', [])
    expected = manifest.get('expected_diagnostic_groups', 17)
    if len(groups) != expected:
        raise ValueError(f'Expected {expected} standard diagnostic groups')
    if any(not group.get('figures') for group in groups):
        raise ValueError('Diagnostic groups cannot omit their figures')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path, output = args.manifest.resolve(), args.output.resolve()
    if output.suffix.lower() != '.pdf':
        parser.error('--output must name a PDF')
    if 'JMP_DS_draft' in output.parts:
        parser.error('The author-owned manuscript subtree is protected')
    if output.exists() or output.with_suffix('.qa.json').exists():
        raise FileExistsError('Refusing to overwrite an existing report or QA sidecar')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    validate(manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    qa = Builder(manifest, manifest_path).build(output)
    output.with_suffix('.qa.json').write_text(json.dumps(qa, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(output=str(output), page_count=qa['page_count'], status=qa['status'])))


if __name__ == '__main__':
    main()
