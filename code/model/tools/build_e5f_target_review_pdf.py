#!/usr/bin/env python3
"""Render the reviewed target catalogue; never alter a calibration contract."""
import argparse
import json
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, PageBreak


def build(packet, output):
    packet = Path(packet)
    review = json.loads((packet / 'review.json').read_text())
    rows = review['rows']
    assert len(rows) == 13 and len({r['id'] for r in rows}) == 13
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('SmallReview', fontName='Helvetica', fontSize=10, leading=13, spaceAfter=6))
    styles.add(ParagraphStyle('ReviewTitle', fontName='Helvetica-Bold', fontSize=19, leading=23, textColor=colors.HexColor('#17324d'), spaceAfter=12))
    styles.add(ParagraphStyle('ReviewHeading', fontName='Helvetica-Bold', fontSize=11, leading=14, textColor=colors.HexColor('#17324d'), spaceBefore=9, spaceAfter=5))
    def p(text, style='SmallReview'):
        return Paragraph(escape(str(text)), styles[style])
    def table(data, widths):
        t = Table([[p(x) for x in row] for row in data], colWidths=widths, repeatRows=1, hAlign='LEFT')
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#dce8f0')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f3f6f8')]),
            ('GRID', (0, 0), (-1, -1), .3, colors.HexColor('#bcc9d3')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 5), ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        return t
    story = [p('Calibration targets: what is verified?', 'ReviewTitle'), p(review['date'])]
    for line in review['summary']:
        story.append(p(line))
    story += [p('Complete current target list', 'ReviewHeading'), table(
        [['Moment', 'Current target', 'Remaining decision']] +
        [[r['label'], r['display_target'], r['short_issue']] for r in rows], [205, 83, 236])]
    for group, heading in [('housing', 'Housing: source and model measurement'), ('fertility', 'Fertility: source and model measurement'), ('wealth', 'Wealth: source and model measurement')]:
        story += [PageBreak(), p(heading, 'ReviewTitle')]
        for r in [x for x in rows if x['group'] == group]:
            story += [p(r['label'] + ' | ' + r['display_target'], 'ReviewHeading'),
                      p('Data: ' + r['data']), p('Model: ' + r['model']),
                      p('Assessment: ' + r['assessment'])]
    story += [PageBreak(), p('Weights and decisions before launch', 'ReviewTitle'),
              p('Current objective: sum of twelve squared gaps divided by fixed working scales squared. The separate 2.1 normalization is unscored. These scales are not all sampling standard errors; cross-moment covariance is not used.'),
              table([['Moment', 'Working scale', 'Basis']] + [[r['label'], '-' if r['working_scale'] is None else f"{r['working_scale']:.6g}", r['weight_basis']] for r in rows], [207, 87, 230])]
    story.append(p('Recommended closure order', 'ReviewHeading'))
    for i, line in enumerate(review['next_steps'], 1):
        story.append(p(f'{i}. {line}'))
    story += [PageBreak(), p('Evidence and limits of this review', 'ReviewTitle')]
    for label, source in review['sources']:
        story += [p(label, 'ReviewHeading'), p(source)]
    story.append(p('This review reads saved empirical receipts and the actual scored observer/source. It does not rerun microdata, change estimates or weights, certify structural identification, or submit a model job. Detailed per-row receipts and the complete target catalogue accompany this PDF.'))
    def footer(canvas, doc):
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#52616b'))
        canvas.drawString(44, 25, 'Target review | no target changes or calibration submission')
        canvas.drawRightString(568, 25, str(doc.page))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), pagesize=letter, leftMargin=44, rightMargin=44, topMargin=40, bottomMargin=42,
                      title='Calibration targets: verification and decisions').build(story, onFirstPage=footer, onLaterPages=footer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--packet', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    build(args.packet, args.output)
