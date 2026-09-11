"""Build the utility/pension visual review from saved results only; no solver imports.

Run with the project/system Python with matplotlib and reportlab. Original diagnostics are verified
and embedded unchanged. Supplemental figures use saved CSV columns directly.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import time

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/fertility-visual-review-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader

ROOT = Path(__file__).resolve().parents[3]
DEFAULT = ROOT / 'output/model/e5f_matched_pf_20260909a'


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    base = args.base.resolve()
    out = (args.output or base / 'visual_review').resolve()
    out.mkdir(parents=True, exist_ok=True)
    path_root = base / 'social_security_repair/history_probe/root_round_02/completed_17370186'
    candidate = base / 'initial_calibration_contract/exact_loop_smoke/completed_17370427'
    decomp = base / 'utility_fiscal_decomposition'
    summary = json.loads((path_root / 'summary.json').read_text())
    if not summary['finite_horizon_market_fiscal_converged']:
        raise ValueError('Expected converged finite diagnostic path')
    checked = {}
    galleries = {}
    def register(label, folder, pins):
        images = []
        for relative, pin in pins.items():
            path = folder / relative
            if digest(path) != pin:
                raise ValueError(f'Graph hash mismatch: {path}')
            checked[str(path)] = pin
            images.append(path)
        if len(images) != 17:
            raise ValueError(f'Expected all 17 original graphs: {label}')
        galleries[label] = images
    register('2023 on the six-date transition (new utility; illustrative shock)', path_root,
             json.loads((path_root / 'diagnostics_receipt.json').read_text())['graph_sha256'])
    saved = json.loads((candidate / 'summary.json').read_text())
    register('Pre-2007 stationary candidate (new utility; two exact repetitions)', candidate / 'raw',
             {v['path']: v['sha256'] for v in saved['original_graphs']})
    pins = json.loads((decomp / 'original_graphs_receipt.json').read_text())
    for case, label in [('old_balanced', 'Old utility, balanced pensions'),
                        ('new_balanced', 'New utility, balanced pensions'),
                        ('old_old', 'Old utility, old pensions (fiscal diagnostic only)'),
                        ('new_old', 'New utility, old pensions (fiscal diagnostic only)')]:
        register('Pre-2007 fixed-parameter comparison: ' + label,
                 decomp / 'collected' / case / 'repetition_01',
                 {'standard_diagnostics/' + Path(v['path']).name: v['sha256'] for v in pins[case]})

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.grid': True,
                         'grid.alpha': .18, 'figure.dpi': 130})
    path = read_csv(path_root / 'transition_path.csv')
    years = [int(r['calendar_year']) for r in path]
    def values(key):
        return [float(r[key]) for r in path]
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 7.8), constrained_layout=True)
    plots = [
        ('payroll_tax_revenue', 'pension_outlays', 'Social Security: revenue and spending', 'Model period units'),
        ('pension_period_units', None, 'Pension benefit per retiree', 'Model period units'),
        ('housing_demand', 'housing_supply', 'Housing demand and supply', 'Housing-service units'),
        ('asset_price', None, 'House asset price', 'Model price units'),
        ('birth_children_topcode_adjusted', None, 'Birth flow (top-parity adjustment)', 'Children per model period'),
        ('owner_rate', None, 'Homeownership', 'Share of household heads'),
    ]
    export = []
    for ax, (key, second, title, ylabel) in zip(axes.flat, plots):
        first_label = {'payroll_tax_revenue': 'Payroll revenue', 'housing_demand': 'Demand'}.get(key, title)
        ax.plot(years, values(key), 'o-', color='#176b94', label=first_label, lw=2)
        if second:
            ax.plot(years, values(second), 'x--', color='#d77a28',
                    label={'pension_outlays': 'Pension spending', 'housing_supply': 'Supply'}[second], lw=1.8)
            ax.legend(fontsize=9)
        ax.set(title=title, ylabel=ylabel, xlabel='Year', xticks=years[::2])
        export.append({'figure': 'transition', 'series': key, 'x': years, 'y': values(key)})
        if second:
            export.append({'figure': 'transition', 'series': second, 'x': years, 'y': values(second)})
    fig.suptitle('New utility and balanced pensions: solved six-date transition', fontsize=17)
    fig.get_layout_engine().set(rect=(0, .075, 1, .94))
    fig.text(.025, .018, 'Illustrative announced preference shock; no fitted historical path or horizon certification.\n'
             'These are model outcomes, not policy effects or population forecasts. Benefits and quantities use four-year model periods.', fontsize=10)
    transition = out / 'transition_overview.png'
    fig.savefig(transition); plt.close(fig)

    old = read_csv(decomp / 'collected/old_balanced/repetition_01/lifecycle_2023.csv')
    new = read_csv(decomp / 'collected/new_balanced/repetition_01/lifecycle_2023.csv')
    if [r['age_node'] for r in old] != [r['age_node'] for r in new]:
        raise ValueError('Static comparison age grids differ')
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 7.8), constrained_layout=True)
    for ax, (key, title, unit) in zip(axes.flat, [
            ('owner_rate', 'Ownership by age', 'Ownership share'),
            ('mean_rooms', 'Housing by age', 'Mean occupied rooms'),
            ('mean_liquid_wealth', 'Liquid wealth by age', 'Model wealth units'),
            ('childless_rate', 'Childlessness by age', 'Childless share')]):
        for rows, label, color, style in [(old, 'Old utility', '#176b94', '-'),
                                          (new, 'New utility', '#d77a28', '--')]:
            xx = [float(r['age_node']) for r in rows]
            yy = [float(r[key]) for r in rows]
            ax.plot(xx, yy, style, color=color, lw=2, label=label)
            export.append({'figure': 'static', 'series': label + ':' + key, 'x': xx, 'y': yy})
        ax.set(title=title, xlabel='Age', ylabel=unit)
        ax.legend()
    fig.suptitle('Direct utility comparison: both pension budgets balanced', fontsize=17)
    fig.get_layout_engine().set(rect=(0, .075, 1, .94))
    fig.text(.025, .018, 'Same structural parameters and fixed fertility preference; separate stationary housing equilibria.\n'
             'Pre-2007 diagnostic. A matched 2023 recalibration is a separate experiment.', fontsize=10)
    static = out / 'utility_comparison.png'
    fig.savefig(static); plt.close(fig)

    late = decomp / 'static_2023/results'
    late_figure = None
    if (late / 'collection_manifest.json').exists():
        late_pins = json.loads((late / 'collection_manifest.json').read_text())['sha256']
        for name, pin in late_pins.items():
            if digest(late / name) != pin:
                raise ValueError('Late comparison artifact changed: ' + name)
        fig, axes = plt.subplots(2, 2, figsize=(13.8, 7.8), constrained_layout=True)
        for case, label, color, style in [('old_balanced', 'Old utility', '#176b94', '-'),
                                          ('new_balanced', 'New utility', '#d77a28', '--')]:
            prefix = case + '/repetition_02/'
            register('Stationary counterpart at the retained 2023 preference: ' + label,
                     late / case / 'repetition_02',
                     {k[len(prefix):]: v for k,v in late_pins.items()
                      if k.startswith(prefix + 'standard_diagnostics/') and k.endswith('.png')})
            source = late / case / 'repetition_02/lifecycle_2023.csv'
            rows = read_csv(source)
            for ax, (key, title, unit) in zip(axes.flat, [
                    ('owner_rate', 'Ownership by age', 'Ownership share'),
                    ('mean_rooms', 'Housing by age', 'Mean occupied rooms'),
                    ('mean_liquid_wealth', 'Liquid wealth by age', 'Model wealth units'),
                    ('childless_rate', 'Childlessness by age', 'Childless share')]):
                xx = [float(r['age_node']) for r in rows]
                yy = [float(r[key]) for r in rows]
                ax.plot(xx, yy, style, color=color, lw=2, label=label)
                ax.set(title=title, xlabel='Age', ylabel=unit)
                ax.legend()
                export.append({'figure': 'static_2023', 'series': label + ':' + key,
                               'x': xx, 'y': yy, 'source_csv': str(source)})
        fig.suptitle('Old and new utility at the retained 2023 preference', fontsize=17)
        fig.get_layout_engine().set(rect=(0, .075, 1, .94))
        fig.text(.025, .018, 'Balanced pensions in both arms; same structural parameters and supply curve; two exact repetitions each.\n'
                 'Stationary counterparts, not recalibrated 2023 economies. They do not reproduce the historical 2023 population.', fontsize=10)
        late_figure = out / 'utility_comparison_2023_preference.png'
        fig.savefig(late_figure); plt.close(fig)

    pages = [(transition, '1. Transition and pension accounts', 'Six dates, 2007-2027; new utility; illustrative shock; terminal horizon remains unchecked.'),
             (static, '2. Direct effect of the utility change', 'Pre-2007 stationary diagnostics; same parameters and preference; balanced pensions in both cases.')]
    for fname, title in [('policy_childless_renter_age30.png', '3. 2023 policies: childless renter, age 30'),
                         ('policy_childless_renter_age42.png', '4. 2023 policies: childless renter, age 42'),
                         ('housing_by_age_income_state.png', '5. Housing across ages and income states'),
                         ('fertility_policy_by_age_income_state.png', '6. Fertility across ages and income states')]:
        pages.append((path_root / 'standard_diagnostics' / fname, title,
                      '2023 slice of the six-date transition; original saved diagnostic, unchanged.'))
    for case, title in [('old_balanced', '7. Fixed-parameter stationary policies: old utility'),
                        ('new_balanced', '8. Fixed-parameter stationary policies: new utility')]:
        pages.append((decomp / 'collected' / case / 'repetition_01/standard_diagnostics/policy_childless_renter_age30.png',
                      title, 'Pre-2007 stationary comparison; balanced pensions; original saved diagnostic, unchanged.'))
    if late_figure:
        pages.append((late_figure, '9. Stationary comparison at the retained 2023 preference',
                      'Both budgets balanced; unchanged structural parameters; no re-estimation or historical-population reproduction.'))
    pdf = ROOT / 'output/pdf/e5f_utility_pension_visual_review.pdf'
    pdf.parent.mkdir(parents=True, exist_ok=True)
    size = landscape(A4)
    doc = canvas.Canvas(str(pdf), pagesize=size)
    doc.setTitle('Utility and pension visual review')
    for number, (p, title, note) in enumerate(pages, 1):
        w, h = size
        doc.setFont('Helvetica-Bold', 15); doc.drawString(26, h - 29, title)
        doc.setFont('Helvetica', 9); doc.drawString(26, h - 45, note)
        img = ImageReader(str(p)); iw, ih = img.getSize()
        scale = min((w - 52) / iw, (h - 83) / ih)
        doc.drawImage(img, (w - iw * scale) / 2, 24 + (h - 83 - ih * scale) / 2,
                      width=iw * scale, height=ih * scale)
        doc.setFont('Helvetica', 8); doc.drawString(26, 12, f'Saved-result review | {number}/{len(pages)} | Full 17-graph sets in the accompanying gallery')
        doc.showPage()
    doc.save()

    def link(p):
        return html.escape(os.path.relpath(p, out))
    sections = []
    for label, images in galleries.items():
        cards = ''.join(f'<figure><figcaption>{html.escape(p.stem.replace("_", " "))}</figcaption>'
                        f'<a href="{link(p)}"><img loading="lazy" src="{link(p)}"></a></figure>' for p in images)
        sections.append(f'<details><summary>{html.escape(label)} - all 17 original graphs</summary>{cards}</details>')
    doc_html = '''<!doctype html><html><meta charset="utf-8"><title>Utility and pension visual review</title>
<style>body{font:17px system-ui;max-width:1180px;margin:32px auto;padding:0 20px;color:#172b3a;background:#f5f7f9}h1{font-size:32px}p{line-height:1.5}img{width:100%;background:white}figure{margin:20px 0;background:white;padding:12px;border:1px solid #dde3e8}figcaption,summary{font-weight:650;padding:12px}details{background:white;margin:16px 0;border:1px solid #ccd6df}summary{cursor:pointer}a{color:#176b94}.note{padding:16px;border-left:4px solid #d77a28;background:#fff8ed}</style>
<h1>Utility, pensions and model behavior</h1>
<p class="note">The six-date transition clears housing and pension accounts. Its preference path is illustrative and its terminal horizon is not certified. The old/new stationary comparison fixes parameters and preferences; it is not a matched 2023 recalibration.</p>
''' + f'<p><a href="{link(pdf)}">Open the {len(pages)}-page PDF</a> · Click any original graph to enlarge it.</p>'
    doc_html += f'<h2>Transition and pension accounts</h2><img src="{transition.name}"><h2>Direct utility comparison</h2><img src="{static.name}">'
    if late_figure:
        doc_html += f'<h2>Stationary counterparts at the 2023 preference</h2><img src="{late_figure.name}"><p>These are fixed-parameter diagnostics. Complete old target/weight tables and parameter restrictions:</p><ul>'
        for case in ['old_balanced','new_balanced']:
            doc_html += f'<li>{case}: <a href="{link(late/case/"target_fit.csv")}">All targets, gaps, weights and contributions</a> · <a href="{link(late/case/"repetition_02/parameters.csv")}">All parameters and restrictions</a></li>'
        doc_html += '</ul>'
    doc_html += '<h2>Policy functions in 2023</h2>' + ''.join(
        f'<figure><figcaption>{html.escape(title)}</figcaption><img src="{link(p)}"></figure>' for p, title, _ in pages[2:4])
    doc_html += '<h2>Complete saved diagnostic sets</h2>' + ''.join(sections)
    doc_html += '<p>All original plots retain their original filenames, data and bytes. No model solve is performed to build this review.</p></html>'
    (out / 'index.html').write_text(doc_html)
    (out / 'plotted_series.json').write_text(json.dumps(export, indent=2))
    receipt = {'source_graphs_verified': len(checked), 'pdf_pages': len(pages),
               'model_solves': 0, 'elapsed_seconds': time.monotonic() - started,
               'source_graph_sha256': checked,
               'source_tables_sha256': {str(p): digest(p) for p in [path_root / 'transition_path.csv',
                    decomp / 'collected/old_balanced/repetition_01/lifecycle_2023.csv',
                    decomp / 'collected/new_balanced/repetition_01/lifecycle_2023.csv']},
               'pdf_sha256': digest(pdf), 'pdf': str(pdf)}
    (out / 'build_receipt.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_graph_sha256','source_tables_sha256']}))


if __name__ == '__main__':
    main()
