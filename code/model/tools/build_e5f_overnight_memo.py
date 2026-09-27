#!/usr/bin/env python3
"""Read saved calibration cases and build an exactly two-page overnight memo.

No model imports or solves. The index has ``primary_objective``, optional
``controllers`` and ``cases`` arrays (each entry has path and weighting), plus
optional status_note/review_notes/next_steps/expected_workers. Relative paths
are relative to the index. Controllers are read through complete.json or
checkpoint.json; explicit cases require status='success'. The primary system
must have 14 target rows and nine outer-search parameter restrictions. The
tenth fitted parameter is psi_child, normalized to completed fertility 2.1.

Example:
  python build_e5f_overnight_memo.py --index run/report_index.json --output run/memo

Writes memo.pdf, summary.json, monitor.json, candidates.csv, target_fit.csv,
parameters.csv, and a small supplemental fit_overview.svg. The tables and
receipts are checked; this collector does not recertify checkpoint arrays.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape


MOMENTS = {
    'initial_normalization': 'Completed fertility (normalization)',
    'cps_childlessness': 'Childlessness, ages 40-44',
    'cps_exactly_one': 'Exactly one child, among mothers',
    'nchs_mean_age': 'Mean age at first birth',
    'nchs_share30': 'First births at age 30 or older',
    'early_fertility': 'Children ever born at 25, capped at 3',
    'wealth_earnings': 'Aggregate wealth / annual earnings',
    'bequest_wealth': 'Annual bequest flow / wealth',
    'old_dispersion': 'Older wealth/income: p90 / median',
    'mean_rooms': 'Mean occupied rooms (AHS)',
    'ownership_30_55': 'Ownership, ages 30-55',
    'first_birth_rooms': 'Housing response to first birth',
    'family_rooms': 'Rooms: 3+ vs 1-2 children at home',
    'recent_parent_ownership': 'Recent-parent ownership gap',
}
PARAMETERS = {
    'H0': 'Housing supply level', 'beta_annual': 'Annual discount factor',
    'chi': 'Owner housing-service premium',
    'first_birth_fixed_cost': 'First-birth fixed cost',
    'kappa_fert': 'First-birth taste scale',
    'kappa_fert_continuation': 'Later-birth taste scale',
    'theta0': 'Bequest strength', 'delta_alpha_jump': 'First-child housing loading',
    'child_benefit_curvature': 'Child-benefit curvature',
    'psi_child': 'Child-benefit level (normalized)',
}
SUCCESS = 'verified_provisional_calibration_point'


def read(path):
    return json.loads(Path(path).read_text())


def resolve(value, base):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Nonfinite number')
    return value


def close(a, b):
    return math.isclose(number(a), number(b), rel_tol=1e-10, abs_tol=1e-12)


def rows(path, key):
    with Path(path).open(newline='') as stream:
        values = list(csv.DictReader(stream))
    result = {r[key]: r for r in values}
    if len(result) != len(values):
        raise ValueError(f'Duplicate {key} in {path}')
    return result


def fmt(value):
    if value is None or value == '':
        return '--'
    value = float(value)
    if value and (abs(value) < .001 or abs(value) >= 10000):
        return f'{value:.2e}'
    return f'{value:.3f}'.rstrip('0').rstrip('.') if value else '0'


def target_spec(objective):
    targets = objective['target_rows']
    names = [r['restriction_id'] for r in targets]
    if len(names) != 14 or set(names) != set(MOMENTS):
        raise ValueError('Expected the complete fourteen-row target registry')
    for row in targets:
        number(row['target'])
        if row['restriction_id'] == 'initial_normalization':
            if row['actual_weight'] is not None or not close(row['target'], 2.1):
                raise ValueError('Expected unscored completed-fertility normalization 2.1')
        elif number(row['actual_weight']) <= 0:
            raise ValueError('Primary weights must be positive')
    restrictions = objective['parameter_restrictions']
    if len(restrictions) != 9 or {r['parameter'] for r in restrictions} != set(PARAMETERS) - {'psi_child'}:
        raise ValueError('Expected all nine outer-search parameter restrictions')
    return targets, restrictions


def gather(index, base):
    """Read explicitly indexed controllers only; never recursively scan outputs."""
    entries, states, errors = [], [], []
    for item in index.get('controllers', []):
        folder = resolve(item['path'], base)
        state = {'path': str(folder), 'weighting': item.get('weighting', 'primary')}
        try:
            source = next((folder / f for f in ('complete.json', 'checkpoint.json') if (folder / f).exists()), None)
            payload = read(source) if source else {}
            state.update(status=payload.get('status', 'running_or_pending'),
                         error=payload.get('error'), records=len(payload.get('records', [])))
            if payload.get('error'):
                errors.append({'path': str(folder), 'status': 'controller_stop', 'error': str(payload['error'])})
            heartbeat = folder / 'heartbeat.json'
            if heartbeat.exists():
                hb = read(heartbeat)
                state['heartbeat'] = hb
                state['heartbeat_age_seconds'] = max(0, datetime.now(timezone.utc).timestamp() - number(hb['epoch']))
                state['stale_30_minutes'] = state['heartbeat_age_seconds'] > 1800 and not (folder / 'complete.json').exists()
            for record in payload.get('records', []):
                candidate = dict(record, weighting=item.get('weighting', 'primary'))
                candidate['path'] = record.get('case_path', str(folder / record['case'] / 'case'))
                candidate['controller'] = str(folder)
                entries.append(candidate)
            states.append(state)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            errors.append({'path': str(folder), 'error': str(exc), 'status': 'collector_error'})
    entries.extend(dict(item) for item in index.get('cases', []))
    unique = {}
    for entry in entries:
        path = resolve(entry['path'], base)
        if (path / 'case' / 'receipt.json').exists():
            path = path / 'case'
        entry['path'] = str(path)
        entry.setdefault('weighting', 'primary')
        key = str(path)
        if key in unique:
            prior = unique[key]
            if prior.get('status') != entry.get('status') or prior['weighting'] != entry['weighting']:
                errors.append({'path': key, 'error': 'Conflicting duplicate index entries', 'status': 'collector_error'})
            continue
        unique[key] = entry
    return list(unique.values()), states, errors


def validate_case(entry, targets, restrictions):
    """Verify complete finite tables and rescore with common primary weights."""
    if entry['weighting'] not in ('primary', 'identity', 'early_fertility_3000'):
        raise ValueError('Unknown weighting experiment')
    path = Path(entry['path'])
    receipt = read(path / 'receipt.json')
    if receipt.get('status') != SUCCESS:
        raise ValueError('Receipt is not a completed provisional calibration point')
    if entry.get('receipt_sha256'):
        actual = hashlib.sha256((path / 'receipt.json').read_bytes()).hexdigest()
        if actual != entry['receipt_sha256']:
            raise ValueError('Receipt differs from controller completion record')
    fits, parameters = rows(path / 'target_fit.csv', 'moment'), rows(path / 'parameters.csv', 'parameter')
    if set(fits) != set(MOMENTS) or not set(PARAMETERS) <= set(parameters):
        raise ValueError('Incomplete target or fitted-parameter table')
    primary_fit, raw_loss = [], 0.
    for target in targets:
        name = target['restriction_id']; row = fits[name]
        value = number(row['model']); gap = value - number(target['target'])
        if not close(row['target'], target['target']) or not close(row['gap'], gap):
            raise ValueError(f'Incompatible target or incorrect gap: {name}')
        weight = target['actual_weight']
        if weight is None:
            if row['weight'] or row['loss_contribution']:
                raise ValueError('Normalization must remain unscored')
            if abs(gap) > .0005 + 1e-12:
                raise ValueError('Completed fertility misses the retained normalization gate')
        else:
            own_weight = number(row['weight'])
            if own_weight <= 0 or not close(row['loss_contribution'], own_weight * gap * gap):
                raise ValueError(f'Invalid own-system loss contribution: {name}')
            if entry['weighting'] == 'primary' and not close(own_weight, weight):
                raise ValueError('Primary-labelled case has different weights')
            if entry['weighting'] == 'identity' and own_weight != 1:
                raise ValueError('Identity-labelled case does not have unit weights')
            if entry['weighting'] == 'early_fertility_3000':
                expected = 3000. if name == 'early_fertility' else weight
                if not close(own_weight, expected):
                    raise ValueError('Early-fertility experiment has unexpected weights')
            raw_loss += number(row['loss_contribution'])
        primary_fit.append(dict(moment=name, label=MOMENTS[name], target=number(target['target']),
                                model=value, gap=gap, weight=weight,
                                loss_contribution=None if weight is None else number(weight) * gap * gap))
    if not close(raw_loss, receipt['loss']) or (entry.get('loss') is not None and not close(raw_loss, entry['loss'])):
        raise ValueError('Saved scalar loss differs from complete target table')
    fitted = []
    for restriction in restrictions:
        name = restriction['parameter']; row = parameters[name]
        estimate, lower, upper = map(number, (row['estimate'], restriction['lower'], restriction['upper']))
        if not lower < upper or not close(row['lower'], lower) or not close(row['upper'], upper):
            raise ValueError(f'Incompatible search bounds: {name}')
        if estimate < lower - 1e-12 or estimate > upper + 1e-12:
            raise ValueError(f'Estimate outside search bounds: {name}')
        near = min(estimate - lower, upper - estimate) <= .01 * (upper - lower)
        if row['near_bound'].lower() != str(near).lower():
            raise ValueError(f'Incorrect near-bound flag: {name}')
        fitted.append(dict(parameter=name, label=PARAMETERS[name], estimate=estimate,
                           lower=lower, upper=upper, near_bound=near, restriction='searched'))
    benefit = number(parameters['psi_child']['estimate'])
    if benefit <= 0 or not close(benefit, receipt['normalization']['psi_child']):
        raise ValueError('Invalid normalized child-benefit level')
    fitted.append(dict(parameter='psi_child', label=PARAMETERS['psi_child'], estimate=benefit,
                       lower=0., upper=None, near_bound=None, restriction='positive; fertility = 2.1'))
    solves = read(path / 'stationary_solves.json')
    if len(solves) != receipt['objective_stationary_solves'] or any(s.get('status') != 'completed' for s in solves):
        raise ValueError('Stationary-solve ledger is incomplete')
    # A weight experiment must retain fixed model objects and the entrant/input
    # checkpoint. Endogenous fiscal values and the derived benefit coefficient
    # may legitimately vary across estimated preference candidates.
    endogenous = {'child_benefit_CRRA_coefficient', 'payroll_tax', 'pension_period'}
    fixed = {name: number(row['estimate']) for name, row in parameters.items()
             if name not in PARAMETERS and name not in endogenous}
    comparison_inputs = dict(fixed_parameters=fixed,
                             inherited_checkpoint=receipt.get('selected_checkpoint_sha256'),
                             owner_grid=receipt.get('retained_owner_grid'),
                             conception_schedule=receipt.get('retained_conception_schedule'))
    diagnostic_folder = path / 'standard_diagnostics'
    if entry.get('controller'):
        export = Path(entry['controller']) / 'selected_export'
        if (export / 'export_receipt.json').exists():
            exported = read(export / 'export_receipt.json')
            selected = exported.get('selected', {})
            if selected.get('case_path') and Path(selected['case_path']).resolve() == path.resolve():
                if (exported.get('status') not in ('verified_selected_export', 'verified_selected_export_from_interrupted_search')
                        or selected.get('receipt_sha256') != hashlib.sha256((path / 'receipt.json').read_bytes()).hexdigest()
                        or selected.get('checkpoint_sha256') != receipt.get('case_checkpoint_sha256')):
                    raise ValueError('Selected diagnostic export does not authenticate this case')
                diagnostic_folder = export / 'standard_diagnostics'
    case = entry.get('case', '')
    phase = ('search' if case.startswith(('initial_', 'de_')) else
             'repeat' if case.startswith('repeat_') else
             'acceptance' if case.startswith('smoke_') else 'other')
    return dict(path=str(path), weighting=entry['weighting'], phase=phase, own_loss=raw_loss,
                primary_loss=sum(r['loss_contribution'] or 0. for r in primary_fit),
                fits=primary_fit, parameters=fitted,
                source_manifest_sha256=receipt.get('source_manifest_sha256'),
                comparison_inputs=comparison_inputs,
                target_weight_fingerprint=receipt.get('target_weight_fingerprint'),
                receipt_sha256=hashlib.sha256((path / 'receipt.json').read_bytes()).hexdigest(),
                stationary_solves=len(solves), solve_seconds=receipt.get('objective_stationary_solve_seconds'),
                diagnostic_pngs=len(list(diagnostic_folder.glob('*.png'))),
                diagnostic_source=str(diagnostic_folder),
                market_residual=receipt.get('market_residual'), estate_funding=receipt.get('estate_funding'),
                observer_warnings=receipt.get('model_observer_warnings', []),
                checkpoint_sha256=receipt.get('case_checkpoint_sha256'),
                status='table_checked_provisional_point')


def collect(index_path):
    index_path = Path(index_path).resolve(); index = read(index_path)
    objective_path = resolve(index['primary_objective'], index_path.parent)
    objective = read(objective_path); targets, restrictions = target_spec(objective)
    entries, controllers, errors = gather(index, index_path.parent)
    counts = Counter(); accepted = []
    for entry in entries:
        status = entry.get('status', 'pending')
        if status != 'success':
            counts[status] += 1
            if entry.get('error'):
                errors.append({'path': entry['path'], 'error': entry['error'], 'status': status})
            continue
        try:
            accepted.append(validate_case(entry, targets, restrictions))
            counts['success'] += 1
        except (OSError, ValueError, TypeError, KeyError) as exc:
            counts['collector_rejected'] += 1
            errors.append({'path': entry['path'], 'status': 'collector_rejected', 'error': str(exc)})
    expected_source = index.get('source_manifest_sha256')
    primary = [x for x in accepted if x['weighting'] == 'primary']
    if not expected_source and primary:
        expected_source = primary[0]['source_manifest_sha256']
    expected_inputs = primary[0]['comparison_inputs'] if primary else None
    # Different code cannot be compared silently, even if target names coincide.
    compatible = []
    for candidate in accepted:
        if (not expected_source or candidate['source_manifest_sha256'] != expected_source
                or expected_inputs is not None and candidate['comparison_inputs'] != expected_inputs):
            errors.append({'path': candidate['path'], 'status': 'comparison_excluded',
                           'error': 'Common primary source/fixed-input identity is missing or differs'})
            counts['comparison_excluded'] += 1
        else:
            compatible.append(candidate)
    primary = [x for x in compatible if x['weighting'] == 'primary']
    best = min(primary, key=lambda x: x['primary_loss']) if primary else None
    groups = {}
    for name in sorted({e['weighting'] for e in entries} |
                       {state['weighting'] for state in controllers} | {'primary'}):
        candidates = [x for x in compatible if x['weighting'] == name]
        winner = min(candidates, key=lambda x: x['primary_loss']) if candidates else None
        own_winner = min(candidates, key=lambda x: x['own_loss']) if candidates else None
        searched = sum(c['phase'] == 'search' for c in candidates)
        groups[name] = dict(completed=len(candidates), search_completed=searched,
                            verification_or_other=len(candidates)-searched,
                            best_under_primary_weights=winner['primary_loss'] if winner else None,
                            best_path=winner['path'] if winner else None,
                            own_winner_primary_loss=own_winner['primary_loss'] if own_winner else None,
                            own_winner_path=own_winner['path'] if own_winner else None)
    if best:
        fits, parameters = best['fits'], best['parameters']
    else:
        fits = [dict(moment=r['restriction_id'], label=MOMENTS[r['restriction_id']], target=number(r['target']),
                     model=None, gap=None, weight=r['actual_weight'], loss_contribution=None) for r in targets]
        parameters = [dict(parameter=r['parameter'], label=PARAMETERS[r['parameter']], estimate=None,
                           lower=number(r['lower']), upper=number(r['upper']), near_bound=None, restriction='searched') for r in restrictions]
        parameters.append(dict(parameter='psi_child', label=PARAMETERS['psi_child'], estimate=None,
                               lower=0., upper=None, near_bound=None, restriction='positive; fertility = 2.1'))
    misses = sorted([r for r in fits if r['loss_contribution'] is not None], key=lambda x: x['loss_contribution'], reverse=True)[:3]
    timings = [number(c['solve_seconds']) for c in compatible if c['solve_seconds'] is not None]
    return dict(generated_utc=datetime.now(timezone.utc).isoformat(), index_path=str(index_path),
                primary_objective_path=str(objective_path), primary_objective_sha256=hashlib.sha256(objective_path.read_bytes()).hexdigest(),
                status='provisional_results' if best else 'no_completed_primary_point',
                status_note=index.get('status_note', ''), lead_note=index.get('lead_note', ''),
                expected_workers=index.get('expected_workers', {}),
                counts=dict(counts), controller_states=controllers, errors=errors,
                groups=groups, selected=best, candidates=compatible, target_fit=fits, parameters=parameters,
                largest_weighted_misses=misses, median_stationary_solve_seconds=statistics.median(timings) if timings else None,
                review_notes=index.get('review_notes', []), next_steps=index.get('next_steps', []),
                repeats=index.get('repeats', {}), benchmark=index.get('benchmark', {}),
                limitations=['Table/receipt checks do not re-certify checkpoint arrays or establish identification.',
                             'All cross-weight comparisons use primary weights on the same fourteen targets.',
                             'Best observed primary point is provisional; global or local optimality is not certified.'])


def write_csv(path, data, fields):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(data)


def write_overview(path, summary):
    """A supplemental weighted-gap view; never replaces the seventeen diagnostics."""
    fits = [r for r in summary['target_fit'] if r['weight'] is not None and r['gap'] is not None]
    width, height = 820, 90 + 29 * max(1, len(fits)); center = 555
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           '<text x="20" y="28" font-family="sans-serif" font-size="17" fill="#172b43">Supplemental: signed gaps under primary weights</text>',
           '<text x="20" y="50" font-family="sans-serif" font-size="12">Bar = (model - target) x sqrt(weight). These are objective units, not statistical z-scores.</text>']
    if fits:
        extent = max(1., max(abs(r['gap'] * math.sqrt(r['weight'])) for r in fits))
        svg.append(f'<line x1="{center}" x2="{center}" y1="68" y2="{height-15}" stroke="#777"/>')
        for i, row in enumerate(fits):
            y = 82 + i * 29; value = row['gap'] * math.sqrt(row['weight']); length = 190 * value / extent
            svg.extend([f'<text x="20" y="{y}" font-family="sans-serif" font-size="12">{escape(row["label"])}</text>',
                        f'<rect x="{min(center, center+length):.3f}" y="{y-12}" width="{max(.3, abs(length)):.3f}" height="16" fill="{"#b75837" if value>0 else "#307e8c"}"/>',
                        f'<text x="775" y="{y}" text-anchor="end" font-family="sans-serif" font-size="11">{fmt(value)}</text>'])
    else:
        svg.append('<text x="20" y="85" font-family="sans-serif" font-size="14">No completed primary point. No model values are imputed.</text>')
    path.write_text('\n'.join(svg + ['</svg>']))


def render_pdf(path, summary):
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.pdfgen import canvas
    from reportlab.platypus import Paragraph, Table, TableStyle

    width, height = 612, 792; left, content_width = 38, 536
    navy, teal, muted = colors.HexColor('#172b43'), colors.HexColor('#167884'), colors.HexColor('#5e6c7b')
    pdf = canvas.Canvas(str(path), pagesize=(width, height))
    pdf.setTitle('Fertility and housing: overnight calibration memo')
    body = ParagraphStyle('body', fontName='Helvetica', fontSize=9, leading=12, textColor=navy, alignment=TA_LEFT)
    small = ParagraphStyle('small', parent=body, fontSize=8, leading=10)
    cell = ParagraphStyle('cell', parent=body, fontSize=8, leading=10)

    def paragraph(text, y, style=body, max_height=100):
        p = Paragraph(text, style); _, h = p.wrap(content_width, max_height)
        if h > max_height:
            raise ValueError('Memo paragraph exceeds its allocated space; shorten index commentary')
        p.drawOn(pdf, left, y-h)
        return y-h-8

    def heading(text, y):
        pdf.setFillColor(teal); pdf.setFont('Helvetica-Bold', 10)
        pdf.drawString(left, y-11, text)
        return y-20

    def table(data, widths, y):
        rendered = [[Paragraph(escape(str(c)), cell) for c in row] for row in data]
        t = Table(rendered, colWidths=widths, hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f0f3')),
                               ('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),
                               ('BOTTOMPADDING',(0,0),(-1,-1),4),('LEFTPADDING',(0,0),(-1,-1),5),
                               ('RIGHTPADDING',(0,0),(-1,-1),5),
                               ('LINEBELOW',(0,0),(-1,0),.6,teal),
                               ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')])]))
        _, h = t.wrap(content_width, height)
        if y-h < 46:
            raise ValueError('Memo table would overflow two-page layout')
        t.drawOn(pdf,left,y-h)
        return y-h-9

    def header(page, subtitle):
        pdf.setFillColor(navy); pdf.setFont('Helvetica-Bold',17)
        pdf.drawString(left,752,'Fertility and housing: overnight memo')
        pdf.setFont('Helvetica',9);pdf.setFillColor(muted)
        pdf.drawString(left,735,subtitle)
        pdf.setStrokeColor(teal);pdf.line(left,724,left+content_width,724)
        pdf.setFont('Helvetica',7)
        pdf.drawString(left,25,'As of '+summary['generated_utc'][:19].replace('T',' ')+' UTC | provisional exploration')
        pdf.drawRightString(left+content_width,25,f'{page} / 2')
        return 710

    def checked_bottom(y):
        if y < 43:
            raise ValueError('Memo commentary exceeds two pages; shorten review notes')

    best = summary['selected']; counts=summary['counts']
    y=header(1,'Fit, progress and the economic misses')
    lead = ('Best completed primary-weight point: loss <b>'+fmt(best['primary_loss'])+'</b>. '
            'This is a provisional searched point, not a certified final calibration.') if best else (
            '<b>No completed primary-weight point is available.</b> Targets and restrictions below are known; '
            'model values remain blank. This memo does not invent results.')
    y=paragraph(lead,y,max_height=43)
    if summary.get('lead_note'):
        y=paragraph(escape(summary['lead_note']),y,small,max_height=35)
    status = ', '.join(f'{k.replace("_", " ")}: {v}' for k,v in sorted(counts.items())) or 'No completed candidate records yet'
    allocation=summary['expected_workers']
    allocation_text=(' Requested capacity: '+', '.join(f'{v} {k}' for k,v in allocation.items())+'.') if allocation else ''
    y=paragraph('<b>Progress:</b> '+escape(status)+'.'+escape(allocation_text)+' '+escape(summary['status_note'][:330]),y,small,max_height=55)
    group_labels={'primary':'Main weights','identity':'Unit weights',
                  'early_fertility_3000':'Early fertility weight 3,000'}
    group_rows=[['Weight system','Search / checks','Loss under main weights']]
    for name, group in summary['groups'].items():
        group_rows.append([group_labels[name],f"{group['search_completed']} / {group['verification_or_other']}",fmt(group['own_winner_primary_loss'])])
    if len(group_rows)>5:
        raise ValueError('At most four weighting systems fit the two-page memo; consolidate the index')
    y=table(group_rows,[225,110,201],y)
    y=paragraph('Checks include acceptance and final-repeat evaluations. Each row scores the point chosen by that weight system using the main weights. '
                'Full re-ranking is saved separately. The target table below reports the main-weight winner.',y,small,max_height=31)
    y=heading('All 14 restrictions: 13 scored moments and one normalization',y)
    data=[['Moment','Target','Model','Gap','Weight','Loss']]
    for row in summary['target_fit']:
        data.append([row['label'],fmt(row['target']),fmt(row['model']),fmt(row['gap']),fmt(row['weight']),fmt(row['loss_contribution'])])
    y=table(data,[225,55,55,55,68,78],y)
    y=paragraph('Gap = model - target. Weight and loss use the primary objective. Normalization is unscored. '
                'Full-precision source tables are saved beside this memo.',y,small,max_height=23)
    if best:
        misses=[]
        for row in summary['largest_weighted_misses']:
            share=100*row['loss_contribution']/best['primary_loss'] if best['primary_loss'] else 0
            misses.append(f'{escape(row["label"])}: {"above" if row["gap"]>0 else "below" if row["gap"]<0 else "at"} target by {fmt(abs(row["gap"]))} ({share:.1f}% of loss)')
        y=paragraph('<b>Largest weighted misses:</b> '+'; '.join(misses)+'.',y,small,max_height=46)
    else:
        y=paragraph('<b>Economics:</b> no complete fit yet supports a diagnosis of which mechanisms are missing.',y,small,max_height=23)
    checked_bottom(y);pdf.showPage()

    y=header(2,'Parameters, reliability and the next decisions')
    y=paragraph('Ten fitted parameters: nine searched jointly and one benefit level solved internally to match completed fertility at 2.1. '
                'Four additional restrictions do not by themselves establish identification.',y,max_height=34)
    data=[['Parameter','Estimate','Lower','Upper','Bound / restriction']]
    for row in summary['parameters']:
        flag=('near bound' if row['near_bound'] else 'interior') if row['estimate'] is not None and row['near_bound'] is not None else '--'
        if row['parameter']=='psi_child':flag='positive; normalized'
        data.append([row['label'],fmt(row['estimate']),fmt(row['lower']),fmt(row['upper']),flag])
    y=table(data,[233,66,62,62,113],y)
    y=paragraph('Near bound means within 1% of the permitted interval. The benefit level has a positivity restriction and is fixed by the fertility equation; '
                'it has no separately adopted upper search bound.',y,small,max_height=31)
    y=heading('Reliability, speed and errors',y)
    timings=summary['median_stationary_solve_seconds']
    timing_text=('Median accumulated stationary-solve time per completed objective: '+fmt(timings/60)+' minutes.') if timings is not None else 'Full-objective timing is not available.'
    diagnostic_text=(f'Selected point has {best["diagnostic_pngs"]}/17 standard diagnostic PNGs.') if best else 'No selected-point diagnostic packet is available.'
    repeats=summary['repeats'].get('status','Exact repeated-result acceptance is not recorded in this report index.')
    benchmark=summary['benchmark'].get('status','No matched warm/cold speed claim is established here.')
    y=paragraph(escape(timing_text+' '+diagnostic_text+' '+str(repeats)+' '+str(benchmark)),y,small,max_height=54)
    errors=summary['errors'];stale=sum(bool(s.get('stale_30_minutes')) for s in summary['controller_states'])
    error_text=f'{len(errors)} recorded error/exclusion entries; {stale} controller heartbeats older than 30 minutes.'
    if errors:
        error_text+=' Latest: '+str(errors[-1]['error'])[:240]
    y=paragraph(escape(error_text),y,small,max_height=44)
    y=heading('What remains uncertain',y)
    notes=summary['review_notes'] or [
        'The first-birth housing observer and bequest/older-wealth data counterparts retain documented approximations.',
        'Interest/credit, tenure scale, rental support, income/entry compatibility and weight choices remain substantive review items; no overnight economic change is implied.',
        'The finer housing grid and conception schedule remain separate tests.']
    if len(notes)>3:
        raise ValueError('Use at most three concise review_notes; full details belong in supporting records')
    for note in notes:
        y=paragraph('&#8226; '+escape(str(note)),y,small,max_height=35)
    y=heading('Next steps',y)
    steps=summary['next_steps'] or [
        'Complete or verify the selected-point repeats, numerical gates and the stable diagnostic packet before treating the fit as a baseline.',
        'Use the largest economic misses and parameters near bounds to choose the next targeted refinement; keep weighting experiments separate.']
    if len(steps)>3:
        raise ValueError('Use at most three concise next_steps')
    for i,step in enumerate(steps,1):
        y=paragraph(f'{i}. '+escape(str(step)),y,small,max_height=35)
    y=paragraph('Audit trail: summary.json, monitor.json and full-precision CSVs accompany this memo. '
                'The supplemental weighted-gap SVG is separate from the unchanged 17 standard diagnostics. '
                'This report checks saved tables and receipts; it does not independently re-solve or re-certify checkpoint arrays.',y,small,max_height=34)
    checked_bottom(y);pdf.showPage();pdf.save()


def build(index, output, pdf=True):
    summary=collect(index);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    monitor={k:summary[k] for k in ('generated_utc','status','counts','groups','controller_states','errors','median_stationary_solve_seconds')}
    monitor['selected_path']=summary['selected']['path'] if summary['selected'] else None
    monitor['largest_weighted_misses']=summary['largest_weighted_misses']
    (output/'monitor.json').write_text(json.dumps(monitor,indent=2,allow_nan=False)+'\n')
    write_csv(output/'target_fit.csv',summary['target_fit'],['moment','label','target','model','gap','weight','loss_contribution'])
    write_csv(output/'parameters.csv',summary['parameters'],['parameter','label','estimate','lower','upper','near_bound','restriction'])
    write_csv(output/'candidates.csv',summary['candidates'],['path','weighting','phase','own_loss','primary_loss','stationary_solves','solve_seconds','diagnostic_pngs','source_manifest_sha256','receipt_sha256','status'])
    write_overview(output/'fit_overview.svg',summary)
    if pdf:render_pdf(output/'memo.pdf',summary)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--json-only',action='store_true',help='Write the complete evidence bundle without the PDF')
    args=parser.parse_args();result=build(args.index,args.output,pdf=not args.json_only)
    print(json.dumps({'status':result['status'],'counts':result['counts'],'groups':result['groups'],'output':str(args.output.resolve())}))


if __name__=='__main__':
    main()
