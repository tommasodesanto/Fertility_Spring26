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
import math
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
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle

ROOT = Path(__file__).resolve().parents[3]
DEFAULT = ROOT / 'output/model/e5f_matched_pf_20260909a'


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assessment(base, out, path_root, late, export):
    """Source-checked tables and calendar-aligned diagnostics, without model calls."""
    source_pins = {}
    def checked_csv(p):
        source_pins[str(p)] = digest(p)
        return read_csv(p)
    late_fits = {c: [r for r in checked_csv(late/c/'target_fit.csv') if r['repetition']=='1']
                for c in ['old_balanced','new_balanced']}
    for c, rows in late_fits.items():
        for r in rows:
            gap = float(r['model'])-float(r['target'])
            if not math.isclose(gap,float(r['gap']),abs_tol=1e-12): raise ValueError('Fit gap mismatch')
            if not math.isclose(gap*gap*float(r['weight']),float(r['loss_contribution']),rel_tol=1e-11,abs_tol=1e-11):
                raise ValueError('Fit contribution mismatch')
        receipt=json.loads((late/c/'comparison_receipt.json').read_text())
        if not math.isclose(sum(float(r['loss_contribution']) for r in rows),receipt['loss'],rel_tol=1e-12):
            raise ValueError('Complete late objective mismatch')
    early_root=base/'initial_calibration_contract/saved_case_scores'
    early=checked_csv(early_root/'joint_12_analysis_28_fit.csv')
    early_params=checked_csv(early_root/'joint_12_analysis_28_parameters.csv')
    for r in early:
        gap=float(r['model'])-float(r['target'])
        if not math.isclose(gap,float(r['gap']),abs_tol=1e-12): raise ValueError('Early gap mismatch')
        if r['scored']=='True' and not math.isclose(gap*gap*float(r['actual_weight']),float(r['loss_contribution']),rel_tol=1e-11,abs_tol=1e-11):
            raise ValueError('Early contribution mismatch')
    early_receipt=json.loads((base/'initial_calibration_contract/exact_loop_smoke/completed_17370427/summary.json').read_text())
    if not math.isclose(sum(float(r['loss_contribution']) for r in early if r['scored']=='True'),early_receipt['loss'],rel_tol=1e-12):
        raise ValueError('Early complete objective disagrees with exact replay')
    labels = {'tfr':'Completed-fertility measure', 'childless_rate':'Childlessness',
              'mean_age_first_birth':'Mean age at first birth', 'share_first_births_age30plus':'First births at 30+',
              'housing_increment_0to1':'First-birth room response',
              'prime30_55_parent_3plus_minus_1to2_mean_rooms':'Rooms: 3+ vs 1-2 children',
              'own_family_gap':'Parent ownership gap (old definition)', 'own_rate':'Ownership (old definition)',
              'aggregate_mean_occupied_rooms_18_85':'Mean occupied rooms',
              'aggregate_wealth_to_annual_gross_labor_earnings':'Wealth / annual labor earnings',
              'annual_bequest_flow_to_aggregate_wealth':'Annual bequests / wealth',
              'old_total_wealth_to_annual_income_p90_p50_7684':'Old wealth/income: p90 / median'}
    def number(v): return '--' if v in ('',None) else f'{float(v):.5g}'
    late_table=[['Moment','Target','Old','New','Old gap','New gap','Weight','Old loss','New loss']]
    for a,b in zip(late_fits['old_balanced'],late_fits['new_balanced']):
        if any(a[k]!=b[k] for k in ['moment','target','weight']): raise ValueError('Mixed late targets')
        late_table.append([labels[a['moment']]]+[number(v) for v in [a['target'],a['model'],b['model'],a['gap'],b['gap'],a['weight'],a['loss_contribution'],b['loss_contribution']]])
    early_table=[['Moment','Target','Model','Gap','Weight','Loss']]
    for r in early:
        early_table.append([r['label'].replace('\u2013','-').replace('\u2212','-')]+[number(r[k]) for k in ['target','model','gap','actual_weight','loss_contribution']])
    late_params={c:{r['parameter']:r for r in checked_csv(late/c/'repetition_02/parameters.csv')} for c in late_fits}
    ptable=[['Parameter','Static old','Static new','Early candidate','Range / restriction','Near bound?']]
    for r in early_params:
        name=r['parameter'];a=late_params['old_balanced'][name];b=late_params['new_balanced'][name]
        bounds=f"[{number(r['lower'])}, {number(r['upper'])}]" if r['lower'] else r['status']
        if name=='psi_child':bounds='Fixed in static tests; early normalized to 2.1'
        if name=='hbar_child_rooms':bounds='Old [0.10,1.80]; new fixed at zero'
        if name=='h_P':bounds='New [0.10,2.30]; old = jump + slope'
        if name=='pension_period':bounds='Derived from payroll budget'
        ptable.append([name,number(a['estimate']),number(b['estimate']),number(r['estimate']),bounds,
                       'Yes: theta1, all three' if name=='theta1' else 'No / fixed'])
    a=late_params['old_balanced']
    old_jump=float(a['h_P']['estimate'])-float(a['hbar_child_rooms']['estimate'])
    ptable.append(['Old first-child jump',number(old_jump),'--','--','Old [0,0.5]; new uses h_P','No'])

    data=base/'path_pilot_20260910/fertility_data'
    declared=json.loads((data/'artifact_sha256.json').read_text())
    for name in ['empirical_blocks.csv','annual_household_stocks_2007_2023.csv']:
        if digest(data/name)!=declared[name]: raise ValueError('Empirical source changed')
    blocks=checked_csv(data/'empirical_blocks.csv')
    path={int(r['calendar_year']):r for r in checked_csv(path_root/'transition_path.csv')}
    first=float(path[2007]['birth_children_topcode_adjusted'])
    births=[]
    for d in blocks:
        year=int(d['decision_year']);r=path[year]
        if int(d['birth_year_start'])!=year+1 or int(d['birth_year_end'])!=year+4:
            raise ValueError('Birth-block clock changed')
        births.append(dict(decision_year=year,birth_year_start=year+1,birth_year_end=year+4,
            data_births=int(d['live_births_total']),data_index=100*float(d['live_births_index_2008_2011']),
            model_births=float(r['birth_children_topcode_adjusted']),
            model_index=100*float(r['birth_children_topcode_adjusted'])/first))
    hh={int(r['year']):float(r['households_stock_proxy']) for r in checked_csv(data/'annual_household_stocks_2007_2023.csv')}
    yy=[2007,2011,2015,2019,2023]
    model_hh=[float(path[y]['adult_population'] or path[y]['household_heads'])/float(path[2007]['adult_population'])*100 for y in yy]
    data_hh=[hh[y]/hh[2007]*100 for y in yy]
    fig,axes=plt.subplots(1,2,figsize=(13.8,6.7),constrained_layout=True)
    x=list(range(4))
    for key,label,color in [('data_index','NCHS observed births','#222222'),('model_index','Current short transition','#176b94')]:
        vals=[r[key] for r in births];line,=axes[0].plot(x,vals,'o-',label=label,color=color,lw=2)
        if list(line.get_ydata())!=vals: raise ValueError('Birth artist changed')
        export.append(dict(figure='birth_data',series=key,x=x,y=vals))
    axes[0].set(title='Birth counts: model versus data',ylabel='Index: 2008-2011 = 100',
                xticks=x,xticklabels=['2008-11','2012-15','2016-19','2020-23']);axes[0].legend()
    for vals,label,color,style in [(data_hh,'Census household stocks','#222222','o-'),(model_hh,'Model (historically conditioned)','#d77a28','x--')]:
        axes[1].plot(yy,vals,style,label=label,color=color,lw=2)
        export.append(dict(figure='household_data',series=label,x=yy,y=vals))
    axes[1].set(title='Household counts: imposed historical information',ylabel='Index: 2007 = 100',xticks=yy)
    axes[1].legend(fontsize=9)
    fig.suptitle('What the short transition matches - and what it does not',fontsize=17,y=.98)
    fig.get_layout_engine().set(rect=(0,.15,1,.74))
    fig.text(.025,.035,'Birth blocks follow the established decision-year + 1 to + 4 clock. Each birth series is normalized to its first block.\n'
             'This compares shapes, not levels or female TFR; national data versus model geography remains an approximation.\n'
             'Household counts enter the historical conditioning and are not independent validation. Sources: NCHS final reports; Census HH-3.',fontsize=10)
    birth_figure=out/'transition_vs_data.png';fig.savefig(birth_figure);plt.close(fig)

    reduced=out/'current_policy_extraction/results'
    receipt=json.loads((reduced/'receipt.json').read_text())
    for name,pin in receipt['output_sha256'].items():
        if digest(reduced/name)!=pin:raise ValueError('Policy extraction changed')
    policy=checked_csv(reduced/'current_transition_conditional_policy.csv')
    pools={}
    for r in policy:pools.setdefault((int(r['age']),int(r['z_index'])),[]).append(r)
    def score(rows):
        return sum(float(a['post_fertility_pre_tenure_mass'])*max(0,float(a['owner_entry_probability'])-float(b['owner_entry_probability']))
                   for a,b in zip(rows,rows[1:]) if a['graph_value_mask']=='True' and b['graph_value_mask']=='True')
    key=max(pools,key=lambda k:score(pools[k]));rows=pools[key]
    products=receipt['cases'][0]['owner_products'];wealth=[float(r['wealth']) for r in rows]
    valid=[r['graph_value_mask']=='True' for r in rows]
    mass=[float(r['post_fertility_pre_tenure_mass']) for r in rows];total=sum(mass)
    expected=[]
    for r,ok in zip(rows,valid):
        value=float(r['tenure_probability_0'])*float(r['housing_conditional_renter'])
        value+=sum(float(r[f'tenure_probability_{i+1}'])*h for i,h in enumerate(products))
        if ok and not math.isclose(float(r['tenure_probability_sum']),1,abs_tol=2e-6):raise ValueError('Tenure probability normalization')
        expected.append(value if ok else float('nan'))
    def masked(name):return [float(r[name]) if ok else float('nan') for r,ok in zip(rows,valid)]
    owner=masked('owner_entry_probability')
    occupied_dips=[i for i in range(len(rows)-1) if valid[i] and valid[i+1] and mass[i]>1e-12 and owner[i]-owner[i+1]>.01]
    fig,axes=plt.subplots(2,2,figsize=(13.8,7.8),constrained_layout=True)
    axes[0,0].plot(wealth,masked('consumption_conditional_renter'),color='#176b94');axes[0,0].set(title='Consumption conditional on renting',ylabel='Model consumption units')
    for vals,label,color,style in [(masked('housing_graph_selected_tenure'),'Selected tenure/product','#888888','--'),
                                  (masked('housing_conditional_renter'),'Conditional rental housing','#176b94',':'),
                                  (expected,'Probability-weighted housing','#d77a28','-')]:
        axes[0,1].plot(wealth,vals,label=label,color=color,ls=style,lw=2)
    axes[0,1].set(title='Three different housing objects',ylabel='Housing units');axes[0,1].legend(fontsize=9)
    axes[1,0].plot(wealth,owner,color='#176b94')
    axes[1,0].scatter([wealth[i] for i in occupied_dips],[owner[i] for i in occupied_dips],color='#b82828',s=25,label='Occupied node before >1pp decline')
    axes[1,0].set(title='Ownership probability: reversals remain',ylabel='Probability');axes[1,0].legend(fontsize=9)
    axes[1,1].plot(wealth,[100*m/total for m in mass],'o-',color='#222222',ms=3)
    axes[1,1].set(title='Where households in this slice actually are',ylabel='Percent of slice mass per wealth node')
    for ax in axes.flat:ax.set_xlabel('Liquid wealth')
    fig.suptitle(f'2023 policy diagnosis: age {key[0]}, income state {key[1]} (z={float(rows[0]["z"]):.3f})',fontsize=17,y=.98)
    fig.get_layout_engine().set(rect=(0,.105,1,.79))
    fig.text(.025,.025,'Slice chosen by the largest mass-weighted ownership decline across the inspected age-30/42 income slices.\n'
             'No smoothing or new solve. Probability-weighted housing corrects the object being viewed; it does not certify the underlying numerical solution.',fontsize=10)
    policy_figure=out/'policy_objects_and_occupancy.png';fig.savefig(policy_figure);plt.close(fig)
    finite=lambda values:[v if math.isfinite(v) else None for v in values]
    export.append(dict(figure='policy_objects',series='expected_housing',x=wealth,y=finite(expected)))
    export.append(dict(figure='policy_objects',series='owner_probability',x=wealth,y=finite(owner)))
    all_fits={'early_candidate':early, **late_fits}
    fit_screen=[]
    for name, fit_rows in all_fits.items():
        for r in fit_rows:
            weight=r.get('actual_weight',r.get('weight'))
            if weight:
                scaled_gap=float(r['gap'])*math.sqrt(float(weight))
                fit_screen.append(dict(case=name,moment=r.get('label',r.get('moment')),
                    scaled_gap=scaled_gap,review_flag=abs(scaled_gap)>2))
    packet=dict(late_target_table=late_table,early_target_table=early_table,parameter_table=ptable,
                birth_comparison=births,selected_policy_slice=list(key),selection_score=score(rows),
                inspected_slice_mass=total,occupied_large_dips=len(occupied_dips),source_sha256=source_pins,
                fit_screen=fit_screen,fit_screen_definition='Flag abs(gap)*sqrt(weight)>2 for review only; objective scales include synthetic scales, not confidence intervals. No target or weight is changed.',
                losses={c:sum(float(r['loss_contribution']) for r in rows if r['loss_contribution']) for c,rows in all_fits.items()})
    (out/'assessment.json').write_text(json.dumps(packet,indent=2))
    return packet,birth_figure,policy_figure


def review_tables(packet, summary):
    """The same explicit questions appear in the PDF, HTML and machine-readable receipt."""
    checks=[['Check','Success means','Current result'],
        ['Pension accounts','Balance at every solved date under the pinned fiscal gate.',
         ('PASS' if summary['social_security_gate'] else 'FAIL')+f"; max scaled residual {summary['maximum_fiscal_residual']:.3g}"],
        ['Housing equilibrium','Demand equals supply under the pinned market gate.',
         ('PASS' if summary['market_gate'] else 'FAIL')+f"; max relative residual {summary['maximum_market_residual']:.3g}"],
        ['Reproducibility','Saved mapping replays and checkpoint reloads correctly.',
         'PASS' if summary['mapping_replay_verified'] and summary['checkpoint_reload_verified'] else 'FAIL'],
        ['Calibration fit','Every target and parameter visible; economically material misses resolved or explicitly assessed.',
         'NEEDS WORK; complete tables follow. A small total loss alone is insufficient.'],
        ['2023 utility comparison','Old and new each recalibrated to the same 2023 targets, fiscal closure and numerical gates.',
         'MISSING; current stationary tests freeze parameters and final preference.'],
        ['Policy functions','Feasible choices and probabilities; occupied reversals explained and sensitive states checked under numerical refinement.',
         'OPEN; plotting objects clarified. Global policy accuracy is not certified.'],
        ['Historical transition','Correctly timed, comparable data overlays; shock estimated; observed and imposed moments distinguished.',
         'NEEDS WORK; birth-count shape misses. Household path is conditioned on data.'],
        ['Terminal horizon','Approved tail checks pass and the dates of interest remain stable when the horizon is extended.',
         'NOT PASSED' if not summary['horizon_verified'] else 'PASS'],
        ['Policy readiness','All above resolved, with approved population, entry, fiscal and geographic closure.',
         'NOT READY' if not summary['production_eligible'] else 'PASS']]
    models=[['Experiment','What was varied / targeted','Interpretation'],
        ['Early stationary candidate','New utility; early-period working targets; fertility normalized to 2.1.',
         f"Loss {packet['losses']['early_candidate']:.3f}. Candidate, not certified SMM."],
        ['Static old / new, late preference','Same structural parameters, final preference and housing supply; balanced pensions.',
         f"Losses {packet['losses']['old_balanced']:.3f} / {packet['losses']['new_balanced']:.3f}. Neither recalibrated."],
        ['Six-date historical diagnostic','New utility; illustrative announced preference path, 2007-2027.',
         'Converged finite path. Neither fitted history nor horizon-certified policy benchmark.']]
    return checks, models


def draw_table_page(doc, size, title, note, rows, widths, footer, font=9):
    w,h=size
    doc.setFont('Helvetica-Bold',17);doc.drawString(26,h-30,title)
    style=ParagraphStyle('cell',fontName='Helvetica',fontSize=font,leading=font+2,textColor=colors.HexColor('#172b3a'))
    n=Paragraph(html.escape(note),ParagraphStyle('note',fontSize=10,leading=13))
    _,nh=n.wrap(w-52,h);n.drawOn(doc,26,h-46-nh)
    cells=[[Paragraph(html.escape(str(v)),style) for v in row] for row in rows]
    table=Table(cells,colWidths=widths,repeatRows=1)
    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dce9ef')),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f3f6f8')]),
        ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),
        ('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),5),
        ('BOTTOMPADDING',(0,0),(-1,-1),5)]))
    for j,heading in enumerate(rows[0]):
        if heading.lower().endswith('loss'):
            for i,row in enumerate(rows[1:],1):
                if row[j]!='--' and float(row[j])>4:
                    table.setStyle(TableStyle([('BACKGROUND',(j,i),(j,i),colors.HexColor('#ffe4d6'))]))
    _,th=table.wrap(w-52,h)
    bottom=h-60-nh-th
    if bottom<35:raise ValueError(f'Table overflows page: {title}; bottom={bottom}')
    table.drawOn(doc,26,bottom)
    doc.setFont('Helvetica',8);doc.drawString(26,12,footer);doc.showPage()


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
    fig.suptitle('New utility and balanced pensions: solved six-date transition', fontsize=17,y=.98)
    fig.get_layout_engine().set(rect=(0, .075, 1, .82))
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
    fig.suptitle('Direct utility comparison: both pension budgets balanced', fontsize=17,y=.98)
    fig.get_layout_engine().set(rect=(0, .075, 1, .82))
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
        fig.suptitle('Old and new utility at the retained 2023 preference', fontsize=17,y=.98)
        fig.get_layout_engine().set(rect=(0, .075, 1, .82))
        fig.text(.025, .018, 'Balanced pensions in both arms; same structural parameters and supply curve; two exact repetitions each.\n'
                 'Stationary counterparts, not recalibrated 2023 economies. They do not reproduce the historical 2023 population.', fontsize=10)
        late_figure = out / 'utility_comparison_2023_preference.png'
        fig.savefig(late_figure); plt.close(fig)

    packet,birth_figure,policy_figure=assessment(base,out,path_root,late,export)
    checks,models=review_tables(packet,summary)
    (out/'scorecard.json').write_text(json.dumps(dict(checks=checks,models=models,
        production_eligible=summary['production_eligible'],model_solves=0,
        fit_screen=packet['fit_screen'],fit_screen_definition=packet['fit_screen_definition']),indent=2))
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
    doc.setTitle('Model scorecard: utility, pensions and fit')
    table_pages=[
        ('The model tester: what counts as success',
         'Read this page first on every update. PASS refers only to the named check and experiment; unresolved checks cannot become green because a solver converged.',
         checks,[120,350,319],9),
        ('Which models are being assessed?',
         'Maintained sequential fertility architecture. Old utility has a per-child housing floor; new utility has an equivalence scale and a first-dependent housing floor. Both retain linear child reward. All cases below balance pensions.',
         models,[155,350,284],10),
        ('Early stationary candidate: every target',
         f"Working objective = {packet['losses']['early_candidate']:.3f}; 12 scored restrictions plus separate fertility normalization. Nine structural coordinates; the preference intercept normalizes fertility. Shaded losses exceed 4 (gap exceeds 2 objective scales): a review screen, not a confidence test; some scales are synthetic. Not certified SMM.",
         packet['early_target_table'],[325,91,91,91,95,96],9),
        ('Static old / new at the final preference: every target',
         f"Diagnostic losses = {packet['losses']['old_balanced']:.3f} / {packet['losses']['new_balanced']:.3f}. Fertility near 1.40 is the stationary outcome at the retained final preference, not a recalibrated 2023 fit. Shading flags gaps exceeding 2 objective scales, not statistical rejection. Targets differ from the early system; do not compare these losses to its loss.",
         packet['late_target_table'],[205]+[73]*8,8),
        ('Parameters: estimates, restrictions and bounds',
         'Static columns are inherited, fixed diagnostic inputs. Early estimates are provisional. Bound proximity follows the saved search flags. Pension amount is an endogenous budget outcome; annual beta is raised to the fourth power in four-year periods.',
         packet['parameter_table'],[145,90,90,90,264,110],8)]
    front_images=[(birth_figure,'Transition versus observed data','Birth-count shapes are compared on aligned four-year blocks; household counts are historically conditioned.'),
                  (policy_figure,'Policy objects and household occupancy','Current 2023 checkpoint: read-only extraction, zero new Bellman solves; numerical reversals remain under investigation.')]
    total_pages=len(table_pages)+len(front_images)+len(pages)
    for number,(title,note,rows,widths,font) in enumerate(table_pages,1):
        draw_table_page(doc,size,title,note,rows,widths,f'Model scorecard | {number}/{total_pages} | Saved inputs; zero model solves',font)
    for number, (p, title, note) in enumerate(front_images+pages, len(table_pages)+1):
        title=title.split('. ',1)[1] if title[:1].isdigit() else title
        w, h = size
        doc.setFont('Helvetica-Bold', 15); doc.drawString(26, h - 29, title)
        doc.setFont('Helvetica', 9); doc.drawString(26, h - 45, note)
        img = ImageReader(str(p)); iw, ih = img.getSize()
        scale = min((w - 52) / iw, (h - 83) / ih)
        doc.drawImage(img, (w - iw * scale) / 2, 24 + (h - 83 - ih * scale) / 2,
                      width=iw * scale, height=ih * scale)
        doc.setFont('Helvetica', 8); doc.drawString(26, 12, f'Saved-result review | {number}/{total_pages} | Full 17-graph sets in the accompanying gallery')
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
''' + f'<p><a href="{link(pdf)}">Open the {total_pages}-page PDF</a> · First page: the tester. Full target and parameter tables follow.</p>'
    doc_html += '<style>table{width:100%;border-collapse:collapse;background:white;font-size:13px;margin:18px 0}th,td{padding:9px;border:1px solid #dde3e8;text-align:left;vertical-align:top}th{background:#dce9ef}tr:nth-child(even){background:#f3f6f8}.scroll{overflow-x:auto}</style>'
    for title,note,rows,_,_ in table_pages:
        def html_cell(i,j,v):
            tag='th' if i==0 else 'td'
            flag=i>0 and rows[0][j].lower().endswith('loss') and v!='--' and float(v)>4
            style=' style="background:#ffe4d6"' if flag else ''
            return f'<{tag}{style}>{html.escape(str(v))}</{tag}>'
        table_html=''.join('<tr>'+''.join(html_cell(i,j,v) for j,v in enumerate(row))+'</tr>' for i,row in enumerate(rows))
        doc_html+=f'<h2>{html.escape(title)}</h2><p>{html.escape(note)}</p><div class="scroll"><table>{table_html}</table></div>'
    doc_html+=f'<h2>Transition versus data</h2><img src="{birth_figure.name}"><p>Housing prices, rents, ownership and housing-quantity historical overlays with aligned geography and samples remain outstanding.</p>'
    doc_html+=f'<h2>What the policy plots mean</h2><img src="{policy_figure.name}"><p>Conditional renter consumption and selected tenure housing are different conditional objects. Probability-weighted physical housing clarifies the graph. This does not resolve all ownership reversals. A prior independent audit verified one stationary reversal only; it is not certification of these transition policies.</p>'
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
    receipt = {'source_graphs_verified': len(checked), 'pdf_pages': total_pages,
               'model_solves': 0, 'elapsed_seconds': time.monotonic() - started,
               'source_graph_sha256': checked,
               'source_tables_sha256': {str(p): digest(p) for p in [path_root / 'transition_path.csv',
                    decomp / 'collected/old_balanced/repetition_01/lifecycle_2023.csv',
                    decomp / 'collected/new_balanced/repetition_01/lifecycle_2023.csv']},
               'assessment_source_sha256':packet['source_sha256'],
               'pdf_sha256': digest(pdf), 'pdf': str(pdf)}
    (out / 'build_receipt.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_graph_sha256','source_tables_sha256','assessment_source_sha256']}))
    print('\n'.join(f'{r[0]}: {r[2]}' for r in checks[1:]))


if __name__ == '__main__':
    main()
