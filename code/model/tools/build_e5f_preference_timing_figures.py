"""Three supplemental timing figures and a complete review PDF; no model solve.

Run with --pilot-root, --output-pdf and --scratch. Uses verified collected results.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import zipfile

# Optional pure-Python PDF libraries, without changing the model environment.
if os.environ.get('E5F_PDF_SITE_PACKAGES'):
    sys.path.append(os.environ['E5F_PDF_SITE_PACKAGES'])
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.platypus import Table, TableStyle

COLORS = {0.: '#24578C', -.5: '#008477', .5: '#C67728'}
LABELS = {0.: 'Linear baseline', -.5: 'Later decline', .5: 'Earlier decline'}
ORDER = (0., -.5, .5)
W, H = 960, 540


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#AAAEB3')
    ax.tick_params(length=0, pad=9, labelsize=12)
    ax.grid(axis='y', color='#E7E9EC', linewidth=.8)
    ax.set_axisbelow(True)


def frame(title, subtitle, page):
    fig = plt.figure(figsize=(W/72,H/72), facecolor='white')
    fig.text(.065,.923,title,fontsize=23,color='#192C40')
    fig.text(.065,.871,subtitle,fontsize=13,color='#4C5663')
    fig.text(.945,.028,str(page),ha='right',fontsize=9,color='#737B84')
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot-root',type=Path,required=True)
    p.add_argument('--output-pdf',type=Path,required=True)
    p.add_argument('--scratch',type=Path,required=True)
    args = p.parse_args()
    root, work = args.pilot_root, args.scratch
    figures = root/'figures'
    figures.mkdir(exist_ok=True); work.mkdir(parents=True,exist_ok=True)
    args.output_pdf.parent.mkdir(parents=True,exist_ok=True)
    inputs = [root/'computation'/name for name in ('all_main_birth_comparisons.csv',
        'all_main_target_fits.csv','all_inherited_parameters.csv','verified_receipts.json')]
    births, fits, parameters = [read_csv(path) for path in inputs[:3]]
    receipt = json.loads(inputs[3].read_text())
    assert not receipt['pending'] and len(receipt['verified']) == 6
    cases = {r['shape']:r for r in receipt['verified'] if r['phase']=='main'}
    assert set(cases)==set(ORDER) and len(births)==12 and len(fits)==36 and len(parameters)==15
    by_shape = {s:[r for r in births if float(r['shape'])==s] for s in ORDER}
    fit_by_shape = {s:[r for r in fits if float(r['shape'])==s] for s in ORDER}
    for s in ORDER:
        assert [int(r['decision_year']) for r in by_shape[s]] == [2007,2011,2015,2019]
        for r in fit_by_shape[s]:
            gap = float(r['model'])-float(r['target'])
            assert math.isclose(gap,float(r['gap']),rel_tol=0,abs_tol=1e-12)
            assert math.isclose(gap*gap*float(r['weight']),float(r['loss_contribution']),rel_tol=1e-12)
        assert math.isclose(sum(float(r['loss_contribution']) for r in fit_by_shape[s]),
                            cases[s]['inherited_objective'],rel_tol=1e-12)
    observed = np.array([float(r['observed_birth_index']) for r in by_shape[0.]])*100
    indices = {s:np.array([float(r['model_shape_index']) for r in by_shape[s]])*100 for s in ORDER}
    levels = {s:np.array([float(r['model_adjusted_births']) for r in by_shape[s]]) for s in ORDER}
    effects = {s:100*(levels[s]/levels[0.]-1) for s in ORDER}
    psi = {}
    for i,s in enumerate((-.5,0.,.5)):
        path = root/'computation/main'/f'case_{i}'/'contract.json'; inputs.append(path)
        psi[s] = np.array(json.loads(path.read_text())['psi_path'],dtype=float)
        assert psi[s][0]==psi[0.][0] if 0. in psi else True
        assert np.all(psi[s][4:]==psi[s][4])
    decline = {s:100*(psi[s][0]-psi[s][:7])/(psi[s][0]-psi[s][4]) for s in ORDER}
    plots, artist_checks = [], []
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.labelsize':13,
        'pdf.fonttype':42,'svg.fonttype':'none','font.size':12})

    def save(fig, name):
        path = work/(name+'.pdf')
        fig.savefig(path)
        fig.savefig(figures/(name+'.png'),dpi=180)
        fig.savefig(figures/(name+'.svg'))
        plots.append(path); plt.close(fig)

    blocks = ['2008-2011','2012-2015','2016-2019','2020-2023']
    x = np.arange(4)
    fig = frame('Aggregate births across four-year periods',
        'Observed birth counts and the model\'s linear preference decline',1)
    ax = fig.add_axes([.09,.255,.80,.535]); style(ax)
    for y,color,label,marker in [(observed,'#222222','Observed US births','o'),
            (indices[0.],COLORS[0.],'Model adjusted births','s')]:
        line, = ax.plot(x,y,color=color,label=label,marker=marker,lw=2.7,ms=7)
        assert np.array_equal(line.get_ydata(),y)
        artist_checks.append(dict(figure=1,series=label,y=y.tolist()))
    ax.set_xticks(x,blocks); ax.set_xlim(-.1,3.9); ax.set_ylim(75,103)
    ax.set_yticks([75,80,85,90,95,100])
    ax.set_ylabel('Birth count index (2008-2011 = 100)')
    ax.legend(loc='upper right',frameon=False,fontsize=12)
    ax.text(3.13,observed[-1],f'Data: {observed[-1]:.1f}\n({observed[-1]-100:.1f}%)',va='center',fontsize=13)
    ax.text(3.13,indices[0.][-1],f'Model: {indices[0.][-1]:.1f}\n({indices[0.][-1]-100:.1f}%)',
            va='center',fontsize=13,color=COLORS[0.])
    fig.text(.065,.145,'The 20.9% compares the last block with the first; it is not a cumulative shortfall.',fontsize=15)
    fig.text(.065,.089,r'$100\,(B_{2020-2023}/B_{2008-2011}-1)=-20.9\%$ for the model; '+
             r'$-11.0\%$ for the data.',fontsize=15)
    fig.text(.065,.028,'Sources: NCHS, Births: Final Data for 2015 and 2023; model calculations.',fontsize=9)
    save(fig,'birth_counts_by_period')

    fig = frame('Preference timing and the birth response',
        'The full path is announced in 2007; preference endpoints and baseline prices are held fixed.',2)
    left = fig.add_axes([.09,.285,.365,.46]); right = fig.add_axes([.59,.285,.355,.46])
    for ax in (left,right): style(ax)
    years = np.arange(2007,2032,4)
    for s in ORDER:
        line, = left.plot(years,decline[s],color=COLORS[s],lw=2.4,marker='o',ms=5,label=LABELS[s])
        assert np.array_equal(line.get_ydata(),decline[s])
        artist_checks.append(dict(figure=2,series=LABELS[s]+' preference decline',y=decline[s].tolist()))
    left.set_title('Timing of the preference decline',fontsize=15,pad=17)
    left.set_xlabel('Decision year'); left.set_ylabel('Decline completed (%)')
    left.set_xticks([2007,2015,2023,2031]); left.set_ylim(-5,110)
    left.set_yticks([0,25,50,75,100]); left.axvline(2023,color='#888888',ls=':',lw=1)
    left.legend(frameon=False,fontsize=10,loc='lower right')
    right.axhline(0,color=COLORS[0.],lw=1.4)
    for s,marker in [(-.5,'o'),(.5,'s')]:
        line, = right.plot(x,effects[s],color=COLORS[s],lw=2.4,marker=marker,ms=6,label=LABELS[s])
        assert np.array_equal(line.get_ydata(),effects[s])
        artist_checks.append(dict(figure=2,series=LABELS[s]+' birth effect',y=effects[s].tolist()))
    right.set_title('Births relative to the linear baseline',fontsize=15,pad=17)
    right.set_xticks(x,blocks,rotation=20,ha='right'); right.set_ylim(-2.8,2.8)
    right.set_yticks([-2,-1,0,1,2]); right.set_ylabel('Difference in birth counts (%)')
    fig.text(.065,.145,'Timing shifts births by roughly 1-2% relative to the baseline, including the first block.',fontsize=14)
    fig.text(.065,.093,'All three paths leave a first-to-last-block birth decline of about 21%, versus 11% in the data.',fontsize=14)
    fig.text(.065,.028,'Source: model calculations. Household responses evaluated at the same inherited price path.',fontsize=9)
    save(fig,'preference_timing_and_births')

    fig = frame('Conditional fit and housing-market clearing',
        'Existing targets and weights; all structural parameters and baseline prices held fixed.',3)
    left = fig.add_axes([.19,.295,.30,.445]); right = fig.add_axes([.66,.295,.285,.445])
    y = np.arange(3)
    for ax in (left,right):
        style(ax); ax.grid(False); ax.set_yticks(y,[LABELS[s] for s in ORDER]); ax.invert_yaxis()
    loss = [cases[s]['inherited_objective'] for s in ORDER]
    gaps = [100*cases[s]['maximum_market_residual'] for s in ORDER]
    left.barh(y,loss,color=[COLORS[s] for s in ORDER],height=.48)
    left.set_xlim(0,125); left.set_xlabel('Existing 12-moment objective')
    left.set_title('Target fit',fontsize=15,pad=20)
    right.barh(y,gaps,color=[COLORS[s] for s in ORDER],height=.48)
    right.set_xlim(0,.062); right.set_xticks([0,.02,.04,.06],['0.00','0.02','0.04','0.06'])
    right.set_xlabel('Maximum market discrepancy (%)')
    right.set_title('Market clearing',fontsize=15,pad=20)
    right.axvline(.02,color='#444444',lw=1.5,ls='--')
    right.text(.021,.98,'Tolerance',fontsize=10,color='#444444',
               transform=right.get_xaxis_transform(),va='top')
    for i,(a,b) in enumerate(zip(loss,gaps)):
        left.text(a+2,i,f'{a:.1f}',va='center',fontsize=13)
        right.text(b+.0015,i,f'{b:.4f}',va='center',fontsize=11)
    improvement=100*(1-cases[-.5]['inherited_objective']/cases[0.]['inherited_objective'])
    fig.text(.065,.148,f'Later timing lowers the objective by {improvement:.1f}%; both alternatives fail market clearing.',fontsize=14)
    fig.text(.065,.096,'The birth-path gap improves only slightly; a full equilibrium comparison requires new prices.',fontsize=14)
    fig.text(.065,.028,'Source: model calculations using the unchanged 12-moment target system.',fontsize=9)
    artist_checks.append(dict(figure=3,loss=loss,market_discrepancy_percent=gaps,market_tolerance_percent=.02))
    save(fig,'conditional_fit_and_market_clearing')

    labels = ['Completed fertility (legacy tfr row)','Childlessness','Mean age at first birth (years)',
        'Share of first births at age 30+','Rooms added after first birth',
        'Parent room gap: 3+ versus 1-2 children','Parent ownership gap','Ownership rate',
        'Mean occupied rooms, ages 18-85','Wealth / annual gross labor earnings',
        'Annual bequest flow / aggregate wealth','Old wealth/income: p90 / median']
    appendix=work/'appendix.pdf'; pdf=canvas.Canvas(str(appendix),pagesize=(W,H))
    formatted=[]
    def head(title,subtitle,page):
        pdf.setFillColor(colors.HexColor('#192C40')); pdf.setFont('Helvetica',23)
        pdf.drawString(62,494,title); pdf.setFillColor(colors.HexColor('#4C5663'))
        pdf.setFont('Helvetica',12); pdf.drawString(62,465,subtitle)
        pdf.setFont('Helvetica',9); pdf.drawRightString(906,18,str(page))
    def table(data,widths,fontsize=11):
        t=Table(data,colWidths=widths,rowHeights=23)
        t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
            ('FONTNAME',(0,1),(-1,-1),'Helvetica'),('FONTSIZE',(0,0),(-1,-1),fontsize),
            ('TEXTCOLOR',(0,0),(-1,-1),colors.HexColor('#263442')),
            ('LINEABOVE',(0,0),(-1,0),1,colors.HexColor('#263442')),
            ('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#90969D')),
            ('LINEBELOW',(0,-1),(-1,-1),1,colors.HexColor('#263442')),
            ('ALIGN',(1,0),(-1,-1),'RIGHT'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('LEFTPADDING',(0,0),(-1,-1),3),('RIGHTPADDING',(0,0),(-1,-1),4)]))
        _,height=t.wrap(W-124,400); t.drawOn(pdf,62,430-height)
    for page,s in enumerate(ORDER,start=4):
        head('Full target fit: '+LABELS[s].lower(),
            'Existing target system; parameters held fixed. All 12 rows are retained.',page)
        data=[['Moment','Target','Model','Gap','Weight','Contribution']]
        for label,row in zip(labels,fit_by_shape[s]):
            values=[f"{float(row[k]):.6f}" for k in ('target','model','gap')]
            values += [f"{float(row['weight']):.6g}",f"{float(row['loss_contribution']):.6f}"]
            data.append([label]+values); formatted.extend(values)
        table(data,[370,92,92,92,92,98],10.5)
        pdf.setFillColor(colors.HexColor('#263442'));pdf.setFont('Helvetica',12)
        pdf.drawString(62,105,f"Total objective: {cases[s]['inherited_objective']:.6f}.  Each contribution equals weight x gap squared.")
        pdf.setFont('Helvetica',10)
        pdf.drawString(62,72,'Shares are fractions; housing outcomes are rooms. Measurement definitions remain those of the inherited target system.')
        pdf.drawString(62,53,'These are conditional diagnostics. Terminal unit-rent distance is 1.088827% against a 1% tolerance for all three cases.')
        pdf.drawString(62,34,'Birth-count figures include imputed additional births at entry into the 3+ child group; they do not measure female TFR.')
        pdf.showPage()
    head('Inherited parameters and restrictions','All 11 previously estimated coordinates are held fixed in the timing comparison.',7)
    data=[['Parameter','Value','Lower bound','Upper bound','Role','Near bound']]
    for row in parameters:
        values=[f"{float(row['value']):.6f}"]
        values += [f"{float(row[k]):.6g}" if row[k] else '-' for k in ('lower_bound','upper_bound')]
        role='Held fixed' if row['is_free_parameter']=='True' else (
            'Normalized' if row['parameter']=='psi_child_2007' else 'Derived' if row['parameter']=='psi_child_2023' else 'External')
        data.append([row['parameter']]+values+[role,'Yes' if row['near_bound']=='True' else 'No'])
        formatted.extend([v for v in values if v!='-'])
    table(data,[305,110,110,110,120,81],11)
    pdf.setFont('Helvetica',10);pdf.setFillColor(colors.HexColor('#263442'))
    pdf.drawString(62,42,'Bounds and near-bound flags reproduce the inherited search domain. No parameter estimate changed in this exercise.')
    pdf.showPage();pdf.save()
    writer=PdfWriter()
    for path in plots+[appendix]:
        writer.append(str(path))
    with args.output_pdf.open('wb') as stream:writer.write(stream)
    reader=PdfReader(args.output_pdf)
    assert len(reader.pages)==7
    appendix_text='\n'.join(page.extract_text() for page in reader.pages[3:])
    assert all(value in appendix_text for value in formatted)
    qa=dict(status='numeric_and_text_checks_passed_visual_review_pending',pages=7,
        input_sha256={str(path):sha(path) for path in inputs},
        pdf_sha256=sha(args.output_pdf),artist_data_checks=artist_checks,
        target_numeric_cells_verified=36*5,parameter_numeric_cells_verified=sum(
            1+bool(r['lower_bound'])+bool(r['upper_bound']) for r in parameters),
        endpoint_decline_percent={str(s):100-indices[s][-1] for s in ORDER},
        observed_endpoint_decline_percent=100-observed[-1],
        meaning='2020-2023 total versus 2008-2011 total; not cumulative and not single-year 2023',
        no_model_solve=True,standard_diagnostic_set_unchanged=True)
    (figures/'verification.json').write_text(json.dumps(qa,indent=2)+'\n')
    with zipfile.ZipFile(figures/'preference_timing_plots.zip','w',zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(figures.iterdir()):
            if path.suffix in ('.png','.svg'):bundle.write(path,path.name)
    print(json.dumps(dict(pdf=str(args.output_pdf),pages=7,figures=str(figures),sha256=qa['pdf_sha256']),indent=2))


if __name__=='__main__':
    main()
