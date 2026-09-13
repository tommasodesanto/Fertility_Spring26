"""Two historical fertility panels motivating the initial-state approximation."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'


def project_cps_age_profile(source,period_level,anchor):
    """Advance observed five-year CPS cohorts to ages40–44 at fixed rates.

    Published single ages are represented at their midpoints, with uniform
    weights inside each five-year band. Future births are a mean-stock
    approximation, without reapplying the CPS initial five-child topcode.
    """
    profile=json.loads((source/'cps_2024_age_profile_scenario.json').read_text())
    means={r['age_start']:r['mean_ceb'] for r in profile['rows']}
    assert means[40]==anchor
    raw=next(r for r in csv.DictReader((source/'nchs_asfr_1990_2023.csv').open()) if r['year']=='2023')
    rates=[float(v)/1000 for k,v in raw.items() if k.startswith('age')]
    scale=period_level/(5*sum(rates))
    rates=[r*scale for r in rates]
    assert abs(5*sum(rates)-period_level)<1e-12
    def integral(lo,hi):
        return sum(max(0.,min(hi,15+5*g)-max(lo,10+5*g))*rate for g,rate in enumerate(rates))
    steady=sum(integral(10,a+.5) for a in range(40,45))/5
    rows=[]
    for year in range(2024,2065,5):
        elapsed=year-2024;initial_low=40-elapsed
        if initial_low>=15:
            initial=means[initial_low]
            added=sum(integral(a+.5,a+.5+elapsed) for a in range(initial_low,initial_low+5))/5
            value=initial+added
            # Independent exposure counting on half-year steps.
            direct=0.
            for a in range(initial_low,initial_low+5):
                for step in range(2*elapsed):
                    midpoint=a+.5+.5*(step+.5)
                    if 10<=midpoint<50:direct+=.5*rates[int((midpoint-10)//5)]/5
            assert abs(direct-added)<1e-12
        else:
            # Cohorts younger than15 in2024: negligible early births imputed
            # at the same fixed age schedule. No arbitrary convergence path.
            initial=None;added=None;value=steady
        rows.append(dict(year=year,initial_age_start=initial_low,initial_mean_ceb=initial,
                         additional_births_to_age40_44=added,projected_ceb40_44=value))
    assert rows[0]['projected_ceb40_44']==anchor
    return rows,dict(asfr_scale_to_existing_wdi_tfr=scale,constant_period_tfr=period_level,
                     limiting_mean_ceb40_44=steady,projection='CPS2024 age means plus subsequent births under fixed2023 age pattern, scaled to the unchanged WDI2023 TFR.',
                     limitations=profile['projection_limitation'],half_year_exposure_check=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--projection',action='store_true',help='Write a separate dotted-extension variant; preserve the historical slide figure.')
    parser.add_argument('--reveal',action='store_true',help='Write source-free history/projection panels with identical axes for a Beamer reveal.')
    args=parser.parse_args()
    if args.reveal:args.projection=True
    source=BASE/'source/stock_forecast'
    p=source/'us_period_fertility_wdi.json'
    c=source/'cps_history_40_44.json'
    period={int(r['date']):float(r['value']) for r in json.loads(p.read_text())[1]}
    completed={int(y):float(v) for y,v in json.loads(c.read_text()).items() if int(y)>=1980}
    assert sorted(period)==list(range(1980,2024))
    assert period[2007]==2.12 and completed[2024]==1.918 and completed[1980]==2.988
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.5,5.1 if args.projection and not args.reveal else 4.5),sharey=True)
    projection_artists=[]
    scenario,scenario_meta=project_cps_age_profile(source,period[2023],completed[2024]) if args.projection else ([],{})
    plotted={}
    for ax,values,title,xlabel in zip(axes,[period,completed],
                                     ['Period fertility','Completed fertility'],
                                     ['Calendar year','Survey year (women aged 40–44)']):
        years=sorted(values)
        line,=ax.plot(years,[values[y] for y in years],color='#245f99',lw=2.2,label='Data',
                      marker='o' if ax is axes[1] else None,ms=3.5)
        assert list(line.get_xdata())==years
        assert list(line.get_ydata())==[values[y] for y in years]
        ax.axhline(2.1,color='#a77a49',lw=1.1,ls='--')
        ax.axvline(2007,color='.55',lw=1,ls=':')
        ax.text(2007,3.10,'2007',ha='center',va='top',color='.4',fontsize=11)
        ax.set(title=title,xlabel=xlabel,xlim=(1979,2025.5),ylim=(1.4,3.15),
               xticks=[1980,1990,2000,2010,2020],yticks=[1.5,2.,2.5,3.])
        ax.grid(axis='y',alpha=.16)
        ax.tick_params(labelleft=True)
        last=years[-1]
        ax.annotate(f'{values[last]:.2f}',(last,values[last]),xytext=(-4,-18),
                    textcoords='offset points',ha='right',color='#245f99',fontsize=11)
        plotted[title]=dict(years=years,values=[values[y] for y in years])
        if args.projection:
            future_years=[2023,2064] if ax is axes[0] else [r['year'] for r in scenario]
            future_values=[period[2023]]*2 if ax is axes[0] else [r['projected_ceb40_44'] for r in scenario]
            assert future_values[0]==values[years[-1]]
            future_line,=ax.plot(future_years,future_values,color='#245f99',lw=2.2,ls=':',label='Constant-rate scenario')
            projection_artists.append(future_line)
            assert list(future_line.get_ydata())==future_values
            ax.set_xlim(1979,2067);ax.set_xticks([1980,2000,2020,2040,2060])
            ax.legend(loc='upper right',frameon=False,fontsize=8.5)
            projection_artists.append(ax.annotate(f'{future_values[-1]:.2f}',(future_years[-1],future_values[-1]),xytext=(0,-18),
                        textcoords='offset points',ha='center',color='#245f99',fontsize=11))
            plotted[title]['scenario']=dict(years=future_years,values=future_values)
    axes[0].set_ylabel('Births per woman')
    axes[1].set_ylabel('Children ever born per woman')
    axes[0].text(1980,2.15,'Approx. replacement: 2.1',color='#906437',fontsize=10)
    axes[1].annotate(f'{completed[1980]:.2f}',(1980,completed[1980]),xytext=(8,1),
                     textcoords='offset points',ha='left',color='#245f99',fontsize=11)
    fig.subplots_adjust(left=.065,right=.98,bottom=.26 if args.projection and not args.reveal else .20,top=.88,wspace=.25)
    if args.projection and not args.reveal:
        fig.text(.065,.035,
                 'Dotted lines: illustrative constant-rate scenario. Period fertility stays at its last observed level; the 2023 age pattern is fixed.\n'
                 'Completed fertility advances the 2024 CPS age profile to ages 40–44. Its limit is slightly lower because some births occur after those ages.\n'
                 'Sources: World Bank WDI; Census CPS Tables 1, 3a and Historical Table 2; NCHS age-specific birth rates.',
                 fontsize=8.4,color='#555555',linespacing=1.4)
    out=BASE/'figures'
    name='fertility_introduction_with_projection' if args.projection else 'fertility_introduction'
    if args.reveal:name='fertility_introduction_reveal_projection'
    for ext in ('png','pdf'):
        fig.savefig(out/f'{name}.{ext}',dpi=170,facecolor='white')
    if args.reveal:
        positions=[list(ax.get_position().bounds) for ax in axes]
        for artist in projection_artists:artist.set_visible(False)
        for ax in axes:
            handles,labels=ax.get_legend_handles_labels()
            keep=[i for i,label in enumerate(labels) if label=='Data']
            ax.legend([handles[i] for i in keep],[labels[i] for i in keep],loc='upper right',frameon=False,fontsize=8.5)
        assert [list(ax.get_position().bounds) for ax in axes]==positions
        for ext in ('png','pdf'):
            fig.savefig(out/f'fertility_introduction_reveal_history.{ext}',dpi=170,facecolor='white')
    plt.close(fig)
    source_paths=[p,c]
    if args.projection:
        source_paths += [source/'cps_2024_age_profile_scenario.json',source/'nchs_asfr_1990_2023.csv']
        original=json.loads((out/'fertility_introduction_verification.json').read_text())
        for title,values in plotted.items():
            assert values['years']==original['plotted'][title]['years']
            assert values['values']==original['plotted'][title]['values']
    (out/f'{name}_verification.json').write_text(json.dumps(dict(
        status='PASS',purpose='Separate illustrative dotted continuation from the CPS2024 age profile at constant period fertility.' if args.projection else 'Historical empirical motivation for approximating2007 by an initial steady state; not evidence of exact stationarity.',
        model_lines=False,housing_panels=False,plotted=plotted,
        completed_definition='CPS children ever born ages40–44, a near-completed measure;2022/2024 counts capped atfive.',
        period_source='World Bank WDI SP.DYN.TFRT.IN, USA,1980–2023; extended September13,2026; all prior1990–2023 values unchanged.',
        completed_source='US Census CPS Historical Table2, through2024.',
        source_sha256={str(q.relative_to(BASE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in source_paths},
        scenario=scenario,scenario_metadata=scenario_meta,
        reveal_identical_axes=args.reveal,visible_sources=False if args.reveal else None,
        artist_data_exact=True),indent=2)+'\n')
    if args.projection:
        with (out/f'{name}.csv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(scenario[0]),lineterminator='\n');writer.writeheader();writer.writerows(scenario)
        print(json.dumps(scenario,indent=2))
    print('PASS: two historical panels, exact source/artist equality; no model or housing series.')


if __name__=='__main__':
    main()
