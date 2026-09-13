"""Two historical fertility panels motivating the initial-state approximation."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'


def main():
    source=BASE/'source/stock_forecast'
    p=source/'us_period_fertility_wdi.json'
    c=source/'cps_history_40_44.json'
    period={int(r['date']):float(r['value']) for r in json.loads(p.read_text())[1]}
    completed={int(y):float(v) for y,v in json.loads(c.read_text()).items() if int(y)>=1990}
    assert sorted(period)==list(range(1990,2024))
    assert period[2007]==2.12 and completed[2024]==1.918
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.5),sharey=True)
    plotted={}
    for ax,values,title,xlabel in zip(axes,[period,completed],
                                     ['Period fertility','Completed fertility'],
                                     ['Calendar year','Survey year (women aged 40–44)']):
        years=sorted(values)
        line,=ax.plot(years,[values[y] for y in years],color='#245f99',lw=2.2,
                      marker='o' if ax is axes[1] else None,ms=3.5)
        assert list(line.get_xdata())==years
        assert list(line.get_ydata())==[values[y] for y in years]
        ax.axhline(2.1,color='#a77a49',lw=1.1,ls='--')
        ax.axvline(2007,color='.55',lw=1,ls=':')
        ax.text(2007,2.28,'2007',ha='center',va='top',color='.4',fontsize=11)
        ax.set(title=title,xlabel=xlabel,xlim=(1989,2025.5),ylim=(1.4,2.3),xticks=[1990,2000,2010,2020])
        ax.grid(axis='y',alpha=.16)
        ax.tick_params(labelleft=True)
        last=years[-1]
        ax.annotate(f'{values[last]:.2f}',(last,values[last]),xytext=(-4,-18),
                    textcoords='offset points',ha='right',color='#245f99',fontsize=11)
        plotted[title]=dict(years=years,values=[values[y] for y in years])
    axes[0].set_ylabel('Births per woman')
    axes[1].set_ylabel('Children ever born per woman')
    axes[0].text(1990,2.125,'Approx. replacement: 2.1',color='#906437',fontsize=10)
    fig.subplots_adjust(left=.065,right=.98,bottom=.20,top=.88,wspace=.25)
    out=BASE/'figures'
    for ext in ('png','pdf'):
        fig.savefig(out/f'fertility_introduction.{ext}',dpi=170,facecolor='white')
    plt.close(fig)
    (out/'fertility_introduction_verification.json').write_text(json.dumps(dict(
        status='PASS',purpose='Historical empirical motivation for approximating2007 by an initial steady state; not evidence of exact stationarity.',
        model_lines=False,housing_panels=False,plotted=plotted,
        completed_definition='CPS children ever born ages40–44, a near-completed measure;2022/2024 counts capped atfive.',
        period_source='World Bank WDI SP.DYN.TFRT.IN, USA,1990–2023; previously retrieved September12,2026.',
        completed_source='US Census CPS Historical Table2, through2024.',
        source_sha256={str(q.relative_to(BASE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (p,c)},
        artist_data_exact=True),indent=2)+'\n')
    print('PASS: two historical panels, exact source/artist equality; no model or housing series.')


if __name__=='__main__':
    main()
