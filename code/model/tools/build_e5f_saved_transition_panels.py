"""Supplemental seminar panels from a verified saved forecast; no model solves."""
from pathlib import Path
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912'
SOURCE=BASE/'matched_short'
OUT=BASE/'seminar_transition_panels'

def main():
    OUT.mkdir(exist_ok=True)
    rows=list(csv.DictReader((SOURCE/'expected_transition.csv').open()))
    fert=json.loads((SOURCE/'fertility.json').read_text())
    receipt=json.loads((SOURCE/'root_receipt.json').read_text())
    assert receipt['finite_horizon_market_fiscal_converged']
    assert not receipt['terminal_distance_passed']
    assert len(rows)==len(fert)==6
    series=[]
    targets={2011:1.974875,2015:1.861,2019:1.755375,2023:1.64575}
    for r,f in zip(rows,fert):
        year=int(r['calendar_year']);assert year==f['calendar_year']
        heads=float(r['household_heads'] or r['adult_population'])
        persons=float(r['resident_persons']) if r['resident_persons'] else float('nan')
        series.append(dict(decision_year=year,birth_window_end=year+4,
            psi=float(r['psi_child']),fertility_household_rate=f['period_tfr_topcode_adjusted'],
            published_tfr_window_average=targets.get(year+4,float('nan')),
            births_raw=float(r['birth_children']),births_topcode_adjusted=float(r['birth_children_topcode_adjusted']),
            heads=heads,resident_persons=persons,asset_price=float(r['asset_price']),
            rent=float(r['renter_price']),ownership_percent=100*float(r['owner_rate']),
            rooms_per_head=float(r['housing_demand'])/heads,pension=float(r['pension_period_units'])))
    with (OUT/'series.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(series[0]));writer.writeheader();writer.writerows(series)
    x=np.array([r['decision_year'] for r in series]);get=lambda k:np.array([r[k] for r in series])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def layout():
        fig,axs=plt.subplots(2,2,figsize=(11.8,7.3))
        for a in axs.flat:a.grid(alpha=.18);a.set_xticks(x);a.set_xlabel('Decision year')
        return fig,axs
    def finish(fig,name):
        fig.suptitle('Expected transition after one permanent preference surprise',fontsize=15,y=.99)
        fig.text(.5,.925,'Sequential model | historical workers’ pinned initial point | preference change −0.01414',ha='center',fontsize=10)
        fig.text(.5,.025,'Finite six-date housing and pension equilibrium; terminal horizon NOT certified.\nLevels conditional on one permanent shock—not a fitted history or a baseline-subtracted IRF.',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.08,1,.90));fig.savefig(OUT/f'{name}.pdf');fig.savefig(OUT/f'{name}.png',dpi=150);plt.close(fig)
    fig,ax=layout()
    a=ax[0,0];a.plot(x,get('psi'),'o-',label='Current preference, believed permanent')
    a.axhline(.16122715651861665,color='gray',ls=':',label='Pre-shock preference')
    a.set(title='Preference input',ylabel='Child-preference coefficient');a.legend(fontsize=8)
    a=ax[0,1];a.plot(x+4,get('fertility_household_rate'),'o-',label='Expected household-rate analogue')
    a.plot(x+4,get('published_tfr_window_average'),'ks--',label='Published female TFR: window average')
    a.set(title='Fertility: model forecast and historical data',xlabel='End of four-year birth window',ylabel='Period fertility');a.set_xticks(x+4);a.legend(fontsize=8)
    a=ax[1,0];a.plot(x,get('births_topcode_adjusted'),'o-',label='Topcode-adjusted birth flow')
    a.plot(x,get('births_raw'),'s--',label='Explicit-state birth flow')
    a.set(title='Births over each four-year interval',ylabel='Births per initial model household');a.legend(fontsize=8)
    a=ax[1,1];a.plot(x,100*get('heads')/get('heads')[4],'o-',label='Household heads')
    a.plot(x,100*get('resident_persons')/get('resident_persons')[4],'s-',label='Resident persons: available from2023')
    a.axvspan(2007,2023,alpha=.07,color='gray');a.set(title='Population and heads (2023 = 100)',ylabel='Index')
    a.legend(fontsize=8);a.text(.02,.04,'Pre-2023 heads externally conditioned;\n2023 persons externally anchored.',transform=a.transAxes,fontsize=8)
    finish(fig,'fertility_demography')
    fig,ax=layout()
    for a,key,title,ylabel in zip(ax.flat,['asset_price','rent','ownership_percent','rooms_per_head'],
        ['House asset price','Rent per physical room','Ownership','Occupied physical rooms'],
        ['Model asset-price units','Model rent units','Percent of household heads','Rooms per household head']):
        a.plot(x,get(key),'o-');a.set(title=title,ylabel=ylabel)
    finish(fig,'housing')
    provenance=dict(source_files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [SOURCE/'expected_transition.csv',SOURCE/'fertility.json',SOURCE/'root_receipt.json']},
        rows=6,finite_market_fiscal_converged=True,terminal_distance_passed=False,
        historical_fit_complete=False,irf_against_no_shock=False,
        initial_calibration='Pinned historical worker initial; NOT verified_final overnight refinement',
        fertility_definition='Sum of four-year age-specific topcode-adjusted birth flows divided by adult-household age mass; approximate analogue of female TFR, no extra factor of four.',
        demographics='Historical heads conditioned externally; resident persons available only from the externally fixed2023anchor. Post2023 person/head accounting is model demographic propagation.',
        shock='Permanent step from .16122715651861665 to .14708715651861665 at2007; no subsequent shocks.',
        figures=['fertility_demography.pdf','housing.pdf'])
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(OUT)

if __name__=='__main__':main()
