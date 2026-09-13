"""Compare saved housing choices at observed ACS household/population counts.

This is demographic standardization, not a replay or a population prediction.
"""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'


def main():
    path = BASE/'source/historical_stock/housing_population_data.json'
    data = json.loads(path.read_text())
    assert data['status'] == 'PASS'
    model_path = BASE/'figures/stock_forecast.csv'
    model = {int(r['year']):r for r in csv.DictReader(model_path.open())}
    rows = []
    for r in data['rows']:
        m = model[r['year']]
        mh = float(m['mean_capped_rooms'])
        p, h, q = r['person_weighted_residents'], r['households'], r['occupied_rooms_capped9']
        mq = mh*h
        assert abs(q/p - r['rooms_per_household']/(p/h)) < 1e-12
        assert abs(mq/p - mh/(p/h)) < 1e-12
        rows.append(dict(year=r['year'],data_persons=p,data_households=h,data_rooms=q,
                         data_persons_per_household=p/h,data_rooms_per_household=q/h,
                         data_rooms_per_person=q/p,model_rooms_per_household=mh,
                         model_rooms_at_data_households=mq,
                         model_rooms_per_person_at_data_demographics=mq/p,
                         data_rooms_per_person_household_weights=r['rooms_per_person_household_weighted']))
    assert [r['year'] for r in rows] == [2007,2011,2015,2019,2023]
    x = [r['year'] for r in rows]
    def vals(key):
        return np.array([r[key] for r in rows])
    def growth(key):
        v=vals(key)
        return float(100*(v[-1]/v[0]-1))
    growths = {k:growth(k) for k in rows[0] if k!='year'}
    # Multiplicative population x household-formation x rooms decomposition.
    parts = [np.log(vals(k)[-1]/vals(k)[0]) for k in
             ('data_persons','data_persons_per_household','data_rooms_per_household')]
    assert abs(parts[0]-parts[1]+parts[2]-np.log(vals('data_rooms')[-1]/vals('data_rooms')[0])) < 1e-12
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.2,8.2))
    red,blue='#bd433c','#245f99'
    a,b,c,d=axes.flat
    model_years=sorted(model)
    assert model_years == list(range(2007,2040,4))
    plotted={}
    for ax,dk,mk,title,ylabel,indexed in [
            (a,'data_households','households','Households','2007 = 100',True),
            (b,'data_rooms','capped_rooms','Total occupied housing','2007 = 100',True),
            (c,'data_rooms_per_household','mean_capped_rooms','Housing per household','Occupied rooms per household',False)]:
        mv=np.array([float(model[y][mk]) for y in model_years]);dv=vals(dk)
        if indexed:
            mv=100*mv/mv[0];dv=100*dv/dv[0]
        ax.plot(model_years[:5],mv[:5],'o-',color=blue,ms=4,label='Model')
        ax.plot(model_years[4:],mv[4:],'o:',color=blue,ms=4,label='Model continuation')
        ax.plot(x,dv,'s--',color=red,ms=4,label='Data')
        ax.set(title=title,ylabel=ylabel)
        ax.axvline(2023,color='.65',lw=.8,ls='--')
        ax.axvspan(2023,2040,color=blue,alpha=.035)
        ax.set_xlim(2005.7,2040.5);ax.set_xticks([2007,2015,2023,2031,2039])
        ax.legend(frameon=False,fontsize=8.5)
        plotted[title]=dict(model_years=model_years,model_values=mv.tolist(),data_years=x,data_values=dv.tolist())
    for key,label,color in [('data_rooms_per_person','ACS',red),('model_rooms_per_person_at_data_demographics','Model at observed demographics',blue)]:
        v=vals(key); d.plot(x,v,'s--' if color==red else 'o-',color=color,ms=4,label=label)
    d.set(title='Housing per resident',ylabel='Occupied rooms per resident')
    d.text(.03,.04,'Observed demographics for both series',transform=d.transAxes,fontsize=9,color='#555555')
    d.set_xlim(2005.7,2024.3);d.set_xticks(x)
    d.legend(frameon=False,fontsize=8.5)
    for ax in axes.flat:
        ax.grid(alpha=.16); ax.margins(y=.18)
        labels=ax.get_legend_handles_labels()[1]
        assert any('Model' in label for label in labels) and any(label in ('Data','ACS') for label in labels)
    fig.suptitle('Housing: model and data',x=.08,y=.97,ha='left',fontsize=19)
    fig.text(.08,.055,
             'Data: national ACS occupied households, heads ages 18–85; rooms capped at nine. Household weights for housing, person weights for residents.\n'
             'Top row: native model aggregates; historical household counts use Census inputs, whose coverage differs from the restricted ACS sample.\n'
             'Bottom right: model rooms per household × ACS households / ACS residents; a comparison at common demographics, not a population prediction.\n'
             'Model resident-person history is unavailable before 2023. Dotted continuation holds preferences fixed; no tax rebate, horizon sensitivity unresolved.',
             fontsize=9,color='#555555',linespacing=1.5)
    fig.subplots_adjust(left=.08,right=.97,top=.89,bottom=.23,hspace=.35,wspace=.26)
    out=BASE/'figures'
    for ext in ('png','pdf'):
        fig.savefig(out/f'housing_population_comparison.{ext}',dpi=170,facecolor='white')
    plt.close(fig)
    with (out/'housing_population_comparison.csv').open('w') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    verification=dict(status='PASS',data_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      model_csv_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
                      growth_2007_2023_percent=growths,
                      standardization='Q_model_standardized = model mean capped rooms per household * ACS HHWT household total. Divide by ACS PERWT residents in those same households for per-person comparison.',
                      model_person_history_before_2023_available=False,
                      every_panel_compares_model_and_data=True,
                      native_model_history_and_continuation_panels=plotted,
                      source_scope='National housing validation. Initial model calibration retains its42-metro targets; no model rerun, reweighting of states, or target change.',
                      checks=['Five household and person ratio identities','Multiplicative growth decomposition','Complete2007–2023model date support'],
                      rows=rows)
    (out/'housing_population_comparison_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(growths,indent=2))


if __name__=='__main__':
    main()
