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
    red,blue,green='#bd433c','#245f99','#3b8065'
    a,b,c,d=axes.flat
    for k,label,color in [('data_rooms','Occupied rooms',blue),('data_households','Households',green),('data_persons','Residents',red)]:
        v=vals(k); a.plot(x,100*v/v[0],'o-',color=color,label=f'{label}: +{growth(k):.1f}%',ms=4)
    a.set(title='What grew in the data?',ylabel='2007 = 100')
    a.legend(frameon=False,fontsize=9)
    b.plot(x,vals('data_persons_per_household'),'s-',color=red,ms=4)
    b.set(title='Residents per household',ylabel='Persons')
    for year,value in [(x[0],vals('data_persons_per_household')[0]),(x[-1],vals('data_persons_per_household')[-1])]:
        b.annotate(f'{value:.2f}',(year,value),xytext=(0,8),textcoords='offset points',ha='center',color=red)
    for key,label,color in [('data_rooms_per_household','ACS',red),('model_rooms_per_household','Model',blue)]:
        v=vals(key); c.plot(x,v,'o-',color=color,ms=4,label=f'{label}: {growth(key):+.1f}%')
    c.set(title='Housing per household',ylabel='Occupied rooms per household')
    c.legend(frameon=False,fontsize=9)
    for key,label,color in [('data_rooms_per_person','ACS',red),('model_rooms_per_person_at_data_demographics','Model at observed demographics',blue)]:
        v=vals(key); d.plot(x,v,'o-',color=color,ms=4,label=f'{label}: {growth(key):+.1f}%')
    d.set(title='Housing per resident: common demographics',ylabel='Occupied rooms per resident')
    d.legend(frameon=False,fontsize=8.5)
    for ax in axes.flat:
        ax.set_xlim(2005.7,2024.3); ax.set_xticks(x); ax.grid(alpha=.16); ax.margins(y=.18)
    fig.suptitle('Housing and population, 2007–2023',x=.08,y=.97,ha='left',fontsize=19)
    fig.text(.08,.055,
             'National ACS: the same occupied households and their residents in every measure; heads ages 18–85, rooms capped at nine.\n'
             'Rooms and households use household weights; residents use person weights. Group quarters are excluded.\n'
             'Bottom right: model rooms per household × observed households / observed residents. This holds demographics common;\n'
             'it is not a model population prediction. The retained model has no separate resident-person history before 2023.',
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
                      source_scope='National housing validation. Initial model calibration retains its42-metro targets; no model rerun, reweighting of states, or target change.',
                      checks=['Five household and person ratio identities','Multiplicative growth decomposition','Complete2007–2023model date support'],
                      rows=rows)
    (out/'housing_population_comparison_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(growths,indent=2))


if __name__=='__main__':
    main()
