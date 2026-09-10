"""Summarize saved ACS components; never re-read raw data or adopt targets."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
SOURCE=OUT/'early_housing_metro_components.csv'
d=pd.read_csv(SOURCE)
ids=list(map(int,(OUT/'active_metros.txt').read_text().split(',')))
N=1000
SEED=20260910
rng=np.random.default_rng(SEED)
frequency=np.array([np.bincount(rng.integers(0,42,size=42),minlength=42) for _ in range(N)])
names=['aggregate_mean_occupied_rooms_capped9_18_85','own_rate_30_55','recent_parent_minus_no_resident_child_ownership_30_55','prime30_55_resident_3plus_minus_1to2_rooms_capped9','own_rate_25_34']

def moments(a):
    return np.column_stack([
        a['capped_rooms_sum']/a['weight'],
        a['own3055_owner_weight']/a['own3055_weight'],
        a['newparent_owner_weight']/a['newparent_weight']-a['nochild_owner_weight']/a['nochild_weight'],
        a['rooms_large_capped_rooms_sum']/a['rooms_large_weight']-a['rooms_small_capped_rooms_sum']/a['rooms_small_weight'],
        a['own2534_owner_weight']/a['own2534_weight'],
    ])

results=[]
for window,years in [('2005_2006',[2005,2006]),('2005_2007',[2005,2006,2007]),('2007',[2007]),('2012',[2012]),('2023',[2023])]:
    c=d[(d.active_metro==1)&d.year.isin(years)].groupby('met2013').sum(numeric_only=True).reindex(ids)
    assert not c.isna().any().any() and len(c)==42
    totals=c.sum(axis=0)
    point=moments(totals.to_dict())[0]
    b=pd.DataFrame(frequency@c.to_numpy(),columns=c.columns)
    draws=moments(b)
    assert draws.shape==(N,5) and np.all(np.isfinite(draws))
    se=draws.std(axis=0,ddof=1)
    cov=np.cov(draws,rowvar=False)
    assert np.allclose(np.diag(cov),se**2,atol=1e-14,rtol=0)
    pd.DataFrame(cov,index=names,columns=names).to_csv(OUT/f'metro_bootstrap_covariance_{window}.csv')
    pd.DataFrame(draws,columns=names).to_csv(OUT/f'metro_bootstrap_draws_{window}.csv',index=False)
    for j,name in enumerate(names):
        results.append({'window':window,'moment':name,'point':point[j],'metro_bootstrap_se':se[j],'p025':np.quantile(draws[:,j],.025),'p975':np.quantile(draws[:,j],.975),'B':N,'seed':SEED,'metro_clusters':42,'status':'diagnostic_new_contract_candidate_not_adopted','uncertainty_interpretation':'empirical metro-cluster resampling; not ACS official design/replicate-weight SE'})
frame=pd.DataFrame(results)
frame.to_csv(OUT/'early_housing_target_candidates.csv',index=False)
meta={'point_source':str(SOURCE),'point_source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'bootstrap':{'B':N,'seed':SEED,'clusters':42,'method':'Resample 42 MET2013 cities with replacement; preserve all household and group totals within each drawn metro; paired same draws across dates','official_ACS_design_SE':False},'new_calibration_contract_adopted':False,'model_observation_rule_required':'For mean-room and family-room targets use min(occupied rooms,9) BEFORE aggregation. Ownership unchanged. Preserve current room response empirical target separately.'}
(OUT/'early_housing_candidate_receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
print(frame[frame.window.isin(['2005_2007','2023'])].to_string(index=False),flush=True)
