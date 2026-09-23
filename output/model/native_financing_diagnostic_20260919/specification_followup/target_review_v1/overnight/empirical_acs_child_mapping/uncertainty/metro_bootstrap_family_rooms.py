"""Paired metro bootstrap for existing and under-18 family_rooms contrasts.

Uses the original early-housing builder's 42-metro bootstrap exactly and one
sequential memmap pass over ACS 2005-06 to form candidate per-metro totals.
No microdata are written.
"""
from pathlib import Path
import hashlib, json, time
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[8]
OUT=Path(__file__).resolve().parent
SOURCE=ROOT/'code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta'
BUILDER=ROOT/'output/model/e5f_matched_pf_20260909a/design_research/housing'
HOUSING=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/full'
PREVIOUS=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_acs_child_mapping'
TARGET_REVIEW=PREVIOUS.parent.parent
YEARS=(2005,2006)
N=1000
SEED=20260910
START=time.monotonic()


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def main():
    ids=[int(x) for x in (BUILDER/'active_metros.txt').read_text().strip().split(',')]
    assert len(ids)==42 and len(set(ids))==42
    idset=set(ids)
    # This matches the authoritative builder's saved old-group aggregate rows.
    metro=pd.read_csv(BUILDER/'early_housing_metro_components.csv')
    old=(metro[(metro.active_metro==1)&metro.year.isin(YEARS)]
         .groupby('met2013').sum(numeric_only=True).reindex(ids))
    assert old.notna().all().all() and len(old)==42

    reader=pd.io.stata.StataReader(str(SOURCE),chunksize=1)
    reader._path_or_buf=SOURCE.open('rb'); reader._read_header(); reader._setup_dtype()
    dtype=reader._dtype; fields=dict(zip(reader._varlist,dtype.names))
    required=('year','sample','serial','met2013','gq','pernum','relate','hhwt','age','ownershp','rooms','unitsstr','nchild','yngch','momloc','poploc')
    assert all(x in fields for x in required)
    data_loc=reader._data_location; nobs=reader._nobs
    reader._path_or_buf.close()
    data=np.memmap(SOURCE,dtype=dtype,mode='r',offset=data_loc,shape=(nobs,))
    yearcol=data[fields['year']]
    def lower(y):
        lo,hi=0,nobs
        while lo<hi:
            mid=(lo+hi)//2
            if int(yearcol[mid])<y:lo=mid+1
            else:hi=mid
        return lo

    # One-pass candidate extraction: retain only eligible head summaries and
    # counts of linked own children by household key; never persist person rows.
    heads={}; linked={}; ranges=[]
    for y in YEARS:
        lo,hi=lower(y),lower(y+1); ranges.append({'year':y,'start_record':lo,'end_record':hi,'records':hi-lo})
        for start in range(lo,hi,250000):
            b=data[start:min(start+250000,hi)]
            a={x:np.asarray(b[fields[x]]) for x in required}
            sample=(a['sample']==y*100+1)
            metro_ok=np.isin(a['met2013'],ids)
            common=sample&metro_ok&np.isin(a['gq'],(1,2))
            ishead=common&(a['pernum']==1)&(a['relate']==1)&(a['hhwt']>0)&(a['age']>=30)&(a['age']<=55)&np.isin(a['ownershp'],(1,2))&(a['rooms']>0)&(a['nchild']>0)&(a['yngch']<18)
            for i in np.flatnonzero(ishead):
                key=(int(a['year'][i]),int(a['sample'][i]),int(a['serial'][i]))
                if key in heads: raise AssertionError(f'Duplicate selected head key {key}')
                heads[key]={'metro':int(a['met2013'][i]),'weight':float(a['hhwt'][i]),'rooms9':float(min(a['rooms'][i],9)),'nchild':int(a['nchild'][i])}
            ischild=common&(a['pernum']>1)&np.isin(a['age'],np.arange(0,100))&((a['momloc']==1)|(a['poploc']==1))
            for i in np.flatnonzero(ischild):
                key=(int(a['year'][i]),int(a['sample'][i]),int(a['serial'][i]))
                if key not in linked: linked[key]=[0,0]
                linked[key][0]+=1
                if int(a['age'][i])<18: linked[key][1]+=1

    # Join only in-memory aggregate summaries, then validate household-level
    # all-age child bins against the field used by the original moment.
    fields_by_metro={m:{'old_hi_n':0,'old_hi_w':0.,'old_hi_r':0.,'old_lo_n':0,'old_lo_w':0.,'old_lo_r':0.,'new_hi_n':0,'new_hi_w':0.,'new_hi_r':0.,'new_lo_n':0,'new_lo_w':0.,'new_lo_r':0.,'family_n':0,'linked_bin_mismatch_n':0} for m in ids}
    for key,head in heads.items():
        total,minor=linked.get(key,(0,0))
        if min(total,3)!=min(head['nchild'],3):
            fields_by_metro[head['metro']]['linked_bin_mismatch_n']+=1
        g=fields_by_metro[head['metro']]; g['family_n']+=1; w=head['weight']; r=head['rooms9']
        if head['nchild']>=3: g['old_hi_n']+=1; g['old_hi_w']+=w; g['old_hi_r']+=w*r
        elif 1<=head['nchild']<=2: g['old_lo_n']+=1; g['old_lo_w']+=w; g['old_lo_r']+=w*r
        if minor>=3: g['new_hi_n']+=1; g['new_hi_w']+=w; g['new_hi_r']+=w*r
        elif 1<=minor<=2: g['new_lo_n']+=1; g['new_lo_w']+=w; g['new_lo_r']+=w*r

    rows=[]
    for m in ids:
        x=fields_by_metro[m]; o=old.loc[m]
        # Original sufficient stats are sourced from the authoritative builder;
        # the direct old grouping above checks that its sample/group totals agree.
        for group,weight_col,room_col,ncol,countcol in [
          ('old_3plus','rooms_large_weight','rooms_large_capped_rooms_sum','rooms_large_n','old_hi_n'),
          ('old_1to2','rooms_small_weight','rooms_small_capped_rooms_sum','rooms_small_n','old_lo_n')]:
            suffix='hi' if group.endswith('3plus') else 'lo'
            assert abs(x[f'old_{suffix}_w']-float(o[weight_col]))<1e-6
            assert abs(x[f'old_{suffix}_r']-float(o[room_col]))<1e-6
            assert int(x[countcol])==int(o[ncol])
        rows.append({'met2013':m,**x})
    stats=pd.DataFrame(rows)
    assert stats.linked_bin_mismatch_n.sum()==0
    assert stats.family_n.sum()==275009
    # Reproduce the earlier aggregate-only candidate receipt by group before
    # using these metro totals for resampling.
    prior_groups=pd.read_csv(PREVIOUS/'family_rooms_group_aggregates.csv').set_index('group')
    candidate_checks=[('under18_children_3plus','new_hi_n','new_hi_w','new_hi_r'),('under18_children_1to2','new_lo_n','new_lo_w','new_lo_r')]
    for label,ncol,wcol,rcol in candidate_checks:
        n=int(stats[ncol].sum()); z=stats[[wcol,rcol]].sum(); prior=prior_groups.loc[label]
        assert n==int(prior.households)
        assert abs(float(z[wcol])-float(prior.household_weight))<1e-6
        assert abs(float(z[rcol]/z[wcol])-float(prior.mean_rooms_cap9))<1e-12
    stats.to_csv(OUT/'metro_sufficient_stats.csv',index=False)

    # Exact builder design: one multinomial-equivalent frequency vector per
    # draw, sampling 42 of 42 metro IDs with replacement; seed and draw order
    # reproduce summarize_early_housing.py's 2005_2006 window.
    rng=np.random.default_rng(SEED)
    frequency=np.array([np.bincount(rng.integers(0,42,size=42),minlength=42) for _ in range(N)])
    c=stats.set_index('met2013').reindex(ids)
    def contrast(arr,hiw,hir,loww,lowr):
        return arr[:,hir]/arr[:,hiw]-arr[:,lowr]/arr[:,loww]
    old_sums=c[['old_hi_w','old_hi_r','old_lo_w','old_lo_r']].to_numpy()
    new_sums=c[['new_hi_w','new_hi_r','new_lo_w','new_lo_r']].to_numpy()
    old_boot=frequency@old_sums; new_boot=frequency@new_sums
    old_draws=contrast(old_boot,0,1,2,3); new_draws=contrast(new_boot,0,1,2,3)
    difference=new_draws-old_draws
    old_point=contrast(old_sums.sum(axis=0,keepdims=True),0,1,2,3)[0]
    new_point=contrast(new_sums.sum(axis=0,keepdims=True),0,1,2,3)[0]
    old_se=float(old_draws.std(ddof=1)); new_se=float(new_draws.std(ddof=1)); diff_se=float(difference.std(ddof=1))
    covariance=np.cov(np.column_stack([old_draws,new_draws]),rowvar=False,ddof=1)
    old_target=float(pd.read_csv(BUILDER/'early_housing_target_candidates.csv').query("window == '2005_2006' and moment == 'prime30_55_resident_3plus_minus_1to2_rooms_capped9'").iloc[0].point)
    old_saved_se=float(pd.read_csv(BUILDER/'early_housing_target_candidates.csv').query("window == '2005_2006' and moment == 'prime30_55_resident_3plus_minus_1to2_rooms_capped9'").iloc[0].metro_bootstrap_se)
    family_row=next(x for x in json.loads((TARGET_REVIEW/'housing_receipt.json').read_text())['findings'] if x['row']=='family_rooms')
    old_saved_weight=float(family_row['weight'])
    expected_candidate=0.336219623885607
    assert abs(old_point-old_target)<1e-12
    assert abs(old_se-old_saved_se)<1e-12
    assert abs(1/old_se**2-old_saved_weight)<1e-9
    assert abs(new_point-expected_candidate)<1e-9
    pd.DataFrame({'draw':np.arange(1,N+1),'existing_family_rooms':old_draws,'under18_candidate':new_draws,'candidate_minus_existing':difference}).to_csv(OUT/'paired_bootstrap_draws.csv',index=False)
    names=['existing_family_rooms','under18_candidate']
    summary={
      'status':'pass',
      'settings':{'bootstrap_draws':N,'seed':SEED,'clusters':42,'cluster':'MET2013 metro','sampling':'Resample 42 of 42 metro IDs with replacement; draw-specific frequency counts preserve every within-metro household/group weight and room-sum total.','pairing':'Identical frequency draws applied to existing and candidate contrasts.','variance':'Sample standard deviation/covariance with ddof=1.','official_ACS_design_SE':False,'membership':'Fixed ACS 2005-06 head and parent-linked under-18 group membership within each metro; only metro multiplicities vary across draws.'},
      'points':{'existing_family_rooms':float(old_point),'under18_candidate':float(new_point),'candidate_minus_existing':float(new_point-old_point)},
      'paired_uncertainty':{'existing_family_rooms_se':old_se,'under18_candidate_se':new_se,'candidate_minus_existing_se':diff_se,'candidate_minus_existing_percentile_95_interval':[float(np.quantile(difference,.025)),float(np.quantile(difference,.975))],'covariance_matrix':covariance.tolist(),'covariance_order':names,'correlation':float(covariance[0,1]/(old_se*new_se)),'inverse_variance_weights':{'existing_family_rooms':1/old_se**2,'under18_candidate':1/new_se**2}},
      'saved_reproduction':{'existing_point_source':old_target,'existing_point_gap':float(old_point-old_target),'existing_saved_se':old_saved_se,'existing_saved_se_gap':float(old_se-old_saved_se),'existing_saved_weight':old_saved_weight,'existing_inverse_variance_weight':float(1/old_se**2),'existing_weight_gap':float(1/old_se**2-old_saved_weight),'assert_tolerances':{'existing_point':1e-12,'existing_se':1e-12,'existing_weight':1e-9,'candidate_point':1e-9}},
      'sample':{'family_households':int(stats.family_n.sum()),'parent_linked_bin_disagreements':int(stats.linked_bin_mismatch_n.sum()),'active_metros':ids,'year_ranges':ranges},
      'input_hashes':{},
      'interpretation':'Metro-cluster empirical resampling uncertainty for fixed observed sample/group membership; not an ACS official design/replicate-weight standard error. Candidate remains diagnostic only.',
      'elapsed_seconds':time.monotonic()-START,
    }
    hashes={
      'raw_source_canonical_sha256':json.loads((HOUSING/'provenance.json').read_text())['source_sha256_from_existing_canonical_receipt'],
      'raw_source_rehashed_this_pass':False,
      'upstream_builder_script_sha256':sha(BUILDER/'summarize_early_housing.py'),
      'upstream_metro_component_sha256':sha(BUILDER/'early_housing_metro_components.csv'),
      'upstream_builder_receipt_sha256':sha(BUILDER/'early_housing_candidate_receipt.json'),
      'parent_link_aggregate_input_sha256':sha(PREVIOUS/'family_rooms_group_aggregates.csv'),
      'parent_link_aggregate_json_sha256':sha(PREVIOUS/'housing_profiles_v1_child_mapping_recomputed.json'),
      'local_script_sha256':sha(Path(__file__).resolve()),
      'sufficient_stats_sha256':None,
    }
    summary['input_hashes']=hashes
    (OUT/'bootstrap_uncertainty.json').write_text(json.dumps(summary,indent=2)+'\n')
    hashes['sufficient_stats_sha256']=sha(OUT/'metro_sufficient_stats.csv')
    summary['input_hashes']=hashes
    (OUT/'bootstrap_uncertainty.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'status':'PASS','old_point':old_point,'old_se':old_se,'candidate_point':new_point,'candidate_se':new_se,'difference_se':diff_se,'seconds':time.monotonic()-START}),flush=True)

if __name__=='__main__':main()
