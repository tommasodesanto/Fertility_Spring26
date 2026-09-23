"""ACS family_rooms child-mapping check; no model or target edits.

Reuses the housing_profiles_v1 raw extract and exact person/household filters.
Counts resident children linked to the head by MOMLOC/POPLOC, split by age.
RELATE==3 is retained only as a secondary relationship-category comparison.
"""
from pathlib import Path
import csv, hashlib, json, time
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[7]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta'
METROS = [12060,12420,12580,14460,15380,16740,16980,17140,18140,19100,19740,19820,26420,26900,27260,28140,29820,31080,33100,33340,33460,34980,35380,35620,36420,36740,37980,38060,38300,38900,39300,39580,40900,41180,41620,41700,41740,41860,42660,45300,47260,47900]
YEARS=(2005,2006)
MAX_SECONDS=2700
START=time.monotonic()


def fail_time():
    if time.monotonic()-START > MAX_SECONDS:
        raise TimeoutError(f'Bounded read exceeded {MAX_SECONDS} seconds')


def main():
    r=pd.io.stata.StataReader(str(SOURCE),chunksize=1)
    r._path_or_buf=SOURCE.open('rb')
    r._read_header(); r._setup_dtype()
    names=r._varlist
    dtype=r._dtype
    fields=dict(zip(names,dtype.names))
    needed=('year','sample','serial','statefip','puma','met2013','gq','pernum','relate','hhwt','age','ownershp','rooms','unitsstr','nchild','yngch','eldch','momloc','poploc')
    missing=[x for x in needed if x not in fields]
    if missing: raise RuntimeError(f'Required source fields unavailable: {missing}')
    data_loc=r._data_location; nobs=r._nobs
    byteorder=r._byteorder
    r._path_or_buf.close()
    data=np.memmap(SOURCE,dtype=dtype,mode='r',offset=data_loc,shape=(nobs,))

    def lower(y):
        lo,hi=0,nobs
        yr=data[fields['year']]
        while lo<hi:
            mid=(lo+hi)//2
            if int(yr[mid])<y:lo=mid+1
            else:hi=mid
        return lo

    all_heads=[]
    year_records=[]
    for y in YEARS:
        fail_time(); lo,hi=lower(y),lower(y+1)
        year_records.append({'year':y,'start_record':lo,'end_record':hi,'raw_records':hi-lo})
        accum=[]
        for start in range(lo,hi,250000):
            block=data[start:min(start+250000,hi)]
            a={x:np.asarray(block[fields[x]]) for x in needed}
            keep=(a['sample']==y*100+1)&np.isin(a['gq'],(1,2))&(a['pernum']==1)&(a['relate']==1)&(a['hhwt']>0)&(a['age']>=18)&(a['age']<=85)&np.isin(a['ownershp'],(1,2))&(a['rooms']>0)&np.isin(a['met2013'],METROS)
            if not keep.any(): continue
            h=pd.DataFrame({x:a[x][keep] for x in needed})
            accum.append(h)
        heads=pd.concat(accum,ignore_index=True) if accum else pd.DataFrame(columns=needed)
        heads['year']=y
        all_heads.append(heads)
    h=pd.concat(all_heads,ignore_index=True)
    head_key_cols=['year','sample','serial']
    assert not h.duplicated(head_key_cols).any(), 'Head household keys must be unique.'
    # Full four-moment target sample exactly follows the canonical raw reader.
    h['w']=h.hhwt.astype(float); h['r9']=np.minimum(h.rooms.astype(float),9)
    base=h
    age_owner=base[base.age.between(30,55)&base.unitsstr.between(3,10)]
    fam=base[base.age.between(30,55)&(base.nchild>0)&(base.yngch<18)].copy()
    recent_yes=age_owner[(age_owner.nchild>0)&(age_owner.eldch<4)]
    recent_no=age_owner[age_owner.nchild==0]
    def wmean(d,x): return float((d.w*d[x]).sum()/d.w.sum()) if d.w.sum()>0 else None
    def share(d,flag): return float(d.loc[flag(d),'w'].sum()/d.w.sum()) if d.w.sum()>0 else None
    standard={
      'mean_rooms':wmean(base,'r9'),
      'ownership_30_55':share(age_owner,lambda d:d.ownershp==1),
      'family_rooms':wmean(fam[fam.nchild>=3],'r9')-wmean(fam[fam.nchild.between(1,2)],'r9'),
      'recent_parent_ownership':share(recent_yes,lambda d:d.ownershp==1)-share(recent_no,lambda d:d.ownershp==1),
    }

    # In the exact family_rooms estimation universe, link person records by
    # household ID and count children whose MOMLOC or POPLOC points to the
    # head. RELATE=3 is retained as a secondary category comparison.
    fam['hhkey']=list(zip(fam.year.astype(int),fam['sample'].astype(int),fam.serial.astype(int)))
    eligible=set(fam.hhkey)
    counts={k:{'minor':0,'adult':0,'all':0,'relate_minor':0,'relate_all':0} for k in eligible}
    for y in YEARS:
        fail_time(); lo,hi=lower(y),lower(y+1)
        for start in range(lo,hi,250000):
            block=data[start:min(start+250000,hi)]
            yy=np.asarray(block[fields['year']]); ss=np.asarray(block[fields['sample']]); ser=np.asarray(block[fields['serial']]); per=np.asarray(block[fields['pernum']]); rel=np.asarray(block[fields['relate']]); age=np.asarray(block[fields['age']]); mom=np.asarray(block[fields['momloc']]); pop=np.asarray(block[fields['poploc']])
            keep=(yy==y)&(ss==y*100+1)&(per>1)&np.isin(age,np.arange(0,100))&((mom==1)|(pop==1)|(rel==3))
            ix=np.flatnonzero(keep)
            for i in ix:
                k=(int(yy[i]),int(ss[i]),int(ser[i]))
                if k not in counts: continue
                
                if int(rel[i])==3:
                    counts[k]['relate_all']+=1
                    if int(age[i])<18: counts[k]['relate_minor']+=1
                if int(mom[i])==1 or int(pop[i])==1:
                    counts[k]['all']+=1
                    if int(age[i])<18: counts[k]['minor']+=1
                    else: counts[k]['adult']+=1
    fam=fam.copy()
    fam['resident_own_children_under18']=[counts[k]['minor'] for k in fam.hhkey]
    fam['resident_own_children_18plus']=[counts[k]['adult'] for k in fam.hhkey]
    fam['resident_own_children_all_ages']=[counts[k]['all'] for k in fam.hhkey]
    fam['relate3_children_under18']=[counts[k]['relate_minor'] for k in fam.hhkey]
    fam['relate3_children_all_ages']=[counts[k]['relate_all'] for k in fam.hhkey]
    # Verify parent-linked all-age counts against NCHILD categories capped at
    # 3. This checks bins, not exact counts above 3 because NCHILD is top-coded.
    # RELATE=3 agreement is a secondary diagnostic.
    fam['relation_nchild_bin']=np.minimum(fam.relate3_children_all_ages,3)
    fam['linked_nchild_bin']=np.minimum(fam.resident_own_children_all_ages,3)
    fam['source_nchild_bin']=np.minimum(fam.nchild.astype(int),3)
    concordance={str(i):int((fam.loc[fam.source_nchild_bin==i,'relation_nchild_bin']==i).sum()) for i in range(4)}
    denom={str(i):int((fam.source_nchild_bin==i).sum()) for i in range(4)}
    mismatch=fam[fam.relation_nchild_bin!=fam.source_nchild_bin]
    linked_mismatch_count=int((fam.linked_nchild_bin!=fam.source_nchild_bin).sum())
    assert linked_mismatch_count==0, 'Parent-linked child bins must match NCHILD capped at 3.'

    actual_hi=fam[fam.resident_own_children_under18>=3]
    actual_lo=fam[fam.resident_own_children_under18.between(1,2)]
    original_hi=fam[fam.nchild>=3]; original_lo=fam[fam.nchild.between(1,2)]
    def weighted_group(d):
        return {'households':int(len(d)),'household_weight':float(d.w.sum()),'mean_rooms_cap9':wmean(d,'r9')}
    existing_moment=weighted_group(original_hi)['mean_rooms_cap9']-weighted_group(original_lo)['mean_rooms_cap9']
    candidate_moment=weighted_group(actual_hi)['mean_rooms_cap9']-weighted_group(actual_lo)['mean_rooms_cap9']
    assert len(actual_hi)+len(actual_lo)==len(fam), 'Candidate groups must span the family_rooms sample.'
    adult_family=fam[fam.resident_own_children_18plus>0]
    target_gaps={k:standard[k]-v for k,v in {'mean_rooms':5.561097376118652,'ownership_30_55':0.6483340343191493,'family_rooms':0.3470669317962507,'recent_parent_ownership':0.16289550916123285}.items()}
    assert all(abs(gap)<1e-9 for gap in target_gaps.values()), f'Existing ACS target reproduction failed: {target_gaps}'
    result={
      'status':'computed',
      'period':'ACS 1-year 2005 and 2006 pooled, household-weighted',
      'sample':'Same 42 active MET2013 geographies and exact head age 30-55, positive HHWT, occupied household and ROOMS>0 filters; family universe additionally NCHILD>0 and YNGCH<18.',
      'relationship_measure':'Counts household person rows whose MOMLOC or POPLOC equals the head PERNUM (1), split AGE<18 versus AGE>=18. These are resident children with an observed ACS parent pointer to the head. A secondary RELATE=3 count is retained as a relationship-category comparison.',
      'existing_four_moment_reproduction':standard,
      'expected_existing_four_moments':{'mean_rooms':5.561097376118652,'ownership_30_55':0.6483340343191493,'family_rooms':0.3470669317962507,'recent_parent_ownership':0.16289550916123285},
      'target_reproduction_max_abs_gap':float(max(abs(standard[k]-v) for k,v in {'mean_rooms':5.561097376118652,'ownership_30_55':0.6483340343191493,'family_rooms':0.3470669317962507,'recent_parent_ownership':0.16289550916123285}.items())),
      'family_sample':{'households':int(len(fam)),'household_weight':float(fam.w.sum()),'families_with_adult_resident_own_child':int(len(adult_family)),'weighted_share_with_adult_resident_own_child':float(adult_family.w.sum()/fam.w.sum()),'families_with_any_minor_linked_child':int((fam.resident_own_children_under18>0).sum()),'weighted_share_with_zero_linked_minor_children':float(fam.loc[fam.resident_own_children_under18==0,'w'].sum()/fam.w.sum())},
      'existing_NCHILD_grouping':{'high_3plus':weighted_group(original_hi),'low_1to2':weighted_group(original_lo),'rooms_difference_high_minus_low':existing_moment},
      'under18_resident_child_grouping_candidate':{'high_3plus':weighted_group(actual_hi),'low_1to2':weighted_group(actual_lo),'rooms_difference_high_minus_low':candidate_moment,'households_excluded_from_two_groups':int(len(fam)-len(actual_hi)-len(actual_lo))},
      'paired_moment_shift_candidate_minus_existing':float(candidate_moment-existing_moment),
      'nchild_relationship_count_concordance':{'meaning':'Secondary check: RELATE=3 all-age resident children, capped at 3, compared with NCHILD categories capped at 3 in the family_rooms universe. Primary parent-linked counts also match NCHILD after this cap. This verifies bins only; exact counts above 3 are unavailable because NCHILD is top-coded.','denominator_by_NCHILD_bin':denom,'exact_bin_agreement_by_NCHILD_bin':concordance,'mismatch_households':int(len(mismatch)),'mismatch_weighted_share':float(mismatch.w.sum()/fam.w.sum()) if len(fam) else None,
        'parent_link_count_bin_agreement_households':int((fam.linked_nchild_bin==fam.source_nchild_bin).sum()),
        'parent_link_count_mismatch_households':linked_mismatch_count,
        'parent_link_count_mismatch_weighted_share':float(fam.loc[fam.linked_nchild_bin!=fam.source_nchild_bin,'w'].sum()/fam.w.sum())},
      'age_count_distribution_weighted':{str(k):float(fam.loc[fam.resident_own_children_under18==k,'w'].sum()/fam.w.sum()) for k in sorted(fam.resident_own_children_under18.unique())},
      'RELATE3_only_candidate_for_comparison':{'adult_child_weighted_share':float(fam.loc[fam.relate3_children_all_ages>fam.relate3_children_under18,'w'].sum()/fam.w.sum()),'rooms_difference_high_minus_low':float(wmean(fam[fam.relate3_children_under18>=3],'r9')-wmean(fam[fam.relate3_children_under18.between(1,2)],'r9'))},
      'year_record_ranges':year_records,
      'active_metros':METROS,
      'elapsed_seconds':time.monotonic()-START,
    }
    (OUT/'housing_profiles_v1_child_mapping_recomputed.json').write_text(json.dumps(result,indent=2)+'\n')
    group_rows=[]
    for label,group in [('existing_NCHILD_3plus',original_hi),('existing_NCHILD_1to2',original_lo),('under18_children_3plus',actual_hi),('under18_children_1to2',actual_lo)]:
      group_rows.append({'group':label,'households':int(len(group)),'household_weight':float(group.w.sum()),'mean_rooms_cap9':wmean(group,'r9')})
    pd.DataFrame(group_rows).to_csv(OUT/'family_rooms_group_aggregates.csv',index=False)
    target_rows=[]
    for name,expected in result['expected_existing_four_moments'].items():
      target_rows.append({'target':name,'expected_target':expected,'recomputed':standard[name],'gap':standard[name]-expected,'status':'reproduced' if abs(standard[name]-expected)<1e-9 else 'mismatch'})
    target_rows.append({'target':'family_rooms_under18_resident_children_candidate','expected_target':'','recomputed':candidate_moment,'gap':candidate_moment-existing_moment,'status':'candidate_only'})
    pd.DataFrame(target_rows).to_csv(OUT/'target_reproduction.csv',index=False)
    print(json.dumps({'status':'PASS','existing_gaps':{x['target']:x['gap'] for x in target_rows[:4]},'adult_family_weighted_share':result['family_sample']['weighted_share_with_adult_resident_own_child'],'old':existing_moment,'new':candidate_moment,'shift':result['paired_moment_shift_candidate_minus_existing'],'seconds':result['elapsed_seconds']}),flush=True)

if __name__=='__main__': main()
