"""Extract comparable national and 42-metro occupied-room stocks at model dates.

Run with pandas >=2.2 (seekable Stata reader), no whole-source copy or bootstrap.
Defaults write only /tmp/full_housing_stock_history.json; no target changes.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[5]
BASE = Path(__file__).resolve().parent
SOURCE = ROOT/'code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta'
YEARS = (2007, 2011, 2012, 2015, 2019, 2023)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('/tmp/full_housing_stock_history.json'))
    args=parser.parse_args();start_time=time.monotonic()
    ids=np.array([int(x) for x in (BASE/'active_metros.txt').read_text().strip().split(',')])
    assert len(set(ids))==42
    fields_needed=('year','sample','met2013','gq','pernum','relate','hhwt','age','ownershp','rooms')
    output=[];checks=[]
    old=list(csv.DictReader((BASE/'early_housing_metro_components.csv').open()))
    oldnational=list(csv.DictReader((ROOT/'code/data/Spatial_aggregate_withmicrodata/output/national_householder_housing_path/national_householder_housing_path.csv').open()))
    with pd.io.stata.StataReader(SOURCE,convert_dates=False,convert_categoricals=False) as reader:
        reader._ensure_open()
        assert hasattr(reader._path_or_buf,'fileno'), 'Reader must remain file backed'
        data=np.memmap(SOURCE,dtype=reader._dtype,mode='r',offset=reader._data_location,shape=(reader._nobs,))
        fields=dict(zip(reader._varlist,reader._dtype.names));assert all(x in fields for x in fields_needed)
        def lower(year):
            lo,hi=0,reader._nobs
            while lo<hi:
                mid=(lo+hi)//2
                if int(data[fields['year']][mid])<year:lo=mid+1
                else:hi=mid
            return lo
        for year in YEARS:
            lo,hi=lower(year),lower(year+1);assert 0<hi-lo<4000000
            totals={g:dict(records=0,households=0.,total_occupied_rooms_literal=0.,total_occupied_rooms_capped9=0.) for g in ('national','42_metros')}
            for index in range(lo,hi,250000):
                if time.monotonic()-start_time>480:raise TimeoutError('Eight-minute extraction limit')
                block=data[index:min(index+250000,hi)]
                a={x:np.asarray(block[fields[x]]) for x in fields_needed}
                assert np.all(a['year']==year)
                keep=(a['sample']==year*100+1)&np.isin(a['gq'],(1,2))&(a['pernum']==1)&(a['relate']==1)&(a['hhwt']>0)&(a['age']>=18)&(a['age']<=85)&np.isin(a['ownershp'],(1,2))&(a['rooms']>0)
                for group,mask in [('national',keep),('42_metros',keep&np.isin(a['met2013'],ids))]:
                    w=a['hhwt'][mask].astype(float);rooms=a['rooms'][mask].astype(float);t=totals[group]
                    t['records']+=int(mask.sum());t['households']+=float(w.sum());t['total_occupied_rooms_literal']+=float((w*rooms).sum());t['total_occupied_rooms_capped9']+=float((w*np.minimum(rooms,9)).sum())
            for group,t in totals.items():
                t.update(calendar_year=year,geography=group,mean_occupied_rooms_capped9=t['total_occupied_rooms_capped9']/t['households'])
                output.append(t)
            ref=[r for r in old if int(r['year'])==year and r['active_metro']=='1']
            if ref:
                for key,column in [('records','n'),('households','weight'),('total_occupied_rooms_capped9','capped_rooms_sum')]:
                    gap=totals['42_metros'][key]-sum(float(r[column]) for r in ref);assert abs(gap)<1e-9,(year,key,gap)
                    checks.append(dict(year=year,geography='42_metros',field=key,absolute_gap=abs(gap)))
            ref=[r for r in oldnational if int(r['calendar_year'])==year]
            if ref:
                for key,column in [('records','base_head_records'),('households','base_head_weight'),('total_occupied_rooms_literal','rooms_weighted_sum')]:
                    gap=totals['national'][key]-float(ref[0][column]);assert abs(gap)<1e-9,(year,key,gap)
                    checks.append(dict(year=year,geography='national',field=key,absolute_gap=abs(gap)))
            print(json.dumps({'year':year,'elapsed_seconds':round(time.monotonic()-start_time,2),'totals':totals}),flush=True)
    for group in ('national','42_metros'):
        base=next(r for r in output if r['calendar_year']==2007 and r['geography']==group)
        for r in output:
            if r['geography']==group:
                r['total_rooms_index_2007_100']=100*r['total_occupied_rooms_capped9']/base['total_occupied_rooms_capped9']
                r['households_index_2007_100']=100*r['households']/base['households']
    receipt=json.loads((BASE/'early_housing_source_receipt.json').read_text())
    assert SOURCE.stat().st_size==receipt['source_size'] and SOURCE.stat().st_mtime_ns==receipt['source_mtime_ns']
    result={'status':'PASS','builder':str(Path(__file__).resolve()),'builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source':str(SOURCE),'source_size':SOURCE.stat().st_size,'source_mtime_ns':SOURCE.stat().st_mtime_ns,'source_sha256_from_existing_receipt':receipt['source_sha256_from_existing_canonical_receipt'],'source_rehashed_this_pass':False,'sample':'ACS1-year SAMPLE=year*100+1; GQ1or2; PERNUM1; RELATE1; HHWT>0; ages18–85; OWNERSHP1or2; ROOMS>0. HHWT already scaled in source.42-metro subset uses fixed MET2013 membership, not admitted PUMA.','room_definition':'Physical rooms occupied; min(ROOMS,9) per household before HHWT aggregation. Covers occupied renter/owner housing, excludes vacancies. Comparable cap9 across2007 and later coding vintages.','metro_ids':[int(x) for x in ids],'series':output,'verification':checks,'elapsed_seconds':time.monotonic()-start_time,'interpretation':'Total stock is weighted sum of occupied capped rooms, not housing-unit count or mean rooms. National household weights are comparable geography to nationally conditioned model heads; historical national household-count growth is an externally conditioned component, not fully predicted.42metro series retains initial calibration geography. No target or model changes.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(str(args.output),flush=True)
if __name__=='__main__':main()
