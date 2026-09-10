"""Bounded empirical diagnostic; no active target edits or model calls.

Uses pandas' DTA header parser with an open file, avoiding pandas 1.5.3's
constructor copy of the entire 9.2 GiB source into BytesIO. Fixed-record,
year-sorted seeking follows the existing national housing builder.
"""
from pathlib import Path
import csv
import json
import sys
import time
import hashlib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta'
YEARS = (2005, 2006, 2007, 2012, 2023)
START = time.monotonic()

class HeaderReader(pd.io.stata.StataReader):
    def __init__(self, path):
        pd.io.stata.StataParser.__init__(self)
        self.col_sizes=[]
        self._convert_dates=False
        self._convert_categoricals=False
        self._index_col=None
        self._convert_missing=False
        self._preserve_dtypes=True
        self._columns=None
        self._order_categoricals=True
        self._encoding=''
        self._chunksize=1
        self._using_iterator=False
        self._has_string_data=False
        self._missing_values=False
        self._can_read_value_labels=False
        self._column_selector_set=False
        self._value_labels_read=False
        self._data_read=False
        self._dtype=None
        self._lines_read=0
        self._native_byteorder='<' if sys.byteorder=='little' else '>'
        self.path_or_buf=path.open('rb')
        self._read_header()
        self._setup_dtype()

def main():
    r=HeaderReader(SOURCE)
    dtype=r._dtype
    fields=dict(zip(r.varlist,dtype.names))
    assert r.byteorder==r._native_byteorder
    data=np.memmap(SOURCE,dtype=dtype,mode='r',offset=r.data_location,shape=(r.nobs,))
    def lower(y):
        lo,hi=0,r.nobs
        while lo<hi:
            mid=(lo+hi)//2
            if int(data[fields['year']][mid])<y:lo=mid+1
            else:hi=mid
        return lo
    ids=[int(x) for x in (OUT/'active_metros.txt').read_text().strip().split(',')]
    assert len(ids)==42
    lookup={}
    for vintage in (2010,2020):
        tab=pd.read_csv(ROOT/f'code/data/mms_center_periphery/data/puma_mms_lookup_{vintage}.csv')
        lookup[vintage]=set(zip(tab.statefip.astype(int),tab.puma.astype(int),tab.cbsacode.astype(int)))
    required=('year','sample','statefip','puma','met2013','gq','pernum','relate','hhwt','age','ownershp','rooms','unitsstr','nchild','yngch','eldch')
    assert all(x in fields for x in required)
    rows=[]
    year_meta=[]
    for y in YEARS:
        lo,hi=lower(y),lower(y+1)
        assert hi-lo<=4000000
        if time.monotonic()-START>300:raise TimeoutError('Five-minute empirical bound exceeded')
        accum={}
        samples=set()
        for start in range(lo,hi,250000):
            block=data[start:min(start+250000,hi)]
            a={x:np.asarray(block[fields[x]]) for x in required}
            assert np.all(a['year']==y)
            samples.update(int(x) for x in np.unique(a['sample']))
            keep=(a['sample']==y*100+1)&np.isin(a['gq'],(1,2))&(a['pernum']==1)&(a['relate']==1)&(a['hhwt']>0)&(a['age']>=18)&(a['age']<=85)&np.isin(a['ownershp'],(1,2))&(a['rooms']>0)
            h=pd.DataFrame({x:a[x][keep].copy() for x in required})
            if h.empty:continue
            # Stata extract contains already-scaled HHWT; validate national receipt below.
            if y in (2012,2023):
                accepted=lookup[2010 if y==2012 else 2020]
                h['active_puma']=[int((int(s),int(p),int(m)) in accepted) for s,p,m in zip(h.statefip,h.puma,h.met2013)]
            else:h['active_puma']=-1
            for (metro,admitted),g in h.groupby(['met2013','active_puma'],sort=False):
                key=(int(metro),int(admitted))
                v=accum.setdefault(key, {'n':0,'weight':0.,'rooms_sum':0.,'capped_rooms_sum':0.,'owner_weight':0.,'rooms9_weight':0.,'rooms_gt9_weight':0.})
                w=g.hhwt.to_numpy(float);rooms=g.rooms.to_numpy(float)
                v['n']+=len(g);v['weight']+=w.sum();v['rooms_sum']+=(w*rooms).sum();v['capped_rooms_sum']+=(w*np.minimum(rooms,9)).sum();v['owner_weight']+=w[g.ownershp.to_numpy()==1].sum();v['rooms9_weight']+=w[rooms==9].sum();v['rooms_gt9_weight']+=w[rooms>9].sum()
                due=(g.age>=30)&(g.age<=55)&g.unitsstr.between(3,10)
                groups={'own3055':due,'own2534':g.age.between(25,34)&g.unitsstr.between(3,10),'newparent':due&(g.nchild>0)&(g.eldch<4),'nochild':due&(g.nchild==0),'rooms_large':g.age.between(30,55)&(g.nchild>=3)&(g.yngch<18),'rooms_small':g.age.between(30,55)&g.nchild.between(1,2)&(g.yngch<18)}
                for name,mask in groups.items():
                    b=mask.to_numpy();gg=g.loc[mask];ww=w[b];rr=rooms[b]
                    for suffix,value in [('n',len(gg)),('weight',ww.sum()),('owner_weight',ww[gg.ownershp.to_numpy()==1].sum()),('rooms_sum',(ww*rr).sum()),('capped_rooms_sum',(ww*np.minimum(rr,9)).sum())]:
                        col=f'{name}_{suffix}';v[col]=v.get(col,0)+float(value)
        for (metro,admitted),values in accum.items():rows.append({'year':y,'met2013':metro,'active_metro':int(metro in ids),'active_puma':admitted,**values})
        year_meta.append({'year':y,'start_record':lo,'end_record':hi,'raw_records':hi-lo,'sample_codes':sorted(samples),'metro_count_nonzero':len({m for m,a in accum if m>0})})
        print(json.dumps(year_meta[-1]),flush=True)
    frame=pd.DataFrame(rows).fillna(0)
    # Independently match existing national raw-data receipts for selected overlap years.
    canonical=pd.read_csv(ROOT/'code/data/Spatial_aggregate_withmicrodata/output/national_householder_housing_path/national_householder_housing_path.csv').set_index('calendar_year')
    replay=[]
    for y in (2007,2023):
        d=frame[frame.year==y].sum(numeric_only=True);old=canonical.loc[y]
        checks={'head_records':int(d['n'])==int(old['base_head_records']),'head_weight':abs(d['weight']-old['base_head_weight'])<1e-6,'rooms_sum':abs(d['rooms_sum']-old['rooms_weighted_sum'])<1e-6,'owner_weight':abs(d['owner_weight']-old['owner_weight'])<1e-6}
        assert all(checks.values()), (y,checks)
        replay.append({'year':y,**{k:bool(v) for k,v in checks.items()}})
    frame.to_csv(OUT/'early_housing_metro_components.csv',index=False)
    metadata={'purpose':'Diagnostic common-metro feasibility; not new calibration targets','source':str(SOURCE),'source_size':SOURCE.stat().st_size,'source_mtime_ns':SOURCE.stat().st_mtime_ns,'source_hash_from_existing_canonical_receipt':'edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e','raw_source_rehashed_this_pass':False,'reader':'File-backed pandas1.5.3 header parser plus numpy memmap; no whole-source copy','years':year_meta,'active_metros':ids,'canonical_overlap_replay':replay,'elapsed_seconds':time.monotonic()-START,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'early_housing_source_receipt.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('PASS',len(frame),'metro/year/admission records; seconds',time.monotonic()-START,flush=True)

if __name__=='__main__':main()
