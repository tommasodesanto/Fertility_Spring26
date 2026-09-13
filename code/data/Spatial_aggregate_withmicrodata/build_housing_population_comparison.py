"""Same-household ACS room and resident counts; file-backed, no model solves.

Run with /opt/anaconda3/bin/python (pandas StataReader with file-backed header).
People use PERWT; households and rooms use HHWT, on the same resident roster.
HHWT times roster size is retained as a separate weighting sensitivity.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'
SOURCE = ROOT/'code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta'


def main():
    start = time.monotonic()
    prior = json.loads((BASE/'source/historical_stock/housing_data.json').read_text())
    assert SOURCE.stat().st_size == prior['source_size']
    assert SOURCE.stat().st_mtime_ns == prior['source_mtime_ns']
    rows = []
    with pd.io.stata.StataReader(SOURCE, convert_dates=False, convert_categoricals=False) as reader:
        reader._ensure_open()
        assert hasattr(reader._path_or_buf, 'fileno')
        raw = np.memmap(SOURCE, dtype=reader._dtype, mode='r', offset=reader._data_location, shape=(reader._nobs,))
        fields = dict(zip(reader._varlist, reader._dtype.names))
        def lower(year):
            lo, hi = 0, reader._nobs
            while lo < hi:
                mid = (lo+hi)//2
                if raw[fields['year']][mid] < year:
                    lo = mid+1
                else:
                    hi = mid
            return lo
        for year in (2007, 2011, 2015, 2019, 2023):
            block = raw[lower(year):lower(year+1)]
            a = {k: np.array(block[fields[k]]) for k in
                 ('sample','serial','pernum','relate','gq','age','ownershp','hhwt','perwt','rooms')}
            assert np.all(a['sample'] == year*100+1)
            assert np.all(np.diff(a['serial']) >= 0)
            first = np.r_[0, np.flatnonzero(np.diff(a['serial']))+1]
            size = np.diff(np.r_[first,len(block)])
            # Check complete, ordered person rosters, not a head-only extract.
            expected = np.arange(len(block))-np.repeat(first,size)+1
            h = {k:v[first] for k,v in a.items()}
            keep = np.isin(h['gq'],(1,2)) & (h['pernum']==1) & (h['relate']==1) & (h['hhwt']>0) & (h['age']>=18) & (h['age']<=85) & np.isin(h['ownershp'],(1,2)) & (h['rooms']>0)
            person_keep = np.repeat(keep,size)
            assert np.array_equal(a['pernum'][person_keep],expected[person_keep]), 'Incomplete selected household roster'
            assert np.all(a['hhwt'][person_keep] == np.repeat(h['hhwt'],size)[person_keep])
            w = h['hhwt'][keep].astype(float)
            households = float(w.sum())
            rooms = float((w*np.minimum(h['rooms'][keep],9)).sum())
            persons = float((w*size[keep]).sum())
            perwt_persons = float(a['perwt'][person_keep].astype(float).sum())
            ref = next(r for r in prior['series'] if r['calendar_year']==year and r['geography']=='national')
            assert households == ref['households'] and rooms == ref['total_occupied_rooms_capped9']
            assert persons == float(a['hhwt'][person_keep].astype(float).sum())
            row = dict(year=year,households=households,occupied_rooms_capped9=rooms,
                       household_weighted_residents=persons,person_weighted_residents=perwt_persons,
                       persons_per_household=perwt_persons/households,rooms_per_household=rooms/households,
                       rooms_per_person=rooms/perwt_persons,
                       rooms_per_person_household_weighted=rooms/persons,
                       mean_household_size_household_weighted=persons/households,
                       sample_households=int(keep.sum()),sample_persons=int(person_keep.sum()))
            assert abs(row['rooms_per_person']*row['persons_per_household']-row['rooms_per_household']) < 1e-12
            rows.append(row)
            print(json.dumps(row),flush=True)
            if time.monotonic()-start>240:
                raise TimeoutError('Four-minute extraction limit')
    result = dict(status='PASS',rows=rows,source=str(SOURCE),source_sha256_from_prior_receipt=prior['source_sha256_from_existing_receipt'],
                  builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sample=prior['sample'],geography='national',
                  population_definition='All resident person records in the exact selected housing households, including children; excludes group quarters and households whose head falls outside ages18–85.',
                  weight_definition='Main count: PERWT sum over residents of selected households. Housing and household totals use HHWT. HHWT times complete household roster size is retained as a weighting sensitivity. This restricted-sample population is not total US population.',
                  sources=['https://usa.ipums.org/usa-action/variables/HHWT','https://usa.ipums.org/usa-action/variables/PERNUM','https://usa.ipums.org/usa-action/variables/PERWT'],
                  checks=['Complete consecutive person rosters in every selected household','HHWT constant within every selected household','All five prior household/room totals reproduce exactly','Person count independently sums from person records','All five room/person accounting identities'],
                  elapsed_seconds=time.monotonic()-start)
    out = BASE/'source/historical_stock/housing_population_data.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(out)


if __name__ == '__main__':
    main()
