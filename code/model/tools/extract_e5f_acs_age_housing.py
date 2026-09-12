"""Extract active-metro ACS housing allocation by household-head age."""
from pathlib import Path
import json, sys, importlib.util
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
OUT = ROOT / "output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/age_housing_allocation"
YEARS = (2005, 2006)

def main():
    spec = importlib.util.spec_from_file_location("early", ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    class CompatReader(mod.HeaderReader):
        def __init__(self, path):
            mod.pd.io.stata.StataParser.__init__(self)
            self.path_or_buf = path.open('rb'); self._path_or_buf = self.path_or_buf
            self.col_sizes=[]; self._convert_dates=False; self._convert_categoricals=False
            self._index_col=None; self._convert_missing=False; self._preserve_dtypes=True
            self._order_categoricals=True; self._encoding=''; self._chunksize=1; self._using_iterator=False
            self._has_string_data=False; self._missing_values=False; self._can_read_value_labels=False
            self._column_selector_set=False; self._value_labels_read=False; self._data_read=False
            self._dtype=None; self._lines_read=0; self._native_byteorder='<' if sys.byteorder=='little' else '>'
            self._read_header(); self._setup_dtype()
    reader = CompatReader(SRC)
    dtype = reader._dtype
    fields = dict(zip(reader.varlist if hasattr(reader,"varlist") else reader._varlist, dtype.names))
    data = np.memmap(SRC, dtype=dtype, mode="r", offset=(reader.data_location if hasattr(reader,"data_location") else reader._data_location), shape=((reader.nobs if hasattr(reader,"nobs") else reader._nobs),))
    nobs = (reader.nobs if hasattr(reader,"nobs") else reader._nobs)
    ids = set(map(int, (ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/active_metros.txt").read_text().split(",")))
    needed = ("year", "sample", "met2013", "gq", "pernum", "relate", "hhwt", "age", "ownershp", "rooms")
    def lower(y):
        lo, hi = 0, nobs
        while lo < hi:
            mid = (lo + hi) // 2
            if int(data[fields["year"]][mid]) < y: lo = mid + 1
            else: hi = mid
        return lo
    rows=[]
    for year in YEARS:
        lo, hi = lower(year), lower(year+1)
        for start in range(lo, hi, 250000):
            block = data[start:min(start+250000, hi)]
            a = {x: np.asarray(block[fields[x]]) for x in needed}
            keep = ((a["sample"] == year*100+1) & np.isin(a["met2013"], list(ids)) & np.isin(a["gq"], (1,2)) &
                    (a["pernum"] == 1) & (a["relate"] == 1) & (a["hhwt"] > 0) & (a["age"] >= 18) & (a["age"] <= 85) &
                    np.isin(a["ownershp"], (1,2)) & (a["rooms"] > 0))
            if not keep.any(): continue
            age = a["age"][keep].astype(int); w = a["hhwt"][keep].astype(float); rooms = a["rooms"][keep].astype(float)
            bins = 18 + 4 * ((age - 18) // 4)
            owner = a["ownershp"][keep] == 1
            for b in np.unique(bins):
                q = bins == b
                rows.append(dict(year=year, age_lower=int(b), age_upper=min(int(b)+3,85), n=int(q.sum()), hhwt=float(w[q].sum()), rooms_capped9_weighted_sum=float((w[q]*np.minimum(rooms[q],9)).sum()), owner_hhwt=float(w[q & owner].sum())))
    frame = pd.DataFrame(rows).groupby(["year","age_lower","age_upper"], as_index=False).sum().sort_values(["year","age_lower"])
    frame["age"] = frame["age_lower"]
    frame["households"] = frame["hhwt"]
    frame["capped_rooms"] = frame["rooms_capped9_weighted_sum"]
    OUT.mkdir(parents=True, exist_ok=True)
    comp = pd.read_csv(ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/early_housing_metro_components.csv")
    c = comp[(comp.active_metro==1) & comp.year.isin(YEARS)].sum(numeric_only=True)
    checks = {"n": bool(int(frame.n.sum()) == int(c.n)), "hhwt": bool(abs(frame.hhwt.sum()-c.weight)<1e-6), "rooms_capped9": bool(abs(frame.rooms_capped9_weighted_sum.sum()-c.capped_rooms_sum)<1e-6), "owner_hhwt": bool(abs(frame.owner_hhwt.sum()-c.owner_weight)<1e-6)}
    payload={"status":"PASS" if all(checks.values()) else "FAIL", "years":list(YEARS), "age_bins":"18-21,...,82-85", "weight":"HHWT", "unit":"household heads (PERNUM=1, RELATE=1)", "checks":checks, "extracted_totals":{"n":int(frame.n.sum()),"hhwt":float(frame.hhwt.sum()),"rooms_capped9":float(frame.rooms_capped9_weighted_sum.sum()),"owner_hhwt":float(frame.owner_hhwt.sum())}, "reference_totals":{"n":int(c.n),"hhwt":float(c.weight),"rooms_capped9":float(c.capped_rooms_sum),"owner_hhwt":float(c.owner_weight)}}
    payload.update(source=str(SRC),source_size=SRC.stat().st_size,source_mtime_ns=SRC.stat().st_mtime_ns,active_metros=sorted(ids))
    assert all(checks.values()), payload
    frame.to_csv(OUT/"data_age_housing_by_year.csv", index=False)
    pooled=frame.drop(columns="year").groupby(["age","age_lower","age_upper"],as_index=False).sum()
    assert len(pooled)==17
    pooled.to_csv(OUT/"data_age_housing.csv",index=False)
    (OUT/"data_verification.json").write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps(payload))
if __name__ == "__main__": main()
