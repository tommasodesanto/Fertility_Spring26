"""Extract six-cell large-owner-home shares for the early ACS sample."""
from pathlib import Path
import importlib.util, json, sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
OUT = ROOT / "output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/age_housing_allocation"
YEARS = (2005, 2006)
AGE_BINS = ((22, 39), (40, 59), (60, 85))

def main():
    spec = importlib.util.spec_from_file_location("early", ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    class CompatReader(mod.HeaderReader):
        def __init__(self, path):
            mod.pd.io.stata.StataParser.__init__(self); self.path_or_buf = path.open('rb'); self._path_or_buf = self.path_or_buf
            self.col_sizes=[]; self._convert_dates=False; self._convert_categoricals=False; self._index_col=None; self._convert_missing=False; self._preserve_dtypes=True; self._order_categoricals=True; self._encoding=''; self._chunksize=1; self._using_iterator=False; self._has_string_data=False; self._missing_values=False; self._can_read_value_labels=False; self._column_selector_set=False; self._value_labels_read=False; self._data_read=False; self._dtype=None; self._lines_read=0; self._native_byteorder='<' if sys.byteorder=='little' else '>'; self._read_header(); self._setup_dtype()
    reader = CompatReader(SRC); dtype = reader._dtype
    names = reader.varlist if hasattr(reader, "varlist") else reader._varlist
    fields = dict(zip(names, dtype.names)); nobs = reader.nobs if hasattr(reader, "nobs") else reader._nobs
    offset = reader.data_location if hasattr(reader, "data_location") else reader._data_location
    data = np.memmap(SRC, dtype=dtype, mode="r", offset=offset, shape=(nobs,))
    metros = set(map(int, (ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/active_metros.txt").read_text().split(",")))
    needed = ("year", "sample", "met2013", "gq", "pernum", "relate", "hhwt", "age", "ownershp", "rooms", "nchild", "yngch")
    totals = {(lo, hi, minor): {"hhwt": 0.0, "records": 0} for lo, hi in AGE_BINS for minor in (False, True)}
    for year in YEARS:
        lo, hi = 0, nobs
        while lo < hi:
            mid = (lo + hi) // 2
            if int(data[fields["year"]][mid]) < year: lo = mid + 1
            else: hi = mid
        start, end = lo, 0
        lo, hi = 0, nobs
        while lo < hi:
            mid = (lo + hi) // 2
            if int(data[fields["year"]][mid]) < year + 1: lo = mid + 1
            else: hi = mid
        end = lo
        for start in range(start, end, 250000):
            block = data[start:min(start + 250000, end)]
            a = {x: np.asarray(block[fields[x]]) for x in needed}
            keep = ((a["sample"] == year * 100 + 1) & np.isin(a["met2013"], list(metros)) & np.isin(a["gq"], (1, 2)) & (a["pernum"] == 1) & (a["relate"] == 1) & (a["hhwt"] > 0) & (a["age"] >= 22) & (a["age"] <= 85) & (a["ownershp"] == 1) & (a["rooms"] >= 6) & (a["nchild"] >= 0) & (a["rooms"] < 99) & (a["nchild"] <= 9))
            for age_lo, age_hi in AGE_BINS:
                q = keep & (a["age"] >= age_lo) & (a["age"] <= age_hi)
                minor = q & (a["nchild"] > 0) & (a["yngch"] < 18)
                nonminor = q & ~((a["nchild"] > 0) & (a["yngch"] < 18))
                for flag, mask in ((False, nonminor), (True, minor)):
                    totals[(age_lo, age_hi, flag)]["hhwt"] += float(a["hhwt"][mask].sum()); totals[(age_lo, age_hi, flag)]["records"] += int(mask.sum())
    denom = sum(v["hhwt"] for v in totals.values())
    assert denom>0
    rows = []
    for lo, hi in AGE_BINS:
        for minor in (False, True):
            v = totals[(lo, hi, minor)]; rows.append({"age_group": f"{lo}-{hi}", "children_under_18": minor, "hhwt": v["hhwt"], "records": v["records"], "share": v["hhwt"] / denom})
    assert abs(sum(r["share"] for r in rows)-1)<1e-12
    OUT.mkdir(parents=True, exist_ok=True); pd.DataFrame(rows).to_csv(OUT / "large_owner_data.csv", index=False)
    prov = {"status": "PASS", "years": list(YEARS), "metros": sorted(metros), "weight": "HHWT", "unit": "PERNUM=1, RELATE=1 household heads", "filters": "sample year*100+1; GQ 1/2; age 22-85; OWNERSHP=1; literal 6<=ROOMS<99; observed NCHILD 0-9", "minor_definition": "nchild > 0 and yngch < 18, YNGCH=99 stays in the denominator as no resident minor", "denominator": "all six HHWT-weighted cells", "may_difference": "May used PERWT and birth-year generation cohorts; this reproduction uses HHWT and actual-age groups 22-39, 40-59, 60-85 for current target comparability", "total_hhwt": denom, "total_records": sum(v["records"] for v in totals.values()), "source": str(SRC), "source_size": SRC.stat().st_size, "source_mtime_ns": SRC.stat().st_mtime_ns}
    (OUT / "large_owner_data_provenance.json").write_text(json.dumps(prov, indent=2) + "\n"); print(pd.DataFrame(rows).to_string(index=False))
if __name__ == "__main__": main()
