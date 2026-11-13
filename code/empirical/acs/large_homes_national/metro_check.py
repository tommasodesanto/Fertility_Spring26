"""ACS 2023, all US households: who holds the owner-occupied homes with 6+ rooms.
Same definitions as code/model/tools/extract_e5f_large_owner_acs.py, without the metro restriction.
Also reports the share of all households (heads 22-85) with no own child under 18 at home."""
from pathlib import Path
import json, sys
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
OUT = ROOT / "output/empirical_acs_large_homes_national_20261112/metro_check"
YEAR = int(sys.argv[1]) if len(sys.argv) > 1 else 2023
MET = list(map(int, (ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/active_metros.txt").read_text().split(",")))
AGE_BINS = ((22, 39), (40, 59), (60, 85))
r = pd.io.stata.StataReader(str(SRC)); r._ensure_open(); r._setup_dtype()
dtype = r._dtype; names = r._varlist; fields = dict(zip(names, dtype.names))
data = np.memmap(SRC, dtype=dtype, mode="r", offset=r._data_location, shape=(r._nobs,))
yr = data[fields["year"]]
def bound(y):
    lo, hi = 0, len(data)
    while lo < hi:
        m = (lo + hi) // 2
        if int(yr[m]) < y: lo = m + 1
        else: hi = m
    return lo
s, e = bound(YEAR), bound(YEAR + 1)
need = ("sample", "gq", "pernum", "relate", "hhwt", "age", "ownershp", "rooms", "nchild", "yngch")
cells = {(lo, hi, k): [0.0, 0] for lo, hi in AGE_BINS for k in (False, True)}
allhh = {"w": 0.0, "w_nochild": 0.0, "n": 0}
for b in range(s, e, 250000):
    blk = data[b:min(b + 250000, e)]
    a = {x: np.asarray(blk[fields[x]]).astype(float) for x in need}
    base = np.isin(np.asarray(blk[fields["met2013"]]), MET) & (a["sample"] == YEAR * 100 + 1) & np.isin(a["gq"], (1, 2)) & (a["pernum"] == 1) & (a["relate"] == 1) & (a["hhwt"] > 0) & (a["age"] >= 22) & (a["age"] <= 85) & (a["nchild"] >= 0) & (a["nchild"] <= 9)
    minor = (a["nchild"] > 0) & (a["yngch"] < 18)
    allhh["w"] += a["hhwt"][base].sum(); allhh["w_nochild"] += a["hhwt"][base & ~minor].sum(); allhh["n"] += int(base.sum())
    keep = base & (a["ownershp"] == 1) & (a["rooms"] >= 6) & (a["rooms"] < 99)
    for lo, hi in AGE_BINS:
        q = keep & (a["age"] >= lo) & (a["age"] <= hi)
        for k, m in ((False, q & ~minor), (True, q & minor)):
            cells[(lo, hi, k)][0] += a["hhwt"][m].sum(); cells[(lo, hi, k)][1] += int(m.sum())
den = sum(v[0] for v in cells.values())
rows = [{"age_group": f"{lo}-{hi}", "children_under_18": k, "hhwt": cells[(lo, hi, k)][0], "records": cells[(lo, hi, k)][1], "share": cells[(lo, hi, k)][0] / den} for lo, hi in AGE_BINS for k in (False, True)]
OUT.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUT / f"large_owner_national_{YEAR}.csv", index=False)
summ = {"year": YEAR, "geography": "all US, no MET2013 restriction", "large_owner_records": sum(v[1] for v in cells.values()),
        "large_owner_no_child_share": sum(v[0] for (lo, hi, k), v in cells.items() if not k) / den,
        "all_households_no_child_share": allhh["w_nochild"] / allhh["w"], "all_households_records": allhh["n"],
        "filters": "sample year*100+1; GQ 1/2; PERNUM=1 RELATE=1 heads 22-85; HHWT>0; large owner: OWNERSHP=1, 6<=ROOMS<99; child = NCHILD>0 & YNGCH<18"}
(OUT / f"large_owner_national_{YEAR}_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
print(pd.DataFrame(rows).to_string(index=False)); print(json.dumps(summ, indent=2))
