#!/usr/bin/env python3
"""Produce national outside-entry measures from staged ACS/IPUMS block files."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
import pandas as pd

WORKDIR = Path(__file__).resolve().parent
STAGED = WORKDIR / "staged_blocks"
OUT = WORKDIR / "national_outside_entry_measures.csv"
parts = sorted(STAGED.glob("block_????????_????????.csv"))
if not parts:
    raise RuntimeError("No staged blocks found. Run stage_acs_blocks.py first.")
diag = pd.concat([pd.read_csv(p) for p in sorted(STAGED.glob("*.diag.csv"))]).groupby("metric")["count"].sum()
if diag.get("pernum1_not_relate1", 0) or diag.get("relate1_not_pernum1", 0):
    raise RuntimeError(f"Ambiguous head identification: {diag.to_dict()}")
d = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
if set(d.year.unique()) != set(range(2012, 2024)):
    raise RuntimeError(f"Years staged are not exactly 2012--2023: {sorted(d.year.unique())}")
# IPUMS CITIZEN: code 1 is born abroad of American parents; codes 2--5 are foreign-born.
d["foreign"] = d.citizen.isin([2, 3, 4, 5])
d["years_since_arrival"] = d.year - d.yrimmig
d["recent4"] = d.foreign & d.years_since_arrival.between(0, 4)
d["recent8"] = d.foreign & d.years_since_arrival.between(0, 8)
rows=[]
defs = [("foreign_born_recent_arrival_4y_head_hhwt", "recent4"), ("foreign_born_recent_arrival_8y_head_hhwt", "recent8"), ("foreign_born_stock_head_hhwt", "foreign")]
for age_label, lo, hi in [("18_22",18,22),("18_25",18,25)]:
    x=d[(d.pernum==1)&(d.relate==1)&d.gq.isin([1,2,5])&d.hhwt.gt(0)&d.age.between(lo,hi)]
    for measure, indicator in defs:
        for period, y in [("pooled",x)]+[(str(year),x[x.year==year]) for year in range(2012,2024)]:
            share=(y.hhwt*y[indicator]).sum()/y.hhwt.sum()
            rows.append(dict(period=period,age_band=age_label,measure=measure,unweighted_n=len(y),weighted_denominator=y.hhwt.sum(),weighted_numerator=(y.hhwt*y[indicator]).sum(),share=share,head_sample_below_5000=len(y)<5000))
x=d[d.perwt.gt(0)&d.age.between(18,22)]
for period,y in [("pooled",x)]+[(str(year),x[x.year==year]) for year in range(2012,2024)]:
    share=(y.perwt*y.recent4).sum()/y.perwt.sum()
    rows.append(dict(period=period,age_band="18_22",measure="foreign_born_recent_arrival_4y_person_perwt",unweighted_n=len(y),weighted_denominator=y.perwt.sum(),weighted_numerator=(y.perwt*y.recent4).sum(),share=share,head_sample_below_5000=False))
o=pd.DataFrame(rows).sort_values(["measure","age_band","period"], key=lambda s: s.ne("pooled") if s.name=="period" else s)
if not o.share.between(0,1, inclusive="neither").all(): raise RuntimeError("A share is outside (0,1).")
o.to_csv(OUT,index=False,float_format="%.10f")
print(f"Wrote {OUT}; head diagnostic={diag.to_dict()}")
