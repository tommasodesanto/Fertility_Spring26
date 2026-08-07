#!/usr/bin/env python3
"""Stage the small national entry samples from fixed record blocks of extract27.

Usage: run one non-overlapping block at a time, then run
``summarize_staged_entry.py``.  Direct record seeking avoids repeatedly
reading the 9.2-GB source extract when execution time is bounded.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

WORKDIR = Path(__file__).resolve().parent
DATA = WORKDIR.parents[2] / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
OUTDIR = WORKDIR / "staged_blocks"
COLUMNS = ["year", "sample", "statefip", "gq", "pernum", "relate", "age", "hhwt", "perwt", "citizen", "yrimmig"]

parser = argparse.ArgumentParser()
parser.add_argument("--start", type=int, required=True)
parser.add_argument("--count", type=int, required=True)
args = parser.parse_args()

OUTDIR.mkdir(exist_ok=True)
reader = pd.read_stata(DATA, columns=COLUMNS, iterator=True, chunksize=250_000, convert_categoricals=False)
if args.start < 0 or args.start >= reader.nobs:
    raise ValueError(f"start must be in [0, {reader.nobs})")
reader._lines_read = args.start  # pandas StataReader: seek directly to this row.
remaining = min(args.count, reader.nobs - args.start)
frames = []
diag = {"pernum1_relate1": 0, "pernum1_not_relate1": 0, "relate1_not_pernum1": 0}
while remaining:
    chunk = reader.read(min(250_000, remaining))
    remaining -= len(chunk)
    chunk = chunk.loc[chunk.year.between(2012, 2023) & chunk.statefip.between(1, 56)].copy()
    if chunk.empty:
        continue
    valid_hhwt = chunk.hhwt.notna() & (chunk.hhwt > 0)
    household = chunk.gq.isin([1, 2, 5])
    p1, r1 = chunk.pernum.eq(1), chunk.relate.eq(1)
    diag["pernum1_relate1"] += int((p1 & r1 & household & valid_hhwt).sum())
    diag["pernum1_not_relate1"] += int((p1 & ~r1 & household & valid_hhwt).sum())
    diag["relate1_not_pernum1"] += int((r1 & ~p1 & household & valid_hhwt).sum())
    keep_head = p1 & r1 & household & valid_hhwt & chunk.age.between(18, 25)
    keep_person = chunk.age.between(18, 22) & chunk.perwt.notna() & (chunk.perwt > 0)
    keep = keep_head | keep_person
    frames.append(chunk.loc[keep])
reader.close()
out = OUTDIR / f"block_{args.start:08d}_{args.count:08d}.csv"
pd.concat(frames, ignore_index=True).to_csv(out, index=False)
(OUTDIR / f"block_{args.start:08d}_{args.count:08d}.diag.csv").write_text(
    "metric,count\n" + "\n".join(f"{k},{v}" for k, v in diag.items()) + "\n"
)
print(f"Wrote {out} ({sum(len(x) for x in frames):,} selected records); diag={diag}")
