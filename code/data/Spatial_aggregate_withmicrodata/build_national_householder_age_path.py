#!/usr/bin/env python3
"""Build national four-year ACS householder-age shares for model transitions.

The raw IPUMS Stata file is sorted by year. This builder uses the file's
fixed-width metadata to seek directly to each requested ACS year, so it does
not materialize the 9 GB extract in memory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_EXTRACT = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
)
DEFAULT_OUTDIR = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/output/national_householder_age_path"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, default=DEFAULT_EXTRACT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--years", default="2007,2011,2015,2019,2023")
    parser.add_argument("--age-start", type=int, default=18)
    parser.add_argument("--period-years", type=int, default=4)
    parser.add_argument("--age-cells", type=int, default=17)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Refusing to write an empty age path")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    extract = args.extract.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    years = tuple(int(item) for item in args.years.split(",") if item.strip())
    if not years:
        raise ValueError("At least one ACS year is required")
    age_end_exclusive = args.age_start + args.period_years * args.age_cells

    reader = pd.read_stata(
        extract,
        iterator=True,
        convert_categoricals=False,
        convert_dates=False,
    )
    try:
        dtype = getattr(reader, "_dtype", None)
        data_location = getattr(reader, "data_location", None)
        raw_names = getattr(dtype, "names", None)
        if dtype is None or data_location is None or raw_names is None:
            raise RuntimeError("Unsupported pandas StataReader metadata")
        field = dict(zip(reader.varlist, raw_names))
        required = {"year", "sample", "gq", "pernum", "relate", "hhwt", "age"}
        missing = sorted(required - set(field))
        if missing:
            raise ValueError(f"ACS extract is missing fields: {missing}")
        record_size = int(dtype.itemsize)
        stream = reader.path_or_buf

        def read_records(start: int, count: int) -> np.ndarray:
            stream.seek(int(data_location) + int(start) * record_size)
            raw = np.frombuffer(
                stream.read(int(count) * record_size),
                dtype=dtype,
                count=int(count),
            )
            if reader.byteorder != reader._native_byteorder:
                raw = raw.byteswap().newbyteorder()
            return raw

        def lower_year_bound(target_year: int) -> int:
            lower, upper = 0, int(reader.nobs)
            while lower < upper:
                midpoint = (lower + upper) // 2
                observed = int(read_records(midpoint, 1)[field["year"]][0])
                if observed < target_year:
                    lower = midpoint + 1
                else:
                    upper = midpoint
            return lower

        rows: list[dict[str, Any]] = []
        receipts: dict[str, Any] = {}
        for year in years:
            first = lower_year_bound(year)
            last = lower_year_bound(year + 1)
            if first >= last:
                raise ValueError(f"ACS extract has no rows for {year}")
            weights = np.zeros(args.age_cells, dtype=float)
            records = np.zeros(args.age_cells, dtype=np.int64)
            all_head_weight = 0.0
            all_head_records = 0
            sample_codes: set[int] = set()
            for start in range(first, last, 250_000):
                raw = read_records(start, min(250_000, last - start))
                sample_codes.update(
                    int(value) for value in np.unique(raw[field["sample"]])
                )
                head = (
                    (raw[field["year"]] == year)
                    & np.isin(raw[field["gq"]], [1, 2])
                    & (raw[field["pernum"]] == 1)
                    & (raw[field["relate"]] == 1)
                    & np.isfinite(raw[field["hhwt"]])
                    & (raw[field["hhwt"]] > 0.0)
                )
                all_head_records += int(np.sum(head))
                all_head_weight += float(np.sum(raw[field["hhwt"]][head]))
                in_range = (
                    head
                    & (raw[field["age"]] >= args.age_start)
                    & (raw[field["age"]] < age_end_exclusive)
                )
                ages = raw[field["age"]][in_range].astype(int)
                hhwt = raw[field["hhwt"]][in_range].astype(float)
                bins = (ages - args.age_start) // args.period_years
                weights += np.bincount(
                    bins, weights=hhwt, minlength=args.age_cells
                )
                records += np.bincount(bins, minlength=args.age_cells)
            in_range_weight = float(np.sum(weights))
            if in_range_weight <= 0.0 or np.any(weights <= 0.0):
                raise RuntimeError(f"Nonpositive age-bin mass in {year}")
            shares = weights / in_range_weight
            if not np.isclose(float(np.sum(shares)), 1.0, atol=1e-14, rtol=0.0):
                raise RuntimeError(f"Age shares fail to sum to one in {year}")
            for index, (record_count, weight, share) in enumerate(
                zip(records, weights, shares, strict=True)
            ):
                age_lower = args.age_start + args.period_years * index
                rows.append(
                    {
                        "year": year,
                        "age_lower": age_lower,
                        "age_upper": age_lower + args.period_years - 1,
                        "unweighted_head_records": int(record_count),
                        "hhwt_head_mass": float(weight),
                        "share_conditional_age_range": float(share),
                    }
                )
            receipts[str(year)] = {
                "raw_year_records": int(last - first),
                "sample_codes": sorted(sample_codes),
                "head_records_all_ages": all_head_records,
                "head_weight_all_ages": all_head_weight,
                "head_weight_age_range": in_range_weight,
                "age_range_share_of_all_heads": in_range_weight / all_head_weight,
            }
            print(
                f"ACS_HEAD_AGE year={year} all_heads={all_head_weight:.0f} "
                f"age_coverage={in_range_weight / all_head_weight:.6f}",
                flush=True,
            )
    finally:
        reader.close()

    write_csv(outdir / "national_householder_age_path.csv", rows)
    write_json(
        outdir / "metadata.json",
        {
            "status": "complete",
            "source": str(extract),
            "source_sha256": sha256(extract),
            "sample": (
                "national IPUMS ACS household heads; GQ in {1,2}; "
                "PERNUM=1; RELATE=1; HHWT>0"
            ),
            "weight": "HHWT",
            "years": years,
            "age_start": args.age_start,
            "age_end": age_end_exclusive - 1,
            "period_years": args.period_years,
            "age_cells": args.age_cells,
            "receipts": receipts,
        },
    )
    print(f"ACS_HEAD_AGE_COMPLETE outdir={outdir}", flush=True)


if __name__ == "__main__":
    main()
