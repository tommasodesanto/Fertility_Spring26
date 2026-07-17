"""
Merge ATTOM assessor-history parquet shards into a single Stata .dta file.

Source: 64 snappy parquet shards on Desktop (assessor-history_*.parquet),
sharing an identical 279-column schema, ~368,751 rows total.

The script:
  1. Reads all shards via pyarrow into a single Arrow table.
  2. Casts decimal128 columns to float64 and date32 to pandas datetime.
  3. Lowercases column names and truncates to Stata's 32-char limit.
     (Verified at write time: no collisions.)
  4. Saves merged.parquet for Python/DuckDB use and assessor_history.dta
     for Stata. Original column names are preserved as Stata variable labels.

Run from repo root:
    python code/data/attom_sample/merge_to_stata.py
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd

SRC_DIR = "/Users/tommasodesanto/Desktop/019e6b85-e6d0-79f3-ac21-a5dc416c8dfa"
OUT_DIR = Path(__file__).resolve().parent
PARQUET_OUT = OUT_DIR / "assessor_history.parquet"
DTA_OUT = OUT_DIR / "assessor_history.dta"


def stata_name(orig: str) -> str:
    """Lowercase + truncate to 32 chars (Stata limit)."""
    return orig.lower()[:32]


def main() -> None:
    files = sorted(glob.glob(os.path.join(SRC_DIR, "*.parquet")))
    if not files:
        sys.exit(f"No parquet files found in {SRC_DIR}")
    print(f"Found {len(files)} shards.")

    tables = []
    for i, fp in enumerate(files, 1):
        tables.append(pq.read_table(fp))
        if i % 16 == 0 or i == len(files):
            print(f"  read {i}/{len(files)}")
    table = pa.concat_tables(tables)
    print(f"Combined table: {table.num_rows:,} rows x {table.num_columns} cols")

    # Cast decimal -> float64 (Stata-friendly). Keep date32 as-is for pandas conversion.
    new_fields = []
    for f in table.schema:
        if pa.types.is_decimal(f.type):
            new_fields.append(pa.field(f.name, pa.float64()))
        else:
            new_fields.append(f)
    casts = []
    for f, nf in zip(table.schema, new_fields):
        col = table.column(f.name)
        if f.type != nf.type:
            col = col.cast(nf.type, safe=False)
        casts.append(col)
    table = pa.Table.from_arrays(casts, schema=pa.schema(new_fields))

    # Write merged parquet (for Python/DuckDB users).
    pq.write_table(table, PARQUET_OUT, compression="snappy")
    print(f"Wrote {PARQUET_OUT}  ({PARQUET_OUT.stat().st_size/1e6:.1f} MB)")

    # Convert to pandas.
    df = table.to_pandas(types_mapper=pd.ArrowDtype)
    # Stata write needs plain numpy dtypes for floats/strings/dates, not Arrow extension dtypes.
    # Easier path: re-convert with default (numpy-backed) dtypes.
    df = table.to_pandas()
    print(f"Pandas frame: {df.shape}, mem={df.memory_usage(deep=True).sum()/1e9:.2f} GB")

    # Build rename map (original -> stata-safe) and variable_labels (stata-safe -> original).
    orig_cols = list(df.columns)
    rename = {c: stata_name(c) for c in orig_cols}
    new_cols = list(rename.values())
    if len(set(new_cols)) != len(new_cols):
        sys.exit("Column-name collision after truncation; need a smarter renamer.")
    df = df.rename(columns=rename)

    # Stata labels are capped at 80 chars; original names are <=42, safe.
    variable_labels = {new: orig for orig, new in rename.items()}

    # pandas to_stata: integer-flagged decimal columns are already float64. Strings can be NA.
    # Use version=118 (Stata 14+, supports strL and UTF-8). convert_dates handles datetimes.
    date_cols = {c: "td" for c in df.columns if pd.api.types.is_datetime64_any_dtype(df[c])}

    # Object columns that are entirely null break to_stata; coerce to "" and stringify mixed.
    obj_cols = [c for c in df.columns if df[c].dtype == object]
    fixed_allnull = 0
    for c in obj_cols:
        s = df[c]
        if s.isna().all():
            df[c] = ""
            fixed_allnull += 1
        else:
            df[c] = s.where(s.notna(), "").astype(str)
    print(f"Coerced {fixed_allnull} all-null object cols to empty string; "
          f"stringified {len(obj_cols)} object cols total.")

    print(f"Writing {DTA_OUT} ...")
    df.to_stata(
        DTA_OUT,
        write_index=False,
        version=118,
        variable_labels=variable_labels,
        convert_dates=date_cols,
    )
    print(f"Wrote {DTA_OUT}  ({DTA_OUT.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
