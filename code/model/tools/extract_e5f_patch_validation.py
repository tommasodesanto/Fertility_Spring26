"""Extract the actual 2023 ACS household-head lifecycle housing levels.

The source is the sorted IPUMS ACS ``extract27.dta`` file.  This deliberately
reads only the 2023 year range through a header-compatible memmap reader, so a
patch validation cannot silently scan the full nine-gigabyte file.  ``with_minor``
counts all valid household heads with a resident minor child under the same NCHILD/YNGCH rule
used by ``extract_e5f_large_owner_acs.py``; YNGCH=99 remains in the denominator
as having no resident minor.
"""
from pathlib import Path
import importlib.util
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
DEFAULT_OUT = ROOT / (
    "output/model/e5f_matched_pf_20260909a/current_candidate_transition/"
    "overnight_20260912/patch_readout/data"
)
YEAR = 2023
CHUNK = 250_000


def _reader(path):
    spec = importlib.util.spec_from_file_location(
        "early_header", ROOT / "output/model/e5f_matched_pf_20260909a/design_research/"
        "housing/inspect_early_housing.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    class CompatReader(mod.HeaderReader):
        def __init__(self, source):
            mod.pd.io.stata.StataParser.__init__(self)
            self.path_or_buf = source.open("rb")
            self._path_or_buf = self.path_or_buf
            self.col_sizes = []
            self._convert_dates = False
            self._convert_categoricals = False
            self._index_col = None
            self._convert_missing = False
            self._preserve_dtypes = True
            self._order_categoricals = True
            self._encoding = ""
            self._chunksize = 1
            self._using_iterator = False
            self._has_string_data = False
            self._missing_values = False
            self._can_read_value_labels = False
            self._column_selector_set = False
            self._value_labels_read = False
            self._data_read = False
            self._dtype = None
            self._lines_read = 0
            self._native_byteorder = "<" if sys.byteorder == "little" else ">"
            self._read_header()
            self._setup_dtype()

    reader = CompatReader(path)
    dtype = reader._dtype
    names = reader.varlist if hasattr(reader, "varlist") else reader._varlist
    fields = dict(zip(names, dtype.names))
    nobs = reader.nobs if hasattr(reader, "nobs") else reader._nobs
    offset = reader.data_location if hasattr(reader, "data_location") else reader._data_location
    data = np.memmap(path, dtype=dtype, mode="r", offset=offset, shape=(nobs,))
    return data, fields, nobs


def _lower_bound(data, field, value):
    lo, hi = 0, len(data)
    while lo < hi:
        mid = (lo + hi) // 2
        if int(data[field][mid]) < value:
            lo = mid + 1
        else:
            hi = mid
    return lo


def extract(source=SRC):
    data, f, nobs = _reader(source)
    metro_file = ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing/active_metros.txt"
    metros = set(map(int, metro_file.read_text().split(",")))
    needed = ("year", "sample", "met2013", "gq", "pernum", "relate", "hhwt", "age",
              "ownershp", "rooms", "nchild", "yngch", "sex", "fertyr", "perwt")
    for field in needed:
        if field not in f:
            raise KeyError(f"extract27.dta is missing required field {field!r}")

    start = _lower_bound(data, f["year"], YEAR)
    stop = _lower_bound(data, f["year"], YEAR + 1)
    totals = {lo: {"households": 0, "hhwt": 0.0, "capped_rooms": 0.0,
                   "owners": 0.0, "with_minor": 0.0, "owner_with_minor": 0.0}
              for lo in range(18, 86, 4)}
    invalid_room_records = 0
    invalid_room_hhwt = 0.0
    fert = {lo: {"women": 0.0, "recent_birth": 0.0, "birth_event": 0.0, "birth_event_records": 0}
            for lo in range(18, 50, 4)}
    for pos in range(start, stop, CHUNK):
        block = data[pos:min(pos + CHUNK, stop)]
        a = {name: np.asarray(block[f[name]]) for name in needed}
        keep = ((a["sample"] == YEAR * 100 + 1) & np.isin(a["met2013"], list(metros)) &
                np.isin(a["gq"], (1, 2)) & (a["pernum"] == 1) & (a["relate"] == 1) &
                (a["hhwt"] > 0) & (a["age"] >= 18) & (a["age"] <= 85) &
                np.isin(a["ownershp"], (1, 2)) & (a["rooms"] > 0))
        female = ((a["sample"] == YEAR*100+1) & np.isin(a["met2013"],list(metros)) & np.isin(a["gq"],(1,2)) &
                  (a["sex"]==2) & (a["perwt"]>0) & np.isin(a["fertyr"],(1,2)) & (a["age"]>=18) & (a["age"]<=49))
        bins=18+4*((a["age"].astype(int)-18)//4)
        for lo in np.unique(bins[female]):
            q=female & (bins==lo);event=q & (a["fertyr"]==2)
            fert[int(lo)]["women"]+=float(a["perwt"][q].sum())
            fert[int(lo)]["recent_birth"]+=float(a["perwt"][event].sum())
            fert[int(lo)]["birth_event"]+=float(a["perwt"][event].sum())
            fert[int(lo)]["birth_event_records"]+=int(event.sum())
        if not keep.any():
            continue
        age = a["age"][keep].astype(int)
        weight = a["hhwt"][keep].astype(float)
        rooms = a["rooms"][keep].astype(float)
        owner = a["ownershp"][keep] == 1
        minor = (a["nchild"][keep] > 0) & (a["yngch"][keep] < 18)
        invalid_room_records += int((rooms >= 99).sum())
        invalid_room_hhwt += float(weight[rooms >= 99].sum())
        for lo in np.unique(18 + 4 * ((age - 18) // 4)):
            q = (18 + 4 * ((age - 18) // 4)) == lo
            totals[int(lo)]["households"] += int(q.sum())
            totals[int(lo)]["hhwt"] += float(weight[q].sum())
            totals[int(lo)]["capped_rooms"] += float((weight[q] * np.minimum(rooms[q], 9)).sum())
            totals[int(lo)]["owners"] += float(weight[q & owner].sum())
            totals[int(lo)]["with_minor"] += float(weight[q & minor].sum())
            totals[int(lo)]["owner_with_minor"] += float(weight[q & owner & minor].sum())

    rows = []
    for lo in range(18, 86, 4):
        hi = min(lo + 3, 85)
        row = {"age_lower": lo, "age_upper": hi, **totals[lo]}
        row["ownership_rate"] = row["owners"] / row["hhwt"]
        row["with_minor_rate"] = row["with_minor"] / row["hhwt"]
        row["mean_capped_rooms"] = row["capped_rooms"] / row["hhwt"]
        rows.append(row)
    frame = pd.DataFrame(rows)
    pooled = {"age_lower": 18, "age_upper": 85}
    for col in ("households", "hhwt", "capped_rooms", "owners", "with_minor", "owner_with_minor"):
        pooled[col] = frame[col].sum()
    pooled["ownership_rate"] = pooled["owners"] / pooled["hhwt"]
    pooled["with_minor_rate"] = pooled["with_minor"] / pooled["hhwt"]
    pooled["owner_with_minor_rate"] = pooled["owner_with_minor"] / pooled["hhwt"]
    pooled["mean_capped_rooms"] = pooled["capped_rooms"] / pooled["hhwt"]
    fertility = pd.DataFrame([{
        "age_lower": lo, "age_upper": min(lo + 3, 49), **fert[lo],
        "recent_birth_rate": fert[lo]["recent_birth"] / fert[lo]["women"],
        "birth_event_rate": fert[lo]["birth_event"] / fert[lo]["women"],
    } for lo in range(18, 50, 4)])
    return frame, pd.DataFrame([pooled]), fertility, {
        "source": str(source), "source_size": source.stat().st_size,
        "source_mtime_ns": source.stat().st_mtime_ns, "year": YEAR,
        "year_row_range": [start, stop], "nobs_source": nobs, "active_metros": sorted(metros),
        "sample": 202301, "weight": "HHWT", "unit": "active matched metros; GQ 1/2; PERNUM=1; RELATE=1 household heads",
        "room_measure": "HHWT-weighted physical ROOMS capped at 9; positive rooms only",
        "ownership_measure": "OWNERSHP=1 among valid household heads",
        "minor_measure": "all valid household heads with NCHILD>0 and YNGCH<18; YNGCH=99 remains denominator as no resident minor",
        "owner_with_minor_measure": "optional separate subset: OWNERSHP=1 and NCHILD>0 and YNGCH<18",
        "room_code_audit": {"rule": "canonical positive ROOMS rule retained; values >=99 are flagged, not silently dropped", "records_ge99": invalid_room_records, "hhwt_ge99": invalid_room_hhwt},
        "aggregate_validation": "pooled ages 18-85 and separately 22-85, sums of four-year bins",
        "fertility_2023_source": {
            "available": True,
            "path": str(source),
            "units": "ACS person records, source-year 2023; SEX identifies female records, AGE is completed age, FERTYR is birth in the last 12 months, and PERWT is the person weight",
            "available_fields": ["SEX", "AGE", "FERTYR", "PERWT", "HHWT", "NCHILD", "YNGCH", "MET2013", "OWNERSHP"],
            "flow_measure_note": "Existing local builders describe FERTYR as births in the last 12 months and construct age-specific rates from female exposure; this patch remeasures the female recent-birth fraction from extract27.dta.",
            "patch_measure": "All women ages18-49, same42metros/GQ1-2/ACS2023; FERTYR in(1,2), eventFERTYR2, PERWT denominator; no household-head, ownership or room restrictions",
            "model_equivalence": "not established by this housing extraction; no model-equivalent female exposure or age-specific flow is invented here",
            "related_existing_aggregate": str(ROOT / "code/data/Spatial_aggregate_withmicrodata/calibration_targets_output/asfr_by_location.csv"),
        },
    }


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output if args.output.is_absolute() else ROOT / args.output
    frame, pooled, fertility, metadata = extract()
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "actual2023_age_housing_levels.csv", index=False)
    pooled.to_csv(out / "actual2023_aggregate_validation.csv", index=False)
    fertility.to_csv(out / "actual2023_female_recent_birth_rates.csv", index=False)
    frame[frame.age_lower >= 22].assign(age_lower=22, age_upper=85).groupby(["age_lower", "age_upper"], as_index=False).sum().assign(
        ownership_rate=lambda x: x.owners / x.hhwt,
        with_minor_rate=lambda x: x.with_minor / x.hhwt,
        owner_with_minor_rate=lambda x: x.owner_with_minor / x.hhwt,
        mean_capped_rooms=lambda x: x.capped_rooms / x.hhwt,
    ).to_csv(out / "actual2023_aggregate_22_85.csv", index=False)
    (out / "actual2023_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(frame.to_string(index=False))
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
