#!/usr/bin/env python3
"""Bounded ACS lifecycle housing diagnostic for the frozen 2005--06 contract.

This is a diagnostic extractor, not a target builder.  It reads only sorted
2005/2006 records from the authoritative ``extract27.dta`` with the existing
file-backed HeaderReader, and writes separate active-42-metro and national
tables.  The national table is descriptive and is never used for the scalar
target gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
HOUSING = ROOT / "output/model/e5f_matched_pf_20260909a/design_research/housing"
SOURCE = ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
METROS = HOUSING / "active_metros.txt"
TARGETS = ROOT / "output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_weights.csv"
SOURCE_RECEIPT = HOUSING / "early_housing_source_receipt.json"
REUSED_READER = HOUSING / "inspect_early_housing.py"
YEARS = (2005, 2006)
CHUNK = 250_000
TIME_LIMIT = 300.0


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _reader_module():
    spec = importlib.util.spec_from_file_location("authoritative_early_housing_reader", REUSED_READER)
    if spec is None or spec.loader is None:
        raise ImportError(REUSED_READER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _open_sorted_source():
    module = _reader_module()
    reader = module.HeaderReader(SOURCE)
    dtype = reader._dtype
    fields = dict(zip(reader.varlist, dtype.names))
    data = np.memmap(SOURCE, dtype=dtype, mode="r", offset=reader.data_location, shape=(reader.nobs,))
    return data, fields, reader.nobs


def _lower_bound(data: np.ndarray, field: str, value: int) -> int:
    lo, hi = 0, len(data)
    while lo < hi:
        mid = (lo + hi) // 2
        if int(data[field][mid]) < value:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _empty_accumulator() -> dict[tuple[Any, ...], dict[str, float]]:
    return {}


def _add_group(acc: dict, key: tuple, frame: pd.DataFrame) -> None:
    if frame.empty:
        return
    w = frame["hhwt"].to_numpy(float)
    rooms = frame["rooms"].to_numpy(float)
    owner = frame["ownershp"].to_numpy(int) == 1
    row = acc.setdefault(key, {"n_records": 0, "hhwt": 0.0, "rooms_sum": 0.0,
                               "rooms_capped9_sum": 0.0, "owner_hhwt": 0.0})
    row["n_records"] += int(len(frame))
    row["hhwt"] += float(w.sum())
    row["rooms_sum"] += float((w * rooms).sum())
    row["rooms_capped9_sum"] += float((w * np.minimum(rooms, 9.0)).sum())
    row["owner_hhwt"] += float(w[owner].sum())


def _add_scalar(stats: dict, frame: pd.DataFrame) -> None:
    if frame.empty:
        return
    w = frame["hhwt"].to_numpy(float)
    rooms = frame["rooms"].to_numpy(float)
    owner = frame["ownershp"].to_numpy(int) == 1
    stats["hhwt"] += float(w.sum())
    stats["rooms_capped9_sum"] += float((w * np.minimum(rooms, 9.0)).sum())
    stats["owner_hhwt"] += float(w[owner].sum())


def _aggregate_chunk(frame: pd.DataFrame, active_metros: set[int], acc: dict,
                     scalar_stats: dict[str, dict[str, float]] | None = None) -> None:
    """Aggregate one already-filtered household-head chunk (fixture-testable)."""
    frame = frame.copy()
    frame["child_bin"] = np.select(
        [frame["nchild"].isna(), frame["nchild"].le(0), frame["nchild"].eq(1),
         frame["nchild"].eq(2), frame["nchild"].gt(2)],
        ["unknown", "0", "1", "2", "3+"], default="unknown"
    )
    frame["young_child"] = np.where(
        frame["nchild"].gt(0) & frame["yngch"].ne(99) & frame["yngch"].lt(18), "with_young_child", "without_young_child"
    )
    frame["recent_parent"] = frame["nchild"].gt(0) & frame["eldch"].ne(99) & frame["eldch"].lt(4)
    frame["tenure"] = np.where(frame["ownershp"].eq(1), "owner", "renter")
    frame["due"] = frame["unitsstr"].between(3, 10)
    frame["age_bin"] = 18 + 4 * ((frame["age"].astype(int) - 18) // 4)
    frame["active42"] = frame["met2013"].isin(active_metros)
    if scalar_stats is not None:
        for scope, scope_mask in (("national", np.ones(len(frame), dtype=bool)),
                                  ("active42", frame["active42"].to_numpy(bool))):
            s = frame.loc[scope_mask]
            _add_scalar(scalar_stats[scope]["mean_rooms"], s)
            due = s.loc[s["due"] & s["age"].between(30, 55)]
            _add_scalar(scalar_stats[scope]["ownership"], due)
            _add_scalar(scalar_stats[scope]["newparent"], due.loc[due["recent_parent"]])
            _add_scalar(scalar_stats[scope]["nochild"], due.loc[due["nchild"].eq(0)])
            family = s.loc[s["age"].between(30, 55) & s["young_child"] .eq("with_young_child")]
            _add_scalar(scalar_stats[scope]["family_large"], family.loc[family["nchild"].ge(3)])
            _add_scalar(scalar_stats[scope]["family_small"], family.loc[family["nchild"].between(1, 2)])
    for scope, scope_mask in (("national", np.ones(len(frame), dtype=bool)),
                              ("active42", frame["active42"].to_numpy(bool))):
        sub_scope = frame.loc[scope_mask]
        for sample, sample_mask in (("all_structures", np.ones(len(sub_scope), dtype=bool)),
                                    ("DUE", sub_scope["due"].to_numpy(bool))):
            sub = sub_scope.loc[sample_mask]
            for age_kind, age_col in (("annual", "age"), ("four_year", "age_bin")):
                for age_value, age_group in sub.groupby(age_col, sort=True):
                    for tenure, tenure_group in age_group.groupby("tenure", sort=True):
                        for child_bin, child_group in tenure_group.groupby("child_bin", sort=True):
                            for young, cell in child_group.groupby("young_child", sort=True):
                                key = (scope, sample, age_kind, int(age_value), tenure, child_bin, young)
                                _add_group(acc, key, cell)


def _target_values() -> dict[str, float]:
    wanted = {"mean_rooms": None, "ownership_30_55": None, "family_rooms": None, "recent_parent_ownership": None}
    with TARGETS.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["restriction_id"] in wanted:
                wanted[row["restriction_id"]] = float(row["target"])
    if any(value is None for value in wanted.values()):
        raise ValueError("working_weights.csv is missing one or more active housing rows")
    return wanted


def _target_from_stats(stats: dict[str, dict[str, float]]) -> dict[str, float]:
    def own(name): return stats[name]["owner_hhwt"] / stats[name]["hhwt"]
    def rooms(name): return stats[name]["rooms_capped9_sum"] / stats[name]["hhwt"]
    return {"mean_rooms": rooms("mean_rooms"),
            "ownership_30_55": own("ownership"),
            "recent_parent_ownership": own("newparent") - own("nochild"),
            "family_rooms": rooms("family_large") - rooms("family_small")}


def _write_outputs(acc: dict, out: Path, metadata: dict) -> None:
    rows = []
    for key, values in sorted(acc.items(), key=lambda item: tuple(map(str, item[0]))):
        scope, sample, age_kind, age_value, tenure, child_bin, young = key
        rows.append({"geography_scope": scope, "sample": sample, "age_kind": age_kind,
                     "age_lower": age_value, "age_upper": age_value if age_kind == "annual" else min(age_value + 3, 85),
                     "tenure": tenure, "current_children": child_bin, "young_child": young, **values})
    pd.DataFrame(rows).to_csv(out / "housing_profile_by_age.csv", index=False)
    (out / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")


def run(output: Path, *, smoke: bool = False) -> dict:
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory must be new and empty: {output}")
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    active_metros = {int(x) for x in METROS.read_text().strip().split(",")}
    data, fields, nobs = _open_sorted_source()
    required = ("year", "sample", "met2013", "gq", "pernum", "relate", "hhwt", "age", "ownershp", "rooms", "unitsstr", "nchild", "yngch", "eldch")
    missing = [name for name in required if name not in fields]
    if missing: raise KeyError(f"missing source fields: {missing}")
    acc = _empty_accumulator(); progress = []
    scalar_names = ("mean_rooms", "ownership", "newparent", "nochild", "family_large", "family_small")
    scalar_stats = {scope: {name: {"hhwt": 0.0, "rooms_capped9_sum": 0.0, "owner_hhwt": 0.0} for name in scalar_names} for scope in ("active42", "national")}
    years = (2005,) if smoke else YEARS
    for year in years:
        lo, hi = _lower_bound(data, fields["year"], year), _lower_bound(data, fields["year"], year + 1)
        stop = min(lo + CHUNK, hi) if smoke else hi
        for start in range(lo, stop, CHUNK):
            if time.monotonic() - started > TIME_LIMIT: raise TimeoutError("five-minute budget exceeded")
            block = data[start:min(start + CHUNK, stop)]
            arr = {name: np.asarray(block[fields[name]]) for name in required}
            keep = ((arr["sample"] == year * 100 + 1) & np.isin(arr["gq"], (1, 2)) & (arr["pernum"] == 1) &
                    (arr["relate"] == 1) & (arr["hhwt"] > 0) & (arr["age"] >= 18) & (arr["age"] <= 85) &
                    np.isin(arr["ownershp"], (1, 2)) & (arr["rooms"] > 0))
            frame = pd.DataFrame({name: arr[name][keep] for name in required})
            _aggregate_chunk(frame, active_metros, acc, scalar_stats)
            progress.append({"year": year, "chunk_start": int(start), "chunk_end": int(min(start + CHUNK, stop)), "records_kept": int(keep.sum()), "elapsed_seconds": time.monotonic() - started})
            (output / "progress.json").write_text(json.dumps(progress, indent=2) + "\n")
    target_values = None if smoke else _target_values()
    target_recomputed = None if smoke else {scope: _target_from_stats(scalar_stats[scope]) for scope in ("active42", "national")}
    gate = None if smoke else {key: {"target": target_values[key], "recomputed": target_recomputed["active42"][key], "pass": bool(abs(target_values[key] - target_recomputed["active42"][key]) <= 1e-10)} for key in target_values}
    if gate and not all(row["pass"] for row in gate.values()): raise AssertionError(f"active42 target gate failed: {gate}")
    receipt = json.loads(SOURCE_RECEIPT.read_text())
    receipt_match = (SOURCE.stat().st_size == receipt.get("source_size") and SOURCE.stat().st_mtime_ns == receipt.get("source_mtime_ns"))
    reused_hash = receipt.get("source_sha256_from_existing_canonical_receipt") if receipt_match else None
    metadata = {"status": "smoke_exact_loop_only" if smoke else "full_pass", "smoke": smoke, "years": list(years), "active_metros": sorted(active_metros), "source": str(SOURCE), "source_size": SOURCE.stat().st_size, "source_mtime_ns": SOURCE.stat().st_mtime_ns, "source_sha256_from_existing_canonical_receipt": reused_hash, "source_hash_reuse_size_mtime_match": receipt_match, "source_rehashed_this_pass": False, "reused_reader": str(REUSED_READER), "reader_sha256": _sha(REUSED_READER), "metro_file_sha256": _sha(METROS), "target_file_sha256": _sha(TARGETS), "target_gate": gate, "chunk_size": CHUNK, "time_limit_seconds": TIME_LIMIT, "progress_file": str(output / "progress.json"), "elapsed_seconds": time.monotonic() - started, "definitions": {"rooms": "positive literal ROOMS; capped at 9 before aggregation", "current_children": "NCHILD bins 0, 1, 2, 3+, unknown preserved", "young_child": "NCHILD>0 and YNGCH<18, with YNGCH=99 in without-young-child denominator", "recent_parent": "NCHILD>0 and ELDCH<4, ELDCH=99 excluded", "DUE": "UNITSSTR in 3:10 for ownership/recent-parent only", "family_rooms": "age 30-55, YNGCH<18, NCHILD>=3 versus NCHILD 1-2, no DUE restriction", "weights": "HHWT", "model_comparison": "empirical current children, not model CEB"}}
    _write_outputs(acc, output, metadata)
    (output / "target_recomputed.json").write_text(json.dumps({"target_values": target_values, "recomputed": target_recomputed, "gate": gate}, indent=2) + "\n")
    (output / "final_receipt.json").write_text(json.dumps({"status": metadata["status"], "rows": sum(1 for _ in acc), "elapsed_seconds": metadata["elapsed_seconds"], "target_gate": gate, "progress_entries": len(progress)}, indent=2) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    run(args.output, smoke=args.smoke)


if __name__ == "__main__":
    main()
