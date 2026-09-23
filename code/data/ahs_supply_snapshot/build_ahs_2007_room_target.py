#!/usr/bin/env python3
"""Reproduce the 2007 national AHS occupied-room target for ages 18--85.

The estimator is the housing-unit-weighted mean of literal PUF ROOMS, with
STATUS=1, positive WGT90GEO and ROOMS, and HHAGE in [18, 85]. The 2007 PUF
topcodes ROOMS at 21; no tail correction is imposed. Fay-BRR uncertainty uses
REPWGT1--REPWGT160 and the Census/HUD factor 4/160.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import urllib.request
import zipfile
from pathlib import Path


SOURCE_URL = (
    "https://www2.census.gov/programs-surveys/ahs/2007/"
    "AHS%202007%20National%20PUF%20v2.0%20Flat%20CSV.zip"
)
SOURCE_SHA256 = "413324add01c65d178091162f5cee4fd11bd82ca1fbbe0194d3c92a727be598f"
ROOT = Path(__file__).resolve().parent
DEFAULT_ZIP = ROOT / "raw/ahs_2007_national_puf_v2_flat_csv.zip"
DEFAULT_OUTPUT = ROOT / "output/ahs_2007_room_target.json"
REPLICATES = 160


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def literal_number(value: str) -> float:
    return float(value.strip().strip("'"))


def calculate(path: Path) -> dict[str, object]:
    totals = [0.0] * (REPLICATES + 1)
    rooms_totals = [0.0] * (REPLICATES + 1)
    capped_rooms_total = top_tail_weight = topcode_weight = 0.0
    n_occupied = n_age = n_valid = 0
    max_rooms = 0
    with zipfile.ZipFile(path) as archive:
        if archive.namelist() != ["ahs2007n.csv"]:
            raise ValueError("Unexpected AHS flat-file archive member")
        with archive.open("ahs2007n.csv") as raw:
            reader = csv.reader(io.TextIOWrapper(raw, encoding="utf-8-sig", newline=""))
            header = next(reader)
            names = ("STATUS", "HHAGE", "ROOMS", "WGT90GEO", "REPWGT0")
            names += tuple(f"REPWGT{i}" for i in range(1, REPLICATES + 1))
            missing = set(names) - set(header)
            if missing:
                raise ValueError(f"Missing AHS fields: {sorted(missing)}")
            ix = {name: header.index(name) for name in names}
            for row in reader:
                if row[ix["STATUS"]].strip("'") != "1":
                    continue
                n_occupied += 1
                age = int(literal_number(row[ix["HHAGE"]]))
                if not 18 <= age <= 85:
                    continue
                n_age += 1
                rooms = int(literal_number(row[ix["ROOMS"]]))
                weight = literal_number(row[ix["WGT90GEO"]])
                if rooms <= 0 or weight <= 0:
                    continue
                if rooms > 21 or not math.isclose(
                    weight, literal_number(row[ix["REPWGT0"]]), rel_tol=0, abs_tol=1e-7
                ):
                    raise ValueError("Unexpected rooms topcode or full-sample replicate weight")
                n_valid += 1
                max_rooms = max(max_rooms, rooms)
                totals[0] += weight
                rooms_totals[0] += weight * rooms
                capped_rooms_total += weight * min(rooms, 9)
                top_tail_weight += weight * (rooms > 9)
                topcode_weight += weight * (rooms == 21)
                for i in range(1, REPLICATES + 1):
                    replicate_weight = literal_number(row[ix[f"REPWGT{i}"]])
                    if replicate_weight < 0:
                        raise ValueError(f"Negative replicate weight {i}")
                    totals[i] += replicate_weight
                    rooms_totals[i] += replicate_weight * rooms
    if n_valid == 0 or any(total <= 0 for total in totals):
        raise ValueError("No eligible AHS records or empty replicate")
    estimate = rooms_totals[0] / totals[0]
    replicate_means = [rooms_totals[i] / totals[i] for i in range(1, REPLICATES + 1)]
    variance = 4.0 / REPLICATES * sum((value - estimate) ** 2 for value in replicate_means)
    return {
        "source_url": SOURCE_URL,
        "source_sha256": SOURCE_SHA256,
        "survey": "2007 AHS National PUF v2.0 flat CSV",
        "sample": "STATUS=1, HHAGE 18--85, WGT90GEO>0, ROOMS>0",
        "weight": "WGT90GEO",
        "room_definition": "Literal ROOMS; public-use topcode 21, no tail imputation",
        "occupied_records": n_occupied,
        "age_eligible_records": n_age,
        "analysis_records": n_valid,
        "weighted_households": totals[0],
        "weighted_rooms": rooms_totals[0],
        "mean_rooms": estimate,
        "standard_error": math.sqrt(variance),
        "variance_method": "Fay BRR: (4/160) * sum((replicate_mean - full_mean)^2)",
        "replicate_count": REPLICATES,
        "mean_rooms_if_capped_at_9": capped_rooms_total / totals[0],
        "share_above_9": top_tail_weight / totals[0],
        "share_at_puf_topcode_21": topcode_weight / totals[0],
        "max_observed_rooms": max_rooms,
        "status": "empirical_target_only; numerical_calibration_contract_unchanged",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not args.source_zip.exists():
        args.source_zip.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(SOURCE_URL, args.source_zip)
    observed_hash = sha256(args.source_zip)
    if observed_hash != SOURCE_SHA256:
        raise ValueError(f"AHS source SHA-256 mismatch: {observed_hash}")
    result = calculate(args.source_zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"AHS 2007 rooms: {result['mean_rooms']:.12f}; SE {result['standard_error']:.12f}")
    print(f"Receipt: {args.output}")


if __name__ == "__main__":
    main()
