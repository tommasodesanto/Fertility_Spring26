#!/usr/bin/env python3
"""Build national ACS housing holdouts and female exposure denominators.

The source IPUMS Stata extract is sorted by year.  This script uses the
fixed-width record metadata exposed by pandas to seek directly to the five
requested annual ACS samples rather than materializing the roughly 10 GB file.

The canonical contract is deliberately narrow: household heads aged 18--85,
weighted by HHWT; ownership is conditional on OWNERSHP in {1, 2}; and mean
rooms is conditional on positive literal ROOMS.  The latter applies no
top-code correction (in particular, the 2007 value 9 is used literally even
though its IPUMS label is "9+").
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_EXTRACT = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
)
DEFAULT_OUTDIR = (
    ROOT
    / "code/data/Spatial_aggregate_withmicrodata/output"
    / "national_householder_housing_path"
)
EXPECTED_SOURCE_SHA256 = (
    "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e"
)
YEARS = (2007, 2011, 2015, 2019, 2023)
EXPECTED_SAMPLE_CODES = {year: year * 100 + 1 for year in YEARS}
AGE_MIN = 18
AGE_MAX = 85
CHUNK_SIZE = 250_000
CSV_NAME = "national_householder_housing_path.csv"
FEMALE_EXPOSURE_CSV_NAME = "national_female_exposure_path.csv"
METADATA_NAME = "metadata.json"
VERIFICATION_NAME = "verification_report.json"
CHECKSUM_NAME = "checksums.sha256"
CSV_COLUMNS = (
    "calendar_year",
    "sample_code",
    "base_head_records",
    "base_head_weight",
    "tenure_valid_records",
    "tenure_valid_weight",
    "owner_records",
    "owner_weight",
    "ownership_rate",
    "rooms_valid_records",
    "rooms_valid_weight",
    "rooms_weighted_sum",
    "mean_rooms_literal",
    "ownership_rate_index_2007",
    "mean_rooms_literal_index_2007",
)
FEMALE_EXPOSURE_COLUMNS = (
    "calendar_year",
    "sample_code",
    "age_lower",
    "age_upper",
    "female_records",
    "female_perwt_mass",
    "female_exposure_share_18_45",
)
FEMALE_CODE = 2
FEMALE_AGE_MIN = 18
FEMALE_AGE_MAX = 45
FEMALE_AGE_BIN_YEARS = 4
FEMALE_AGE_BINS = tuple(
    (lower, lower + FEMALE_AGE_BIN_YEARS - 1)
    for lower in range(
        FEMALE_AGE_MIN, FEMALE_AGE_MAX + 1, FEMALE_AGE_BIN_YEARS
    )
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, default=DEFAULT_EXTRACT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--expected-source-sha256",
        default=EXPECTED_SOURCE_SHA256,
        help="Fail unless the source extract has this SHA-256 fingerprint.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help=(
            "Recompute from the raw extract, require byte-identical canonical "
            "CSV and metadata, and refresh only the verification report."
        ),
    )
    return parser.parse_args()


def json_bytes(payload: Any) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def csv_bytes(rows: list[dict[str, Any]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=list(CSV_COLUMNS),
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        rendered = dict(row)
        for name in (
            "ownership_rate",
            "mean_rooms_literal",
            "ownership_rate_index_2007",
            "mean_rooms_literal_index_2007",
        ):
            rendered[name] = f"{float(row[name]):.15f}"
        for name in (
            "base_head_weight",
            "tenure_valid_weight",
            "owner_weight",
            "rooms_valid_weight",
            "rooms_weighted_sum",
        ):
            value = float(row[name])
            if not np.isclose(value, round(value), atol=1e-7, rtol=0.0):
                raise RuntimeError(f"Expected integer-valued ACS total for {name}")
            rendered[name] = str(int(round(value)))
        writer.writerow(rendered)
    return stream.getvalue().encode("utf-8")


def female_exposure_csv_bytes(rows: list[dict[str, Any]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=list(FEMALE_EXPOSURE_COLUMNS),
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        rendered = dict(row)
        mass = float(row["female_perwt_mass"])
        if not np.isclose(mass, round(mass), atol=1e-7, rtol=0.0):
            raise RuntimeError("Expected integer-valued ACS PERWT mass")
        rendered["female_perwt_mass"] = str(int(round(mass)))
        rendered["female_exposure_share_18_45"] = (
            f"{float(row['female_exposure_share_18_45']):.15f}"
        )
        writer.writerow(rendered)
    return stream.getvalue().encode("utf-8")


def selected_value_labels(
    reader: pd.io.stata.StataReader, variable: str, codes: tuple[int, ...]
) -> dict[str, str]:
    index = reader.varlist.index(variable)
    label_name = reader.lbllist[index]
    labels = reader.value_labels().get(label_name, {}) if label_name else {}
    return {str(code): str(labels.get(code, "")) for code in codes}


def build_contract() -> dict[str, Any]:
    return {
        "artifact_id": "national_householder_housing_path_v1",
        "calendar_years": list(YEARS),
        "annual_sample_codes": {
            str(year): EXPECTED_SAMPLE_CODES[year] for year in YEARS
        },
        "base_sample": [
            "YEAR equals calendar year",
            "GQ in {1,2}",
            "PERNUM equals 1",
            "RELATE equals 1",
            "HHWT finite and strictly positive",
            f"AGE in [{AGE_MIN},{AGE_MAX}]",
        ],
        "weight": "HHWT",
        "ownership_estimand": (
            "HHWT-weighted share with OWNERSHP=1 among base-sample heads "
            "with OWNERSHP in {1,2}"
        ),
        "rooms_estimand": (
            "HHWT-weighted arithmetic mean of the literal numeric ROOMS value "
            "among base-sample heads with finite ROOMS>0; no top-code adjustment"
        ),
        "rooms_comparability_warning": (
            "IPUMS labels ROOMS=9 as 9+ through 2007, whereas later annual ACS "
            "samples retain values above 9. The requested literal series is "
            "therefore not harmonized across the 2007 coding break; do not "
            "interpret the 2007--2011 change as a pure change in housing quantity."
        ),
        "csv_columns": list(CSV_COLUMNS),
        "female_exposure_artifact": {
            "purpose": (
                "separate denominator for standardizing model age-specific "
                "birth rates; never an adult-household stock"
            ),
            "sample": [
                "YEAR equals calendar year",
                "GQ in {1,2}",
                "SEX equals 2 (Female, verified from the source value label)",
                "PERWT finite and strictly positive",
                f"AGE in [{FEMALE_AGE_MIN},{FEMALE_AGE_MAX}]",
            ],
            "weight": "PERWT",
            "age_bins": [list(pair) for pair in FEMALE_AGE_BINS],
            "csv_columns": list(FEMALE_EXPOSURE_COLUMNS),
        },
        "expected_source_sha256": EXPECTED_SOURCE_SHA256,
    }


def compute(
    extract: Path, expected_source_sha256: str
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
    list[dict[str, Any]],
]:
    if not extract.is_file():
        raise FileNotFoundError(extract)
    source_sha256 = sha256_file(extract)
    if source_sha256 != expected_source_sha256:
        raise RuntimeError(
            "Source fingerprint mismatch: "
            f"expected {expected_source_sha256}, observed {source_sha256}"
        )

    contract = build_contract()
    contract_sha256 = sha256_bytes(
        json.dumps(contract, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )
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
            raise RuntimeError("Unsupported pandas StataReader fixed-width metadata")
        field = dict(zip(reader.varlist, raw_names, strict=True))
        required = {
            "year",
            "sample",
            "gq",
            "pernum",
            "relate",
            "hhwt",
            "age",
            "ownershp",
            "rooms",
            "sex",
            "perwt",
        }
        missing = sorted(required - set(field))
        if missing:
            raise ValueError(f"ACS extract is missing required variables: {missing}")

        record_size = int(dtype.itemsize)
        stream = reader.path_or_buf

        def read_records(start: int, count: int) -> np.ndarray:
            stream.seek(int(data_location) + int(start) * record_size)
            raw = np.frombuffer(
                stream.read(int(count) * record_size),
                dtype=dtype,
                count=int(count),
            )
            if raw.size != count:
                raise RuntimeError(
                    f"Short Stata read at record {start}: {raw.size} != {count}"
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

        variable_labels = reader.variable_labels()
        data_dictionary = {
            variable.upper(): {
                "stata_name": variable,
                "variable_label": str(variable_labels.get(variable, "")),
                "storage_dtype": str(dtype.fields[field[variable]][0]),
            }
            for variable in sorted(required)
        }
        data_dictionary["GQ"]["selected_codes"] = selected_value_labels(
            reader, "gq", (1, 2)
        )
        data_dictionary["RELATE"]["selected_codes"] = selected_value_labels(
            reader, "relate", (1,)
        )
        data_dictionary["OWNERSHP"]["code_definitions"] = selected_value_labels(
            reader, "ownershp", (0, 1, 2)
        )
        data_dictionary["ROOMS"]["code_definitions"] = selected_value_labels(
            reader, "rooms", tuple(range(0, 28)) + (30,)
        )
        data_dictionary["SEX"]["code_definitions"] = selected_value_labels(
            reader, "sex", (1, 2, 9)
        )
        if data_dictionary["SEX"]["code_definitions"].get(str(FEMALE_CODE)) != "Female":
            raise RuntimeError(
                "Source SEX coding does not identify code 2 as Female; refusing to guess"
            )

        rows: list[dict[str, Any]] = []
        female_exposure_rows: list[dict[str, Any]] = []
        householder_receipts: dict[str, Any] = {}
        female_exposure_receipts: dict[str, Any] = {}
        checks: list[dict[str, Any]] = []
        for year in YEARS:
            first = lower_year_bound(year)
            last = lower_year_bound(year + 1)
            if first >= last:
                raise RuntimeError(f"No records found for calendar year {year}")

            sample_codes: set[int] = set()
            base_records = 0
            base_weight = 0.0
            tenure_valid_records = 0
            tenure_valid_weight = 0.0
            owner_records = 0
            owner_weight = 0.0
            rooms_valid_records = 0
            rooms_valid_weight = 0.0
            rooms_weighted_sum = 0.0
            female_bin_records = np.zeros(len(FEMALE_AGE_BINS), dtype=np.int64)
            female_bin_weights = np.zeros(len(FEMALE_AGE_BINS), dtype=float)
            household_person_records_18_45 = 0
            household_person_perwt_18_45 = 0.0

            for start in range(first, last, CHUNK_SIZE):
                raw = read_records(start, min(CHUNK_SIZE, last - start))
                sample_codes.update(
                    int(value) for value in np.unique(raw[field["sample"]])
                )
                hhwt_all = raw[field["hhwt"]].astype(float)
                age = raw[field["age"]]
                base = (
                    (raw[field["year"]] == year)
                    & np.isin(raw[field["gq"]], (1, 2))
                    & (raw[field["pernum"]] == 1)
                    & (raw[field["relate"]] == 1)
                    & np.isfinite(hhwt_all)
                    & (hhwt_all > 0.0)
                    & (age >= AGE_MIN)
                    & (age <= AGE_MAX)
                )
                hhwt = hhwt_all[base]
                tenure = raw[field["ownershp"]][base]
                rooms = raw[field["rooms"]][base].astype(float)

                base_records += int(np.sum(base))
                base_weight += float(np.sum(hhwt))

                tenure_valid = np.isin(tenure, (1, 2))
                owner = tenure_valid & (tenure == 1)
                tenure_valid_records += int(np.sum(tenure_valid))
                tenure_valid_weight += float(np.sum(hhwt[tenure_valid]))
                owner_records += int(np.sum(owner))
                owner_weight += float(np.sum(hhwt[owner]))

                rooms_valid = np.isfinite(rooms) & (rooms > 0.0)
                rooms_valid_records += int(np.sum(rooms_valid))
                rooms_valid_weight += float(np.sum(hhwt[rooms_valid]))
                rooms_weighted_sum += float(
                    np.sum(hhwt[rooms_valid] * rooms[rooms_valid])
                )

                perwt_all = raw[field["perwt"]].astype(float)
                person_18_45 = (
                    (raw[field["year"]] == year)
                    & np.isin(raw[field["gq"]], (1, 2))
                    & np.isfinite(perwt_all)
                    & (perwt_all > 0.0)
                    & (age >= FEMALE_AGE_MIN)
                    & (age <= FEMALE_AGE_MAX)
                )
                household_person_records_18_45 += int(np.sum(person_18_45))
                household_person_perwt_18_45 += float(
                    np.sum(perwt_all[person_18_45])
                )
                female = person_18_45 & (raw[field["sex"]] == FEMALE_CODE)
                female_ages = age[female].astype(int)
                female_perwt = perwt_all[female]
                female_bins = (
                    (female_ages - FEMALE_AGE_MIN) // FEMALE_AGE_BIN_YEARS
                )
                female_bin_records += np.bincount(
                    female_bins, minlength=len(FEMALE_AGE_BINS)
                )
                female_bin_weights += np.bincount(
                    female_bins,
                    weights=female_perwt,
                    minlength=len(FEMALE_AGE_BINS),
                )

            observed_codes = sorted(sample_codes)
            expected_code = EXPECTED_SAMPLE_CODES[year]
            if observed_codes != [expected_code]:
                raise RuntimeError(
                    f"Unexpected sample code(s) in {year}: {observed_codes}; "
                    f"expected [{expected_code}]"
                )
            if base_weight <= 0 or tenure_valid_weight <= 0 or rooms_valid_weight <= 0:
                raise RuntimeError(f"Nonpositive estimation denominator in {year}")

            ownership_rate = owner_weight / tenure_valid_weight
            mean_rooms = rooms_weighted_sum / rooms_valid_weight
            row = {
                "calendar_year": year,
                "sample_code": expected_code,
                "base_head_records": base_records,
                "base_head_weight": base_weight,
                "tenure_valid_records": tenure_valid_records,
                "tenure_valid_weight": tenure_valid_weight,
                "owner_records": owner_records,
                "owner_weight": owner_weight,
                "ownership_rate": ownership_rate,
                "rooms_valid_records": rooms_valid_records,
                "rooms_valid_weight": rooms_valid_weight,
                "rooms_weighted_sum": rooms_weighted_sum,
                "mean_rooms_literal": mean_rooms,
            }
            rows.append(row)
            householder_receipts[str(year)] = {
                "raw_year_records": int(last - first),
                "sample_codes_observed": observed_codes,
                "base_head_records": base_records,
                "base_head_weight": base_weight,
                "tenure_invalid_records": base_records - tenure_valid_records,
                "tenure_invalid_weight": base_weight - tenure_valid_weight,
                "rooms_invalid_records": base_records - rooms_valid_records,
                "rooms_invalid_weight": base_weight - rooms_valid_weight,
            }
            female_total_records = int(np.sum(female_bin_records))
            female_total_weight = float(np.sum(female_bin_weights))
            if female_total_records <= 0 or female_total_weight <= 0.0:
                raise RuntimeError(f"Nonpositive female exposure denominator in {year}")
            if np.any(female_bin_records <= 0) or np.any(female_bin_weights <= 0.0):
                raise RuntimeError(f"Empty female exposure age bin in {year}")
            for (age_lower, age_upper), record_count, weight in zip(
                FEMALE_AGE_BINS,
                female_bin_records,
                female_bin_weights,
                strict=True,
            ):
                female_exposure_rows.append(
                    {
                        "calendar_year": year,
                        "sample_code": expected_code,
                        "age_lower": age_lower,
                        "age_upper": age_upper,
                        "female_records": int(record_count),
                        "female_perwt_mass": float(weight),
                        "female_exposure_share_18_45": (
                            float(weight) / female_total_weight
                        ),
                    }
                )
            female_exposure_receipts[str(year)] = {
                "raw_year_records": int(last - first),
                "sample_codes_observed": observed_codes,
                "household_person_records_age_18_45_with_valid_perwt": (
                    household_person_records_18_45
                ),
                "household_person_perwt_mass_age_18_45": (
                    household_person_perwt_18_45
                ),
                "female_records_age_18_45": female_total_records,
                "female_perwt_mass_age_18_45": female_total_weight,
            }
            print(
                "ACS_HOUSING_PATH "
                f"year={year} heads={base_records} head_weight={base_weight:.0f} "
                f"ownership={ownership_rate:.9f} rooms={mean_rooms:.9f}",
                flush=True,
            )

        ownership_2007 = float(rows[0]["ownership_rate"])
        rooms_2007 = float(rows[0]["mean_rooms_literal"])
        for row in rows:
            row["ownership_rate_index_2007"] = (
                float(row["ownership_rate"]) / ownership_2007
            )
            row["mean_rooms_literal_index_2007"] = (
                float(row["mean_rooms_literal"]) / rooms_2007
            )

        def record_check(name: str, passed: bool, detail: str) -> None:
            checks.append({"check": name, "passed": bool(passed), "detail": detail})
            if not passed:
                raise RuntimeError(f"Self-check failed: {name}: {detail}")

        record_check(
            "canonical_years_and_order",
            [row["calendar_year"] for row in rows] == list(YEARS),
            str([row["calendar_year"] for row in rows]),
        )
        record_check(
            "sample_codes",
            all(
                row["sample_code"] == EXPECTED_SAMPLE_CODES[row["calendar_year"]]
                for row in rows
            ),
            str([row["sample_code"] for row in rows]),
        )
        record_check(
            "valid_denominators_nested",
            all(
                0 < row["tenure_valid_records"] <= row["base_head_records"]
                and 0 < row["rooms_valid_records"] <= row["base_head_records"]
                and 0 < row["tenure_valid_weight"] <= row["base_head_weight"]
                and 0 < row["rooms_valid_weight"] <= row["base_head_weight"]
                for row in rows
            ),
            "valid tenure and room samples lie within the base sample",
        )
        record_check(
            "owner_numerator_nested",
            all(
                0 <= row["owner_records"] <= row["tenure_valid_records"]
                and 0 <= row["owner_weight"] <= row["tenure_valid_weight"]
                for row in rows
            ),
            "owner numerator lies within valid-tenure denominator",
        )
        record_check(
            "outcome_ranges",
            all(
                0.0 <= row["ownership_rate"] <= 1.0
                and 1.0 <= row["mean_rooms_literal"] <= 30.0
                for row in rows
            ),
            "ownership in [0,1], literal mean rooms in [1,30]",
        )
        record_check(
            "base_indices",
            rows[0]["ownership_rate_index_2007"] == 1.0
            and rows[0]["mean_rooms_literal_index_2007"] == 1.0,
            "both 2007 indices equal one exactly",
        )
        record_check(
            "female_exposure_rows",
            len(female_exposure_rows) == len(YEARS) * len(FEMALE_AGE_BINS),
            f"{len(female_exposure_rows)} year-by-age-bin rows",
        )
        female_share_sums = {
            year: sum(
                row["female_exposure_share_18_45"]
                for row in female_exposure_rows
                if row["calendar_year"] == year
            )
            for year in YEARS
        }
        record_check(
            "female_exposure_shares_sum_to_one",
            all(
                np.isclose(value, 1.0, atol=1e-14, rtol=0.0)
                for value in female_share_sums.values()
            ),
            str(female_share_sums),
        )

        sample_label_definitions = {
            str(year): selected_value_labels(
                reader, "sample", (EXPECTED_SAMPLE_CODES[year],)
            )[str(EXPECTED_SAMPLE_CODES[year])]
            for year in YEARS
        }
    finally:
        reader.close()

    builder_path = Path(__file__).resolve()
    metadata = {
        "artifact_id": contract["artifact_id"],
        "status": "complete",
        "contract": contract,
        "contract_sha256": contract_sha256,
        "source": {
            "path": str(extract),
            "repository_relative_path": str(extract.relative_to(ROOT)),
            "sha256": source_sha256,
            "size_bytes": extract.stat().st_size,
        },
        "builder": {
            "path": str(builder_path),
            "repository_relative_path": str(builder_path.relative_to(ROOT)),
            "sha256": sha256_file(builder_path),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
        },
        "data_dictionary": data_dictionary,
        "annual_sample_labels": sample_label_definitions,
        "receipts": {
            "householder_housing": householder_receipts,
            "female_exposure": female_exposure_receipts,
        },
        "self_checks": checks,
    }
    return rows, female_exposure_rows, metadata, checks


def write_verification_and_checksums(
    outdir: Path,
    report: dict[str, Any],
) -> None:
    report_path = outdir / VERIFICATION_NAME
    atomic_write(report_path, json_bytes(report))
    names = (CSV_NAME, FEMALE_EXPOSURE_CSV_NAME, METADATA_NAME, VERIFICATION_NAME)
    checksum_lines = [
        f"{sha256_file(outdir / name)}  {name}" for name in names
    ]
    atomic_write(
        outdir / CHECKSUM_NAME,
        ("\n".join(checksum_lines) + "\n").encode("utf-8"),
    )


def main() -> None:
    args = parse_args()
    extract = args.extract.resolve()
    outdir = args.outdir.resolve()
    rows, female_exposure_rows, metadata, checks = compute(
        extract, args.expected_source_sha256
    )
    expected_csv = csv_bytes(rows)
    expected_female_exposure_csv = female_exposure_csv_bytes(female_exposure_rows)
    metadata["output_files"] = {
        CSV_NAME: {
            "rows": len(rows),
            "bytes": len(expected_csv),
            "sha256": sha256_bytes(expected_csv),
        },
        FEMALE_EXPOSURE_CSV_NAME: {
            "rows": len(female_exposure_rows),
            "bytes": len(expected_female_exposure_csv),
            "sha256": sha256_bytes(expected_female_exposure_csv),
        },
    }
    expected_metadata = json_bytes(metadata)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / CSV_NAME
    female_exposure_csv_path = outdir / FEMALE_EXPOSURE_CSV_NAME
    metadata_path = outdir / METADATA_NAME
    if args.verify:
        missing = [
            str(path)
            for path in (csv_path, female_exposure_csv_path, metadata_path)
            if not path.is_file()
        ]
        if missing:
            raise FileNotFoundError(
                "Verification requires canonical artifacts: " + ", ".join(missing)
            )
        observed_csv = csv_path.read_bytes()
        observed_female_exposure_csv = female_exposure_csv_path.read_bytes()
        observed_metadata = metadata_path.read_bytes()
        csv_exact = observed_csv == expected_csv
        female_exposure_csv_exact = (
            observed_female_exposure_csv == expected_female_exposure_csv
        )
        metadata_exact = observed_metadata == expected_metadata
        report = {
            "artifact_id": metadata["artifact_id"],
            "status": (
                "pass"
                if csv_exact and female_exposure_csv_exact and metadata_exact
                else "fail"
            ),
            "mode": "independent_recompute_and_byte_compare",
            "contract_sha256": metadata["contract_sha256"],
            "source_sha256": metadata["source"]["sha256"],
            "canonical_csv_byte_exact": csv_exact,
            "canonical_female_exposure_csv_byte_exact": female_exposure_csv_exact,
            "canonical_metadata_byte_exact": metadata_exact,
            "expected_csv_sha256": sha256_bytes(expected_csv),
            "observed_csv_sha256": sha256_bytes(observed_csv),
            "expected_female_exposure_csv_sha256": sha256_bytes(
                expected_female_exposure_csv
            ),
            "observed_female_exposure_csv_sha256": sha256_bytes(
                observed_female_exposure_csv
            ),
            "expected_metadata_sha256": sha256_bytes(expected_metadata),
            "observed_metadata_sha256": sha256_bytes(observed_metadata),
            "self_checks": checks,
        }
        write_verification_and_checksums(outdir, report)
        if not csv_exact or not female_exposure_csv_exact or not metadata_exact:
            raise RuntimeError("Canonical artifacts differ from exact recomputation")
        print(
            "ACS_HOUSING_PATH_VERIFY_PASS "
            f"contract_sha256={metadata['contract_sha256']} "
            f"csv_sha256={sha256_bytes(expected_csv)}",
            flush=True,
        )
        return

    atomic_write(csv_path, expected_csv)
    atomic_write(female_exposure_csv_path, expected_female_exposure_csv)
    atomic_write(metadata_path, expected_metadata)
    csv_disk_exact = csv_path.read_bytes() == expected_csv
    female_exposure_csv_disk_exact = (
        female_exposure_csv_path.read_bytes() == expected_female_exposure_csv
    )
    metadata_disk_exact = metadata_path.read_bytes() == expected_metadata
    report = {
        "artifact_id": metadata["artifact_id"],
        "status": (
            "pass"
            if csv_disk_exact
            and female_exposure_csv_disk_exact
            and metadata_disk_exact
            else "fail"
        ),
        "mode": "build_write_readback_self_check",
        "contract_sha256": metadata["contract_sha256"],
        "source_sha256": metadata["source"]["sha256"],
        "canonical_csv_byte_exact": csv_disk_exact,
        "canonical_female_exposure_csv_byte_exact": (
            female_exposure_csv_disk_exact
        ),
        "canonical_metadata_byte_exact": metadata_disk_exact,
        "canonical_csv_sha256": sha256_file(csv_path),
        "canonical_female_exposure_csv_sha256": sha256_file(
            female_exposure_csv_path
        ),
        "canonical_metadata_sha256": sha256_file(metadata_path),
        "self_checks": checks,
    }
    write_verification_and_checksums(outdir, report)
    if (
        not csv_disk_exact
        or not female_exposure_csv_disk_exact
        or not metadata_disk_exact
    ):
        raise RuntimeError("Atomic write/readback verification failed")
    print(
        "ACS_HOUSING_PATH_BUILD_COMPLETE "
        f"outdir={outdir} contract_sha256={metadata['contract_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
