#!/usr/bin/env python3
"""Build June 2024 CPS completed-fertility targets and bootstrap covariance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = SCRIPT_DIR / "source_manifest.json"
DEFAULT_CACHE_PATH = SCRIPT_DIR / "cache" / "jun24pub.csv"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output"

MOMENT_KEYS = (
    "tfr",
    "childless_rate",
    "parity_share_0",
    "parity_share_1",
    "parity_share_2",
    "parity_share_3plus",
    "capped_top_bin_mean",
)

MOMENT_METADATA = {
    "tfr": {
        "model_key": "tfr",
        "label": "Completed fertility",
        "unit": "children_per_woman",
        "definition": (
            "Weighted mean children ever born among women ages 40-44; "
            "PTSF1 is capped at 5 in the public-use file."
        ),
        "published_reference": 1.918,
        "source_variable": "PTSF1",
    },
    "childless_rate": {
        "model_key": "childless_rate",
        "label": "Childlessness",
        "unit": "share",
        "definition": "Weighted share of women ages 40-44 with PTSF1 equal to 0.",
        "published_reference": 0.188,
        "source_variable": "PTSF1",
    },
    "parity_share_0": {
        "model_key": "",
        "label": "Parity share: 0",
        "unit": "share",
        "definition": "Weighted share of women ages 40-44 with PTSF1 equal to 0.",
        "published_reference": None,
        "source_variable": "PTSF1",
    },
    "parity_share_1": {
        "model_key": "",
        "label": "Parity share: 1",
        "unit": "share",
        "definition": "Weighted share of women ages 40-44 with PTSF1 equal to 1.",
        "published_reference": None,
        "source_variable": "PTSF1",
    },
    "parity_share_2": {
        "model_key": "",
        "label": "Parity share: 2",
        "unit": "share",
        "definition": "Weighted share of women ages 40-44 with PTSF1 equal to 2.",
        "published_reference": None,
        "source_variable": "PTSF1",
    },
    "parity_share_3plus": {
        "model_key": "",
        "label": "Parity share: 3+",
        "unit": "share",
        "definition": "Weighted share of women ages 40-44 with PTSF1 at least 3.",
        "published_reference": None,
        "source_variable": "PTSF1",
    },
    "capped_top_bin_mean": {
        "model_key": "tfr_top_bin_weight",
        "label": "Capped top-bin mean",
        "unit": "children_per_woman_in_3plus_bin",
        "definition": (
            "Weighted mean PTSF1 among women with PTSF1 at least 3. "
            "Because public-use PTSF1 codes five or more births as 5, this is "
            "E[min(N,5) | N>=3], not an uncapped conditional mean."
        ),
        "published_reference": None,
        "source_variable": "PTSF1",
    },
}

WINDOWS = (
    ("primary_40_44", 40, 44),
    ("shift_41_45", 41, 45),
    ("shift_43_47", 43, 47),
)

EXPECTED_PRIMARY = {
    "tfr": 1.91842485035,
    "childless_rate": 0.188179593021,
    "capped_top_bin_mean": 3.602359422009,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-file",
        type=Path,
        help="Use this local jun24pub.csv instead of the managed cache.",
    )
    parser.add_argument(
        "--refresh-download",
        action="store_true",
        help="Redownload the managed cache and replace it only after verification.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help="Number of person-level bootstrap replicates; must be at least 1000.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260724,
        help="Seed for NumPy's PCG64 random-number generator.",
    )
    args = parser.parse_args()
    if args.bootstrap_replicates < 1000:
        parser.error("--bootstrap-replicates must be at least 1000")
    return args


def load_manifest() -> dict:
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source(path: Path, manifest: dict) -> None:
    actual_size = path.stat().st_size
    expected_size = int(manifest["byte_size"])
    if actual_size != expected_size:
        raise RuntimeError(
            f"Source byte-size mismatch for {path}: "
            f"expected {expected_size}, found {actual_size}"
        )
    actual_hash = sha256_file(path)
    expected_hash = manifest["sha256"]
    if actual_hash != expected_hash:
        raise RuntimeError(
            f"Source SHA-256 mismatch for {path}: "
            f"expected {expected_hash}, found {actual_hash}"
        )


def download_source(url: str, destination: Path, manifest: dict) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".download")
    if temporary.exists():
        temporary.unlink()

    print(f"Downloading {url}")
    started = time.monotonic()
    downloaded = 0
    next_report = 32 * 1024 * 1024
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Fertility-Spring26-CPS-builder/1.0"},
    )
    try:
        with urllib.request.urlopen(request) as response, temporary.open("wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)
                if downloaded >= next_report:
                    print(f"  downloaded {downloaded / (1024**2):.0f} MiB")
                    next_report += 32 * 1024 * 1024
        verify_source(temporary, manifest)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    elapsed = time.monotonic() - started
    print(f"Verified cache at {destination} ({elapsed:.1f}s)")


def resolve_source(args: argparse.Namespace, manifest: dict) -> Path:
    if args.raw_file is not None:
        source = args.raw_file.expanduser().resolve()
        if not source.exists():
            raise FileNotFoundError(f"Missing --raw-file: {source}")
        verify_source(source, manifest)
        return source

    source = DEFAULT_CACHE_PATH
    if args.refresh_download or not source.exists():
        download_source(manifest["source_url"], source, manifest)
    else:
        verify_source(source, manifest)
        print(f"Using verified cache: {source}")
    return source


def parse_int(row: dict[str, str], key: str) -> int:
    return int(row[key])


def read_analysis_records(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    required = {"PESEX", "PRTAGE", "PTSF1", "PWSSWGT"}
    timing_columns = ("PTSF2", "PTSAYFC")
    ages: list[int] = []
    births: list[int] = []
    weights: list[float] = []
    raw_rows = 0
    female_40_47_rows = 0
    excluded_not_in_universe = 0
    excluded_nonpositive_weight = 0

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in {path}")
        header = set(reader.fieldnames)
        missing = sorted(required - header)
        if missing:
            raise RuntimeError(f"Missing required CPS columns: {', '.join(missing)}")
        timing_presence = {column: column in header for column in timing_columns}

        for row in reader:
            raw_rows += 1
            if parse_int(row, "PESEX") != 2:
                continue
            age = parse_int(row, "PRTAGE")
            if not 40 <= age <= 47:
                continue
            female_40_47_rows += 1
            children_ever_born = parse_int(row, "PTSF1")
            if not 0 <= children_ever_born <= 5:
                excluded_not_in_universe += 1
                continue
            weight = float(row["PWSSWGT"])
            if not np.isfinite(weight) or weight <= 0:
                excluded_nonpositive_weight += 1
                continue
            ages.append(age)
            births.append(children_ever_born)
            weights.append(weight)

    diagnostics = {
        "raw_rows": raw_rows,
        "female_age_40_47_rows": female_40_47_rows,
        "excluded_not_in_fertility_universe": excluded_not_in_universe,
        "excluded_nonpositive_weight": excluded_nonpositive_weight,
        "analysis_rows_age_40_47": len(ages),
        "timing_variable_presence": timing_presence,
    }
    return (
        np.asarray(ages, dtype=np.int16),
        np.asarray(births, dtype=np.int8),
        np.asarray(weights, dtype=np.float64),
        diagnostics,
    )


def estimate_vector(births: np.ndarray, weights: np.ndarray) -> np.ndarray:
    total_weight = weights.sum()
    if total_weight <= 0:
        raise RuntimeError("Nonpositive total analysis weight")
    shares = np.asarray(
        [(weights[births == parity].sum() / total_weight) for parity in range(3)]
        + [weights[births >= 3].sum() / total_weight],
        dtype=np.float64,
    )
    top = births >= 3
    top_weight = weights[top].sum()
    if top_weight <= 0:
        raise RuntimeError("No weighted observations in the 3+ parity bin")
    capped_top_bin_mean = np.dot(weights[top], births[top]) / top_weight
    completed_fertility = np.dot(weights, births) / total_weight
    return np.asarray(
        [
            completed_fertility,
            shares[0],
            shares[0],
            shares[1],
            shares[2],
            shares[3],
            capped_top_bin_mean,
        ],
        dtype=np.float64,
    )


def subset_for_window(
    ages: np.ndarray,
    births: np.ndarray,
    weights: np.ndarray,
    age_min: int,
    age_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    keep = (ages >= age_min) & (ages <= age_max)
    return births[keep], weights[keep]


def bootstrap_primary(
    births: np.ndarray,
    weights: np.ndarray,
    replicates: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(births)
    draws = np.empty((replicates, len(MOMENT_KEYS)), dtype=np.float64)
    for draw in range(replicates):
        indices = rng.integers(0, n, size=n)
        draws[draw, :] = estimate_vector(births[indices], weights[indices])
        if (draw + 1) % 100 == 0 or draw + 1 == replicates:
            print(f"  bootstrap {draw + 1}/{replicates}")
    return draws


def float_text(value: float | None) -> str:
    if value is None:
        return ""
    return f"{float(value):.12f}"


def population_text(value: float) -> str:
    return f"{float(value):.4f}"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_targets(
    output_dir: Path,
    estimates: np.ndarray,
    standard_errors: np.ndarray,
    n_persons: int,
    weighted_population: float,
    args: argparse.Namespace,
    manifest: dict,
) -> None:
    rows = []
    for index, key in enumerate(MOMENT_KEYS):
        metadata = MOMENT_METADATA[key]
        reference = metadata["published_reference"]
        rows.append(
            {
                "moment_key": key,
                "model_key": metadata["model_key"],
                "label": metadata["label"],
                "estimate": float_text(estimates[index]),
                "standard_error": float_text(standard_errors[index]),
                "published_reference": float_text(reference),
                "gap_from_published": float_text(
                    None if reference is None else estimates[index] - reference
                ),
                "unit": metadata["unit"],
                "definition": metadata["definition"],
                "source_variable": metadata["source_variable"],
                "weighting": "PWSSWGT supplement final weight",
                "age_min": 40,
                "age_max": 44,
                "n_persons": n_persons,
                "weighted_population": population_text(weighted_population),
                "bootstrap_method": "person-level resampling with replacement",
                "bootstrap_replicates": args.bootstrap_replicates,
                "bootstrap_seed": args.seed,
                "source_url": manifest["source_url"],
            }
        )
    write_csv(
        output_dir / "cps_fertility_targets.csv",
        [
            "moment_key",
            "model_key",
            "label",
            "estimate",
            "standard_error",
            "published_reference",
            "gap_from_published",
            "unit",
            "definition",
            "source_variable",
            "weighting",
            "age_min",
            "age_max",
            "n_persons",
            "weighted_population",
            "bootstrap_method",
            "bootstrap_replicates",
            "bootstrap_seed",
            "source_url",
        ],
        rows,
    )


def write_covariance(
    output_dir: Path,
    covariance: np.ndarray,
    correlation: np.ndarray,
    args: argparse.Namespace,
) -> None:
    rows = []
    for row_index, row_key in enumerate(MOMENT_KEYS):
        for column_index, column_key in enumerate(MOMENT_KEYS):
            rows.append(
                {
                    "row_moment": row_key,
                    "column_moment": column_key,
                    "covariance": float_text(covariance[row_index, column_index]),
                    "correlation": float_text(correlation[row_index, column_index]),
                    "bootstrap_replicates": args.bootstrap_replicates,
                    "bootstrap_seed": args.seed,
                }
            )
    write_csv(
        output_dir / "cps_fertility_covariance.csv",
        [
            "row_moment",
            "column_moment",
            "covariance",
            "correlation",
            "bootstrap_replicates",
            "bootstrap_seed",
        ],
        rows,
    )


def write_sensitivities(
    output_dir: Path,
    ages: np.ndarray,
    births: np.ndarray,
    weights: np.ndarray,
) -> None:
    rows = []
    for window_label, age_min, age_max in WINDOWS:
        window_births, window_weights = subset_for_window(
            ages, births, weights, age_min, age_max
        )
        for weighting, analysis_weights in (
            ("weighted_PWSSWGT", window_weights),
            ("unweighted", np.ones(len(window_births), dtype=np.float64)),
        ):
            estimates = estimate_vector(window_births, analysis_weights)
            weighted_population = int(np.rint(window_weights).sum()) / 10_000
            for index, key in enumerate(MOMENT_KEYS):
                rows.append(
                    {
                        "window": window_label,
                        "age_min": age_min,
                        "age_max": age_max,
                        "weighting": weighting,
                        "moment_key": key,
                        "estimate": float_text(estimates[index]),
                        "n_persons": len(window_births),
                        "weighted_population": population_text(weighted_population),
                    }
                )
    write_csv(
        output_dir / "cps_fertility_sensitivities.csv",
        [
            "window",
            "age_min",
            "age_max",
            "weighting",
            "moment_key",
            "estimate",
            "n_persons",
            "weighted_population",
        ],
        rows,
    )


def write_build_diagnostics(
    output_dir: Path,
    diagnostics: dict,
    args: argparse.Namespace,
    source_path: Path,
    manifest: dict,
) -> None:
    try:
        source_file_used = str(source_path.relative_to(SCRIPT_DIR))
    except ValueError:
        source_file_used = str(source_path)
    payload = {
        "source_file_used": source_file_used,
        "source_sha256": manifest["sha256"],
        "source_byte_size": manifest["byte_size"],
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_seed": args.seed,
        "sample": {
            "sex": "women (PESEX=2)",
            "primary_age_window": "40-44",
            "fertility_universe": "PTSF1 in {0,1,2,3,4,5}",
            "positive_weight": "PWSSWGT > 0",
        },
        **diagnostics,
        "timing_variables_used_for_targets": False,
    }
    with (output_dir / "build_diagnostics.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def verify_outputs(
    estimates: np.ndarray,
    covariance: np.ndarray,
    correlation: np.ndarray,
) -> None:
    by_key = dict(zip(MOMENT_KEYS, estimates))
    for key, expected in EXPECTED_PRIMARY.items():
        difference = abs(by_key[key] - expected)
        if difference > 5e-10:
            raise RuntimeError(
                f"Reproduction gate failed for {key}: "
                f"expected {expected:.12f}, found {by_key[key]:.12f}"
            )
    parity_sum = sum(by_key[key] for key in MOMENT_KEYS[2:6])
    if abs(parity_sum - 1.0) > 1e-12:
        raise RuntimeError(f"Parity shares do not sum to one: {parity_sum}")
    if not np.allclose(covariance, covariance.T, atol=1e-15, rtol=1e-12):
        raise RuntimeError("Bootstrap covariance is not symmetric")
    if not np.isfinite(covariance).all() or not np.isfinite(correlation).all():
        raise RuntimeError("Non-finite value in bootstrap covariance/correlation")


def main() -> int:
    args = parse_args()
    manifest = load_manifest()
    source_path = resolve_source(args, manifest)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Reading CPS fertility-universe records for ages 40-47...")
    ages, births, weights, diagnostics = read_analysis_records(source_path)
    primary_births, primary_weights = subset_for_window(
        ages, births, weights, 40, 44
    )
    primary_estimates = estimate_vector(primary_births, primary_weights)
    weighted_population = int(np.rint(primary_weights).sum()) / 10_000
    print(
        f"Primary sample: n={len(primary_births):,}, "
        f"weighted population={weighted_population:,.0f}"
    )
    print(
        f"Completed fertility={primary_estimates[0]:.12f}; "
        f"childlessness={primary_estimates[1]:.12f}; "
        f"capped top-bin mean={primary_estimates[-1]:.12f}"
    )

    print(
        f"Running seeded person-level bootstrap "
        f"(B={args.bootstrap_replicates}, seed={args.seed})..."
    )
    bootstrap_draws = bootstrap_primary(
        primary_births,
        primary_weights,
        args.bootstrap_replicates,
        args.seed,
    )
    standard_errors = bootstrap_draws.std(axis=0, ddof=1)
    covariance = np.cov(bootstrap_draws, rowvar=False, ddof=1)
    correlation = np.corrcoef(bootstrap_draws, rowvar=False)
    verify_outputs(primary_estimates, covariance, correlation)

    write_targets(
        output_dir,
        primary_estimates,
        standard_errors,
        len(primary_births),
        weighted_population,
        args,
        manifest,
    )
    write_covariance(output_dir, covariance, correlation, args)
    write_sensitivities(output_dir, ages, births, weights)
    write_build_diagnostics(
        output_dir, diagnostics, args, source_path, manifest
    )

    timing_presence = diagnostics["timing_variable_presence"]
    print(
        "First-birth timing variable presence: "
        + ", ".join(
            f"{key}={'present' if value else 'absent'}"
            for key, value in timing_presence.items()
        )
    )
    print("Timing variables were not used to construct any target.")
    print(f"Wrote outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
