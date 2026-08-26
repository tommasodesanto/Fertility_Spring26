#!/usr/bin/env python3
"""Diagnose the asymptotic convergence rate of the terminal person law."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from demographic_transition.person_cohort_law import (  # noqa: E402
    _advance_existing_linear,
)
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402


DEFAULT_SOURCE_DIR = person_pf.DEFAULT_SOURCE_DIR
DEFAULT_HEADSHIP_DIR = person_pf.DEFAULT_HEADSHIP_DIR
MODEL_AGE_START = 18
MODEL_AGE_CELL_WIDTH = 4
ATTENUATION_THRESHOLDS = (0.10, 0.02, 0.01)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terminal-summaries", type=Path, nargs="+", required=True)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot-periods", type=int, default=160)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def build_existing_cohort_operator(survival: np.ndarray) -> np.ndarray:
    shape = tuple(np.asarray(survival).shape)
    number_of_cells = int(np.prod(shape))
    operator = np.zeros((number_of_cells, number_of_cells))
    for column in range(number_of_cells):
        basis = np.zeros(shape)
        basis.flat[column] = 1.0
        operator[:, column] = _advance_existing_linear(
            basis,
            survival,
            topcoded_last_age=True,
        ).ravel()
    return operator


def build_four_year_operator(
    *,
    annual_existing_operator: np.ndarray,
    newborn_survivor_vector: np.ndarray,
    headship_rates: np.ndarray,
    annual_births_per_head: float,
) -> np.ndarray:
    annual = np.asarray(annual_existing_operator, dtype=float)
    newborn = np.asarray(newborn_survivor_vector, dtype=float).reshape(-1)
    headship = np.asarray(headship_rates, dtype=float).reshape(-1)
    if annual.shape != (newborn.size, newborn.size) or headship.size != newborn.size:
        raise ValueError("Terminal demographic operator dimensions differ")
    four_year_existing = np.linalg.matrix_power(annual, 4)
    four_year_birth_response = sum(
        np.linalg.matrix_power(annual, power) @ newborn for power in range(4)
    )
    return four_year_existing + float(annual_births_per_head) * np.outer(
        four_year_birth_response, headship
    )


def validate_block_operator(
    *,
    four_year_operator: np.ndarray,
    annual_existing_operator: np.ndarray,
    newborn_survivor_vector: np.ndarray,
    headship_rates: np.ndarray,
    annual_births_per_head: float,
) -> float:
    initial = np.linspace(0.2, 1.2, four_year_operator.shape[0])
    annual_births = float(annual_births_per_head) * float(
        np.dot(np.asarray(headship_rates).reshape(-1), initial)
    )
    direct = initial.copy()
    for _ in range(4):
        direct = (
            np.asarray(annual_existing_operator) @ direct
            + np.asarray(newborn_survivor_vector).reshape(-1) * annual_births
        )
    matrix_result = np.asarray(four_year_operator) @ initial
    return float(np.max(np.abs(direct - matrix_result)))


def minimum_periods(spectral_radius: float, threshold: float) -> int:
    return int(math.ceil(math.log(float(threshold)) / math.log(spectral_radius)))


def main() -> None:
    args = parse_args()
    if int(args.plot_periods) < 2:
        raise ValueError("--plot-periods must be at least two")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    source_dir = args.source_dir.resolve()
    headship_dir = args.headship_dir.resolve()
    source_paths = {
        "births_mid": source_dir / "births_mid.csv",
        "survival": source_dir / "survival.csv",
        "acs_headship": headship_dir / "acs_headship_profiles.csv",
    }
    demographic_source = person_pf.demographic_source
    births = demographic_source.births_path(
        demographic_source.read_csv(source_paths["births_mid"])
    )
    survival = demographic_source.survival_path(
        demographic_source.read_csv(source_paths["survival"]), nativity=0
    )[2100]
    raw_headship = demographic_source.fixed_2023_headship(
        demographic_source.read_csv(source_paths["acs_headship"])
    )
    birth_sex_shares = np.asarray(births[2100], dtype=float)
    birth_sex_shares /= float(np.sum(birth_sex_shares))
    newborn_survivor = np.zeros_like(survival)
    newborn_survivor[:, 0] = birth_sex_shares * survival[:, 0]
    annual_existing = build_existing_cohort_operator(survival)

    terminal_summaries = []
    for path in args.terminal_summaries:
        resolved = path.resolve()
        summary = json.loads(resolved.read_text(encoding="utf-8"))
        if summary.get("status") != "complete_corrected_terminal_root":
            raise RuntimeError(f"Incomplete terminal root: {resolved}")
        if not all((summary.get("gates") or {}).values()):
            raise RuntimeError(f"Terminal root gates failed: {resolved}")
        expected_inputs = {
            "input_births_mid": sha256(source_paths["births_mid"]),
            "input_survival": sha256(source_paths["survival"]),
            "input_acs_headship": sha256(source_paths["acs_headship"]),
        }
        for name, digest in expected_inputs.items():
            if (summary.get("source_hashes") or {}).get(name) != digest:
                raise RuntimeError(f"{resolved} uses a different {name} fingerprint")
        terminal_summaries.append((resolved, summary))
    if len({summary["case"] for _, summary in terminal_summaries}) != len(
        terminal_summaries
    ):
        raise RuntimeError("Duplicate terminal policy cases were supplied")

    rows: list[dict[str, Any]] = []
    decay_rows: list[dict[str, Any]] = []
    for path, summary in terminal_summaries:
        aligned_headship = raw_headship.copy()
        factors = np.asarray(summary["headship_alignment_factors"], dtype=float)
        for index, factor in enumerate(factors):
            lower = MODEL_AGE_START + MODEL_AGE_CELL_WIDTH * index
            upper = lower + MODEL_AGE_CELL_WIDTH
            aligned_headship[:, lower:upper] *= float(factor)
        operator = build_four_year_operator(
            annual_existing_operator=annual_existing,
            newborn_survivor_vector=newborn_survivor,
            headship_rates=aligned_headship,
            annual_births_per_head=float(summary["annual_births_per_head"]),
        )
        validation = validate_block_operator(
            four_year_operator=operator,
            annual_existing_operator=annual_existing,
            newborn_survivor_vector=newborn_survivor,
            headship_rates=aligned_headship,
            annual_births_per_head=float(summary["annual_births_per_head"]),
        )
        if validation > 1e-12:
            raise RuntimeError(f"Four-year block-operator identity failed: {validation}")
        eigenvalue_moduli = np.sort(np.abs(np.linalg.eigvals(operator)))[::-1]
        spectral_radius = float(eigenvalue_moduli[0])
        if not 0.0 < spectral_radius < 1.0:
            raise RuntimeError(
                f"Terminal demographic operator is not stable: {spectral_radius}"
            )
        row: dict[str, Any] = {
            "case": summary["case"],
            "terminal_summary": str(path),
            "terminal_summary_sha256": sha256(path),
            "annual_births_per_head": float(summary["annual_births_per_head"]),
            "renewal_ratio": float(summary["renewal_ratio"]),
            "four_year_spectral_radius": spectral_radius,
            "spectral_half_life_years": 4.0 * math.log(0.5) / math.log(spectral_radius),
            "block_operator_validation_max_abs": validation,
            "second_eigenvalue_modulus": float(eigenvalue_moduli[1]),
        }
        for threshold in ATTENUATION_THRESHOLDS:
            label = f"{100.0 * threshold:g}pct"
            periods = minimum_periods(spectral_radius, threshold)
            row[f"periods_to_{label}_modal_attenuation"] = periods
            row[f"years_to_{label}_modal_attenuation"] = 4 * periods
        rows.append(row)
        for period in range(int(args.plot_periods) + 1):
            decay_rows.append(
                {
                    "case": summary["case"],
                    "period": period,
                    "years_after_2023": 4 * period,
                    "dominant_mode_remaining": spectral_radius**period,
                }
            )

    write_csv(output_dir / "terminal_spectral_summary.csv", rows)
    write_csv(output_dir / "terminal_modal_decay.csv", decay_rows)
    figure, axis = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    for case in [summary["case"] for _, summary in terminal_summaries]:
        selected = [row for row in decay_rows if row["case"] == case]
        axis.semilogy(
            [row["years_after_2023"] for row in selected],
            [row["dominant_mode_remaining"] for row in selected],
            label=case,
        )
    for threshold in ATTENUATION_THRESHOLDS:
        axis.axhline(threshold, color="0.65", linestyle=":", linewidth=1)
    axis.set(
        title="Terminal demographic convergence rate",
        xlabel="Years after 2023",
        ylabel="Dominant-mode share remaining",
    )
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    figure.savefig(output_dir / "terminal_modal_decay.png", dpi=200)
    figure.savefig(output_dir / "terminal_modal_decay.pdf")
    plt.close(figure)

    summary_payload = {
        "status": "complete_person_demography_terminal_spectral_diagnostic",
        "promotion_status": "diagnostic_only",
        "interpretation": (
            "The spectral calculation holds the terminal survival, migration, "
            "birth-sex shares, aligned headship schedule, and births-per-head "
            "rule fixed. It gives the local asymptotic convergence rate of the "
            "person block, not a guarantee that a full nonlinear equilibrium "
            "path or household distribution has reached its terminal state."
        ),
        "period_years": 4,
        "terminal_cases": rows,
        "source_hashes": {name: sha256(path) for name, path in source_paths.items()},
        "driver_sha256": sha256(Path(__file__).resolve()),
    }
    write_json(output_dir / "summary.json", summary_payload)


if __name__ == "__main__":
    main()
