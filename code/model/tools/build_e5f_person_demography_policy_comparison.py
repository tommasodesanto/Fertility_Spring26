#!/usr/bin/env python3
"""Compare two coherent-person-law rebated property-tax transitions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


MARKET_TOLERANCE = 2.0e-4
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--reform-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_case(directory: Path, expected_case: str) -> tuple[dict[str, Any], list[dict[str, str]]]:
    root = directory.resolve()
    summary_path = root / "summary.json"
    path_file = root / "transition_path.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
        raise RuntimeError(f"Incomplete path packet: {root}")
    if summary.get("case") != expected_case:
        raise RuntimeError(f"Wrong case in {root}")
    if not all((summary.get("accounting_gates") or {}).values()):
        raise RuntimeError(f"Accounting gates failed in {root}")
    rows = read_csv(path_file)
    if len(rows) != int(summary["horizon"]):
        raise RuntimeError(f"Path row count differs from horizon in {root}")
    return summary, rows


def verify_path_gates(
    summary: dict[str, Any], rows: list[dict[str, str]]
) -> dict[str, Any]:
    market_residuals = np.asarray(
        [abs(float(row["relative_market_residual"])) for row in rows], dtype=float
    )
    fiscal_residuals = np.asarray(
        [abs(float(row["government_budget_residual"])) for row in rows], dtype=float
    )
    finite = bool(
        np.all(np.isfinite(market_residuals))
        and np.all(np.isfinite(fiscal_residuals))
    )
    maximum_market = float(np.max(market_residuals)) if finite else float("inf")
    maximum_fiscal = float(np.max(fiscal_residuals)) if finite else float("inf")
    market_pass = finite and maximum_market <= MARKET_TOLERANCE
    fiscal_pass = finite and maximum_fiscal <= FISCAL_ABSOLUTE_TOLERANCE
    recomputed_path_converged = market_pass and fiscal_pass
    summary_market = float(summary["maximum_market_residual"])
    summary_fiscal = float(summary["maximum_fiscal_residual"])
    summary_matches_rows = bool(
        np.isclose(summary_market, maximum_market, rtol=1e-12, atol=1e-14)
        and np.isclose(summary_fiscal, maximum_fiscal, rtol=1e-12, atol=1e-14)
    )
    summary_flag_consistent = bool(summary.get("path_converged")) == bool(
        recomputed_path_converged
    )
    return {
        "all_path_values_finite": finite,
        "maximum_market_residual_recomputed": maximum_market,
        "maximum_fiscal_residual_recomputed": maximum_fiscal,
        "market_tolerance": MARKET_TOLERANCE,
        "fiscal_absolute_tolerance": FISCAL_ABSOLUTE_TOLERANCE,
        "market_pass": market_pass,
        "fiscal_pass": fiscal_pass,
        "summary_matches_saved_path": summary_matches_rows,
        "summary_path_flag_consistent": summary_flag_consistent,
        "all_checks_pass": bool(
            recomputed_path_converged
            and summary_matches_rows
            and summary_flag_consistent
        ),
    }


def pct(reform: float, baseline: float) -> float:
    return 100.0 * (reform / baseline - 1.0)


def nearest_row(rows: list[dict[str, Any]], year: int) -> dict[str, Any]:
    return min(rows, key=lambda row: abs(int(row["calendar_year"]) - int(year)))


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline_dir = args.baseline_dir.resolve()
    reform_dir = args.reform_dir.resolve()
    baseline_summary, baseline = load_case(
        baseline_dir, "rebated-tax1-baseline"
    )
    reform_summary, reform = load_case(reform_dir, "rebated-tax2-reform")
    if int(baseline_summary["horizon"]) != int(reform_summary["horizon"]):
        raise RuntimeError("Policy paths have different horizons")
    if baseline_summary["source_hashes"] != reform_summary["source_hashes"]:
        raise RuntimeError("Policy paths use different source or input fingerprints")
    if baseline_summary["terminal_root"]["source_hashes"] != reform_summary[
        "terminal_root"
    ]["source_hashes"]:
        raise RuntimeError("Terminal roots use different source or input fingerprints")

    baseline_path_gates = verify_path_gates(baseline_summary, baseline)
    reform_path_gates = verify_path_gates(reform_summary, reform)

    effect_rows: list[dict[str, Any]] = []
    for left, right in zip(baseline, reform, strict=True):
        if int(left["calendar_year"]) != int(right["calendar_year"]):
            raise RuntimeError("Policy path calendars differ")
        left_people = float(left["resident_persons_actual_scale"])
        right_people = float(right["resident_persons_actual_scale"])
        left_heads = float(left["household_heads"])
        right_heads = float(right["household_heads"])
        left_birth_rate = float(left["birth_children_topcode_adjusted"]) / max(
            left_heads * 4.0, 1e-15
        )
        right_birth_rate = float(right["birth_children_topcode_adjusted"]) / max(
            right_heads * 4.0, 1e-15
        )
        effect_rows.append(
            {
                "period": int(left["period"]),
                "calendar_year": int(left["calendar_year"]),
                "asset_price_percent_change": pct(
                    float(right["asset_price"]), float(left["asset_price"])
                ),
                "renter_price_percent_change": pct(
                    float(right["renter_price"]), float(left["renter_price"])
                ),
                "equal_transfer_difference": float(
                    right["equal_transfer_period_units"]
                )
                - float(left["equal_transfer_period_units"]),
                "equal_transfer_percent_change": pct(
                    float(right["equal_transfer_period_units"]),
                    float(left["equal_transfer_period_units"]),
                ),
                "resident_person_difference": right_people - left_people,
                "resident_person_percent_change": pct(right_people, left_people),
                "household_head_difference_model_units": right_heads - left_heads,
                "household_head_percent_change": pct(right_heads, left_heads),
                "annual_births_per_head_baseline": left_birth_rate,
                "annual_births_per_head_reform": right_birth_rate,
                "annual_births_per_head_percent_change": pct(
                    right_birth_rate, left_birth_rate
                ),
                "owner_rate_percentage_point_change": 100.0
                * (float(right["owner_rate"]) - float(left["owner_rate"])),
            }
        )
    with (output_dir / "policy_effects.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(effect_rows[0]))
        writer.writeheader()
        writer.writerows(effect_rows)

    years = np.asarray([row["calendar_year"] for row in effect_rows])
    figure, axes = plt.subplots(3, 2, figsize=(11, 10), constrained_layout=True)
    panels = [
        ("asset_price_percent_change", "Asset-price effect", "Percent"),
        ("renter_price_percent_change", "Rent-equivalent effect", "Percent"),
        (
            "annual_births_per_head_percent_change",
            "Annual-birth-rate effect",
            "Percent",
        ),
        ("resident_person_difference", "Resident-person effect", "Persons"),
        ("household_head_percent_change", "Household-head effect", "Percent"),
        (
            "owner_rate_percentage_point_change",
            "Ownership effect",
            "Percentage points",
        ),
    ]
    for axis, (key, title, ylabel) in zip(axes.flat, panels, strict=True):
        axis.plot(years, [float(row[key]) for row in effect_rows])
        axis.axhline(0.0, color="0.5", linewidth=1)
        axis.grid(alpha=0.2)
        axis.set(title=title, xlabel="Year", ylabel=ylabel)
    figure.savefig(output_dir / "policy_effects.png", dpi=200)
    figure.savefig(output_dir / "policy_effects.pdf")
    plt.close(figure)

    accepted = bool(
        baseline_path_gates["all_checks_pass"]
        and reform_path_gates["all_checks_pass"]
        and (baseline_summary.get("terminal_convergence") or {}).get(
            "all_checks_pass"
        )
        and (reform_summary.get("terminal_convergence") or {}).get(
            "all_checks_pass"
        )
    )
    selected_years = {}
    for requested_year in (2023, 2050, 2080):
        row = nearest_row(effect_rows, requested_year)
        selected_years[str(requested_year)] = {
            "requested_year": requested_year,
            "model_calendar_year": int(row["calendar_year"]),
            **row,
        }
    summary = {
        "status": (
            "complete_accepted_person_demography_policy_comparison"
            if accepted
            else "complete_unaccepted_person_demography_policy_comparison"
        ),
        "accepted": accepted,
        "interpretation": (
            "The comparison is accepted only when both policy paths clear all "
            "date-by-date markets and budgets and both simulated endpoints pass "
            "the declared terminal-state gates."
        ),
        "horizon": int(baseline_summary["horizon"]),
        "baseline_dir": str(baseline_dir),
        "baseline_summary_sha256": sha256(baseline_dir / "summary.json"),
        "baseline_path_sha256": sha256(baseline_dir / "transition_path.csv"),
        "reform_dir": str(reform_dir),
        "reform_summary_sha256": sha256(reform_dir / "summary.json"),
        "reform_path_sha256": sha256(reform_dir / "transition_path.csv"),
        "baseline_path_converged": bool(baseline_summary.get("path_converged")),
        "reform_path_converged": bool(reform_summary.get("path_converged")),
        "baseline_path_gate_audit": baseline_path_gates,
        "reform_path_gate_audit": reform_path_gates,
        "baseline_terminal_converged": bool(
            (baseline_summary.get("terminal_convergence") or {}).get(
                "all_checks_pass"
            )
        ),
        "reform_terminal_converged": bool(
            (reform_summary.get("terminal_convergence") or {}).get(
                "all_checks_pass"
            )
        ),
        "selected_year_effects": selected_years,
        "terminal_root_effects": {
            "resident_person_difference": float(
                reform_summary["terminal_root"]["resident_persons_actual_scale"]
            )
            - float(
                baseline_summary["terminal_root"]["resident_persons_actual_scale"]
            ),
            "resident_person_percent_change": pct(
                float(
                    reform_summary["terminal_root"]["resident_persons_actual_scale"]
                ),
                float(
                    baseline_summary["terminal_root"][
                        "resident_persons_actual_scale"
                    ]
                ),
            ),
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
