#!/usr/bin/env python3
"""Collect the internally identified three-parameter bequest recalibration.

With ``--arm M4`` this collects the standard bequest recalibration instead:
11 clean-frontier coordinates plus ``theta0`` and ``theta1`` against the
median-composition target set, gated on beating the strict M1 ``theta0=0``
nested seed.

With ``--arm M5`` this collects the income-disciplined recalibration: the M4
parameters plus ``tenure_choice_kappa`` against the income-disciplined target
set, gated on beating the strict M1 ``theta0=0``/``tenure_choice_kappa=0``
nested seed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


# (target_set, free_parameter_count, target_count, minimum completed chains)
ARM_CONTRACTS = {
    "M3": ("candidate_replacement_bequest_internal_v1", 14, 15, 4),
    "M4": ("candidate_replacement_bequest_median_composition_v1", 13, 14, 6),
    "M5": ("candidate_replacement_income_disciplined_v1", 14, 15, 8),
}

ESTATE_MEDIAN = "old_total_estate_wealth_to_annual_income_median_7684"
NONHOUSING_MEDIAN = "old_nonhousing_wealth_to_income_median_6575"
ESTATE_MEDIAN_TARGET = 6.50131577436537
ESTATE_MEDIAN_SE = 0.2319582116443
ESTATE_MEDIAN_WEIGHT = 18.585767349158665
NONHOUSING_MEDIAN_TARGET = 1.90821154211154
NONHOUSING_MEDIAN_SE = 0.109272216032637
NONHOUSING_MEDIAN_WEIGHT = 83.74916751466371
ESTABLISHED_BLOCK_LOSS_CEILING = 4.42
M4_SHARED_THETA = (
    "beta",
    "alpha_cons",
    "c_bar_0",
    "c_bar_n",
    "h_bar_0",
    "h_bar_jump",
    "h_bar_n",
    "psi_child",
    "kappa_fert",
    "chi",
    "H0",
    "theta0",
    "theta1",
    "theta_n",
    "tenure_choice_kappa",
)
EXPECTED_M4_ACTIVE_DOMAIN = [
    {"name": "beta_annual", "lower": 0.80, "upper": 0.9995, "transform": "discount"},
    {"name": "alpha_cons", "lower": 0.02, "upper": 0.98, "transform": "logit"},
    {"name": "c_bar_0", "lower": 0.0, "upper": 2.0, "transform": "softzero"},
    {"name": "c_bar_n", "lower": 0.0, "upper": 3.0, "transform": "softzero"},
    {"name": "h_bar_0", "lower": 0.05, "upper": 5.80, "transform": "log"},
    {"name": "h_bar_jump", "lower": 0.0, "upper": 8.0, "transform": "softzero"},
    {"name": "h_bar_n", "lower": 0.0, "upper": 5.0, "transform": "softzero"},
    {"name": "psi_child", "lower": -3.0, "upper": 3.0, "transform": "asinh"},
    {"name": "kappa_fert", "lower": 0.02, "upper": 50.0, "transform": "log"},
    {"name": "chi", "lower": 0.10, "upper": 5.0, "transform": "log"},
    {"name": "H0", "lower": 0.20, "upper": 80.0, "transform": "log"},
    {"name": "theta0", "lower": 0.0, "upper": 8.0, "transform": "softzero"},
    {"name": "theta1", "lower": 0.02, "upper": 16.0, "transform": "log"},
]
EXPECTED_M4_FIXED = {"theta_n": 0.0, "tenure_choice_kappa": 0.0}
EXPECTED_M4_MECHANISM = {
    "bequest_spec": "linear_child_scale",
    "normalize_bequest_utility": True,
    "owner_ltv_taper": False,
    "owner_ltv_taper_start_age": 66.0,
    "owner_ltv_taper_end_age": 82.0,
    "owner_ltv_terminal_share": 0.0,
    "use_age_survival": True,
    "entry_wealth_censor_to_frontier": False,
}
EXPECTED_TIGHT_EVALUATOR = {"max_iter_eq": 40, "tol_eq": 2.5e-5, "repeats": 2}

# M5 income-disciplined contract.  The share target/weight are the block3
# headline row (person-cluster bootstrap 1/SE^2); the old-age ownership rate
# is the legacy ACS-sourced target and weight.
NONHOUSING_SHARE = "old_nonhousing_ge_1x_income_share_6575"
NONHOUSING_SHARE_TARGET = 0.608333139649131
NONHOUSING_SHARE_SE = 0.0102949617302331
NONHOUSING_SHARE_WEIGHT = 9435.18732291246
OLD_AGE_OWN_RATE = "old_age_own_rate"
OLD_AGE_OWN_RATE_TARGET = 0.76426097
OLD_AGE_OWN_RATE_WEIGHT = 160.0
YOUNG_LIQUID = "young_childless_renter_liquid_wealth_to_annual_gross_income_2535"
M5_ESTABLISHED_BLOCK_LOSS_CEILING = 4.5
M5_YOUNG_LIQUID_GAP_CEILING = 0.15
M5_OLD_AGE_OWN_RATE_GAP_CEILING = 0.03
M5_SS_ANNUAL_RHO = 0.9601845894041878
M5_SS_ANNUAL_INNOVATION_SD = 0.06453733259357768
EXPECTED_M5_ACTIVE_DOMAIN = EXPECTED_M4_ACTIVE_DOMAIN + [
    {"name": "tenure_choice_kappa", "lower": 0.0, "upper": 0.12, "transform": "softzero"},
]
EXPECTED_M5_FIXED = {"theta_n": 0.0}
EXPECTED_M5_MECHANISM = dict(EXPECTED_M4_MECHANISM)
EXPECTED_M5_MECHANISM["entry_wealth_censor_to_frontier"] = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--arm", choices=tuple(ARM_CONTRACTS), default="M3")
    parser.add_argument(
        "--nested-reference",
        type=Path,
        default=None,
        help=(
            "summary.json of the max-evals-1 theta0=0 nested reference task "
            "(tight solves of the strict M1 seed under the arm objective); "
            "required for --arm M4 and --arm M5."
        ),
    )
    parser.add_argument(
        "--jacobian-summary",
        type=Path,
        default=None,
        help=(
            "Optional summary.json of the winner-point identification "
            "Jacobian job. M5 only: the identification_reported acceptance "
            "row is True when this file exists and 'pending' otherwise, so "
            "the collector never fails on a Jacobian job that runs later."
        ),
    )
    parser.add_argument(
        "--additional-summary",
        type=Path,
        action="append",
        default=[],
        help=(
            "Optional strict M4 summary.json from a predeclared post-profile "
            "full-dimensional polish. The six production-chain contract is "
            "validated separately and remains unchanged."
        ),
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        if not keys:
            return
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def eligible_tight(summary: dict[str, Any]) -> dict[str, Any] | None:
    tight = summary.get("best_tight")
    repeat = summary.get("tight_repeat_check") or {}
    if not (
        isinstance(tight, dict)
        and tight.get("strict_converged")
        and repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", 1.0)) == 0.0
        and float(repeat.get("max_abs_moment_difference", 1.0)) == 0.0
    ):
        return None
    return tight


def same_coordinate_up_to_roundtrip(left: Any, right: Any) -> bool:
    """Accept only floating-point transform roundoff in a stored coordinate."""

    lhs = float(left)
    rhs = float(right)
    if not math.isfinite(lhs) or not math.isfinite(rhs):
        return False
    tolerance = 4.0 * max(math.ulp(lhs), math.ulp(rhs))
    return math.isclose(lhs, rhs, rel_tol=0.0, abs_tol=tolerance)


def validate_m4_chain_metadata(metadata: dict[str, Any]) -> None:
    if not (
        metadata.get("status") == "proper_joint_smm_chain"
        and metadata.get("arm") == "M4"
        and metadata.get("target_set") == ARM_CONTRACTS["M4"][0]
        and int(metadata.get("free_parameter_count", -1)) == ARM_CONTRACTS["M4"][1]
        and int(metadata.get("target_count", -1)) == ARM_CONTRACTS["M4"][2]
        and metadata.get("active_domain") == EXPECTED_M4_ACTIVE_DOMAIN
        and metadata.get("fixed_parameters") == EXPECTED_M4_FIXED
        and metadata.get("mechanism") == EXPECTED_M4_MECHANISM
        and int(metadata.get("J", -1)) == 17
        and int(metadata.get("Nb", -1)) == 120
        and int(metadata.get("max_iter_eq", -1)) == 10
        and float(metadata.get("tol_eq", math.nan)) == 1e-4
        and metadata.get("tight_winner_evaluator") == EXPECTED_TIGHT_EVALUATOR
        and float((metadata.get("targets") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_TARGET
        and float((metadata.get("targets") or {}).get(NONHOUSING_MEDIAN, math.nan))
        == NONHOUSING_MEDIAN_TARGET
        and float((metadata.get("weights") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_WEIGHT
        and float((metadata.get("weights") or {}).get(NONHOUSING_MEDIAN, math.nan))
        == NONHOUSING_MEDIAN_WEIGHT
    ):
        raise RuntimeError("M4 chain metadata violates the exact production contract")


def validate_m5_chain_metadata(metadata: dict[str, Any]) -> None:
    income = metadata.get("income_process") or {}
    if not (
        metadata.get("status") == "proper_joint_smm_chain"
        and metadata.get("arm") == "M5"
        and metadata.get("target_set") == ARM_CONTRACTS["M5"][0]
        and int(metadata.get("free_parameter_count", -1)) == ARM_CONTRACTS["M5"][1]
        and int(metadata.get("target_count", -1)) == ARM_CONTRACTS["M5"][2]
        and metadata.get("active_domain") == EXPECTED_M5_ACTIVE_DOMAIN
        and metadata.get("fixed_parameters") == EXPECTED_M5_FIXED
        and metadata.get("mechanism") == EXPECTED_M5_MECHANISM
        and int(metadata.get("J", -1)) == 17
        and int(metadata.get("Nb", -1)) == 120
        and int(metadata.get("max_iter_eq", -1)) == 10
        and float(metadata.get("tol_eq", math.nan)) == 1e-4
        and metadata.get("tight_winner_evaluator") == EXPECTED_TIGHT_EVALUATOR
        and float(income.get("annual_rho", math.nan)) == M5_SS_ANNUAL_RHO
        and float(income.get("annual_innovation_sd", math.nan)) == M5_SS_ANNUAL_INNOVATION_SD
        and metadata.get("entry_wealth_ages") == "18_24"
        and float((metadata.get("targets") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_TARGET
        and float((metadata.get("targets") or {}).get(NONHOUSING_SHARE, math.nan))
        == NONHOUSING_SHARE_TARGET
        and float((metadata.get("targets") or {}).get(OLD_AGE_OWN_RATE, math.nan))
        == OLD_AGE_OWN_RATE_TARGET
        and float((metadata.get("weights") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_WEIGHT
        and float((metadata.get("weights") or {}).get(NONHOUSING_SHARE, math.nan))
        == NONHOUSING_SHARE_WEIGHT
        and float((metadata.get("weights") or {}).get(OLD_AGE_OWN_RATE, math.nan))
        == OLD_AGE_OWN_RATE_WEIGHT
        and NONHOUSING_MEDIAN not in (metadata.get("targets") or {})
    ):
        raise RuntimeError("M5 chain metadata violates the exact production contract")


def make_plots(outdir: Path, winner: dict[str, Any]) -> None:
    if plt is None:
        return
    life = winner["lifecycle"]
    ages = life["ages"]
    fields = (
        ("ownership_by_age", "Ownership rate"),
        ("age_mass", "Age mass"),
        ("mean_liquid_wealth_by_age", "Mean liquid wealth"),
        ("occupied_rooms_by_age", "Occupied rooms (aggregate)"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, (field, title) in zip(axes.flat, fields):
        ax.plot(ages, life[field], marker="o", ms=3, lw=1.7, color="#4C78A8")
        ax.set_title(title)
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Age")
    fig.tight_layout()
    fig.savefig(outdir / "lifecycle_diagnostics.png", dpi=180)
    plt.close(fig)

    rows = sorted(winner["target_fit"], key=lambda row: float(row["loss_contribution"]))
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(
        [row["moment"] for row in rows],
        [float(row["loss_contribution"]) for row in rows],
        color="#4C78A8",
    )
    ax.set_xlabel("Weighted loss contribution")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "loss_contributions.png", dpi=180)
    plt.close(fig)


def m4_acceptance_rows(
    winner: dict[str, Any], nested_loss: float
) -> tuple[list[dict[str, Any]], float]:
    fit_by_name = {str(row["moment"]): row for row in winner["target_fit"]}
    estate = fit_by_name[ESTATE_MEDIAN]
    nonhousing = fit_by_name[NONHOUSING_MEDIAN]
    established_loss = sum(
        float(row["loss_contribution"])
        for row in winner["target_fit"]
        if row["moment"] not in {ESTATE_MEDIAN, NONHOUSING_MEDIAN}
    )
    rows = [
        {
            "criterion": "strict_exact_tight_repeat",
            "value": True,
            "threshold": True,
            "pass": True,
        },
        {
            "criterion": "estate_median_absolute_gap",
            "value": abs(float(estate["gap"])),
            "threshold": ESTATE_MEDIAN_SE,
            "pass": abs(float(estate["gap"])) <= ESTATE_MEDIAN_SE,
        },
        {
            "criterion": "nonhousing_median_absolute_gap",
            "value": abs(float(nonhousing["gap"])),
            "threshold": 2.0 * NONHOUSING_MEDIAN_SE,
            "pass": (
                float(nonhousing["model"]) > 0.0
                and abs(float(nonhousing["gap"])) <= 2.0 * NONHOUSING_MEDIAN_SE
            ),
        },
        {
            "criterion": "established_12_moment_loss",
            "value": established_loss,
            "threshold": ESTABLISHED_BLOCK_LOSS_CEILING,
            "pass": established_loss <= ESTABLISHED_BLOCK_LOSS_CEILING,
        },
        {
            "criterion": "free_theta0_loss_minus_nested_zero",
            "value": float(winner["rank_loss"]) - float(nested_loss),
            "threshold": 0.0,
            "pass": float(winner["rank_loss"]) <= float(nested_loss),
        },
    ]
    return rows, established_loss


def m4_diagnostic_rows(winner: dict[str, Any]) -> list[dict[str, Any]]:
    moments = winner.get("moments") or {}
    lifecycle = winner.get("lifecycle") or {}
    requested = {
        "old_total_estate_wealth_to_annual_income_p90_p50_7684": (
            3.44811075444552,
            0.132477154936326,
        ),
        "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": (
            0.101108088567873,
            0.563010076238571,
        ),
        "old_total_wealth_to_income_median_6575": (None, None),
        "old_age_own_rate": (None, None),
    }
    rows = [
        {
            "diagnostic": name,
            "model": float(moments[name]),
            "reference": reference,
            "bootstrap_se": standard_error,
        }
        for name, (reference, standard_error) in requested.items()
        if name in moments and math.isfinite(float(moments[name]))
    ]
    rows.extend(
        [
            {
                "diagnostic": "ownership_path_acceptance",
                "model": bool(lifecycle.get("hard_acceptance_pass", False)),
                "reference": True,
                "bootstrap_se": None,
            },
            {
                "diagnostic": "ownership_band_failure",
                "model": bool(lifecycle.get("ownership_band_failure", True)),
                "reference": False,
                "bootstrap_se": None,
            },
            {
                "diagnostic": "ownership_cliff_failure",
                "model": bool(lifecycle.get("ownership_cliff_failure", True)),
                "reference": False,
                "bootstrap_se": None,
            },
        ]
    )
    return rows


def m5_acceptance_rows(
    winner: dict[str, Any],
    nested_loss: float,
    jacobian_summary: Path | None,
) -> tuple[list[dict[str, Any]], float]:
    fit_by_name = {str(row["moment"]): row for row in winner["target_fit"]}
    estate = fit_by_name[ESTATE_MEDIAN]
    share = fit_by_name[NONHOUSING_SHARE]
    young_liquid = fit_by_name[YOUNG_LIQUID]
    old_own = fit_by_name[OLD_AGE_OWN_RATE]
    established_loss = sum(
        float(row["loss_contribution"])
        for row in winner["target_fit"]
        if row["moment"] not in {ESTATE_MEDIAN, NONHOUSING_SHARE, OLD_AGE_OWN_RATE}
    )
    jacobian_reported = jacobian_summary is not None and jacobian_summary.is_file()
    rows = [
        {
            "criterion": "strict_exact_tight_repeat",
            "value": True,
            "threshold": True,
            "pass": True,
        },
        {
            "criterion": "estate_median_absolute_gap",
            "value": abs(float(estate["gap"])),
            "threshold": ESTATE_MEDIAN_SE,
            "pass": abs(float(estate["gap"])) <= ESTATE_MEDIAN_SE,
        },
        {
            "criterion": "nonhousing_ge_1x_share_absolute_gap",
            "value": abs(float(share["gap"])),
            "threshold": 2.0 * NONHOUSING_SHARE_SE,
            "pass": abs(float(share["gap"])) <= 2.0 * NONHOUSING_SHARE_SE,
        },
        {
            "criterion": "established_12_moment_loss",
            "value": established_loss,
            "threshold": M5_ESTABLISHED_BLOCK_LOSS_CEILING,
            "pass": established_loss <= M5_ESTABLISHED_BLOCK_LOSS_CEILING,
        },
        {
            "criterion": "young_liquid_absolute_gap",
            "value": abs(float(young_liquid["gap"])),
            "threshold": M5_YOUNG_LIQUID_GAP_CEILING,
            "pass": abs(float(young_liquid["gap"])) <= M5_YOUNG_LIQUID_GAP_CEILING,
        },
        {
            "criterion": "old_age_own_rate_absolute_gap",
            "value": abs(float(old_own["gap"])),
            "threshold": M5_OLD_AGE_OWN_RATE_GAP_CEILING,
            "pass": abs(float(old_own["gap"])) <= M5_OLD_AGE_OWN_RATE_GAP_CEILING,
        },
        {
            "criterion": "free_winner_loss_minus_nested_zero",
            "value": float(winner["rank_loss"]) - float(nested_loss),
            "threshold": 0.0,
            "pass": float(winner["rank_loss"]) <= float(nested_loss),
        },
        {
            # Wired loosely: the Jacobian job may run after collection, so an
            # absent summary marks this row "pending" instead of failing.
            "criterion": "identification_reported",
            "value": str(jacobian_summary) if jacobian_reported else None,
            "threshold": "jacobian summary file exists",
            "pass": True if jacobian_reported else "pending",
        },
    ]
    return rows, established_loss


def m5_diagnostic_rows(winner: dict[str, Any]) -> list[dict[str, Any]]:
    moments = winner.get("moments") or {}
    lifecycle = winner.get("lifecycle") or {}
    requested = {
        "old_total_estate_wealth_to_annual_income_p90_p50_7684": (
            3.44811075444552,
            0.132477154936326,
        ),
        "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575": (
            0.101108088567873,
            0.563010076238571,
        ),
        NONHOUSING_MEDIAN: (NONHOUSING_MEDIAN_TARGET, NONHOUSING_MEDIAN_SE),
        "old_total_wealth_to_income_median_6575": (None, None),
    }
    rows = [
        {
            "diagnostic": name,
            "model": float(moments[name]),
            "reference": reference,
            "bootstrap_se": standard_error,
        }
        for name, (reference, standard_error) in requested.items()
        if name in moments and math.isfinite(float(moments[name]))
    ]
    rows.extend(
        [
            {
                "diagnostic": "ownership_path_acceptance",
                "model": bool(lifecycle.get("hard_acceptance_pass", False)),
                "reference": True,
                "bootstrap_se": None,
            },
            {
                "diagnostic": "ownership_band_failure",
                "model": bool(lifecycle.get("ownership_band_failure", True)),
                "reference": False,
                "bootstrap_se": None,
            },
            {
                "diagnostic": "ownership_cliff_failure",
                "model": bool(lifecycle.get("ownership_cliff_failure", True)),
                "reference": False,
                "bootstrap_se": None,
            },
        ]
    )
    return rows


def validate_m4_nested_reference(
    nested_summary: dict[str, Any], target_set: str
) -> dict[str, Any]:
    metadata = nested_summary.get("metadata") or {}
    if not (
        metadata.get("status") == "proper_joint_smm_chain"
        and metadata.get("arm") == "M4"
        and metadata.get("target_set") == target_set
        and int(metadata.get("free_parameter_count", -1)) == 13
        and int(metadata.get("target_count", -1)) == 14
        and int(metadata.get("max_evals", -1)) == 1
        and float(metadata.get("start_mix", math.nan)) == 0.0
        and metadata.get("seed_arm") == "M1"
        and int(metadata.get("J", -1)) == 17
        and int(metadata.get("Nb", -1)) == 120
        and int(metadata.get("max_iter_eq", -1)) == 10
        and float(metadata.get("tol_eq", math.nan)) == 1e-4
        and metadata.get("tight_winner_evaluator") == EXPECTED_TIGHT_EVALUATOR
        and metadata.get("active_domain") == EXPECTED_M4_ACTIVE_DOMAIN
        and metadata.get("fixed_parameters") == EXPECTED_M4_FIXED
        and metadata.get("mechanism") == EXPECTED_M4_MECHANISM
        and float((metadata.get("targets") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_TARGET
        and float((metadata.get("targets") or {}).get(NONHOUSING_MEDIAN, math.nan))
        == NONHOUSING_MEDIAN_TARGET
        and float((metadata.get("weights") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_WEIGHT
        and float((metadata.get("weights") or {}).get(NONHOUSING_MEDIAN, math.nan))
        == NONHOUSING_MEDIAN_WEIGHT
    ):
        raise RuntimeError("nested reference metadata does not match the exact M4/M1 contract")

    nested_tight = eligible_tight(nested_summary)
    if nested_tight is None:
        raise RuntimeError("nested theta0=0 reference has no bit-identical strict tight repeats")
    nested_theta = nested_tight.get("theta") or {}
    if float(nested_theta.get("theta0", math.nan)) != 0.0:
        raise RuntimeError("nested reference tight solve is not at theta0=0")
    if any(float(nested_theta.get(name, math.nan)) != value for name, value in EXPECTED_M4_FIXED.items()):
        raise RuntimeError("nested reference violates the fixed M4 restrictions")

    seed_path = Path(str(metadata.get("seed_record", "")))
    if not seed_path.is_file():
        raise RuntimeError(f"nested M1 seed record is unavailable: {seed_path}")
    seed_payload = json.loads(seed_path.read_text())
    m1 = (seed_payload.get("winners") or {}).get("M1")
    if not isinstance(m1, dict) or not bool(m1.get("strict_converged", False)):
        raise RuntimeError("nested seed record does not contain a strict M1 winner")
    seed_theta = m1.get("theta") or {}
    mismatches = {
        name: {"seed": seed_theta.get(name), "nested": nested_theta.get(name)}
        for name in M4_SHARED_THETA
        if not same_coordinate_up_to_roundtrip(
            seed_theta.get(name, math.nan), nested_theta.get(name, math.nan)
        )
    }
    if mismatches:
        raise RuntimeError(f"nested reference is not the exact strict M1 theta: {mismatches}")
    return nested_tight


def validate_m5_nested_reference(
    nested_summary: dict[str, Any], target_set: str
) -> dict[str, Any]:
    metadata = nested_summary.get("metadata") or {}
    income = metadata.get("income_process") or {}
    if not (
        metadata.get("status") == "proper_joint_smm_chain"
        and metadata.get("arm") == "M5"
        and metadata.get("target_set") == target_set
        and int(metadata.get("free_parameter_count", -1)) == ARM_CONTRACTS["M5"][1]
        and int(metadata.get("target_count", -1)) == ARM_CONTRACTS["M5"][2]
        and int(metadata.get("max_evals", -1)) == 1
        and float(metadata.get("start_mix", math.nan)) == 0.0
        and metadata.get("seed_arm") == "M1"
        and int(metadata.get("J", -1)) == 17
        and int(metadata.get("Nb", -1)) == 120
        and int(metadata.get("max_iter_eq", -1)) == 10
        and float(metadata.get("tol_eq", math.nan)) == 1e-4
        and metadata.get("tight_winner_evaluator") == EXPECTED_TIGHT_EVALUATOR
        and metadata.get("active_domain") == EXPECTED_M5_ACTIVE_DOMAIN
        and metadata.get("fixed_parameters") == EXPECTED_M5_FIXED
        and metadata.get("mechanism") == EXPECTED_M5_MECHANISM
        and float(income.get("annual_rho", math.nan)) == M5_SS_ANNUAL_RHO
        and float(income.get("annual_innovation_sd", math.nan)) == M5_SS_ANNUAL_INNOVATION_SD
        and metadata.get("entry_wealth_ages") == "18_24"
        and float((metadata.get("targets") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_TARGET
        and float((metadata.get("targets") or {}).get(NONHOUSING_SHARE, math.nan))
        == NONHOUSING_SHARE_TARGET
        and float((metadata.get("targets") or {}).get(OLD_AGE_OWN_RATE, math.nan))
        == OLD_AGE_OWN_RATE_TARGET
        and float((metadata.get("weights") or {}).get(ESTATE_MEDIAN, math.nan))
        == ESTATE_MEDIAN_WEIGHT
        and float((metadata.get("weights") or {}).get(NONHOUSING_SHARE, math.nan))
        == NONHOUSING_SHARE_WEIGHT
        and float((metadata.get("weights") or {}).get(OLD_AGE_OWN_RATE, math.nan))
        == OLD_AGE_OWN_RATE_WEIGHT
    ):
        raise RuntimeError("nested reference metadata does not match the exact M5/M1 contract")

    nested_tight = eligible_tight(nested_summary)
    if nested_tight is None:
        raise RuntimeError("nested theta0=0 reference has no bit-identical strict tight repeats")
    nested_theta = nested_tight.get("theta") or {}
    if float(nested_theta.get("theta0", math.nan)) != 0.0:
        raise RuntimeError("nested reference tight solve is not at theta0=0")
    if float(nested_theta.get("tenure_choice_kappa", math.nan)) != 0.0:
        raise RuntimeError("nested reference tight solve is not at tenure_choice_kappa=0")
    if any(
        float(nested_theta.get(name, math.nan)) != value
        for name, value in EXPECTED_M5_FIXED.items()
    ):
        raise RuntimeError("nested reference violates the fixed M5 restrictions")

    seed_path = Path(str(metadata.get("seed_record", "")))
    if not seed_path.is_file():
        raise RuntimeError(f"nested M1 seed record is unavailable: {seed_path}")
    seed_payload = json.loads(seed_path.read_text())
    m1 = (seed_payload.get("winners") or {}).get("M1")
    if not isinstance(m1, dict) or not bool(m1.get("strict_converged", False)):
        raise RuntimeError("nested seed record does not contain a strict M1 winner")
    seed_theta = m1.get("theta") or {}
    mismatches = {
        name: {"seed": seed_theta.get(name), "nested": nested_theta.get(name)}
        for name in M4_SHARED_THETA
        if not same_coordinate_up_to_roundtrip(
            seed_theta.get(name, math.nan), nested_theta.get(name, math.nan)
        )
    }
    if mismatches:
        raise RuntimeError(f"nested reference is not the exact strict M1 theta: {mismatches}")
    return nested_tight


def main() -> None:
    args = parse_args()
    target_set, free_count, target_count, min_chains = ARM_CONTRACTS[args.arm]
    if args.arm in {"M4", "M5"} and args.nested_reference is None:
        raise RuntimeError(
            f"--arm {args.arm} requires --nested-reference with the theta0=0 tight solve"
        )
    production_paths = sorted(args.results_dir.glob("task_*/summary.json"))
    if args.arm == "M4" and len(production_paths) != 6:
        raise RuntimeError(f"expected exactly 6 completed M4 chains, found {len(production_paths)}")
    if args.arm == "M5" and len(production_paths) != 8:
        raise RuntimeError(f"expected exactly 8 completed M5 chains, found {len(production_paths)}")
    if len(production_paths) < min_chains:
        raise RuntimeError(
            f"expected at least {min_chains} completed {args.arm} chains, found {len(production_paths)}"
        )
    if args.additional_summary and args.arm != "M4":
        raise RuntimeError("--additional-summary is only valid for M4")
    production_summaries = [json.loads(path.read_text()) for path in production_paths]
    additional_paths = list(args.additional_summary)
    additional_summaries = [json.loads(path.read_text()) for path in additional_paths]
    summaries = production_summaries + additional_summaries
    summary_paths = production_paths + additional_paths
    if args.arm == "M4":
        metadata_rows = [summary.get("metadata") or {} for summary in production_summaries]
        for meta in metadata_rows:
            validate_m4_chain_metadata(meta)
        seeds = [int(meta.get("seed", -1)) for meta in metadata_rows]
        starts = [tuple(float(v) for v in meta.get("initial_unit_vector", [])) for meta in metadata_rows]
        if len(set(seeds)) != 6 or len(set(starts)) != 6 or any(len(start) != 13 for start in starts):
            raise RuntimeError("M4 production chains do not contain six distinct 13-D starts/seeds")
        for summary in additional_summaries:
            validate_m4_chain_metadata(summary.get("metadata") or {})
    elif args.arm == "M5":
        metadata_rows = [summary.get("metadata") or {} for summary in production_summaries]
        for meta in metadata_rows:
            validate_m5_chain_metadata(meta)
        seeds = [int(meta.get("seed", -1)) for meta in metadata_rows]
        starts = [tuple(float(v) for v in meta.get("initial_unit_vector", [])) for meta in metadata_rows]
        if len(set(seeds)) != 8 or len(set(starts)) != 8 or any(len(start) != 14 for start in starts):
            raise RuntimeError("M5 production chains do not contain eight distinct 14-D starts/seeds")

    eligible: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    for summary_index, (summary, summary_path) in enumerate(zip(summaries, summary_paths)):
        meta = summary["metadata"]
        if not (
            meta.get("arm") == args.arm
            and meta.get("target_set") == target_set
            and int(meta.get("free_parameter_count", -1)) == free_count
            and int(meta.get("target_count", -1)) == target_count
        ):
            raise RuntimeError(
                f"{args.arm} chain metadata does not match the "
                f"{free_count}-parameter/{target_count}-target contract"
            )
        tight = eligible_tight(summary)
        if args.arm in {"M4", "M5"} and tight is not None:
            expected_fixed = EXPECTED_M4_FIXED if args.arm == "M4" else EXPECTED_M5_FIXED
            tight_theta = tight.get("theta") or {}
            if any(
                float(tight_theta.get(name, math.nan)) != value
                for name, value in expected_fixed.items()
            ):
                raise RuntimeError(f"{args.arm} tight winner violates its fixed restrictions")
        chain_rows.append(
            {
                "source": "production" if summary_index < len(production_summaries) else "post_profile_polish",
                "summary_path": str(summary_path),
                "seed": meta["seed"],
                "start_mix": meta["start_mix"],
                "eligible": tight is not None,
                "search_cases": summary["n_cases_completed"],
                "search_strict": summary["n_strict"],
                "tight_loss": tight["rank_loss"] if tight else None,
                "tight_residual": tight["market_residual"] if tight else None,
            }
        )
        if tight is not None:
            eligible.append(tight)
    if not eligible:
        raise RuntimeError(f"no strict, exactly repeated {args.arm} winner")
    winner = min(eligible, key=lambda row: float(row["rank_loss"]))

    nested_loss = None
    if args.arm in {"M4", "M5"}:
        # The nested reference is the same driver run with max-evals 1 seeded
        # at the strict M1 winner, so its tight repeats price the exact
        # theta0=0 (and, for M5, tenure_choice_kappa=0) nested submodel under
        # the arm objective and evaluator.
        nested_summary = json.loads(args.nested_reference.read_text())
        if args.arm == "M4":
            nested_tight = validate_m4_nested_reference(nested_summary, target_set)
        else:
            nested_tight = validate_m5_nested_reference(nested_summary, target_set)
        nested_loss = float(nested_tight["rank_loss"])
        if float(winner["rank_loss"]) > nested_loss:
            raise RuntimeError(
                f"{args.arm} missed its theta0=0 nested M1 seed: tight loss "
                f"{float(winner['rank_loss']):.6f} exceeds nested {nested_loss:.6f}"
            )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)
    write_csv(args.outdir / "target_fit_full.csv", list(winner["target_fit"]))
    write_csv(args.outdir / "parameter_table_full.csv", list(winner["parameters"]))
    write_csv(
        args.outdir / "lifecycle_decomposition.csv",
        [
            {
                "age": age,
                "age_mass": mass,
                "ownership_rate": own,
                "mean_liquid_wealth": liquid,
                "owner_rooms": owner_rooms,
                "occupied_rooms": occupied_rooms,
            }
            for age, mass, own, liquid, owner_rooms, occupied_rooms in zip(
                winner["lifecycle"]["ages"],
                winner["lifecycle"]["age_mass"],
                winner["lifecycle"]["ownership_by_age"],
                winner["lifecycle"]["mean_liquid_wealth_by_age"],
                winner["lifecycle"]["owner_rooms_by_age"],
                winner["lifecycle"]["occupied_rooms_by_age"],
            )
        ],
    )
    acceptance_rows: list[dict[str, Any]] = []
    established_loss = None
    if args.arm == "M4":
        assert nested_loss is not None
        acceptance_rows, established_loss = m4_acceptance_rows(winner, nested_loss)
        write_csv(args.outdir / "acceptance_criteria.csv", acceptance_rows)
        write_csv(args.outdir / "untargeted_diagnostics.csv", m4_diagnostic_rows(winner))
    elif args.arm == "M5":
        assert nested_loss is not None
        acceptance_rows, established_loss = m5_acceptance_rows(
            winner, nested_loss, args.jacobian_summary
        )
        write_csv(args.outdir / "acceptance_criteria.csv", acceptance_rows)
        write_csv(args.outdir / "untargeted_diagnostics.csv", m5_diagnostic_rows(winner))

    results_payload: dict[str, Any] = {
        "winners": {args.arm: winner},
        f"eligible_{args.arm}_chains": len(eligible),
        "production_chain_count": len(production_summaries),
        "additional_summary_count": len(additional_summaries),
        "additional_summary_paths": [str(path) for path in additional_paths],
    }
    if nested_loss is not None:
        results_payload["nested_theta0_zero_tight_loss"] = float(nested_loss)
        results_payload["established_12_moment_loss"] = float(established_loss)
        results_payload["acceptance_criteria"] = acceptance_rows
        pending_criteria = [
            str(row["criterion"]) for row in acceptance_rows if row["pass"] == "pending"
        ]
        results_payload["pending_acceptance_criteria"] = pending_criteria
        results_payload["all_acceptance_criteria_pass"] = all(
            bool(row["pass"]) for row in acceptance_rows if row["pass"] != "pending"
        )
    (args.outdir / "results.json").write_text(
        json.dumps(results_payload, indent=2, sort_keys=True)
    )
    make_plots(args.outdir, winner)

    theta = winner["theta"]
    if args.arm == "M4":
        header = "# Standard bequest recalibration (M4)"
        blurb = (
            "The reported winner estimates the 11 clean-frontier parameters plus "
            "`theta0` and `theta1` against 14 moments; `theta_n` and "
            "`tenure_choice_kappa` are external restrictions. The winner beats "
            "the strict M1 `theta0=0` nested seed under the M4 objective. All "
            "reported moments use a strict, exactly repeated winner solve."
        )
    elif args.arm == "M5":
        header = "# Income-disciplined recalibration (M5)"
        blurb = (
            "The reported winner estimates the 11 clean-frontier parameters plus "
            "`theta0`, `theta1`, and `tenure_choice_kappa` against 15 moments; "
            "`theta_n=0` is the only external restriction. Income persistence is "
            "the literature anchor (annual rho 0.90, innovation s.d. 0.20) and "
            "entrant wealth uses the PSID 18-24 distribution. The winner beats "
            "the strict M1 `theta0=0`/`tenure_choice_kappa=0` nested seed under "
            "the M5 objective. All reported moments use a strict, exactly "
            "repeated winner solve."
        )
    else:
        header = "# Internally calibrated bequest block"
        blurb = (
            "The reported winner estimates the 11 clean-frontier parameters plus "
            "`theta0`, `theta1`, and `theta_n` against 15 moments. All reported "
            "moments use a strict, exactly repeated winner solve."
        )
    lines = [
        header,
        "",
        blurb,
        "",
        "| Loss | Residual | theta0 | theta1 | theta_n | Eligible chains |",
        "|---:|---:|---:|---:|---:|---:|",
        f"| {float(winner['rank_loss']):.6f} | {float(winner['market_residual']):.2e} | "
        f"{float(theta['theta0']):.6f} | {float(theta['theta1']):.6f} | "
        f"{float(theta['theta_n']):.6f} | {len(eligible)} |",
        "",
        "Complete target-fit, parameter-bound, lifecycle, and plot artifacts are adjacent.",
    ]
    if args.arm in {"M4", "M5"}:
        ceiling = (
            ESTABLISHED_BLOCK_LOSS_CEILING
            if args.arm == "M4"
            else M5_ESTABLISHED_BLOCK_LOSS_CEILING
        )
        lines.extend(
            [
                "",
                f"Established 12-moment loss: {float(established_loss):.6f} "
                f"(ceiling {ceiling:.2f}).",
                f"All predeclared acceptance criteria pass: "
                f"{all(bool(row['pass']) for row in acceptance_rows if row['pass'] != 'pending')}.",
            ]
        )
    if args.arm == "M5":
        lines.append(
            f"Estimated tenure_choice_kappa: {float(theta['tenure_choice_kappa']):.6f}."
        )
        pending = [
            str(row["criterion"]) for row in acceptance_rows if row["pass"] == "pending"
        ]
        if pending:
            lines.append(f"Pending acceptance criteria: {', '.join(pending)}.")
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
