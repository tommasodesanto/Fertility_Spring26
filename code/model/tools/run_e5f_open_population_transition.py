#!/usr/bin/env python3
"""Cohort-accounting transition for the sequential child-room-floor model.

The script makes no change to the household state space or stationary solver.
It starts from a stationary sequential-model distribution, applies a permanent
change in the child-preference shifter, and carries the full age, wealth,
income, tenure, parity, and children-at-home distribution through calendar
time. At each date, households treat the current price and primitives as
permanent and the contemporaneous housing market clears against the initial
housing stock. This is therefore a transparent temporary-equilibrium
transition diagnostic, not a perfect-foresight asset-pricing exercise.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
import math
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import audit_closed_reproductive_closure as closure_audit
import run_dynamic_population_transition as calendar


DEFAULT_SOURCE = closure_audit.E5F_FLOOR_SOURCE
DEFAULT_OUTDIR = ROOT / "output/model/e5f_floor_open_population_transition"
CENSUS_HH3_URL = (
    "https://www2.census.gov/programs-surveys/demo/tables/families/"
    "time-series/households/hh3.xls"
)
CENSUS_HH3_HOUSEHOLDS_THOUSANDS = {
    2007: 116011.0,
    2011: 119927.0,
    2015: 124587.0,
    2019: 128579.0,
    2023: 131434.0,
}
ACS_NATIONAL_HEAD_AGE_SOURCE = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
)
ACS_NATIONAL_HEAD_AGE_SOURCE_SHA256 = (
    "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e"
)
ACS_NATIONAL_HEAD_AGE_SHARES_18_85 = {
    2007: (
        0.01662266202182001, 0.04529676719571153, 0.06359694782099141,
        0.06696559869635438, 0.07884514768615880, 0.08072893564675242,
        0.08968295548613242, 0.09115500500154539, 0.08684345526254360,
        0.07812470544072851, 0.07104778040230673, 0.05661696446968036,
        0.04505232136554861, 0.03889108604181299, 0.03548381308904295,
        0.03128817505209684, 0.02375767932077307,
    ),
    2011: (
        0.01423640057870632, 0.04015597541289629, 0.05944231037589900,
        0.06829743892768564, 0.06969370947726934, 0.07775638667739104,
        0.07863758615588470, 0.08591835898663702, 0.08895820986694107,
        0.08466006473759148, 0.07633661295406242, 0.06975907935230928,
        0.05367679625793782, 0.04331476984142344, 0.03516810939479943,
        0.02998756910537614, 0.02400062189718957,
    ),
    2015: (
        0.01237826965969259, 0.03850574185466547, 0.05703943274201573,
        0.06801644249983209, 0.07105444888041931, 0.06998319536695549,
        0.07633662912233560, 0.07639420525997405, 0.08415207247359373,
        0.08484820293023589, 0.08099594893356551, 0.07322657866496309,
        0.06532914452875623, 0.05029551243914543, 0.03838511193415870,
        0.02999548094268535, 0.02306358176700575,
    ),
    2019: (
        0.01204386212620710, 0.03672383072129470, 0.05910831022458462,
        0.06727556723592544, 0.07113343807388117, 0.07127626076428463,
        0.06827752862343339, 0.07300446221666929, 0.07298973943669219,
        0.07877541510471971, 0.08171188110840960, 0.07839923886409916,
        0.06788946863719994, 0.06041097476218209, 0.04486839220718734,
        0.03290898088489887, 0.02320264900833074,
    ),
    2023: (
        0.01269617145483895, 0.03962038556026601, 0.05765140549135367,
        0.06890576989592977, 0.07222893075363992, 0.07333366867966330,
        0.07215408288932150, 0.06656022972977625, 0.07142826081255550,
        0.06959182196947911, 0.07656515710357487, 0.07798978045201599,
        0.07085160357313883, 0.06165837774909073, 0.05063149245400102,
        0.03500355974319133, 0.02312930168816325,
    ),
}
ACS_NATIONAL_HEAD_AGE_RECEIPTS = {
    2007: {
        "sample": 200701,
        "head_records_all_ages": 1_174_259,
        "head_weight_all_ages": 112_376_192.0,
        "age_18_85_share_of_all_heads": 0.9785258874050475,
    },
    2011: {
        "sample": 201101,
        "head_records_all_ages": 1_204_817,
        "head_weight_all_ages": 114_990_599.0,
        "age_18_85_share_of_all_heads": 0.9746011584825295,
    },
    2015: {
        "sample": 201501,
        "head_records_all_ages": 1_226_681,
        "head_weight_all_ages": 118_204_832.0,
        "age_18_85_share_of_all_heads": 0.9729968145464646,
    },
    2019: {
        "sample": 201901,
        "head_records_all_ages": 1_276_656,
        "head_weight_all_ages": 122_797_400.0,
        "age_18_85_share_of_all_heads": 0.9723886417790605,
    },
    2023: {
        "sample": 202301,
        "head_records_all_ages": 1_342_971,
        "head_weight_all_ages": 131_325_888.0,
        "age_18_85_share_of_all_heads": 0.9753331117776260,
    },
}
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--periods", type=int, default=25)
    parser.add_argument(
        "--old-psi-child",
        type=float,
        default=0.10,
        help="Child-preference shifter in the initial stationary equilibrium.",
    )
    parser.add_argument(
        "--new-psi-child",
        type=float,
        default=0.0,
        help="Permanent post-shock child-preference shifter.",
    )
    parser.add_argument(
        "--preference-transition-periods",
        type=int,
        default=0,
        help=(
            "Number of model periods over which psi_child moves linearly from "
            "the old to the new value. Zero preserves the immediate permanent shock."
        ),
    )
    parser.add_argument(
        "--historical-start-year",
        type=int,
        default=None,
        help=(
            "If supplied, download the paper-facing FRED comparison series and "
            "index the model's date 0 to this calendar year."
        ),
    )
    parser.add_argument(
        "--housing-supply-mode",
        choices=("fixed-stock", "static-elastic"),
        default="fixed-stock",
        help=(
            "Fixed-stock is the short-run benchmark. Static-elastic is the "
            "Oswald-style sequence of contemporaneous housing-supply equilibria."
        ),
    )
    parser.add_argument(
        "--housing-supply-elasticity",
        type=float,
        default=None,
        help=(
            "Optional diagnostic override for the static housing-supply "
            "elasticity; the retained model value is used when omitted."
        ),
    )
    parser.add_argument(
        "--outside-origin-entry-share",
        type=float,
        default=0.169,
        help="Initial steady-state share of entrant households supplied from outside.",
    )
    parser.add_argument(
        "--market-tol",
        type=float,
        default=2e-4,
        help=(
            "Relative clearing tolerance. The transition distribution is atomic "
            "on the deterministic tenure/wealth grid, so its demand schedule can "
            "jump across the exact root."
        ),
    )
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument(
        "--renewal-clock",
        choices=("birth-vintage", "household-state"),
        default="birth-vintage",
        help=(
            "Birth-vintage is the paper transition: births wait four complete "
            "model periods before supplying entrant households. Household-state "
            "reuses the model's memoryless child-exit flow as a diagnostic."
        ),
    )
    parser.add_argument(
        "--birth-to-entry-delay-periods",
        type=int,
        default=4,
        help="Completed four-year child stages between birth and adult entry.",
    )
    parser.add_argument(
        "--initial-birth-pipeline-multiplier",
        type=float,
        default=1.0,
        help=(
            "Diagnostic scale factor on the locally born cohorts already in the "
            "birth-to-entry pipeline at date 0. Future cohorts remain model-generated."
        ),
    )
    parser.add_argument(
        "--initial-age-profile",
        choices=("stationary", "acs-national-heads-2007"),
        default="stationary",
        help=(
            "Optionally reweight the stationary household distribution to the "
            "national 2007 ACS four-year householder-age shares while preserving "
            "every model state conditional on age."
        ),
    )
    parser.add_argument(
        "--household-count-path",
        choices=("none", "census-hh3-2007-2023"),
        default="none",
        help=(
            "Apply a dated age-specific household-formation/migration residual "
            "that matches Census HH-3 household totals and national ACS four-year "
            "householder-age shares through 2023."
        ),
    )
    parser.add_argument(
        "--policy-case",
        choices=("none", "dependent-child-ltv95"),
        default="none",
        help="Optional tenure policy activated at --policy-start-period.",
    )
    parser.add_argument(
        "--policy-start-period",
        type=int,
        default=4,
        help="First model date at which the optional tenure policy is active.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run at most four calendar periods after the exact baseline gates.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def preference_shifter_at_date(
    period: int,
    old_value: float,
    new_value: float,
    transition_periods: int,
) -> float:
    """Return the dated fertility-preference shifter.

    With a positive transition length, date 0 remains the old steady state and
    the new value is reached exactly at ``period == transition_periods``.
    """
    if transition_periods <= 0:
        return float(new_value)
    weight = min(max(float(period) / float(transition_periods), 0.0), 1.0)
    return float(old_value + weight * (new_value - old_value))


def acs_2007_age_reweight_diagnostic(
    stationary_age_mass: np.ndarray,
    ages: np.ndarray,
    stationary_entry_flow: float,
    periods: int = 4,
    period_years: float = 4.0,
) -> dict[str, Any]:
    """Advance the ACS-reweighted age stock at the stationary entrant flow."""
    stationary = np.asarray(stationary_age_mass, dtype=float).reshape(-1)
    age_grid = np.asarray(ages, dtype=float).reshape(-1)
    if stationary.shape != age_grid.shape or stationary.size < 2:
        raise ValueError("Age masses and age labels must be conformable vectors")
    if np.any(~np.isfinite(stationary)) or np.any(stationary <= 0.0):
        raise ValueError("Stationary age masses must be finite and strictly positive")
    if periods < 1:
        raise ValueError("The age-reweight horizon must be positive")
    if period_years <= 0.0:
        raise ValueError("The model period length must be positive")

    initial_total = float(np.sum(stationary))
    target_shares = np.asarray(
        ACS_NATIONAL_HEAD_AGE_SHARES_18_85[2007], dtype=float
    )
    if target_shares.shape != stationary.shape:
        raise RuntimeError(
            "The national ACS age bins do not match the model age grid"
        )
    if not math.isclose(float(np.sum(target_shares)), 1.0, rel_tol=0.0, abs_tol=2e-14):
        raise RuntimeError("National ACS age shares do not sum to one")
    reweighted = initial_total * target_shares
    if not math.isclose(float(np.sum(reweighted)), initial_total, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("ACS age reweighting changed date-0 total mass")

    survival = stationary[1:] / stationary[:-1]
    if np.any(survival < -1e-14) or np.any(survival > 1.0 + 1e-10):
        raise RuntimeError("Stationary age masses imply an invalid survival schedule")
    if not math.isclose(
        float(stationary_entry_flow),
        float(stationary[0]),
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        raise RuntimeError("Stationary entry flow does not match the youngest age mass")

    current = reweighted.copy()
    path = [
        {
            "period": 0,
            "adult_household_mass": float(np.sum(current)),
            "population_index": 1.0,
        }
    ]
    for period in range(1, periods + 1):
        next_mass = np.zeros_like(current)
        next_mass[0] = float(stationary_entry_flow)
        next_mass[1:] = survival * current[:-1]
        current = next_mass
        path.append(
            {
                "period": period,
                "adult_household_mass": float(np.sum(current)),
                "population_index": float(np.sum(current)) / initial_total,
            }
        )

    return {
        "status": "diagnostic_only",
        "source": "IPUMS USA ACS microdata extract",
        "source_path": str(ACS_NATIONAL_HEAD_AGE_SOURCE),
        "source_sha256": ACS_NATIONAL_HEAD_AGE_SOURCE_SHA256,
        "source_year": 2007,
        "units": "HHWT-weighted household heads ages 18--85",
        "sample": (
            "national ACS; GQ in {1,2}; PERNUM=1; RELATE=1; HHWT>0; "
            "age 18--85"
        ),
        "mapping": (
            "exact closed integer bins [18+4j, 21+4j], j=0,...,16; "
            "every non-age model state is preserved conditional on age"
        ),
        "source_receipt": ACS_NATIONAL_HEAD_AGE_RECEIPTS[2007],
        "future_entry_rule": "old steady-state entrant flow held fixed",
        "horizon_periods": periods,
        "period_years": period_years,
        "horizon_calendar_year": 2007 + periods * period_years,
        "model_age_grid": age_grid.tolist(),
        "stationary_age_mass": stationary.tolist(),
        "reweighted_initial_age_mass": reweighted.tolist(),
        "age_survival": survival.tolist(),
        "path": path,
        "population_index_at_horizon": path[-1]["population_index"],
    }


def reweight_distribution_to_acs_2007_ages(
    g_pre: np.ndarray,
    diagnostic: dict[str, Any],
) -> np.ndarray:
    """Apply the diagnostic's age weights, preserving every within-age state."""
    out = np.asarray(g_pre, dtype=float).copy()
    stationary = np.asarray(diagnostic["stationary_age_mass"], dtype=float)
    target = np.asarray(diagnostic["reweighted_initial_age_mass"], dtype=float)
    if out.shape[3] != stationary.size or stationary.shape != target.shape:
        raise ValueError("The ACS age diagnostic does not match the KFE age axis")
    factors = target / stationary
    out *= factors.reshape((1, 1, 1, -1, 1, 1, 1))
    if not math.isclose(
        float(np.sum(out)),
        float(np.sum(g_pre)),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("ACS age reweighting changed total household mass")
    return out


def census_household_model_age_levels() -> dict[int, float]:
    """Convert published HH-3 totals to heads in the model's ages 18--85."""
    return {
        year: float(total)
        * float(
            ACS_NATIONAL_HEAD_AGE_RECEIPTS[year][
                "age_18_85_share_of_all_heads"
            ]
        )
        for year, total in CENSUS_HH3_HOUSEHOLDS_THOUSANDS.items()
    }


def census_household_target_indices() -> dict[int, float]:
    """Return model-age HH-3 household totals indexed to 2007."""
    levels = census_household_model_age_levels()
    base = levels[2007]
    return {
        int((year - 2007) / 4): count / base
        for year, count in levels.items()
    }


def census_household_age_target_years() -> dict[int, int]:
    """Map transition dates to the HH-3 benchmark year at that date."""
    return {
        int((year - 2007) / 4): year
        for year in CENSUS_HH3_HOUSEHOLDS_THOUSANDS
    }


def reweight_distribution_to_observed_age_path(
    g_pre: np.ndarray,
    ages: np.ndarray,
    *,
    year: int,
    initial_mass: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Match HH-3 total mass and national ACS four-year age shares."""
    if year not in CENSUS_HH3_HOUSEHOLDS_THOUSANDS:
        raise ValueError(f"No Census HH-3 target is available for {year}")
    out = np.asarray(g_pre, dtype=float).copy()
    age_grid = np.asarray(ages, dtype=float).reshape(-1)
    if out.shape[3] != age_grid.size:
        raise ValueError("The HH-3 bridge does not match the model age axis")
    published_total = float(CENSUS_HH3_HOUSEHOLDS_THOUSANDS[year])
    model_age_coverage = float(
        ACS_NATIONAL_HEAD_AGE_RECEIPTS[year][
            "age_18_85_share_of_all_heads"
        ]
    )
    model_age_levels = census_household_model_age_levels()
    model_age_total = float(model_age_levels[year])
    base_model_age_total = float(model_age_levels[2007])
    target_total = float(initial_mass) * model_age_total / base_model_age_total
    target_shares = np.asarray(
        ACS_NATIONAL_HEAD_AGE_SHARES_18_85[year], dtype=float
    )
    if target_shares.shape != age_grid.shape:
        raise RuntimeError(f"The {year} ACS age shares do not match model ages")
    if not math.isclose(float(np.sum(target_shares)), 1.0, rel_tol=0.0, abs_tol=2e-14):
        raise RuntimeError(f"The {year} ACS age shares do not sum to one")
    target_age_masses = target_total * target_shares
    before_age_masses = np.sum(out, axis=(0, 1, 2, 4, 5, 6))
    group_rows: list[dict[str, Any]] = []
    for index, (age, before, target_mass, target_share) in enumerate(
        zip(age_grid, before_age_masses, target_age_masses, target_shares, strict=True)
    ):
        before = float(before)
        if before <= 0.0:
            raise RuntimeError(f"Model age {age:g} has no mass in {year}")
        factor = float(target_mass) / before
        out[:, :, :, index, :, :, :] *= factor
        after = float(np.sum(out[:, :, :, index, :, :, :]))
        group_rows.append(
            {
                "model_age": float(age),
                "closed_integer_age_bin": [int(age), int(age + 3)],
                "acs_share_conditional_18_85": float(target_share),
                "mass_before_bridge": before,
                "target_mass": float(target_mass),
                "mass_after_bridge": after,
                "net_residual": after - before,
                "multiplicative_factor": factor,
                "target_gap": after - float(target_mass),
            }
        )
    before_total = float(np.sum(before_age_masses))
    after_total = float(np.sum(out))
    if not math.isclose(after_total, target_total, rel_tol=0.0, abs_tol=2e-12):
        raise RuntimeError(
            f"Observed household age bridge missed the {year} total: "
            f"model={after_total:.12g}, target={target_total:.12g}"
        )
    return out, {
        "year": int(year),
        "hh3_published_total_thousands": published_total,
        "acs_age_18_85_share_of_all_heads": model_age_coverage,
        "hh3_total_for_model_ages_thousands": model_age_total,
        "acs_age_source": str(ACS_NATIONAL_HEAD_AGE_SOURCE),
        "acs_age_source_sha256": ACS_NATIONAL_HEAD_AGE_SOURCE_SHA256,
        "acs_source_receipt": ACS_NATIONAL_HEAD_AGE_RECEIPTS[year],
        "acs_age_normalization": "shares conditional on household-head age 18--85",
        "model_mass_before_bridge": before_total,
        "model_target_mass": target_total,
        "model_mass_after_bridge": after_total,
        "net_residual": after_total - before_total,
        "absolute_group_residual": float(
            np.sum([abs(float(row["net_residual"])) for row in group_rows])
        ),
        "maximum_absolute_target_gap": max(
            abs(float(row["target_gap"])) for row in group_rows
        ),
        "groups": group_rows,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def configure_sequential_model() -> Any:
    chain = closure_audit.load_chain(profile="e5f-floor")
    from intergen_eqscale_seq_optimized import solver as sequential_model

    # Reuse the already-tested calendar-time scaffolding, but route every
    # household and transition operation through the sequential model.
    calendar.model = sequential_model
    return chain, sequential_model


def apply_sequential_fertility(
    g_pre: np.ndarray,
    fert_probs: np.ndarray,
    P: SimpleNamespace,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Apply at most one sequential birth per household during the period."""
    model = calendar.model
    if not bool(getattr(P, "sequential_births", False)):
        raise RuntimeError("The E5F transition requires sequential_births=True")
    if not model.independent_child_maturation_active(P):
        raise RuntimeError("The E5F transition requires independent child maturation")

    out = np.asarray(g_pre, dtype=float).copy()
    fecundity = model.get_fecundity_by_age(P)
    continuation = getattr(P, "_fert2_probs", None)
    if continuation is None:
        raise RuntimeError("Sequential continuation-birth probabilities are unavailable")
    births = 0.0
    births_by_loc = np.zeros(int(P.I))
    for j in range(int(P.J)):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        pi_j = float(fecundity[j])
        for zz in range(g_pre.shape[4]):
            # Take every continuation at-risk pool before any birth lands, as
            # in the stationary KFE. This prevents two births in one period.
            at_risk = {
                (parity, child_state): g_pre[
                    :, :, :, j, zz, parity, child_state
                ].copy()
                for parity in range(1, int(P.n_parity) - 1)
                for child_state in range(parity + 1)
            }

            settled = model.readiness_settled_state(P)
            childless = g_pre[:, :, :, j, zz, 0, settled]
            try_mass = childless * fert_probs[:, :, :, j, zz, 1]
            realized = pi_j * try_mass
            out[:, :, :, j, zz, 0, settled] -= realized
            out[:, :, :, j, zz, 1, 1] += realized
            births += float(np.sum(realized))
            births_by_loc += np.sum(realized, axis=(0, 1))

            for parity in range(1, int(P.n_parity) - 1):
                for child_state in range(parity + 1):
                    pool = at_risk[(parity, child_state)]
                    attempt_probability = continuation[
                        :, :, :, j, zz, 1, parity - 1, child_state
                    ]
                    realized = pi_j * pool * attempt_probability
                    destination = model.birth_destination_child_state(P, child_state)
                    out[:, :, :, j, zz, parity, child_state] -= realized
                    out[:, :, :, j, zz, parity + 1, destination] += realized
                    births += float(np.sum(realized))
                    births_by_loc += np.sum(realized, axis=(0, 1))

    minimum = float(np.min(out))
    if minimum < -1e-13:
        raise RuntimeError(f"Sequential fertility produced negative mass: {minimum:.3e}")
    out[out < 0.0] = 0.0
    mass_error = float(np.sum(out) - np.sum(g_pre))
    if abs(mass_error) > 2e-10:
        raise RuntimeError(f"Sequential fertility failed mass conservation: {mass_error:.3e}")
    return out, births, births_by_loc


def children_at_home_units(distribution: np.ndarray, P: SimpleNamespace) -> float:
    total = 0.0
    for parity in range(1, int(P.n_parity)):
        for child_state in range(parity + 1):
            total += child_state * float(
                np.sum(distribution[:, :, :, :, parity, child_state])
            )
    return total


def advance_sequential_calendar_distribution(
    evaluation: calendar.PeriodEvaluation,
    entry_by_loc_next: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Advance incumbents one age cell and measure independently matured children."""
    model = calendar.model
    policy = evaluation.policy
    g_post = evaluation.g_post_fertility
    next_pre = np.zeros_like(g_post)
    mature_by_loc = np.zeros(int(P.I))
    deaths = 0.0
    _, _, Pi_z = model.income_transition_values(P)
    stochastic = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    if int(P.I) != 1:
        raise NotImplementedError("The sequential transition is currently one-market only")

    for j in range(int(P.J) - 1):
        survival = float(P.survival_probs[j]) if bool(P.use_age_survival) else 1.0
        cohort = g_post[:, :, :, j, :, :, :]
        deaths += (1.0 - survival) * float(np.sum(cohort))
        survivors = survival * cohort
        advanced = model.advance_cohort_one_period_markov_income(
            survivors,
            j,
            policy.loc_probs,
            policy.tenure_choice,
            policy.tenure_probs,
            policy.bp_pol,
            P,
            b_grid,
            shared,
            policy.maps.lmm_idx,
            policy.maps.lmm_wt,
            policy.maps.tmx_idx,
            policy.maps.tmx_wt,
            stochastic,
            P.Pi_child if stochastic else None,
            Pi_z,
        )
        next_pre[:, :, :, j + 1, :, :, :] = advanced
        matured_children = children_at_home_units(survivors, P) - children_at_home_units(
            advanced, P
        )
        if matured_children < -2e-10:
            raise RuntimeError(
                f"Children-at-home units rose during maturation at age index {j}: "
                f"{matured_children:.3e}"
            )
        mature_by_loc[0] += max(matured_children, 0.0)

    deaths += float(np.sum(g_post[:, :, :, int(P.J) - 1, :, :, :]))
    next_pre[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        np.asarray(entry_by_loc_next, dtype=float), P, b_grid
    )
    expected_mass = float(np.sum(g_post)) - deaths + float(np.sum(entry_by_loc_next))
    mass_residual = float(np.sum(next_pre)) - expected_mass
    return next_pre, mature_by_loc, deaths, mass_residual


def independent_child_distribution_rows(
    scenario: str,
    period: int,
    g_post: np.ndarray,
    P: SimpleNamespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    total = float(np.sum(g_post))
    age_rows: list[dict[str, Any]] = []
    for j in range(int(P.J)):
        mass = float(np.sum(g_post[:, :, :, j, :, :, :]))
        age_rows.append(
            {
                "scenario": scenario,
                "period": period,
                "years_from_start": period * float(P.period_years),
                "age": float(P.age_start + j * P.da),
                "adult_mass": mass,
                "adult_share": mass / max(total, 1e-15),
            }
        )
    child_rows: list[dict[str, Any]] = []
    for parity in range(int(P.n_parity)):
        for child_state in range(int(P.n_child_states)):
            mass = float(np.sum(g_post[:, :, :, :, :, parity, child_state]))
            valid_child_count = child_state if parity >= 1 and child_state <= parity else 0
            child_rows.append(
                {
                    "scenario": scenario,
                    "period": period,
                    "years_from_start": period * float(P.period_years),
                    "parity_state": parity,
                    "child_state": child_state,
                    "child_state_label": (
                        "no_children_at_home"
                        if valid_child_count == 0
                        else f"active_child_stage_{valid_child_count}"
                    ),
                    "household_mass": mass,
                    "household_share": mass / max(total, 1e-15),
                    "child_units": valid_child_count * mass,
                }
            )
    return age_rows, child_rows


def operator_gates(
    solution: SimpleNamespace,
    policy: calendar.PolicyBundle,
    g_pre: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> dict[str, Any]:
    P._fert2_probs = np.asarray(solution.fert2_probs, dtype=float).copy()
    evaluation = calendar.evaluate_period(
        policy.price,
        g_pre,
        P,
        b_grid,
        shared,
        calendar.SolveCounter(),
        supplied_policy=policy,
    )
    empty_next, mature_raw_by_loc, _, zero_entry_mass_residual = (
        advance_sequential_calendar_distribution(
            evaluation, np.zeros(int(P.I)), P, b_grid, shared
        )
    )
    empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        np.asarray(solution.entry_by_loc, dtype=float), P, b_grid
    )
    effective_mature = float(P.entrant_conversion_factor) * float(
        np.sum(mature_raw_by_loc)
    )
    return {
        "stationary_post_fertility_nesting_l1": float(
            np.sum(
                np.abs(
                    evaluation.g_post_fertility
                    - np.asarray(solution.g_beginning_distribution, dtype=float)
                )
            )
        ),
        "one_step_constant_path_nesting_l1": float(np.sum(np.abs(empty_next - g_pre))),
        "zero_entry_mass_accounting_residual": float(zero_entry_mass_residual),
        "entry_flow_E": float(solution.entry_rate),
        "effective_mature_flow_B": float(solution.entrants_mature_total),
        "reconstructed_effective_mature_flow_B": effective_mature,
        "mature_flow_abs_error": abs(effective_mature - float(solution.entrants_mature_total)),
        "birth_flow": float(evaluation.births),
        "stationary_birth_flow": float(solution.total_births_kfe),
        "birth_flow_abs_error": abs(float(evaluation.births) - float(solution.total_births_kfe)),
    }


def run_birth_vintage_scenario(
    *,
    label: str,
    initial_g_pre: np.ndarray,
    baseline_price: float,
    baseline_B: float,
    baseline_births: float,
    parameters: SimpleNamespace,
    b_grid: np.ndarray,
    periods: int,
    outside_flow: float,
    retention: float,
    conversion: float,
    delay_periods: int,
    initial_birth_pipeline_multiplier: float,
    population_target_indices: dict[int, float] | None,
    population_age_target_years: dict[int, int] | None,
    policy_case: str,
    policy_start_period: int,
    old_psi_child: float,
    new_psi_child: float,
    preference_transition_periods: int,
    supply_rule: calendar.HousingSupplyRule,
    market_tol: float,
    market_max_iter: int,
    outdir: Path,
    counter: calendar.SolveCounter,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Run the open transition with a dated birth-to-entry queue.

    The household child state continues to govern housing demand exactly as in
    the calibrated model. Population renewal instead uses a separate aggregate
    birth-vintage queue, so a newborn cannot become an entrant household in the
    same four-year period.
    """
    if delay_periods < 1:
        raise ValueError("The birth-to-entry delay must be at least one period")
    P = copy.deepcopy(parameters)
    shared = calendar.model.precompute_shared(P, b_grid)
    g_pre = np.asarray(initial_g_pre, dtype=float).copy()
    outside_shares = np.asarray(P.entry_shares, dtype=float).reshape(-1)
    outside_shares = outside_shares / np.sum(outside_shares)
    if int(P.I) != 1:
        raise NotImplementedError("The birth-vintage transition is currently one-market only")

    # At the old steady state this yield maps births exactly into the model's
    # effective mature-entrant flow. It combines child survival and the
    # children-to-households conversion, while the queue supplies the timing.
    maturation_survival_yield = baseline_B / max(conversion * baseline_births, 1e-15)
    if not 0.0 < maturation_survival_yield <= 1.0 + 1e-10:
        raise RuntimeError(
            "The stationary birth-to-entry yield is invalid: "
            f"{maturation_survival_yield:.12g}"
        )
    scheduled_entries = [
        float(baseline_B) * float(initial_birth_pipeline_multiplier)
    ] * int(delay_periods)
    target_indices = dict(population_target_indices or {})
    target_age_years = dict(population_age_target_years or {})
    initial_mass = float(np.sum(g_pre))
    if 0 in target_indices and not math.isclose(
        float(target_indices[0]), 1.0, rel_tol=0.0, abs_tol=1e-12
    ):
        raise ValueError("The date-0 population target must equal one")
    if policy_case not in {"none", "dependent-child-ltv95"}:
        raise ValueError(f"Unknown transition policy: {policy_case}")

    path_rows: list[dict[str, Any]] = []
    age_rows: list[dict[str, Any]] = []
    child_rows: list[dict[str, Any]] = []
    price_guess = float(baseline_price)
    initial_policy = None
    started = time.perf_counter()
    start_solves = counter.total
    max_market_residual = 0.0
    max_mass_residual = 0.0
    max_projection = 0.0
    min_mass = math.inf
    max_nonfinite = 0
    max_age_target_gap = 0.0
    census_age_bridge_audits: list[dict[str, Any]] = []
    grid_resolution_fallback_periods: list[int] = []
    previous_psi_child: float | None = None
    previous_policy_active: bool | None = None

    for period in range(int(periods)):
        dated_psi_child = preference_shifter_at_date(
            period,
            old_psi_child,
            new_psi_child,
            preference_transition_periods,
        )
        P.psi_child = dated_psi_child
        policy_active = policy_case != "none" and period >= int(policy_start_period)
        P.parent_dp_waiver = bool(policy_active)
        P.parent_dp_waiver_phi = 0.95 if policy_active else 1.0
        P.parent_dp_waiver_locations = np.array([], dtype=int)
        P.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
        P.parent_dp_waiver_birth_state_only = False
        shared = calendar.model.precompute_shared(P, b_grid)
        if previous_psi_child is not None and not math.isclose(
            dated_psi_child, previous_psi_child, rel_tol=0.0, abs_tol=1e-14
        ):
            # A policy bundle is indexed by both price and preferences.  It is
            # a valid warm start only while the preference shifter is unchanged.
            initial_policy = None
        if previous_policy_active is not None and policy_active != previous_policy_active:
            initial_policy = None
        previous_psi_child = dated_psi_child
        previous_policy_active = policy_active
        target_tolerance = min(float(market_tol), 5e-5)
        try:
            evaluation = calendar.clear_scalar_housing_market(
                g_pre,
                price_guess,
                P,
                b_grid,
                shared,
                counter,
                target_tolerance,
                min(int(market_max_iter), 18),
                supply_rule,
                initial_policy=initial_policy,
            )
        except RuntimeError as error:
            if (
                float(market_tol) <= target_tolerance
                or "Housing market did not clear" not in str(error)
            ):
                raise
            # The deterministic tenure/wealth grid can jump across the exact
            # market root. Retry only the failed date at the declared grid-
            # resolution acceptance tolerance; all smooth dates retain the
            # tighter 5e-5 target.
            evaluation = calendar.clear_scalar_housing_market(
                g_pre,
                price_guess,
                P,
                b_grid,
                shared,
                counter,
                float(market_tol),
                int(market_max_iter),
                supply_rule,
                initial_policy=initial_policy,
            )
            grid_resolution_fallback_periods.append(period)
        initial_policy = evaluation.policy
        price_guess = float(evaluation.policy.price[0])
        entry_flow = float(np.sum(evaluation.g_pre[:, :, :, 0, :, :, :]))

        empty_next, model_mature_raw_by_loc, deaths, _ = (
            advance_sequential_calendar_distribution(
                evaluation, np.zeros(int(P.I)), P, b_grid, shared
            )
        )
        scheduled_B = float(scheduled_entries.pop(0))
        new_potential_B = (
            conversion * maturation_survival_yield * float(evaluation.births)
        )
        scheduled_entries.append(new_potential_B)
        retained_B = retention * scheduled_B
        next_target_index = target_indices.get(period + 1)
        dated_outside_flow = float(outside_flow)
        entrants_next_by_loc = dated_outside_flow * outside_shares
        entrants_next_by_loc[0] += retained_B
        empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
            entrants_next_by_loc, P, b_grid
        )
        provisional_next = empty_next
        bridge_audit = None
        target_next_year = target_age_years.get(period + 1)
        if target_next_year is not None:
            ages = P.age_start + np.arange(int(P.J)) * P.da
            next_pre, bridge_audit = reweight_distribution_to_observed_age_path(
                provisional_next,
                ages,
                year=int(target_next_year),
                initial_mass=initial_mass,
            )
            bridge_audit["transition_from_period"] = int(period)
            bridge_audit["transition_to_period"] = int(period + 1)
            census_age_bridge_audits.append(bridge_audit)
            max_age_target_gap = max(
                max_age_target_gap,
                float(bridge_audit["maximum_absolute_target_gap"]),
            )
        else:
            next_pre = provisional_next

        health = calendar.distribution_health(
            {
                "pre_fertility": evaluation.g_pre,
                "post_fertility": evaluation.g_post_fertility,
                "current_choice": evaluation.g_current,
                "next_pre_fertility": next_pre,
            }
        )
        period_min = health["min_distribution_mass"]
        nonfinite = int(health["nonfinite_distribution_count"])
        if period_min is None or nonfinite > 0 or float(period_min) < -1e-14:
            raise RuntimeError(
                f"{label} t={period}: distribution-health gate failed: {health}"
            )
        min_mass = min(min_mass, float(period_min))
        max_nonfinite = max(max_nonfinite, nonfinite)

        expected_next_mass = (
            float(np.sum(evaluation.g_post_fertility))
            - deaths
            + float(np.sum(entrants_next_by_loc))
            + (
                float(bridge_audit["net_residual"])
                if bridge_audit is not None
                else 0.0
            )
        )
        mass_residual = float(np.sum(next_pre)) - expected_next_mass
        max_mass_residual = max(max_mass_residual, abs(mass_residual))
        max_market_residual = max(
            max_market_residual, float(evaluation.relative_market_residual)
        )
        max_projection = max(
            max_projection, float(evaluation.feasibility_projection_mass)
        )

        adult_population = float(np.sum(evaluation.g_post_fertility))
        population_index = adult_population / max(initial_mass, 1e-15)
        population_target_index = target_indices.get(period)
        population_target_gap = (
            population_index - float(population_target_index)
            if population_target_index is not None
            else math.nan
        )
        current_mass = float(np.sum(evaluation.g_current))
        owner_rate = float(
            np.sum(evaluation.g_current[:, 1:, :, :, :, :, :])
        ) / max(current_mass, 1e-15)
        dependent_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, :, 1:]))
        dependent_owner_rate = float(
            np.sum(evaluation.g_current[:, 1:, :, :, :, :, 1:])
        ) / max(dependent_mass, 1e-15)
        parent_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, 1:, :]))
        parent_owner_rate = float(
            np.sum(evaluation.g_current[:, 1:, :, :, :, 1:, :])
        ) / max(parent_mass, 1e-15)
        childless_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, 0, :]))
        childless_owner_rate = float(
            np.sum(evaluation.g_current[:, 1:, :, :, :, 0, :])
        ) / max(childless_mass, 1e-15)
        ages = P.age_start + np.arange(int(P.J)) * P.da
        age_mass = np.sum(
            evaluation.g_post_fertility, axis=(0, 1, 2, 4, 5, 6)
        )
        mean_age = float(np.sum(ages * age_mass) / max(np.sum(age_mass), 1e-15))
        demand = float(np.sum(evaluation.demand_by_loc))
        supply = float(np.sum(evaluation.supply_by_loc))
        model_mature_effective = conversion * float(np.sum(model_mature_raw_by_loc))
        row = {
            "scenario": label,
            "period": period,
            "years_from_start": period * float(P.period_years),
            "psi_child": dated_psi_child,
            "asset_price": price_guess,
            "asset_price_index": price_guess / baseline_price,
            "housing_user_cost": float(P.user_cost_rate * price_guess),
            "adult_population": adult_population,
            "population_index": population_index,
            "population_target_index": population_target_index,
            "population_target_gap": population_target_gap,
            "next_population_target_index": next_target_index,
            "next_population_target_year": target_next_year,
            "mean_adult_age": mean_age,
            "entry_flow_E": entry_flow,
            "birth_children": float(evaluation.births),
            "births_per_adult": float(evaluation.births) / max(adult_population, 1e-15),
            "births_over_entry": float(evaluation.births) / max(entry_flow, 1e-15),
            "mature_children_raw": scheduled_B / max(conversion, 1e-15),
            "effective_mature_entrant_flow_B": scheduled_B,
            "mature_to_current_entry_flow_ratio_diagnostic": scheduled_B / max(entry_flow, 1e-15),
            "model_state_same_period_mature_flow_B": model_mature_effective,
            "birth_queue_new_potential_flow_B": new_potential_B,
            "birth_queue_scheduled_flows": list(scheduled_entries),
            "closure": "open_birth_vintage",
            "housing_supply_mode": supply_rule.mode,
            "outside_entry_flow_M": dated_outside_flow,
            "outside_entry_flow_share_of_initial_mass": dated_outside_flow
            / max(initial_mass, 1e-15),
            "retained_mature_entrants": retained_B,
            "entrant_flow_next": float(np.sum(entrants_next_by_loc)),
            "youngest_age_mass_next": float(
                np.sum(next_pre[:, :, :, 0, :, :, :])
            ),
            "census_age_bridge_net_residual": (
                float(bridge_audit["net_residual"])
                if bridge_audit is not None
                else 0.0
            ),
            "census_age_bridge_absolute_group_residual": (
                float(bridge_audit["absolute_group_residual"])
                if bridge_audit is not None
                else 0.0
            ),
            "census_age_bridge_maximum_target_gap": (
                float(bridge_audit["maximum_absolute_target_gap"])
                if bridge_audit is not None
                else 0.0
            ),
            "adult_deaths": deaths,
            "housing_demand": demand,
            "housing_demand_per_adult": demand / max(adult_population, 1e-15),
            "housing_supply": supply,
            "relative_market_residual": float(evaluation.relative_market_residual),
            "owner_rate": owner_rate,
            "dependent_child_owner_rate": dependent_owner_rate,
            "parent_owner_rate": parent_owner_rate,
            "childless_owner_rate": childless_owner_rate,
            "policy_case": policy_case,
            "policy_active": policy_active,
            "mass_accounting_residual": mass_residual,
            "feasibility_frontier_projection_mass": float(
                evaluation.feasibility_projection_mass
            ),
            "min_distribution_mass": period_min,
            "nonfinite_distribution_count": nonfinite,
        }
        path_rows.append(row)
        period_ages, period_children = independent_child_distribution_rows(
            label, period, evaluation.g_post_fertility, P
        )
        age_rows.extend(period_ages)
        child_rows.extend(period_children)
        calendar.write_json_atomic(
            outdir / "latest_completed_period.json",
            {
                "status": "running",
                "scenario": label,
                "completed_period": period,
                "periods_requested": periods,
                "model_solve_count": counter.total,
                "latest": row,
            },
        )
        print(
            f"{label} t={period:02d} p={price_guess:.6f} "
            f"pop={adult_population:.6f} births={evaluation.births:.6f} "
            f"queue_B={scheduled_B:.6f} market={evaluation.relative_market_residual:.2e}",
            flush=True,
        )
        g_pre = next_pre

    scenario_summary = {
        "scenario": label,
        "closure": "open_birth_vintage",
        "housing_supply_mode": supply_rule.mode,
        "post_target_outside_entry_flow_M": outside_flow,
        "population_target_indices": target_indices,
        "population_age_target_years": target_age_years,
        "maximum_absolute_population_target_gap": max(
            (
                abs(float(row["population_target_gap"]))
                for row in path_rows
                if math.isfinite(float(row["population_target_gap"]))
            ),
            default=0.0,
        ),
        "census_age_bridge_audits": census_age_bridge_audits,
        "maximum_absolute_census_age_target_gap": max_age_target_gap,
        "policy_case": policy_case,
        "policy_start_period": policy_start_period,
        "renewal_retention": retention,
        "entrant_conversion_factor": conversion,
        "maturation_survival_yield": maturation_survival_yield,
        "birth_to_entry_delay_periods": delay_periods,
        "birth_to_entry_delay_years": delay_periods * float(P.period_years) + float(P.period_years),
        "initial_birth_pipeline_multiplier": initial_birth_pipeline_multiplier,
        "old_psi_child": old_psi_child,
        "new_psi_child": new_psi_child,
        "preference_transition_periods": preference_transition_periods,
        "preference_transition_years": preference_transition_periods * float(P.period_years),
        "periods": periods,
        "model_solve_count": counter.total - start_solves,
        "elapsed_seconds": time.perf_counter() - started,
        "max_market_residual": max_market_residual,
        "market_tolerance": market_tol,
        "smooth_date_target_tolerance": min(float(market_tol), 5e-5),
        "grid_resolution_fallback_periods": grid_resolution_fallback_periods,
        "market_resolution_note": (
            "Deterministic tenure choices on the finite wealth grid make the "
            "off-stationary demand schedule locally discontinuous; the reported "
            "price is the first bracketed grid equilibrium within this tolerance."
        ),
        "max_abs_mass_accounting_residual": max_mass_residual,
        "max_feasibility_frontier_projection_mass": max_projection,
        "min_distribution_mass": min_mass if math.isfinite(min_mass) else None,
        "max_nonfinite_distribution_count": max_nonfinite,
    }
    if max_market_residual > market_tol:
        raise RuntimeError(
            f"{label}: market gate failed, max residual={max_market_residual:.3e}"
        )
    if max_mass_residual > 2e-8:
        raise RuntimeError(
            f"{label}: mass gate failed, max residual={max_mass_residual:.3e}"
        )
    if max_projection > 1e-6:
        raise RuntimeError(
            f"{label}: feasibility projection gate failed, mass={max_projection:.3e}"
        )
    if scenario_summary["maximum_absolute_population_target_gap"] > 2e-10:
        raise RuntimeError(
            f"{label}: Census household-count target gap is "
            f"{scenario_summary['maximum_absolute_population_target_gap']:.3e}"
        )
    if max_age_target_gap > 2e-10:
        raise RuntimeError(
            f"{label}: Census householder-age target gap is "
            f"{max_age_target_gap:.3e}"
        )
    return path_rows, age_rows, child_rows, scenario_summary


def solve_stationary_open_endpoint(
    *,
    chain: Any,
    base_overrides: dict[str, Any],
    old_price: float,
    new_psi_child: float,
    outside_flow: float,
    retention: float,
    conversion: float,
    maturation_survival_yield: float,
    supply_rule: calendar.HousingSupplyRule,
    policy_case: str,
    outdir: Path,
) -> dict[str, Any]:
    """Solve the stationary endpoint consistent with the birth-vintage closure."""
    cache: dict[float, tuple[float, dict[str, Any]]] = {}
    endpoint_overrides = dict(base_overrides)
    if policy_case == "dependent-child-ltv95":
        endpoint_overrides.update(
            parent_dp_waiver=True,
            parent_dp_waiver_phi=0.95,
            parent_dp_waiver_locations=np.array([], dtype=int),
            parent_dp_waiver_owner_rungs=np.array([], dtype=int),
            parent_dp_waiver_birth_state_only=False,
        )
    elif policy_case != "none":
        raise ValueError(f"Unknown stationary policy: {policy_case}")

    def evaluate(asset_price: float) -> tuple[float, dict[str, Any]]:
        key = round(float(asset_price), 12)
        if key not in cache:
            solution, parameters, price, elapsed = closure_audit.solve_pe(
                chain,
                endpoint_overrides,
                asset_price=float(asset_price),
                psi_child=float(new_psi_child),
            )
            readout = closure_audit.readout(
                chain,
                solution,
                parameters,
                price,
                label="postshock_open_endpoint_search",
                price_ratio=float(asset_price / old_price),
                psi_child=float(new_psi_child),
                elapsed=elapsed,
            )
            queue_B = (
                conversion
                * maturation_survival_yield
                * float(solution.total_births_kfe)
            )
            denominator = float(solution.entry_rate) - retention * queue_B
            if denominator <= 0.0:
                population_scale = math.inf
                market_gap = math.inf
            else:
                population_scale = outside_flow / denominator
                stationary_supply = float(
                    supply_rule.quantity(np.array([float(asset_price)]))[0]
                )
                market_gap = (
                    population_scale * float(readout["housing_demand_per_adult"])
                    - stationary_supply
                )
            row = dict(readout)
            row.update(
                policy_case=policy_case,
                queue_mature_entrant_flow_B=queue_B,
                renewal_denominator=denominator,
                stationary_population_scale=population_scale,
                housing_supply_mode=supply_rule.mode,
                stationary_housing_supply=(
                    float(supply_rule.quantity(np.array([float(asset_price)]))[0])
                    if math.isfinite(population_scale)
                    else math.nan
                ),
                fixed_stock_market_gap=market_gap,
                fixed_stock_relative_market_gap=(
                    market_gap / max(
                        float(supply_rule.quantity(np.array([float(asset_price)]))[0]),
                        1e-15,
                    )
                    if math.isfinite(market_gap)
                    else math.inf
                ),
            )
            cache[key] = float(market_gap), row
            calendar.write_csv(outdir / "stationary_endpoint_search.csv", [item[1] for item in cache.values()])
        return cache[key]

    # The fixed-stock case can approach a vacancy corner.  Start essentially
    # at zero so a failed bracket is evidence against a positive-price root,
    # not merely evidence against a root above five percent of the old price.
    lower = 1e-4
    upper = 1.25 * float(old_price)
    f_lower, _ = evaluate(lower)
    f_upper, _ = evaluate(upper)
    audited_prices = [lower, upper]
    if math.isfinite(f_lower) and math.isfinite(f_upper) and f_lower * f_upper > 0.0:
        audited_prices = list(np.geomspace(lower, upper, 25))
        audited = [(price, *evaluate(float(price))) for price in audited_prices]
        bracket = None
        for left_row, right_row in zip(audited[:-1], audited[1:]):
            left_price, left_gap, _ = left_row
            right_price, right_gap, _ = right_row
            if (
                math.isfinite(left_gap)
                and math.isfinite(right_gap)
                and left_gap * right_gap <= 0.0
            ):
                bracket = (float(left_price), float(left_gap), float(right_price), float(right_gap))
                break
        if bracket is not None:
            lower, f_lower, upper, f_upper = bracket
        else:
            finite_rows = [row for row in audited if math.isfinite(row[1])]
            max_row = max(finite_rows, key=lambda row: row[1]) if finite_rows else None
            result = {
                "status": "no_positive_price_root_on_audited_grid",
                "housing_supply_mode": supply_rule.mode,
                "lower_asset_price": lower,
                "lower_market_gap": f_lower,
                "upper_asset_price": upper,
                "upper_market_gap": f_upper,
                "audit_grid_points": len(audited_prices),
                "maximum_market_gap_on_grid": max_row[1] if max_row is not None else math.nan,
                "price_at_maximum_market_gap": max_row[0] if max_row is not None else math.nan,
                "interpretation": (
                    "At the post-shock renewal scale, households demand less than "
                    "the inherited fixed stock throughout the audited positive-price "
                    "grid; a terminal equilibrium needs vacancies, scrappage, or "
                    "stock adjustment."
                    if supply_rule.mode == "fixed-stock"
                    else "The post-shock stationary equilibrium has no sign-changing "
                    "market root on the audited price grid."
                ),
            }
            calendar.write_csv(outdir / "stationary_open_endpoint.csv", [result])
            return result
    if not math.isfinite(f_lower) or not math.isfinite(f_upper):
        result = {
            "status": "nonfinite_stationary_endpoint_bracket",
            "housing_supply_mode": supply_rule.mode,
            "lower_asset_price": lower,
            "lower_market_gap": f_lower,
            "upper_asset_price": upper,
            "upper_market_gap": f_upper,
            "interpretation": "The stationary endpoint bracket contains a nonfinite market gap.",
        }
        calendar.write_csv(outdir / "stationary_open_endpoint.csv", [result])
        return result
    midpoint = 0.5 * (lower + upper)
    for _ in range(28):
        midpoint = 0.5 * (lower + upper)
        f_midpoint, row = evaluate(midpoint)
        if abs(float(row["fixed_stock_relative_market_gap"])) <= 2.5e-5:
            break
        if f_lower * f_midpoint <= 0.0:
            upper = midpoint
            f_upper = f_midpoint
        else:
            lower = midpoint
            f_lower = f_midpoint
    _, root = evaluate(midpoint)
    root["status"] = "complete"
    calendar.write_csv(outdir / "stationary_open_endpoint.csv", [root])
    return root


def make_paper_transition_plot(
    paths: list[dict[str, Any]],
    *,
    old_births: float,
    old_B: float,
    stationary_endpoint: dict[str, Any] | None,
    delay_years: float,
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([row["years_from_start"] for row in paths], dtype=float)
    population = np.asarray([row["population_index"] for row in paths], dtype=float)
    prices = np.asarray([row["asset_price_index"] for row in paths], dtype=float)
    births = np.asarray([row["birth_children"] for row in paths], dtype=float) / old_births
    mature = np.asarray(
        [row["effective_mature_entrant_flow_B"] for row in paths], dtype=float
    ) / old_B
    period_years = float(years[1] - years[0]) if len(years) > 1 else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), constrained_layout=True)
    axes[0].plot(years, population, lw=2.2, color="#21476b", label="Adult households")
    axes[0].plot(years, prices, lw=2.2, color="#b14f32", label="House price")
    if stationary_endpoint is not None and stationary_endpoint.get("status") == "complete":
        axes[0].axhline(
            float(stationary_endpoint["stationary_population_scale"]),
            color="#21476b",
            lw=1.1,
            ls=":",
            label="New population steady state",
        )
        axes[0].axhline(
            float(stationary_endpoint["price_ratio"]),
            color="#b14f32",
            lw=1.1,
            ls=":",
            label="New price steady state",
        )
    axes[0].axhline(1.0, color="black", lw=0.8, alpha=0.45)
    axes[0].set(
        title="Movement toward the new steady state",
        xlabel="Years after the preference change",
        ylabel="Index (old steady state = 1)",
    )
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].plot(years, births, lw=2.2, color="#4f7f3f", label="Births")
    axes[1].plot(
        years + period_years,
        mature,
        lw=2.2,
        color="#6a4c93",
        label="Locally born entrants",
    )
    axes[1].axvline(delay_years, color="black", lw=1.0, ls="--", alpha=0.65)
    axes[1].axhline(1.0, color="black", lw=0.8, alpha=0.45)
    axes[1].set(
        title="Births reach entry with a delay",
        xlabel="Years after the preference change",
        ylabel="Flow index (old steady state = 1)",
    )
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2)
    fig.savefig(outdir / "steady_state_transition.png", dpi=220)
    fig.savefig(outdir / "steady_state_transition.pdf")
    plt.close(fig)


def make_population_target_accounting_plot(
    paths: list[dict[str, Any]],
    *,
    post_target_outside_flow: float,
    outdir: Path,
) -> None:
    """Plot the Census household target and age-specific bridge residual."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([row["years_from_start"] for row in paths], dtype=float)
    model = np.asarray([row["population_index"] for row in paths], dtype=float)
    target = np.asarray(
        [
            float(row["population_target_index"])
            if row["population_target_index"] is not None
            else math.nan
            for row in paths
        ],
        dtype=float,
    )
    net_residual = np.asarray(
        [row["census_age_bridge_net_residual"] for row in paths],
        dtype=float,
    ) / max(float(paths[0]["adult_population"]), 1e-15)
    absolute_residual = np.asarray(
        [row["census_age_bridge_absolute_group_residual"] for row in paths],
        dtype=float,
    ) / max(float(paths[0]["adult_population"]), 1e-15)
    initial_mass = float(paths[0]["adult_population"])

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.8), constrained_layout=True)
    axes[0].plot(years, model, color="#21476b", lw=2.2, label="Model households")
    axes[0].plot(
        years,
        target,
        color="#222222",
        marker="o",
        lw=1.8,
        label="Census HH-3 target",
    )
    axes[0].set(
        title="Household-count bridge",
        xlabel="Years after 2007",
        ylabel="Index (2007 = 1)",
    )
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].plot(
        years,
        net_residual,
        color="#6a4c93",
        marker="o",
        lw=2.0,
        label="Net bridge residual",
    )
    axes[1].plot(
        years,
        absolute_residual,
        color="#b14f32",
        marker="s",
        lw=1.6,
        label="Absolute age-group reallocation",
    )
    axes[1].axhline(
        float(post_target_outside_flow) / max(initial_mass, 1e-15),
        color="#555555",
        lw=1.1,
        ls="--",
        label="Maintained outside flow",
    )
    axes[1].set(
        title="Reduced-form age bridge",
        xlabel="Transition beginning at model date",
        ylabel="Flow / 2007 household mass",
    )
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2)
    fig.savefig(outdir / "population_target_accounting.png", dpi=220)
    fig.savefig(outdir / "population_target_accounting.pdf")
    plt.close(fig)


FRED_COMPARISON_SERIES = {
    "house_price": {
        "id": "CSUSHPINSA",
        "title": "S&P Cotality Case-Shiller U.S. National Home Price Index",
    },
    "rent": {
        "id": "CUUR0000SEHA",
        "title": "Consumer Price Index for All Urban Consumers: Rent of Primary Residence",
    },
    "consumer_price": {
        "id": "CPIAUCSL",
        "title": "Consumer Price Index for All Urban Consumers: All Items",
    },
    "population": {
        "id": "POPTHM",
        "title": "Population, Total for United States",
    },
    "fertility": {
        "id": "SPDYNTFRTINUSA",
        "title": "Fertility Rate, Total for the United States",
    },
}


def download_fred_annual_comparison(outdir: Path) -> tuple[dict[str, dict[int, float]], dict[str, Any]]:
    """Download and archive annual averages of the comparison series."""
    raw_dir = outdir / "observed_fred_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    annual: dict[str, dict[int, float]] = {}
    metadata: dict[str, Any] = {
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "aggregation": "arithmetic mean of nonmissing observations within calendar year",
        "series": {},
    }
    for name, specification in FRED_COMPARISON_SERIES.items():
        series_id = str(specification["id"])
        url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={series_id}"
        completed = subprocess.run(
            [
                "curl",
                "--location",
                "--fail",
                "--silent",
                "--show-error",
                "--max-time",
                "45",
                url,
            ],
            check=True,
            capture_output=True,
            timeout=50,
        )
        raw = completed.stdout
        if not raw:
            raise RuntimeError(f"FRED returned an empty file for {series_id}")
        raw_path = raw_dir / f"{series_id}.csv"
        raw_path.write_bytes(raw)
        values_by_year: dict[int, list[float]] = {}
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
        for row in reader:
            date = str(row.get("observation_date", ""))
            value_text = str(row.get(series_id, "")).strip()
            if len(date) < 4 or value_text in {"", "."}:
                continue
            try:
                year = int(date[:4])
                value = float(value_text)
            except ValueError:
                continue
            if math.isfinite(value):
                values_by_year.setdefault(year, []).append(value)
        annual[name] = {
            year: float(sum(values) / len(values))
            for year, values in values_by_year.items()
            if values
        }
        metadata["series"][name] = {
            **specification,
            "fred_page": f"https://fred.stlouisfed.org/series/{series_id}",
            "download_url": url,
            "raw_file": str(raw_path),
            "raw_sha256": hashlib.sha256(raw).hexdigest(),
            "first_year": min(annual[name]),
            "last_year": max(annual[name]),
        }
    write_json(outdir / "observed_data_metadata.json", metadata)
    return annual, metadata


def make_historical_comparison(
    paths: list[dict[str, Any]],
    *,
    start_year: int,
    old_user_cost: float,
    household_counts_by_year: dict[int, float] | None,
    outdir: Path,
) -> dict[str, Any]:
    """Compare the deliberately sparse transition with broad U.S. aggregates."""
    annual, source_metadata = download_fred_annual_comparison(outdir)
    household_counts = dict(household_counts_by_year or {})
    required = (
        ("house_price", "rent", "consumer_price")
        if household_counts
        else ("house_price", "rent", "consumer_price", "population")
    )
    missing_base = [name for name in required if start_year not in annual[name]]
    if missing_base:
        raise RuntimeError(
            f"FRED comparison has no {start_year} normalization for {missing_base}"
        )
    base_house_price = annual["house_price"][start_year]
    base_rent = annual["rent"][start_year]
    base_cpi = annual["consumer_price"][start_year]
    if household_counts and start_year not in household_counts:
        raise RuntimeError(f"Census household path has no {start_year} normalization")
    base_population = (
        household_counts[start_year]
        if household_counts
        else annual["population"][start_year]
    )
    base_real_house_price = base_house_price / base_cpi
    base_real_rent = base_rent / base_cpi
    base_price_to_rent = base_house_price / base_rent

    rows: list[dict[str, Any]] = []
    for model_row in paths:
        calendar_year_float = start_year + float(model_row["years_from_start"])
        calendar_year = int(round(calendar_year_float))
        if abs(calendar_year_float - calendar_year) > 1e-8:
            continue
        if household_counts and calendar_year not in household_counts:
            continue
        if any(calendar_year not in annual[name] for name in required):
            continue
        house_price = annual["house_price"][calendar_year]
        rent = annual["rent"][calendar_year]
        cpi = annual["consumer_price"][calendar_year]
        population = (
            household_counts[calendar_year]
            if household_counts
            else annual["population"][calendar_year]
        )
        model_price = float(model_row["asset_price_index"])
        model_rent = float(model_row["housing_user_cost"]) / old_user_cost
        rows.append(
            {
                "calendar_year": calendar_year,
                "model_period": int(model_row["period"]),
                "model_psi_child": float(model_row["psi_child"]),
                "model_population_index": float(model_row["population_index"]),
                "observed_population_index": population / base_population,
                "model_asset_price_index": model_price,
                "model_rent_index": model_rent,
                "model_price_to_rent_index": model_price / max(model_rent, 1e-15),
                "observed_nominal_house_price_index": house_price / base_house_price,
                "observed_nominal_rent_index": rent / base_rent,
                "observed_real_house_price_index": (house_price / cpi) / base_real_house_price,
                "observed_real_rent_index": (rent / cpi) / base_real_rent,
                "observed_price_to_rent_index": (house_price / rent) / base_price_to_rent,
                "observed_tfr": annual["fertility"].get(calendar_year),
            }
        )
    if not rows:
        raise RuntimeError("No model dates overlap the downloaded historical series")
    calendar.write_csv(outdir / "historical_comparison.csv", rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([row["calendar_year"] for row in rows], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9), constrained_layout=True)
    axes[0].plot(
        years,
        [row["observed_population_index"] for row in rows],
        color="#222222",
        marker="o",
        lw=2.0,
        label="U.S. households" if household_counts else "U.S. population",
    )
    axes[0].plot(
        years,
        [row["model_population_index"] for row in rows],
        color="#21476b",
        marker="o",
        lw=2.0,
        label="Model",
    )
    axes[0].set(
        title="Adult households" if household_counts else "Population",
        ylabel=f"Index ({start_year} = 1)",
    )
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(
        years,
        [row["observed_real_house_price_index"] for row in rows],
        color="#b14f32",
        marker="o",
        lw=2.0,
        label="Real house price",
    )
    axes[1].plot(
        years,
        [row["observed_real_rent_index"] for row in rows],
        color="#7f7f7f",
        marker="o",
        lw=2.0,
        label="Real rent",
    )
    axes[1].plot(
        years,
        [row["model_asset_price_index"] for row in rows],
        color="#21476b",
        marker="o",
        lw=2.0,
        ls="--",
        label="Model price and rent",
    )
    axes[1].set(title="Housing costs", ylabel=f"Index ({start_year} = 1)")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(
        years,
        [row["observed_price_to_rent_index"] for row in rows],
        color="#222222",
        marker="o",
        lw=2.0,
        label="U.S. data",
    )
    axes[2].plot(
        years,
        [row["model_price_to_rent_index"] for row in rows],
        color="#6a4c93",
        marker="o",
        lw=2.0,
        label="Model user-cost restriction",
    )
    axes[2].set(title="House-price-to-rent ratio", ylabel=f"Index ({start_year} = 1)")
    axes[2].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xlabel("Year")
        axis.axhline(1.0, color="black", lw=0.7, alpha=0.35)
        axis.grid(alpha=0.18)
    fig.savefig(outdir / "historical_comparison.png", dpi=220)
    fig.savefig(outdir / "historical_comparison.pdf")
    plt.close(fig)

    first = rows[0]
    last = rows[-1]
    summary = {
        "start_year": start_year,
        "last_overlapping_year": int(last["calendar_year"]),
        "number_of_model_dates": len(rows),
        "normalization": f"annual averages indexed to {start_year}=1",
        "population_units": (
            "Census HH-3 households with heads ages 18--85 versus model adult-household mass"
            if household_counts
            else "FRED total persons versus model adult-household mass"
        ),
        "house_prices_and_rents": "deflated by CPIAUCSL; price-to-rent uses nominal ratios",
        "first_row": first,
        "last_row": last,
        "sources": source_metadata["series"],
    }
    if household_counts:
        summary["sources"]["households"] = {
            "title": (
                "Table HH-3 total households, restricted to head ages 18--85 "
                "using national ACS age coverage"
            ),
            "url": CENSUS_HH3_URL,
            "published_hh3_totals_thousands": CENSUS_HH3_HOUSEHOLDS_THOUSANDS,
            "model_age_totals_thousands": household_counts,
            "revision_note": "The revised 2011 row is used.",
        }
    write_json(outdir / "historical_comparison_summary.json", summary)
    return summary


def main() -> None:
    args = parse_args()
    if args.periods < 1:
        raise ValueError("--periods must be positive")
    if args.preference_transition_periods < 0:
        raise ValueError("--preference-transition-periods cannot be negative")
    if args.initial_birth_pipeline_multiplier <= 0.0:
        raise ValueError("--initial-birth-pipeline-multiplier must be positive")
    if (
        args.housing_supply_elasticity is not None
        and args.housing_supply_elasticity < 0.0
    ):
        raise ValueError("--housing-supply-elasticity cannot be negative")
    if (
        args.housing_supply_elasticity is not None
        and args.housing_supply_mode != "static-elastic"
    ):
        raise ValueError(
            "--housing-supply-elasticity requires --housing-supply-mode static-elastic"
        )
    if args.preference_transition_periods > 0 and args.renewal_clock != "birth-vintage":
        raise ValueError(
            "A gradual preference path is implemented only for the paper-facing "
            "birth-vintage renewal clock."
        )
    if args.policy_start_period < 0:
        raise ValueError("--policy-start-period cannot be negative")
    if (
        args.initial_age_profile != "stationary"
        or args.household_count_path != "none"
    ) and args.historical_start_year != 2007:
        raise ValueError(
            "The ACS age and Census household-count options require "
            "--historical-start-year 2007"
        )
    if args.household_count_path != "none" and args.renewal_clock != "birth-vintage":
        raise ValueError("The Census household-count bridge requires the birth-vintage clock")
    if (
        args.household_count_path == "census-hh3-2007-2023"
        and args.initial_age_profile != "acs-national-heads-2007"
    ):
        raise ValueError(
            "The observed age path requires "
            "--initial-age-profile acs-national-heads-2007"
        )
    if args.policy_case != "none" and args.renewal_clock != "birth-vintage":
        raise ValueError("The dated tenure policy requires the birth-vintage clock")
    if not 0.0 < args.outside_origin_entry_share < 1.0:
        raise ValueError("--outside-origin-entry-share must lie strictly between zero and one")
    if args.smoke:
        smoke_horizon = 5 if args.household_count_path != "none" else 4
        args.periods = min(int(args.periods), smoke_horizon)

    started = time.perf_counter()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve()
    chain, model = configure_sequential_model()
    winner, metadata = closure_audit.load_winner(source)
    theta = closure_audit.theta_from_winner(winner)
    base = closure_audit.make_overrides(
        chain, theta, nb=120, profile="e5f-floor"
    )
    old_overrides = dict(base)
    old_overrides["psi_child"] = float(args.old_psi_child)

    print("OLD_STEADY_STATE_START", flush=True)
    old_solution, old_parameters, old_price, old_seconds = closure_audit.solve_ge(
        chain, old_overrides
    )
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError(
            f"Production dimension gate failed: J={old_parameters.J}, Nb={old_parameters.Nb}"
        )
    if not model.independent_child_maturation_active(old_parameters):
        raise RuntimeError("Source routing did not activate independent child maturation")
    if not bool(getattr(old_parameters, "child_room_floor", False)):
        raise RuntimeError("Source routing did not activate the child-room floor")

    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    old_shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price, old_parameters, b_grid, old_shared
    )

    # Install the sequential operators before calling any generic calendar helper.
    calendar.apply_fertility = apply_sequential_fertility
    calendar.advance_calendar_distribution = advance_sequential_calendar_distribution
    calendar.distribution_rows = independent_child_distribution_rows

    initial_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution,
        old_policy,
        old_parameters,
        b_grid,
        old_shared,
    )
    gates = operator_gates(
        old_solution,
        old_policy,
        initial_g_pre,
        old_parameters,
        b_grid,
        old_shared,
    )
    gates.update(reconstruction)
    gate_tolerance = 5e-9
    gates["passed"] = bool(
        gates["stationary_post_fertility_nesting_l1"] <= gate_tolerance
        and gates["one_step_constant_path_nesting_l1"] <= gate_tolerance
        and gates["mature_flow_abs_error"] <= gate_tolerance
        and gates["birth_flow_abs_error"] <= gate_tolerance
    )
    if not gates["passed"]:
        write_json(
            outdir / "summary.json",
            {
                "status": "failed_operator_gates",
                "source": source,
                "old_psi_child": args.old_psi_child,
                "gates": gates,
            },
        )
        raise RuntimeError(f"Sequential calendar operator failed exact nesting: {gates}")
    print(
        "OPERATOR_GATES_PASSED "
        f"one_step_l1={gates['one_step_constant_path_nesting_l1']:.3e} "
        f"mature_error={gates['mature_flow_abs_error']:.3e}",
        flush=True,
    )

    E_old = float(old_solution.entry_rate)
    B_old = float(old_solution.entrants_mature_total)
    outside_share = float(args.outside_origin_entry_share)
    retention = (1.0 - outside_share) * E_old / B_old
    if not 0.0 <= retention <= 1.0:
        raise RuntimeError(
            "Initial outside-origin share implies an invalid retention rate: "
            f"rho={retention:.8g}"
        )
    conversion = float(old_parameters.entrant_conversion_factor)
    B_old_raw = B_old / conversion

    acs_age_reweight = None
    if args.historical_start_year == 2007:
        stationary_age_mass = np.sum(
            initial_g_pre, axis=(0, 1, 2, 4, 5, 6)
        )
        model_ages = (
            float(old_parameters.age_start)
            + np.arange(int(old_parameters.J), dtype=float)
            * float(old_parameters.da)
        )
        acs_age_reweight = acs_2007_age_reweight_diagnostic(
            stationary_age_mass,
            model_ages,
            E_old,
            periods=4,
            period_years=float(old_parameters.period_years),
        )
        write_json(
            outdir / "acs_2007_head_age_reweight_diagnostic.json",
            acs_age_reweight,
        )
        if args.initial_age_profile == "acs-national-heads-2007":
            initial_g_pre = reweight_distribution_to_acs_2007_ages(
                initial_g_pre,
                acs_age_reweight,
            )

    supply_rule, supply_normalization = calendar.normalize_date0_housing_supply(
        initial_g_pre,
        old_policy,
        old_parameters,
        b_grid,
        old_shared,
        str(args.housing_supply_mode),
    )
    retained_supply_elasticity = float(supply_rule.elasticity)
    if args.housing_supply_elasticity is not None:
        supply_rule = calendar.HousingSupplyRule(
            mode=supply_rule.mode,
            initial_price=supply_rule.initial_price,
            initial_stock=supply_rule.initial_stock,
            elasticity=float(args.housing_supply_elasticity),
        )
    supply_normalization["retained_housing_supply_elasticity"] = (
        retained_supply_elasticity
    )
    supply_normalization["transition_housing_supply_elasticity"] = float(
        supply_rule.elasticity
    )
    supply_normalization["elasticity_status"] = (
        "diagnostic_override"
        if args.housing_supply_elasticity is not None
        else "retained_model_value"
    )

    population_target_indices = (
        census_household_target_indices()
        if args.household_count_path == "census-hh3-2007-2023"
        else None
    )
    population_age_target_years = (
        census_household_age_target_years()
        if args.household_count_path == "census-hh3-2007-2023"
        else None
    )

    postshock_parameters = copy.deepcopy(old_parameters)
    postshock_parameters.psi_child = float(args.new_psi_child)
    postshock_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    counter = calendar.SolveCounter()
    if args.renewal_clock == "birth-vintage":
        momentum_suffix = (
            ""
            if math.isclose(
                float(args.initial_birth_pipeline_multiplier),
                1.0,
                rel_tol=0.0,
                abs_tol=1e-14,
            )
            else "_inherited_momentum"
        )
        census_suffix = (
            "_census_household_bridge"
            if args.household_count_path == "census-hh3-2007-2023"
            else ""
        )
        age_suffix = (
            "_acs_age_initial"
            if args.initial_age_profile == "acs-national-heads-2007"
            else ""
        )
        policy_suffix = (
            "_dependent_child_ltv95"
            if args.policy_case == "dependent-child-ltv95"
            else ""
        )
        scenario = (
            "preference_decline_open_birth_vintage_"
            + str(args.housing_supply_mode).replace("-", "_")
            + momentum_suffix
            + census_suffix
            + age_suffix
            + policy_suffix
        )
        paths, ages, children, scenario_summary = run_birth_vintage_scenario(
            label=scenario,
            initial_g_pre=initial_g_pre,
            baseline_price=float(np.asarray(old_price).reshape(-1)[0]),
            baseline_B=B_old,
            baseline_births=float(old_solution.total_births_kfe),
            parameters=postshock_parameters,
            b_grid=b_grid,
            periods=int(args.periods),
            outside_flow=outside_share * E_old,
            retention=retention,
            conversion=conversion,
            delay_periods=int(args.birth_to_entry_delay_periods),
            initial_birth_pipeline_multiplier=float(
                args.initial_birth_pipeline_multiplier
            ),
            population_target_indices=population_target_indices,
            population_age_target_years=population_age_target_years,
            policy_case=str(args.policy_case),
            policy_start_period=int(args.policy_start_period),
            old_psi_child=float(args.old_psi_child),
            new_psi_child=float(args.new_psi_child),
            preference_transition_periods=int(args.preference_transition_periods),
            supply_rule=supply_rule,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            outdir=outdir,
            counter=counter,
        )
    else:
        scenario = "preference_decline_open_household_state_fixed_stock"
        paths, ages, children, scenario_summary = calendar.run_scenario(
            scenario,
            "nesting",
            initial_g_pre,
            None,
            float(np.asarray(old_price).reshape(-1)[0]),
            E_old,
            B_old_raw,
            postshock_parameters,
            b_grid,
            int(args.periods),
            retention,
            conversion,
            supply_rule,
            float(args.market_tol),
            int(args.market_max_iter),
            outdir,
            counter,
        )
    calendar.write_csv(outdir / "transition_path.csv", paths)
    calendar.write_csv(outdir / "adult_age_distribution.csv", ages)
    calendar.write_csv(outdir / "child_pipeline.csv", children)
    calendar.make_plots(paths, ages, children, outdir)

    stationary_endpoint = None
    if args.renewal_clock == "birth-vintage" and not args.smoke:
        print("STATIONARY_ENDPOINT_START", flush=True)
        stationary_endpoint = solve_stationary_open_endpoint(
            chain=chain,
            base_overrides=base,
            old_price=float(np.asarray(old_price).reshape(-1)[0]),
            new_psi_child=float(args.new_psi_child),
            outside_flow=outside_share * E_old,
            retention=retention,
            conversion=conversion,
            maturation_survival_yield=float(
                scenario_summary["maturation_survival_yield"]
            ),
            supply_rule=supply_rule,
            policy_case=str(args.policy_case),
            outdir=outdir,
        )
        if stationary_endpoint.get("status") == "complete":
            print(
                "STATIONARY_ENDPOINT_DONE "
                f"p={stationary_endpoint['asset_price']:.6f} "
                f"S={stationary_endpoint['stationary_population_scale']:.6f}",
                flush=True,
            )
        else:
            print(
                "STATIONARY_ENDPOINT_NO_POSITIVE_ROOT "
                f"mode={supply_rule.mode}",
                flush=True,
            )
    make_paper_transition_plot(
        paths,
        old_births=float(old_solution.total_births_kfe),
        old_B=B_old,
        stationary_endpoint=stationary_endpoint,
        delay_years=(
            float(scenario_summary.get("birth_to_entry_delay_years", 0.0))
            if args.renewal_clock == "birth-vintage"
            else 0.0
        ),
        outdir=outdir,
    )
    if population_target_indices:
        make_population_target_accounting_plot(
            paths,
            post_target_outside_flow=outside_share * E_old,
            outdir=outdir,
        )

    historical_comparison = None
    if args.historical_start_year is not None:
        historical_comparison = make_historical_comparison(
            paths,
            start_year=int(args.historical_start_year),
            old_user_cost=float(old_parameters.user_cost_rate)
            * float(np.asarray(old_price).reshape(-1)[0]),
            household_counts_by_year=(
                census_household_model_age_levels()
                if args.household_count_path == "census-hh3-2007-2023"
                else None
            ),
            outdir=outdir,
        )

    old_row = closure_audit.readout(
        chain,
        old_solution,
        old_parameters,
        old_price,
        label="old_stationary_equilibrium",
        price_ratio=1.0,
        psi_child=float(args.old_psi_child),
        elapsed=old_seconds,
    )
    final_row = paths[-1]
    peak_population = max(paths, key=lambda row: float(row["adult_population"]))
    peak_price = max(paths, key=lambda row: float(row["asset_price"]))
    outside_flow = outside_share * E_old
    summary = {
        "status": "complete",
        "method": "exact cohort accounting with temporary-equilibrium household policies",
        "not_perfect_foresight": True,
        "source": source,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_arm": metadata.get("arm", metadata.get("winner_metadata", {}).get("arm")),
        "theta": theta,
        "old_psi_child": float(args.old_psi_child),
        "new_psi_child": float(args.new_psi_child),
        "preference_transition_periods": int(args.preference_transition_periods),
        "preference_transition_years": int(args.preference_transition_periods)
        * float(old_parameters.period_years),
        "renewal_clock": args.renewal_clock,
        "birth_to_entry_delay_periods": int(args.birth_to_entry_delay_periods),
        "initial_birth_pipeline_multiplier": float(
            args.initial_birth_pipeline_multiplier
        ),
        "initial_age_profile": str(args.initial_age_profile),
        "household_count_path": str(args.household_count_path),
        "housing_supply_elasticity_override": args.housing_supply_elasticity,
        "policy_case": str(args.policy_case),
        "policy_start_period": int(args.policy_start_period),
        "period_years": float(old_parameters.period_years),
        "periods": int(args.periods),
        "old_steady_state": old_row,
        "renewal_closure": {
            "equation": (
                "G_(t+1) is reweighted to HH-3 totals and national ACS "
                "four-year head-age shares through 2023; "
                "E_(t+1) = M_bar + rho * B_t afterward"
                if population_target_indices
                else "E_(t+1) = M_bar + rho * B_t"
            ),
            "outside_origin_entry_share_at_old_steady_state": outside_share,
            "outside_entry_flow_M": outside_flow,
            "retention_rho": retention,
            "entrant_conversion_factor": conversion,
            "old_entry_flow_E": E_old,
            "old_effective_mature_flow_B": B_old,
            "identity_residual": E_old - outside_flow - retention * B_old,
        },
        "operator_gates": gates,
        "supply_normalization": supply_normalization,
        "scenario_summary": scenario_summary,
        "stationary_open_endpoint": stationary_endpoint,
        "historical_comparison": historical_comparison,
        "acs_2007_head_age_reweight_diagnostic": acs_age_reweight,
        "peak_population": peak_population,
        "peak_asset_price": peak_price,
        "last_simulated_period": final_row,
        "model_solve_count": counter.total + 1,
        "elapsed_seconds": time.perf_counter() - started,
        "interpretation": {
            "household_expectations": "current price and post-shock primitives are treated as permanent at each date",
            "housing_supply": (
                "fixed at the old steady-state occupied stock"
                if args.housing_supply_mode == "fixed-stock"
                else "date-0-normalized contemporaneous constant-elastic supply"
            ),
            "population_normalization": "none after the old steady state",
            "initial_age_profile": str(args.initial_age_profile),
            "household_count_bridge": (
                "Census HH-3 totals and national ACS four-year householder ages "
                "matched through a dated reduced-form formation/migration residual"
                if population_target_indices
                else "none"
            ),
            "policy": str(args.policy_case),
            "initial_child_pipeline": (
                "old steady-state pipeline"
                if math.isclose(
                    float(args.initial_birth_pipeline_multiplier),
                    1.0,
                    rel_tol=0.0,
                    abs_tol=1e-14,
                )
                else "diagnostically scaled inherited cohorts; future cohorts are model-generated"
            ),
            "purpose": (
                "measure housing and policy incidence after matching the observed "
                "household-count and four-year head-age path with an accounting residual"
                if population_target_indices
                else "test whether inherited cohorts can generate short-run population and price momentum after a fertility-preference decline"
            ),
        },
    }
    write_json(outdir / "summary.json", summary)

    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    readme = f"""# Sequential-model open population transition

This packet starts from the paper's sequential child-room-floor model at
`psi_child={args.old_psi_child:g}` and permanently moves to
`psi_child={args.new_psi_child:g}`. The preference transition lasts
`{args.preference_transition_periods}` model periods; zero means an immediate
shift. The initial age profile is `{args.initial_age_profile}` and the
household-count path is `{args.household_count_path}`. With the Census bridge,
an age-specific reduced-form household-formation/migration residual matches
HH-3 household totals and national ACS four-year head-age shares through 2023;
afterward
the model returns to the old steady-state outside flow. The bridge preserves
the model's full conditional state distribution within each four-year age cell.
The observed bridge is used only over the reported historical window; it is not
a structural post-2023 headship law.
The retention rate is fixed throughout, and the full household distribution is
carried forward without population renormalization.

Housing supply is `{args.housing_supply_mode}`. The static-elastic option is a
sequence of contemporaneous supply equilibria; it is not a construction-sector
transition with a predetermined stock. Its transition elasticity is
`{supply_rule.elasticity:g}`; the retained model value is
`{retained_supply_elasticity:g}`.

The renewal clock is `{args.renewal_clock}`. The paper-facing default uses a
separate birth-vintage queue: household children-at-home states continue to
govern housing choices, but locally born entrant households arrive only after
the declared child-to-entry delay. The date-0 inherited pipeline multiplier is
`{args.initial_birth_pipeline_multiplier:g}`; values other than one are
diagnostic initial conditions, not re-estimated household parameters.
When the historical start is 2007, the packet also writes a national ACS
head-age reweighting check with the old steady-state entrant flow held fixed.
The optional tenure policy is `{args.policy_case}` and starts in period
`{args.policy_start_period}`.

The calculation is a temporary-equilibrium transition diagnostic: households
treat the current price as permanent at each date. It is not yet the paper's
perfect-foresight asset-price transition and should not be used for welfare.
At every date, the maintained stationary user-cost identity keeps the model
house-price-to-rent ratio fixed. The optional historical comparison reports
that restriction directly instead of introducing an unestimated wedge.

Exact command:

```bash
{command}
```

Key acceptance checks are in `summary.json`; the exact path is in
`transition_path.csv` and the standard visual packet is regenerated by the
same command.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print("TRANSITION_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
