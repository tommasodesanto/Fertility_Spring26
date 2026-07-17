#!/usr/bin/env python3
"""One fully re-estimated chain for the intergenerational bequest/exit battery.

Every calibration arm re-estimates the eleven clean-frontier coordinates. Most
headline bequest arms additionally estimate ``theta0`` while ``theta1`` and the
age-dependent owner borrowing schedule are externally fixed. A5 is retained
only as an explicitly underidentified diagnostic with both bequest parameters
free. M2 combines externally pinned age survival with a normalized child-blind
warm glow and no owner-LTV taper. M4 uses the same mechanism but estimates both
remaining De Nardi parameters, ``theta0`` and ``theta1``, with ``theta_n=0``.
M4_PROFILE is a diagnostic conditional profile that fixes ``theta1`` while
re-estimating the eleven clean-frontier coordinates and ``theta0``.
M5 is the income-disciplined variant of M4: it re-estimates the eleven
clean-frontier coordinates plus ``theta0``, ``theta1``, and
``tenure_choice_kappa`` with ``theta_n=0``, under literature-anchored income
persistence, the model-entry-age (18-24) entrant wealth distribution, and the
income-disciplined target system.
The script writes a recoverable record after every tight model solve.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from intergen_housing_fertility.calibration import (
    base_overrides,
    diagnostic_loss,
    external_entry_wealth_overrides_1824,
    extract_moments,
    get_target_set,
)
from intergen_housing_fertility.local_panel import income_process_overrides, jsonable
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
    validate_production_profile,
)
from intergen_housing_fertility.solver import InfeasibleThetaError, run_model_cp_dt
from tools.run_intergen_combined_wide_discovery import inverse, transform


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SEED = (
    ROOT
    / "output/model/full_calibration_audit_20260713_claude/round3/remote_mirror/sweep"
    / "clean_k0_frontier__tight_rep1.json"
)
DEFAULT_ACCEPTANCE = (
    ROOT
    / "code/data/mms_center_periphery/output_ownership_audit"
    / "acs_ownership_4year_acceptance_bins_6284.csv"
)
MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
# Arm M5 literature anchor in the Sommer-Sullivan tradition: annual AR(1)
# persistence 0.90 with annual innovation s.d. 0.20.  The in-repo PSID
# re-estimates (block1_income_process_estimates.csv in code/data/
# psid_followup_mar2026/output/intergen_income_entry_targets_20260716/)
# bracket these values -- AR(1)-only rho 0.9436 / s.d. 0.2846 and
# persistent-plus-transitory rho 0.9749 / s.d. 0.1769 -- and are reported as
# appendix evidence, not used directly.
# July-16 feasibility probes (m5 contract, section: income) rejected the
# SS-range values (0.90, sd 0.12-0.20): interior dead-node mass at age 22
# under the hard-debt/no-default structure. M5 therefore retains the matched
# process; the income-risk upgrade is deferred to M6 with a designed
# forbearance/default margin.
SS_ANNUAL_RHO = MATCHED_ANNUAL_RHO
SS_ANNUAL_INNOVATION_SD = MATCHED_ANNUAL_INNOVATION_SD
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0
PERIOD_YEARS = 4.0
SSA_2023_LIFE_TABLE_URL = "https://www.ssa.gov/oact/STATS/table4c6.html"
# Four-year survival from the combined male+female survivor counts in the
# 2023 SSA period life table.  Mortality starts at 66 so all children born by
# the model's final fertile age have matured before household death is active.
POSTRETIREMENT_SURVIVAL = {
    66.0: 0.9391263063710125,
    70.0: 0.9184976343249724,
    74.0: 0.8849521927812863,
    78.0: 0.8300468061015381,
}

# This is the audited relaxed discovery domain, not the narrower production
# box.  All eleven coordinates are re-estimated in every arm.
BASE_DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    ("beta_annual", 0.80, 0.9995, "discount"),
    ("alpha_cons", 0.02, 0.98, "logit"),
    ("c_bar_0", 0.0, 2.0, "softzero"),
    ("c_bar_n", 0.0, 3.0, "softzero"),
    ("h_bar_0", 0.05, 5.80, "log"),
    ("h_bar_jump", 0.0, 8.0, "softzero"),
    ("h_bar_n", 0.0, 5.0, "softzero"),
    ("psi_child", -3.0, 3.0, "asinh"),
    ("kappa_fert", 0.02, 50.0, "log"),
    ("chi", 0.10, 5.0, "log"),
    ("H0", 0.20, 80.0, "log"),
)
THETA0_DOMAIN = ("theta0", 0.0, 1.5, "softzero")
THETA1_DOMAIN = ("theta1", 0.10, 2.0, "log")
THETA_N_DOMAIN = ("theta_n", 0.0, 1.5, "softzero")
M4_THETA0_DOMAIN = ("theta0", 0.0, 8.0, "softzero")
M4_THETA1_DOMAIN = ("theta1", 0.02, 16.0, "log")
M5_KAPPA_DOMAIN = ("tenure_choice_kappa", 0.0, 0.12, "softzero")
SCHEDULE_ARMS = {"A2", "A3", "A4", "A5", "A3_PROFILE"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-record", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--seed-arm", type=str, default=None)
    parser.add_argument("--acceptance-csv", type=Path, default=DEFAULT_ACCEPTANCE)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--arm",
        choices=(
            "A0",
            "A1",
            "A2",
            "A3",
            "A4",
            "A5",
            "A3_PROFILE",
            "M0",
            "M1",
            "M2",
            "M3",
            "M4",
            "M4_PROFILE",
            "M5",
        ),
        required=True,
    )
    parser.add_argument("--ltv-terminal", type=float, default=0.4)
    parser.add_argument("--theta1", type=float, default=0.25)
    parser.add_argument("--seed-theta0", type=float, default=0.30)
    parser.add_argument("--seed-theta-n", type=float, default=0.75)
    parser.add_argument("--seed-kappa", type=float, default=0.0)
    parser.add_argument("--fixed-theta0", type=float, default=None)
    parser.add_argument("--seed", type=int, default=2026071401)
    parser.add_argument("--start-mix", type=float, default=0.0)
    parser.add_argument("--max-evals", type=int, default=2000)
    parser.add_argument("--minutes", type=float, default=230.0)
    parser.add_argument("--initial-step", type=float, default=0.015)
    parser.add_argument("--shrink", type=float, default=0.5)
    parser.add_argument("--J", type=int, default=PRODUCTION_J)
    parser.add_argument("--Nb", type=int, default=PRODUCTION_SEARCH_NB)
    parser.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    parser.add_argument("--tol-eq", type=float, default=1e-4)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def load_theta(path: Path, seed_arm: str | None = None) -> dict[str, float]:
    payload: Any = json.loads(path.read_text())
    if seed_arm is not None:
        winners = payload.get("winners") if isinstance(payload, dict) else None
        if not isinstance(winners, dict) or not isinstance(winners.get(seed_arm), dict):
            raise ValueError(f"{path} does not contain winner arm {seed_arm!r}")
        payload = winners[seed_arm]
    for key in ("record", "best"):
        if isinstance(payload, dict) and key in payload and isinstance(payload[key], dict):
            payload = payload[key]
    if isinstance(payload, dict) and "theta_used" in payload:
        payload = payload["theta_used"]
    elif isinstance(payload, dict) and "theta" in payload:
        payload = payload["theta"]
    if not isinstance(payload, dict):
        raise ValueError(f"{path} does not contain a theta object")
    theta = {str(k): float(v) for k, v in payload.items() if isinstance(v, (int, float))}
    required = {"beta" if name == "beta_annual" else name for name, *_ in BASE_DOMAIN}
    missing = required - set(theta)
    if missing:
        raise ValueError(f"seed record missing coordinates: {sorted(missing)}")
    return theta


def target_system(target_set: str = PRODUCTION_TARGET_SET) -> tuple[dict[str, float], dict[str, float]]:
    targets, weights = get_target_set(target_set)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = ROOMS_WEIGHT
    expected = 14 if target_set == "candidate_replacement_bequest_median_composition_v1" else 15
    if len(targets) != expected or set(targets) != set(weights):
        raise ValueError(f"battery requires the audited {expected}-moment target system")
    return targets, weights


def arm_contract(args: argparse.Namespace) -> tuple[list[tuple[str, float, float, str]], dict[str, float], dict[str, Any]]:
    if not 0.0 <= float(args.ltv_terminal) <= 1.0:
        raise ValueError("ltv-terminal multiplier must lie in [0,1]")
    if float(args.theta1) <= 0.0:
        raise ValueError("theta1 must be positive")
    if args.arm == "M4" and not (
        M4_THETA1_DOMAIN[1] <= float(args.theta1) <= M4_THETA1_DOMAIN[2]
    ):
        raise ValueError("M4 theta1 start must lie inside [0.02,16]")
    if args.arm == "M4_PROFILE" and not (
        M4_THETA1_DOMAIN[1] <= float(args.theta1) <= M4_THETA1_DOMAIN[2]
    ):
        raise ValueError("M4 profile theta1 value must lie inside [0.02,16]")
    if args.arm == "M4" and not (
        M4_THETA0_DOMAIN[1] <= float(args.seed_theta0) <= M4_THETA0_DOMAIN[2]
    ):
        raise ValueError("M4 theta0 start must lie inside [0,8]")
    if args.arm == "M5":
        if not M4_THETA1_DOMAIN[1] <= float(args.theta1) <= M4_THETA1_DOMAIN[2]:
            raise ValueError("M5 theta1 start must lie inside [0.02,16]")
        if not M4_THETA0_DOMAIN[1] <= float(args.seed_theta0) <= M4_THETA0_DOMAIN[2]:
            raise ValueError("M5 theta0 start must lie inside [0,8]")
        if not M5_KAPPA_DOMAIN[1] <= float(args.seed_kappa) <= M5_KAPPA_DOMAIN[2]:
            raise ValueError("M5 tenure_choice_kappa start must lie inside [0,0.12]")
    active = list(BASE_DOMAIN)
    fixed: dict[str, float] = {"theta_n": 0.0, "tenure_choice_kappa": 0.0}
    if args.arm == "M3":
        active.extend((THETA0_DOMAIN, THETA1_DOMAIN, THETA_N_DOMAIN))
        fixed.pop("theta_n")
    elif args.arm == "M4":
        active.extend((M4_THETA0_DOMAIN, M4_THETA1_DOMAIN))
    elif args.arm == "M5":
        active.extend((M4_THETA0_DOMAIN, M4_THETA1_DOMAIN, M5_KAPPA_DOMAIN))
        fixed.pop("tenure_choice_kappa")
    elif args.arm == "M4_PROFILE":
        active.append(M4_THETA0_DOMAIN)
        fixed["theta1"] = float(args.theta1)
    elif args.arm in {"A1", "A3", "A4", "M2"}:
        active.append(THETA0_DOMAIN)
        fixed["theta1"] = float(args.theta1)
    elif args.arm == "A5":
        active.extend((THETA0_DOMAIN, THETA1_DOMAIN))
    elif args.arm == "A3_PROFILE":
        if args.fixed_theta0 is None or not 0.0 <= float(args.fixed_theta0) <= THETA0_DOMAIN[2]:
            raise ValueError("A3_PROFILE requires fixed-theta0 inside the theta0 domain")
        fixed.update(theta0=float(args.fixed_theta0), theta1=float(args.theta1))
    else:
        fixed.update(theta0=0.0, theta1=float(args.theta1))

    bequest_spec = "linear_child_scale"
    if args.arm in {"A1", "A3", "A5", "A3_PROFILE"}:
        bequest_spec = "parent_gated_luxury"
    elif args.arm == "A4":
        bequest_spec = "equal_division_luxury"
    mechanism = {
        "bequest_spec": bequest_spec,
        "normalize_bequest_utility": True,
        "owner_ltv_taper": args.arm in SCHEDULE_ARMS,
        "owner_ltv_taper_start_age": 66.0,
        "owner_ltv_taper_end_age": 82.0,
        "owner_ltv_terminal_share": float(args.ltv_terminal),
        "use_age_survival": args.arm in {"M1", "M2", "M3", "M4", "M4_PROFILE", "M5"},
        # M5 only: entrant debt beyond the state-conditional feasible frontier
        # is censored to the nearest alive node (relocated, never left dead);
        # the July-11 entry gate stays active as the backstop.
        "entry_wealth_censor_to_frontier": args.arm == "M5",
    }
    if args.arm in {"M0", "M1", "M2", "M3", "M4", "M4_PROFILE", "M5"}:
        mechanism.update(
            owner_ltv_taper=False,
            owner_ltv_terminal_share=0.0,
            bequest_spec="linear_child_scale",
            normalize_bequest_utility=True,
        )
    return active, fixed, mechanism


def survival_schedule(args: argparse.Namespace) -> np.ndarray:
    schedule = np.ones(int(args.J) - 1, dtype=float)
    if args.arm not in {"M1", "M2", "M3", "M4", "M4_PROFILE", "M5"}:
        return schedule
    for j in range(int(args.J) - 1):
        age = float(18.0 + j * PERIOD_YEARS)
        if age in POSTRETIREMENT_SURVIVAL:
            schedule[j] = POSTRETIREMENT_SURVIVAL[age]
    return schedule


def theta_from_unit(
    unit: np.ndarray,
    active: list[tuple[str, float, float, str]],
    fixed: dict[str, float],
) -> dict[str, float]:
    theta = dict(fixed)
    for u, (name, lo, hi, kind) in zip(np.asarray(unit, dtype=float), active):
        value = transform(float(u), lo, hi, kind)
        theta["beta" if name == "beta_annual" else name] = value**PERIOD_YEARS if name == "beta_annual" else value
    return theta


def unit_from_theta(theta: dict[str, float], active: list[tuple[str, float, float, str]]) -> np.ndarray:
    unit = []
    for name, lo, hi, kind in active:
        key = "beta" if name == "beta_annual" else name
        value = float(theta.get(key, 0.0))
        if name == "beta_annual":
            value = value ** (1.0 / PERIOD_YEARS)
        unit.append(inverse(value, lo, hi, kind))
    return np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)


def common_overrides(args: argparse.Namespace, mechanism: dict[str, Any]) -> dict[str, Any]:
    if args.arm == "M5":
        income = income_process_overrides(
            5, "rouwenhorst", SS_ANNUAL_INNOVATION_SD, SS_ANNUAL_RHO
        )
    else:
        income = income_process_overrides(
            5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
        )
    overrides = {
        **base_overrides(J=args.J, Nb=args.Nb, n_house=5, max_iter_eq=args.max_iter_eq),
        **production_profile_overrides(),
        **income,
        **mechanism,
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
        "eta_supply": np.array([1.75]),
        "lambda_d": 0.0,
        "debt_taper_start_age": 42.0,
        "debt_taper_end_age": 62.0,
        "survival_probs": survival_schedule(args),
    }
    if args.arm == "M5":
        overrides.update(external_entry_wealth_overrides_1824())
    return overrides


def load_acceptance(path: Path) -> list[dict[str, Any]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty acceptance file: {path}")
    return rows


def lifecycle_diagnostics(sol: Any, P: Any, acceptance: list[dict[str, Any]]) -> dict[str, Any]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float)
    mass_age = g.sum(axis=(0, 1, 2, 4, 5, 6))
    own_age = g[:, 1:].sum(axis=(0, 1, 2, 4, 5, 6)) / np.maximum(mass_age, 1e-300)
    liquid_age = np.einsum("b,btijznc->j", b_grid, g) / np.maximum(mass_age, 1e-300)
    ages = np.asarray([float(P.age_start + j * P.da) for j in range(g.shape[3])])
    tenure_rooms = np.concatenate(([0.0], np.asarray(P.H_own, dtype=float)))
    owner_rooms_age = np.einsum("t,btijznc->j", tenure_rooms, g)
    renter_rooms_age = np.einsum("bijznc,bijznc->j", g[:, 0], np.asarray(sol.hR_pol)[:, 0])
    occupied_rooms_age = owner_rooms_age + renter_rooms_age

    checks: list[dict[str, Any]] = []
    for row in acceptance:
        state_age = float(row["age_state"])
        idx = int(np.argmin(np.abs(ages - state_age)))
        model = float(own_age[idx])
        lower, upper = float(row["acceptance_lower"]), float(row["acceptance_upper"])
        checks.append(
            {
                "age_state": state_age,
                "model": model,
                "empirical": float(row["owner_rate"]),
                "lower": lower,
                "upper": upper,
                "within_band": bool(lower <= model <= upper),
            }
        )
    adjacent = [checks[i + 1]["model"] - checks[i]["model"] for i in range(len(checks) - 1)]
    cliff = any(delta < -0.15 for delta in adjacent)

    def window(values: np.ndarray, lo: float, hi: float) -> float:
        selected = (ages >= lo) & (ages <= hi)
        weights = mass_age[selected]
        return float(np.sum(values[selected] * weights) / max(float(np.sum(weights)), 1e-300))

    early = window(liquid_age, 62.0, 74.0)
    late = window(liquid_age, 74.0, 82.0)
    old = ages >= 62.0
    return {
        "ages": ages.tolist(),
        "ownership_by_age": own_age.tolist(),
        "mean_liquid_wealth_by_age": liquid_age.tolist(),
        "ownership_acceptance_checks": checks,
        "adjacent_ownership_changes": adjacent,
        "ownership_cliff_failure": cliff,
        "ownership_band_failure": not all(row["within_band"] for row in checks),
        "hard_acceptance_pass": bool(not cliff and all(row["within_band"] for row in checks)),
        "wealth_62_74": early,
        "wealth_74_82": late,
        "decum_ratio_wealth_74plus_over_62_74": late / max(early, 1e-12),
        "age_mass": mass_age.tolist(),
        "owner_rooms_by_age": owner_rooms_age.tolist(),
        "occupied_rooms_by_age": occupied_rooms_age.tolist(),
        "old_household_mass_share_62plus": float(np.sum(mass_age[old]) / max(np.sum(mass_age), 1e-300)),
        "old_owner_rooms_share_62plus": float(
            np.sum(owner_rooms_age[old]) / max(np.sum(owner_rooms_age), 1e-300)
        ),
        "old_occupied_rooms_share_62plus": float(
            np.sum(occupied_rooms_age[old]) / max(np.sum(occupied_rooms_age), 1e-300)
        ),
        "survival_probs": np.asarray(P.survival_probs, dtype=float).tolist(),
    }


def target_fit(
    moments: dict[str, float], targets: dict[str, float], weights: dict[str, float]
) -> list[dict[str, float | str]]:
    rows = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights[name])
        rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )
    return rows


def parameter_table(
    theta: dict[str, float],
    active: list[tuple[str, float, float, str]],
    fixed: dict[str, float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, lo, hi, kind in active:
        key = "beta" if name == "beta_annual" else name
        estimate = float(theta[key]) ** (1.0 / PERIOD_YEARS) if name == "beta_annual" else float(theta[key])
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower": lo,
                "upper": hi,
                "transform": kind,
                "role": "estimated",
                "near_bound": min(estimate - lo, hi - estimate) <= 0.02 * (hi - lo),
            }
        )
    for name, estimate in sorted(fixed.items()):
        rows.append(
            {
                "parameter": name,
                "estimate": float(estimate),
                "lower": None,
                "upper": None,
                "transform": "external restriction",
                "role": "fixed",
                "near_bound": False,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25
        # Complete the transformed initial simplex so smoke exercises theta0
        # and every active bequest coordinate, not only the first core coordinate.
        args.max_evals = min(int(args.max_evals), 15)
        args.minutes = min(float(args.minutes), 8.0)
    else:
        validate_production_profile(
            PRODUCTION_PROFILE_NAME,
            J=args.J,
            Nb=args.Nb,
            n_house=5,
            income_states=5,
            target_set=PRODUCTION_TARGET_SET,
            max_iter_eq=args.max_iter_eq,
            stage="search",
        )
        if not math.isclose(float(args.tol_eq), 1e-4):
            raise ValueError(
                "proper bequest/exit searches require tol_eq=1e-4; each winner "
                "is re-solved twice at max_iter_eq=40 and tol_eq=2.5e-5"
            )

    seed_theta = load_theta(args.seed_record, args.seed_arm)
    active, fixed, mechanism = arm_contract(args)
    seed_theta.update(fixed)
    if args.arm == "M3":
        seed_theta.update(
            theta0=float(args.seed_theta0),
            theta1=float(args.theta1),
            theta_n=float(args.seed_theta_n),
        )
    elif args.arm == "M4":
        # For M4 these are dispersed optimizer starts, never external values.
        seed_theta.update(theta0=float(args.seed_theta0), theta1=float(args.theta1))
    elif args.arm == "M5":
        # Dispersed optimizer starts; the nested reference sets theta0=0 and
        # tenure_choice_kappa=0 at the strict M1 winner coordinates.
        seed_theta.update(
            theta0=float(args.seed_theta0),
            theta1=float(args.theta1),
            tenure_choice_kappa=float(args.seed_kappa),
        )
    x0 = unit_from_theta(seed_theta, active)
    rng = np.random.default_rng(int(args.seed))
    start_mix = float(np.clip(args.start_mix, 0.0, 0.25))
    if start_mix > 0.0:
        x0 = np.clip((1.0 - start_mix) * x0 + start_mix * rng.random(len(active)), 0.0, 1.0)

    if args.arm == "M3":
        target_set = "candidate_replacement_bequest_internal_v1"
    elif args.arm in {"M4", "M4_PROFILE"}:
        target_set = "candidate_replacement_bequest_median_composition_v1"
    elif args.arm == "M5":
        target_set = "candidate_replacement_income_disciplined_v1"
    else:
        target_set = PRODUCTION_TARGET_SET
    targets, weights = target_system(target_set)
    acceptance = load_acceptance(args.acceptance_csv)
    overrides = common_overrides(args, mechanism)
    max_evals = max(1, int(args.max_evals))
    initial_step = float(np.clip(args.initial_step, 1e-5, 0.25))
    shrink = float(np.clip(args.shrink, 0.1, 0.9))
    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path = args.outdir / "cases.jsonl"
    best_path = args.outdir / "best.json"
    cases_path.write_text("")
    metadata = {
        "status": "smoke" if args.smoke else "proper_joint_smm_chain",
        "arm": args.arm,
        "identification_status": (
            "diagnostic_underidentified"
            if args.arm == "A5"
            else "conditional_theta1_profile"
            if args.arm == "M4_PROFILE"
            else "rank_to_be_verified_at_winner"
        ),
        "free_parameter_count": len(active),
        "target_count": len(targets),
        "active_domain": [
            {"name": name, "lower": lo, "upper": hi, "transform": kind}
            for name, lo, hi, kind in active
        ],
        "fixed_parameters": fixed,
        "mechanism": mechanism,
        "survival_source": (
            SSA_2023_LIFE_TABLE_URL
            if args.arm in {"M1", "M2", "M3", "M4", "M4_PROFILE", "M5"}
            else None
        ),
        "income_process": {
            "states": 5,
            "process": "rouwenhorst",
            "annual_rho": SS_ANNUAL_RHO if args.arm == "M5" else MATCHED_ANNUAL_RHO,
            "annual_innovation_sd": (
                SS_ANNUAL_INNOVATION_SD if args.arm == "M5" else MATCHED_ANNUAL_INNOVATION_SD
            ),
        },
        "entry_wealth_ages": "18_24" if args.arm == "M5" else "25_35",
        "seed_record": str(args.seed_record),
        "seed_arm": args.seed_arm,
        "seed": int(args.seed),
        "start_mix": start_mix,
        "initial_unit_vector": x0,
        "algorithm": "bounded_nelder_mead_in_transformed_coordinates",
        "minutes": float(args.minutes),
        "max_evals": max_evals,
        "J": int(args.J),
        "Nb": int(args.Nb),
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "tight_winner_evaluator": {"max_iter_eq": 40, "tol_eq": 2.5e-5, "repeats": 2},
        "target_set": target_set,
        "targets": targets,
        "weights": weights,
        "acceptance_file": str(args.acceptance_csv),
    }
    (args.outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True))

    started = time.perf_counter()
    # Preserve five minutes inside the declared chain budget for two fresh
    # tight winner solves.  Smoke tests exercise only the search loop.
    tight_reserve_sec = 0.0 if args.smoke else min(300.0, max(0.0, float(args.minutes) * 60.0 - 1.0))
    deadline = started + max(1.0, float(args.minutes) * 60.0 - tight_reserve_sec)
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    eval_idx = 0

    def score(record: dict[str, Any] | None) -> float:
        if record is None or not bool(record.get("strict_converged", False)):
            return math.inf
        return float(record.get("rank_loss", math.inf))

    def can_eval() -> bool:
        return eval_idx < max_evals and time.perf_counter() < deadline

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal eval_idx, best
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = theta_from_unit(clipped, active, fixed)
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt({**overrides, **theta}, verbose=False)
            moments = extract_moments(sol, P)
            loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
            residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
            timings = dict(getattr(sol, "timings", {}))
            strict = bool(
                timings.get("strict_converged", getattr(sol, "converged", False))
                and math.isfinite(residual)
                and residual <= float(P.tol_eq)
            )
            lifecycle = lifecycle_diagnostics(sol, P, acceptance)
            fit = target_fit(moments, targets, weights)
            status, error, census = "ok", "", []
            price = float(np.asarray(p_eq).reshape(-1)[0])
        except InfeasibleThetaError as exc:
            moments, fit, lifecycle, timings = {}, [], {}, {}
            loss, residual, price, strict = math.inf, math.inf, math.nan, False
            status, error, census = "infeasible_theta", str(exc), list(exc.census)
        except Exception as exc:  # noqa: BLE001 - each failed proposal is a checkpoint.
            moments, fit, lifecycle, timings = {}, [], {}, {}
            loss, residual, price, strict = math.inf, math.inf, math.nan, False
            status, error, census = f"failed:{type(exc).__name__}", str(exc), []
        record = {
            "case": eval_idx,
            "label": label,
            "arm": args.arm,
            "status": status,
            "strict_converged": strict,
            "rank_loss": loss,
            "market_residual": residual,
            "price": price,
            "theta": theta,
            "moments": moments,
            "target_fit": fit,
            "parameters": parameter_table(theta, active, fixed),
            "lifecycle": lifecycle,
            "timings": timings,
            "origin": origin,
            "unit_vector": clipped,
            "feasibility_census": census,
            "error": error,
            "elapsed_sec": time.perf_counter() - t0,
        }
        with cases_path.open("a") as handle:
            handle.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
        records.append(record)
        (args.outdir / "latest_completed_case.json").write_text(
            json.dumps(jsonable(record), indent=2, sort_keys=True)
        )
        if score(record) < score(best):
            best = record
            best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
        elapsed = time.perf_counter() - started
        print(
            f"eval {eval_idx + 1}/{max_evals} {label}: status={status} "
            f"loss={score(record):.7g} resid={residual:.2e} best={score(best):.7g} "
            f"elapsed={elapsed / 60.0:.1f}m",
            flush=True,
        )
        eval_idx += 1
        return record

    simplex: list[tuple[np.ndarray, dict[str, Any]]] = []
    simplex.append((x0.copy(), evaluate(x0, "seed", {"phase": "seed"})))
    for j in range(len(active)):
        if not can_eval():
            break
        step_j = initial_step * (0.75 + 0.5 * rng.random())
        trial = x0.copy()
        trial[j] = trial[j] + step_j if trial[j] + step_j <= 1.0 else max(0.0, trial[j] - step_j)
        simplex.append(
            (
                trial,
                evaluate(trial, f"nm_init_{j:02d}", {"phase": "initial_simplex", "dim": j, "step": step_j}),
            )
        )

    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, shrink
    iteration = 0
    while can_eval() and len(simplex) >= 2:
        simplex.sort(key=lambda item: score(item[1]))
        centroid = np.mean([unit for unit, _ in simplex[:-1]], axis=0)
        worst_unit, worst_record = simplex[-1]
        second_worst = score(simplex[-2][1])
        best_loss = score(simplex[0][1])
        reflected = np.clip(centroid + alpha * (centroid - worst_unit), 0.0, 1.0)
        reflected_record = evaluate(reflected, f"nm_reflect_{iteration:04d}", {"phase": "reflect", "iteration": iteration})
        reflected_loss = score(reflected_record)
        if reflected_loss < best_loss and can_eval():
            expanded = np.clip(centroid + gamma * (reflected - centroid), 0.0, 1.0)
            expanded_record = evaluate(expanded, f"nm_expand_{iteration:04d}", {"phase": "expand", "iteration": iteration})
            simplex[-1] = (expanded, expanded_record) if score(expanded_record) < reflected_loss else (reflected, reflected_record)
        elif reflected_loss < second_worst:
            simplex[-1] = (reflected, reflected_record)
        else:
            if reflected_loss < score(worst_record):
                contracted = np.clip(centroid + rho * (reflected - centroid), 0.0, 1.0)
                threshold, phase = reflected_loss, "outside_contract"
            else:
                contracted = np.clip(centroid - rho * (centroid - worst_unit), 0.0, 1.0)
                threshold, phase = score(worst_record), "inside_contract"
            contract_record = evaluate(contracted, f"nm_contract_{iteration:04d}", {"phase": phase, "iteration": iteration}) if can_eval() else worst_record
            if score(contract_record) < threshold:
                simplex[-1] = (contracted, contract_record)
            else:
                best_unit = simplex[0][0].copy()
                new_simplex = [simplex[0]]
                for slot, (unit, old_record) in enumerate(simplex[1:], start=1):
                    if not can_eval():
                        new_simplex.append((unit, old_record))
                        continue
                    shrunk = np.clip(best_unit + sigma * (unit - best_unit), 0.0, 1.0)
                    shrunk_record = evaluate(shrunk, f"nm_shrink_{iteration:04d}_{slot:02d}", {"phase": "shrink", "iteration": iteration, "slot": slot})
                    new_simplex.append((shrunk, shrunk_record))
                simplex = new_simplex
        iteration += 1

    tight_records: list[dict[str, Any]] = []
    if not args.smoke and best is not None:
        tight_path = args.outdir / "tight_cases.jsonl"
        tight_path.write_text("")
        tight_overrides = {**overrides, "max_iter_eq": 40, "tol_eq": 2.5e-5}
        winning_unit = np.asarray(best["unit_vector"], dtype=float)
        winning_theta = theta_from_unit(winning_unit, active, fixed)
        for repeat_idx in range(2):
            t0 = time.perf_counter()
            try:
                sol, P, p_eq = run_model_cp_dt({**tight_overrides, **winning_theta}, verbose=False)
                moments = extract_moments(sol, P)
                loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
                residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
                timings = dict(getattr(sol, "timings", {}))
                strict = bool(
                    timings.get("strict_converged", getattr(sol, "converged", False))
                    and math.isfinite(residual)
                    and residual <= 2.5e-5
                )
                lifecycle = lifecycle_diagnostics(sol, P, acceptance)
                fit = target_fit(moments, targets, weights)
                status, error, census = "ok", "", []
                price = float(np.asarray(p_eq).reshape(-1)[0])
            except InfeasibleThetaError as exc:
                moments, fit, lifecycle, timings = {}, [], {}, {}
                loss, residual, price, strict = math.inf, math.inf, math.nan, False
                status, error, census = "infeasible_theta", str(exc), list(exc.census)
            except Exception as exc:  # noqa: BLE001 - persist verification failures.
                moments, fit, lifecycle, timings = {}, [], {}, {}
                loss, residual, price, strict = math.inf, math.inf, math.nan, False
                status, error, census = f"failed:{type(exc).__name__}", str(exc), []
            tight_record = {
                "case": repeat_idx,
                "label": f"tight_winner_repeat_{repeat_idx + 1}",
                "arm": args.arm,
                "status": status,
                "strict_converged": strict,
                "rank_loss": loss,
                "market_residual": residual,
                "price": price,
                "theta": winning_theta,
                "moments": moments,
                "target_fit": fit,
                "parameters": parameter_table(winning_theta, active, fixed),
                "lifecycle": lifecycle,
                "timings": timings,
                "origin": {"phase": "tight_winner_verification", "search_case": best["case"]},
                "unit_vector": winning_unit,
                "feasibility_census": census,
                "error": error,
                "evaluator": {"max_iter_eq": 40, "tol_eq": 2.5e-5},
                "elapsed_sec": time.perf_counter() - t0,
            }
            with tight_path.open("a") as handle:
                handle.write(json.dumps(jsonable(tight_record), sort_keys=True) + "\n")
            tight_records.append(tight_record)
            (args.outdir / "latest_tight_case.json").write_text(
                json.dumps(jsonable(tight_record), indent=2, sort_keys=True)
            )
            print(
                f"tight repeat {repeat_idx + 1}/2: strict={strict} loss={loss:.7g} "
                f"resid={residual:.2e}",
                flush=True,
            )

    strict_tight = [row for row in tight_records if bool(row.get("strict_converged", False))]
    best_tight = min(strict_tight, key=lambda row: float(row["rank_loss"])) if strict_tight else None
    if best_tight is not None:
        (args.outdir / "best_tight.json").write_text(
            json.dumps(jsonable(best_tight), indent=2, sort_keys=True)
        )
    repeat_check = None
    if len(tight_records) == 2 and tight_records[0].get("moments") and tight_records[1].get("moments"):
        repeat_check = {
            "loss_abs_difference": abs(float(tight_records[0]["rank_loss"]) - float(tight_records[1]["rank_loss"])),
            "max_abs_moment_difference": max(
                abs(float(tight_records[0]["moments"][name]) - float(tight_records[1]["moments"][name]))
                for name in targets
            ),
            "both_strict": all(bool(row.get("strict_converged", False)) for row in tight_records),
        }

    summary = {
        "status": "complete" if eval_idx >= max_evals else "time_budget_reached",
        "best": best,
        "best_search": best,
        "best_tight": best_tight,
        "tight_repeats": tight_records,
        "tight_repeat_check": repeat_check,
        "n_cases_completed": len(records),
        "n_strict": sum(bool(row.get("strict_converged", False)) for row in records),
        "n_hard_acceptance": sum(bool(row.get("lifecycle", {}).get("hard_acceptance_pass", False)) for row in records),
        "elapsed_sec": time.perf_counter() - started,
        "search_elapsed_budget_sec": max(1.0, float(args.minutes) * 60.0 - tight_reserve_sec),
        "tight_reserve_sec": tight_reserve_sec,
        "iterations_completed": iteration,
        "metadata": metadata,
    }
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
