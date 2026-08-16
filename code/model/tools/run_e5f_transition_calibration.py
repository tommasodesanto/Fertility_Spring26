#!/usr/bin/env python3
"""Transition-based calibration readout for the sequential E5F model.

The stationary E5F estimate supplies the search center. Each panel candidate
derives an old steady-state preference intercept that matches completed
fertility, reweights the benchmark to the observed 2007 householder-age
distribution, moves preferences over four model periods, propagates the
unnormalized calendar distribution, and measures the existing E5 target system
at the 2023 date.
The Census/ACS household path is an explicit external bridge. The bounded
panel and coordinate designs are diagnostics, not claims of global optimality.
The 2007 trend origin is a timing normalization, not an estimated shock date.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
import time
import traceback
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

import audit_closed_reproductive_closure as closure
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system
from intergen_eqscale_seq_optimized.e5f_floor_profile import E5F_DOMAIN
from intergen_eqscale_seq_optimized.e5f_income_entry_profile import (
    E5F_INCOME_ENTRY_DOMAIN,
    E5F_INCOME_ENTRY_PROFILE_NAME,
    e5f_income_entry_overrides,
)


DEFAULT_OUTDIR = ROOT / "output/model/e5f_transition_calibration"
TARGET_AGE = 42.0
START_YEAR = 2007
TERMINAL_YEAR = 2023
TRANSITION_PERIODS = 4
TRANSITION_SEARCH_DOMAIN: tuple[tuple[str, float, float, str], ...] = tuple(
    row for row in E5F_DOMAIN if row[0] != "psi_child"
) + (("psi_child_2023", -1.25, 0.20, "asinh"),)
BASELINE_MODEL_PROFILE = "e5f-floor"
REPAIRED_MODEL_PROFILE = "e5f-income-entry"
DEFAULT_REPAIRED_PREFERENCE_CHANGE = -0.2902


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=closure.E5F_FLOOR_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--model-profile",
        choices=(BASELINE_MODEL_PROFILE, REPAIRED_MODEL_PROFILE),
        default=BASELINE_MODEL_PROFILE,
        help=(
            "Household profile. The repair keeps the E5F room floor, adds the "
            "externally measured permanent-income process, and estimates one "
            "first-birth fixed cost."
        ),
    )
    parser.add_argument("--old-psi-child", type=float, default=0.1062)
    parser.add_argument("--old-completed-fertility-target", type=float, default=2.12)
    parser.add_argument("--old-completed-fertility-tol", type=float, default=5e-4)
    parser.add_argument(
        "--candidate-new-psis",
        default="-0.2419,-0.10,0.0",
        help="Comma-separated terminal child-preference intercepts.",
    )
    parser.add_argument(
        "--candidate-psi-changes",
        default=None,
        help=(
            "Optional comma-separated changes from the normalized 2007 intercept. "
            "This is the invariant coordinate for the repaired profile."
        ),
    )
    parser.add_argument(
        "--fixed-first-birth-cost",
        type=float,
        default=None,
        help="Diagnostic repaired-profile cost held fixed outside a joint panel.",
    )
    parser.add_argument("--outside-origin-entry-share", type=float, default=0.169)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument("--nb", type=int, default=120)
    parser.add_argument("--post-2023-periods", type=int, default=0)
    parser.add_argument(
        "--policy-case",
        choices=("none", "dependent-child-ltv95"),
        default="none",
    )
    parser.add_argument("--expected-source-sha256", default=None)
    parser.add_argument("--expected-target-set", default=None)
    parser.add_argument("--expected-target-fingerprint", default=None)
    parser.add_argument("--panel-task-id", type=int, default=None)
    parser.add_argument("--panel-size", type=int, default=64)
    parser.add_argument("--panel-seed", type=int, default=2026081601)
    parser.add_argument(
        "--panel-design",
        choices=("mixed", "coordinate"),
        default="mixed",
    )
    parser.add_argument(
        "--panel-local-radius",
        type=float,
        default=0.18,
        help="Radius in transformed unit space for the local half of the panel.",
    )
    parser.add_argument(
        "--panel-center-json",
        type=Path,
        default=None,
        help="Optional collected best-candidate JSON used to center a refinement panel.",
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


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def parse_candidates(raw: str) -> list[float]:
    values = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("--candidate-new-psis must contain finite numbers")
    return list(dict.fromkeys(values))


def transform_unit(u: float, lower: float, upper: float, kind: str) -> float:
    u = float(np.clip(u, 0.0, 1.0))
    if kind == "log":
        return lower * (upper / lower) ** u
    if kind == "discount":
        return lower + (upper - lower) * (1.0 - (1.0 - u) ** 2)
    if kind == "asinh":
        return math.sinh(
            (1.0 - u) * math.asinh(lower) + u * math.asinh(upper)
        )
    if kind == "softzero":
        return lower + (upper - lower) * u**2
    raise ValueError(f"Unsupported transition-search transform: {kind}")


def inverse_unit(value: float, lower: float, upper: float, kind: str) -> float:
    value = float(np.clip(value, lower, upper))
    if kind == "log":
        return math.log(value / lower) / math.log(upper / lower)
    if kind == "discount":
        return 1.0 - math.sqrt(
            max(0.0, 1.0 - (value - lower) / (upper - lower))
        )
    if kind == "asinh":
        return (math.asinh(value) - math.asinh(lower)) / (
            math.asinh(upper) - math.asinh(lower)
        )
    if kind == "softzero":
        return math.sqrt(max(0.0, (value - lower) / (upper - lower)))
    raise ValueError(f"Unsupported transition-search transform: {kind}")


def latin_hypercube(size: int, dimensions: int, seed: int) -> np.ndarray:
    if size < 1 or dimensions < 1:
        raise ValueError("Latin-hypercube dimensions must be positive")
    rng = np.random.default_rng(int(seed))
    out = np.empty((size, dimensions), dtype=float)
    for column in range(dimensions):
        out[:, column] = (rng.permutation(size) + rng.random(size)) / size
    return out


def transition_unit_from_candidate(
    theta: dict[str, float], terminal_preference_coordinate: float
) -> np.ndarray:
    values: list[float] = []
    for name, lower, upper, kind in TRANSITION_SEARCH_DOMAIN:
        if name in {"psi_child_2023", "psi_child_change_2023"}:
            value = float(terminal_preference_coordinate)
        elif name == "beta_annual":
            value = float(theta["beta"]) ** 0.25
        else:
            value = float(theta[name])
        values.append(inverse_unit(value, lower, upper, kind))
    return np.clip(np.asarray(values, dtype=float), 0.0, 1.0)


def candidate_from_transition_unit(
    source_theta: dict[str, float], unit: np.ndarray
) -> tuple[dict[str, float], float]:
    theta = dict(source_theta)
    terminal_preference_coordinate = math.nan
    for u, (name, lower, upper, kind) in zip(
        np.asarray(unit, dtype=float), TRANSITION_SEARCH_DOMAIN, strict=True
    ):
        value = transform_unit(float(u), lower, upper, kind)
        if name in {"psi_child_2023", "psi_child_change_2023"}:
            terminal_preference_coordinate = value
        elif name == "beta_annual":
            theta["beta"] = value**4
        else:
            theta[name] = value
    if not math.isfinite(terminal_preference_coordinate):
        raise RuntimeError("Transition candidate omitted the terminal preference intercept")
    return theta, float(terminal_preference_coordinate)


def panel_candidate(
    source_theta: dict[str, float], args: argparse.Namespace
) -> tuple[dict[str, float], float, dict[str, Any]]:
    task = int(args.panel_task_id)
    size = int(args.panel_size)
    if not 1 <= task <= size:
        raise ValueError(f"--panel-task-id must lie in [1,{size}]")
    if not 0.0 < float(args.panel_local_radius) <= 0.5:
        raise ValueError("--panel-local-radius must lie in (0,0.5]")

    center_theta = dict(source_theta)
    terminal_name = TRANSITION_SEARCH_DOMAIN[-1][0]
    center_terminal_coordinate = (
        DEFAULT_REPAIRED_PREFERENCE_CHANGE
        if terminal_name == "psi_child_change_2023"
        else -0.2419
    )
    center_status = "working_stationary_estimate"
    center_sha256 = None
    if args.panel_center_json is not None:
        center_path = args.panel_center_json.resolve()
        payload = json.loads(center_path.read_text(encoding="utf-8"))
        candidate_payload = payload.get("best_candidate", payload)
        center_theta.update(candidate_payload["theta"])
        if terminal_name == "psi_child_change_2023":
            center_terminal_coordinate = float(candidate_payload["new_psi_child"]) - float(
                candidate_payload["old_psi_child"]
            )
        else:
            center_terminal_coordinate = float(candidate_payload["new_psi_child"])
        center_status = str(center_path)
        center_sha256 = source_sha256(center_path)
    center = transition_unit_from_candidate(center_theta, center_terminal_coordinate)

    if task == 1:
        unit = center
        design = "anchor"
    elif str(args.panel_design) == "coordinate":
        expected_size = 1 + 2 * len(TRANSITION_SEARCH_DOMAIN)
        if size != expected_size:
            raise ValueError(
                f"Coordinate panel requires --panel-size={expected_size}"
            )
        offset = task - 2
        dimension = offset // 2
        sign = -1.0 if offset % 2 == 0 else 1.0
        unit = center.copy()
        unit[dimension] = np.clip(
            unit[dimension] + sign * float(args.panel_local_radius), 0.0, 1.0
        )
        design = f"coordinate_{dimension}_{'minus' if sign < 0 else 'plus'}"
    else:
        local_count = max(1, (size - 1) // 2)
        if task - 1 <= local_count:
            lhs = latin_hypercube(
                local_count, len(TRANSITION_SEARCH_DOMAIN), int(args.panel_seed)
            )
            perturbation = 2.0 * lhs[task - 2] - 1.0
            unit = np.clip(
                center + float(args.panel_local_radius) * perturbation, 0.0, 1.0
            )
            design = "local_latin_hypercube"
        else:
            wide_count = size - 1 - local_count
            lhs = latin_hypercube(
                wide_count,
                len(TRANSITION_SEARCH_DOMAIN),
                int(args.panel_seed) + 1,
            )
            unit = lhs[task - 2 - local_count]
            design = "wide_latin_hypercube"
    theta, terminal_coordinate = candidate_from_transition_unit(source_theta, unit)
    return theta, terminal_coordinate, {
        "task_id": task,
        "panel_size": size,
        "panel_seed": int(args.panel_seed),
        "design": design,
        "panel_design": str(args.panel_design),
        "local_radius": float(args.panel_local_radius),
        "center": center_status,
        "center_sha256": center_sha256,
        "terminal_preference_coordinate": terminal_name,
        "unit_vector": unit.tolist(),
        "domain": [
            {"name": name, "lower": lower, "upper": upper, "transform": kind}
            for name, lower, upper, kind in TRANSITION_SEARCH_DOMAIN
        ],
    }


def source_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def activate_model_profile(
    name: str,
    theta: dict[str, float],
) -> tuple[tuple[tuple[str, float, float, str], ...], dict[str, Any], dict[str, Any]]:
    """Return the active domain, external overrides, and serializable metadata."""
    if name == BASELINE_MODEL_PROFILE:
        return tuple(E5F_DOMAIN), {}, {
            "name": BASELINE_MODEL_PROFILE,
            "permanent_income_levels": False,
            "first_birth_fixed_cost": False,
        }
    if name == REPAIRED_MODEL_PROFILE:
        theta.setdefault("first_birth_fixed_cost", 0.0)
        overrides = e5f_income_entry_overrides()
        return tuple(E5F_INCOME_ENTRY_DOMAIN), overrides, {
            "name": REPAIRED_MODEL_PROFILE,
            "profile_id": E5F_INCOME_ENTRY_PROFILE_NAME,
            "permanent_income_levels": "externally measured PSID distribution",
            "permanent_income_log_variance": float(
                overrides["permanent_income_log_variance"]
            ),
            "income_state_count": int(len(np.asarray(overrides["z_grid"]))),
            "first_birth_fixed_cost": "one additional free parameter",
            "first_birth_fixed_cost_semantics": (
                "one-time utility cost paid only when the first child arrives"
            ),
        }
    raise ValueError(f"Unknown model profile: {name}")


def solve_old_steady_state(
    chain: Any,
    base_overrides: dict[str, Any],
    *,
    initial_psi: float,
    completed_fertility_target: float,
    completed_fertility_tolerance: float,
    normalize: bool,
) -> tuple[Any, SimpleNamespace, np.ndarray, float, dict[str, Any]]:
    """Solve the initial steady state, optionally deriving its fertility intercept."""
    cache: dict[float, tuple[Any, SimpleNamespace, np.ndarray, float, float]] = {}

    def evaluate(psi: float) -> tuple[Any, SimpleNamespace, np.ndarray, float, float]:
        key = float(psi)
        if key not in cache:
            overrides = dict(base_overrides)
            overrides["psi_child"] = key
            solution, parameters, price, seconds = closure.solve_ge(chain, overrides)
            moments = chain.extract_moments(solution, parameters)
            cache[key] = (
                solution,
                parameters,
                price,
                float(seconds),
                float(moments["tfr"]),
            )
            print(
                "OLD_SS_NORMALIZATION_EVAL "
                f"psi={key:.8f} tfr={cache[key][4]:.8f} "
                f"seconds={float(seconds):.2f}",
                flush=True,
            )
        return cache[key]

    initial = evaluate(float(initial_psi))
    target = float(completed_fertility_target)
    tolerance = float(completed_fertility_tolerance)
    if not normalize or abs(initial[4] - target) <= tolerance:
        return (*initial[:4], {
            "status": "normalized_at_initial_guess" if normalize else "fixed_intercept",
            "psi_child": float(initial_psi),
            "completed_fertility": initial[4],
            "target": target,
            "absolute_gap": abs(initial[4] - target),
            "stationary_solves": len(cache),
            "stationary_solve_seconds": sum(row[3] for row in cache.values()),
        })

    lower, upper = -3.0, 3.0
    low = evaluate(lower)
    high = evaluate(upper)
    low_gap, high_gap = low[4] - target, high[4] - target
    while low_gap * high_gap > 0.0 and max(abs(lower), abs(upper)) < 24.0:
        if low_gap > 0.0 and high_gap > 0.0:
            lower *= 2.0
            low = evaluate(lower)
            low_gap = low[4] - target
        else:
            upper *= 2.0
            high = evaluate(upper)
            high_gap = high[4] - target
    if low_gap == 0.0:
        chosen_psi, chosen = lower, low
    elif high_gap == 0.0:
        chosen_psi, chosen = upper, high
    elif low_gap * high_gap > 0.0:
        raise RuntimeError(
            "Old-steady-state fertility normalization is not bracketed: "
            f"TFR({lower:g})={low[4]:.6g}, TFR({upper:g})={high[4]:.6g}, "
            f"target={target:.6g}"
        )
    else:
        chosen_psi, chosen = float(initial_psi), initial
        for _ in range(18):
            midpoint = 0.5 * (lower + upper)
            mid = evaluate(midpoint)
            mid_gap = mid[4] - target
            if abs(mid_gap) < abs(chosen[4] - target):
                chosen_psi, chosen = midpoint, mid
            if abs(mid_gap) <= tolerance:
                break
            if low_gap * mid_gap <= 0.0:
                upper, high_gap = midpoint, mid_gap
            else:
                lower, low_gap = midpoint, mid_gap
        if abs(chosen[4] - target) > tolerance:
            raise RuntimeError(
                "Old-steady-state fertility normalization missed tolerance: "
                f"best={chosen[4]:.8g}, target={target:.8g}, "
                f"gap={chosen[4] - target:.3e}"
            )
    return (*chosen[:4], {
        "status": "derived_intercept",
        "psi_child": float(chosen_psi),
        "completed_fertility": chosen[4],
        "target": target,
        "absolute_gap": abs(chosen[4] - target),
        "stationary_solves": len(cache),
        "stationary_solve_seconds": sum(row[3] for row in cache.values()),
    })


def first_birth_flows_by_age(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> np.ndarray:
    """Realized first-birth flows by model age in the current period."""
    model = calendar.model
    flows = np.zeros(int(P.J))
    fecundity = model.get_fecundity_by_age(P)
    settled = model.readiness_settled_state(P)
    for j in range(int(P.J)):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        for zz in range(evaluation.g_pre.shape[4]):
            childless = evaluation.g_pre[:, :, :, j, zz, 0, settled]
            attempt = evaluation.policy.fert_probs[:, :, :, j, zz, 1]
            flows[j] += float(np.sum(float(fecundity[j]) * childless * attempt))
    return flows


def first_birth_housing_response(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> float:
    """Horizon-zero birth-versus-no-birth housing response on the dated risk set."""
    if int(getattr(P, "housing_event_horizon", 0)) != 0:
        raise ValueError("The transition calibration currently requires housing_event_horizon=0")
    model = calendar.model
    policy = evaluation.policy
    fecundity = model.get_fecundity_by_age(P)
    settled = model.readiness_settled_state(P)
    shape = (
        len(policy.bp_pol),
        1 + int(P.n_house),
        int(P.I),
        evaluation.g_pre.shape[4],
        int(P.n_parity),
        int(P.n_child_states),
    )
    response_sum = 0.0
    response_mass = 0.0
    for j in range(int(P.J)):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        for zz in range(evaluation.g_pre.shape[4]):
            childless = evaluation.g_pre[:, :, :, j, zz, 0, settled]
            attempt = policy.fert_probs[:, :, :, j, zz, 1]
            realized = float(fecundity[j]) * childless * attempt
            mass = float(np.sum(realized))
            if mass <= 1e-14:
                continue
            birth_cohort = np.zeros(shape)
            birth_cohort[:, :, :, zz, 1, 1] = realized
            control_cohort = np.zeros(shape)
            control_cohort[:, :, :, zz, 0, settled] = realized
            observed_birth = model.realize_current_choices_markov_income(
                birth_cohort,
                j,
                policy.loc_probs,
                policy.tenure_choice,
                policy.tenure_probs,
                policy.maps.lmm_idx,
                policy.maps.lmm_wt,
                policy.maps.tmx_idx,
                policy.maps.tmx_wt,
                use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
            )
            observed_control = model.realize_current_choices_markov_income(
                control_cohort,
                j,
                policy.loc_probs,
                policy.tenure_choice,
                policy.tenure_probs,
                policy.maps.lmm_idx,
                policy.maps.lmm_wt,
                policy.maps.tmx_idx,
                policy.maps.tmx_wt,
                use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
            )
            birth_housing = model.mean_housing_distribution_markov(
                observed_birth, j, policy.hR_pol, P
            )
            control_housing = model.mean_housing_distribution_markov(
                observed_control, j, policy.hR_pol, P
            )
            response_sum += mass * (birth_housing - control_housing)
            response_mass += mass
    return response_sum / max(response_mass, 1e-15)


def target_age_fertility_moments(
    distribution: np.ndarray,
    P: SimpleNamespace,
    target_age: float = TARGET_AGE,
) -> dict[str, float]:
    """Completed-fertility stock for the one model cohort centered at age 42."""
    model = calendar.model
    age_index = model.age_to_index(P, target_age)
    block = np.asarray(distribution[:, :, :, age_index, :, :, :], dtype=float)
    total = float(np.sum(block))
    parity = np.asarray(
        [float(np.sum(block[:, :, :, :, n, :])) / max(total, 1e-15) for n in range(int(P.n_parity))]
    )
    weights = np.arange(int(P.n_parity), dtype=float)
    weights[-1] = float(getattr(P, "tfr_top_bin_weight", weights[-1]))
    return {
        "target_age": float(P.age_start + age_index * P.da),
        "target_age_index": int(age_index),
        "target_age_mass": total,
        "tfr": float(np.sum(weights * parity)),
        "childless_rate": float(parity[0]),
    }


def transition_cross_section_moments(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    target_names: tuple[str, ...],
) -> dict[str, float]:
    """Measure the E5 cross-sectional targets on a dated transition distribution."""
    model = calendar.model
    from intergen_eqscale_seq_optimized.calibration import extract_moments

    stats = model.compute_markov_statistics(
        evaluation.g_current,
        evaluation.policy.fert_probs,
        evaluation.policy.loc_probs,
        P,
        b_grid,
        evaluation.policy.price,
        evaluation.policy.hR_pol,
        asset_g=evaluation.g_post_fertility,
        bequest_g=evaluation.g_current,
        bp_pol=evaluation.policy.bp_pol,
    )
    stats.total_mass = float(np.sum(evaluation.g_current))
    stats.aggregate_housing_demand = float(np.sum(evaluation.demand_by_loc))
    stats.best_max_abs_rel_excess = float(evaluation.relative_market_residual)
    stats.g = evaluation.g_current
    stats.hR_pol = evaluation.policy.hR_pol
    stats.owner_user_cost = np.asarray(P.user_cost_rate * evaluation.policy.price)
    stats.housing_increment_0to1_eventstudy_t3 = first_birth_housing_response(
        evaluation, P
    )
    extracted = extract_moments(stats, P)
    extracted.update(target_age_fertility_moments(evaluation.g_current, P))
    moments = {name: float(extracted[name]) for name in target_names}
    nonfinite = [name for name, value in moments.items() if not math.isfinite(value)]
    if nonfinite:
        raise RuntimeError(f"Transition target moments are nonfinite: {nonfinite}")
    return moments


def cohort_timing_moments(
    old_stationary_flows: np.ndarray,
    period_records: dict[int, dict[str, Any]],
    P: SimpleNamespace,
    terminal_period: int = TRANSITION_PERIODS,
    target_age: float = TARGET_AGE,
) -> dict[str, Any]:
    """First-birth timing for the cohort aged 42 at the 2023 model date."""
    model = calendar.model
    target_index = model.age_to_index(P, target_age)
    initial_index = target_index - int(terminal_period)
    if initial_index < 0:
        raise ValueError("The target cohort does not exist at the initial date")
    ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    rows: list[dict[str, Any]] = []
    for age_index in range(initial_index):
        rows.append(
            {
                "source": "old_steady_state_prehistory",
                "period": age_index - initial_index,
                "age_index": age_index,
                "age": float(ages[age_index]),
                "first_birth_flow": float(old_stationary_flows[age_index]),
            }
        )
    for period in range(int(terminal_period) + 1):
        age_index = initial_index + period
        flows = np.asarray(period_records[period]["first_birth_flows_by_age"], dtype=float)
        rows.append(
            {
                "source": "simulated_transition",
                "period": period,
                "age_index": age_index,
                "age": float(ages[age_index]),
                "first_birth_flow": float(flows[age_index]),
            }
        )
    total = sum(float(row["first_birth_flow"]) for row in rows)
    mean_age = sum(
        float(row["age"]) * float(row["first_birth_flow"]) for row in rows
    ) / max(total, 1e-15)
    late_share = sum(
        float(row["first_birth_flow"]) for row in rows if float(row["age"]) >= 30.0
    ) / max(total, 1e-15)
    return {
        "target_cohort_initial_age": float(ages[initial_index]),
        "target_cohort_terminal_age": float(ages[target_index]),
        "first_birth_flow_total": total,
        "mean_age_first_birth": mean_age,
        "share_first_births_age30plus": late_share,
        "flow_rows": rows,
    }


def target_fit_rows(
    candidate: str,
    moments: dict[str, float],
    target_system: Any,
) -> tuple[list[dict[str, Any]], float]:
    rows: list[dict[str, Any]] = []
    loss = 0.0
    for name, target, weight in zip(
        target_system.moment_names,
        target_system.target_values,
        target_system.weights,
        strict=True,
    ):
        model = float(moments[name])
        gap = model - float(target)
        contribution = float(weight) * gap**2
        loss += contribution
        rows.append(
            {
                "candidate": candidate,
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": float(weight),
                "loss_contribution": contribution,
                "standardized_gap": math.copysign(math.sqrt(contribution), gap),
            }
        )
    return rows, loss


def parameter_rows(theta: dict[str, float], candidates: list[float]) -> list[dict[str, Any]]:
    bounds = {name: (lower, upper) for name, lower, upper, _ in E5F_DOMAIN}
    rows: list[dict[str, Any]] = []
    for name, value in theta.items():
        domain_name = "beta_annual" if name == "beta" else name
        lower, upper = bounds.get(domain_name, (math.nan, math.nan))
        estimate = float(value) ** 0.25 if name == "beta" else float(value)
        rows.append(
            {
                "parameter": domain_name,
                "value": estimate,
                "lower_bound": lower,
                "upper_bound": upper,
                "status": (
                    "replaced_by_dated_preference_path"
                    if name == "psi_child"
                    else "candidate_value_in_bounded_transition_pilot"
                ),
                "near_bound": bool(
                    math.isfinite(lower)
                    and math.isfinite(upper)
                    and min(abs(estimate - lower), abs(upper - estimate))
                    <= 0.02 * max(upper - lower, 1e-12)
                ),
            }
        )
    rows.extend(
        [
            {
                "parameter": "psi_child_2007",
                "value": math.nan,
                "lower_bound": -3.0,
                "upper_bound": 3.0,
                "status": "externally_normalized_to_old_completed_fertility",
                "near_bound": False,
            },
            {
                "parameter": "psi_child_2023",
                "value": math.nan,
                "lower_bound": min(candidates),
                "upper_bound": max(candidates),
                "status": "one_parameter_bounded_pilot_grid",
                "near_bound": False,
            },
        ]
    )
    return rows


def make_diagnostic_plot(
    summaries: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = sorted(summaries, key=lambda row: float(row["new_psi_child"]))
    x = np.asarray([float(row["new_psi_child"]) for row in ordered])
    fertility = np.asarray([float(row["terminal_tfr"]) for row in ordered])
    loss = np.asarray([float(row["transition_loss"]) for row in ordered])
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), constrained_layout=True)
    axes[0].plot(x, fertility, marker="o")
    axes[0].axhline(1.918, color="black", ls="--", lw=1.0, label="CPS target")
    axes[0].set(
        xlabel=r"Terminal child-preference intercept $\psi_{2023}$",
        ylabel="Completed fertility, age-42 cohort",
        title="Transition fertility target",
    )
    axes[0].legend(frameon=False)
    axes[1].plot(x, loss, marker="o")
    axes[1].set(
        xlabel=r"Terminal child-preference intercept $\psi_{2023}$",
        ylabel="12-moment transition loss",
        title="Fixed-parameter pilot objective",
        yscale="log",
    )
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(outdir / "transition_calibration_pilot.png", dpi=220)
    fig.savefig(outdir / "transition_calibration_pilot.pdf")
    plt.close(fig)

    moments = list(dict.fromkeys(row["moment"] for row in fit_rows))
    labels = [str(row["candidate"]) for row in summaries]
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    offsets = np.linspace(-0.30, 0.30, max(len(labels), 1))
    positions = np.arange(len(moments), dtype=float)
    for offset, label in zip(offsets, labels, strict=True):
        rows = [row for row in fit_rows if row["candidate"] == label]
        lookup = {row["moment"]: float(row["standardized_gap"]) for row in rows}
        ax.scatter(
            positions + offset,
            [lookup[name] for name in moments],
            s=24,
            label=label,
        )
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xticks(positions)
    ax.set_xticklabels(moments, rotation=48, ha="right")
    ax.set_ylabel("Signed square root of loss contribution")
    ax.set_title("2023 target fit measured on the transition distribution")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "transition_target_fit.png", dpi=220)
    fig.savefig(outdir / "transition_target_fit.pdf")
    plt.close(fig)


def main() -> None:
    global E5F_DOMAIN, TRANSITION_SEARCH_DOMAIN
    args = parse_args()
    candidates = parse_candidates(args.candidate_new_psis)
    candidate_psi_changes = (
        None
        if args.candidate_psi_changes is None
        else parse_candidates(args.candidate_psi_changes)
    )
    if int(args.nb) != 120:
        raise ValueError("The transition calibration pilot is hard-gated to Nb=120")
    if int(args.post_2023_periods) < 0:
        raise ValueError("--post-2023-periods must be nonnegative")
    if not 0.0 < float(args.outside_origin_entry_share) < 1.0:
        raise ValueError("--outside-origin-entry-share must lie in (0,1)")

    started = time.perf_counter()
    source = args.source.resolve()
    actual_source_sha256 = source_sha256(source)
    if (
        args.expected_source_sha256 is not None
        and actual_source_sha256 != str(args.expected_source_sha256)
    ):
        raise RuntimeError(
            "Source fingerprint gate failed: "
            f"actual={actual_source_sha256}, expected={args.expected_source_sha256}"
        )
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    chain, model = transition.configure_sequential_model()
    winner, source_metadata = closure.load_winner(source)
    theta = closure.theta_from_winner(winner)
    active_domain, profile_overrides, model_profile = activate_model_profile(
        str(args.model_profile), theta
    )
    if args.fixed_first_birth_cost is not None:
        if str(args.model_profile) != REPAIRED_MODEL_PROFILE:
            raise ValueError("--fixed-first-birth-cost requires the repaired profile")
        if args.panel_task_id is not None:
            raise ValueError("A fixed first-birth cost cannot be combined with a joint panel")
        theta["first_birth_fixed_cost"] = float(args.fixed_first_birth_cost)
    E5F_DOMAIN = active_domain
    terminal_preference_row = (
        ("psi_child_change_2023", -1.50, 0.20, "asinh")
        if str(args.model_profile) == REPAIRED_MODEL_PROFILE
        else ("psi_child_2023", -1.25, 0.20, "asinh")
    )
    TRANSITION_SEARCH_DOMAIN = tuple(
        row for row in active_domain if row[0] != "psi_child"
    ) + (terminal_preference_row,)
    panel_metadata = None
    panel_terminal_preference_coordinate = math.nan
    if args.panel_task_id is not None:
        theta, panel_terminal_preference_coordinate, panel_metadata = panel_candidate(
            theta, args
        )
    target_system = e5_target_system()
    if (
        args.expected_target_set is not None
        and target_system.name != str(args.expected_target_set)
    ):
        raise RuntimeError(
            "Target-set gate failed: "
            f"actual={target_system.name}, expected={args.expected_target_set}"
        )
    if (
        args.expected_target_fingerprint is not None
        and target_system.fingerprint != str(args.expected_target_fingerprint)
    ):
        raise RuntimeError(
            "Target fingerprint gate failed: "
            f"actual={target_system.fingerprint}, "
            f"expected={args.expected_target_fingerprint}"
        )
    target_system.require_identified(len(E5F_DOMAIN))
    base = closure.make_overrides(chain, theta, nb=int(args.nb), profile="e5f-floor")
    base.update(profile_overrides)
    base.update(theta)

    print("TRANSITION_CALIBRATION_OLD_STEADY_STATE", flush=True)
    (
        old_solution,
        old_parameters,
        old_price,
        old_seconds,
        old_fertility_normalization,
    ) = solve_old_steady_state(
        chain,
        base,
        initial_psi=float(args.old_psi_child),
        completed_fertility_target=float(args.old_completed_fertility_target),
        completed_fertility_tolerance=float(args.old_completed_fertility_tol),
        normalize=(
            panel_metadata is not None
            or str(args.model_profile) == REPAIRED_MODEL_PROFILE
        ),
    )
    old_psi_child = float(old_fertility_normalization["psi_child"])
    if panel_metadata is not None:
        if str(args.model_profile) == REPAIRED_MODEL_PROFILE:
            new_psi_child = old_psi_child + float(panel_terminal_preference_coordinate)
            panel_metadata["psi_child_change_2023"] = float(
                panel_terminal_preference_coordinate
            )
        else:
            new_psi_child = float(panel_terminal_preference_coordinate)
        panel_metadata["old_psi_child"] = old_psi_child
        panel_metadata["new_psi_child"] = new_psi_child
        candidates = [new_psi_child]
    elif candidate_psi_changes is not None:
        candidates = [old_psi_child + change for change in candidate_psi_changes]
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError("Production dimension gate failed")
    income_profile_gates: dict[str, Any] = {
        "permanent_income_levels_enabled": bool(
            getattr(old_parameters, "permanent_income_levels_enabled", False)
        ),
        "income_state_count": int(old_parameters.Nz),
        "stationary_weight_max_abs_gap": float(
            np.max(
                np.abs(
                    np.asarray(old_parameters.z_weights, dtype=float)
                    @ np.asarray(old_parameters.Pi_z, dtype=float)
                    - np.asarray(old_parameters.z_weights, dtype=float)
                )
            )
        ),
    }
    if str(args.model_profile) == REPAIRED_MODEL_PROFILE:
        if (
            not income_profile_gates["permanent_income_levels_enabled"]
            or income_profile_gates["income_state_count"] != 15
            or income_profile_gates["stationary_weight_max_abs_gap"] > 1.0e-12
        ):
            raise RuntimeError(f"Measured-income profile gates failed: {income_profile_gates}")
    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    old_shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price, old_parameters, b_grid, old_shared
    )
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    initial_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution, old_policy, old_parameters, b_grid, old_shared
    )
    operator_gates = transition.operator_gates(
        old_solution,
        old_policy,
        initial_g_pre,
        old_parameters,
        b_grid,
        old_shared,
    )
    operator_gates.update(reconstruction)
    if max(
        float(operator_gates["stationary_post_fertility_nesting_l1"]),
        float(operator_gates["one_step_constant_path_nesting_l1"]),
        float(operator_gates["mature_flow_abs_error"]),
        float(operator_gates["birth_flow_abs_error"]),
        float(operator_gates["topcode_adjusted_birth_flow_abs_error"]),
    ) > 5e-9:
        raise RuntimeError(f"Stationary operator gates failed: {operator_gates}")

    old_evaluation = calendar.evaluate_period(
        old_price,
        initial_g_pre,
        old_parameters,
        b_grid,
        old_shared,
        calendar.SolveCounter(),
        supplied_policy=old_policy,
    )
    old_transition_moments = transition_cross_section_moments(
        old_evaluation,
        old_parameters,
        b_grid,
        target_system.moment_names,
    )
    old_stationary_moments = chain.extract_moments(old_solution, old_parameters)
    stationary_measurement_gaps = {
        name: float(old_transition_moments[name]) - float(old_stationary_moments[name])
        for name in target_system.moment_names
    }
    max_stationary_measurement_gap = max(
        abs(value) for value in stationary_measurement_gaps.values()
    )
    if max_stationary_measurement_gap > 2e-8:
        raise RuntimeError(
            "Transition measurement does not nest stationary measurement: "
            f"max gap={max_stationary_measurement_gap:.3e}, "
            f"gaps={stationary_measurement_gaps}"
        )

    old_renewal = closure.topcode_consistent_renewal_accounting(
        old_solution, old_parameters
    )
    E_old = float(old_solution.entry_rate)
    B_old_raw = float(old_solution.entrants_mature_total)
    B_old = float(old_renewal["topcode_adjusted_mature_entrant_households"])
    births_old_raw = float(old_solution.total_births_kfe)
    births_old = float(old_renewal["topcode_adjusted_birth_children"])
    outside_share = float(args.outside_origin_entry_share)
    retention = (1.0 - outside_share) * E_old / B_old
    if not 0.0 <= retention <= 1.0:
        raise RuntimeError(f"Outside share implies invalid retention: {retention}")
    conversion = float(old_parameters.entrant_conversion_factor)
    stationary_initial_g_pre = initial_g_pre
    model_ages = (
        float(old_parameters.age_start)
        + np.arange(int(old_parameters.J), dtype=float) * float(old_parameters.da)
    )
    stationary_age_mass = np.sum(
        stationary_initial_g_pre, axis=(0, 1, 2, 4, 5, 6)
    )
    age_reweight = transition.acs_2007_age_reweight_diagnostic(
        stationary_age_mass,
        model_ages,
        E_old,
        periods=TRANSITION_PERIODS,
        period_years=float(old_parameters.period_years),
    )
    initial_g_pre = transition.reweight_distribution_to_acs_2007_ages(
        stationary_initial_g_pre,
        age_reweight,
    )
    supply_rule, supply_normalization = calendar.normalize_date0_housing_supply(
        initial_g_pre,
        old_policy,
        old_parameters,
        b_grid,
        old_shared,
        "static-elastic",
    )
    old_first_birth_flows = first_birth_flows_by_age(old_evaluation, old_parameters)
    population_target_indices = transition.census_household_target_indices()
    population_age_target_years = transition.census_household_age_target_years()
    population_target = float(population_target_indices[TRANSITION_PERIODS])

    all_fit_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    target_measurements = {
        "tfr": "completed fertility of the model cohort centered at age 42 in 2023",
        "childless_rate": "childless share of the model cohort centered at age 42 in 2023",
        "mean_age_first_birth": "cohort first-birth flows: old-steady-state prehistory plus dated 2007--2023 transition",
        "share_first_births_age30plus": "same dated target-cohort first-birth ledger",
        "housing_increment_0to1": "horizon-zero birth-versus-no-birth housing response on the 2023 risk set",
        "remaining_targets": "2023 transition cross-section with beginning-of-period wealth timing",
    }
    for index, new_psi in enumerate(candidates, start=1):
        label = (
            f"task_{int(args.panel_task_id):03d}"
            if args.panel_task_id is not None
            else f"psi_{new_psi:+.6f}".replace("+", "p").replace("-", "m")
        )
        candidate_dir = outdir / "cases" / label
        candidate_dir.mkdir(parents=True, exist_ok=True)
        period_records: dict[int, dict[str, Any]] = {}

        def observer(
            period: int,
            evaluation: calendar.PeriodEvaluation,
            parameters: SimpleNamespace,
            grid: np.ndarray,
        ) -> None:
            period_records[int(period)] = {
                "moments": transition_cross_section_moments(
                    evaluation, parameters, grid, target_system.moment_names
                ),
                "first_birth_flows_by_age": first_birth_flows_by_age(
                    evaluation, parameters
                ),
            }

        print(
            f"TRANSITION_CALIBRATION_CASE {index}/{len(candidates)} "
            f"new_psi={new_psi:.6f}",
            flush=True,
        )
        counter = calendar.SolveCounter()
        case_started = time.perf_counter()
        total_periods = TRANSITION_PERIODS + 1 + int(args.post_2023_periods)
        paths, ages, children, scenario_summary = transition.run_birth_vintage_scenario(
            label=label,
            initial_g_pre=initial_g_pre,
            baseline_price=float(np.asarray(old_price).reshape(-1)[0]),
            baseline_B=B_old,
            baseline_births=births_old,
            baseline_raw_B=B_old_raw,
            baseline_raw_births=births_old_raw,
            parameters=copy.deepcopy(old_parameters),
            b_grid=b_grid,
            periods=total_periods,
            outside_flow=outside_share * E_old,
            retention=retention,
            conversion=conversion,
            delay_periods=4,
            initial_birth_pipeline_multiplier=1.0,
            population_target_indices=population_target_indices,
            population_age_target_years=population_age_target_years,
            policy_case=str(args.policy_case),
            policy_start_period=TRANSITION_PERIODS,
            old_psi_child=old_psi_child,
            new_psi_child=float(new_psi),
            preference_transition_periods=TRANSITION_PERIODS,
            supply_rule=supply_rule,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            outdir=candidate_dir,
            counter=counter,
            period_observer=observer,
        )
        calendar.write_csv(candidate_dir / "transition_path.csv", paths)
        calendar.write_csv(candidate_dir / "adult_age_distribution.csv", ages)
        calendar.write_csv(candidate_dir / "child_pipeline.csv", children)
        if sorted(period_records) != list(range(total_periods)):
            raise RuntimeError(f"Incomplete transition moment ledger: {sorted(period_records)}")
        timing = cohort_timing_moments(
            old_first_birth_flows,
            period_records,
            old_parameters,
        )
        terminal_moments = dict(period_records[TRANSITION_PERIODS]["moments"])
        terminal_moments["mean_age_first_birth"] = float(timing["mean_age_first_birth"])
        terminal_moments["share_first_births_age30plus"] = float(
            timing["share_first_births_age30plus"]
        )
        fit_rows, loss = target_fit_rows(label, terminal_moments, target_system)
        all_fit_rows.extend(fit_rows)
        calendar.write_csv(candidate_dir / "target_fit.csv", fit_rows)
        write_json(candidate_dir / "cohort_timing_ledger.json", timing)
        terminal = paths[TRANSITION_PERIODS]
        endpoint = paths[-1]
        summary = {
            "candidate": label,
            "theta": theta,
            "new_psi_child": float(new_psi),
            "transition_loss": float(loss),
            "terminal_tfr": float(terminal_moments["tfr"]),
            "terminal_childless_rate": float(terminal_moments["childless_rate"]),
            "terminal_mean_age_first_birth": float(terminal_moments["mean_age_first_birth"]),
            "terminal_share_first_births_age30plus": float(
                terminal_moments["share_first_births_age30plus"]
            ),
            "terminal_population_index": float(terminal["population_index"]),
            "census_household_population_target_index": population_target,
            "population_target_gap": float(terminal["population_index"]) - population_target,
            "terminal_housing_cost_index": float(terminal["asset_price_index"]),
            "policy_case": str(args.policy_case),
            "post_2023_periods": int(args.post_2023_periods),
            "endpoint_year": START_YEAR + int(endpoint["years_from_start"]),
            "endpoint_population_index": float(endpoint["population_index"]),
            "endpoint_housing_cost_index": float(endpoint["asset_price_index"]),
            "endpoint_owner_rate": float(endpoint["owner_rate"]),
            "endpoint_dependent_child_owner_rate": float(
                endpoint["dependent_child_owner_rate"]
            ),
            "endpoint_topcode_adjusted_births": float(
                endpoint["birth_children_topcode_adjusted"]
            ),
            "model_solve_count": int(counter.total),
            "elapsed_seconds": time.perf_counter() - case_started,
            "max_market_residual": float(scenario_summary["max_market_residual"]),
            "max_mass_residual": float(scenario_summary["max_abs_mass_accounting_residual"]),
        }
        summaries.append(summary)
        calendar.write_csv(outdir / "candidate_summary.csv", summaries)
        calendar.write_csv(outdir / "target_fit_long.csv", all_fit_rows)
        write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "running",
                "completed_cases": index,
                "total_cases": len(candidates),
                "latest": summary,
            },
        )
        best = min(summaries, key=lambda row: float(row["transition_loss"]))
        write_json(outdir / "best_so_far.json", best)

    best = min(summaries, key=lambda row: float(row["transition_loss"]))
    params = parameter_rows(theta, candidates)
    for row in params:
        if row["parameter"] == "psi_child_2007":
            row["value"] = old_psi_child
        elif row["parameter"] == "psi_child_2023":
            row["value"] = float(best["new_psi_child"])
    calendar.write_csv(outdir / "parameter_table.csv", params)
    make_diagnostic_plot(summaries, all_fit_rows, outdir)
    summary = {
        "status": (
            "complete_transition_calibration_panel_task"
            if panel_metadata is not None
            else "complete_bounded_fixed_parameter_pilot"
        ),
        "source": str(source),
        "source_sha256": actual_source_sha256,
        "source_metadata": source_metadata,
        "model_profile": model_profile,
        "income_profile_gates": income_profile_gates,
        "target_set": target_system.name,
        "target_fingerprint": target_system.fingerprint,
        "target_count": target_system.count,
        "stationary_free_parameter_count": len(E5F_DOMAIN),
        "transition_free_parameter_count": len(TRANSITION_SEARCH_DOMAIN),
        "identification_status": (
            f"{len(TRANSITION_SEARCH_DOMAIN)} transition parameters against "
            f"{target_system.count} dated target moments"
            if panel_metadata is not None
            else "fixed-parameter transition frontier; not a joint transition calibration"
        ),
        "panel_design": panel_metadata,
        "fixed_first_birth_cost": args.fixed_first_birth_cost,
        "candidate_psi_changes": candidate_psi_changes,
        "policy_case": str(args.policy_case),
        "post_2023_periods": int(args.post_2023_periods),
        "old_psi_child": old_psi_child,
        "old_fertility_normalization": old_fertility_normalization,
        "old_model_completed_fertility": float(old_transition_moments["tfr"]),
        "old_completed_fertility_reference": float(
            args.old_completed_fertility_target
        ),
        "outside_origin_entry_share": outside_share,
        "outside_origin_entry_status": "provisional_across_market_anchor_not_national_migration",
        "renewal_retention": retention,
        "housing_supply": supply_normalization,
        "target_measurements": target_measurements,
        "population_bridge": {
            "status": "externally_normalized_not_estimated",
            "total_households": "Census HH-3, ages 18--85 after ACS coverage adjustment",
            "age_distribution": "national ACS household-head age shares at each four-year date",
            "initial_age_reweight": age_reweight,
            "target_indices": population_target_indices,
        },
        "population_validation_status": (
            "matched by the declared Census/ACS household-formation and migration bridge; "
            "not generated by post-2007 fertility"
        ),
        "census_household_population_target_index": population_target,
        "stationary_operator_gates": operator_gates,
        "stationary_measurement_max_abs_gap": max_stationary_measurement_gap,
        "stationary_measurement_gaps": stationary_measurement_gaps,
        "candidates": summaries,
        "best_candidate": best,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json(outdir / "summary.json", summary)
    (outdir / "README.md").write_text(
        "\n".join(
            [
                "# E5F transition-calibration pilot",
                "",
                "This packet measures the existing twelve-row E5 target system on the",
                "simulated 2023 distribution. Date 0 preserves conditional states from",
                "an old steady-state benchmark but reweights age masses to the observed",
                "2007 householder-age profile. In panel mode the existing eight structural parameters",
                "and the terminal child-preference intercept form a nine-dimensional",
                "transition-calibration candidate. This remains a search packet, not a",
                "promoted estimate.",
                "",
                "The fertility stock targets use the model cohort centered at age 42.",
                "First-birth timing combines old-steady-state prehistory with realized",
                "first-birth flows during 2007--2023. The previous period-TFR",
                "normalization is therefore not treated as the completed-fertility",
                "calibration target.",
                "",
                "The 2007--2023 household-count and age paths are matched by an explicit",
                "Census/ACS bridge. This is an externally normalized migration and",
                "household-formation input: the code does not attribute that growth to",
                "post-2007 births, which cannot mature before the terminal date.",
                "",
                "Primary files: `summary.json`, `candidate_summary.csv`,",
                "`target_fit_long.csv`, `parameter_table.csv`, and the two diagnostic",
                "figures. Each case folder contains the full transition and cohort timing",
                "ledger.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(
        "TRANSITION_CALIBRATION_COMPLETE "
        f"best={best['candidate']} loss={best['transition_loss']:.6f} "
        f"tfr={best['terminal_tfr']:.6f}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        try:
            failed_args = parse_args()
            write_json(
                failed_args.outdir.resolve() / "failure.json",
                {
                    "status": "invalid_candidate",
                    "model_profile": str(failed_args.model_profile),
                    "panel_task_id": failed_args.panel_task_id,
                    "error_type": type(error).__name__,
                    "error": str(error),
                    "traceback": traceback.format_exc(),
                },
            )
        except Exception:
            pass
        raise
