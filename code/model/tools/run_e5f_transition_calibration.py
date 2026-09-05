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

from intergen_eqscale_seq_optimized.e5_profile import (
    E5_HOUSING_EVENT_HORIZON,
)
from intergen_eqscale_seq_optimized.e5_young_ownership_profile import (
    E5_BASELINE_TARGET_PROFILE,
    E5_YOUNG_OWNERSHIP_PROFILE,
    E5_YOUNG_OWNERSHIP_PROVENANCE,
    e5_target_system_for_profile,
)
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
    parser.add_argument(
        "--replacement-fertility",
        type=float,
        default=transition.REPLACEMENT_FERTILITY,
        help=(
            "Births per woman needed for replacement; its inverse maps "
            "top-code-adjusted births into future entrant households."
        ),
    )
    parser.add_argument(
        "--old-completed-fertility-target",
        type=float,
        default=None,
        help="Old-steady-state completed fertility; defaults to replacement fertility.",
    )
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
    parser.add_argument(
        "--fixed-first-child-room-jump",
        type=float,
        default=None,
        help=(
            "Diagnostic nonnegative extra housing-service floor applied whenever "
            "at least one child is at home. It is external to the active "
            "ten-parameter calibration and nests the current model at zero."
        ),
    )
    parser.add_argument(
        "--estimate-first-child-room-jump",
        action="store_true",
        help=(
            "Add the nonnegative first-child housing-service intercept to the "
            "joint transition-search domain. This adds one free parameter and "
            "requires the repaired income-entry profile."
        ),
    )
    parser.add_argument("--outside-origin-entry-share", type=float, default=0.169)
    parser.add_argument(
        "--housing-supply-elasticity",
        type=float,
        default=None,
        help=(
            "Externally fixed dated-transition housing-supply elasticity. "
            "Omission retains the calibrated stationary value exactly."
        ),
    )
    parser.add_argument(
        "--fixed-tenure-choice-kappa",
        type=float,
        default=None,
        help=(
            "Externally fixed nonnegative tenure-choice dispersion. It is not "
            "added to the free-parameter count or the dated target system."
        ),
    )
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
    parser.add_argument(
        "--target-profile",
        choices=(E5_BASELINE_TARGET_PROFILE, E5_YOUNG_OWNERSHIP_PROFILE),
        default=E5_BASELINE_TARGET_PROFILE,
        help=(
            "Explicit dated-target contract. The default is byte-for-byte the "
            "live E5 target system; the young-ownership profile appends one "
            "diagnostic overidentifying row."
        ),
    )
    parser.add_argument("--expected-target-set", default=None)
    parser.add_argument("--expected-target-fingerprint", default=None)
    parser.add_argument("--expected-code-bundle-sha256", default=None)
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
            old_psi = candidate_payload.get(
                "old_psi_child", payload.get("old_psi_child")
            )
            new_psi = candidate_payload.get(
                "new_psi_child", payload.get("new_psi_child")
            )
            if old_psi is None or new_psi is None:
                raise ValueError(
                    "A repaired-profile center must save old_psi_child and "
                    "new_psi_child either in best_candidate or at the top level"
                )
            center_terminal_coordinate = float(new_psi) - float(old_psi)
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


def code_fingerprint_contract(model: Any) -> dict[str, Any]:
    """Pin every active file that can change the transition experiment."""
    profile_dir = Path(model.__file__).resolve().parent
    paths = {
        "transition_calibration_driver": Path(__file__).resolve(),
        "calendar_operator": Path(calendar.__file__).resolve(),
        "e5f_transition_operator": Path(transition.__file__).resolve(),
        "closure_audit": Path(closure.__file__).resolve(),
    }
    paths.update(
        {
            f"sequential_package/{path.name}": path
            for path in sorted(profile_dir.glob("*.py"))
        }
    )
    files = {name: source_sha256(path) for name, path in paths.items()}
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
    return {
        "files": files,
        "bundle_sha256": hashlib.sha256(encoded).hexdigest(),
    }


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


def configure_first_child_room_jump(
    *,
    model_profile_name: str,
    theta: dict[str, float],
    active_domain: tuple[tuple[str, float, float, str], ...],
    profile_overrides: dict[str, Any],
    model_profile: dict[str, Any],
    fixed_jump: float | None,
    estimate_jump: bool,
) -> tuple[
    tuple[tuple[str, float, float, str], ...],
    dict[str, Any],
    dict[str, Any],
]:
    """Apply the nested first-child housing-intercept contract."""
    if estimate_jump and fixed_jump is not None:
        raise ValueError("The first-child room jump cannot be both fixed and estimated")
    if not estimate_jump and fixed_jump is None:
        return active_domain, profile_overrides, model_profile
    if model_profile_name != REPAIRED_MODEL_PROFILE:
        raise ValueError("The first-child room jump requires the repaired profile")

    if estimate_jump:
        theta.setdefault("hbar_first_child_jump", 0.0)
        profile_overrides["hbar_first_child_jump"] = 0.0
        active_domain = tuple(active_domain) + (
            ("hbar_first_child_jump", 0.0, 0.5, "softzero"),
        )
        model_profile["first_child_room_jump"] = "one additional free parameter"
        model_profile["first_child_room_jump_semantics"] = (
            "extra housing-service floor whenever at least one child is at home"
        )
        model_profile["first_child_room_jump_status"] = (
            "jointly estimated transition parameter"
        )
        return active_domain, profile_overrides, model_profile

    jump = float(fixed_jump)
    if not math.isfinite(jump) or jump < 0.0:
        raise ValueError(
            "--fixed-first-child-room-jump must be finite and nonnegative"
        )
    theta["hbar_first_child_jump"] = jump
    profile_overrides["hbar_first_child_jump"] = jump
    model_profile["first_child_room_jump"] = jump
    model_profile["first_child_room_jump_status"] = (
        "externally fixed diagnostic; not included in the free-parameter count"
    )
    return active_domain, profile_overrides, model_profile


def validated_structural_stationary_age_mass(
    observed_age_mass: np.ndarray,
    *,
    entry_flow: float,
    structural_survival: np.ndarray,
    relative_tolerance: float = 1.0e-8,
) -> tuple[np.ndarray, dict[str, float]]:
    """Remove only KFE roundoff from a stationary cohort-age profile.

    A stationary cohort's age mass is the entrant flow multiplied by cumulative
    structural survival. Positive tenure smoothing can leave solver-scale
    roundoff of order 1e-9 in otherwise unit survival cells. The ACS bridge must
    not reinterpret that noise as survival above one, so this helper checks the
    full observed profile against the exact structural recursion before using
    the recursion. Any economically meaningful mismatch still fails closed.
    """
    observed = np.asarray(observed_age_mass, dtype=float).reshape(-1)
    survival = np.asarray(structural_survival, dtype=float).reshape(-1)
    entry = float(entry_flow)
    tolerance = float(relative_tolerance)
    if observed.size < 2 or survival.size != observed.size - 1:
        raise ValueError("Age mass and structural survival dimensions do not conform")
    if (
        not math.isfinite(entry)
        or entry <= 0.0
        or not math.isfinite(tolerance)
        or tolerance <= 0.0
        or np.any(~np.isfinite(observed))
        or np.any(observed <= 0.0)
        or np.any(~np.isfinite(survival))
        or np.any(survival < 0.0)
        or np.any(survival > 1.0)
    ):
        raise ValueError("Age-mass structural recursion inputs are invalid")
    structural = np.empty_like(observed)
    structural[0] = entry
    structural[1:] = entry * np.cumprod(survival)
    scale = max(entry, 1.0e-15)
    relative_gap = float(np.max(np.abs(observed - structural)) / scale)
    if relative_gap > tolerance:
        raise RuntimeError(
            "Stationary age mass disagrees with structural survival: "
            f"relative_gap={relative_gap:.3e}, tolerance={tolerance:.3e}"
        )
    return structural, {
        "status": "structural_survival_recursion_after_kfe_roundoff_gate",
        "max_abs_mass_gap": float(np.max(np.abs(observed - structural))),
        "max_relative_mass_gap_to_entry": relative_gap,
        "relative_tolerance": tolerance,
    }


def normalize_distribution_mass_roundoff(
    distribution: np.ndarray,
    *,
    expected_mass: float,
    stage: str,
    relative_tolerance: float = 5.0e-9,
) -> tuple[np.ndarray, dict[str, float | str]]:
    """Enforce exact mass for a pure redistribution after a fail-closed gate."""
    values = np.asarray(distribution, dtype=float)
    expected = float(expected_mass)
    actual = float(np.sum(values))
    tolerance = float(relative_tolerance)
    scale = max(abs(expected), 1.0e-15)
    relative_gap = abs(actual - expected) / scale
    if (
        not math.isfinite(expected)
        or expected <= 0.0
        or not math.isfinite(actual)
        or actual <= 0.0
        or not math.isfinite(tolerance)
        or tolerance <= 0.0
        or relative_gap > tolerance
    ):
        raise RuntimeError(
            f"{stage} mass gate failed: actual={actual:.17g}, "
            f"expected={expected:.17g}, relative_gap={relative_gap:.3e}, "
            f"tolerance={tolerance:.3e}"
        )
    normalized = values * (expected / actual)
    return normalized, {
        "stage": stage,
        "actual_mass_before_normalization": actual,
        "expected_mass": expected,
        "relative_gap": relative_gap,
        "relative_tolerance": tolerance,
    }


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

    initial_gap = initial[4] - target
    step = 0.25 if initial_gap < 0.0 else -0.25
    probe_psi = float(initial_psi) + step
    probe = evaluate(probe_psi)
    probe_gap = probe[4] - target
    while initial_gap * probe_gap > 0.0 and abs(step) < 24.0:
        step *= 2.0
        probe_psi = float(initial_psi) + step
        probe = evaluate(probe_psi)
        probe_gap = probe[4] - target
    if initial_gap * probe_gap > 0.0:
        raise RuntimeError(
            "Old-steady-state fertility normalization is not bracketed: "
            f"TFR({float(initial_psi):g})={initial[4]:.6g}, "
            f"TFR({probe_psi:g})={probe[4]:.6g}, "
            f"target={target:.6g}"
        )
    if float(initial_psi) < probe_psi:
        lower, low, low_gap = float(initial_psi), initial, initial_gap
        upper, high, high_gap = probe_psi, probe, probe_gap
    else:
        lower, low, low_gap = probe_psi, probe, probe_gap
        upper, high, high_gap = float(initial_psi), initial, initial_gap
    chosen_psi, chosen = min(
        ((lower, low), (upper, high)), key=lambda item: abs(item[1][4] - target)
    )
    for _ in range(14):
        width = upper - lower
        denominator = high_gap - low_gap
        if denominator == 0.0:
            candidate_psi = 0.5 * (lower + upper)
        else:
            candidate_psi = lower - low_gap * width / denominator
            candidate_psi = float(
                np.clip(candidate_psi, lower + 0.05 * width, upper - 0.05 * width)
            )
        candidate = evaluate(candidate_psi)
        candidate_gap = candidate[4] - target
        if abs(candidate_gap) < abs(chosen[4] - target):
            chosen_psi, chosen = candidate_psi, candidate
        if abs(candidate_gap) <= tolerance:
            break
        if low_gap * candidate_gap <= 0.0:
            upper, high, high_gap = candidate_psi, candidate, candidate_gap
        else:
            lower, low, low_gap = candidate_psi, candidate, candidate_gap
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


def first_birth_accounting_by_age(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> dict[str, np.ndarray]:
    """Return first-birth risk sets, realized flows, and hazards by age."""
    model = calendar.model
    flows = np.zeros(int(P.J))
    at_risk = np.zeros(int(P.J))
    fecundity = model.get_fecundity_by_age(P)
    settled = model.readiness_settled_state(P)
    for j in range(int(P.J)):
        childless_by_income = evaluation.g_pre[:, :, :, j, :, 0, settled]
        at_risk[j] = float(np.sum(childless_by_income))
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        for zz in range(evaluation.g_pre.shape[4]):
            childless = evaluation.g_pre[:, :, :, j, zz, 0, settled]
            attempt = evaluation.policy.fert_probs[:, :, :, j, zz, 1]
            flows[j] += float(np.sum(float(fecundity[j]) * childless * attempt))
    hazards = np.divide(
        flows,
        at_risk,
        out=np.zeros_like(flows),
        where=at_risk > 1e-15,
    )
    if np.any(hazards < -1e-14) or np.any(hazards > 1.0 + 1e-12):
        raise RuntimeError(f"Invalid first-birth hazards: {hazards}")
    return {
        "at_risk": at_risk,
        "flow": flows,
        "hazard": np.clip(hazards, 0.0, 1.0),
    }


def first_birth_flows_by_age(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> np.ndarray:
    """Backward-compatible realized first-birth flows by model age."""
    return first_birth_accounting_by_age(evaluation, P)["flow"]


def period_fertility_diagnostics(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> dict[str, Any]:
    """Measure current-rate fertility from dated age-specific birth flows.

    Each model age cell spans one four-year period. Dividing its four-year
    birth flow by the corresponding adult-household mass and summing over
    ages therefore gives the model analogue of period TFR; no extra factor of
    four is required. This is a diagnostic because model households, not all
    U.S. women, form the denominator.
    """
    model = calendar.model
    n_progressions = int(P.n_parity) - 1
    if n_progressions != 3:
        raise RuntimeError("Period-fertility diagnostics require parity states 0/1/2/3+")
    flows = np.zeros((int(P.J), n_progressions), dtype=float)
    age_mass = np.sum(evaluation.g_pre, axis=(0, 1, 2, 4, 5, 6))
    fecundity = model.get_fecundity_by_age(P)
    continuation = calendar.policy_continuation_birth_probs(evaluation.policy, P)
    if continuation is None:
        raise RuntimeError("Sequential continuation-birth probabilities are unavailable")
    settled = model.readiness_settled_state(P)
    for j in range(int(P.J)):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        pi_j = float(fecundity[j])
        for zz in range(evaluation.g_pre.shape[4]):
            childless = evaluation.g_pre[:, :, :, j, zz, 0, settled]
            first_attempt = evaluation.policy.fert_probs[:, :, :, j, zz, 1]
            flows[j, 0] += float(np.sum(pi_j * childless * first_attempt))
            for parity in range(1, int(P.n_parity) - 1):
                for child_state in range(parity + 1):
                    pool = evaluation.g_pre[
                        :, :, :, j, zz, parity, child_state
                    ]
                    attempt = continuation[
                        :, :, :, j, zz, 1, parity - 1, child_state
                    ]
                    flows[j, parity] += float(np.sum(pi_j * pool * attempt))
    explicit_by_age = np.sum(flows, axis=1)
    top_weight = float(P.tfr_top_bin_weight)
    topcode_by_age = (
        flows[:, 0] + flows[:, 1] + (top_weight - 2.0) * flows[:, 2]
    )
    explicit_total = float(np.sum(explicit_by_age))
    if abs(explicit_total - float(evaluation.births)) > 2e-10:
        raise RuntimeError(
            "Dated parity birth flows do not reproduce aggregate births: "
            f"flows={explicit_total:.12g}, aggregate={evaluation.births:.12g}"
        )
    topcode_accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        P,
    )
    topcode_total = float(np.sum(topcode_by_age))
    if abs(
        topcode_total
        - float(topcode_accounting["topcode_adjusted_birth_children"])
    ) > 2e-10:
        raise RuntimeError("Dated parity flows do not reproduce top-code-adjusted births")
    explicit_rates = np.divide(
        explicit_by_age,
        age_mass,
        out=np.zeros_like(explicit_by_age),
        where=age_mass > 1e-15,
    )
    topcode_rates = np.divide(
        topcode_by_age,
        age_mass,
        out=np.zeros_like(topcode_by_age),
        where=age_mass > 1e-15,
    )
    age_cell_starts = (
        float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    )
    ages = age_cell_starts + 0.5 * float(P.da)
    first_birth_accounting = first_birth_accounting_by_age(evaluation, P)
    if np.max(np.abs(first_birth_accounting["flow"] - flows[:, 0])) > 2e-12:
        raise RuntimeError("First-birth flow diagnostics disagree")
    first_birth_rates = np.divide(
        flows[:, 0],
        age_mass,
        out=np.zeros(int(P.J), dtype=float),
        where=age_mass > 1e-15,
    )
    first_total = float(np.sum(flows[:, 0]))
    return {
        "period_tfr_explicit_states": float(np.sum(explicit_rates)),
        "period_tfr_topcode_adjusted": float(np.sum(topcode_rates)),
        "period_first_birth_mean_age": float(
            np.sum(ages * flows[:, 0]) / max(first_total, 1e-15)
        ),
        "period_first_birth_share_age30plus": float(
            np.sum(flows[ages >= 30.0, 0]) / max(first_total, 1e-15)
        ),
        "age_cell_start": age_cell_starts,
        "age_cell_midpoint": ages,
        # Backward-compatible display field, now explicitly midpoint-coded.
        "age": ages,
        "age_mass": age_mass,
        "birth_flow_first": flows[:, 0],
        "birth_flow_second": flows[:, 1],
        "birth_flow_third_bin_entry": flows[:, 2],
        "birth_flow_explicit": explicit_by_age,
        "birth_flow_topcode_adjusted": topcode_by_age,
        "childless_at_risk_mass": first_birth_accounting["at_risk"],
        "first_birth_hazard": first_birth_accounting["hazard"],
        "first_birth_rate_all_households": first_birth_rates,
        "age_specific_birth_rate_explicit": explicit_rates,
        "age_specific_birth_rate_topcode_adjusted": topcode_rates,
    }


def first_birth_housing_response(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> float:
    """One-period birth-versus-no-birth housing response on the dated risk set.

    The PSID row is the change from event time -1 to +3.  One model period is
    four years, so both branches are advanced once from the same pre-birth
    state before their housing choices are compared.  The branch holds the
    current temporary-equilibrium policy fixed and deliberately does not add a
    second birth inside the four-year window.
    """
    horizon = int(getattr(P, "housing_event_horizon", 0))
    if horizon != E5_HOUSING_EVENT_HORIZON:
        raise ValueError(
            "The E5 transition calibration requires a one-period housing-event horizon"
        )
    model = calendar.model
    policy = evaluation.policy
    fecundity = model.get_fecundity_by_age(P)
    settled = model.readiness_settled_state(P)
    _, _, income_transition = model.income_transition_values(P)
    stochastic_child_aging = bool(
        getattr(P, "use_stochastic_aging", False) and hasattr(P, "Pi_child")
    )
    child_transition = P.Pi_child if stochastic_child_aging else None
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
        if j + horizon >= int(P.J):
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
            birth_cohort = model.advance_cohort_horizon_markov_income(
                birth_cohort,
                j,
                horizon,
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
                stochastic_child_aging,
                child_transition,
                income_transition,
            )
            control_cohort = model.advance_cohort_horizon_markov_income(
                control_cohort,
                j,
                horizon,
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
                stochastic_child_aging,
                child_transition,
                income_transition,
            )
            observed_birth = model.realize_current_choices_markov_income(
                birth_cohort,
                j + horizon,
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
                j + horizon,
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
                observed_birth, j + horizon, policy.hR_pol, P
            )
            control_housing = model.mean_housing_distribution_markov(
                observed_control, j + horizon, policy.hR_pol, P
            )
            response_sum += mass * (birth_housing - control_housing)
            response_mass += mass
    return response_sum / max(response_mass, 1e-15)


def begin_dated_first_birth_housing_branch(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    *,
    origin_period: int,
) -> dict[str, Any]:
    """Create equal birth/control branches and carry them one dated period.

    The branches are selected from successful first births in the origin-date
    pre-fertility risk set.  Survival and the origin-date location, tenure,
    saving, income, and child-maturation kernels carry both branches to the
    next date.  The aggregate Census age bridge is intentionally not applied
    to this fixed cohort.
    """
    model = calendar.model
    policy = evaluation.policy
    fecundity = model.get_fecundity_by_age(P)
    settled = model.readiness_settled_state(P)
    treated = np.zeros_like(evaluation.g_pre)
    control = np.zeros_like(evaluation.g_pre)
    for j in range(int(P.J) - 1):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        for zz in range(evaluation.g_pre.shape[4]):
            childless = evaluation.g_pre[:, :, :, j, zz, 0, settled]
            attempt = policy.fert_probs[:, :, :, j, zz, 1]
            realized = float(fecundity[j]) * childless * attempt
            treated[:, :, :, j, zz, 1, 1] = realized
            control[:, :, :, j, zz, 0, settled] = realized

    origin_mass = float(np.sum(treated))
    if origin_mass <= 1e-14:
        raise RuntimeError("The dated first-birth branch has zero treated mass")
    if not math.isclose(
        float(np.sum(control)), origin_mass, rel_tol=0.0, abs_tol=2e-13
    ):
        raise RuntimeError("Birth and control branches do not have equal origin mass")

    _, _, income_transition = model.income_transition_values(P)
    stochastic_child_aging = bool(
        getattr(P, "use_stochastic_aging", False) and hasattr(P, "Pi_child")
    )
    child_transition = P.Pi_child if stochastic_child_aging else None

    def advance(
        branch: np.ndarray, label: str
    ) -> tuple[np.ndarray, float, dict[str, float | str]]:
        next_pre = np.zeros_like(branch)
        deaths = 0.0
        for j in range(int(P.J) - 1):
            cohort = branch[:, :, :, j, :, :, :]
            survival = (
                float(P.survival_probs[j])
                if bool(getattr(P, "use_age_survival", False))
                else 1.0
            )
            deaths += (1.0 - survival) * float(np.sum(cohort))
            next_pre[:, :, :, j + 1, :, :, :] = (
                model.advance_cohort_one_period_markov_income(
                    survival * cohort,
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
                    stochastic_child_aging,
                    child_transition,
                    income_transition,
                )
            )
        expected_survivor_mass = float(np.sum(branch)) - deaths
        next_pre, mass_gate = normalize_distribution_mass_roundoff(
            next_pre,
            expected_mass=expected_survivor_mass,
            stage=f"{label}_branch_advancement",
        )
        return next_pre, deaths, mass_gate

    treated_next, treated_deaths, treated_advance_gate = advance(treated, "treated")
    control_next, control_deaths, control_advance_gate = advance(control, "control")
    treated_next_mass = float(np.sum(treated_next))
    control_next_mass = float(np.sum(control_next))
    if not math.isclose(
        treated_next_mass, control_next_mass, rel_tol=0.0, abs_tol=2e-11
    ):
        raise RuntimeError(
            "Birth and control branch masses differ after advancement: "
            f"treated={treated_next_mass:.17g}, "
            f"control={control_next_mass:.17g}, "
            f"gap={treated_next_mass - control_next_mass:.3e}, "
            f"origin={origin_mass:.17g}"
        )
    if not math.isclose(
        treated_deaths, control_deaths, rel_tol=0.0, abs_tol=2e-13
    ):
        raise RuntimeError("Birth and control branches have different survival losses")
    return {
        "status": "pending_one_period_dated_branch",
        "origin_period": int(origin_period),
        "origin_price": float(policy.price[0]),
        "origin_mass": origin_mass,
        "survivor_mass_before_destination_gating": treated_next_mass,
        "survival_loss": treated_deaths,
        "treated_next_pre": treated_next,
        "control_next_pre": control_next,
        "branch_mass_roundoff_gates": [
            treated_advance_gate,
            control_advance_gate,
        ],
    }


def finish_dated_first_birth_housing_branch(
    branch: dict[str, Any],
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    *,
    destination_period: int,
) -> dict[str, Any]:
    """Measure the prior date's branch with the destination policy and price."""
    origin_period = int(branch["origin_period"])
    if int(destination_period) != origin_period + 1:
        raise ValueError("The first-birth housing branch must span exactly one period")
    treated_pre, treated_projection = calendar.gate_pre_fertility_distribution(
        np.asarray(branch["treated_next_pre"], dtype=float),
        evaluation.policy,
        P,
        b_grid,
        shared,
    )
    control_pre, control_projection = calendar.gate_pre_fertility_distribution(
        np.asarray(branch["control_next_pre"], dtype=float),
        evaluation.policy,
        P,
        b_grid,
        shared,
    )
    expected_destination_mass = float(
        branch["survivor_mass_before_destination_gating"]
    )
    treated_pre, treated_gate = normalize_distribution_mass_roundoff(
        treated_pre,
        expected_mass=expected_destination_mass,
        stage="treated_destination_feasibility_gate",
    )
    control_pre, control_gate = normalize_distribution_mass_roundoff(
        control_pre,
        expected_mass=expected_destination_mass,
        stage="control_destination_feasibility_gate",
    )
    # Only the treated branch can take a continuation birth at the destination
    # date.  The empirical comparison group is confirmed childless, so the
    # otherwise identical control branch is held childless by construction.
    treated_post, continuation_births, _ = transition.apply_sequential_fertility(
        treated_pre,
        evaluation.policy.fert_probs,
        P,
        calendar.policy_continuation_birth_probs(evaluation.policy, P),
    )
    treated_post, treated_fertility_gate = normalize_distribution_mass_roundoff(
        treated_post,
        expected_mass=expected_destination_mass,
        stage="treated_destination_fertility",
    )
    control_parent_mass = float(np.sum(control_pre[..., 1:, :]))
    if control_parent_mass > 2e-13:
        raise RuntimeError("The confirmed-childless control branch acquired parent mass")
    treated_current = calendar.model.realize_current_cross_section(
        treated_post,
        evaluation.policy.loc_probs,
        evaluation.policy.tenure_choice,
        evaluation.policy.tenure_probs,
        evaluation.policy.maps.lmm_idx,
        evaluation.policy.maps.lmm_wt,
        evaluation.policy.maps.tmx_idx,
        evaluation.policy.maps.tmx_wt,
        use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
    )
    control_current = calendar.model.realize_current_cross_section(
        control_pre,
        evaluation.policy.loc_probs,
        evaluation.policy.tenure_choice,
        evaluation.policy.tenure_probs,
        evaluation.policy.maps.lmm_idx,
        evaluation.policy.maps.lmm_wt,
        evaluation.policy.maps.tmx_idx,
        evaluation.policy.maps.tmx_wt,
        use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
    )
    treated_current, treated_current_gate = normalize_distribution_mass_roundoff(
        treated_current,
        expected_mass=expected_destination_mass,
        stage="treated_destination_current_choice",
    )
    control_current, control_current_gate = normalize_distribution_mass_roundoff(
        control_current,
        expected_mass=expected_destination_mass,
        stage="control_destination_current_choice",
    )
    treated_mass = float(np.sum(treated_current))
    control_mass = float(np.sum(control_current))
    if not math.isclose(treated_mass, control_mass, rel_tol=0.0, abs_tol=2e-11):
        raise RuntimeError("Birth and control branch masses differ at destination")
    if treated_mass <= 1e-14:
        raise RuntimeError("The dated first-birth branch has zero destination mass")
    treated_housing = float(
        np.sum(calendar.housing_demand_by_location(treated_current, evaluation.policy.hR_pol, P))
        / treated_mass
    )
    control_housing = float(
        np.sum(calendar.housing_demand_by_location(control_current, evaluation.policy.hR_pol, P))
        / control_mass
    )
    return {
        "status": "complete_one_period_dated_matched_branch_did",
        "origin_period": origin_period,
        "destination_period": int(destination_period),
        "origin_price": float(branch["origin_price"]),
        "destination_price": float(evaluation.policy.price[0]),
        "origin_mass": float(branch["origin_mass"]),
        "destination_mass": treated_mass,
        "survival_loss": float(branch["survival_loss"]),
        "treated_feasibility_projection_mass": float(treated_projection),
        "control_feasibility_projection_mass": float(control_projection),
        "treated_continuation_births": float(continuation_births),
        "treated_mean_housing": treated_housing,
        "control_mean_housing": control_housing,
        "housing_response": treated_housing - control_housing,
        "census_age_bridge_applied": False,
        "branch_mass_roundoff_gates": [
            *list(branch.get("branch_mass_roundoff_gates") or []),
            treated_gate,
            control_gate,
            treated_fertility_gate,
            treated_current_gate,
            control_current_gate,
        ],
    }


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
    shared: SimpleNamespace,
    target_names: tuple[str, ...],
    *,
    housing_increment_override: float | None = None,
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
    housing_response = (
        first_birth_housing_response(evaluation, P, b_grid, shared)
        if housing_increment_override is None
        else float(housing_increment_override)
    )
    if not math.isfinite(housing_response):
        raise RuntimeError("The first-birth housing response is nonfinite")
    stats.housing_increment_0to1_eventstudy_t3 = housing_response
    extracted = extract_moments(stats, P)
    extracted.update(target_age_fertility_moments(evaluation.g_current, P))
    moments = {name: float(extracted[name]) for name in target_names}
    nonfinite = [name for name, value in moments.items() if not math.isfinite(value)]
    if nonfinite:
        raise RuntimeError(f"Transition target moments are nonfinite: {nonfinite}")
    return moments


def cohort_timing_moments(
    old_stationary_accounting: dict[str, np.ndarray],
    period_records: dict[int, dict[str, Any]],
    P: SimpleNamespace,
    terminal_period: int = TRANSITION_PERIODS,
    target_age: float = TARGET_AGE,
    terminal_childless_probability: float | None = None,
    terminal_identity_tolerance: float = 2.0e-9,
) -> dict[str, Any]:
    """First-birth timing for a fixed synthetic cohort reaching age 42.

    Dated Census age reweighting changes the mass of an age cell. We therefore
    splice age-specific hazards, not raw birth counts, and propagate one common
    childless risk set through the cohort's life.
    """
    model = calendar.model
    target_index = model.age_to_index(P, target_age)
    initial_period = max(0, int(terminal_period) - target_index)
    initial_index = max(0, target_index - int(terminal_period))
    age_cell_starts = (
        float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
    )
    ages = age_cell_starts + 0.5 * float(P.da)
    rows: list[dict[str, Any]] = []
    hazards: list[float] = []
    for age_index in range(initial_index):
        hazard = float(old_stationary_accounting["hazard"][age_index])
        hazards.append(hazard)
        rows.append(
            {
                "source": "old_steady_state_prehistory",
                "period": age_index - initial_index,
                "age_index": age_index,
                "age_cell_start": float(age_cell_starts[age_index]),
                "age": float(ages[age_index]),
                "source_at_risk_mass": float(
                    old_stationary_accounting["at_risk"][age_index]
                ),
                "source_first_birth_flow": float(
                    old_stationary_accounting["flow"][age_index]
                ),
                "first_birth_hazard": hazard,
            }
        )
    for period in range(initial_period, int(terminal_period) + 1):
        age_index = initial_index + period - initial_period
        accounting = period_records[period]["first_birth_accounting_by_age"]
        hazard = float(np.asarray(accounting["hazard"], dtype=float)[age_index])
        hazards.append(hazard)
        rows.append(
            {
                "source": "simulated_transition",
                "period": period,
                "age_index": age_index,
                "age_cell_start": float(age_cell_starts[age_index]),
                "age": float(ages[age_index]),
                "source_at_risk_mass": float(
                    np.asarray(accounting["at_risk"], dtype=float)[age_index]
                ),
                "source_first_birth_flow": float(
                    np.asarray(accounting["flow"], dtype=float)[age_index]
                ),
                "first_birth_hazard": hazard,
            }
        )
    synthetic_childless = 1.0
    for row, hazard in zip(rows, hazards, strict=True):
        synthetic_flow = synthetic_childless * hazard
        row["synthetic_childless_risk_set"] = synthetic_childless
        row["synthetic_first_birth_probability"] = synthetic_flow
        synthetic_childless *= 1.0 - hazard
    total = sum(float(row["synthetic_first_birth_probability"]) for row in rows)
    raw_total = total
    raw_childless = synthetic_childless
    identity_normalization = None
    if terminal_childless_probability is not None:
        terminal_childless = float(terminal_childless_probability)
        tolerance = float(terminal_identity_tolerance)
        gap = raw_childless - terminal_childless
        if (
            not math.isfinite(terminal_childless)
            or not 0.0 <= terminal_childless <= 1.0
            or not math.isfinite(tolerance)
            or tolerance <= 0.0
            or abs(gap) > tolerance
            or total <= 0.0
        ):
            raise RuntimeError(
                "Hazard-spliced cohort childlessness does not match the terminal "
                "age-42 stock before numerical identity normalization: "
                f"gap={gap:.3e}, tolerance={tolerance:.3e}"
            )
        target_total = 1.0 - terminal_childless
        probability_scale = target_total / total
        for row in rows:
            raw_probability = float(row["synthetic_first_birth_probability"])
            row[
                "synthetic_first_birth_probability_before_terminal_identity_normalization"
            ] = raw_probability
            row["synthetic_first_birth_probability"] = (
                probability_scale * raw_probability
            )
        total = target_total
        synthetic_childless = terminal_childless
        identity_normalization = {
            "status": "terminal_stock_identity_after_2e-9_roundoff_gate",
            "raw_synthetic_childless_probability": raw_childless,
            "terminal_childless_probability": terminal_childless,
            "raw_gap": gap,
            "absolute_tolerance": tolerance,
            "first_birth_probability_scale": probability_scale,
        }
    mean_age = sum(
        float(row["age"]) * float(row["synthetic_first_birth_probability"])
        for row in rows
    ) / max(total, 1e-15)
    late_share = sum(
        float(row["synthetic_first_birth_probability"])
        for row in rows
        if float(row["age"]) >= 30.0
    ) / max(total, 1e-15)
    return {
        "target_cohort_initial_age": float(age_cell_starts[initial_index]),
        "target_cohort_terminal_age": float(age_cell_starts[target_index]),
        "target_cohort_initial_age_cell_midpoint": float(ages[initial_index]),
        "target_cohort_terminal_age_cell_midpoint": float(ages[target_index]),
        "first_birth_age_measurement": (
            "four-year model cells labeled by interval midpoints"
        ),
        "synthetic_ever_first_birth_probability": total,
        "synthetic_childless_probability": synthetic_childless,
        "synthetic_ever_first_birth_probability_before_terminal_identity_normalization": raw_total,
        "synthetic_childless_probability_before_terminal_identity_normalization": raw_childless,
        "terminal_childless_identity_normalization": identity_normalization,
        "mean_age_first_birth": mean_age,
        "share_first_births_age30plus": late_share,
        "hazard_rows": rows,
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


def parameter_rows(
    theta: dict[str, float],
    best: dict[str, Any],
    *,
    old_psi_child: float,
    joint_panel: bool,
) -> list[dict[str, Any]]:
    """Report the actual dated-transition coordinates and their true bounds."""
    rows: list[dict[str, Any]] = []
    for name, lower, upper, transform in TRANSITION_SEARCH_DOMAIN:
        if name == "beta_annual":
            estimate = float(theta["beta"]) ** 0.25
        elif name == "psi_child_change_2023":
            estimate = float(best["new_psi_child"]) - float(old_psi_child)
        elif name == "psi_child_2023":
            estimate = float(best["new_psi_child"])
        else:
            estimate = float(theta[name])
        rows.append(
            {
                "parameter": name,
                "value": estimate,
                "lower_bound": lower,
                "upper_bound": upper,
                "transform": transform,
                "is_free_parameter": True,
                "status": (
                    "estimated_free_transition_parameter"
                    if joint_panel
                    else "fixed_or_bounded_pilot_transition_parameter"
                ),
                "near_bound": bool(
                    min(abs(estimate - lower), abs(upper - estimate))
                    <= 0.02 * max(upper - lower, 1e-12)
                ),
            }
        )
    rows.extend(
        [
            {
                "parameter": "psi_child_2007",
                "value": math.nan,
                "lower_bound": math.nan,
                "upper_bound": math.nan,
                "transform": "derived",
                "is_free_parameter": False,
                "status": "externally_normalized_to_old_completed_fertility",
                "near_bound": False,
            },
            {
                "parameter": "psi_child_2023",
                "value": math.nan,
                "lower_bound": math.nan,
                "upper_bound": math.nan,
                "transform": "derived",
                "is_free_parameter": False,
                "status": "derived_from_old_intercept_and_transition_coordinate",
                "near_bound": False,
            },
        ]
    )
    return rows


def make_diagnostic_plot(
    summaries: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    outdir: Path,
    target_count: int,
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
        ylabel=f"{int(target_count)}-moment transition loss",
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
    if args.housing_supply_elasticity is not None and (
        not math.isfinite(float(args.housing_supply_elasticity))
        or float(args.housing_supply_elasticity) < 0.0
    ):
        raise ValueError("--housing-supply-elasticity must be finite and nonnegative")
    if args.fixed_tenure_choice_kappa is not None and (
        not math.isfinite(float(args.fixed_tenure_choice_kappa))
        or float(args.fixed_tenure_choice_kappa) < 0.0
    ):
        raise ValueError("--fixed-tenure-choice-kappa must be finite and nonnegative")
    replacement_fertility = float(args.replacement_fertility)
    renewal_conversion = transition.effective_birth_to_household_conversion(
        replacement_fertility
    )
    if not math.isclose(
        replacement_fertility,
        transition.REPLACEMENT_FERTILITY,
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        raise ValueError("Production calibration is hard-gated to replacement fertility 2.1")
    old_completed_fertility_target = (
        replacement_fertility
        if args.old_completed_fertility_target is None
        else float(args.old_completed_fertility_target)
    )
    if not math.isclose(
        old_completed_fertility_target,
        replacement_fertility,
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        raise ValueError(
            "The old steady state must be normalized to the declared replacement fertility"
        )

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
    code_contract = code_fingerprint_contract(model)
    if (
        args.expected_code_bundle_sha256 is not None
        and str(code_contract["bundle_sha256"])
        != str(args.expected_code_bundle_sha256)
    ):
        raise RuntimeError(
            "Code-bundle fingerprint gate failed: "
            f"actual={code_contract['bundle_sha256']}, "
            f"expected={args.expected_code_bundle_sha256}"
        )
    winner, source_metadata = closure.load_winner(source)
    theta = closure.theta_from_winner(winner)
    active_domain, profile_overrides, model_profile = activate_model_profile(
        str(args.model_profile), theta
    )
    active_domain, profile_overrides, model_profile = configure_first_child_room_jump(
        model_profile_name=str(args.model_profile),
        theta=theta,
        active_domain=active_domain,
        profile_overrides=profile_overrides,
        model_profile=model_profile,
        fixed_jump=args.fixed_first_child_room_jump,
        estimate_jump=bool(args.estimate_first_child_room_jump),
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
    retained_tenure_choice_kappa = float(theta.get("tenure_choice_kappa", 0.0))
    if args.fixed_tenure_choice_kappa is not None:
        theta["tenure_choice_kappa"] = float(args.fixed_tenure_choice_kappa)
        profile_overrides["tenure_choice_kappa"] = float(
            args.fixed_tenure_choice_kappa
        )
    model_profile["tenure_choice_kappa"] = {
        "value": float(theta.get("tenure_choice_kappa", retained_tenure_choice_kappa)),
        "retained_value": retained_tenure_choice_kappa,
        "status": (
            "externally_fixed_profile_not_estimated"
            if args.fixed_tenure_choice_kappa is not None
            else "retained_model_value"
        ),
    }
    target_system = e5_target_system_for_profile(str(args.target_profile))
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
    if (
        args.fixed_tenure_choice_kappa is not None
        and float(args.fixed_tenure_choice_kappa) > 0.0
    ):
        base["normalize_transition_mass_roundoff"] = True

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
        completed_fertility_target=old_completed_fertility_target,
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
        old_shared,
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
    kfe_B_old_raw = float(old_solution.entrants_mature_total)
    kfe_B_old_adjusted = float(
        old_renewal["topcode_adjusted_mature_entrant_households"]
    )
    births_old_raw = float(old_solution.total_births_kfe)
    births_old = float(old_renewal["topcode_adjusted_birth_children"])
    B_old_raw = renewal_conversion * births_old_raw
    outside_share = float(args.outside_origin_entry_share)
    renewal_old_state = transition.stationary_renewal_from_births(
        E_old,
        births_old,
        outside_share,
        renewal_conversion,
    )
    B_old = float(renewal_old_state["queue_mature_flow_B"])
    outside_flow = float(renewal_old_state["outside_flow_M"])
    retention = float(renewal_old_state["retention_rho"])
    old_queue_B_over_E = float(renewal_old_state["queue_B_over_E"])
    old_renewal_residual = float(renewal_old_state["identity_residual"])
    old_fertility_ratio = float(old_transition_moments["tfr"]) / replacement_fertility
    if abs(old_queue_B_over_E - old_fertility_ratio) > 2e-8:
        raise RuntimeError(
            "Old renewal flow does not match the stationary fertility ratio: "
            f"B/E={old_queue_B_over_E:.12g}, "
            f"completed/replacement={old_fertility_ratio:.12g}"
        )
    if abs(old_queue_B_over_E - 1.0) > float(args.old_completed_fertility_tol):
        raise RuntimeError(
            "Old renewal state is not at replacement: "
            f"B/E={old_queue_B_over_E:.12g}"
        )
    if abs(old_renewal_residual) > 2e-12:
        raise RuntimeError(
            f"Old renewal accounting does not nest exactly: {old_renewal_residual:.3e}"
        )
    stationary_initial_g_pre = initial_g_pre
    model_ages = (
        float(old_parameters.age_start)
        + np.arange(int(old_parameters.J), dtype=float) * float(old_parameters.da)
    )
    stationary_age_mass = np.sum(
        stationary_initial_g_pre, axis=(0, 1, 2, 4, 5, 6)
    )
    implied_age_survival = stationary_age_mass[1:] / stationary_age_mass[:-1]
    structural_age_survival = (
        np.asarray(old_parameters.survival_probs, dtype=float)[: int(old_parameters.J) - 1]
        if bool(old_parameters.use_age_survival)
        else np.ones(int(old_parameters.J) - 1, dtype=float)
    )
    age_survival_diagnostic = {
        "stationary_age_mass": stationary_age_mass.tolist(),
        "implied_age_survival": implied_age_survival.tolist(),
        "structural_age_survival": structural_age_survival.tolist(),
        "implied_minus_structural": (
            implied_age_survival - structural_age_survival
        ).tolist(),
        "max_implied_age_survival": float(np.max(implied_age_survival)),
        "max_abs_implied_minus_structural": float(
            np.max(np.abs(implied_age_survival - structural_age_survival))
        ),
        "stationary_entry_flow": E_old,
        "youngest_age_mass": float(stationary_age_mass[0]),
    }
    bridge_structural_survival = None
    if args.fixed_tenure_choice_kappa is not None:
        (
            structural_stationary_age_mass,
            structural_roundoff_gate,
        ) = validated_structural_stationary_age_mass(
            stationary_age_mass,
            entry_flow=E_old,
            structural_survival=structural_age_survival,
        )
        age_survival_diagnostic["structural_roundoff_gate"] = (
            structural_roundoff_gate
        )
        age_survival_diagnostic["bridge_stationary_age_mass"] = (
            stationary_age_mass.tolist()
        )
        age_survival_diagnostic["exact_structural_stationary_age_mass"] = (
            structural_stationary_age_mass.tolist()
        )
        bridge_structural_survival = structural_age_survival
    write_json(outdir / "stationary_age_survival_diagnostic.json", age_survival_diagnostic)
    print(
        "STATIONARY_AGE_SURVIVAL_DIAGNOSTIC "
        f"max_implied={age_survival_diagnostic['max_implied_age_survival']:.12g} "
        "max_gap="
        f"{age_survival_diagnostic['max_abs_implied_minus_structural']:.3e}",
        flush=True,
    )
    age_reweight = transition.acs_2007_age_reweight_diagnostic(
        stationary_age_mass,
        model_ages,
        E_old,
        periods=TRANSITION_PERIODS,
        period_years=float(old_parameters.period_years),
        structural_survival=bridge_structural_survival,
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
        "externally_fixed_profile_not_estimated"
        if args.housing_supply_elasticity is not None
        else "retained_model_value"
    )
    old_first_birth_accounting = first_birth_accounting_by_age(
        old_evaluation, old_parameters
    )
    population_target_indices = transition.census_household_target_indices()
    population_age_target_years = transition.census_household_age_target_years()
    population_target = float(population_target_indices[TRANSITION_PERIODS])

    all_fit_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    target_measurements = {
        "tfr": "completed fertility of the model cohort centered at age 42 in 2023",
        "childless_rate": "childless share of the model cohort centered at age 42 in 2023",
        "mean_age_first_birth": (
            "four-year age-cell midpoint mean for a fixed synthetic cohort built "
            "from dated first-birth hazards: old-steady-state prehistory plus "
            "the 2007--2023 transition"
        ),
        "share_first_births_age30plus": "same hazard-based fixed-cohort ledger",
        "housing_increment_0to1": (
            "2019--2023 matched birth-versus-confirmed-childless branch: 2019 "
            "policies advance the fixed cohort, 2023 policies govern continuation "
            "fertility and housing, and the Census age bridge is not applied"
        ),
        "remaining_targets": "2023 transition cross-section with beginning-of-period wealth timing",
    }
    if "own_rate_2534" in target_system.moment_names:
        target_measurements["own_rate_2534"] = (
            "2023 transition cross-section: owner mass divided by total household "
            "mass over model ages 25--34"
        )
    for index, new_psi in enumerate(candidates, start=1):
        label = (
            f"task_{int(args.panel_task_id):03d}"
            if args.panel_task_id is not None
            else f"psi_{new_psi:+.6f}".replace("+", "p").replace("-", "m")
        )
        candidate_dir = outdir / "cases" / label
        candidate_dir.mkdir(parents=True, exist_ok=True)
        period_records: dict[int, dict[str, Any]] = {}
        pending_housing_branch: dict[str, Any] | None = None

        def observer(
            period: int,
            evaluation: calendar.PeriodEvaluation,
            parameters: SimpleNamespace,
            grid: np.ndarray,
            shared: SimpleNamespace,
        ) -> None:
            nonlocal pending_housing_branch
            if pending_housing_branch is None:
                same_policy_response = first_birth_housing_response(
                    evaluation, parameters, grid, shared
                )
                housing_measurement = {
                    "status": "same_policy_one_period_diagnostic_no_prior_date",
                    "origin_period": None,
                    "destination_period": int(period),
                    "housing_response": same_policy_response,
                    "census_age_bridge_applied": False,
                }
            else:
                housing_measurement = finish_dated_first_birth_housing_branch(
                    pending_housing_branch,
                    evaluation,
                    parameters,
                    grid,
                    shared,
                    destination_period=int(period),
                )
            period_records[int(period)] = {
                "moments": transition_cross_section_moments(
                    evaluation,
                    parameters,
                    grid,
                    shared,
                    target_system.moment_names,
                    housing_increment_override=float(
                        housing_measurement["housing_response"]
                    ),
                ),
                "first_birth_housing_did": housing_measurement,
                "first_birth_accounting_by_age": first_birth_accounting_by_age(
                    evaluation, parameters
                ),
                "period_fertility_diagnostics": period_fertility_diagnostics(
                    evaluation, parameters
                ),
            }
            pending_housing_branch = (
                begin_dated_first_birth_housing_branch(
                    evaluation,
                    parameters,
                    grid,
                    shared,
                    origin_period=int(period),
                )
                if int(period) < total_periods - 1
                else None
            )

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
            baseline_births=births_old,
            baseline_raw_births=births_old_raw,
            parameters=copy.deepcopy(old_parameters),
            b_grid=b_grid,
            periods=total_periods,
            outside_flow=outside_flow,
            retention=retention,
            effective_birth_to_household_conversion=renewal_conversion,
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
        dated_moment_rows: list[dict[str, Any]] = []
        dated_fertility_rows: list[dict[str, Any]] = []
        dated_fertility_age_rows: list[dict[str, Any]] = []
        dated_timing_ledgers: list[dict[str, Any]] = []
        dated_housing_did_rows: list[dict[str, Any]] = []
        target_lookup = {
            name: (float(target), float(weight))
            for name, target, weight in zip(
                target_system.moment_names,
                target_system.target_values,
                target_system.weights,
                strict=True,
            )
        }
        for period in range(total_periods):
            period_timing = cohort_timing_moments(
                old_first_birth_accounting,
                period_records,
                old_parameters,
                terminal_period=period,
                terminal_childless_probability=float(
                    period_records[period]["moments"]["childless_rate"]
                ),
            )
            dated_timing_ledgers.append(
                {
                    "period": period,
                    "calendar_year": START_YEAR
                    + period * float(old_parameters.period_years),
                    **period_timing,
                }
            )
            dated_moments = dict(period_records[period]["moments"])
            housing_did = dict(period_records[period]["first_birth_housing_did"])
            dated_housing_did_rows.append(
                {
                    "candidate": label,
                    "period": period,
                    "calendar_year": START_YEAR
                    + period * float(old_parameters.period_years),
                    "status": housing_did["status"],
                    "origin_period": housing_did.get("origin_period"),
                    "destination_period": housing_did.get("destination_period"),
                    "origin_price": housing_did.get("origin_price"),
                    "destination_price": housing_did.get("destination_price"),
                    "origin_mass": housing_did.get("origin_mass"),
                    "destination_mass": housing_did.get("destination_mass"),
                    "survival_loss": housing_did.get("survival_loss"),
                    "treated_feasibility_projection_mass": housing_did.get(
                        "treated_feasibility_projection_mass"
                    ),
                    "control_feasibility_projection_mass": housing_did.get(
                        "control_feasibility_projection_mass"
                    ),
                    "treated_continuation_births": housing_did.get(
                        "treated_continuation_births"
                    ),
                    "treated_mean_housing": housing_did.get(
                        "treated_mean_housing"
                    ),
                    "control_mean_housing": housing_did.get(
                        "control_mean_housing"
                    ),
                    "housing_response": housing_did["housing_response"],
                    "census_age_bridge_applied": housing_did[
                        "census_age_bridge_applied"
                    ],
                }
            )
            dated_moments["mean_age_first_birth"] = float(
                period_timing["mean_age_first_birth"]
            )
            dated_moments["share_first_births_age30plus"] = float(
                period_timing["share_first_births_age30plus"]
            )
            for name in target_system.moment_names:
                target, weight = target_lookup[name]
                is_terminal = period == TRANSITION_PERIODS
                model_value = float(dated_moments[name])
                gap = model_value - target if is_terminal else math.nan
                dated_moment_rows.append(
                    {
                        "candidate": label,
                        "period": period,
                        "calendar_year": START_YEAR
                        + period * float(old_parameters.period_years),
                        "moment": name,
                        "model": model_value,
                        "role": (
                            "calibration_moment_at_2023_model_state"
                            if is_terminal
                            else "untargeted_transition_diagnostic"
                        ),
                        "target": target if is_terminal else None,
                        "weight": weight if is_terminal else None,
                        "gap": gap if is_terminal else None,
                        "loss_contribution": weight * gap**2 if is_terminal else None,
                    }
                )
            fertility = period_records[period]["period_fertility_diagnostics"]
            year = START_YEAR + period * float(old_parameters.period_years)
            dated_fertility_rows.append(
                {
                    "candidate": label,
                    "period": period,
                    "calendar_year": year,
                    "period_tfr_explicit_states": float(
                        fertility["period_tfr_explicit_states"]
                    ),
                    "period_tfr_topcode_adjusted": float(
                        fertility["period_tfr_topcode_adjusted"]
                    ),
                    "period_first_birth_mean_age": float(
                        fertility["period_first_birth_mean_age"]
                    ),
                    "period_first_birth_share_age30plus": float(
                        fertility["period_first_birth_share_age30plus"]
                    ),
                    "role": "untargeted_current-rate_diagnostic",
                    "denominator": "model pre-fertility adult-household mass by age",
                }
            )
            ages_diagnostic = np.asarray(fertility["age"], dtype=float)
            age_starts_diagnostic = np.asarray(
                fertility["age_cell_start"], dtype=float
            )
            for age_index, age in enumerate(ages_diagnostic):
                dated_fertility_age_rows.append(
                    {
                        "candidate": label,
                        "period": period,
                        "calendar_year": year,
                        "age_index": age_index,
                        "age_cell_start": float(
                            age_starts_diagnostic[age_index]
                        ),
                        "age_cell_midpoint": float(age),
                        "age": float(age),
                        "age_mass": float(fertility["age_mass"][age_index]),
                        "birth_flow_first": float(
                            fertility["birth_flow_first"][age_index]
                        ),
                        "birth_flow_second": float(
                            fertility["birth_flow_second"][age_index]
                        ),
                        "birth_flow_third_bin_entry": float(
                            fertility["birth_flow_third_bin_entry"][age_index]
                        ),
                        "birth_flow_explicit": float(
                            fertility["birth_flow_explicit"][age_index]
                        ),
                        "birth_flow_topcode_adjusted": float(
                            fertility["birth_flow_topcode_adjusted"][age_index]
                        ),
                        "childless_at_risk_mass": float(
                            fertility["childless_at_risk_mass"][age_index]
                        ),
                        "first_birth_hazard": float(
                            fertility["first_birth_hazard"][age_index]
                        ),
                        "first_birth_rate_all_households": float(
                            fertility["first_birth_rate_all_households"][age_index]
                        ),
                        "age_specific_birth_rate_explicit": float(
                            fertility["age_specific_birth_rate_explicit"][age_index]
                        ),
                        "age_specific_birth_rate_topcode_adjusted": float(
                            fertility["age_specific_birth_rate_topcode_adjusted"][age_index]
                        ),
                    }
                )
        calendar.write_csv(candidate_dir / "dated_model_moments.csv", dated_moment_rows)
        calendar.write_csv(
            candidate_dir / "dated_first_birth_housing_did.csv",
            dated_housing_did_rows,
        )
        calendar.write_csv(
            candidate_dir / "dated_period_fertility.csv", dated_fertility_rows
        )
        calendar.write_csv(
            candidate_dir / "dated_period_fertility_by_age.csv",
            dated_fertility_age_rows,
        )
        write_json(candidate_dir / "dated_cohort_timing_ledgers.json", dated_timing_ledgers)
        timing = dated_timing_ledgers[TRANSITION_PERIODS]
        terminal_moments = dict(period_records[TRANSITION_PERIODS]["moments"])
        terminal_moments["mean_age_first_birth"] = float(timing["mean_age_first_birth"])
        terminal_moments["share_first_births_age30plus"] = float(
            timing["share_first_births_age30plus"]
        )
        synthetic_childless_gap = float(timing["synthetic_childless_probability"]) - float(
            terminal_moments["childless_rate"]
        )
        if abs(synthetic_childless_gap) > 2e-10:
            raise RuntimeError(
                "Hazard-spliced cohort childlessness does not match the terminal "
                "age-42 stock: "
                f"gap={synthetic_childless_gap:.3e}"
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
            "terminal_first_birth_housing_response": float(
                terminal_moments["housing_increment_0to1"]
            ),
            "terminal_synthetic_childless_consistency_gap": synthetic_childless_gap,
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
    params = parameter_rows(
        theta,
        best,
        old_psi_child=old_psi_child,
        joint_panel=panel_metadata is not None,
    )
    params.extend(
        [
            {
                "parameter": "tenure_choice_kappa",
                "value": float(theta.get("tenure_choice_kappa", 0.0)),
                "lower_bound": math.nan,
                "upper_bound": math.nan,
                "transform": "external_profile",
                "is_free_parameter": False,
                "status": model_profile["tenure_choice_kappa"]["status"],
                "near_bound": False,
            },
            {
                "parameter": "housing_supply_elasticity",
                "value": float(supply_rule.elasticity),
                "lower_bound": math.nan,
                "upper_bound": math.nan,
                "transform": "external_profile",
                "is_free_parameter": False,
                "status": supply_normalization["elasticity_status"],
                "near_bound": False,
            },
        ]
    )
    for row in params:
        if row["parameter"] == "psi_child_2007":
            row["value"] = old_psi_child
        elif row["parameter"] == "psi_child_2023":
            row["value"] = float(best["new_psi_child"])
    calendar.write_csv(outdir / "parameter_table.csv", params)
    make_diagnostic_plot(
        summaries,
        all_fit_rows,
        outdir,
        target_count=target_system.count,
    )
    summary = {
        "status": (
            "complete_transition_calibration_panel_task"
            if panel_metadata is not None
            else "complete_bounded_fixed_parameter_pilot"
        ),
        "source": str(source),
        "source_sha256": actual_source_sha256,
        "source_metadata": source_metadata,
        "code_fingerprints": code_contract,
        "model_profile": model_profile,
        "income_profile_gates": income_profile_gates,
        "target_set": target_system.name,
        "target_profile": str(args.target_profile),
        "target_fingerprint": target_system.fingerprint,
        "target_count": target_system.count,
        "target_profile_metadata": (
            E5_YOUNG_OWNERSHIP_PROVENANCE
            if str(args.target_profile) == E5_YOUNG_OWNERSHIP_PROFILE
            else {}
        ),
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
        "fixed_first_child_room_jump": args.fixed_first_child_room_jump,
        "estimate_first_child_room_jump": bool(
            args.estimate_first_child_room_jump
        ),
        "external_closure_contract": {
            "housing_supply_elasticity": float(supply_rule.elasticity),
            "housing_supply_elasticity_status": supply_normalization[
                "elasticity_status"
            ],
            "tenure_choice_kappa": float(theta.get("tenure_choice_kappa", 0.0)),
            "tenure_choice_kappa_status": model_profile[
                "tenure_choice_kappa"
            ]["status"],
            "identification": (
                "both objects are externally fixed; the transition parameters "
                f"are disciplined by {target_system.count} dated moments"
            ),
        },
        "candidate_psi_changes": candidate_psi_changes,
        "policy_case": str(args.policy_case),
        "post_2023_periods": int(args.post_2023_periods),
        "old_psi_child": old_psi_child,
        "old_fertility_normalization": old_fertility_normalization,
        "old_model_completed_fertility": float(old_transition_moments["tfr"]),
        "old_completed_fertility_reference": old_completed_fertility_target,
        "outside_origin_entry_share": outside_share,
        "outside_origin_entry_status": "provisional_across_market_anchor_not_national_migration",
        "renewal_retention": retention,
        "renewal_accounting_contract": {
            "replacement_fertility": replacement_fertility,
            "effective_birth_to_household_conversion": renewal_conversion,
            "birth_measure": "topcode_adjusted_birth_children",
            "top_bin_mean_children": float(old_parameters.tfr_top_bin_weight),
            "birth_vintage_queue_waiting_slots": 4,
            "birth_to_entry_effect_lag_dates": 5,
            "birth_to_entry_effect_lag_years": 20.0,
            "outside_flow_formula": "M = outside_origin_entry_share * E_old",
            "retention_formula": "rho = (1 - outside_origin_entry_share) * E_old / B_old",
            "status": (
                "external replacement normalization; household KFE maturation "
                "is diagnostic and does not determine population renewal"
            ),
        },
        "renewal_accounting_old_state": {
            "old_entry_flow_E": E_old,
            "old_queue_mature_flow_B": B_old,
            "old_raw_birth_queue_flow_B": B_old_raw,
            "old_queue_B_over_E": old_queue_B_over_E,
            "old_renewal_residual": old_renewal_residual,
            "outside_flow_M": outside_flow,
            "outside_origin_entry_share": outside_share,
            "retention_rho": retention,
            "kfe_entrant_conversion_factor_diagnostic": float(
                old_parameters.entrant_conversion_factor
            ),
            "kfe_old_topcode_adjusted_mature_flow_diagnostic": kfe_B_old_adjusted,
            "kfe_old_raw_mature_flow_diagnostic": kfe_B_old_raw,
        },
        "housing_supply": supply_normalization,
        "target_measurements": target_measurements,
        "dated_measurement_contract": {
            "first_birth_housing": (
                "one-period matched branch DID; the origin-date policy carries "
                "equal birth and confirmed-childless branches through survival, "
                "saving, tenure, location, income, and child maturation; the "
                "destination policy governs continuation fertility and housing; "
                "the Census age bridge is excluded"
            ),
            "first_birth_housing_horizon_periods": 1,
            "first_birth_housing_horizon_years": 4.0,
            "first_birth_housing_terminal_origin_year": 2019,
            "first_birth_housing_terminal_destination_year": 2023,
            "cohort_first_birth_timing": (
                "age-specific first-birth hazards spliced into a fixed synthetic "
                "cohort; raw age-cell masses never weight timing"
            ),
            "first_birth_age_cells": (
                "four-year intervals [18,22),...,[42,46) labeled by midpoints "
                "20,24,...,44; the NCHS target applies the identical binned-age "
                "operator with boundary collapse"
            ),
            "period_fertility": (
                "sum over four-year age cells of dated birth flows divided by "
                "all pre-fertility adult-household mass in the age cell"
            ),
            "period_fertility_topcode_timing": (
                "the extra 0.602 children represented by the 3+ bin are imputed "
                "to the date and age of entry into parity 3; explicit-state TFR "
                "is the main historical diagnostic and the top-coded version is "
                "a sensitivity"
            ),
            "period_fertility_status": (
                "untargeted diagnostic; model household denominator is not literal "
                "national female exposure"
            ),
            "dated_maintained_moments": (
                "all twelve moments saved at every model date; only the 2023 "
                "model-state measurements enter the calibration loss"
            ),
        },
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
                "2007 householder-age profile. In panel mode the active structural parameters",
                "and the terminal child-preference change are estimated jointly on the",
                "dated 2023 target system. This remains a search packet, not a",
                "promoted estimate.",
                "",
                "The fertility stock targets use the model cohort centered at age 42.",
                "First-birth timing combines old-steady-state prehistory with dated",
                "first-birth hazards during 2007--2023 and labels four-year age cells",
                "by their midpoints. The previous period-TFR",
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
