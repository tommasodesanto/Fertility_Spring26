"""Direct no-inversion calibration objective.

This module keeps the live April 2026 target definitions, but treats the
center amenity and center rent shifter as ordinary calibrated parameters. The
geography moments therefore enter the same SMM loss as the fertility, tenure,
housing, wealth, and migration moments.

The current model uses the benchmark-normalized outside-option population
closure. Potential young adults choose between the modeled locations and an
outside economy; the benchmark size target normalizes the baseline city to one,
while the residual outside-origin pool is held fixed in counterfactuals.
"""

from __future__ import annotations

import copy
import time
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Mapping, Sequence

import numpy as np

from .objective import extract_moments
from .parameters import asdict, build_calibration_setup
from .solver import run_model_cp_dt
from .theta import apply_theta


DIRECT_GEOMETRY_NAMES = ["E_C", "r_bar_C"]
DIRECT_GEOMETRY_LB = np.array([-1.50, 0.010])
DIRECT_GEOMETRY_UB = np.array([2.50, 0.300])
OUTSIDE_VALUE_NAME = "outside_value"
OUTSIDE_VALUE_LB = -120.0
OUTSIDE_VALUE_UB = 20.0
DEFAULT_OUTSIDE_VALUE_X0 = -41.95
RENEWAL_FLOW_NAME = "outside_entry_flow"
RENEWAL_FLOW_LB = 1e-4
RENEWAL_FLOW_UB = 0.080
DEFAULT_RENEWAL_FLOW_X0 = 0.003
POPULATION_SCALE_CLOSURES = {"accounting_scale_prices", "scaled_housing", "scaled_housing_accounting"}
BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES = {
    "outside_option_benchmark_normalized",
    "benchmark_outside_option_normalized",
}
RENEWAL_VALVE_FIXED_CLOSURES = {"renewal_valve", "renewal_scale", "demographic_valve"}
RENEWAL_VALVE_CALIBRATED_CLOSURES = {
    "renewal_valve_calibrated",
    "renewal_calibrated",
    "demographic_valve_calibrated",
    "benchmark_renewal_valve",
}
RENEWAL_VALVE_CLOSURES = RENEWAL_VALVE_FIXED_CLOSURES | RENEWAL_VALVE_CALIBRATED_CLOSURES
OPEN_CITY_CLOSURES = (
    POPULATION_SCALE_CLOSURES
    | BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES
    | RENEWAL_VALVE_CLOSURES
    | {"outside_option", "outside_option_local_births", "local_births_outside", "open_city"}
)


@dataclass
class DirectCalibrationSetup:
    """Container for the no-inversion SMM setup."""

    targets: dict[str, float]
    weights: dict[str, float]
    P_base: SimpleNamespace
    names: list[str]
    lb: np.ndarray
    ub: np.ndarray
    x0: np.ndarray
    setup_mode: str
    population_closure: str


@dataclass
class DirectObjectiveResult:
    """One evaluated parameter vector."""

    loss: float
    theta: np.ndarray
    moments: dict[str, float]
    elapsed_sec: float
    solve_ok: bool
    converged: bool
    p_eq: list[float]
    timings: dict[str, Any]
    error: str | None = None


def build_direct_calibration_setup(
    setup_mode: str = "benchmark",
    geo_weight: float = 100.0,
    population_closure: str = "outside_option_benchmark_normalized",
    scale_target: float = 1.0,
    scale_weight: float = 100.0,
    outside_value_x0: float | None = None,
    outside_flow_x0: float | None = None,
    renewal_retention: float = 1.0,
    hR_max: float | None = None,
    owner_h_bar_scale: float | None = None,
    weight_overrides: Mapping[str, float] | None = None,
    parent_dp_waiver: bool | None = None,
    parent_dp_waiver_phi: float | None = None,
    H_own: Sequence[float] | None = None,
) -> DirectCalibrationSetup:
    """Build the direct-geometry counterpart of ``build_calibration_setup``.

    The first 13 parameters are unchanged. The next two parameters are:

    ``E_C``
        Center location amenity, with periphery normalized to zero.

    ``r_bar_C``
        Center housing supply rent shifter, with periphery left at the live
        baseline value.

    Under the default benchmark-normalized outside-option closure, no
    population-scale parameter is estimated. The outside value is chosen to
    hit the target city-entry probability, and the outside-born potential
    entrant flow is computed mechanically from the benchmark normalized
    distribution:

    ``M = S^* [E_0 / q^E - B_0]``.

    Under the fixed renewal-valve closure, the final parameter is:

    ``outside_entry_flow``
        Exogenous outside-born young-adult entry flow, M. It is a nuisance
        closure parameter calibrated to benchmark city scale and held fixed in
        counterfactuals.

    Under the older value-based scale closure, the final parameter is:

    ``outside_value``
        The value of the outside option. Kept for diagnostics, not the default
        benchmark closure.
    """

    base = build_calibration_setup(setup_mode)
    if hR_max is not None:
        base.P_base.hR_max = float(hR_max)
    if owner_h_bar_scale is not None:
        base.P_base.owner_h_bar_scale = float(owner_h_bar_scale)
    if parent_dp_waiver is not None:
        base.P_base.parent_dp_waiver = bool(parent_dp_waiver)
    if parent_dp_waiver_phi is not None:
        base.P_base.parent_dp_waiver_phi = float(parent_dp_waiver_phi)
    if H_own is not None:
        ladder = np.asarray(H_own, dtype=float).reshape(-1)
        if ladder.size == 0 or np.any(~np.isfinite(ladder)) or np.any(np.diff(ladder) <= 0):
            raise ValueError("H_own override must be a strictly increasing finite sequence")
        base.P_base.H_own = ladder
        base.P_base.n_house = int(ladder.size)
        base.P_base.h_own_min = float(ladder[0])
        base.P_base.h_own_max = float(ladder[-1])
    targets = dict(base.targets)
    weights = dict(base.weights)
    if weight_overrides:
        for key, value in weight_overrides.items():
            if key in targets:
                weights[key] = float(value)
    closure = str(population_closure or "outside_option_benchmark_normalized").lower()

    targets["inv_pop_share_C"] = float(base.inversion_targets["pop_share_C"])
    targets["inv_rent_ratio_C_over_P"] = float(base.inversion_targets["rent_ratio"])
    weights["inv_pop_share_C"] = float(geo_weight)
    weights["inv_rent_ratio_C_over_P"] = float(geo_weight)

    names = list(base.names) + DIRECT_GEOMETRY_NAMES
    lb_parts = [np.asarray(base.lb, dtype=float), DIRECT_GEOMETRY_LB]
    ub_parts = [np.asarray(base.ub, dtype=float), DIRECT_GEOMETRY_UB]
    x0_parts = [
        np.asarray(base.x0, dtype=float),
        np.array([float(base.P_base.E_loc[1]), float(base.P_base.r_bar[1])]),
    ]

    if closure in POPULATION_SCALE_CLOSURES:
        targets["implied_total_population"] = float(scale_target)
        weights["implied_total_population"] = float(scale_weight)

    if closure in BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES:
        base.P_base.population_closure = closure
        base.P_base.target_city_entry_prob = float(getattr(base.P_base, "target_city_entry_prob", 0.89))
        base.P_base.calibrate_outside_value_to_entry_prob = True
        base.P_base.outside_benchmark_target_total_population = float(scale_target)
    elif closure in RENEWAL_VALVE_CALIBRATED_CLOSURES:
        base.P_base.population_closure = closure
        base.P_base.renewal_calibrate_outside_flow = True
        base.P_base.renewal_target_total_population = float(scale_target)
        base.P_base.renewal_retention = float(renewal_retention)
        base.P_base.outside_entry_shares = np.ones(base.P_base.I) / base.P_base.I
    elif closure in RENEWAL_VALVE_FIXED_CLOSURES:
        names.append(RENEWAL_FLOW_NAME)
        lb_parts.append(np.array([RENEWAL_FLOW_LB]))
        ub_parts.append(np.array([RENEWAL_FLOW_UB]))
        x0_parts.append(np.array([DEFAULT_RENEWAL_FLOW_X0 if outside_flow_x0 is None else float(outside_flow_x0)]))
        base.P_base.population_closure = closure
        base.P_base.renewal_calibrate_outside_flow = False
        base.P_base.outside_entry_flow = float(x0_parts[-1][0])
        base.P_base.renewal_retention = float(renewal_retention)
        base.P_base.outside_entry_shares = np.ones(base.P_base.I) / base.P_base.I
    elif closure in POPULATION_SCALE_CLOSURES:
        names.append(OUTSIDE_VALUE_NAME)
        lb_parts.append(np.array([OUTSIDE_VALUE_LB]))
        ub_parts.append(np.array([OUTSIDE_VALUE_UB]))
        x0_parts.append(np.array([DEFAULT_OUTSIDE_VALUE_X0 if outside_value_x0 is None else float(outside_value_x0)]))
        base.P_base.population_closure = closure
        base.P_base.outside_value = float(x0_parts[-1][0])
        base.P_base.outside_value_is_calibrated = True
        base.P_base.allow_uncalibrated_outside_value = False
    elif closure in OPEN_CITY_CLOSURES:
        base.P_base.population_closure = closure
    else:
        base.P_base.population_closure = closure

    lb = np.concatenate(lb_parts)
    ub = np.concatenate(ub_parts)
    x0 = np.concatenate(x0_parts)

    return DirectCalibrationSetup(
        targets=targets,
        weights=weights,
        P_base=base.P_base,
        names=names,
        lb=lb,
        ub=ub,
        x0=x0,
        setup_mode=str(setup_mode),
        population_closure=closure,
    )


def evaluate_direct_theta(
    theta: Sequence[float],
    setup: DirectCalibrationSetup,
    *,
    max_iter_eq: int | None = None,
    force_full: bool = False,
    eq_penalty_weight: float = 0.0,
    verbose: bool = False,
) -> DirectObjectiveResult:
    """Evaluate one no-inversion SMM objective point."""

    theta_arr = np.asarray(theta, dtype=float).reshape(-1)
    t0 = time.perf_counter()
    if len(theta_arr) != len(setup.names):
        return _failed_result(theta_arr, t0, f"theta has {len(theta_arr)} entries; expected {len(setup.names)}")
    if np.any(theta_arr < setup.lb) or np.any(theta_arr > setup.ub) or np.any(~np.isfinite(theta_arr)):
        return _failed_result(theta_arr, t0, "theta outside direct calibration bounds")

    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta_arr)}
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    structural_theta = [theta_dict[name] for name in structural_names]
    P = apply_theta(P, structural_theta, structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if OUTSIDE_VALUE_NAME in theta_dict:
        P.outside_value = theta_dict[OUTSIDE_VALUE_NAME]
        P.outside_value_is_calibrated = True
        P.allow_uncalibrated_outside_value = False
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    if max_iter_eq is not None:
        P.max_iter_eq = int(max_iter_eq)
    if force_full:
        P.force_full_bellman = True

    try:
        sol, P_out, p_eq = run_model_cp_dt(P, verbose=verbose)
    except Exception as exc:  # pragma: no cover - exercised by cluster jobs
        return _failed_result(theta_arr, t0, repr(exc))

    rent_ratio = float((P_out.user_cost_rate * p_eq[1]) / (P_out.user_cost_rate * p_eq[0]))
    moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
    moments_ns.inv_pop_share_C = float(sol.pop_share[1])
    moments_ns.inv_rent_ratio_C_over_P = rent_ratio
    moments_ns.pop_share_C = float(sol.pop_share[1])
    moments_ns.rent_ratio_C_over_P = rent_ratio
    moments_ns.E_C = float(P_out.E_loc[1])
    moments_ns.r_bar_C = float(P_out.r_bar[1])
    if OUTSIDE_VALUE_NAME in theta_dict:
        moments_ns.outside_value = float(P_out.outside_value)
    if RENEWAL_FLOW_NAME in theta_dict:
        moments_ns.outside_entry_flow = float(P_out.outside_entry_flow)
        moments_ns.renewal_retention = float(getattr(P_out, "renewal_retention", np.nan))
    _attach_population_scale_moments(moments_ns, sol, P_out)

    moments = _namespace_to_float_dict(moments_ns)
    if setup.population_closure in RENEWAL_VALVE_CLOSURES or setup.population_closure in BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES:
        scale = getattr(sol, "accounting_scale", None)
        scale_error = _population_scale_error(scale)
        if scale_error is not None:
            return _failed_result(theta_arr, t0, scale_error)

    loss = compute_smm_loss(moments, setup.targets, setup.weights)
    timings = _jsonable_mapping(getattr(sol, "timings", {}))
    converged = bool(timings.get("accepted", False))

    if eq_penalty_weight > 0:
        best_eq_error = float(timings.get("best_eq_error", np.inf))
        if np.isfinite(best_eq_error):
            excess = max(0.0, best_eq_error - 0.02)
            loss += float(eq_penalty_weight) * (excess / 0.02) ** 2
        else:
            loss = 1e6

    if not np.isfinite(loss):
        loss = 1e6

    return DirectObjectiveResult(
        loss=float(loss),
        theta=theta_arr.copy(),
        moments=moments,
        elapsed_sec=float(time.perf_counter() - t0),
        solve_ok=True,
        converged=converged,
        p_eq=[float(x) for x in np.asarray(p_eq).reshape(-1)],
        timings=timings,
    )


def compute_smm_loss(moments: Mapping[str, float], targets: Mapping[str, float], weights: Mapping[str, float]) -> float:
    """Compute weighted relative SMM loss over all requested targets."""

    loss = 0.0
    for name, target_val in targets.items():
        if name not in weights:
            continue
        model_val = float(moments.get(name, np.nan))
        if not np.isfinite(model_val):
            return 1e6
        target_val = float(target_val)
        dev = (model_val - target_val) / abs(target_val) if abs(target_val) > 1e-6 else model_val - target_val
        loss += float(weights[name]) * dev**2
    return float(loss)


def direct_theta_record(theta: Sequence[float], names: Sequence[str]) -> dict[str, float]:
    """Named parameter dictionary for JSON/CSV output."""

    theta_arr = np.asarray(theta, dtype=float).reshape(-1)
    return {name: float(theta_arr[k]) for k, name in enumerate(names)}


def _failed_result(theta: np.ndarray, t0: float, error: str) -> DirectObjectiveResult:
    return DirectObjectiveResult(
        loss=1e6,
        theta=np.asarray(theta, dtype=float).reshape(-1).copy(),
        moments={},
        elapsed_sec=float(time.perf_counter() - t0),
        solve_ok=False,
        converged=False,
        p_eq=[],
        timings={},
        error=error,
    )


def _population_scale_error(scale: Any) -> str | None:
    """Return a failure reason when the population closure has no finite scale."""

    if scale is None:
        return "renewal scale diagnostics missing"
    denominator = float(getattr(scale, "denominator", np.nan))
    outside_flow = float(getattr(scale, "outside_entry_flow", np.nan))
    implied_population = float(getattr(scale, "implied_total_population", np.nan))
    scale_factor = float(getattr(scale, "scale_factor", np.nan))
    finite_flag = bool(getattr(scale, "finite_stationary_scale", False))
    if (
        not finite_flag
        or not np.isfinite(denominator)
        or denominator <= 1e-12
        or not np.isfinite(outside_flow)
        or outside_flow < 0.0
        or not np.isfinite(implied_population)
        or not np.isfinite(scale_factor)
    ):
        return (
            "invalid population scale: "
            f"finite={finite_flag}, denominator={denominator:.6g}, "
            f"outside_entry_flow={outside_flow:.6g}, "
            f"implied_total_population={implied_population:.6g}, "
            f"scale_factor={scale_factor:.6g}"
        )
    return None


def _namespace_to_float_dict(obj: SimpleNamespace) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in vars(obj).items():
        arr = np.asarray(value)
        if arr.size == 1:
            out[key] = float(arr.reshape(-1)[0])
    return out


def _attach_population_scale_moments(moments: SimpleNamespace, sol: SimpleNamespace, P_out: SimpleNamespace) -> None:
    scale = getattr(sol, "accounting_scale", None)
    if scale is not None:
        moments.scale_factor = float(scale.scale_factor)
        moments.implied_total_population = float(scale.implied_total_population)
        moments.implied_entry_total = float(scale.implied_entry_total)
        moments.city_entry_prob_total = float(getattr(scale, "city_entry_prob_total", np.nan))
        moments.outside_entry_prob = float(getattr(scale, "outside_entry_prob", np.nan))
        moments.mature_cityborn_per_entry = float(scale.mature_cityborn_per_entry)
        moments.entry_per_unit_scale = float(getattr(scale, "entry_per_unit_scale", np.nan))
        moments.mature_cityborn_per_unit_scale = float(getattr(scale, "mature_cityborn_per_unit_scale", np.nan))
        moments.outside_entry_flow = float(getattr(scale, "outside_entry_flow", np.nan))
        moments.renewal_retention = float(getattr(scale, "renewal_retention", np.nan))
        moments.population_scale_denominator = float(scale.denominator)
        moments.stationary_entry_relative_residual = float(scale.stationary_entry_relative_residual)
        moments.implied_housing_demand_P = float(scale.implied_housing_demand[0])
        moments.implied_housing_demand_C = float(scale.implied_housing_demand[1])
        moments.reference_housing_demand_P = float(scale.reference_housing_demand[0])
        moments.reference_housing_demand_C = float(scale.reference_housing_demand[1])
        return

    total_mass = float(getattr(sol, "total_mass", np.nan))
    moments.scale_factor = total_mass / max(float(getattr(P_out, "N_target", 1.0)), 1e-12)
    moments.implied_total_population = total_mass
    moments.implied_entry_total = float(getattr(sol, "entry_rate", np.nan))


def _jsonable_mapping(d: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in d.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        elif isinstance(value, np.generic):
            out[key] = value.item()
        else:
            out[key] = value
    return out
