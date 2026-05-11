"""Regression checks for the population closures.

This script is intentionally small and deterministic. It verifies that the
default normalized solver path is unchanged, that the accounting-scale equation
closes in residuals, and that the renewal-valve closure has the intended
stationary-scale accounting.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import accounting_population_scale, renewal_population_scale, run_model_cp_dt
from dt_cp_model.utils import make_grid


BENCHMARK_DIR = ROOT / "benchmarks"
PRECLOSURE_REFERENCE = BENCHMARK_DIR / "fast_ge_iter1_2026_05_05_preclosure.json"

SCALAR_KEYS = [
    "tfr",
    "own_rate",
    "mean_age_first_birth",
    "migration_rate_2245",
    "prime_childless_renter_median_rooms",
    "prime_childless_owner_median_rooms",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "total_mass",
]


def main() -> None:
    sol, P, p_eq = solve_one_iter_fast()
    report = compact_report(sol, p_eq)
    check_default_path_unchanged(report)

    b_grid = make_grid(P)
    baseline_scale = accounting_population_scale(sol, P, b_grid, calibrate_outside_value=True)
    check_scale_residuals(baseline_scale, label="baseline calibrated scale")
    assert_close(baseline_scale.scale_factor, 1.0, 1e-10, "calibrated baseline scale")

    baseline_renewal = calibrated_renewal_scale(sol, P, b_grid)
    check_scale_residuals(baseline_renewal, label="baseline renewal scale")
    assert_close(baseline_renewal.scale_factor, 1.0, 1e-10, "calibrated baseline renewal scale")

    check_uncalibrated_scaled_closure_rejected()
    check_fertility_scale_sensitivity(baseline_scale)
    check_renewal_fertility_scale_sensitivity(baseline_renewal)
    check_scaled_housing_smoke(baseline_scale)
    check_renewal_scaled_housing_smoke(baseline_renewal)

    print("population closure checks passed")
    print(f"  calibrated outside_value = {baseline_scale.outside_value:.12g}")
    print(f"  mature children per entrant = {baseline_scale.mature_cityborn_per_entry:.12g}")


def solve_one_iter_fast(psi_child: float | None = None) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    setup = build_calibration_setup("fast")
    P = setup.P_base
    P.max_iter_eq = 1
    P.force_full_bellman = True
    if psi_child is not None:
        P.psi_child = float(psi_child)
    return run_model_cp_dt(P, verbose=False)


def compact_report(sol: SimpleNamespace, p_eq: np.ndarray) -> dict[str, object]:
    return {
        "p_eq": [float(x) for x in p_eq],
        "tfr": float(2 * sol.mean_parity),
        "own_rate": float(sol.own_rate),
        "pop_share": [float(x) for x in sol.pop_share],
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "migration_rate_2245": float(getattr(sol, "migration_rate_2245", np.nan)),
        "prime_childless_renter_median_rooms": float(getattr(sol, "prime_childless_renter_median_rooms", np.nan)),
        "prime_childless_owner_median_rooms": float(getattr(sol, "prime_childless_owner_median_rooms", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(getattr(sol, "housing_increment_1to2_proxy_t3", np.nan)),
        "total_mass": float(getattr(sol, "total_mass", np.nan)),
    }


def check_default_path_unchanged(report: dict[str, object]) -> None:
    reference = json.loads(PRECLOSURE_REFERENCE.read_text())
    max_abs = max_abs_report_diff(report, reference)
    if max_abs > 1e-9:
        raise AssertionError(f"default normalized path changed: max_abs_diff={max_abs:.6g}")
    print(f"  default normalized path unchanged: max_abs_diff={max_abs:.3g}")


def check_scale_residuals(scale: SimpleNamespace, *, label: str) -> None:
    if not scale.finite_stationary_scale:
        raise AssertionError(f"{label}: stationary scale is not finite")
    if scale.denominator <= 0:
        raise AssertionError(f"{label}: non-positive denominator {scale.denominator:.6g}")
    if not np.isfinite(scale.scale_factor):
        raise AssertionError(f"{label}: non-finite scale factor")
    if abs(scale.stationary_entry_residual) > 1e-12:
        raise AssertionError(f"{label}: entry residual {scale.stationary_entry_residual:.6g}")
    if abs(scale.stationary_scale_residual) > 1e-12:
        raise AssertionError(f"{label}: scale residual {scale.stationary_scale_residual:.6g}")
    if scale.stationary_entry_relative_residual > 1e-10:
        raise AssertionError(f"{label}: relative residual {scale.stationary_entry_relative_residual:.6g}")


def calibrated_renewal_scale(sol: SimpleNamespace, P: SimpleNamespace, b_grid: np.ndarray) -> SimpleNamespace:
    raw = renewal_population_scale(sol, P, b_grid, outside_entry_flow=1.0)
    outside_flow = raw.entry_per_unit_scale - raw.renewal_retention * raw.mature_cityborn_per_unit_scale
    if outside_flow <= 0:
        raise AssertionError(f"baseline renewal flow is non-positive: {outside_flow:.6g}")
    return renewal_population_scale(sol, P, b_grid, outside_entry_flow=outside_flow)


def check_uncalibrated_scaled_closure_rejected() -> None:
    setup = build_calibration_setup("fast")
    P = setup.P_base
    P.max_iter_eq = 1
    P.population_closure = "accounting_scale_prices"
    try:
        run_model_cp_dt(P, verbose=False)
    except ValueError as exc:
        if "requires a calibrated outside value" not in str(exc):
            raise
    else:
        raise AssertionError("uncalibrated accounting_scale_prices closure was accepted")
    print("  uncalibrated scaled-housing closure rejected")


def check_fertility_scale_sensitivity(baseline_scale: SimpleNamespace) -> None:
    low_sol, low_P, _ = solve_one_iter_fast(psi_child=0.03)
    high_sol, high_P, _ = solve_one_iter_fast(psi_child=0.12)
    low_scale = accounting_population_scale(
        low_sol,
        low_P,
        make_grid(low_P),
        outside_value=baseline_scale.outside_value,
        kappa_entry=baseline_scale.kappa_entry,
    )
    high_scale = accounting_population_scale(
        high_sol,
        high_P,
        make_grid(high_P),
        outside_value=baseline_scale.outside_value,
        kappa_entry=baseline_scale.kappa_entry,
    )
    check_scale_residuals(low_scale, label="low-fertility scale")
    check_scale_residuals(high_scale, label="high-fertility scale")
    if not high_scale.scale_factor > low_scale.scale_factor:
        raise AssertionError(
            "fertility sensitivity failed: "
            f"low_scale={low_scale.scale_factor:.6g}, high_scale={high_scale.scale_factor:.6g}"
        )
    print(
        "  fixed-outside-value fertility sensitivity: "
        f"scale {low_scale.scale_factor:.6g} -> {high_scale.scale_factor:.6g}"
    )


def check_renewal_fertility_scale_sensitivity(baseline_scale: SimpleNamespace) -> None:
    low_sol, low_P, _ = solve_one_iter_fast(psi_child=0.03)
    high_sol, high_P, _ = solve_one_iter_fast(psi_child=0.12)
    low_scale = renewal_population_scale(
        low_sol,
        low_P,
        make_grid(low_P),
        outside_entry_flow=baseline_scale.outside_entry_flow,
        renewal_retention=baseline_scale.renewal_retention,
    )
    high_scale = renewal_population_scale(
        high_sol,
        high_P,
        make_grid(high_P),
        outside_entry_flow=baseline_scale.outside_entry_flow,
        renewal_retention=baseline_scale.renewal_retention,
    )
    check_scale_residuals(low_scale, label="low-fertility renewal scale")
    check_scale_residuals(high_scale, label="high-fertility renewal scale")
    if not high_scale.scale_factor > low_scale.scale_factor:
        raise AssertionError(
            "renewal fertility sensitivity failed: "
            f"low_scale={low_scale.scale_factor:.6g}, high_scale={high_scale.scale_factor:.6g}"
        )
    print(
        "  renewal-valve fertility sensitivity: "
        f"scale {low_scale.scale_factor:.6g} -> {high_scale.scale_factor:.6g}"
    )


def check_scaled_housing_smoke(baseline_scale: SimpleNamespace) -> None:
    setup = build_calibration_setup("fast")
    P = setup.P_base
    P.max_iter_eq = 3
    P.force_full_bellman = True
    P.population_closure = "accounting_scale_prices"
    P.outside_value = baseline_scale.outside_value
    P.outside_value_is_calibrated = True
    P.kappa_entry = baseline_scale.kappa_entry
    P.collect_ge_trace = True
    sol, _, _ = run_model_cp_dt(P, verbose=False)
    if not hasattr(sol, "accounting_scale"):
        raise AssertionError("scaled-housing smoke did not attach accounting_scale")
    check_scale_residuals(sol.accounting_scale, label="scaled-housing smoke")
    timings = getattr(sol, "timings", {})
    if timings.get("population_closure") != "accounting_scale_prices":
        raise AssertionError("scaled-housing smoke did not report the expected closure")
    iterations = int(timings.get("iterations_completed", 0))
    if iterations < 1 or iterations > 3:
        raise AssertionError(f"scaled-housing smoke iteration count changed: {timings}")
    print(
        "  scaled-housing smoke finite: "
        f"scale={sol.accounting_scale.scale_factor:.6g}, "
        f"reason={timings.get('convergence_reason')}"
    )


def check_renewal_scaled_housing_smoke(baseline_scale: SimpleNamespace) -> None:
    setup = build_calibration_setup("fast")
    P = setup.P_base
    P.max_iter_eq = 3
    P.force_full_bellman = True
    P.population_closure = "renewal_valve"
    P.outside_entry_flow = baseline_scale.outside_entry_flow
    P.renewal_retention = baseline_scale.renewal_retention
    P.collect_ge_trace = True
    sol, _, _ = run_model_cp_dt(P, verbose=False)
    if not hasattr(sol, "accounting_scale"):
        raise AssertionError("renewal scaled-housing smoke did not attach accounting_scale")
    check_scale_residuals(sol.accounting_scale, label="renewal scaled-housing smoke")
    timings = getattr(sol, "timings", {})
    if timings.get("population_closure") != "renewal_valve":
        raise AssertionError("renewal scaled-housing smoke did not report the expected closure")
    iterations = int(timings.get("iterations_completed", 0))
    if iterations < 1 or iterations > 3:
        raise AssertionError(f"renewal scaled-housing smoke iteration count changed: {timings}")
    print(
        "  renewal scaled-housing smoke finite: "
        f"scale={sol.accounting_scale.scale_factor:.6g}, "
        f"reason={timings.get('convergence_reason')}"
    )


def max_abs_report_diff(a: dict[str, object], b: dict[str, object]) -> float:
    diffs = [abs(float(a[key]) - float(b[key])) for key in SCALAR_KEYS]
    diffs.extend(abs(float(x) - float(y)) for x, y in zip(a["p_eq"], b["p_eq"]))
    diffs.extend(abs(float(x) - float(y)) for x, y in zip(a["pop_share"], b["pop_share"]))
    return max(diffs)


def assert_close(value: float, target: float, tol: float, label: str) -> None:
    if abs(float(value) - float(target)) > tol:
        raise AssertionError(f"{label}: {value:.12g} differs from {target:.12g} by more than {tol:.1e}")


if __name__ == "__main__":
    main()
