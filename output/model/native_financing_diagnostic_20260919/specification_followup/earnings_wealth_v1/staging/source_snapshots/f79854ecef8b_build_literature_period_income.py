"""Pure Sommer (2016) annual-to-period income-process adapter.

The adapter is deliberately independent of the household solver.  It maps the
annual persistent-plus-iid log process into a mean-one four-year period proxy,
then constructs a finite Rouwenhorst/Gauss-Hermite approximation directly at
the fitted period parameters.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def _validate_positive(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _validate_period_years(value: int | float) -> int:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 1.0 or not numeric.is_integer():
        raise ValueError("period_years must be a finite positive integer")
    return int(numeric)


def _stationary_weights(n: int) -> np.ndarray:
    return np.array([math.comb(n - 1, k) for k in range(n)], dtype=float) / 2.0 ** (n - 1)


def _rouwenhorst_period(n: int, rho: float, stationary_variance: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rouwenhorst chain with period-rho and period innovation variance."""
    n = int(n)
    if n < 2:
        raise ValueError("n_persistent must be at least 2")
    rho = float(rho)
    if not np.isfinite(rho) or not 0.0 < rho < 1.0:
        raise ValueError("rho_period must lie strictly between zero and one")
    stationary_variance = _validate_positive("stationary persistent log variance", stationary_variance)
    p = (1.0 + rho) / 2.0
    pi = np.array([[p, 1.0 - p], [1.0 - p, p]], dtype=float)
    for size in range(3, n + 1):
        old = pi
        pi = np.zeros((size, size), dtype=float)
        pi[:-1, :-1] += p * old
        pi[:-1, 1:] += (1.0 - p) * old
        pi[1:, :-1] += (1.0 - p) * old
        pi[1:, 1:] += p * old
        pi[1:-1] *= 0.5
    pi /= pi.sum(axis=1, keepdims=True)
    weights = _stationary_weights(n)
    log_grid = np.linspace(-math.sqrt(stationary_variance) * math.sqrt(n - 1),
                           math.sqrt(stationary_variance) * math.sqrt(n - 1), n)
    levels = np.exp(log_grid)
    levels /= float(weights @ levels)
    return levels, weights, pi


def _iid_lognormal_rule(log_sd: float, n: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Mean-one Gauss-Hermite approximation to exp(N(-v/2,v))."""
    log_sd = float(log_sd)
    n = int(n)
    if not np.isfinite(log_sd) or log_sd < 0.0:
        raise ValueError("iid log standard deviation must be finite and nonnegative")
    if n < 2:
        raise ValueError("n_iid must be at least 2")
    nodes, weights = np.polynomial.hermite.hermgauss(n)
    weights = weights / math.sqrt(math.pi)
    levels = np.exp(-0.5 * log_sd**2 + math.sqrt(2.0) * log_sd * nodes)
    levels /= float(weights @ levels)
    return levels, weights


def annual_block_level_covariances(
    rho_annual: float,
    persistent_innovation_sd_annual: float,
    transitory_sd_annual: float,
    period_years: int = 4,
    max_lag: int = 4,
) -> np.ndarray:
    """Exact stationary level covariances of block-average annual earnings."""
    rho_annual = float(rho_annual)
    if not np.isfinite(rho_annual) or not 0.0 < rho_annual < 1.0:
        raise ValueError("rho_annual must lie strictly between zero and one")
    eta = _validate_positive("persistent annual innovation SD", persistent_innovation_sd_annual)
    eps = _validate_positive("transitory annual SD", transitory_sd_annual)
    period_years = _validate_period_years(period_years)
    max_lag = int(max_lag)
    if period_years < 1 or max_lag < 0:
        raise ValueError("period_years and max_lag must be nonnegative integers (period_years positive)")
    vp = eta**2 / (1.0 - rho_annual**2)
    ve = eps**2
    out = np.empty(max_lag + 1, dtype=float)
    for lag in range(max_lag + 1):
        cov = 0.0
        for i in range(period_years):
            for j in range(period_years):
                annual_lag = abs(lag * period_years + j - i)
                log_cov = vp + ve if annual_lag == 0 else vp * rho_annual**annual_lag
                cov += math.expm1(log_cov)
        out[lag] = cov / period_years**2
    return out


def _chain_level_covariances(z: np.ndarray, weights: np.ndarray, transition: np.ndarray, max_lag: int = 4) -> np.ndarray:
    mean = float(weights @ z)
    centered = z - mean
    out = np.empty(max_lag + 1, dtype=float)
    matrix_power = np.eye(z.size)
    for lag in range(max_lag + 1):
        out[lag] = float((weights * centered) @ matrix_power @ centered)
        matrix_power = matrix_power @ transition
    return out


def _chain_log_covariances(z: np.ndarray, weights: np.ndarray, transition: np.ndarray, max_lag: int = 4) -> np.ndarray:
    """Autocovariances of the actual discrete log nodes."""
    return _chain_level_covariances(np.log(np.asarray(z, dtype=float)), weights, transition, max_lag)


def _fit_period_proxy(covariances: np.ndarray) -> tuple[float, float, float]:
    c0, c1, c2 = map(float, covariances[:3])
    if min(c0, c1, c2) <= -1.0:
        raise ValueError("level covariances must exceed -1")
    den = math.log1p(c1)
    num = math.log1p(c2)
    if den == 0.0:
        raise ValueError("period proxy persistence is undefined when C1 is zero")
    rho = num / den
    vp = den / rho if rho != 0.0 else float("nan")
    ve = math.log1p(c0) - vp
    if not (np.isfinite(rho) and 0.0 < rho < 1.0 and np.isfinite(vp) and vp > 0.0 and np.isfinite(ve) and ve >= 0.0):
        raise ValueError(f"inadmissible period proxy: rho_period={rho}, Vp={vp}, Ve={ve}")
    return rho, vp, ve


def _build_grid(n_persistent: int, iid_states: int, rho: float, vp: float, ve: float) -> tuple[dict[str, Any], dict[str, Any]]:
    persistent, wp, pi_p = _rouwenhorst_period(n_persistent, rho, vp)
    iid, we = _iid_lognormal_rule(math.sqrt(ve), iid_states)
    z = np.multiply.outer(persistent, iid).reshape(-1)
    w = np.multiply.outer(wp, we).reshape(-1)
    pi_eps = np.broadcast_to(we, (we.size, we.size)).copy()
    pi = np.kron(pi_p, pi_eps)
    w /= float(w.sum())
    z /= float(w @ z)
    override = {
        "use_income_types": True,
        "income_type_transition": "markov",
        "income_shock_persistence": float(rho),
        "z_grid": z,
        "z_weights": w,
        "Pi_z": pi,
        "permanent_income_levels_enabled": False,
        "permanent_income_log_variance": 0.0,
    }
    return override, {"persistent_levels": persistent, "persistent_weights": wp, "persistent_transition": pi_p, "iid_levels": iid, "iid_weights": we}


def build_literature_period_income_adapter(
    *,
    rho_annual: float = 0.95,
    persistent_innovation_sd_annual: float = 0.21,
    transitory_sd_annual: float = 0.17,
    period_years: int = 4,
    n_persistent: int = 5,
    iid_states: int = 3,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build deterministic overrides and an auditable moment receipt."""
    annual_cov = annual_block_level_covariances(rho_annual, persistent_innovation_sd_annual, transitory_sd_annual, period_years, 4)
    rho, vp, ve = _fit_period_proxy(annual_cov)
    overrides, components = _build_grid(n_persistent, iid_states, rho, vp, ve)
    actual_cov = _chain_level_covariances(overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"], 4)
    actual_log_cov = _chain_log_covariances(overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"], 4)
    target_log_cov = np.array([vp + ve] + [vp * rho**k for k in range(1, 5)])
    proxy_cov = np.expm1(target_log_cov)
    grid_table = []
    for n in (5, 9, 15):
        small, _ = _build_grid(n, iid_states, rho, vp, ve)
        got = _chain_level_covariances(small["z_grid"], small["z_weights"], small["Pi_z"], 4)
        grid_table.append({"n_persistent": n, "joint_states": int(n * iid_states), "level_covariances": got.tolist(), "errors_vs_exact_block": (got - annual_cov).tolist(), "mean": float(small["z_weights"] @ small["z_grid"])})
    metadata = {
        "referenceannualparams": {"rho": float(rho_annual), "persistent_innovation_sd": float(persistent_innovation_sd_annual), "transitory_sd": float(transitory_sd_annual), "permanent_type": False},
        "period_years": int(period_years), "n_persistent": int(n_persistent), "iid_states": int(iid_states),
        "continuous_covariance_targets": annual_cov.tolist(),
        "proxy_fit": {"rho_period": rho, "Vp": vp, "Ve": ve, "persistent_innovation_variance_period": vp * (1.0 - rho**2), "proxy_covariances": proxy_cov.tolist(), "fit_errors": (proxy_cov - annual_cov).tolist()},
        "actual_chain_level_covariances": actual_cov.tolist(), "actual_chain_log_covariances": actual_log_cov.tolist(),
        "actual_chain_mean": float(overrides["z_weights"] @ overrides["z_grid"]), "actual_chain_probabilities": overrides["z_weights"].tolist(),
        "grid_resolution_table": grid_table,
        "persistent_component": {k: np.asarray(v).tolist() for k, v in components.items()},
        "permanent_income_levels_enabled": False,
    }
    return overrides, metadata


def conventional_endpoint_comparison(
    *, rho_annual: float = 0.95, persistent_innovation_sd_annual: float = 0.21,
    transitory_sd_annual: float = 0.17, period_years: int = 4, max_lag: int = 4,
) -> dict[str, Any]:
    """Diagnostic endpoint conversion, retained separately from exact block matching."""
    period_years = _validate_period_years(period_years)
    rho = float(rho_annual) ** period_years
    vp = float(persistent_innovation_sd_annual) ** 2 / (1.0 - float(rho_annual) ** 2)
    ve = math.log1p(math.expm1(float(transitory_sd_annual) ** 2) / float(period_years))
    log_cov = np.array([vp + ve] + [vp * rho**k for k in range(1, int(max_lag) + 1)])
    exact = annual_block_level_covariances(rho_annual, persistent_innovation_sd_annual, transitory_sd_annual, period_years, max_lag)
    grids = []
    for n in (5, 9, 15):
        grid, _ = _build_grid(n, 3, rho, vp, ve)
        actual = _chain_level_covariances(grid["z_grid"], grid["z_weights"], grid["Pi_z"], max_lag)
        grids.append({"n_persistent": n, "joint_states": 3 * n, "actual_level_covariances": actual.tolist(), "errors_vs_exact_block": (actual - exact).tolist()})
    return {"label": "conventional_endpoint_diagnostic_only", "rho_period": rho, "Vp": vp, "Ve": ve,
            "continuous_level_covariances": np.expm1(log_cov).tolist(),
            "continuous_log_covariances": log_cov.tolist(), "exact_annual_block_covariances": exact.tolist(),
            "grid_resolution_table": grids}


def build_conventional_endpoint_income_adapter(
    *, rho_annual: float = 0.95, persistent_innovation_sd_annual: float = 0.21,
    transitory_sd_annual: float = 0.17, period_years: int = 4,
    n_persistent: int = 5, iid_states: int = 3,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the conventional endpoint diagnostic as an explicit opt-in object."""
    period_years = _validate_period_years(period_years)
    rho_annual = float(rho_annual)
    eta = _validate_positive("persistent annual innovation SD", persistent_innovation_sd_annual)
    eps = _validate_positive("transitory annual SD", transitory_sd_annual)
    if not np.isfinite(rho_annual) or not 0.0 < rho_annual < 1.0:
        raise ValueError("rho_annual must lie strictly between zero and one")
    n_persistent, iid_states = int(n_persistent), int(iid_states)
    if n_persistent < 2 or iid_states < 2:
        raise ValueError("n_persistent and iid_states must be at least two")
    rho = rho_annual**period_years
    vp = eta**2 / (1.0 - rho_annual**2)
    ve = math.log1p(math.expm1(eps**2) / period_years)
    exact = annual_block_level_covariances(rho_annual, eta, eps, period_years, 4)
    continuous_log = np.array([vp + ve] + [vp * rho**k for k in range(1, 5)])
    continuous = np.expm1(continuous_log)
    overrides, components = _build_grid(n_persistent, iid_states, rho, vp, ve)
    discrete_level = _chain_level_covariances(overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"], 4)
    discrete_log = _chain_log_covariances(overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"], 4)
    metadata = {
        "label": "conventional_endpoint_diagnostic_only", "referenceannualparams": {"rho": rho_annual, "persistent_innovation_sd": eta, "transitory_sd": eps, "permanent_type": False},
        "period_years": period_years, "n_persistent": n_persistent, "iid_states": iid_states,
        "rho_period": rho, "Vp": vp, "Ve": ve,
        "continuous_endpoint_level_covariances": continuous.tolist(), "continuous_endpoint_log_covariances": continuous_log.tolist(),
        "exact_annual_block_level_covariances": exact.tolist(),
        "continuous_level_errors_vs_exact_block": (continuous - exact).tolist(),
        "discrete_level_covariances": discrete_level.tolist(), "discrete_log_covariances": discrete_log.tolist(),
        "discrete_level_errors_vs_continuous_endpoint": (discrete_level - continuous).tolist(),
        "discrete_level_errors_vs_exact_block": (discrete_level - exact).tolist(),
        "components": {k: np.asarray(v).tolist() for k, v in components.items()},
        "iid_transition_independent": bool(np.allclose(overrides["Pi_z"].reshape(n_persistent, iid_states, n_persistent, iid_states)[:, 0, :, :], overrides["Pi_z"].reshape(n_persistent, iid_states, n_persistent, iid_states)[:, -1, :, :])),
        "stationary": bool(np.allclose(overrides["z_weights"] @ overrides["Pi_z"], overrides["z_weights"], atol=1e-12)),
    }
    return overrides, metadata


def diagnose_literature_period_income_adapter(
    *, rho_annual: float = 0.95, persistent_innovation_sd_annual: float = 0.21,
    transitory_sd_annual: float = 0.17, period_years: int = 4,
) -> dict[str, Any]:
    """Return the pinned continuous diagnostic, including an admissibility blocker."""
    cov = annual_block_level_covariances(rho_annual, persistent_innovation_sd_annual,
                                         transitory_sd_annual, period_years, 4)
    endpoint = conventional_endpoint_comparison(rho_annual=rho_annual, persistent_innovation_sd_annual=persistent_innovation_sd_annual, transitory_sd_annual=transitory_sd_annual, period_years=period_years)
    try:
        rho, vp, ve = _fit_period_proxy(cov)
        return {"status": "admissible", "continuous_covariance_targets": cov.tolist(), "conventional_endpoint": endpoint,
                "proxy_fit": {"rho_period": rho, "Vp": vp, "Ve": ve}}
    except ValueError as exc:
        den = math.log1p(float(cov[1])); num = math.log1p(float(cov[2]))
        rho = num / den if den else float("nan")
        vp = den / rho if rho else float("nan")
        ve = math.log1p(float(cov[0])) - vp if np.isfinite(vp) else float("nan")
        return {"status": "blocked_inadmissible_proxy", "continuous_covariance_targets": cov.tolist(), "conventional_endpoint": endpoint,
                "proxy_fit": {"rho_period": rho, "Vp": vp, "Ve": ve}, "blocker": str(exc)}


def write_receipt(path: str | Path, overrides: dict[str, Any], metadata: dict[str, Any]) -> None:
    payload = {"metadata": metadata, "overrides": {k: (np.asarray(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in overrides.items()}}
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


__all__ = ["annual_block_level_covariances", "build_literature_period_income_adapter", "build_conventional_endpoint_income_adapter", "conventional_endpoint_comparison", "diagnose_literature_period_income_adapter", "write_receipt"]
