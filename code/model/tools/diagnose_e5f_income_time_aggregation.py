"""Diagnose annual-to-four-year aggregation for the income candidate.

This is a data-free Monte Carlo diagnostic.  It does not estimate parameters,
solve the household model, or modify a calibration contract.  The retained
annual process is lognormal with stationary persistent log variance ``Vp``,
annual persistence ``rho``, and iid transitory log variance ``Ve``.  Annual
levels are normalized to have mean one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CANDIDATE = ROOT / "output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json"
DEFAULT_OUTPUT = ROOT / "output/model/native_financing_diagnostic_20260919/earnings_time_aggregation"
DEFAULT_SEED = 20260920
DEFAULT_BATCHES = {"smoke": 2, "full": 20}
DEFAULT_HOUSEHOLDS = {"smoke": 2_000, "full": 20_000}
DEFAULT_YEARS = {"smoke": 40, "full": 120}
TOTAL_TIME_BUDGET_SECONDS = {"smoke": 60.0, "full": 720.0}
REPORT_LAGS = (0, 1, 2, 4)
QUANTILE_LEVELS = (0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
MC_SIGMA_TOLERANCE = 6.0
MC_ABSOLUTE_FLOOR = 1e-3


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_fingerprint(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def dependency_paths(candidate_path: Path) -> dict[str, Path]:
    source_dir = ROOT / "code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727"
    paths = {
        "driver": Path(__file__).resolve(),
        "candidate_json": candidate_path.resolve(),
        "constructor": (ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py").resolve(),
        "local_panel": (ROOT / "code/model/intergen_eqscale_seq_optimized/local_panel.py").resolve(),
        "source_autocovariance_fit": (source_dir / "md_autocovariance_fit.csv").resolve(),
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"missing diagnostic dependency: {path}")
    return paths


def dependency_hashes(paths: dict[str, Path]) -> dict[str, str]:
    return {name: sha256(path) for name, path in paths.items()}


def prepare_output_directory(output: Path) -> None:
    """Reserve a fresh output directory; never overwrite an existing packet."""
    if output.exists():
        if not output.is_dir() or any(output.iterdir()):
            raise FileExistsError(f"refusing to overwrite nonempty output directory: {output}")
    else:
        output.mkdir(parents=True)


def final_status(
    *,
    requested_batches: int,
    completed_batches: int,
    exact_checks_pass: bool,
    mean_checks_pass: bool,
    plots_present: bool,
    dependencies_unchanged: bool,
    stop_reason: str | None,
) -> str:
    if stop_reason == "time_budget":
        return "stopped_time_budget"
    if completed_batches < requested_batches:
        return "failed"
    if not all((exact_checks_pass, mean_checks_pass, plots_present, dependencies_unchanged)):
        return "failed"
    return "completed"


def annual_level_covariance(v_p: float, v_e: float, rho: float, lag: int) -> float:
    """Covariance of mean-one annual levels at an integer lag."""
    distance = abs(int(lag))
    same_period = v_e if distance == 0 else 0.0
    return float(np.exp(v_p * rho**distance + same_period) - 1.0)


def block_average_level_covariance(
    v_p: float, v_e: float, rho: float, block_size: int, block_lag: int
) -> float:
    """Exact covariance of two block averages from the annual process.

    The 16 terms are retained explicitly for the four-year case.  This is a
    level covariance, not a covariance of log block averages.
    """
    size = int(block_size)
    if size < 1:
        raise ValueError("block_size must be positive")
    pair_sum = 0.0
    for i in range(size):
        for j in range(size):
            pair_sum += annual_level_covariance(v_p, v_e, rho, size * int(block_lag) + j - i)
    return float(pair_sum / (size * size))


def continuous_endpoint_moments(
    v_p: float, v_e_period: float, rho_period: float, lags: Iterable[int] = REPORT_LAGS
) -> dict[str, dict[str, float]]:
    """Moments of endpoint persistent lognormal times iid period shock."""
    out: dict[str, dict[str, float]] = {}
    for lag in lags:
        lag_i = int(lag)
        log_cov = v_p + v_e_period if lag_i == 0 else v_p * rho_period**lag_i
        level_cov = np.exp(log_cov) - 1.0 if lag_i == 0 else np.exp(v_p * rho_period**lag_i) - 1.0
        out[str(lag_i)] = {
            "log_covariance": float(log_cov),
            "level_covariance": float(level_cov),
        }
    return {"mean": {"value": 1.0}, "moments": out}


def markov_payload_moments(
    z_grid: np.ndarray,
    z_weights: np.ndarray,
    transition: np.ndarray,
    lags: Iterable[int] = REPORT_LAGS,
) -> dict[str, dict[str, float]]:
    """Compute exact moments implied by a finite level Markov payload."""
    z = np.asarray(z_grid, dtype=float).reshape(-1)
    weights = np.asarray(z_weights, dtype=float).reshape(-1)
    pi = np.asarray(transition, dtype=float)
    if z.size != weights.size or pi.shape != (z.size, z.size):
        raise ValueError("incompatible Markov payload shapes")
    if np.any(z <= 0.0) or np.any(weights < 0.0):
        raise ValueError("Markov payload must have positive levels and nonnegative weights")
    weights = weights / weights.sum()
    if not np.allclose(weights @ pi, weights, atol=1e-10, rtol=0.0):
        raise ValueError("weights are not stationary for the supplied transition")
    mean = float(weights @ z)
    centered_level = z - mean
    logz = np.log(z)
    centered_log = logz - float(weights @ logz)
    out: dict[str, dict[str, float]] = {}
    for lag in lags:
        lag_i = int(lag)
        if lag_i < 0:
            raise ValueError("lags must be nonnegative")
        pi_lag = np.linalg.matrix_power(pi, lag_i)
        out[str(lag_i)] = {
            "log_covariance": float((weights * centered_log) @ pi_lag @ centered_log),
            "level_covariance": float((weights * centered_level) @ pi_lag @ centered_level),
        }
    return {"mean": {"value": mean}, "moments": out}


def load_candidate(candidate_path: Path) -> dict[str, Any]:
    candidate = read_json(candidate_path)
    period_years = float(candidate["period_years"])
    if not np.isclose(period_years, 4.0, atol=0.0, rtol=0.0):
        raise ValueError("this diagnostic's MC loop is hardcoded to four-year blocks")
    annual = candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
    mapping = candidate["four_year_diagnostic_mapping"]
    rho = float(annual["rho_annual"])
    v_p = float(annual["persistent_variance"])
    v_e = float(annual["transitory_variance"])
    rho_period = rho ** period_years
    v_e_period = float(mapping["transitory_log_variance_period"])
    expected = math.log(1.0 + (math.exp(v_e) - 1.0) / period_years)
    if not np.isclose(v_e_period, expected, atol=1e-12, rtol=0.0):
        raise ValueError("candidate transitory mapping does not match the retained formula")
    if not (0.0 < rho < 1.0 and v_p > 0.0 and v_e >= 0.0):
        raise ValueError("candidate annual primitives are invalid")
    return {
        "candidate": candidate,
        "rho_annual": rho,
        "persistent_variance_annual": v_p,
        "transitory_variance_annual": v_e,
        "period_years": int(round(period_years)),
        "rho_period": rho_period,
        "transitory_log_variance_period": v_e_period,
    }


def payload_from_candidate(candidate_path: Path, primitives: dict[str, Any]) -> dict[str, Any]:
    """Build the existing 15-state payload through its production constructor."""
    model_path = ROOT / "code/model"
    if str(model_path) not in sys.path:
        sys.path.insert(0, str(model_path))
    from intergen_eqscale_seq_optimized.local_panel import income_process_fingerprint

    annual_eta_sd = math.sqrt(
        (1.0 - primitives["rho_annual"] ** 2) * primitives["persistent_variance_annual"]
    )
    from build_persistent_transitory_income_candidate import (
        build_persistent_transitory_income_candidate,
    )

    overrides, metadata = build_persistent_transitory_income_candidate(
        rho_annual=primitives["rho_annual"],
        persistent_innovation_sd_annual=annual_eta_sd,
        transitory_log_sd_period=math.sqrt(primitives["transitory_log_variance_period"]),
        period_years=float(primitives["period_years"]),
    )
    payload_fingerprint = income_process_fingerprint(overrides)
    candidate = primitives["candidate"]
    expected_state_count = int(candidate.get("state_count", overrides["z_grid"].size))
    if overrides["z_grid"].size != expected_state_count:
        raise ValueError("rebuilt payload state count disagrees with candidate metadata")
    for field, expected in (("persistent_states", 5), ("transitory_states", 3), ("joint_states", 15)):
        if field in metadata and int(metadata[field]) != expected:
            raise ValueError(f"rebuilt payload metadata has unexpected {field}")
    expected_fingerprints: dict[str, str] = {}
    expected_candidate_hashes: dict[str, str] = {}
    for plan_name in ("calibration_plan.json", "calibration_plan.remote.json", "search_plan.remote.json"):
        plan_path = candidate_path.parent / plan_name
        if plan_path.is_file():
            plan = read_json(plan_path)
            fingerprint = plan.get("candidate_payload_fingerprint")
            if fingerprint:
                expected_fingerprints[plan_name] = str(fingerprint)
            candidate_hash = plan.get("candidate_json_sha256")
            if candidate_hash:
                expected_candidate_hashes[plan_name] = str(candidate_hash)
    candidate_fingerprint = json_fingerprint(candidate)
    if expected_fingerprints and any(candidate_fingerprint != value for value in expected_fingerprints.values()):
        raise ValueError("candidate JSON fingerprint disagrees with candidate plan metadata")
    candidate_sha = sha256(candidate_path)
    if expected_candidate_hashes and any(candidate_sha != value for value in expected_candidate_hashes.values()):
        raise ValueError("candidate JSON SHA-256 disagrees with candidate plan metadata")
    metadata_payload = candidate.get("payload", candidate.get("income_payload", {}))
    for key in ("z_grid", "z_weights", "Pi_z"):
        if key in metadata_payload:
            expected_array = np.asarray(metadata_payload[key], dtype=float)
            if not np.array_equal(expected_array, np.asarray(overrides[key], dtype=float)):
                raise ValueError(f"rebuilt payload array disagrees with candidate metadata: {key}")
    return {
        "overrides": overrides,
        "metadata": metadata,
        "payload_fingerprint": payload_fingerprint,
        "expected_fingerprints": expected_fingerprints,
        "candidate_fingerprint": candidate_fingerprint,
        "expected_candidate_hashes": expected_candidate_hashes,
        "constructor_sha256": sha256(ROOT / "code/model/tools/build_persistent_transitory_income_candidate.py"),
        "candidate_sha256": candidate_sha,
    }


def _batch_simulation(
    rng: np.random.Generator, v_p: float, v_e: float, rho: float, years: int, households: int
) -> dict[str, Any]:
    """Simulate one independent panel and return only summaries."""
    if years % 4:
        raise ValueError("years must be divisible by four")
    p = rng.normal(0.0, math.sqrt(v_p), size=households)
    eta_sd = math.sqrt(v_p * (1.0 - rho * rho))
    annual = np.empty((households, years), dtype=float)
    for year in range(years):
        if year:
            p = rho * p + rng.normal(0.0, eta_sd, size=households)
        e = rng.normal(0.0, math.sqrt(v_e), size=households)
        annual[:, year] = np.exp(p + e - 0.5 * (v_p + v_e))
    blocks = annual.reshape(households, years // 4, 4).mean(axis=2)
    block_log = np.log(blocks)
    annual_mean = float(annual.mean())
    annual_covariances: dict[str, float] = {}
    annual_centered = annual - annual_mean
    for lag in REPORT_LAGS:
        if lag == 0:
            annual_covariances[str(lag)] = float(np.mean(annual_centered * annual_centered))
        else:
            annual_covariances[str(lag)] = float(
                np.mean(annual_centered[:, :-lag] * annual_centered[:, lag:])
            )
    block_mean = float(blocks.mean())
    block_centered = blocks - block_mean
    log_mean = float(block_log.mean())
    log_centered = block_log - log_mean
    level_moments: dict[str, float] = {}
    log_moments: dict[str, float] = {}
    for lag in REPORT_LAGS:
        if lag == 0:
            level_moments[str(lag)] = float(np.mean(block_centered * block_centered))
            log_moments[str(lag)] = float(np.mean(log_centered * log_centered))
        else:
            level_moments[str(lag)] = float(
                np.mean(block_centered[:, :-lag] * block_centered[:, lag:])
            )
            log_moments[str(lag)] = float(
                np.mean(log_centered[:, :-lag] * log_centered[:, lag:])
            )
    quantiles = {str(q): float(x) for q, x in zip(QUANTILE_LEVELS, np.quantile(blocks, QUANTILE_LEVELS))}
    return {
        "annual_mean": annual_mean,
        "annual_level_covariance": annual_covariances,
        "block_mean": block_mean,
        "block_level_covariance": level_moments,
        "block_log_covariance": log_moments,
        "block_quantiles": quantiles,
    }


def _aggregate_batch_statistics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def summarize(values: list[float]) -> dict[str, float]:
        array = np.asarray(values, dtype=float)
        return {
            "estimate": float(array.mean()),
            "mc_se": float(array.std(ddof=1) / math.sqrt(array.size)) if array.size > 1 else 0.0,
            "batches": int(array.size),
        }

    output: dict[str, Any] = {
        "annual_mean": summarize([row["annual_mean"] for row in rows]),
        "annual_level_covariance": {},
        "block_mean": summarize([row["block_mean"] for row in rows]),
        "block_level_covariance": {},
        "block_log_covariance": {},
        "block_quantiles": {},
    }
    for key in map(str, REPORT_LAGS):
        output["annual_level_covariance"][key] = summarize(
            [row["annual_level_covariance"][key] for row in rows]
        )
        output["block_level_covariance"][key] = summarize(
            [row["block_level_covariance"][key] for row in rows]
        )
        output["block_log_covariance"][key] = summarize(
            [row["block_log_covariance"][key] for row in rows]
        )
    for q in map(str, QUANTILE_LEVELS):
        output["block_quantiles"][q] = summarize([row["block_quantiles"][q] for row in rows])
    return output


def _empty_batch_statistics() -> dict[str, Any]:
    empty = {"estimate": None, "mc_se": None, "batches": 0}
    return {
        "annual_mean": dict(empty),
        "annual_level_covariance": {str(k): dict(empty) for k in REPORT_LAGS},
        "block_mean": dict(empty),
        "block_level_covariance": {str(k): dict(empty) for k in REPORT_LAGS},
        "block_log_covariance": {str(k): dict(empty) for k in REPORT_LAGS},
        "block_quantiles": {str(q): dict(empty) for q in QUANTILE_LEVELS},
    }


def _validity_checks(
    monte_carlo: dict[str, Any], exact_block: dict[str, float], exact_annual: dict[str, float]
) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    mean_checks: dict[str, Any] = {}
    for label, estimate in (("annual", monte_carlo["annual_mean"]), ("block", monte_carlo["block_mean"])):
        gap = abs(float(estimate["estimate"]) - 1.0)
        tolerance = MC_SIGMA_TOLERANCE * float(estimate["mc_se"]) + MC_ABSOLUTE_FLOOR
        mean_checks[f"{label}_mean"] = {
            "passed": bool(gap <= tolerance),
            "estimate": float(estimate["estimate"]),
            "target": 1.0,
            "gap": float(gap),
            "tolerance": float(tolerance),
        }
    for label, estimates, targets in (
        ("annual", monte_carlo["annual_level_covariance"], exact_annual),
        ("block", monte_carlo["block_level_covariance"], exact_block),
    ):
        for key, target in targets.items():
            estimate = estimates[str(key)]
            gap = abs(float(estimate["estimate"]) - float(target))
            tolerance = MC_SIGMA_TOLERANCE * float(estimate["mc_se"]) + MC_ABSOLUTE_FLOOR
            checks[f"{label}_lag_{key}"] = {
                "passed": bool(gap <= tolerance),
                "estimate": float(estimate["estimate"]),
                "target": float(target),
                "gap": float(gap),
                "tolerance": float(tolerance),
            }
    return {
        "statistical_tolerance": f"{MC_SIGMA_TOLERANCE:g} batch standard errors + {MC_ABSOLUTE_FLOOR:g} absolute floor",
        "mean_one_checks": mean_checks,
        "all_mean_one_checks_pass": all(item["passed"] for item in mean_checks.values()),
        "exact_level_moment_checks": checks,
        "all_exact_level_moment_checks_pass": all(item["passed"] for item in checks.values()),
    }


def _write_plot_packet(output: Path, continuous: dict[str, Any], markov: dict[str, Any], exact: dict[str, float], mc: dict[str, Any]) -> dict[str, Any]:
    try:
        return _write_plot_packet_impl(output, continuous, markov, exact, mc)
    except Exception as exc:  # pragma: no cover - environment-specific fallback
        return {"status": "failed", "error": str(exc)}


def _write_plot_packet_impl(output: Path, continuous: dict[str, Any], markov: dict[str, Any], exact: dict[str, float], mc: dict[str, Any]) -> dict[str, Any]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - environment-specific fallback
        return {"status": "unavailable", "error": str(exc)}
    labels = [str(x) for x in REPORT_LAGS]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(x, [continuous["moments"][k]["level_covariance"] for k in labels], "o-", label="continuous endpoint")
    ax.plot(x, [markov["moments"][k]["level_covariance"] for k in labels], "s--", label="15-state Markov")
    estimates = [mc["block_level_covariance"][k]["estimate"] for k in labels]
    errors = [MC_SIGMA_TOLERANCE * mc["block_level_covariance"][k]["mc_se"] for k in labels]
    ax.errorbar(x, estimates, yerr=errors, fmt="d", capsize=3, label="annual-simulated block average")
    ax.plot(x, [exact[k] for k in labels], "k:", label="exact block-average covariance")
    ax.set_xticks(x, labels)
    ax.set_xlabel("four-year block lag")
    ax.set_ylabel("level covariance")
    ax.legend(frameon=False)
    fig.tight_layout()
    moments_path = output / "level_covariance_comparison.png"
    fig.savefig(moments_path, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    q = np.asarray(QUANTILE_LEVELS)
    est = np.asarray([mc["block_quantiles"][str(v)]["estimate"] for v in q])
    err = MC_SIGMA_TOLERANCE * np.asarray([mc["block_quantiles"][str(v)]["mc_se"] for v in q])
    ax.errorbar(q, est, yerr=err, fmt="o", capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel("quantile")
    ax.set_ylabel("four-year average level")
    ax.set_title("Simulated block-average distribution")
    fig.tight_layout()
    quantile_path = output / "block_average_quantiles.png"
    fig.savefig(quantile_path, dpi=150)
    plt.close(fig)
    return {"status": "written", "files": [moments_path.name, quantile_path.name]}


def run_diagnostic(mode: str, output: Path, candidate_path: Path, seed: int | None = None) -> dict[str, Any]:
    mode = str(mode).lower()
    if mode not in DEFAULT_BATCHES:
        raise ValueError("mode must be smoke or full")
    design = {
        "mode": mode,
        "batches": DEFAULT_BATCHES[mode],
        "households_per_batch": DEFAULT_HOUSEHOLDS[mode],
        "annual_years_per_batch": DEFAULT_YEARS[mode],
        "block_years": 4,
        "seed": int(DEFAULT_SEED if seed is None else seed),
        "total_time_budget_seconds": TOTAL_TIME_BUDGET_SECONDS[mode],
        "retain_raw_panels": False,
    }
    if design["annual_years_per_batch"] % design["block_years"]:
        raise ValueError("design years must be divisible by block_years")
    prepare_output_directory(output)
    candidate_path = candidate_path.resolve()
    paths = dependency_paths(candidate_path)
    hashes_before = dependency_hashes(paths)
    plan = {
        "status": "planned",
        "diagnostic": "annual_lognormal_to_four_year_block_average",
        "design": design,
        "source_files": {name: str(path) for name, path in paths.items()},
        "source_hashes_before": hashes_before,
        "stop_rule": "stop before the next batch when the mode-specific total time budget is reached; never retry or change the seed",
    }
    write_json(output / "plan.json", plan)
    write_json(output / "progress.json", {
        "status": "planned",
        "completed_batches": 0,
        "total_batches": design["batches"],
        "total_time_budget_seconds": design["total_time_budget_seconds"],
    })
    started = time.perf_counter()
    deadline = started + float(design["total_time_budget_seconds"])
    primitives = load_candidate(candidate_path)
    if primitives["period_years"] != 4:
        raise ValueError("MC loop requires period_years == 4")
    payload = payload_from_candidate(candidate_path, primitives)
    overrides = payload["overrides"]
    plan.update({
        "status": "running",
        "resolved_period_years": primitives["period_years"],
        "payload_fingerprint": payload["payload_fingerprint"],
        "expected_payload_fingerprints": payload["expected_fingerprints"],
    })
    write_json(output / "plan.json", plan)
    continuous = continuous_endpoint_moments(
        primitives["persistent_variance_annual"],
        primitives["transitory_log_variance_period"],
        primitives["rho_period"],
    )
    markov = markov_payload_moments(overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"])
    exact_annual = {str(k): annual_level_covariance(
        primitives["persistent_variance_annual"],
        primitives["transitory_variance_annual"],
        primitives["rho_annual"],
        k,
    ) for k in REPORT_LAGS}
    exact_block = {str(k): block_average_level_covariance(
        primitives["persistent_variance_annual"],
        primitives["transitory_variance_annual"],
        primitives["rho_annual"],
        primitives["period_years"],
        k,
    ) for k in REPORT_LAGS}
    rows: list[dict[str, Any]] = []
    heartbeat = output / "heartbeat.json"
    seed_sequence = np.random.SeedSequence(design["seed"])
    stop_reason: str | None = None
    for batch_index, child_seed in enumerate(seed_sequence.spawn(design["batches"]), start=1):
        if time.perf_counter() >= deadline:
            stop_reason = "time_budget"
            break
        batch_started = time.perf_counter()
        row = _batch_simulation(
            np.random.default_rng(child_seed),
            primitives["persistent_variance_annual"],
            primitives["transitory_variance_annual"],
            primitives["rho_annual"],
            design["annual_years_per_batch"],
            design["households_per_batch"],
        )
        rows.append(row)
        latest = {
            "batch_index": batch_index,
            "elapsed_seconds": time.perf_counter() - started,
            "batch_seconds": time.perf_counter() - batch_started,
            "summary": row,
        }
        write_json(output / "latest_completed_batch.json", latest)
        heartbeat_payload = {
            "status": "running" if batch_index < design["batches"] else "completed",
            "mode": mode,
            "completed_batches": batch_index,
            "total_batches": design["batches"],
            "last_batch_seconds": time.perf_counter() - batch_started,
            "elapsed_seconds": time.perf_counter() - started,
        }
        write_json(heartbeat, heartbeat_payload)
        write_json(output / "progress.json", {
            **heartbeat_payload,
            "time_budget_seconds": design["total_time_budget_seconds"],
            "stop_reason": None,
        })
        print(json.dumps({"event": "batch_complete", **heartbeat_payload}), flush=True)
        if time.perf_counter() >= deadline:
            stop_reason = "time_budget"
            if batch_index < design["batches"]:
                break
    if rows:
        monte_carlo = _aggregate_batch_statistics(rows)
        validity = _validity_checks(monte_carlo, exact_block, exact_annual)
    else:
        monte_carlo = _empty_batch_statistics()
        validity = {
            "statistical_tolerance": f"{MC_SIGMA_TOLERANCE:g} batch standard errors + {MC_ABSOLUTE_FLOOR:g} absolute floor",
            "mean_one_checks": {},
            "all_mean_one_checks_pass": False,
            "exact_level_moment_checks": {},
            "all_exact_level_moment_checks_pass": False,
        }
    hashes_after = dependency_hashes(paths)
    dependencies_unchanged = hashes_before == hashes_after
    complete = len(rows) == design["batches"] and stop_reason is None
    validity["all_dependencies_unchanged"] = dependencies_unchanged
    validity["all_batches_completed"] = complete
    validity["completed_batches"] = len(rows)
    validity["stop_reason"] = stop_reason
    result = {
        "status": "running",
        "diagnostic": "annual_lognormal_to_four_year_block_average",
        "design": design,
        "source_files": {name: str(path) for name, path in paths.items()},
        "source_hashes_before": hashes_before,
        "source_hashes_after": hashes_after,
        "payload_contract": {
            "rebuilt_15_state_payload_fingerprint": payload["payload_fingerprint"],
            "candidate_json_fingerprint": payload["candidate_fingerprint"],
            "expected_candidate_payload_fingerprints": payload["expected_fingerprints"],
            "expected_candidate_json_sha256": payload["expected_candidate_hashes"],
        },
        "annual_primitives": {
            "rho_annual": primitives["rho_annual"],
            "persistent_log_variance": primitives["persistent_variance_annual"],
            "transitory_log_variance": primitives["transitory_variance_annual"],
        },
        "period_mapping": {
            "rho_period": primitives["rho_period"],
            "transitory_log_variance_period": primitives["transitory_log_variance_period"],
            "mapping_is_transitory_moment_match_only": True,
            "persistent_path_is_not_block_averaged": True,
        },
        "continuous_endpoint_moments": continuous,
        "markov_payload_moments": markov,
        "exact_annual_level_mean": 1.0,
        "exact_annual_level_covariances": exact_annual,
        "exact_block_average_level_mean": 1.0,
        "exact_block_average_level_covariances": exact_block,
        "monte_carlo": monte_carlo,
        "validity": validity,
        "elapsed_seconds": time.perf_counter() - started,
    }
    plots = (
        _write_plot_packet(output, continuous, markov, exact_block, monte_carlo)
        if rows
        else {"status": "failed", "files": [], "error": "no Monte Carlo batch completed"}
    )
    result["plots"] = plots
    required_plots = ("level_covariance_comparison.png", "block_average_quantiles.png")
    plots_present = plots.get("status") == "written" and all((output / name).is_file() for name in required_plots)
    validity["all_required_plots_present"] = plots_present
    result["validity"] = validity
    result["status"] = final_status(
        requested_batches=design["batches"],
        completed_batches=len(rows),
        exact_checks_pass=bool(validity["all_exact_level_moment_checks_pass"]),
        mean_checks_pass=bool(validity["all_mean_one_checks_pass"]),
        plots_present=plots_present,
        dependencies_unchanged=dependencies_unchanged,
        stop_reason=stop_reason,
    )
    write_json(output / "receipt.json", result)
    write_json(output / "progress.json", {
        "status": result["status"],
        "completed_batches": len(rows),
        "total_batches": design["batches"],
        "elapsed_seconds": result["elapsed_seconds"],
        "time_budget_seconds": design["total_time_budget_seconds"],
        "stop_reason": stop_reason,
    })
    summary_lines = [
        "# Income time aggregation diagnostic",
        "",
        f"Mode: `{mode}`; batches: `{design['batches']}`; households per batch: `{design['households_per_batch']}`; annual years: `{design['annual_years_per_batch']}`.",
        "",
        "The annual process is simulated with stationary persistent log variance and iid transitory log variance from the retained no-fixed candidate. Four-year levels are arithmetic averages of four annual mean-one levels. The 15-state and continuous endpoint objects are reported as approximations for comparison; neither is labeled an exact block-average process.",
        "",
        "For annual mean-one levels, the exact covariance is `exp(Vp * rho^|d| + Ve * 1[d=0]) - 1`; block covariance is the average of its 16 annual pair covariances. Block log moments are reported from simulation because the log of an arithmetic average has no matching closed form here.",
        "",
        f"Receipt status: **{result['status']}**. Mean-one checks: **{'PASS' if validity['all_mean_one_checks_pass'] else 'FAIL'}**; exact level-moment checks: **{'PASS' if validity['all_exact_level_moment_checks_pass'] else 'FAIL'}** under {validity['statistical_tolerance']}.",
        "",
        "See `receipt.json` for complete moments, Monte Carlo SEs, source hashes, and validity flags. No raw panels are retained.",
        "",
        f"Supplemental plots: `{', '.join(plots.get('files', []))}`.",
    ]
    (output / "README.md").write_text("\n".join(summary_lines) + "\n")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(DEFAULT_BATCHES), default="smoke")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--candidate-json", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args(argv)
    result = run_diagnostic(args.mode, args.output, args.candidate_json, args.seed)
    return 0 if result["status"] == "completed" else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
