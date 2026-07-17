#!/usr/bin/env python3
"""Run the two fixed-theta acceptance cases for the entrant-feasibility repair.

This is a diagnostic-only driver.  It reproduces the July 10 combined-spec
Nb=120 evaluation environment, runs the current-bound theta (expected to be
inadmissible), then runs the same theta with ``c_bar_0=1.26`` (expected to
solve strictly).  It does not calibrate, run policy counterfactuals, or edit
production configuration.

Run from ``code/model`` with the project virtual environment.  The driver
forces NUMBA/OMP/MKL/OpenBLAS to one thread before importing NumPy so the
acceptance run is deterministic and matches the production worker contract.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any


THREAD_ENVIRONMENT = {
    "NUMBA_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
}
THREAD_ENVIRONMENT_BEFORE = {name: os.environ.get(name) for name in THREAD_ENVIRONMENT}
for _name, _value in THREAD_ENVIRONMENT.items():
    os.environ[_name] = _value

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

import numpy as np  # noqa: E402  (thread environment must be fixed first)

from intergen_housing_fertility.calibration import (  # noqa: E402
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import income_process_overrides  # noqa: E402
from intergen_housing_fertility.production_profile import (  # noqa: E402
    PRODUCTION_INCOME_STATES,
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    production_profile_metadata,
    production_profile_overrides,
    validate_production_profile,
)
from intergen_housing_fertility.solver import (  # noqa: E402
    DEAD_MASS_TOL,
    DEAD_VALUE_CUTOFF,
    InfeasibleThetaError,
    entry_wealth_grid_weights,
    income_at_state,
    income_transition_values,
    run_model_cp_dt,
)


MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0
PERIOD_YEARS = 4.0
DIAGNOSTIC_C_BAR_0 = 1.26
BUDGET_TOL = 1e-9
POSITIVE_MASS_TOL = 1e-12
EXPECTED_TARGET_COUNT = 15

DEFAULT_RECORD = (
    ROOT
    / "output/model/combined_recalibration/overnight_20260710_report/current_bound_best.json"
)
DEFAULT_SPEC = ROOT / "output/model/full_audit_20260711/FEASIBILITY_FIX_SPEC_FOR_CODEX.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--record",
        type=Path,
        default=DEFAULT_RECORD,
        help="current-bound candidate JSON (default: July 10 extracted record)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "output directory (default: output/model/"
            "feasibility_fix_diagnostic_<today YYYYMMDD>)"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace this driver's existing acceptance artifacts in --outdir",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="show equilibrium iteration output (quiet by default)",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_reference(path: Path) -> tuple[dict[str, float], dict[str, Any]]:
    payload = json.loads(path.read_text())
    raw_theta = payload.get("theta")
    if not isinstance(raw_theta, dict):
        raise ValueError(f"{path} does not contain a theta object")
    theta = {str(name): float(value) for name, value in raw_theta.items()}
    if not math.isclose(theta.get("c_bar_0", math.nan), 1.28, rel_tol=0.0, abs_tol=0.0):
        raise ValueError(
            f"Run A requires the exact current-bound c_bar_0=1.28; got {theta.get('c_bar_0')}"
        )
    return theta, payload


def build_environment() -> tuple[dict[str, float], dict[str, float], dict[str, Any], dict[str, Any]]:
    if PRODUCTION_SEARCH_NB != 120:
        raise ValueError(
            f"acceptance contract requires production Nb=120, found {PRODUCTION_SEARCH_NB}"
        )
    validate_production_profile(
        PRODUCTION_PROFILE_NAME,
        J=PRODUCTION_J,
        Nb=PRODUCTION_SEARCH_NB,
        n_house=5,
        income_states=PRODUCTION_INCOME_STATES,
        target_set=PRODUCTION_TARGET_SET,
        max_iter_eq=PRODUCTION_MAX_ITER_EQ,
        stage="search",
    )
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    targets = dict(targets)
    weights = dict(weights)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = ROOMS_WEIGHT
    if len(targets) != EXPECTED_TARGET_COUNT:
        raise ValueError(
            f"acceptance contract requires {EXPECTED_TARGET_COUNT} targets, found {len(targets)}"
        )
    income = income_process_overrides(
        PRODUCTION_INCOME_STATES,
        "rouwenhorst",
        MATCHED_ANNUAL_INNOVATION_SD,
        MATCHED_ANNUAL_RHO,
    )
    fixed = production_profile_overrides()
    fixed.update(
        {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": np.array([1.75]),
            "normalize_bequest_utility": True,
            # Pin the finalized feasibility design explicitly rather than
            # relying silently on future parameter defaults.
            "lambda_d": 0.0,
            "debt_taper_start_age": 42.0,
            "debt_taper_end_age": 62.0,
            "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
        }
    )
    return targets, weights, income, fixed


def validate_reference_fit(
    payload: dict[str, Any], targets: dict[str, float], weights: dict[str, float]
) -> dict[str, dict[str, Any]]:
    rows = payload.get("fit")
    if not isinstance(rows, list):
        raise ValueError("reference record has no fit table")
    by_name = {str(row["moment"]): dict(row) for row in rows}
    if set(by_name) != set(targets):
        raise ValueError(
            "reference and active target systems differ: "
            f"reference_only={sorted(set(by_name) - set(targets))}, "
            f"active_only={sorted(set(targets) - set(by_name))}"
        )
    for name in targets:
        row = by_name[name]
        if not math.isclose(float(row["target"]), float(targets[name]), rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"target drift for {name}: {row['target']} != {targets[name]}")
        if not math.isclose(float(row["weight"]), float(weights[name]), rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"weight drift for {name}: {row['weight']} != {weights[name]}")
    return by_name


def active_budget_slack(sol: Any, P: Any) -> np.ndarray:
    """Budget slack aligned with the realized post-choice cross-section."""

    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float)
    z_grid, _, _ = income_transition_values(P)
    p = np.asarray(sol.p_eq, dtype=float).reshape(-1)
    r_hat = float(P.user_cost_rate) * p
    slack = np.empty_like(g)
    for i in range(int(P.I)):
        for j in range(int(P.J)):
            for zz, z_value in enumerate(z_grid):
                resources = float(P.R_gross) * b_grid + income_at_state(
                    P, i, j, float(z_value)
                )
                for ten in range(1 + int(P.n_house)):
                    c = np.asarray(sol.c_pol[:, ten, i, j, zz, :, :], dtype=float)
                    bp = np.asarray(sol.bp_pol[:, ten, i, j, zz, :, :], dtype=float)
                    if ten == 0:
                        h = np.asarray(sol.hR_pol[:, ten, i, j, zz, :, :], dtype=float)
                        uses = c + float(r_hat[i]) * h + bp
                    else:
                        house = float(P.H_own[ten - 1])
                        size_cost = (
                            float(getattr(P, "owner_size_cost", 0.0))
                            * float(p[i])
                            * max(
                                house - float(getattr(P, "owner_size_cost_ref", 6.0)),
                                0.0,
                            )
                            ** float(getattr(P, "owner_size_cost_power", 2.0))
                        )
                        operating_cost = (
                            (float(P.delta) + float(P.tau_H)) * float(p[i]) * house
                            + size_cost
                        )
                        uses = c + operating_cost + bp
                    slack[:, ten, i, j, zz, :, :] = resources[:, None, None] - uses
    return slack


def entry_diagnostics(sol: Any, P: Any) -> dict[str, Any]:
    """Reconstruct the exact entrant mixture and measure dead/artifact mass."""

    b_grid = np.asarray(sol.b_grid, dtype=float)
    z_grid, z_weights, _ = income_transition_values(P)
    total_mass = 0.0
    dead_mass = 0.0
    dead_uniform_mass = 0.0
    rows: list[dict[str, Any]] = []
    for i in range(int(P.I)):
        loc_mass = float(np.asarray(P.entry_by_loc, dtype=float).reshape(-1)[i])
        for zz, (z_value, z_weight) in enumerate(zip(z_grid, z_weights)):
            indices, node_weights = entry_wealth_grid_weights(
                b_grid, P, i=i, j=0, z_value=float(z_value)
            )
            for raw_index, raw_weight in zip(indices, node_weights):
                b_index = int(raw_index)
                mass = loc_mass * float(z_weight) * float(raw_weight)
                value = float(sol.V[b_index, 0, i, 0, zz, 0, 0])
                probabilities = np.asarray(
                    sol.fert_probs[b_index, 0, i, 0, zz, :], dtype=float
                )
                is_dead = value <= DEAD_VALUE_CUTOFF
                is_uniform = bool(
                    probabilities.size > 0
                    and np.max(np.abs(probabilities - 1.0 / probabilities.size)) <= 1e-8
                )
                total_mass += mass
                dead_mass += mass if is_dead else 0.0
                dead_uniform_mass += mass if is_dead and is_uniform else 0.0
                if is_dead:
                    rows.append(
                        {
                            "location": i,
                            "income_state": zz,
                            "z": float(z_value),
                            "asset_index": b_index,
                            "b": float(b_grid[b_index]),
                            "mass": mass,
                            "value": value,
                            "fertility_probabilities": probabilities.tolist(),
                        }
                    )
    denominator = max(total_mass, 1e-300)
    return {
        "total_entry_flow_before_population_normalization": total_mass,
        "dead_entry_mass": dead_mass,
        "dead_entry_mass_share": dead_mass / denominator,
        "dead_uniform_fertility_mass": dead_uniform_mass,
        "dead_uniform_fertility_mass_share": dead_uniform_mass / denominator,
        "dead_entry_rows": rows[:20],
    }


def solution_diagnostics(sol: Any, P: Any) -> dict[str, Any]:
    g = np.asarray(sol.g, dtype=float)
    values = np.asarray(sol.V, dtype=float)
    fertility = np.asarray(sol.fert_probs, dtype=float)
    total_mass = float(np.sum(g))
    dead = values <= DEAD_VALUE_CUTOFF
    positive_dead_mass = float(np.sum(np.where(dead, g, 0.0)))

    # Fertility is chosen only from the childless state.  A dead childless row
    # must now have the exact all-zero mask, never the old uniform softmax.
    childless_dead = values[..., 0, 0] <= DEAD_VALUE_CUTOFF
    n_options = fertility.shape[-1]
    uniform = np.max(np.abs(fertility - 1.0 / n_options), axis=-1) <= 1e-8
    probability_sums = np.sum(fertility, axis=-1)
    all_zero = np.max(np.abs(fertility), axis=-1) <= 1e-15

    slack = active_budget_slack(sol, P)
    positive_mass = g > POSITIVE_MASS_TOL
    bad_budget = positive_mass & (slack < -BUDGET_TOL)
    positive_mass_total = float(np.sum(g[positive_mass]))
    bad_budget_mass = float(np.sum(g[bad_budget]))
    if np.any(positive_mass):
        minimum_on_path_slack = float(np.min(slack[positive_mass]))
    else:
        minimum_on_path_slack = math.nan

    worst_budget_rows: list[dict[str, Any]] = []
    candidate_indices = np.argwhere(positive_mass)
    if candidate_indices.size:
        ranked = sorted(
            candidate_indices,
            key=lambda idx: float(slack[tuple(idx)]),
        )[:8]
        z_grid, _, _ = income_transition_values(P)
        b_grid = np.asarray(sol.b_grid, dtype=float)
        for idx in ranked:
            b, ten, loc, j, zz, parity, child_state = (int(value) for value in idx)
            key = (b, ten, loc, j, zz, parity, child_state)
            worst_budget_rows.append(
                {
                    "age": float(P.age_start + P.da * j),
                    "b": float(b_grid[b]),
                    "tenure": ten,
                    "location": loc,
                    "z": float(z_grid[zz]),
                    "parity": parity,
                    "child_state": child_state,
                    "mass": float(g[key]),
                    "budget_slack": float(slack[key]),
                }
            )

    return {
        "feasibility_gate_completed_without_exception": True,
        "dead_value_cutoff": DEAD_VALUE_CUTOFF,
        "dead_mass_tolerance": DEAD_MASS_TOL,
        "population_mass": total_mass,
        "positive_mass_on_dead_nodes": positive_dead_mass,
        "positive_mass_on_dead_nodes_share": positive_dead_mass / max(total_mass, 1e-300),
        "fertility": {
            "uniform_probability_rows_all": int(np.sum(uniform)),
            "dead_uniform_probability_rows": int(np.sum(childless_dead & uniform)),
            "dead_nonzero_probability_rows": int(
                np.sum(childless_dead & (probability_sums > 1e-15))
            ),
            "dead_all_zero_probability_rows": int(np.sum(childless_dead & all_zero)),
            "maximum_probability_sum_on_dead_rows": (
                float(np.max(probability_sums[childless_dead]))
                if np.any(childless_dead)
                else 0.0
            ),
        },
        "budget": {
            "positive_mass_threshold": POSITIVE_MASS_TOL,
            "budget_slack_tolerance": BUDGET_TOL,
            "positive_mass_total": positive_mass_total,
            "budget_violation_mass": bad_budget_mass,
            "budget_violation_mass_share": bad_budget_mass
            / max(positive_mass_total, 1e-300),
            "minimum_budget_slack_on_positive_mass": minimum_on_path_slack,
            "worst_positive_mass_rows": worst_budget_rows,
        },
        "entrant": entry_diagnostics(sol, P),
    }


def fit_rows(
    moments: dict[str, float], targets: dict[str, float], weights: dict[str, float]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name in targets:
        model = float(moments.get(name, math.nan))
        target = float(targets[name])
        weight = float(weights[name])
        gap = model - target
        rows.append(
            {
                "moment": name,
                "target": target,
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )
    return rows


def run_case(
    label: str,
    theta: dict[str, float],
    *,
    targets: dict[str, float],
    weights: dict[str, float],
    income: dict[str, Any],
    fixed: dict[str, Any],
    verbose: bool,
) -> dict[str, Any]:
    overrides = {
        **base_overrides(
            J=PRODUCTION_J,
            Nb=PRODUCTION_SEARCH_NB,
            n_house=5,
            max_iter_eq=PRODUCTION_MAX_ITER_EQ,
        ),
        **fixed,
        **income,
        **theta,
    }
    started = time.perf_counter()
    try:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=verbose)
        moments = extract_moments(sol, P)
        market_residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
        timings = dict(getattr(sol, "timings", {}))
        solver_strict = bool(
            timings.get("strict_converged", getattr(sol, "converged", False))
        )
        strict = bool(
            solver_strict
            and np.isfinite(market_residual)
            and market_residual <= float(getattr(P, "tol_eq", 1e-4))
        )
        return {
            "label": label,
            "status": "ok",
            "theta": theta,
            "rank_loss": diagnostic_loss(moments, targets=targets, weights=weights),
            "strict_converged": strict,
            "market_residual": market_residual,
            "tol_eq": float(getattr(P, "tol_eq", 1e-4)),
            "p_eq": np.asarray(p_eq, dtype=float),
            "moments": moments,
            "target_fit": fit_rows(moments, targets, weights),
            "diagnostics": solution_diagnostics(sol, P),
            "timings": timings,
            "feasibility_error": "",
            "feasibility_census": [],
            "elapsed_sec": time.perf_counter() - started,
        }
    except InfeasibleThetaError as exc:
        return {
            "label": label,
            "status": "infeasible_theta",
            "theta": theta,
            "rank_loss": math.inf,
            "strict_converged": False,
            "market_residual": math.inf,
            "tol_eq": math.nan,
            "p_eq": [math.nan],
            "moments": {},
            "target_fit": [],
            "diagnostics": {},
            "timings": {},
            "feasibility_error": str(exc),
            "feasibility_stage": exc.stage,
            "feasibility_dead_mass": exc.dead_mass,
            "feasibility_census": list(exc.census),
            "elapsed_sec": time.perf_counter() - started,
        }
    except Exception as exc:  # noqa: BLE001 - both records must be written on failure.
        return {
            "label": label,
            "status": f"failed: {type(exc).__name__}: {exc}",
            "theta": theta,
            "rank_loss": math.inf,
            "strict_converged": False,
            "market_residual": math.inf,
            "tol_eq": math.nan,
            "p_eq": [math.nan],
            "moments": {},
            "target_fit": [],
            "diagnostics": {},
            "timings": {},
            "feasibility_error": "",
            "feasibility_census": [],
            "elapsed_sec": time.perf_counter() - started,
        }


def comparison_rows(
    reference_fit: dict[str, dict[str, Any]],
    run_b: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
) -> list[dict[str, Any]]:
    repaired_moments = dict(run_b.get("moments") or {})
    rows: list[dict[str, Any]] = []
    for name in targets:
        prior = reference_fit[name]
        repaired_raw = repaired_moments.get(name)
        repaired = float(repaired_raw) if repaired_raw is not None else math.nan
        target = float(targets[name])
        weight = float(weights[name])
        repaired_gap = repaired - target
        contaminated_model = float(prior["model"])
        rows.append(
            {
                "moment": name,
                "target": target,
                "weight": weight,
                "contaminated_14p780_model": contaminated_model,
                "contaminated_gap": float(prior["gap"]),
                "contaminated_loss_contribution": float(prior["loss_contribution"]),
                "run_b_cbar1p26_model": repaired,
                "run_b_gap": repaired_gap,
                "run_b_loss_contribution": weight * repaired_gap * repaired_gap,
                "run_b_minus_contaminated": repaired - contaminated_model,
            }
        )
    return rows


def acceptance_checks(
    run_a: dict[str, Any], run_b: dict[str, Any], targets: dict[str, float]
) -> list[dict[str, Any]]:
    diagnostics = dict(run_b.get("diagnostics") or {})
    entrant = dict(diagnostics.get("entrant") or {})
    fertility = dict(diagnostics.get("fertility") or {})
    budget = dict(diagnostics.get("budget") or {})
    moments = dict(run_b.get("moments") or {})
    run_b_entry_dead_mass = (
        entrant.get("dead_entry_mass")
        if entrant
        else run_b.get("feasibility_dead_mass", math.inf)
    )
    checks = [
        {
            "gate": "run_a_is_infeasible_theta",
            "pass": run_a.get("status") == "infeasible_theta",
            "observed": run_a.get("status"),
        },
        {
            "gate": "run_a_has_entry_census",
            "pass": bool(run_a.get("feasibility_census"))
            and run_a.get("feasibility_stage") == "entry",
            "observed": {
                "stage": run_a.get("feasibility_stage"),
                "rows": len(run_a.get("feasibility_census") or []),
            },
        },
        {
            "gate": "run_b_strict_market_convergence",
            "pass": run_b.get("status") == "ok" and bool(run_b.get("strict_converged")),
            "observed": {
                "status": run_b.get("status"),
                "strict": run_b.get("strict_converged"),
                "residual": run_b.get("market_residual"),
                "tolerance": run_b.get("tol_eq"),
            },
        },
        {
            "gate": "run_b_zero_dead_entry_mass",
            "pass": float(run_b_entry_dead_mass) <= DEAD_MASS_TOL,
            "observed": run_b_entry_dead_mass,
        },
        {
            "gate": "run_b_zero_positive_population_dead_mass",
            "pass": float(diagnostics.get("positive_mass_on_dead_nodes", math.inf))
            <= DEAD_MASS_TOL,
            "observed": diagnostics.get("positive_mass_on_dead_nodes"),
        },
        {
            "gate": "run_b_zero_uniform_fertility_artifact_rows",
            "pass": int(fertility.get("dead_uniform_probability_rows", -1)) == 0
            and int(fertility.get("dead_nonzero_probability_rows", -1)) == 0,
            "observed": {
                "dead_uniform": fertility.get("dead_uniform_probability_rows"),
                "dead_nonzero": fertility.get("dead_nonzero_probability_rows"),
            },
        },
        {
            "gate": "run_b_budget_slack_on_positive_mass",
            "pass": float(budget.get("budget_violation_mass", math.inf)) <= DEAD_MASS_TOL
            and float(budget.get("minimum_budget_slack_on_positive_mass", -math.inf))
            >= -BUDGET_TOL,
            "observed": {
                "violation_mass": budget.get("budget_violation_mass"),
                "minimum_slack": budget.get("minimum_budget_slack_on_positive_mass"),
                "tolerance": BUDGET_TOL,
            },
        },
        {
            "gate": "run_b_all_15_moments_finite",
            "pass": len(targets) == EXPECTED_TARGET_COUNT
            and all(np.isfinite(float(moments.get(name, math.nan))) for name in targets),
            "observed": {
                "target_count": len(targets),
                "finite_target_moments": sum(
                    np.isfinite(float(moments.get(name, math.nan))) for name in targets
                ),
            },
        },
    ]
    return checks


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: Any, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "—"
    return f"{number:.{digits}f}" if np.isfinite(number) else "—"


def summary_markdown(
    run_a: dict[str, Any],
    run_b: dict[str, Any],
    checks: list[dict[str, Any]],
    comparison: list[dict[str, Any]],
    reference: dict[str, Any],
) -> str:
    passed = all(bool(row["pass"]) for row in checks)
    census = run_a.get("feasibility_census") or []
    lines = [
        "# Entrant-feasibility fix: Nb=120 acceptance",
        "",
        f"**Overall verdict: {'PASS' if passed else 'FAIL'}.** This is a fixed-theta diagnostic, not a calibration or policy run.",
        "",
        "Configuration: July 9 production profile, combined fixed specification, five-state matched Rouwenhorst income process, "
        "`lambda_d=0`, debt taper ages 42–62, and one NUMBA/OMP thread.",
        "",
        "## Acceptance gates",
        "",
        "| Gate | Result |",
        "|---|---:|",
    ]
    for row in checks:
        lines.append(f"| `{row['gate']}` | {'PASS' if row['pass'] else 'FAIL'} |")

    lines.extend(
        [
            "",
            "## Run A — current theta, c̄₀=1.28",
            "",
            f"Expected `infeasible_theta`; observed `{run_a.get('status')}`. "
            f"Stage: `{run_a.get('feasibility_stage', '—')}`; dead mass: "
            f"{fmt(run_a.get('feasibility_dead_mass'), 12)}.",
            "",
            "| Age | z | b | Mass | Slack |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in census[:8]:
        lines.append(
            f"| {fmt(row.get('age'), 1)} | {fmt(row.get('z'), 6)} | "
            f"{fmt(row.get('b'), 6)} | {fmt(row.get('mass'), 9)} | "
            f"{fmt(row.get('slack'), 9)} |"
        )
    if not census:
        lines.append("| — | — | — | — | — |")

    diagnostics = dict(run_b.get("diagnostics") or {})
    entrant = dict(diagnostics.get("entrant") or {})
    fertility = dict(diagnostics.get("fertility") or {})
    budget = dict(diagnostics.get("budget") or {})
    run_b_census = run_b.get("feasibility_census") or []
    run_b_configured_cbar = float(
        (run_b.get("theta") or {}).get("c_bar_0", DIAGNOSTIC_C_BAR_0)
    )
    run_b_row_cbar_caps = [
        {
            "z": row.get("z"),
            "b": row.get("b"),
            "maximum_c_bar_0": run_b_configured_cbar + float(row["slack"]),
        }
        for row in run_b_census
        if row.get("slack") is not None and np.isfinite(float(row["slack"]))
    ]
    lines.extend(
        [
            "",
            "## Run B — diagnostic-only c̄₀=1.26",
            "",
            f"Status `{run_b.get('status')}`; strict convergence `{run_b.get('strict_converged')}`; "
            f"residual {fmt(run_b.get('market_residual'), 9)}; loss {fmt(run_b.get('rank_loss'), 6)}.",
            "",
            f"Dead entrant mass: {fmt(entrant.get('dead_entry_mass'), 12)}. "
            f"Positive population mass on dead nodes: {fmt(diagnostics.get('positive_mass_on_dead_nodes'), 12)}. "
            f"Dead uniform fertility rows: {fertility.get('dead_uniform_probability_rows', '—')}. "
            f"Budget-violation mass: {fmt(budget.get('budget_violation_mass'), 12)}; "
            f"minimum on-path slack: {fmt(budget.get('minimum_budget_slack_on_positive_mass'), 12)}.",
            (
                f"Run B rejection stage: `{run_b.get('feasibility_stage', '—')}`; "
                f"dead mass: {fmt(run_b.get('feasibility_dead_mass'), 12)}."
                if run_b.get("status") == "infeasible_theta"
                else ""
            ),
        ]
    )
    if run_b_census:
        lines.extend(
            [
                "",
                "Run B infeasible-entry census:",
                "",
                "| Age | z | b | Mass | Slack |",
                "|---:|---:|---:|---:|---:|",
            ]
        )
        for row in run_b_census[:8]:
            lines.append(
                f"| {fmt(row.get('age'), 1)} | {fmt(row.get('z'), 6)} | "
                f"{fmt(row.get('b'), 6)} | {fmt(row.get('mass'), 9)} | "
                f"{fmt(row.get('slack'), 9)} |"
            )
    if run_b.get("status") == "infeasible_theta" and run_b_row_cbar_caps:
        most_restrictive_cbar = min(
            float(row["maximum_c_bar_0"]) for row in run_b_row_cbar_caps
        )
        row_caps_text = "; ".join(
            f"z={fmt(row['z'], 6)}, b={fmt(row['b'], 6)}: "
            f"c̄₀≤{fmt(row['maximum_c_bar_0'], 9)}"
            for row in run_b_row_cbar_caps
        )
        lines.extend(
            [
                "",
                "**Acceptance-spec inconsistency.** With `lambda_d=0` and "
                "`s_next=1`, the rollover rule gives `b'=b`, so entrant slack is "
                "`y - c̄₀ - r*h̄₀ + q*b`. The reported low-z, negative-b rows "
                "therefore cannot pass at c̄₀=1.26; this is not a code-gate failure. "
                "Each row implies `maximum c̄₀ = configured c̄₀ + slack`: "
                f"{row_caps_text}. The most restrictive value is "
                f"c̄₀≤{fmt(most_restrictive_cbar, 9)}.",
            ]
        )
    lines.extend(
        [
            "",
            "## Headline moment comparison",
            "",
            "| Moment | Contaminated 14.780 record | Run B | Change |",
            "|---|---:|---:|---:|",
        ]
    )
    headline_names = {
        "tfr",
        "childless_rate",
        "own_rate",
        "old_age_own_rate",
        "aggregate_mean_occupied_rooms_18_85",
    }
    for row in comparison:
        if row["moment"] in headline_names:
            lines.append(
                f"| `{row['moment']}` | {fmt(row['contaminated_14p780_model'])} | "
                f"{fmt(row['run_b_cbar1p26_model'])} | {fmt(row['run_b_minus_contaminated'])} |"
            )
    lines.extend(
        [
            "",
            f"Reference loss: {fmt(reference.get('loss'), 9)}. The complete 15-moment comparison is in "
            "`target_fit_comparison.csv` and `target_fit_comparison.json`.",
            "",
        ]
    )
    return "\n".join(lines)


def check_output_paths(outdir: Path, overwrite: bool) -> None:
    managed = [
        outdir / "run_a.json",
        outdir / "run_b.json",
        outdir / "target_fit_comparison.json",
        outdir / "target_fit_comparison.csv",
        outdir / "SUMMARY.md",
    ]
    existing = [path for path in managed if path.exists()]
    if existing and not overwrite:
        joined = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"refusing to overwrite existing artifacts: {joined}")
    outdir.mkdir(parents=True, exist_ok=True)


def main() -> None:
    args = parse_args()
    outdir = args.outdir or (
        ROOT / "output/model" / f"feasibility_fix_diagnostic_{date.today():%Y%m%d}"
    )
    check_output_paths(outdir, args.overwrite)
    theta, reference = load_reference(args.record)
    targets, weights, income, fixed = build_environment()
    reference_fit = validate_reference_fit(reference, targets, weights)

    run_a = run_case(
        "run_a_current_theta_expected_infeasible",
        dict(theta),
        targets=targets,
        weights=weights,
        income=income,
        fixed=fixed,
        verbose=args.verbose,
    )
    print(
        f"Run A: status={run_a['status']} stage={run_a.get('feasibility_stage')} "
        f"census_rows={len(run_a.get('feasibility_census') or [])}",
        flush=True,
    )
    for row in run_a.get("feasibility_census") or []:
        print("  " + json.dumps(jsonable(row), sort_keys=True), flush=True)

    theta_b = dict(theta)
    theta_b["c_bar_0"] = DIAGNOSTIC_C_BAR_0
    run_b = run_case(
        "run_b_cbar1p26_diagnostic_only",
        theta_b,
        targets=targets,
        weights=weights,
        income=income,
        fixed=fixed,
        verbose=args.verbose,
    )
    print(
        f"Run B: status={run_b['status']} strict={run_b['strict_converged']} "
        f"loss={fmt(run_b.get('rank_loss'), 9)} residual={fmt(run_b.get('market_residual'), 9)}",
        flush=True,
    )

    comparison = comparison_rows(reference_fit, run_b, targets, weights)
    checks = acceptance_checks(run_a, run_b, targets)
    passed = all(bool(row["pass"]) for row in checks)
    provenance = {
        "driver": str(Path(__file__).resolve()),
        "record": str(args.record.resolve()),
        "record_sha256": sha256_file(args.record),
        "spec": str(DEFAULT_SPEC.resolve()),
        "spec_sha256": sha256_file(DEFAULT_SPEC) if DEFAULT_SPEC.exists() else None,
        "production_profile": production_profile_metadata(),
        "acceptance_grid": {
            "J": PRODUCTION_J,
            "Nb": PRODUCTION_SEARCH_NB,
            "n_house": 5,
            "income_states": PRODUCTION_INCOME_STATES,
            "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
            "target_set": PRODUCTION_TARGET_SET,
        },
        "combined_fixed_spec": {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": [1.75],
            "normalize_bequest_utility": True,
            "lambda_d": 0.0,
            "debt_taper_start_age": 42.0,
            "debt_taper_end_age": 62.0,
        },
        "income_process": {
            "method": "rouwenhorst",
            "annual_rho": MATCHED_ANNUAL_RHO,
            "annual_innovation_sd": MATCHED_ANNUAL_INNOVATION_SD,
        },
        "thread_environment_before_driver": THREAD_ENVIRONMENT_BEFORE,
        "thread_environment_enforced": THREAD_ENVIRONMENT,
    }
    run_a["acceptance_expectation"] = "infeasible_theta at current theta"
    run_a["provenance"] = provenance
    run_b["acceptance_expectation"] = (
        "diagnostic-only strict solve with zero dead mass, fertility artifact, and budget violations"
    )
    run_b["provenance"] = provenance
    comparison_payload = {
        "status": "pass" if passed else "fail",
        "reference_loss": float(reference.get("loss", math.nan)),
        "reference_market_residual": float(reference.get("market_residual", math.nan)),
        "run_b_loss": run_b.get("rank_loss"),
        "run_b_market_residual": run_b.get("market_residual"),
        "acceptance_checks": checks,
        "rows": comparison,
        "provenance": provenance,
    }

    (outdir / "run_a.json").write_text(
        json.dumps(jsonable(run_a), indent=2, sort_keys=True) + "\n"
    )
    (outdir / "run_b.json").write_text(
        json.dumps(jsonable(run_b), indent=2, sort_keys=True) + "\n"
    )
    (outdir / "target_fit_comparison.json").write_text(
        json.dumps(jsonable(comparison_payload), indent=2, sort_keys=True) + "\n"
    )
    write_csv(outdir / "target_fit_comparison.csv", jsonable(comparison))
    (outdir / "SUMMARY.md").write_text(
        summary_markdown(run_a, run_b, checks, comparison, reference)
    )

    for row in checks:
        print(f"{'PASS' if row['pass'] else 'FAIL'} {row['gate']}", flush=True)
    print(f"Artifacts: {outdir.resolve()}", flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
