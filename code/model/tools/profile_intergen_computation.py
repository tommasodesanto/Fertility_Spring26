#!/usr/bin/env python3
"""Profile the production M5 solver without modifying production code.

The driver reconstructs the exact M5 mechanism and parameter point from the
collected report.  It can profile either the tight general-equilibrium solve or
a single fixed-price solve at the collected equilibrium price.  Results are
written as machine-readable JSON plus standard ``pstats`` artifacts.
"""

from __future__ import annotations

import argparse
import cProfile
import io
import json
import math
import os
import platform
import pstats
import resource
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numba
import numpy as np

from intergen_housing_fertility.calibration import (
    diagnostic_loss,
    extract_moments,
)
from intergen_housing_fertility import solver as solver_module
from intergen_housing_fertility.solver import run_model_cp_dt
from tools.run_intergen_bequest_exit_chain import (
    arm_contract,
    common_overrides,
    load_theta,
    target_system,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULTS = (
    ROOT
    / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/computational_audit_20260719"
TARGET_SET = "candidate_replacement_income_disciplined_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--mode",
        choices=("tight-ge", "loose-ge", "fixed-price", "fixed-price-fast"),
        default="tight-ge",
    )
    parser.add_argument("--top", type=int, default=100)
    parser.add_argument(
        "--numba-scatter",
        action="store_true",
        help="Enable the existing compiled scatter switch for an audit benchmark.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, SimpleNamespace):
        return {str(key): jsonable(item) for key, item in vars(value).items()}
    return value


def load_collected_winner(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    winner = payload.get("winners", {}).get("M5")
    if not isinstance(winner, dict):
        raise ValueError(f"M5 winner is missing from {path}")
    return winner


def build_overrides(args: argparse.Namespace, winner: dict[str, Any]) -> dict[str, Any]:
    tight = args.mode == "tight-ge"
    max_iter_eq = 40 if tight else 10
    tol_eq = 2.5e-5 if tight else 1e-4
    contract_args = SimpleNamespace(
        arm="M5",
        ltv_terminal=0.0,
        theta1=0.25,
        seed_theta0=0.30,
        seed_kappa=0.0,
        fixed_theta0=None,
        J=17,
        Nb=120,
        max_iter_eq=max_iter_eq,
        tol_eq=tol_eq,
    )
    _, fixed, mechanism = arm_contract(contract_args)
    theta = load_theta(args.results, "M5")
    theta.update(fixed)
    overrides = {
        **common_overrides(contract_args, mechanism),
        **theta,
        "max_iter_eq": max_iter_eq,
        "tol_eq": tol_eq,
        "use_numba_scatter": bool(args.numba_scatter),
    }
    if args.mode in {"fixed-price", "fixed-price-fast"}:
        overrides.update(
            solve_mode="pe",
            p_fixed=np.array([float(winner["price"])]),
            w_fixed=np.array([1.0]),
            entry_shares_fixed=np.array([1.0]),
        )
    return overrides


def run_profiled_model(
    mode: str,
    overrides: dict[str, Any],
) -> tuple[Any, Any, np.ndarray]:
    if mode != "fixed-price-fast":
        return run_model_cp_dt(overrides, verbose=False)

    original = solver_module.solve_markov_income_partial_equilibrium

    def fast_partial_equilibrium(
        p_fixed: np.ndarray,
        P: Any,
        b_grid: np.ndarray,
        verbose: bool = False,
    ) -> tuple[Any, Any, np.ndarray]:
        _ = verbose
        price = np.asarray(p_fixed, dtype=float).reshape(-1).copy()
        P.eq_iter = 1
        solution = solver_module.solve_markov_income_at_prices(
            price,
            P,
            b_grid,
            verbose=False,
            fast_stats=True,
        )
        return solution, P, price

    solver_module.solve_markov_income_partial_equilibrium = fast_partial_equilibrium
    try:
        return run_model_cp_dt(overrides, verbose=False)
    finally:
        solver_module.solve_markov_income_partial_equilibrium = original


def probability_summary(
    name: str,
    probs: np.ndarray,
    *,
    choice_axis: int,
) -> dict[str, Any]:
    array = np.asarray(probs, dtype=float)
    totals = np.sum(array, axis=choice_axis)
    active = totals > 1e-12
    return {
        "name": name,
        "shape": list(array.shape),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "active_row_count": int(np.sum(active)),
        "max_active_sum_error": (
            float(np.max(np.abs(totals[active] - 1.0))) if np.any(active) else 0.0
        ),
        "negative_count": int(np.sum(array < -1e-12)),
        "above_one_count": int(np.sum(array > 1.0 + 1e-12)),
    }


def invariant_summary(sol: Any, P: Any) -> dict[str, Any]:
    distribution = np.asarray(sol.g, dtype=float)
    values = np.asarray(sol.V, dtype=float)
    alive = values > -1e9
    adjacent_alive = alive[1:] & alive[:-1]
    value_diff = values[1:] - values[:-1]
    monotonicity_violations = adjacent_alive & (value_diff < -1e-8)
    economically_alive = (values[1:] > -1e6) & (values[:-1] > -1e6)
    material_monotonicity_violations = economically_alive & (value_diff < -1e-8)
    violation_mass = np.where(
        monotonicity_violations,
        distribution[1:] + distribution[:-1],
        0.0,
    )
    bp = np.asarray(sol.bp_pol, dtype=float)
    grid = np.asarray(sol.b_grid, dtype=float)
    return {
        "distribution": {
            "shape": list(distribution.shape),
            "mass": float(np.sum(distribution)),
            "minimum": float(np.min(distribution)),
            "negative_count": int(np.sum(distribution < -1e-14)),
            "nonfinite_count": int(np.sum(~np.isfinite(distribution))),
        },
        "values": {
            "shape": list(values.shape),
            "alive_nonfinite_count": int(np.sum(alive & ~np.isfinite(values))),
            "adjacent_alive_pairs": int(np.sum(adjacent_alive)),
            "wealth_monotonicity_violation_count": int(np.sum(monotonicity_violations)),
            "wealth_monotonicity_violation_share": float(
                np.sum(monotonicity_violations) / max(np.sum(adjacent_alive), 1)
            ),
            "material_wealth_monotonicity_violation_count": int(
                np.sum(material_monotonicity_violations)
            ),
            "mass_adjacent_to_wealth_monotonicity_violations": float(
                np.sum(violation_mass)
            ),
            "worst_wealth_value_drop": (
                float(np.min(value_diff[adjacent_alive])) if np.any(adjacent_alive) else 0.0
            ),
        },
        "savings_policy": {
            "minimum": float(np.min(bp)),
            "maximum": float(np.max(bp)),
            "grid_minimum": float(grid[0]),
            "grid_maximum": float(grid[-1]),
            "below_grid_count": int(np.sum(bp < grid[0] - 1e-10)),
            "above_grid_count": int(np.sum(bp > grid[-1] + 1e-10)),
            "above_grid_by_more_than_golden_tolerance_count": int(
                np.sum(bp > grid[-1] + 1e-3)
            ),
            "nonfinite_count": int(np.sum(~np.isfinite(bp))),
        },
        "probabilities": [
            probability_summary("tenure", sol.tenure_probs, choice_axis=-1),
            probability_summary("location", sol.loc_probs, choice_axis=3),
            probability_summary("fertility", sol.fert_probs, choice_axis=-1),
        ],
        "market": {
            "strict_converged": bool(getattr(sol, "converged", False)),
            "best_max_abs_rel_excess": float(
                getattr(sol, "best_max_abs_rel_excess", math.nan)
            ),
            "housing_demand": jsonable(sol.housing_demand),
            "housing_supply": jsonable(sol.housing_supply),
        },
        "configuration": {
            "J": int(P.J),
            "Nb": int(P.Nb),
            "income_states": int(P.Nz),
            "owner_rungs": int(P.n_house),
        },
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    winner = load_collected_winner(args.results)
    overrides = build_overrides(args, winner)
    tag = args.mode.replace("-", "_")
    if args.numba_scatter:
        tag += "_numba_scatter"
    configured_numba_threads = int(numba.get_num_threads())
    if configured_numba_threads != 1:
        tag += f"_threads{configured_numba_threads}"
    raw_path = args.outdir / f"{tag}.pstats"
    text_path = args.outdir / f"{tag}_profile.txt"
    summary_path = args.outdir / f"{tag}_summary.json"
    invariants_path = args.outdir / f"{tag}_invariants.json"

    load_before = os.getloadavg()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    profiler = cProfile.Profile()
    profiler.enable()
    sol, P, price = run_profiled_model(args.mode, overrides)
    profiler.disable()
    wall_seconds = time.perf_counter() - wall_start
    cpu_seconds = time.process_time() - cpu_start
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    load_after = os.getloadavg()
    profiler.dump_stats(raw_path)

    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats("cumulative")
    stats.print_stats(max(1, int(args.top)))
    text_path.write_text(stream.getvalue())

    targets, weights = target_system(TARGET_SET)
    if args.mode == "fixed-price-fast":
        loss = math.nan
        selected_moment_gaps: dict[str, float] = {}
    else:
        moments = extract_moments(sol, P)
        loss = diagnostic_loss(moments, targets=targets, weights=weights)
        selected_moment_gaps = {
            name: float(moments[name]) - float(row["model"])
            for name, row in ((row["moment"], row) for row in winner["target_fit"])
            if name in moments
        }
    summary = {
        "mode": args.mode,
        "results_source": str(args.results),
        "wall_seconds": wall_seconds,
        "cpu_seconds": cpu_seconds,
        "cpu_to_wall_ratio": cpu_seconds / max(wall_seconds, 1e-12),
        "max_rss_before": int(usage_before.ru_maxrss),
        "max_rss_after": int(usage_after.ru_maxrss),
        "load_average_before": list(load_before),
        "load_average_after": list(load_after),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "numba": numba.__version__,
        "numba_threads": configured_numba_threads,
        "numba_scatter": bool(args.numba_scatter),
        "thread_environment": {
            name: os.environ.get(name)
            for name in (
                "NUMBA_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            )
        },
        "price": jsonable(price),
        "loss": float(loss),
        "collected_loss": float(winner["rank_loss"]),
        "max_abs_collected_target_moment_gap": float(
            max((abs(value) for value in selected_moment_gaps.values()), default=0.0)
        ),
        "timings": jsonable(getattr(sol, "timings", {})),
        "profile_artifact": str(raw_path),
        "profile_text": str(text_path),
    }
    summary_path.write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    if args.mode != "fixed-price-fast":
        invariants_path.write_text(
            json.dumps(jsonable(invariant_summary(sol, P)), indent=2, sort_keys=True)
        )
    print(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
