"""Bounded multicore local diagnostics for the one-market model."""

from __future__ import annotations

import concurrent.futures as cf
import json
import math
import os
import time
from pathlib import Path
from typing import Any

import numpy as np

from .calibration import (
    OLD_NONLOCATION_TARGETS,
    OLD_NONLOCATION_WEIGHTS,
    PERIOD_YEARS,
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    informed_smoke_candidates,
    jsonable,
)
from .diagnostics import write_diagnostics
from .parameters import make_persistent_transition_matrix
from .production_profile import (
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    production_profile_metadata,
    production_profile_overrides,
    validate_production_profile,
)
from .solver import run_model_cp_dt


DEFAULT_INCOME_GRID_5 = np.array([0.60, 0.80, 1.00, 1.20, 1.40])
DEFAULT_INCOME_WEIGHTS_5 = np.array([0.10, 0.20, 0.40, 0.20, 0.10])
DEFAULT_RHO_Z = 0.85
DEFAULT_TENURE_CHOICE_KAPPA = 0.01

GLOBAL_DE_BOUNDS = list(PRODUCTION_SEARCH_BOUNDS)


def run_local_panel(
    outdir: Path,
    *,
    n_cases: int = 144,
    seed: int = 20260608,
    J: int = 17,
    Nb: int = 60,
    n_house: int = 6,
    max_iter_eq: int = 25,
    workers: int = 6,
    minutes: float = 30.0,
    income_states: int = 5,
    diagnostic_best: int = 3,
    target_set: str = "candidate_no_timing_v0",
    include_anchors: bool = True,
    seed_theta: dict[str, Any] | None = None,
    progress: bool = True,
) -> dict[str, Any]:
    """Run a bounded multicore diagnostic panel.

    This is not a formal calibration routine. It is a local economic check
    around the current production parameter ledger. The ranking loss excludes
    mean age at first birth because the current model has no direct timing
    shifter for that object.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    cases_path.write_text("")

    workers = max(1, int(workers))
    os.environ["NUMBA_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    rank_targets, rank_weights = get_target_set(target_set)
    income = income_process_overrides(income_states)
    seed_theta_clean = keep_internal_theta(seed_theta) if seed_theta is not None else None
    candidates = local_panel_candidates(
        n_cases,
        seed,
        include_anchors=include_anchors,
        seed_theta=seed_theta_clean,
    )
    meta = {
        "status": "bounded_multicore_diagnostic_not_formal_calibration",
        "n_cases_requested": int(n_cases),
        "seed": int(seed),
        "include_anchors": bool(include_anchors),
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "workers": int(workers),
        "minutes_budget": float(minutes),
        "income_states": int(income_states),
        "z_grid": jsonable(income["z_grid"]),
        "z_weights": jsonable(income["z_weights"]),
        "income_shock_persistence": float(DEFAULT_RHO_Z),
        "rank_target_set": str(target_set),
        "rank_targets": rank_targets,
        "rank_weights": rank_weights,
        "seed_theta": jsonable(seed_theta_clean),
        "full_old_targets": OLD_NONLOCATION_TARGETS,
        "full_old_weights": OLD_NONLOCATION_WEIGHTS,
        "varied_internal_parameters": [
            "beta",
            "alpha_cons",
            "c_bar_0",
            "c_bar_n",
            "h_bar_0",
            "h_bar_jump",
            "h_bar_n",
            "psi_child",
            "kappa_fert",
            "tenure_choice_kappa",
            "chi",
            "theta0",
            "theta_n",
        ],
        "fixed_external_or_first_stage": [
            "entry_wealth_mode",
            "entry_wealth_ratio_nodes",
            "entry_wealth_ratio_weights",
            "q",
            "delta",
            "tau_H",
            "phi",
            "pti_limit",
            "psi",
            "income_age_profile",
            "H_own",
            "hR_max",
            "H0",
            "r_bar",
            "eta_supply",
            "theta1",
        ],
        "ranking_note": (
            "candidate_no_timing_v0 excludes mean_age_first_birth, keeps parity composition diagnostic, "
            "and adds candidate targets for midlife liquid wealth/income, housing user-cost share, "
            "and childless renter/owner median rooms. This ledger is a documented trial target set, "
            "not a finalized empirical target system."
        ),
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))

    start = time.perf_counter()
    deadline = start + max(1.0, float(minutes) * 60.0)
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    next_idx = 0
    futures: dict[cf.Future[dict[str, Any]], int] = {}

    def submit_next(executor: cf.ProcessPoolExecutor) -> None:
        nonlocal next_idx
        if next_idx >= len(candidates) or time.perf_counter() >= deadline:
            return
        idx = next_idx
        next_idx += 1
        futures[
            executor.submit(
                run_local_panel_case,
                idx,
                candidates[idx],
                J,
                Nb,
                n_house,
                max_iter_eq,
                income,
                rank_targets,
                rank_weights,
            )
        ] = idx

    with cf.ProcessPoolExecutor(max_workers=workers, initializer=init_panel_worker) as executor:
        for _ in range(min(workers, len(candidates))):
            submit_next(executor)
        while futures:
            done, _ = cf.wait(futures, timeout=30.0, return_when=cf.FIRST_COMPLETED)
            if not done:
                if progress:
                    elapsed = time.perf_counter() - start
                    print(
                        f"heartbeat: completed={len(records)}, submitted={next_idx}, "
                        f"elapsed={elapsed / 60:.1f} min",
                        flush=True,
                    )
                continue
            for fut in sorted(done, key=lambda f: futures[f]):
                futures.pop(fut)
                try:
                    record = fut.result()
                except Exception as exc:  # noqa: BLE001 - checkpoint case-level failures.
                    record = {
                        "case": -1,
                        "label": "executor_failure",
                        "status": f"failed: {type(exc).__name__}: {exc}",
                        "rank_loss": math.inf,
                        "full_old_nonlocation_loss": math.inf,
                        "theta": {},
                        "moments": {},
                        "p_eq": [math.nan],
                        "market_residual": math.inf,
                        "elapsed_sec": math.nan,
                        "timings": {},
                    }
                records.append(record)
                append_jsonl(cases_path, record)
                if is_better_record(record, best):
                    best = record
                    best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
                if progress:
                    best_loss = float(best["rank_loss"]) if best is not None else math.inf
                    elapsed = time.perf_counter() - start
                    print(
                        f"case {record.get('case')} {record.get('label')}: "
                        f"rank={float(record.get('rank_loss', math.inf)):.4g}, "
                        f"full={float(record.get('full_old_nonlocation_loss', math.inf)):.4g}, "
                        f"resid={float(record.get('market_residual', math.inf)):.2e}, "
                        f"case_sec={float(record.get('elapsed_sec', math.nan)):.1f}, "
                        f"best={best_loss:.4g}, done={len(records)}, "
                        f"elapsed={elapsed / 60:.1f}m",
                        flush=True,
                    )
                submit_next(executor)

    records_sorted = sorted(records, key=record_sort_key)
    converged_records = strict_records(records_sorted)
    summary = {
        "best": converged_records[0] if converged_records else None,
        "top_records": converged_records[: min(10, len(converged_records))],
        "elapsed_sec": float(time.perf_counter() - start),
        "n_cases_completed": int(len(records)),
        "n_cases_submitted": int(next_idx),
        "stopped_by_time_budget": bool(next_idx < len(candidates)),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    write_panel_summary_tables(outdir, records_sorted, rank_targets)
    write_panel_plots(outdir, records_sorted, rank_targets)
    if diagnostic_best > 0 and converged_records:
        write_best_case_diagnostics(
            outdir,
            converged_records[: int(diagnostic_best)],
            J=J,
            Nb=Nb,
            n_house=n_house,
            max_iter_eq=max_iter_eq,
            income=income,
        )
    return summary


def run_global_de_panel(
    outdir: Path,
    *,
    max_evals: int = 240,
    seed: int = 20260609,
    J: int = 16,
    Nb: int = 60,
    n_house: int = 6,
    max_iter_eq: int = 25,
    minutes: float = 115.0,
    income_states: int = 5,
    pop_size: int = 20,
    mutation: float = 0.85,
    crossover: float = 0.70,
    target_set: str = "candidate_no_timing_v0",
    seed_theta: dict[str, Any] | None = None,
    profile_name: str | None = None,
    progress: bool = True,
) -> dict[str, Any]:
    """Run an independent differential-evolution global search panel.

    Each array task runs a full restart from broad bounds over the 14-parameter
    internal calibration vector. This changes the search algorithm only; it
    uses the same model solution and ranking target system as `run_local_panel`.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    cases_path.write_text("")

    os.environ["NUMBA_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    rng = np.random.default_rng(seed)
    rank_targets, rank_weights = get_target_set(target_set)
    income = income_process_overrides(income_states)
    profile_extra_overrides: dict[str, Any] = {}
    if profile_name is not None:
        validate_production_profile(
            profile_name,
            J=J,
            Nb=Nb,
            n_house=n_house,
            income_states=income_states,
            target_set=target_set,
            max_iter_eq=max_iter_eq,
            stage="search",
        )
        profile_extra_overrides = production_profile_overrides()
    dim = len(GLOBAL_DE_BOUNDS)
    seed_theta_clean = keep_internal_theta(seed_theta) if seed_theta is not None else None
    pop_size = max(4, int(pop_size))
    max_evals = max(1, int(max_evals))
    mutation = float(mutation)
    crossover = min(max(float(crossover), 0.0), 1.0)

    meta = {
        "status": "global_de_panel_not_formal_calibration",
        "algorithm": "independent_latin_hypercube_plus_differential_evolution_rand_1_bin",
        "seed": int(seed),
        "max_evals_requested": int(max_evals),
        "pop_size": int(pop_size),
        "mutation": float(mutation),
        "crossover": float(crossover),
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "minutes_budget": float(minutes),
        "income_states": int(income_states),
        "z_grid": jsonable(income["z_grid"]),
        "z_weights": jsonable(income["z_weights"]),
        "income_shock_persistence": float(DEFAULT_RHO_Z),
        "rank_target_set": str(target_set),
        "rank_targets": rank_targets,
        "rank_weights": rank_weights,
        "production_profile": str(profile_name) if profile_name is not None else None,
        "production_profile_spec": production_profile_metadata() if profile_name is not None else None,
        "seed_theta": jsonable(seed_theta_clean),
        "bounds": [
            {"name": name, "lower": float(lo), "upper": float(hi)}
            for name, lo, hi in GLOBAL_DE_BOUNDS
        ],
        "bounds_note": (
            "beta_annual is transformed to the four-year beta used by the model; "
            "all other bounds are on model-period parameters."
        ),
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))

    start = time.perf_counter()
    deadline = start + max(1.0, float(minutes) * 60.0)
    pop = latin_hypercube(rng, pop_size, dim)
    seeded_unit = global_unit_from_theta(seed_theta_clean) if seed_theta_clean is not None else None
    if seeded_unit is not None:
        pop[0] = seeded_unit
    pop_loss = np.full(pop_size, math.inf)
    pop_records: list[dict[str, Any] | None] = [None] * pop_size
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    eval_idx = 0
    generation = 0

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal eval_idx, best
        theta = theta_from_global_unit(unit)
        record = run_local_panel_case(
            eval_idx,
            {"label": label, "theta": theta},
            J,
            Nb,
            n_house,
            max_iter_eq,
            income,
            rank_targets,
            rank_weights,
            profile_extra_overrides,
        )
        record["algorithm"] = "global_de"
        record["origin"] = jsonable(origin)
        record["unit_vector"] = jsonable(unit)
        append_jsonl(cases_path, record)
        records.append(record)
        if is_better_record(record, best):
            best = record
            best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
        if progress:
            best_loss = float(best["rank_loss"]) if best is not None else math.inf
            elapsed = time.perf_counter() - start
            print(
                f"eval {eval_idx + 1}/{max_evals} {label}: "
                f"rank={float(record.get('rank_loss', math.inf)):.4g}, "
                f"full={float(record.get('full_old_nonlocation_loss', math.inf)):.4g}, "
                f"resid={float(record.get('market_residual', math.inf)):.2e}, "
                f"case_sec={float(record.get('elapsed_sec', math.nan)):.1f}, "
                f"best={best_loss:.4g}, elapsed={elapsed / 60:.1f}m",
                flush=True,
            )
        eval_idx += 1
        return record

    for i in range(pop_size):
        if eval_idx >= max_evals or time.perf_counter() >= deadline:
            break
        label = "warm_start" if i == 0 and seeded_unit is not None else f"init_{i:03d}"
        phase = "seed_theta" if i == 0 and seeded_unit is not None else "latin_hypercube"
        record = evaluate(pop[i], label, {"phase": phase, "slot": i})
        pop_loss[i] = record_selection_loss(record)
        pop_records[i] = record

    while eval_idx < max_evals and time.perf_counter() < deadline:
        order = np.arange(pop_size)
        rng.shuffle(order)
        for i in order:
            if eval_idx >= max_evals or time.perf_counter() >= deadline:
                break
            choices = [j for j in range(pop_size) if j != int(i)]
            if len(choices) < 3:
                break
            r1, r2, r3 = rng.choice(choices, size=3, replace=False)
            mutant = pop[r1] + mutation * (pop[r2] - pop[r3])
            out_of_bounds = (mutant < 0.0) | (mutant > 1.0) | ~np.isfinite(mutant)
            if np.any(out_of_bounds):
                mutant = mutant.copy()
                mutant[out_of_bounds] = rng.random(int(np.sum(out_of_bounds)))
            trial = pop[i].copy()
            mask = rng.random(dim) < crossover
            mask[rng.integers(0, dim)] = True
            trial[mask] = np.clip(mutant[mask], 0.0, 1.0)
            record = evaluate(
                trial,
                f"de_g{generation:03d}_i{int(i):03d}",
                {
                    "phase": "differential_evolution_rand_1_bin",
                    "generation": int(generation),
                    "slot": int(i),
                    "base": int(r1),
                    "diff_a": int(r2),
                    "diff_b": int(r3),
                },
            )
            trial_loss = record_selection_loss(record)
            if trial_loss <= pop_loss[i]:
                pop[i] = trial
                pop_loss[i] = trial_loss
                pop_records[i] = record

        generation += 1
        if eval_idx < max_evals and time.perf_counter() < deadline and generation % 3 == 0:
            worst = int(np.argmax(pop_loss))
            immigrant = rng.random(dim)
            record = evaluate(
                immigrant,
                f"immigrant_g{generation:03d}_i{worst:03d}",
                {"phase": "random_immigrant", "generation": int(generation), "slot": worst},
            )
            immigrant_loss = record_selection_loss(record)
            if immigrant_loss <= pop_loss[worst]:
                pop[worst] = immigrant
                pop_loss[worst] = immigrant_loss
                pop_records[worst] = record

    records_sorted = sorted(records, key=record_sort_key)
    converged_records = strict_records(records_sorted)
    summary = {
        "best": converged_records[0] if converged_records else None,
        "top_records": converged_records[: min(10, len(converged_records))],
        "elapsed_sec": float(time.perf_counter() - start),
        "n_cases_completed": int(len(records)),
        "n_cases_submitted": int(eval_idx),
        "stopped_by_time_budget": bool(eval_idx < max_evals),
        "generations_completed": int(generation),
        "final_population_losses": jsonable(pop_loss),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    write_panel_summary_tables(outdir, records_sorted, rank_targets)
    write_panel_plots(outdir, records_sorted, rank_targets)
    return summary


def run_local_polish(
    outdir: Path,
    *,
    max_evals: int = 240,
    seed: int = 20260628,
    J: int = 17,
    Nb: int = 60,
    n_house: int = 6,
    max_iter_eq: int = PRODUCTION_MAX_ITER_EQ,
    minutes: float = 115.0,
    income_states: int = 5,
    target_set: str = "candidate_no_timing_v0",
    seed_theta: dict[str, Any] | None = None,
    method: str = "nelder-mead",
    initial_step: float = 0.06,
    min_step: float = 0.003,
    shrink: float = 0.5,
    profile_name: str | None = None,
    fixed_theta: dict[str, Any] | None = None,
    progress: bool = True,
) -> dict[str, Any]:
    """Run a true derivative-free local polish around a seed theta.

    The existing `local-panel` routine is a seeded/random panel, not a local
    optimizer. This routine optimizes in the normalized bounded parameter cube
    used by the global-DE code and checkpoints every model evaluation.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    cases_path.write_text("")

    os.environ["NUMBA_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    rng = np.random.default_rng(seed)
    rank_targets, rank_weights = get_target_set(target_set)
    income = income_process_overrides(income_states)
    profile_extra_overrides: dict[str, Any] = {}
    if profile_name is not None:
        validate_production_profile(
            profile_name,
            J=J,
            Nb=Nb,
            n_house=n_house,
            income_states=income_states,
            target_set=target_set,
            max_iter_eq=max_iter_eq,
            stage="search",
        )
        profile_extra_overrides = production_profile_overrides()
    fixed_theta_clean = keep_internal_theta(fixed_theta) if fixed_theta is not None else {}
    fixed_keys = set(fixed_theta_clean or {})
    search_bounds = [
        bound
        for bound in GLOBAL_DE_BOUNDS
        if ("beta" if bound[0] == "beta_annual" else bound[0]) not in fixed_keys
    ]
    seed_theta_clean = keep_internal_theta(seed_theta) if seed_theta is not None else None
    x0 = global_unit_from_theta(seed_theta_clean, bounds=search_bounds)
    if x0 is None:
        raise ValueError("local polish requires a seed theta containing all searched parameters.")

    method = str(method).strip().lower().replace("_", "-")
    if method not in {"nelder-mead", "pattern"}:
        raise ValueError(f"unknown local polish method: {method}")
    max_evals = max(1, int(max_evals))
    initial_step = min(max(float(initial_step), 1e-5), 0.75)
    min_step = min(max(float(min_step), 1e-6), initial_step)
    shrink = min(max(float(shrink), 0.05), 0.95)
    dim = len(search_bounds)

    meta = {
        "status": "local_polish_true_derivative_free_optimizer",
        "algorithm": method,
        "seed": int(seed),
        "max_evals_requested": int(max_evals),
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "minutes_budget": float(minutes),
        "income_states": int(income_states),
        "z_grid": jsonable(income["z_grid"]),
        "z_weights": jsonable(income["z_weights"]),
        "income_shock_persistence": float(DEFAULT_RHO_Z),
        "rank_target_set": str(target_set),
        "rank_targets": rank_targets,
        "rank_weights": rank_weights,
        "production_profile": str(profile_name) if profile_name is not None else None,
        "production_profile_spec": production_profile_metadata() if profile_name is not None else None,
        "fixed_theta": jsonable(fixed_theta_clean),
        "seed_theta": jsonable(seed_theta_clean),
        "initial_unit_vector": jsonable(x0),
        "initial_step": float(initial_step),
        "min_step": float(min_step),
        "shrink": float(shrink),
        "bounds": [
            {"name": name, "lower": float(lo), "upper": float(hi)}
            for name, lo, hi in search_bounds
        ],
        "source_controlled_bounds": [
            {"name": name, "lower": float(lo), "upper": float(hi)}
            for name, lo, hi in GLOBAL_DE_BOUNDS
        ],
        "bounds_note": (
            "Local polish is performed in the same normalized bounded cube as "
            "global-DE; beta_annual is transformed to four-year beta."
        ),
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))

    start = time.perf_counter()
    deadline = start + max(1.0, float(minutes) * 60.0)
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    eval_idx = 0

    def loss_of(record: dict[str, Any]) -> float:
        return record_selection_loss(record)

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal eval_idx, best
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = theta_from_global_unit(clipped, bounds=search_bounds)
        theta.update(fixed_theta_clean or {})
        record = run_local_panel_case(
            eval_idx,
            {"label": label, "theta": theta},
            J,
            Nb,
            n_house,
            max_iter_eq,
            income,
            rank_targets,
            rank_weights,
            profile_extra_overrides,
        )
        record["algorithm"] = f"local_polish_{method}"
        record["origin"] = jsonable(origin)
        record["unit_vector"] = jsonable(clipped)
        append_jsonl(cases_path, record)
        records.append(record)
        if is_better_record(record, best):
            best = record
            best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
        if progress:
            best_loss = loss_of(best) if best is not None else math.inf
            elapsed = time.perf_counter() - start
            print(
                f"eval {eval_idx + 1}/{max_evals} {label}: "
                f"rank={loss_of(record):.4g}, "
                f"full={float(record.get('full_old_nonlocation_loss', math.inf)):.4g}, "
                f"resid={float(record.get('market_residual', math.inf)):.2e}, "
                f"case_sec={float(record.get('elapsed_sec', math.nan)):.1f}, "
                f"best={best_loss:.4g}, elapsed={elapsed / 60:.1f}m",
                flush=True,
            )
        eval_idx += 1
        return record

    def can_eval() -> bool:
        return eval_idx < max_evals and time.perf_counter() < deadline

    if method == "pattern":
        best_unit = x0.copy()
        best_record = evaluate(best_unit, "seed", {"phase": "seed"})
        step = initial_step
        iteration = 0
        while can_eval() and step >= min_step:
            improved = False
            order = np.arange(dim)
            rng.shuffle(order)
            for j in order:
                if not can_eval():
                    break
                candidates: list[tuple[float, np.ndarray, dict[str, Any]]] = [
                    (loss_of(best_record), best_unit, best_record)
                ]
                for direction in (1.0, -1.0):
                    if not can_eval():
                        break
                    trial = best_unit.copy()
                    trial[int(j)] = np.clip(trial[int(j)] + direction * step, 0.0, 1.0)
                    rec = evaluate(
                        trial,
                        f"pattern_i{iteration:03d}_d{int(j):02d}_{'p' if direction > 0 else 'm'}",
                        {"phase": "coordinate_probe", "iteration": iteration, "dim": int(j), "step": step, "direction": direction},
                    )
                    candidates.append((loss_of(rec), trial, rec))
                candidates.sort(key=lambda item: item[0])
                if candidates[0][0] + 1e-10 < loss_of(best_record):
                    best_unit = candidates[0][1].copy()
                    best_record = candidates[0][2]
                    improved = True
            if not improved:
                step *= shrink
            iteration += 1
    else:
        simplex: list[tuple[np.ndarray, dict[str, Any]]] = []
        seed_record = evaluate(x0, "seed", {"phase": "seed"})
        simplex.append((x0.copy(), seed_record))
        for j in range(dim):
            if not can_eval():
                break
            step_j = initial_step * (0.75 + 0.5 * rng.random())
            trial = x0.copy()
            if trial[j] + step_j <= 1.0:
                trial[j] += step_j
            else:
                trial[j] = max(0.0, trial[j] - step_j)
            rec = evaluate(trial, f"nm_init_{j:02d}", {"phase": "initial_simplex", "dim": int(j), "step": float(step_j)})
            simplex.append((trial, rec))

        alpha = 1.0
        gamma = 2.0
        rho = 0.5
        sigma = shrink
        iteration = 0
        while can_eval() and len(simplex) >= 2:
            simplex.sort(key=lambda item: loss_of(item[1]))
            if len(simplex) < dim + 1 and not can_eval():
                break
            n_centroid = max(1, len(simplex) - 1)
            centroid = np.mean([u for u, _ in simplex[:n_centroid]], axis=0)
            worst_unit, worst_record = simplex[-1]
            second_worst_loss = loss_of(simplex[-2][1]) if len(simplex) > 1 else math.inf
            best_loss = loss_of(simplex[0][1])

            reflected = np.clip(centroid + alpha * (centroid - worst_unit), 0.0, 1.0)
            reflected_record = evaluate(reflected, f"nm_reflect_{iteration:03d}", {"phase": "reflect", "iteration": iteration})
            reflected_loss = loss_of(reflected_record)

            if reflected_loss < best_loss and can_eval():
                expanded = np.clip(centroid + gamma * (reflected - centroid), 0.0, 1.0)
                expanded_record = evaluate(expanded, f"nm_expand_{iteration:03d}", {"phase": "expand", "iteration": iteration})
                simplex[-1] = (expanded, expanded_record) if loss_of(expanded_record) < reflected_loss else (reflected, reflected_record)
            elif reflected_loss < second_worst_loss:
                simplex[-1] = (reflected, reflected_record)
            else:
                if reflected_loss < loss_of(worst_record):
                    contracted = np.clip(centroid + rho * (reflected - centroid), 0.0, 1.0)
                    phase = "outside_contract"
                    threshold = reflected_loss
                else:
                    contracted = np.clip(centroid - rho * (centroid - worst_unit), 0.0, 1.0)
                    phase = "inside_contract"
                    threshold = loss_of(worst_record)
                contract_record = evaluate(contracted, f"nm_contract_{iteration:03d}", {"phase": phase, "iteration": iteration}) if can_eval() else worst_record
                if loss_of(contract_record) < threshold:
                    simplex[-1] = (contracted, contract_record)
                else:
                    best_unit = simplex[0][0].copy()
                    new_simplex = [simplex[0]]
                    for k, (unit, _) in enumerate(simplex[1:], start=1):
                        if not can_eval():
                            new_simplex.append((unit, simplex[k][1]))
                            continue
                        shrunk = np.clip(best_unit + sigma * (unit - best_unit), 0.0, 1.0)
                        shrunk_record = evaluate(shrunk, f"nm_shrink_{iteration:03d}_{k:02d}", {"phase": "shrink", "iteration": iteration, "slot": k})
                        new_simplex.append((shrunk, shrunk_record))
                    simplex = new_simplex
            iteration += 1

    records_sorted = sorted(records, key=record_sort_key)
    converged_records = strict_records(records_sorted)
    summary = {
        "best": converged_records[0] if converged_records else None,
        "top_records": converged_records[: min(10, len(converged_records))],
        "elapsed_sec": float(time.perf_counter() - start),
        "n_cases_completed": int(len(records)),
        "n_cases_submitted": int(eval_idx),
        "stopped_by_time_budget": bool(eval_idx < max_evals),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    write_panel_summary_tables(outdir, records_sorted, rank_targets)
    write_panel_plots(outdir, records_sorted, rank_targets)
    return summary


def init_panel_worker() -> None:
    try:
        from numba import set_num_threads

        set_num_threads(1)
    except Exception:
        pass


def run_local_panel_case(
    idx: int,
    candidate: dict[str, Any],
    J: int,
    Nb: int,
    n_house: int,
    max_iter_eq: int,
    income: dict[str, Any],
    rank_targets: dict[str, float],
    rank_weights: dict[str, float],
    extra_overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    theta = dict(candidate["theta"])
    overrides = {
        **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
        **income,
        **(extra_overrides or {}),
        **theta,
    }
    try:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
        moments = extract_moments(sol, P)
        rank_loss = diagnostic_loss(moments, targets=rank_targets, weights=rank_weights)
        full_loss = diagnostic_loss(moments, targets=OLD_NONLOCATION_TARGETS, weights=OLD_NONLOCATION_WEIGHTS)
        status = "ok"
        market_residual = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
        timings = getattr(sol, "timings", {})
        solver_strict = bool(timings.get("strict_converged", getattr(sol, "converged", False)))
        strict_converged = bool(
            solver_strict
            and np.isfinite(market_residual)
            and market_residual <= float(getattr(P, "tol_eq", 1e-4))
        )
        tol_eq = float(getattr(P, "tol_eq", 1e-4))
    except Exception as exc:  # noqa: BLE001 - panel should checkpoint failed parameter vectors.
        moments = {}
        rank_loss = math.inf
        full_loss = math.inf
        status = f"failed: {type(exc).__name__}: {exc}"
        p_eq = np.array([np.nan])
        market_residual = math.inf
        timings = {}
        strict_converged = False
        tol_eq = math.nan
    return {
        "case": int(idx),
        "label": str(candidate["label"]),
        "status": status,
        "rank_loss": float(rank_loss),
        "full_old_nonlocation_loss": float(full_loss),
        "theta": jsonable(theta),
        "moments": jsonable(moments),
        "p_eq": jsonable(p_eq),
        "market_residual": float(market_residual),
        "strict_converged": bool(strict_converged),
        "tol_eq": float(tol_eq),
        "elapsed_sec": float(time.perf_counter() - t0),
        "timings": jsonable(timings),
    }


def local_panel_candidates(
    n_cases: int,
    seed: int,
    *,
    include_anchors: bool = True,
    seed_theta: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    if seed_theta is not None:
        candidates.append({"label": "warm_start", "theta": dict(seed_theta)})
    if include_anchors:
        anchors = [
            keep_internal_candidate(candidate)
            for candidate in informed_smoke_candidates()
            if str(candidate["label"]) != "baseline"
        ]
        candidates.append({"label": "baseline", "theta": {}})
        candidates.extend(anchors)
    rng = np.random.default_rng(seed)
    while len(candidates) < int(n_cases):
        idx = len(candidates)
        candidates.append({"label": f"draw_{idx:04d}", "theta": draw_internal_candidate(rng, idx)})
    return candidates[: int(n_cases)]


def latin_hypercube(rng: np.random.Generator, n: int, dim: int) -> np.ndarray:
    out = np.empty((int(n), int(dim)))
    for j in range(int(dim)):
        out[:, j] = (np.arange(int(n)) + rng.random(int(n))) / float(n)
        rng.shuffle(out[:, j])
    return out


def theta_from_global_unit(
    unit: np.ndarray,
    bounds: list[tuple[str, float, float]] | tuple[tuple[str, float, float], ...] | None = None,
) -> dict[str, float]:
    unit = np.asarray(unit, dtype=float)
    active_bounds = GLOBAL_DE_BOUNDS if bounds is None else list(bounds)
    theta: dict[str, float] = {}
    for u, (name, lo, hi) in zip(unit, active_bounds):
        value = float(lo + np.clip(u, 0.0, 1.0) * (hi - lo))
        if name == "beta_annual":
            theta["beta"] = value**PERIOD_YEARS
        else:
            theta[name] = value
    return theta


def global_unit_from_theta(
    theta: dict[str, Any] | None,
    bounds: list[tuple[str, float, float]] | tuple[tuple[str, float, float], ...] | None = None,
) -> np.ndarray | None:
    if theta is None:
        return None
    active_bounds = GLOBAL_DE_BOUNDS if bounds is None else list(bounds)
    unit = np.full(len(active_bounds), np.nan)
    for idx, (name, lo, hi) in enumerate(active_bounds):
        source_name = "beta" if name == "beta_annual" else name
        if source_name not in theta:
            if source_name == "tenure_choice_kappa":
                value = DEFAULT_TENURE_CHOICE_KAPPA
            else:
                return None
        else:
            value = float(theta[source_name])
        if name == "beta_annual":
            value = value ** (1.0 / PERIOD_YEARS)
        unit[idx] = (value - float(lo)) / max(float(hi - lo), 1e-12)
    return np.clip(unit, 0.0, 1.0)


def keep_internal_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    allowed = {
        "beta",
        "alpha_cons",
        "c_bar_0",
        "c_bar_n",
        "h_bar_0",
        "h_bar_jump",
        "h_bar_n",
        "psi_child",
        "kappa_fert",
        "tenure_choice_kappa",
        "chi",
        "theta0",
        "theta_n",
    }
    theta = {k: v for k, v in dict(candidate["theta"]).items() if k in allowed}
    return {"label": str(candidate["label"]), "theta": theta}


def keep_internal_theta(theta: dict[str, Any] | None) -> dict[str, Any] | None:
    if theta is None:
        return None
    return keep_internal_candidate({"label": "seed_theta", "theta": theta})["theta"]


def draw_internal_candidate(rng: np.random.Generator, idx: int) -> dict[str, Any]:
    beta_annual = rng.uniform(0.94, 0.985)
    mode = idx % 4
    if mode == 0:
        alpha_lo, alpha_hi = 0.56, 0.70
        child_cost_lo, child_cost_hi = 0.22, 0.55
        psi_child_lo, psi_child_hi = 0.07, 0.16
    elif mode == 1:
        alpha_lo, alpha_hi = 0.62, 0.78
        child_cost_lo, child_cost_hi = 0.35, 0.75
        psi_child_lo, psi_child_hi = 0.04, 0.12
    elif mode == 2:
        alpha_lo, alpha_hi = 0.58, 0.75
        child_cost_lo, child_cost_hi = 0.25, 0.65
        psi_child_lo, psi_child_hi = 0.05, 0.14
    else:
        alpha_lo, alpha_hi = 0.60, 0.80
        child_cost_lo, child_cost_hi = 0.25, 0.75
        psi_child_lo, psi_child_hi = 0.03, 0.14
    return {
        "beta": beta_annual**PERIOD_YEARS,
        "alpha_cons": rng.uniform(alpha_lo, alpha_hi),
        "c_bar_0": rng.uniform(0.04, 0.16) * PERIOD_YEARS,
        "c_bar_n": rng.uniform(child_cost_lo, child_cost_hi),
        "h_bar_0": rng.uniform(3.0, 4.8),
        "h_bar_jump": rng.uniform(0.40, 1.35),
        "h_bar_n": rng.uniform(0.30, 1.10),
        "psi_child": rng.uniform(psi_child_lo, psi_child_hi),
        "kappa_fert": rng.uniform(3.5, 7.5),
        "tenure_choice_kappa": rng.uniform(0.0, 0.06),
        "chi": rng.uniform(0.90, 1.35),
        "theta0": rng.uniform(0.20, 1.10),
        "theta_n": rng.uniform(0.00, 0.65),
    }


def income_process_overrides(income_states: int) -> dict[str, Any]:
    if int(income_states) == 5:
        z_grid = DEFAULT_INCOME_GRID_5.copy()
        z_weights = DEFAULT_INCOME_WEIGHTS_5.copy()
    elif int(income_states) == 3:
        z_grid = np.array([0.70, 1.00, 1.30])
        z_weights = np.array([0.30, 0.40, 0.30])
    else:
        n = max(1, int(income_states))
        z_grid = np.linspace(0.60, 1.40, n)
        center = 0.5 * (n - 1)
        z_weights = np.maximum(1.0 - np.abs(np.arange(n) - center) / max(center + 1.0, 1.0), 0.05)
        z_weights = z_weights / np.sum(z_weights)
    return {
        "use_income_types": True,
        "income_type_transition": "markov",
        "z_grid": z_grid,
        "z_weights": z_weights,
        "income_shock_persistence": DEFAULT_RHO_Z,
        "Pi_z": make_persistent_transition_matrix(z_weights, DEFAULT_RHO_Z),
    }


def ranking_targets_without_age() -> tuple[dict[str, float], dict[str, float]]:
    """Backward-compatible alias for the old no-timing diagnostic target set."""

    return get_target_set("old_nonlocation_no_timing")


def is_better_record(record: dict[str, Any], best: dict[str, Any] | None) -> bool:
    candidate_loss = record_selection_loss(record)
    if not np.isfinite(candidate_loss):
        return False
    if best is None:
        return True
    return candidate_loss < record_selection_loss(best)


def record_is_strictly_converged(record: dict[str, Any]) -> bool:
    return bool(record.get("strict_converged", False))


def record_selection_loss(record: dict[str, Any]) -> float:
    if not record_is_strictly_converged(record):
        return math.inf
    loss = float(record.get("rank_loss", math.inf))
    return loss if np.isfinite(loss) else math.inf


def record_sort_key(record: dict[str, Any]) -> tuple[float, float, int]:
    return (
        record_selection_loss(record),
        float(record.get("rank_loss", math.inf)),
        int(record.get("case", -1)),
    )


def strict_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [record for record in records if np.isfinite(record_selection_loss(record))]


def append_jsonl(path: Path, record: dict[str, Any]) -> None:
    with path.open("a") as fh:
        fh.write(json.dumps(jsonable(record), sort_keys=True) + "\n")


def write_panel_summary_tables(outdir: Path, records: list[dict[str, Any]], targets: dict[str, float]) -> None:
    rows = []
    for record in records:
        moments = dict(record.get("moments", {}))
        row = {
            "case": record.get("case"),
            "label": record.get("label"),
            "rank_loss": record.get("rank_loss"),
            "full_old_nonlocation_loss": record.get("full_old_nonlocation_loss"),
            "market_residual": record.get("market_residual"),
            "strict_converged": record.get("strict_converged"),
            "elapsed_sec": record.get("elapsed_sec"),
            "p_eq": first_scalar(record.get("p_eq")),
        }
        for name in targets:
            row[name] = moments.get(name)
            row[f"{name}_target"] = targets[name]
            row[f"{name}_gap"] = none_safe_sub(moments.get(name), targets[name])
        row["mean_age_first_birth"] = moments.get("mean_age_first_birth")
        row["parity_share_0"] = moments.get("parity_share_0")
        row["parity_share_1"] = moments.get("parity_share_1")
        row["parity_share_2plus"] = moments.get("parity_share_2plus")
        row["young_liquid_wealth_to_income"] = moments.get("young_liquid_wealth_to_income")
        row["young_childless_renter_liquid_wealth_to_annual_gross_income_2535"] = moments.get(
            "young_childless_renter_liquid_wealth_to_annual_gross_income_2535"
        )
        rows.append(row)
    write_csv(outdir / "ranked_cases.csv", rows)


def write_panel_plots(outdir: Path, records: list[dict[str, Any]], targets: dict[str, float]) -> None:
    if not records:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    ok = strict_records(records)
    if not ok:
        return
    best = ok[0]
    plot_best_moment_fit(figdir / "best_moment_fit.png", best, targets, plt)
    plot_tradeoffs(figdir / "panel_tradeoffs.png", ok, targets, plt)
    plot_loss_distribution(figdir / "loss_distribution.png", ok, plt)


def plot_best_moment_fit(path: Path, best: dict[str, Any], targets: dict[str, float], plt: Any) -> None:
    names = list(targets)
    values = [float(dict(best.get("moments", {})).get(name, np.nan)) for name in names]
    target_values = [targets[name] for name in names]
    x = np.arange(len(names))
    width = 0.38
    fig, ax = plt.subplots(figsize=(13, 5.2))
    ax.bar(x - width / 2, values, width, label="model")
    ax.bar(x + width / 2, target_values, width, label="target")
    ax.set_xticks(x, names, rotation=35, ha="right")
    ax.set_title(f"Best local-panel fit: {best.get('label')} case {best.get('case')}")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_tradeoffs(path: Path, records: list[dict[str, Any]], targets: dict[str, float], plt: Any) -> None:
    loss = np.array([float(r.get("rank_loss", np.nan)) for r in records])
    colors = np.log10(np.maximum(loss, 1e-8))
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    panels = [
        ("tfr", "childless_rate", "TFR vs childlessness"),
        (
            "own_rate",
            "young_childless_renter_liquid_wealth_to_annual_gross_income_2535"
            if "young_childless_renter_liquid_wealth_to_annual_gross_income_2535" in targets
            else "young_liquid_wealth_to_income",
            "Ownership vs young liquid wealth/income",
        ),
        ("housing_increment_0to1", "housing_increment_1to2", "Housing increments"),
        ("old_age_own_rate", "old_age_parent_childless_gap", "Old ownership and old parent-childless gap"),
    ]
    for ax, (xname, yname, title) in zip(axes.ravel(), panels):
        x = np.array([float(dict(r.get("moments", {})).get(xname, np.nan)) for r in records])
        y = np.array([float(dict(r.get("moments", {})).get(yname, np.nan)) for r in records])
        sc = ax.scatter(x, y, c=colors, s=26, cmap="viridis_r", alpha=0.85)
        if xname in targets:
            ax.axvline(targets[xname], color="0.35", linestyle="--", linewidth=1)
        if yname in targets:
            ax.axhline(targets[yname], color="0.35", linestyle="--", linewidth=1)
        ax.set_xlabel(xname)
        ax.set_ylabel(yname)
        ax.set_title(title)
        ax.grid(alpha=0.2)
    fig.colorbar(sc, ax=axes.ravel().tolist(), label="log10 rank loss")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_loss_distribution(path: Path, records: list[dict[str, Any]], plt: Any) -> None:
    loss = np.array([float(r.get("rank_loss", np.nan)) for r in records])
    loss = loss[np.isfinite(loss)]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(loss, bins=min(30, max(5, len(loss) // 4)), color="tab:blue", alpha=0.85)
    ax.set_xlabel("rank loss")
    ax.set_ylabel("case count")
    ax.set_title("Local-panel loss distribution")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_best_case_diagnostics(
    outdir: Path,
    records: list[dict[str, Any]],
    *,
    J: int,
    Nb: int,
    n_house: int,
    max_iter_eq: int,
    income: dict[str, Any],
) -> None:
    diag_root = outdir / "diagnostics"
    diag_root.mkdir(parents=True, exist_ok=True)
    for rank, record in enumerate(records, start=1):
        theta = dict(record.get("theta", {}))
        overrides = {
            **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
            **income,
            **theta,
        }
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
        case_dir = diag_root / f"rank{rank:02d}_case{int(record.get('case', -1)):04d}_{record.get('label')}"
        write_diagnostics(sol, P, case_dir)
        moments = extract_moments(sol, P)
        payload = {
            "rank": rank,
            "source_record": record,
            "p_eq": jsonable(p_eq),
            "moments": jsonable(moments),
            "timings": jsonable(getattr(sol, "timings", {})),
        }
        (case_dir / "record.json").write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))


def load_panel_records(path: Path) -> list[dict[str, Any]]:
    records = []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    columns = list(rows[0].keys())
    with path.open("w") as fh:
        fh.write(",".join(columns) + "\n")
        for row in rows:
            fh.write(",".join(csv_cell(row.get(col)) for col in columns) + "\n")


def csv_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not np.isfinite(value):
            return ""
        return f"{value:.10g}"
    text = str(value)
    if any(ch in text for ch in [",", "\n", '"']):
        text = '"' + text.replace('"', '""') + '"'
    return text


def first_scalar(value: Any) -> float | None:
    if isinstance(value, list) and value:
        try:
            return float(value[0])
        except Exception:
            return None
    try:
        return float(value)
    except Exception:
        return None


def none_safe_sub(value: Any, target: float) -> float | None:
    try:
        out = float(value) - float(target)
    except Exception:
        return None
    return out if np.isfinite(out) else None
