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
from .solver import run_model_cp_dt


DEFAULT_INCOME_GRID_5 = np.array([0.60, 0.80, 1.00, 1.20, 1.40])
DEFAULT_INCOME_WEIGHTS_5 = np.array([0.10, 0.20, 0.40, 0.20, 0.10])
DEFAULT_RHO_Z = 0.85


def run_local_panel(
    outdir: Path,
    *,
    n_cases: int = 144,
    seed: int = 20260608,
    J: int = 16,
    Nb: int = 60,
    n_house: int = 6,
    max_iter_eq: int = 25,
    workers: int = 6,
    minutes: float = 30.0,
    income_states: int = 5,
    diagnostic_best: int = 3,
    target_set: str = "candidate_no_timing_v0",
    include_anchors: bool = True,
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
    candidates = local_panel_candidates(n_cases, seed, include_anchors=include_anchors)
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
        "full_old_targets": OLD_NONLOCATION_TARGETS,
        "full_old_weights": OLD_NONLOCATION_WEIGHTS,
        "varied_internal_parameters": [
            "beta",
            "alpha_cons",
            "b_entry_fixed",
            "c_bar_0",
            "c_bar_n",
            "h_bar_0",
            "h_bar_jump",
            "h_bar_n",
            "psi_child",
            "kappa_fert",
            "chi",
            "theta0",
            "theta_n",
        ],
        "fixed_external_or_first_stage": [
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
            "tenure_choice_kappa",
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

    records_sorted = sorted(records, key=lambda r: float(r.get("rank_loss", math.inf)))
    summary = {
        "best": records_sorted[0] if records_sorted else None,
        "top_records": records_sorted[: min(10, len(records_sorted))],
        "elapsed_sec": float(time.perf_counter() - start),
        "n_cases_completed": int(len(records)),
        "n_cases_submitted": int(next_idx),
        "stopped_by_time_budget": bool(next_idx < len(candidates)),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    write_panel_summary_tables(outdir, records_sorted, rank_targets)
    write_panel_plots(outdir, records_sorted, rank_targets)
    if diagnostic_best > 0 and records_sorted:
        write_best_case_diagnostics(
            outdir,
            records_sorted[: int(diagnostic_best)],
            J=J,
            Nb=Nb,
            n_house=n_house,
            max_iter_eq=max_iter_eq,
            income=income,
        )
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
) -> dict[str, Any]:
    t0 = time.perf_counter()
    theta = dict(candidate["theta"])
    overrides = {
        **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
        **income,
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
    except Exception as exc:  # noqa: BLE001 - panel should checkpoint failed parameter vectors.
        moments = {}
        rank_loss = math.inf
        full_loss = math.inf
        status = f"failed: {type(exc).__name__}: {exc}"
        p_eq = np.array([np.nan])
        market_residual = math.inf
        timings = {}
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
        "elapsed_sec": float(time.perf_counter() - t0),
        "timings": jsonable(timings),
    }


def local_panel_candidates(n_cases: int, seed: int, *, include_anchors: bool = True) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
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


def keep_internal_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    allowed = {
        "beta",
        "alpha_cons",
        "b_entry_fixed",
        "c_bar_0",
        "c_bar_n",
        "h_bar_0",
        "h_bar_jump",
        "h_bar_n",
        "psi_child",
        "kappa_fert",
        "chi",
        "theta0",
        "theta_n",
    }
    theta = {k: v for k, v in dict(candidate["theta"]).items() if k in allowed}
    return {"label": str(candidate["label"]), "theta": theta}


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
        "b_entry_fixed": rng.uniform(-2.0, 6.0),
        "c_bar_0": rng.uniform(0.04, 0.16) * PERIOD_YEARS,
        "c_bar_n": rng.uniform(child_cost_lo, child_cost_hi),
        "h_bar_0": rng.uniform(3.0, 4.8),
        "h_bar_jump": rng.uniform(0.40, 1.35),
        "h_bar_n": rng.uniform(0.30, 1.10),
        "psi_child": rng.uniform(psi_child_lo, psi_child_hi),
        "kappa_fert": rng.uniform(3.5, 7.5),
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
    if best is None:
        return True
    return float(record.get("rank_loss", math.inf)) < float(best.get("rank_loss", math.inf))


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
    ok = [r for r in records if np.isfinite(float(r.get("rank_loss", math.inf)))]
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
        ("own_rate", "young_liquid_wealth_to_income", "Ownership vs young liquid wealth/income"),
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
