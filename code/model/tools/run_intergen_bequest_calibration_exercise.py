#!/usr/bin/env python3
"""Small diagnostic calibration around the parent-gated luxury bequest proposal.

The exercise has two layers.  First, it profiles fixed positive values of
``theta0`` and ``theta1`` at the verified clean-frontier coordinates, with and
without an externally pinned post-retirement owner-LTV taper.  Second, when
``--polish-evals`` is positive, it performs a bounded coordinate polish of the
remaining eleven parameters inside the best cell of each LTV arm.

This is a diagnostic frontier, not a promoted SMM specification: the proposed
late-life decumulation statistic is reported but has no empirical target yet.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import numpy as np

from intergen_housing_fertility.calibration import (
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
)
from intergen_housing_fertility.local_panel import income_process_overrides, jsonable
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)
from intergen_housing_fertility.solver import InfeasibleThetaError, run_model_cp_dt


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SEED = (
    ROOT
    / "output/model/full_calibration_audit_20260713_claude/round3/remote_mirror/sweep"
    / "clean_k0_frontier__tight_rep1.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_bequest_calibration_exercise"
MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0
PERIOD_YEARS = 4.0

# Eleven coordinates left free when theta0, theta1, theta_n, and tenure noise
# are fixed within a profile cell.  The h_bar_0 domain includes the verified
# relaxed frontier; this is intentionally wider than the production box.
POLISH_BOUNDS: tuple[tuple[str, float, float], ...] = (
    ("beta_annual", 0.940, 0.995),
    ("alpha_cons", 0.400, 0.950),
    ("c_bar_0", 0.080, 1.280),
    ("c_bar_n", 0.050, 1.500),
    ("h_bar_0", 0.250, 6.000),
    ("h_bar_jump", 0.050, 2.500),
    ("h_bar_n", 0.020, 2.000),
    ("psi_child", 0.000, 0.350),
    ("kappa_fert", 1.000, 12.000),
    ("chi", 0.400, 1.500),
    ("H0", 1.000, 20.000),
)
THETA1_PROFILE_BOUND = (0.10, 2.00)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-record", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--theta0-grid", default="0.05,0.30")
    parser.add_argument("--theta1-grid", default="0.25,1.00")
    parser.add_argument("--schedule-arms", default="none,linear_66_82")
    parser.add_argument("--polish-evals", type=int, default=0)
    parser.add_argument("--initial-step", type=float, default=0.04)
    parser.add_argument("--minutes", type=float, default=30.0)
    parser.add_argument("--J", type=int, default=PRODUCTION_J)
    parser.add_argument("--Nb", type=int, default=PRODUCTION_SEARCH_NB)
    parser.add_argument("--max-iter-eq", type=int, default=40)
    parser.add_argument("--tol-eq", type=float, default=2.5e-5)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def parse_grid(raw: str, *, positive: bool = False) -> list[float]:
    values = [float(item.strip()) for item in str(raw).split(",") if item.strip()]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError(f"invalid grid: {raw!r}")
    if positive and any(value <= 0.0 for value in values):
        raise ValueError("grid values must be strictly positive")
    return sorted(set(values))


def load_theta(path: Path) -> dict[str, float]:
    payload: Any = json.loads(path.read_text())
    while isinstance(payload, dict) and "theta" not in payload and "record" in payload:
        payload = payload["record"]
    raw = payload.get("theta", payload) if isinstance(payload, dict) else payload
    if not isinstance(raw, dict):
        raise ValueError(f"{path} does not contain a theta object")
    theta = {str(key): float(value) for key, value in raw.items() if isinstance(value, (int, float))}
    required = {"beta" if name == "beta_annual" else name for name, _, _ in POLISH_BOUNDS}
    required.update({"theta0", "theta_n", "tenure_choice_kappa"})
    missing = required - set(theta)
    if missing:
        raise ValueError(f"seed record is missing coordinates: {sorted(missing)}")
    return theta


def target_system() -> tuple[dict[str, float], dict[str, float]]:
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = ROOMS_WEIGHT
    if set(targets) != set(weights) or len(targets) != 15:
        raise ValueError("exercise requires the active 15-moment target system")
    return targets, weights


def schedule_overrides(arm: str) -> dict[str, Any]:
    if arm == "none":
        return {"owner_ltv_taper": False}
    if arm == "linear_66_82":
        return {
            "owner_ltv_taper": True,
            "owner_ltv_taper_start_age": 66.0,
            "owner_ltv_taper_end_age": 82.0,
            "owner_ltv_terminal_share": 0.0,
        }
    raise ValueError(f"unknown schedule arm: {arm}")


def common_overrides(args: argparse.Namespace) -> dict[str, Any]:
    income = income_process_overrides(
        5,
        "rouwenhorst",
        MATCHED_ANNUAL_INNOVATION_SD,
        MATCHED_ANNUAL_RHO,
    )
    return {
        **base_overrides(J=args.J, Nb=args.Nb, n_house=5, max_iter_eq=args.max_iter_eq),
        **production_profile_overrides(),
        **income,
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
        "eta_supply": np.array([1.75]),
        "lambda_d": 0.0,
        "debt_taper_start_age": 42.0,
        "debt_taper_end_age": 62.0,
        "normalize_bequest_utility": True,
    }


def target_fit_rows(
    label: str,
    moments: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights[name])
        rows.append(
            {
                "label": label,
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )
    return rows


def late_life_stats(sol: Any, P: Any) -> dict[str, Any]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float)
    mass_age = g.sum(axis=(0, 1, 2, 4, 5, 6))
    liquid_age = np.einsum("b,btijznc->j", b_grid, g) / np.maximum(mass_age, 1e-300)
    own_age = g[:, 1:].sum(axis=(0, 1, 2, 4, 5, 6)) / np.maximum(mass_age, 1e-300)
    ages = np.asarray([float(P.age_start + j * P.da) for j in range(g.shape[3])])

    def window(values: np.ndarray, lo: float, hi: float) -> float:
        selected = (ages >= lo) & (ages <= hi)
        weights = mass_age[selected]
        return float(np.sum(values[selected] * weights) / max(float(np.sum(weights)), 1e-300))

    early = window(liquid_age, 62.0, 74.0)
    late = window(liquid_age, 74.0, 82.0)
    return {
        "ages": ages.tolist(),
        "mean_liquid_wealth_by_age": liquid_age.tolist(),
        "ownership_by_age": own_age.tolist(),
        "wealth_62_74": early,
        "wealth_74_82": late,
        "decum_ratio_wealth_74plus_over_62_74": late / max(early, 1e-12),
        "own_62_74": window(own_age, 62.0, 74.0),
        "own_74plus": window(own_age, 74.0, 82.0),
    }


def strict_convergence(sol: Any, P: Any) -> tuple[bool, float]:
    residual = float(getattr(sol, "best_max_abs_rel_excess", math.nan))
    timings = dict(getattr(sol, "timings", {}))
    strict = bool(
        timings.get("strict_converged", getattr(sol, "converged", False))
        and math.isfinite(residual)
        and residual <= float(P.tol_eq)
    )
    return strict, residual


def evaluate(
    *,
    label: str,
    stage: str,
    arm: str,
    theta: dict[str, float],
    proposed: bool,
    args: argparse.Namespace,
    targets: dict[str, float],
    weights: dict[str, float],
) -> dict[str, Any]:
    started = time.perf_counter()
    overrides = {
        **common_overrides(args),
        **schedule_overrides(arm),
        **theta,
        "theta_n": 0.0,
        "tenure_choice_kappa": 0.0,
        "bequest_spec": "parent_gated_luxury" if proposed else "linear_child_scale",
    }
    try:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=not args.quiet)
        moments = extract_moments(sol, P)
        strict, residual = strict_convergence(sol, P)
        status = "ok" if strict else "non_strict"
        profiles = late_life_stats(sol, P)
        fit = target_fit_rows(label, moments, targets, weights)
        loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
        timings = dict(getattr(sol, "timings", {}))
        price = float(np.asarray(p_eq).reshape(-1)[0])
        feasibility_census: list[dict[str, Any]] = []
        error = ""
    except InfeasibleThetaError as exc:
        moments, profiles, fit, timings = {}, {}, [], {}
        status, strict, residual, loss, price = "infeasible_theta", False, math.inf, math.inf, math.nan
        feasibility_census, error = list(exc.census), str(exc)
    except Exception as exc:  # noqa: BLE001 - every failed cell must checkpoint.
        moments, profiles, fit, timings = {}, {}, [], {}
        status, strict, residual, loss, price = f"failed:{type(exc).__name__}", False, math.inf, math.inf, math.nan
        feasibility_census, error = [], str(exc)
    return {
        "label": label,
        "stage": stage,
        "schedule_arm": arm,
        "bequest_spec": "parent_gated_luxury" if proposed else "linear_child_scale",
        "status": status,
        "strict_converged": strict,
        "rank_loss": loss,
        "market_residual": residual,
        "price": price,
        "theta": theta,
        "moments": moments,
        "target_fit": fit,
        "late_life": profiles,
        "timings": timings,
        "feasibility_census": feasibility_census,
        "error": error,
        "elapsed_sec": float(time.perf_counter() - started),
    }


def admissible(record: dict[str, Any]) -> bool:
    return bool(record.get("strict_converged", False)) and math.isfinite(float(record.get("rank_loss", math.inf)))


def best_record(records: list[dict[str, Any]], arm: str, *, positive: bool = True) -> dict[str, Any] | None:
    eligible = [
        record
        for record in records
        if record.get("schedule_arm") == arm
        and admissible(record)
        and (not positive or float(record.get("theta", {}).get("theta0", 0.0)) > 0.0)
    ]
    return min(eligible, key=lambda record: float(record["rank_loss"])) if eligible else None


def unit_from_theta(theta: dict[str, float]) -> np.ndarray:
    out = []
    for name, lo, hi in POLISH_BOUNDS:
        key = "beta" if name == "beta_annual" else name
        value = float(theta[key]) ** (1.0 / PERIOD_YEARS) if name == "beta_annual" else float(theta[key])
        out.append((value - lo) / (hi - lo))
    return np.clip(np.asarray(out, dtype=float), 0.0, 1.0)


def theta_from_unit(unit: np.ndarray, fixed: dict[str, float]) -> dict[str, float]:
    theta = dict(fixed)
    for value, (name, lo, hi) in zip(np.asarray(unit, dtype=float), POLISH_BOUNDS):
        level = float(lo + np.clip(value, 0.0, 1.0) * (hi - lo))
        theta["beta" if name == "beta_annual" else name] = level**PERIOD_YEARS if name == "beta_annual" else level
    return theta


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def parameter_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    bounds = {name: (lo, hi) for name, lo, hi in POLISH_BOUNDS}
    bounds.update({"theta0": (0.0, 2.0), "theta1": THETA1_PROFILE_BOUND, "theta_n": (0.0, 1.5), "tenure_choice_kappa": (0.0, 0.12)})
    rows: list[dict[str, Any]] = []
    for record in records:
        for parameter, (lo, hi) in bounds.items():
            key = "beta" if parameter == "beta_annual" else parameter
            estimate = float(record["theta"][key])
            if parameter == "beta_annual":
                estimate = estimate ** (1.0 / PERIOD_YEARS)
            span = hi - lo
            rows.append(
                {
                    "label": record["label"],
                    "parameter": parameter,
                    "estimate": estimate,
                    "lower_or_external_restriction": lo,
                    "upper_or_external_restriction": hi,
                    "role": "profile_fixed" if parameter in {"theta0", "theta1", "theta_n", "tenure_choice_kappa"} else ("locally_searched" if record["stage"] == "polish" else "held_at_clean_frontier"),
                    "near_bound": min(estimate - lo, hi - estimate) <= 0.02 * span,
                }
            )
    return rows


def write_progress(outdir: Path, records: list[dict[str, Any]]) -> None:
    latest = records[-1]
    (outdir / "latest_completed_case.json").write_text(json.dumps(jsonable(latest), indent=2, sort_keys=True))
    best = {
        arm: best_record(records, arm, positive=False)
        for arm in sorted({str(record["schedule_arm"]) for record in records})
    }
    (outdir / "best_so_far.json").write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
    with (outdir / "cases.jsonl").open("a") as handle:
        handle.write(json.dumps(jsonable(latest), sort_keys=True) + "\n")


def plot_packet(outdir: Path, records: list[dict[str, Any]], arms: list[str]) -> None:
    selected = [record for record in records if record["stage"] == "baseline"]
    selected.extend(filter(None, (best_record(records, arm) for arm in arms)))
    selected_unique = {record["label"]: record for record in selected if record is not None}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for record in selected_unique.values():
        late = record.get("late_life", {})
        if not late:
            continue
        axes[0].plot(late["ages"], late["ownership_by_age"], marker="o", label=record["label"])
        axes[1].plot(late["ages"], late["mean_liquid_wealth_by_age"], marker="o", label=record["label"])
    axes[0].axhline(0.76426097, color="black", linestyle="--", linewidth=1, label="old-age target")
    axes[0].set(ylabel="ownership rate", xlabel="age", ylim=(0, 1.02))
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set(ylabel="mean liquid wealth", xlabel="age")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / "age_profiles.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, max(len(arms), 1), figsize=(5.2 * max(len(arms), 1), 4.2), squeeze=False)
    for axis, arm in zip(axes[0], arms):
        grid = [record for record in records if record["stage"] == "profile" and record["schedule_arm"] == arm and admissible(record)]
        if not grid:
            axis.set_axis_off()
            continue
        scatter = axis.scatter(
            [record["theta"]["theta1"] for record in grid],
            [record["theta"]["theta0"] for record in grid],
            c=[record["rank_loss"] for record in grid],
            s=180,
            cmap="viridis_r",
        )
        for record in grid:
            axis.annotate(f"{record['rank_loss']:.2f}", (record["theta"]["theta1"], record["theta"]["theta0"]), ha="center", va="center", fontsize=8)
        axis.set(title=arm, xlabel=r"luxury shifter $\theta_1$", ylabel=r"bequest strength $\theta_0$")
        fig.colorbar(scatter, ax=axis, label="15-moment loss")
    fig.tight_layout()
    fig.savefig(outdir / "profile_surface.png", dpi=180)
    plt.close(fig)

    fit_records = [record for record in selected_unique.values() if record.get("target_fit")]
    if fit_records:
        moments = [row["moment"] for row in fit_records[0]["target_fit"]]
        x = np.arange(len(moments))
        width = 0.8 / len(fit_records)
        fig, axis = plt.subplots(figsize=(13, 5.2))
        for idx, record in enumerate(fit_records):
            contrib = {row["moment"]: row["loss_contribution"] for row in record["target_fit"]}
            axis.bar(x + (idx - (len(fit_records) - 1) / 2) * width, [contrib[name] for name in moments], width=width, label=record["label"])
        axis.set_xticks(x, moments, rotation=65, ha="right", fontsize=8)
        axis.set_ylabel("loss contribution")
        axis.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / "loss_contributions.png", dpi=180)
        plt.close(fig)


def write_readme(
    outdir: Path,
    args: argparse.Namespace,
    records: list[dict[str, Any]],
    arms: list[str],
    elapsed: float,
) -> None:
    lines = [
        "# Parent-gated luxury bequest calibration exercise",
        "",
        f"Completed {len(records)} solves in {elapsed / 60.0:.1f} minutes.",
        "",
        "This is a diagnostic profile, not a promoted SMM calibration. The active objective has 15 moments. Within each positive-bequest profile cell, `theta0`, `theta1`, `theta_n=0`, and deterministic tenure are externally fixed; the optional polish searches the remaining 11 coordinates. The late-life wealth ratio is diagnostic-only because no empirical target or weight has been approved.",
        "",
        "The proposed terminal utility is normalized as $1\\{n\\ge1\\}\\theta_0[u(\\theta_1+b)-u(\\theta_1)]$. Omitting the subtraction would add a parent-only CRRA level constant and contaminate fertility comparisons.",
        "",
        "## Best strict record by schedule arm",
        "",
        "| arm | loss | theta0 | theta1 | old ownership | own 74+ | decumulation ratio | residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in arms:
        record = best_record(records, arm)
        if record is None:
            lines.append(f"| {arm} | no strict positive-bequest record | | | | | | |")
            continue
        moments, late = record["moments"], record["late_life"]
        lines.append(
            f"| {arm} | {record['rank_loss']:.6f} | {record['theta']['theta0']:.4g} | {record['theta']['theta1']:.4g} | {moments.get('old_age_own_rate', math.nan):.4f} | {late.get('own_74plus', math.nan):.4f} | {late.get('decum_ratio_wealth_74plus_over_62_74', math.nan):.4f} | {record['market_residual']:.2e} |"
        )
    lines.extend(
        [
            "",
            "Complete target fits are in `target_fit_all.csv`; all parameter values, restrictions, search roles, and bound flags are in `parameters_all.csv`. `latest_completed_case.json` and `best_so_far.json` are recoverable run checkpoints.",
            "",
            "Regenerate the packet with:",
            "",
            "```bash",
            "cd code/model",
            f"PYTHONPATH=$PWD .venv/bin/python tools/run_intergen_bequest_calibration_exercise.py --seed-record {args.seed_record} --outdir {outdir} --theta0-grid {args.theta0_grid} --theta1-grid {args.theta1_grid} --schedule-arms {args.schedule_arms} --polish-evals {args.polish_evals}",
            "```",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    if args.smoke:
        # Preserve the full age loop so the 66--82 LTV schedule is exercised.
        args.J = PRODUCTION_J
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25
        args.theta0_grid = "0.05"
        args.theta1_grid = "0.25"
        args.polish_evals = 0
        args.minutes = min(float(args.minutes), 5.0)
    theta0_grid = parse_grid(args.theta0_grid, positive=True)
    theta1_grid = parse_grid(args.theta1_grid, positive=True)
    if any(not THETA1_PROFILE_BOUND[0] <= value <= THETA1_PROFILE_BOUND[1] for value in theta1_grid):
        raise ValueError(f"theta1 grid must lie in {THETA1_PROFILE_BOUND}")
    arms = [item.strip() for item in str(args.schedule_arms).split(",") if item.strip()]
    for arm in arms:
        schedule_overrides(arm)
    if not arms:
        raise ValueError("at least one schedule arm is required")

    seed = load_theta(args.seed_record)
    seed.update(theta0=0.0, theta_n=0.0, tenure_choice_kappa=0.0)
    targets, weights = target_system()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "cases.jsonl").write_text("")
    metadata = {
        "status": "smoke" if args.smoke else "diagnostic_profile_not_promoted_smm",
        "seed_record": str(args.seed_record.resolve()),
        "seed_theta": seed,
        "target_set": PRODUCTION_TARGET_SET,
        "target_count": len(targets),
        "targets": targets,
        "weights": weights,
        "theta0_grid": theta0_grid,
        "theta1_grid": theta1_grid,
        "schedule_arms": arms,
        "polish_bounds": POLISH_BOUNDS,
        "polish_free_parameter_count": len(POLISH_BOUNDS),
        "profile_fixed": ["theta0", "theta1", "theta_n=0", "tenure_choice_kappa=0"],
        "decumulation_moment_status": "diagnostic_only_no_empirical_target",
        "J": args.J,
        "Nb": args.Nb,
        "max_iter_eq": args.max_iter_eq,
        "tol_eq": args.tol_eq,
        "time_budget_minutes": args.minutes,
        "polish_eval_budget_per_arm": args.polish_evals,
    }
    (args.outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True))

    records: list[dict[str, Any]] = []
    started = time.perf_counter()
    deadline = started + max(1.0, float(args.minutes) * 60.0)
    max_cases = math.inf if int(args.max_cases) <= 0 else int(args.max_cases)

    def run_case(**kwargs: Any) -> bool:
        if len(records) >= max_cases or time.perf_counter() >= deadline:
            return False
        record = evaluate(args=args, targets=targets, weights=weights, **kwargs)
        records.append(record)
        write_progress(args.outdir, records)
        print(
            f"case {len(records)} {record['label']}: status={record['status']} "
            f"loss={record['rank_loss']:.6g} residual={record['market_residual']:.2e} "
            f"elapsed={(time.perf_counter() - started) / 60.0:.1f}m",
            flush=True,
        )
        return True

    baseline_theta = dict(seed)
    baseline_theta["theta1"] = 0.25
    run_case(
        label="baseline_no_bequest",
        stage="baseline",
        arm="none",
        theta=baseline_theta,
        proposed=False,
    )

    for arm in arms:
        for theta0 in theta0_grid:
            for theta1 in theta1_grid:
                theta = dict(seed)
                theta.update(theta0=theta0, theta1=theta1)
                if not run_case(
                    label=f"profile_{arm}_th0_{theta0:g}_th1_{theta1:g}",
                    stage="profile",
                    arm=arm,
                    theta=theta,
                    proposed=True,
                ):
                    break

    for arm in arms:
        start_record = best_record(records, arm)
        if start_record is None or int(args.polish_evals) <= 0:
            continue
        current_theta = dict(start_record["theta"])
        fixed = {
            "theta0": current_theta["theta0"],
            "theta1": current_theta["theta1"],
            "theta_n": 0.0,
            "tenure_choice_kappa": 0.0,
        }
        current_unit = unit_from_theta(current_theta)
        current_best = start_record
        step = float(np.clip(args.initial_step, 0.002, 0.20))
        eval_count = 0
        while eval_count < int(args.polish_evals) and time.perf_counter() < deadline:
            improved_in_pass = False
            for dim in range(len(POLISH_BOUNDS)):
                for direction in (1.0, -1.0):
                    if eval_count >= int(args.polish_evals):
                        break
                    candidate_unit = current_unit.copy()
                    candidate_unit[dim] = np.clip(candidate_unit[dim] + direction * step, 0.0, 1.0)
                    theta = theta_from_unit(candidate_unit, fixed)
                    completed = run_case(
                        label=f"polish_{arm}_{eval_count + 1:03d}",
                        stage="polish",
                        arm=arm,
                        theta=theta,
                        proposed=True,
                    )
                    if not completed:
                        break
                    candidate = records[-1]
                    eval_count += 1
                    if admissible(candidate) and float(candidate["rank_loss"]) < float(current_best["rank_loss"]):
                        current_best = candidate
                        current_unit = candidate_unit
                        improved_in_pass = True
                    if time.perf_counter() >= deadline:
                        break
                if eval_count >= int(args.polish_evals) or time.perf_counter() >= deadline:
                    break
            if not improved_in_pass:
                step *= 0.5
            if step < 0.002:
                break

    all_fit = [row for record in records for row in record.get("target_fit", [])]
    write_csv(args.outdir / "target_fit_all.csv", all_fit)
    write_csv(args.outdir / "parameters_all.csv", parameter_rows(records))
    plot_packet(args.outdir, records, arms)
    elapsed = time.perf_counter() - started
    write_readme(args.outdir, args, records, arms, elapsed)
    summary = {
        "status": "complete" if time.perf_counter() < deadline else "time_budget_reached",
        "case_count": len(records),
        "elapsed_sec": elapsed,
        "best_by_arm": {arm: best_record(records, arm, positive=False) for arm in arms},
        "baseline": next((record for record in records if record["stage"] == "baseline"), None),
    }
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
