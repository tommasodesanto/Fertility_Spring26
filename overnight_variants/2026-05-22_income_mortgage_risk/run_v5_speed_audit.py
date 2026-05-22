#!/usr/bin/env python3
"""Speed audit for the V5 HANK-z benchmark-normalized outside closure."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z
from run_income_mortgage_risk_v2_scenarios import BENCHMARK_P
from run_income_mortgage_risk_v5_hank_z_outside_closure import prepare_parameters


def parse_cases(raw: str, default_iter: int) -> list[tuple[int, int, int]]:
    cases: list[tuple[int, int, int]] = []
    for part in raw.split(","):
        token = part.strip().lower()
        if not token:
            continue
        bits = token.split("x")
        if len(bits) not in (2, 3):
            raise ValueError(f"case must have form NBxNZ or NBxNZxITER, got {part!r}")
        nb_raw, nz_raw = bits[:2]
        iter_raw = bits[2] if len(bits) == 3 else default_iter
        cases.append((int(nb_raw), int(nz_raw), int(iter_raw)))
    if not cases:
        raise ValueError("at least one audit case is required")
    return cases


def configure_benchmark_normalized_v5(P: SimpleNamespace, kappa_entry: float, target_q: float) -> SimpleNamespace:
    P.population_closure = "outside_option_benchmark_normalized"
    P.kappa_entry = float(kappa_entry)
    P.local_birth_entry_weight = 1.0
    P.target_city_entry_prob = float(target_q)
    P.calibrate_outside_value_to_entry_prob = True
    P.outside_benchmark_target_total_population = float(getattr(P, "N_target", 1.0))
    P.allow_uncalibrated_outside_value = False
    P.collect_ge_trace = True
    return P


def run_case(
    *,
    nb: int,
    nz: int,
    max_iter_eq: int,
    rho_z: float,
    sigma_z: float,
    kappa_entry: float,
    target_q: float,
    tol_eq: float,
    quiet: bool,
) -> dict[str, object]:
    P, _, _ = prepare_parameters(nb, max_iter_eq, tol_eq)
    P = configure_benchmark_normalized_v5(P, kappa_entry, target_q)
    b_grid = make_grid(P)
    t0 = time.perf_counter()
    sol, P_out, p_eq = solve_equilibrium_hank_z(
        BENCHMARK_P,
        P,
        b_grid,
        nz=nz,
        rho_z=rho_z,
        sigma_z=sigma_z,
        verbose=not quiet,
    )
    elapsed = time.perf_counter() - t0
    timings = dict(getattr(sol, "timings", {}))
    scale = getattr(sol, "accounting_scale", SimpleNamespace())
    bellman = float(timings.get("bellman_hank_z", np.nan))
    dist = float(timings.get("distribution_hank_z", np.nan))
    n_full = max(int(timings.get("n_full", 0)), 1)
    n_dist = max(int(timings.get("n_dist", 0)), 1)
    residual = elapsed - bellman - dist
    return {
        "nb": int(nb),
        "nz": int(nz),
        "max_iter_eq": int(max_iter_eq),
        "iterations_completed": int(timings.get("iterations_completed", 0)),
        "accepted": bool(timings.get("accepted", False)),
        "convergence_reason": timings.get("convergence_reason", ""),
        "best_eq_error": float(timings.get("best_eq_error", np.nan)),
        "final_eq_error": float(timings.get("final_eq_error", np.nan)),
        "elapsed_sec": float(elapsed),
        "bellman_sec": bellman,
        "distribution_sec": dist,
        "residual_final_stats_overhead_sec": float(residual),
        "bellman_sec_per_iter": bellman / n_full,
        "distribution_sec_per_iter": dist / n_dist,
        "elapsed_sec_per_iter": float(elapsed) / max(int(timings.get("iterations_completed", 0)), 1),
        "p_P": float(p_eq[0]),
        "p_C": float(p_eq[1]),
        "scale_factor": float(getattr(scale, "scale_factor", np.nan)),
        "city_entry_prob_total": float(getattr(scale, "city_entry_prob_total", np.nan)),
        "outside_entry_prob": float(getattr(scale, "outside_entry_prob", np.nan)),
        "outside_entry_flow": float(getattr(scale, "outside_entry_flow", np.nan)),
        "tfr_trace_final": float(getattr(sol, "mean_parity", np.nan) * 2.0),
        "own_rate": float(getattr(sol, "own_rate", np.nan)),
        "rho_z": float(rho_z),
        "sigma_z": float(sigma_z),
        "kappa_entry": float(kappa_entry),
        "target_q": float(target_q),
        "tol_eq": float(tol_eq),
        "timings": timings,
        "ge_trace": getattr(sol, "ge_trace", []),
        "P_population_closure": str(getattr(P_out, "population_closure", "")),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "nb",
        "nz",
        "max_iter_eq",
        "iterations_completed",
        "accepted",
        "convergence_reason",
        "best_eq_error",
        "final_eq_error",
        "elapsed_sec",
        "bellman_sec",
        "distribution_sec",
        "residual_final_stats_overhead_sec",
        "bellman_sec_per_iter",
        "distribution_sec_per_iter",
        "elapsed_sec_per_iter",
        "p_P",
        "p_C",
        "scale_factor",
        "city_entry_prob_total",
        "outside_entry_prob",
        "outside_entry_flow",
        "tfr_trace_final",
        "own_rate",
        "rho_z",
        "sigma_z",
        "kappa_entry",
        "target_q",
        "tol_eq",
        "P_population_closure",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def write_report(path: Path, rows: list[dict[str, object]], full_v5: dict[str, object] | None) -> None:
    table = "\n".join(
        "| {nb} | {nz} | {it} | {elapsed:.2f} | {bell:.2f} | {dist:.2f} | {resid:.2f} | {bell_it:.2f} | {dist_it:.2f} |".format(
            nb=int(row["nb"]),
            nz=int(row["nz"]),
            it=int(row["iterations_completed"]),
            elapsed=float(row["elapsed_sec"]),
            bell=float(row["bellman_sec"]),
            dist=float(row["distribution_sec"]),
            resid=float(row["residual_final_stats_overhead_sec"]),
            bell_it=float(row["bellman_sec_per_iter"]),
            dist_it=float(row["distribution_sec_per_iter"]),
        )
        for row in rows
    )
    full_block = "Accepted full V5 run was not found in the JSON log."
    if full_v5 is not None:
        timings = full_v5.get("timings", {})
        full_block = "\n".join(
            [
                f"- elapsed seconds: `{float(full_v5.get('elapsed_sec', np.nan)):.2f}`",
                f"- iterations completed: `{int(timings.get('iterations_completed', 0))}`",
                f"- Bellman seconds: `{float(timings.get('bellman_hank_z', np.nan)):.2f}`",
                f"- forward-distribution seconds: `{float(timings.get('distribution_hank_z', np.nan)):.2f}`",
                f"- residual final-stat / overhead seconds: `{float(full_v5.get('elapsed_sec', np.nan)) - float(timings.get('bellman_hank_z', np.nan)) - float(timings.get('distribution_hank_z', np.nan)):.2f}`",
                f"- Bellman seconds per GE iteration: `{float(timings.get('bellman_hank_z', np.nan)) / max(int(timings.get('n_full', 1)), 1):.2f}`",
                f"- forward seconds per GE iteration: `{float(timings.get('distribution_hank_z', np.nan)) / max(int(timings.get('n_dist', 1)), 1):.2f}`",
            ]
        )

    text = f"""# V5 HANK-z Speed Audit

This audit runs the benchmark-normalized outside-option V5 closure for a fixed
small number of GE iterations. These are timing probes, not accepted
calibrations. The goal is to separate Bellman time, fast forward-distribution
time, and residual final-statistics / overhead time as \(N_z\) and \(N_b\)
change.

## Fixed-Iteration Audit

| \(N_b\) | \(N_z\) | GE iters | elapsed sec | Bellman sec | forward sec | residual sec | Bellman / iter | forward / iter |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

## Accepted Full V5 Reference

{full_block}

## Read

The runtime problem is not the benchmark normalization arithmetic itself. The
expensive objects are the structural HANK-\(z\) Bellman and forward passes,
with an additional final-statistics overhead after the GE loop. The accepted
\(N_b=30,N_z=7\) V5 solve is therefore not calibration-ready without a faster
mode or a smaller audit/calibration grid.
"""
    path.write_text(text)


def load_full_v5_reference(root: Path) -> dict[str, object] | None:
    path = root / "income_mortgage_risk_v5_hank_z_outside_closure.log"
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    if not data.get("ok", False):
        return None
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="30x3x4,30x5x4,30x7x4,20x7x4,30x7x8")
    parser.add_argument("--max-iter-eq", type=int, default=4)
    parser.add_argument("--rho-z", type=float, default=0.95)
    parser.add_argument("--sigma-z", type=float, default=0.35)
    parser.add_argument("--kappa-entry", type=float, default=1_000_000.0)
    parser.add_argument("--target-city-entry-prob", type=float, default=0.9)
    parser.add_argument("--tol-eq", type=float, default=0.0)
    parser.add_argument("--output-stem", default="speed_audit_v5_benchmark_normalized")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    rows = []
    for nb, nz, max_iter_eq in parse_cases(args.cases, args.max_iter_eq):
        rows.append(
            run_case(
                nb=nb,
                nz=nz,
                max_iter_eq=max_iter_eq,
                rho_z=args.rho_z,
                sigma_z=args.sigma_z,
                kappa_entry=args.kappa_entry,
                target_q=args.target_city_entry_prob,
                tol_eq=args.tol_eq,
                quiet=args.quiet,
            )
        )

    output_stem = str(args.output_stem)
    write_csv(root / f"{output_stem}.csv", rows)
    (root / f"{output_stem}.json").write_text(json.dumps(rows, indent=2, default=str))
    report_name = "REPORT_V5_SPEED_AUDIT.md" if output_stem == "speed_audit_v5_benchmark_normalized" else f"{output_stem}.md"
    write_report(root / report_name, rows, load_full_v5_reference(root))


if __name__ == "__main__":
    main()
