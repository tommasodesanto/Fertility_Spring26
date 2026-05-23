#!/usr/bin/env python3
"""Run and plot the V5 benchmark-normalized outside-closure prototype at Nz=5."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z

from plot_income_mortgage_risk_v4_hank_z_ge import (
    plot_aggregate_sorting,
    plot_center_ownership_heatmaps,
    plot_fertility_wealth,
    plot_moment_gaps,
    plot_moments,
    plot_sorting_age_profiles,
    plot_tenure_stacks,
    plot_trace,
    plot_z_process,
    setup_plot_style,
    sorting_tables,
)
from run_income_mortgage_risk_v2_scenarios import (
    BENCHMARK_P,
    namespace_float_dict,
)
from run_income_mortgage_risk_v3_hank_z import write_diagnostics as write_hank_z_diagnostics
from run_income_mortgage_risk_v5_hank_z_outside_closure import (
    prepare_parameters,
    write_closure_diagnostics,
    write_report,
    write_results,
    write_trace,
)


OUTPUT_STEM = "income_mortgage_risk_v5_hank_z_outside_closure_nz5"
FIGURE_DIR = "figures_v5_hank_z_outside_closure_nz5"


def configure_v5_benchmark_normalized(
    P: SimpleNamespace,
    *,
    kappa_entry: float,
    target_city_entry_prob: float,
    local_birth_entry_weight: float,
) -> SimpleNamespace:
    P.population_closure = "outside_option_benchmark_normalized"
    P.kappa_entry = float(kappa_entry)
    P.local_birth_entry_weight = float(local_birth_entry_weight)
    P.target_city_entry_prob = float(target_city_entry_prob)
    P.calibrate_outside_value_to_entry_prob = True
    P.outside_benchmark_target_total_population = float(getattr(P, "N_target", 1.0))
    P.allow_uncalibrated_outside_value = False
    P.collect_ge_trace = True
    return P


def solve_v5(args: argparse.Namespace):
    P, targets, weights = prepare_parameters(args.nb, args.max_iter_eq, args.tol_eq)
    P = configure_v5_benchmark_normalized(
        P,
        kappa_entry=args.kappa_entry,
        target_city_entry_prob=args.target_city_entry_prob,
        local_birth_entry_weight=args.local_birth_entry_weight,
    )
    b_grid = make_grid(P)
    t0 = time.perf_counter()
    sol, P_out, p_eq = solve_equilibrium_hank_z(
        BENCHMARK_P,
        P,
        b_grid,
        nz=args.nz,
        rho_z=args.rho_z,
        sigma_z=args.sigma_z,
        verbose=not args.quiet,
    )
    elapsed = time.perf_counter() - t0
    return sol, P_out, p_eq, b_grid, targets, weights, elapsed


def extract_v5_moments(sol: SimpleNamespace, P: SimpleNamespace, p_eq: np.ndarray) -> dict[str, float]:
    moments_ns = extract_moments(sol, P, p_eq, 0.0, 0.0, 0.0, True)
    moments_ns.inv_pop_share_C = float(sol.pop_share[1])
    moments_ns.inv_rent_ratio_C_over_P = float(
        (P.user_cost_rate * p_eq[1]) / max(P.user_cost_rate * p_eq[0], 1e-12)
    )
    return namespace_float_dict(moments_ns)


def write_plot_data(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    rows: list[dict[str, float | str]] = []
    for records in tables.values():
        rows.extend(records)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    path = out_dir / f"plot_data_{OUTPUT_STEM}.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_index(out_dir: Path, images: list[Path], summary: dict[str, object]) -> Path:
    bullets = "\n".join(f"- `{img.name}`" for img in images)
    text = f"""# V5 HANK-z Outside-Closure Nz=5 Figure Set

Generated from the isolated benchmark-normalized outside-option closure branch.

## Run Summary

- accepted: `{summary.get("accepted")}`
- strict converged: `{summary.get("strict_converged")}`
- convergence reason: `{summary.get("convergence_reason")}`
- iterations: `{summary.get("iterations_completed")}`
- final equilibrium error: `{float(summary.get("final_eq_error", np.nan)):.6g}`
- elapsed seconds for solve and plots: `{float(summary.get("elapsed_sec", np.nan)):.2f}`
- prices: `{summary.get("p_eq")}`
- \(S\): `{float(summary.get("scale_factor", np.nan)):.6g}`
- \(q^E\): `{float(summary.get("city_entry_prob_total", np.nan)):.6g}`
- outside probability: `{float(summary.get("outside_entry_prob", np.nan)):.6g}`
- residual outside-born flow \(M\): `{float(summary.get("outside_entry_flow", np.nan)):.6g}`
- calibrated outside value: `{float(summary.get("outside_value", np.nan)):.6g}`
- SMM loss: `{float(summary.get("loss", np.nan)):.6g}`
- \(b\) states: `{summary.get("Nb")}`
- \(z\) states: `{summary.get("Nz")}`

## Images

{bullets}

## PDF

- `{summary.get("pdf")}`

## Notes

This is a separate `Nz=5` plotting/evaluation run. It does not overwrite the
canonical `Nz=7` V5 outside-closure output files.
"""
    path = out_dir / "README_FIGURES.md"
    path.write_text(text)
    return path


def write_pdf_packet(out_dir: Path, images: list[Path], summary: dict[str, object]) -> Path:
    pdf_path = out_dir / "HANK_Z_OUTSIDE_CLOSURE_NZ5_FIGURE_PACKET.pdf"
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis("off")
        lines = [
            "V5 HANK-z Outside-Option Closure",
            "Benchmark-normalized Nz=5 figure packet",
            "",
            f"accepted: {summary.get('accepted')}",
            f"strict converged: {summary.get('strict_converged')}",
            f"convergence reason: {summary.get('convergence_reason')}",
            f"iterations: {summary.get('iterations_completed')}",
            f"final GE error: {float(summary.get('final_eq_error', np.nan)):.6g}",
            f"elapsed seconds: {float(summary.get('elapsed_sec', np.nan)):.2f}",
            f"prices: {summary.get('p_eq')}",
            f"S: {float(summary.get('scale_factor', np.nan)):.6g}",
            f"qE: {float(summary.get('city_entry_prob_total', np.nan)):.6g}",
            f"outside probability: {float(summary.get('outside_entry_prob', np.nan)):.6g}",
            f"M: {float(summary.get('outside_entry_flow', np.nan)):.6g}",
            f"outside value: {float(summary.get('outside_value', np.nan)):.6g}",
            f"SMM loss: {float(summary.get('loss', np.nan)):.6g}",
            f"b states: {summary.get('Nb')}",
            f"z states: {summary.get('Nz')}",
            "",
            "Structural state in this prototype:",
            r"$(b,d,i,a,n,s,z)$",
            "",
            r"The mortgage-account state $\mu$ is not active in this packet.",
        ]
        ax.text(0.08, 0.90, "\n".join(lines), va="top", ha="left", fontsize=14)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for image_path in images:
            img = plt.imread(image_path)
            fig = plt.figure(figsize=(11, 8.5))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(image_path.stem.replace("_", " "), fontsize=13, pad=12)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    return pdf_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=5)
    parser.add_argument("--rho-z", type=float, default=0.95)
    parser.add_argument("--sigma-z", type=float, default=0.35)
    parser.add_argument("--target-city-entry-prob", type=float, default=0.9)
    parser.add_argument("--kappa-entry", type=float, default=1_000_000.0)
    parser.add_argument("--local-birth-entry-weight", type=float, default=1.0)
    parser.add_argument("--max-iter-eq", type=int, default=60)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    setup_plot_style()
    root = Path(__file__).resolve().parent
    out_dir = root / FIGURE_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = root / f"{OUTPUT_STEM}.log"
    try:
        sol, P, p_eq, b_grid, targets, weights, elapsed = solve_v5(args)
        scale = getattr(sol, "accounting_scale", SimpleNamespace())
        moments = extract_v5_moments(sol, P, p_eq)
        loss = compute_smm_loss(moments, targets, weights)

        results_path = root / f"results_{OUTPUT_STEM}.csv"
        diagnostics_path = root / f"diagnostics_{OUTPUT_STEM}.csv"
        closure_path = root / f"diagnostics_{OUTPUT_STEM}_closure.csv"
        trace_path = root / f"diagnostics_{OUTPUT_STEM}_trace.csv"
        report_path = root / "REPORT_V5_HANK_Z_OUTSIDE_CLOSURE_NZ5.md"

        write_results(results_path, targets, moments)
        write_hank_z_diagnostics(diagnostics_path, sol, P)
        write_closure_diagnostics(closure_path, sol)
        write_trace(trace_path, sol)
        write_report(report_path, sol, P, p_eq, moments, targets, loss, elapsed)

        tables = sorting_tables(sol, P, b_grid)
        plot_data_path = write_plot_data(out_dir, tables)
        images = [
            plot_trace(out_dir, trace_path, sol),
            plot_moments(out_dir, results_path),
            plot_moment_gaps(out_dir, results_path),
            plot_z_process(out_dir, sol, P),
            plot_sorting_age_profiles(out_dir, tables),
            plot_center_ownership_heatmaps(out_dir, tables),
            plot_tenure_stacks(out_dir, tables),
            plot_fertility_wealth(out_dir, tables),
            plot_aggregate_sorting(out_dir, tables),
        ]

        summary = {
            "ok": True,
            "accepted": bool(sol.timings.get("accepted", False)),
            "strict_converged": bool(sol.timings.get("strict_converged", False)),
            "convergence_reason": sol.timings.get("convergence_reason"),
            "iterations_completed": int(sol.timings.get("iterations_completed", 0)),
            "best_eq_error": float(sol.timings.get("best_eq_error", np.nan)),
            "final_eq_error": float(sol.timings.get("final_eq_error", np.nan)),
            "elapsed_sec": float(elapsed),
            "loss": float(loss),
            "p_eq": [float(x) for x in p_eq],
            "Nb": int(P.Nb),
            "Nz": int(len(sol.hank_z.z_grid)),
            "rho_z": float(args.rho_z),
            "sigma_z": float(args.sigma_z),
            "kappa_entry": float(args.kappa_entry),
            "target_city_entry_prob": float(args.target_city_entry_prob),
            "scale_factor": float(getattr(scale, "scale_factor", np.nan)),
            "city_entry_prob_total": float(getattr(scale, "city_entry_prob_total", np.nan)),
            "outside_entry_prob": float(getattr(scale, "outside_entry_prob", np.nan)),
            "outside_entry_flow": float(getattr(scale, "outside_entry_flow", np.nan)),
            "outside_value": float(getattr(scale, "outside_value", np.nan)),
            "moments": moments,
            "timings": dict(sol.timings),
            "results": results_path.name,
            "diagnostics": diagnostics_path.name,
            "closure_diagnostics": closure_path.name,
            "trace": trace_path.name,
            "report": report_path.name,
            "plot_data": str(plot_data_path.relative_to(root)),
            "images": [str(img.relative_to(root)) for img in images],
        }
        pdf_path = write_pdf_packet(out_dir, images, summary)
        summary["pdf"] = str(pdf_path.relative_to(root))
        index_path = write_index(out_dir, images, summary)
        summary["figure_index"] = str(index_path.relative_to(root))
        (out_dir / "plot_run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
        log_path.write_text(json.dumps(summary, indent=2, default=str))

        if not args.quiet:
            print(f"Wrote {len(images)} figures to {out_dir}")
            print(f"Wrote PDF packet to {pdf_path}")
            print(f"Wrote plot data to {plot_data_path}")
            print(f"Wrote figure index to {index_path}")
    except Exception as exc:
        log_path.write_text(json.dumps({"ok": False, "error": repr(exc)}, indent=2))
        raise


if __name__ == "__main__":
    main()
