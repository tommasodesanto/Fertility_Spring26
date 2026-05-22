#!/usr/bin/env python3
"""Generate figure set for the V4 full-GE HANK-z prototype."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from dt_cp_model.parameters import apply_overrides, finalize_location_choice_spec, setup_parameters
from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import BENCHMARK_P, build_base_parameters


MOMENT_PLOT_ORDER = [
    "tfr",
    "childless_rate",
    "mean_age_first_birth",
    "tfr_gradient",
    "own_rate",
    "own_gradient",
    "own_family_gap",
    "prime_childless_renter_median_rooms",
    "prime_childless_owner_median_rooms",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "young_liquid_wealth_to_income",
    "center_share_nonparents",
    "center_share_newparents",
    "migration_rate",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
    "inv_pop_share_C",
    "inv_rent_ratio_C_over_P",
]


def setup_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 180,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def savefig(fig: plt.Figure, out_dir: Path, name: str) -> Path:
    path = out_dir / name
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def age_vector(P) -> np.ndarray:
    return float(getattr(P, "age_start", 18)) + np.arange(P.J) * float(getattr(P, "da", 1.0))


def solve_v4(args: argparse.Namespace):
    P_override, _, _ = build_base_parameters(args.nb, force_full=False)
    P = setup_parameters()
    P = apply_overrides(P, P_override)
    P = finalize_location_choice_spec(P)
    P.solve_mode = "ge"
    P.Nb = int(args.nb)
    P.max_iter_eq = int(args.max_iter_eq)
    P.tol_eq = float(args.tol_eq)
    P.collect_ge_trace = True
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
    return sol, P_out, p_eq, b_grid, elapsed


def plot_trace(out_dir: Path, trace_path: Path, sol=None) -> Path:
    if trace_path.exists():
        trace = pd.read_csv(trace_path)
    elif sol is not None and hasattr(sol, "ge_trace"):
        trace = pd.DataFrame(sol.ge_trace)
        if "p" in trace:
            trace["p_P"] = trace["p"].apply(lambda x: x[0])
            trace["p_C"] = trace["p"].apply(lambda x: x[1])
    else:
        raise FileNotFoundError("No GE trace available for plotting.")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    ax = axes[0, 0]
    ax.plot(trace["iter"], trace["err"], marker="o", label="max error")
    ax.plot(trace["iter"], trace["err_p"], marker="o", label="prices")
    ax.plot(trace["iter"], trace["err_e"], marker="o", label="entry")
    ax.set_yscale("log")
    ax.set_xlabel("GE iteration")
    ax.set_ylabel("relative error")
    ax.set_title("Equilibrium Convergence")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(trace["iter"], trace["p_P"], marker="o", label="Periphery")
    ax.plot(trace["iter"], trace["p_C"], marker="o", label="Center")
    ax.set_xlabel("GE iteration")
    ax.set_ylabel("price")
    ax.set_title("House Prices")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(trace["iter"], trace["own_rate"], marker="o", color="#2f6f4e")
    ax.set_xlabel("GE iteration")
    ax.set_ylabel("ownership")
    ax.set_title("Ownership During Fixed Point")

    ax = axes[1, 1]
    ax.plot(trace["iter"], trace["pop_C"], marker="o", label="Center", color="#7b3f98")
    ax.plot(trace["iter"], trace["pop_P"], marker="o", label="Periphery", color="#c46a27")
    ax.set_xlabel("GE iteration")
    ax.set_ylabel("population share")
    ax.set_title("Spatial Population Shares")
    ax.legend(frameon=False)

    return savefig(fig, out_dir, "01_equilibrium_trace.png")


def plot_moments(out_dir: Path, results_path: Path) -> Path:
    df = pd.read_csv(results_path)
    df = df[df["moment"].isin(MOMENT_PLOT_ORDER)].copy()
    df["moment"] = pd.Categorical(df["moment"], MOMENT_PLOT_ORDER, ordered=True)
    df = df.sort_values("moment")
    x = np.arange(len(df))
    width = 0.27

    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.bar(x - width, df["target"], width=width, label="Target", color="#252525")
    ax.bar(x, df["benchmark_model"], width=width, label="Benchmark", color="#9e9e9e")
    ax.bar(x + width, df["model"], width=width, label="HANK-z GE", color="#2a6fbb")
    ax.axhline(0.0, color="black", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("_", "\n") for m in df["moment"]], rotation=0)
    ax.set_title("Targets, Current Benchmark, and HANK-z Full-GE Prototype")
    ax.set_ylabel("moment value")
    ax.legend(frameon=False, ncol=3)
    return savefig(fig, out_dir, "02_moment_target_benchmark_model.png")


def plot_moment_gaps(out_dir: Path, results_path: Path) -> Path:
    df = pd.read_csv(results_path)
    df = df[df["moment"].isin(MOMENT_PLOT_ORDER)].copy()
    df["moment"] = pd.Categorical(df["moment"], MOMENT_PLOT_ORDER, ordered=True)
    df = df.sort_values("moment")
    df["benchmark_gap"] = df["benchmark_model"] - df["target"]
    df["model_gap"] = df["model"] - df["target"]
    y = np.arange(len(df))
    height = 0.36

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(y - height / 2, df["benchmark_gap"], height=height, label="Benchmark - target", color="#9e9e9e")
    ax.barh(y + height / 2, df["model_gap"], height=height, label="HANK-z GE - target", color="#2a6fbb")
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df["moment"])
    ax.invert_yaxis()
    ax.set_xlabel("raw deviation from target")
    ax.set_title("Moment Deviations From Target")
    ax.legend(frameon=False)
    return savefig(fig, out_dir, "02b_moment_deviations_from_target.png")


def plot_z_process(out_dir: Path, sol, P) -> Path:
    z = np.asarray(sol.hank_z.z_grid)
    Pi = np.asarray(sol.hank_z.Pi_z)
    diag = sol.hank_z.diagnostics
    ages = age_vector(P)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    im = ax.imshow(Pi, vmin=0, vmax=max(1.0, float(np.max(Pi))), cmap="Blues")
    ax.set_xticks(np.arange(len(z)))
    ax.set_yticks(np.arange(len(z)))
    ax.set_xticklabels([f"{v:.2f}" for v in z])
    ax.set_yticklabels([f"{v:.2f}" for v in z])
    ax.set_xlabel("next z")
    ax.set_ylabel("current z")
    ax.set_title("Earnings Transition Matrix")
    for i in range(Pi.shape[0]):
        for j in range(Pi.shape[1]):
            ax.text(j, i, f"{Pi[i, j]:.2f}", ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1]
    for iz, zval in enumerate(z):
        ax.plot(ages, diag.z_by_age[:, iz], label=f"z={zval:.2f}", linewidth=2)
    ax.set_xlabel("age")
    ax.set_ylabel("mass share")
    ax.set_title("Distribution of z by Age")
    ax.legend(frameon=False)
    return savefig(fig, out_dir, "03_z_process_and_age_mass.png")


def expected_fertility_by_age_z(sol, P, g_z: np.ndarray) -> np.ndarray:
    Nb, nt, I, J, _, _, Nz = g_z.shape
    out = np.full((J, Nz), np.nan)
    for j in range(J):
        if not ((j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end)):
            continue
        for iz in range(Nz):
            base = g_z[:, :, :, j, 0, 0, iz]
            mass = float(np.sum(base))
            if mass <= 1e-14:
                continue
            pr = sol.hank_z.fert_probs_z[:, :, :, j, :, iz]
            out[j, iz] = float(np.sum(base * (pr @ np.arange(P.n_parity))) / mass)
    return out


def mean_housing_by_age_z(sol, P, g_z: np.ndarray) -> np.ndarray:
    Nb, nt, I, J, npar, ncs, Nz = g_z.shape
    out = np.full((J, Nz), np.nan)
    hR = sol.hank_z.hR_pol_z
    for j in range(J):
        for iz in range(Nz):
            total_mass = float(np.sum(g_z[:, :, :, j, :, :, iz]))
            if total_mass <= 1e-14:
                continue
            total_h = 0.0
            renter_mass = g_z[:, 0, :, j, :, :, iz]
            total_h += float(np.sum(renter_mass * hR[:, 0, :, j, :, :, iz]))
            for ten in range(1, nt):
                mt = float(np.sum(g_z[:, ten, :, j, :, :, iz]))
                total_h += mt * float(P.H_own[ten - 1])
            out[j, iz] = total_h / total_mass
    return out


def sorting_tables(sol, P, b_grid: np.ndarray) -> dict[str, list[dict[str, float | str]]]:
    g_z = np.asarray(sol.hank_z.g_z)
    z_grid = np.asarray(sol.hank_z.z_grid)
    diag = sol.hank_z.diagnostics
    ages = age_vector(P)
    center_idx = min(1, P.I - 1)
    fert_age_z = expected_fertility_by_age_z(sol, P, g_z)
    housing_age_z = mean_housing_by_age_z(sol, P, g_z)

    age_z_rows: list[dict[str, float | str]] = []
    aggregate_rows: list[dict[str, float | str]] = []
    tenure_rows: list[dict[str, float | str]] = []

    for iz, zval in enumerate(z_grid):
        gz_all = g_z[..., iz]
        mass_all = float(np.sum(gz_all))
        aggregate_rows.append(
            {
                "diagnostic": "aggregate_z",
                "z": float(zval),
                "mass": mass_all,
                "center_share": float(np.sum(gz_all[:, :, center_idx, :, :, :]) / max(mass_all, 1e-14)),
                "own_rate": float(np.sum(gz_all[:, 1:, :, :, :, :]) / max(mass_all, 1e-14)),
                "mean_liquid_wealth": float(np.sum(gz_all * b_grid.reshape(-1, 1, 1, 1, 1, 1)) / max(mass_all, 1e-14)),
            }
        )
        for j, age in enumerate(ages):
            gjz = g_z[:, :, :, j, :, :, iz]
            mass = float(np.sum(gjz))
            if mass <= 1e-14:
                continue
            age_z_rows.append(
                {
                    "diagnostic": "age_z",
                    "age": float(age),
                    "z": float(zval),
                    "mass": mass,
                    "mass_share": float(diag.z_by_age[j, iz]),
                    "center_share": float(np.sum(gjz[:, :, center_idx, :, :]) / mass),
                    "own_rate": float(np.sum(gjz[:, 1:, :, :, :]) / mass),
                    "mean_liquid_wealth": float(diag.mean_liquid_wealth_by_z_age[j, iz]),
                    "mean_income": float(diag.mean_income_by_z_age[j, iz]),
                    "mean_housing": float(housing_age_z[j, iz]),
                    "expected_fertility_choice": float(fert_age_z[j, iz]) if np.isfinite(fert_age_z[j, iz]) else np.nan,
                }
            )

    age_groups = {
        "25-34": (25, 34),
        "35-44": (35, 44),
        "45-54": (45, 54),
        "55-64": (55, 64),
    }
    for group, (lo, hi) in age_groups.items():
        jidx = np.where((ages >= lo) & (ages <= hi))[0]
        for iz, zval in enumerate(z_grid):
            block = np.take(g_z[..., iz], jidx, axis=3)
            mass = float(np.sum(block))
            if mass <= 1e-14:
                continue
            renter = float(np.sum(block[:, 0, :, :, :, :]) / mass)
            tenure_rows.append({"diagnostic": "tenure_agegroup_z", "age_group": group, "z": float(zval), "tenure": "renter", "share": renter})
            for ten in range(1, 1 + P.n_house):
                share = float(np.sum(block[:, ten, :, :, :, :]) / mass)
                tenure_rows.append(
                    {
                        "diagnostic": "tenure_agegroup_z",
                        "age_group": group,
                        "z": float(zval),
                        "tenure": f"owner_H{ten}",
                        "share": share,
                    }
                )

    fert_wealth_rows = []
    for rec in diag.fertility_by_z_wealth_bin:
        row = {"diagnostic": "fertility_by_z_wealth_bin"}
        row.update(rec)
        fert_wealth_rows.append(row)

    return {
        "age_z": age_z_rows,
        "aggregate_z": aggregate_rows,
        "tenure_agegroup_z": tenure_rows,
        "fertility_by_z_wealth_bin": fert_wealth_rows,
    }


def write_plot_data(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    rows = []
    for records in tables.values():
        rows.extend(records)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    path = out_dir / "plot_data_income_mortgage_risk_v4_hank_z_ge.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def plot_sorting_age_profiles(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    age_df = pd.DataFrame(tables["age_z"])
    z_vals = sorted(age_df["z"].unique())
    colors = ["#8c510a", "#4d4d4d", "#01665e", "#2a6fbb", "#7b3f98"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    panels = [
        ("mean_liquid_wealth", "Mean Liquid Wealth", "wealth"),
        ("own_rate", "Ownership by z and Age", "ownership"),
        ("center_share", "Center Share by z and Age", "center share"),
        ("mean_housing", "Mean Housing Services by z and Age", "rooms/services"),
    ]
    for ax, (col, title, ylabel) in zip(axes.ravel(), panels):
        for idx, zval in enumerate(z_vals):
            d = age_df[age_df["z"] == zval]
            ax.plot(d["age"], d[col], linewidth=2, label=f"z={zval:.2f}", color=colors[idx % len(colors)])
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("age")
    axes[0, 0].legend(frameon=False)
    return savefig(fig, out_dir, "04_sorting_age_profiles.png")


def plot_center_ownership_heatmaps(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    age_df = pd.DataFrame(tables["age_z"])
    z_vals = sorted(age_df["z"].unique())
    ages = sorted(age_df["age"].unique())
    labels = [f"{z:.2f}" for z in z_vals]

    center = np.full((len(z_vals), len(ages)), np.nan)
    own = np.full_like(center, np.nan)
    fertility = np.full_like(center, np.nan)
    for i, zval in enumerate(z_vals):
        d = age_df[age_df["z"] == zval].set_index("age")
        center[i, :] = d.loc[ages, "center_share"]
        own[i, :] = d.loc[ages, "own_rate"]
        fertility[i, :] = d.loc[ages, "expected_fertility_choice"]

    fig, axes = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True)
    for ax, arr, title, cmap in [
        (axes[0], center, "Center Share", "Purples"),
        (axes[1], own, "Ownership Rate", "Greens"),
        (axes[2], fertility, "Expected Fertility Choice During Fertile Ages", "Oranges"),
    ]:
        im = ax.imshow(arr, aspect="auto", interpolation="nearest", cmap=cmap)
        ax.set_yticks(np.arange(len(z_vals)))
        ax.set_yticklabels(labels)
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    step = 5
    axes[-1].set_xticks(np.arange(0, len(ages), step))
    axes[-1].set_xticklabels([f"{int(a)}" for a in ages[::step]])
    axes[-1].set_xlabel("age")
    return savefig(fig, out_dir, "05_sorting_heatmaps.png")


def plot_tenure_stacks(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    df = pd.DataFrame(tables["tenure_agegroup_z"])
    groups = list(dict.fromkeys(df["age_group"].tolist()))
    tenures = list(dict.fromkeys(df["tenure"].tolist()))
    z_vals = sorted(df["z"].unique())
    colors = ["#4d4d4d", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharey=True)
    for ax, group in zip(axes.ravel(), groups):
        bottom = np.zeros(len(z_vals))
        for k, tenure in enumerate(tenures):
            vals = []
            for zval in z_vals:
                rec = df[(df["age_group"] == group) & (df["z"] == zval) & (df["tenure"] == tenure)]
                vals.append(float(rec["share"].iloc[0]) if not rec.empty else 0.0)
            ax.bar([f"{z:.2f}" for z in z_vals], vals, bottom=bottom, label=tenure, color=colors[k % len(colors)])
            bottom += np.asarray(vals)
        ax.set_title(f"Tenure/Rung Shares, Ages {group}")
        ax.set_xlabel("z")
        ax.set_ylabel("share")
        ax.set_ylim(0, 1)
    axes[0, 0].legend(frameon=False, ncol=2, bbox_to_anchor=(0.0, 1.42), loc="upper left")
    return savefig(fig, out_dir, "06_tenure_rung_sorting_by_z.png")


def plot_fertility_wealth(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    df = pd.DataFrame(tables["fertility_by_z_wealth_bin"])
    z_vals = sorted(df["z"].unique())
    bins = sorted(df["wealth_bin"].unique())
    arr = np.full((len(z_vals), len(bins)), np.nan)
    mass = np.full_like(arr, np.nan)
    for i, zval in enumerate(z_vals):
        for j, wb in enumerate(bins):
            rec = df[(df["z"] == zval) & (df["wealth_bin"] == wb)]
            if not rec.empty:
                arr[i, j] = float(rec["expected_completed_fertility_choice"].iloc[0])
                mass[i, j] = float(rec["mass"].iloc[0])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    im = axes[0].imshow(arr, aspect="auto", interpolation="nearest", cmap="YlGnBu")
    axes[0].set_title("Expected Fertility Choice")
    axes[0].set_yticks(np.arange(len(z_vals)))
    axes[0].set_yticklabels([f"{z:.2f}" for z in z_vals])
    axes[0].set_xticks(np.arange(len(bins)))
    axes[0].set_xticklabels([f"bin {int(b)}" for b in bins])
    axes[0].set_ylabel("z")
    axes[0].set_xlabel("liquid-wealth bin")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    im = axes[1].imshow(mass, aspect="auto", interpolation="nearest", cmap="YlOrBr")
    axes[1].set_title("Fertile Childless Mass")
    axes[1].set_yticks(np.arange(len(z_vals)))
    axes[1].set_yticklabels([f"{z:.2f}" for z in z_vals])
    axes[1].set_xticks(np.arange(len(bins)))
    axes[1].set_xticklabels([f"bin {int(b)}" for b in bins])
    axes[1].set_xlabel("liquid-wealth bin")
    fig.colorbar(im, ax=axes[1], fraction=0.046)
    return savefig(fig, out_dir, "07_fertility_by_z_and_wealth.png")


def plot_aggregate_sorting(out_dir: Path, tables: dict[str, list[dict[str, float | str]]]) -> Path:
    df = pd.DataFrame(tables["aggregate_z"])
    x = np.arange(len(df))
    labels = [f"{z:.2f}" for z in df["z"]]
    width = 0.25

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.bar(x - width, df["center_share"], width, label="Center share", color="#7b3f98")
    ax.bar(x, df["own_rate"], width, label="Ownership", color="#2f6f4e")
    ax.bar(x + width, df["mass"] / max(float(df["mass"].sum()), 1e-14), width, label="Population mass", color="#777777")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("z")
    ax.set_ylabel("share")
    ax.set_title("Aggregate Sorting by z")
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(df["z"], df["mean_liquid_wealth"], marker="o", linewidth=2, color="#2a6fbb")
    ax.set_xlabel("z")
    ax.set_ylabel("mean liquid wealth")
    ax.set_title("Wealth Gradient by z")
    return savefig(fig, out_dir, "08_aggregate_sorting_by_z.png")


def write_index(out_dir: Path, images: list[Path], log_summary: dict) -> Path:
    bullets = "\n".join(f"- `{img.name}`" for img in images)
    text = f"""# V4 HANK-z GE Figure Set

Generated from the isolated full-equilibrium HANK-\(z\) prototype.

## Run Summary

- accepted: `{log_summary.get("accepted")}`
- strict converged: `{log_summary.get("strict_converged")}`
- iterations: `{log_summary.get("iterations_completed")}`
- final equilibrium error: `{log_summary.get("final_eq_error")}`
- elapsed seconds for plot solve: `{log_summary.get("plot_elapsed_sec"):.2f}`
- \(b\) states: `{log_summary.get("Nb")}`
- \(z\) states: `{log_summary.get("Nz")}`

## Images

{bullets}

## PDF

- `{log_summary.get("pdf", "HANK_Z_GE_FIGURE_PACKET.pdf")}`

## Notes

These figures use the copied prototype only. The structural state is
\((b,d,i,a,n,s,z)\). The mortgage-account state \(\mu\) is not active in this
full-GE figure set.
"""
    path = out_dir / "README_FIGURES.md"
    path.write_text(text)
    return path


def write_pdf_packet(out_dir: Path, images: list[Path], log_summary: dict) -> Path:
    """Compile the PNG figure set into a single PDF packet."""

    pdf_path = out_dir / "HANK_Z_GE_FIGURE_PACKET.pdf"
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis("off")
        lines = [
            "HANK-z Full-Equilibrium Prototype",
            "",
            "Isolated overnight branch figure packet",
            "",
            f"accepted: {log_summary.get('accepted')}",
            f"strict converged: {log_summary.get('strict_converged')}",
            f"iterations: {log_summary.get('iterations_completed')}",
            f"final GE error: {log_summary.get('final_eq_error'):.6g}",
            f"plot solve seconds: {log_summary.get('plot_elapsed_sec'):.2f}",
            f"b states: {log_summary.get('Nb')}",
            f"z states: {log_summary.get('Nz')}",
            "",
            "Structural state in this prototype:",
            r"$(b,d,i,a,n,s,z)$",
            "",
            r"The mortgage-account state $\mu$ is not active in this full-GE packet.",
        ]
        ax.text(0.08, 0.90, "\n".join(lines), va="top", ha="left", fontsize=15)
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
    parser.add_argument("--nz", type=int, default=3, choices=[3, 5])
    parser.add_argument("--rho-z", type=float, default=0.82)
    parser.add_argument("--sigma-z", type=float, default=0.28)
    parser.add_argument("--max-iter-eq", type=int, default=35)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    setup_plot_style()
    root = Path(__file__).resolve().parent
    out_dir = root / "figures_v4_hank_z_ge"
    out_dir.mkdir(parents=True, exist_ok=True)

    sol, P, _, b_grid, elapsed = solve_v4(args)
    tables = sorting_tables(sol, P, b_grid)
    data_path = write_plot_data(out_dir, tables)

    images = [
        plot_trace(out_dir, root / "diagnostics_income_mortgage_risk_v4_hank_z_ge_trace.csv", sol),
        plot_moments(out_dir, root / "results_income_mortgage_risk_v4_hank_z_ge.csv"),
        plot_moment_gaps(out_dir, root / "results_income_mortgage_risk_v4_hank_z_ge.csv"),
        plot_z_process(out_dir, sol, P),
        plot_sorting_age_profiles(out_dir, tables),
        plot_center_ownership_heatmaps(out_dir, tables),
        plot_tenure_stacks(out_dir, tables),
        plot_fertility_wealth(out_dir, tables),
        plot_aggregate_sorting(out_dir, tables),
    ]

    summary = {
        "accepted": bool(sol.timings.get("accepted", False)),
        "strict_converged": bool(sol.timings.get("strict_converged", False)),
        "iterations_completed": int(sol.timings.get("iterations_completed", 0)),
        "final_eq_error": float(sol.timings.get("final_eq_error", np.nan)),
        "plot_elapsed_sec": float(elapsed),
        "Nb": int(P.Nb),
        "Nz": int(len(sol.hank_z.z_grid)),
        "plot_data": str(data_path.name),
        "images": [img.name for img in images],
    }
    pdf_path = write_pdf_packet(out_dir, images, summary)
    summary["pdf"] = pdf_path.name
    (out_dir / "plot_run_summary.json").write_text(json.dumps(summary, indent=2))
    index_path = write_index(out_dir, images, summary)

    if not args.quiet:
        print(f"Wrote {len(images)} figures to {out_dir}")
        print(f"Wrote PDF packet to {pdf_path}")
        print(f"Wrote plot data to {data_path}")
        print(f"Wrote figure index to {index_path}")


if __name__ == "__main__":
    main()
