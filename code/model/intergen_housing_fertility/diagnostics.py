"""Diagnostic output for the intergenerational housing/fertility model."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


def write_diagnostics(sol: SimpleNamespace, P: SimpleNamespace, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "summary.json").write_text(_dumps_json(_summary(sol, P)))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ages = P.age_start + np.arange(P.J) * P.da

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ages, sol.own_by_age, lw=2.0)
    ax.axvspan(
        P.age_start + (P.A_f_start - 1) * P.da,
        P.age_start + P.A_f_end * P.da,
        color="0.9",
        label="fertile window",
    )
    ax.set_xlabel("age")
    ax.set_ylabel("ownership rate")
    ax.set_ylim(0.0, 1.05)
    ax.legend(frameon=False)
    ax.set_title("Ownership by age")
    fig.tight_layout()
    fig.savefig(outdir / "ownership_by_age.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ages, getattr(sol, "fert_by_age", np.zeros(P.J)), lw=2.0, color="tab:green")
    ax.set_xlabel("age")
    ax.set_ylabel("expected births")
    ax.set_title("Fertility by age")
    fig.tight_layout()
    fig.savefig(outdir / "fertility_by_age.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([0, 1], [float(sol.aggregate_housing_demand), float(sol.aggregate_housing_supply)])
    ax.set_xticks([0, 1], ["demand", "supply"])
    ax.set_ylabel("service units per adult")
    rel_excess = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
    ax.set_title(f"Housing-services clearing, max rel. residual={rel_excess:.2e}")
    fig.tight_layout()
    fig.savefig(outdir / "housing_market.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    markets = np.arange(P.I)
    demand = np.asarray(getattr(sol, "housing_demand", np.zeros(P.I)), dtype=float).reshape(-1)
    supply = np.asarray(getattr(sol, "housing_supply", np.zeros(P.I)), dtype=float).reshape(-1)
    width = 0.35
    ax.bar(markets - width / 2, demand, width, label="demand")
    ax.bar(markets + width / 2, supply, width, label="supply")
    ax.set_xticks(markets, [str(i + 1) for i in markets])
    ax.set_xlabel("housing-services market")
    ax.set_ylabel("service units per adult")
    ax.legend(frameon=False)
    ax.set_title("Housing market by market")
    fig.tight_layout()
    fig.savefig(outdir / "market_clearing_by_market.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    residual = (demand - supply) / np.maximum(supply, 1e-12)
    ax.axhline(0.0, color="0.25", lw=1.0)
    ax.bar(markets, residual, color="tab:red")
    ax.set_xticks(markets, [str(i + 1) for i in markets])
    ax.set_xlabel("housing-services market")
    ax.set_ylabel("(demand - supply) / supply")
    ax.set_title("Housing market relative residual")
    fig.tight_layout()
    fig.savefig(outdir / "market_clearing_residuals.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([0, 1], [float(sol.aggregate_rental_demand), float(sol.aggregate_owner_demand)])
    ax.set_xticks([0, 1], ["renters", "owners"])
    ax.set_ylabel("service units per adult")
    ax.set_title("Housing services by tenure")
    fig.tight_layout()
    fig.savefig(outdir / "tenure_services.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    labels = [f"{h:g}" for h in P.H_own]
    ax.bar(np.arange(P.n_house), sol.owner_demand_by_size)
    ax.set_xticks(np.arange(P.n_house), labels)
    ax.set_xlabel("owner housing rung")
    ax.set_ylabel("service units per adult")
    ax.set_title("Owner demand by housing rung")
    fig.tight_layout()
    fig.savefig(outdir / "owner_rungs.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.bar([0], [float(sol.owner_user_cost[0])], width=0.45, label="user cost")
    ax1.set_xticks([0], ["aggregate"])
    ax1.set_ylabel("flow user cost")
    ax2 = ax1.twinx()
    ax2.plot([0], [float(sol.owner_asset_price[0])], marker="s", color="tab:orange", label="asset price")
    ax2.set_ylabel("asset price")
    ax1.set_title("Housing price")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "housing_prices.png", dpi=180)
    plt.close(fig)

    if hasattr(sol, "type_values"):
        z = np.asarray(sol.type_values, dtype=float)
        x = np.arange(len(z))
        width = 0.35
        fig, ax1 = plt.subplots(figsize=(7, 4))
        ax1.bar(x - width / 2, getattr(sol, "own_rate_by_income_type", np.zeros_like(z)), width, label="ownership")
        ax1.set_ylabel("ownership rate")
        ax1.set_ylim(0.0, 1.05)
        ax1.set_xticks(x, [f"{v:g}" for v in z])
        ax1.set_xlabel("income state")
        ax2 = ax1.twinx()
        ax2.bar(x + width / 2, getattr(sol, "mean_fertility_by_income_type", np.zeros_like(z)), width, color="tab:green", label="children")
        ax2.set_ylabel("mean completed children")
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.legend(h1 + h2, l1 + l2, frameon=False)
        ax1.set_title("Outcomes by income state")
        fig.tight_layout()
        fig.savefig(outdir / "income_state_outcomes.png", dpi=180)
        plt.close(fig)

    write_distribution_diagnostics(sol, P, outdir, plt, ages)
    write_policy_diagnostics(sol, P, outdir, plt)


def write_distribution_diagnostics(sol: SimpleNamespace, P: SimpleNamespace, outdir: Path, plt: Any, ages: np.ndarray) -> None:
    if not hasattr(sol, "g") or not hasattr(sol, "b_grid"):
        return
    g = np.asarray(sol.g)
    if g.ndim < 7:
        return
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))), dtype=float).reshape(-1)
    if len(z_grid) != g.shape[4]:
        return

    own_age_z = np.full((P.J, len(z_grid)), np.nan)
    wealth_age_z = np.full((P.J, len(z_grid)), np.nan)
    housing_age_z = np.full((P.J, len(z_grid)), np.nan)
    fert_age_z = np.full((P.J, len(z_grid)), np.nan)
    for j in range(P.J):
        for zz in range(len(z_grid)):
            gjz = g[:, :, :, j, zz, :, :]
            mass = float(np.sum(gjz))
            if mass <= 1e-14:
                continue
            own_age_z[j, zz] = float(np.sum(gjz[:, 1:, :, :, :]) / mass)
            wealth_age_z[j, zz] = weighted_liquid_wealth(gjz, b_grid)
            housing_age_z[j, zz] = mean_housing_for_distribution(gjz, sol, P, j, zz)
            if P.A_f_start - 1 <= j <= P.A_f_end and hasattr(sol, "fert_probs"):
                childless = g[:, :, :, j, zz, 0, 0]
                cmass = float(np.sum(childless))
                if cmass > 1e-14:
                    probs = np.asarray(sol.fert_probs)[:, :, :, j, zz, :]
                    fert_age_z[j, zz] = float(np.sum(childless * (probs @ np.arange(P.n_parity))) / cmass)

    plot_age_by_income_state(outdir, plt, ages, z_grid, own_age_z, "ownership rate", "Ownership by age and income state", "ownership_by_age_income_state.png", ylim=(0.0, 1.05))
    plot_age_by_income_state(outdir, plt, ages, z_grid, wealth_age_z, "mean liquid wealth", "Liquid wealth by age and income state", "liquid_wealth_by_age_income_state.png")
    plot_age_by_income_state(outdir, plt, ages, z_grid, housing_age_z, "mean housing services", "Housing by age and income state", "housing_by_age_income_state.png")
    plot_age_by_income_state(outdir, plt, ages, z_grid, fert_age_z, "expected children", "Fertility policy by age and income state", "fertility_policy_by_age_income_state.png")

    ages_to_plot = np.asarray(getattr(P, "diagnostic_policy_ages", np.array([30.0, 42.0])), dtype=float).reshape(-1)
    for age in ages_to_plot:
        j = age_to_index(P, float(age))
        fig, ax = plt.subplots(figsize=(7, 4))
        for zz, z_value in enumerate(z_grid):
            mass = g[:, 0, 0, j, zz, 0, 0]
            if float(np.sum(mass)) <= 1e-14:
                continue
            ax.plot(b_grid, mass / max(float(np.sum(mass)), 1e-14), lw=1.8, label=f"z={z_value:g}")
        ax.set_xlabel("liquid wealth")
        ax.set_ylabel("share within state")
        ax.set_title(f"Childless-renter wealth distribution at age {P.age_start + j * P.da:g}")
        ax.legend(frameon=False)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(outdir / f"wealth_dist_childless_renter_age{int(round(P.age_start + j * P.da))}.png", dpi=180)
        plt.close(fig)


def weighted_liquid_wealth(g_age_z: np.ndarray, b_grid: np.ndarray) -> float:
    mass_b = np.sum(g_age_z, axis=(1, 2, 3, 4))
    total = float(np.sum(mass_b))
    return float(np.sum(mass_b * b_grid) / max(total, 1e-14))


def mean_housing_for_distribution(g_age_z: np.ndarray, sol: SimpleNamespace, P: SimpleNamespace, j: int, zz: int) -> float:
    nt = 1 + P.n_house
    total_h = 0.0
    total_m = 0.0
    for ten in range(nt):
        for nn in range(P.n_parity):
            for cs in range(P.n_child_states):
                mass = g_age_z[:, ten, :, nn, cs]
                mh = float(np.sum(mass))
                if mh <= 1e-15:
                    continue
                if ten == 0:
                    total_h += float(np.sum(mass[:, 0] * sol.hR_pol[:, ten, 0, j, zz, nn, cs]))
                else:
                    total_h += mh * float(P.H_own[ten - 1])
                total_m += mh
    return total_h / max(total_m, 1e-14)


def plot_age_by_income_state(
    outdir: Path,
    plt: Any,
    ages: np.ndarray,
    z_grid: np.ndarray,
    values: np.ndarray,
    ylabel: str,
    title: str,
    filename: str,
    ylim: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    for zz, z_value in enumerate(z_grid):
        ax.plot(ages, values[:, zz], lw=1.8, label=f"z={z_value:g}")
    ax.set_xlabel("age")
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=180)
    plt.close(fig)


def write_policy_diagnostics(sol: SimpleNamespace, P: SimpleNamespace, outdir: Path, plt: Any) -> None:
    if not hasattr(sol, "b_grid"):
        return
    if not hasattr(sol, "c_pol") or not hasattr(sol, "hR_pol") or not hasattr(sol, "fert_probs"):
        return
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))), dtype=float).reshape(-1)
    Nz = len(z_grid)
    if np.asarray(sol.c_pol).ndim < 7 or Nz <= 0:
        return

    ages_to_plot = getattr(P, "diagnostic_policy_ages", np.array([30.0, 42.0]))
    ages_to_plot = np.asarray(ages_to_plot, dtype=float).reshape(-1)
    for age in ages_to_plot:
        j = age_to_index(P, float(age))
        if j < 0 or j >= P.J:
            continue
        plot_policy_childless_renter(sol, P, outdir, plt, b_grid, z_grid, j)


def age_to_index(P: SimpleNamespace, age: float) -> int:
    idx = int(round((float(age) - float(P.age_start)) / max(float(P.da), 1e-12)))
    return int(np.clip(idx, 0, P.J - 1))


def plot_policy_childless_renter(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    outdir: Path,
    plt: Any,
    b_grid: np.ndarray,
    z_grid: np.ndarray,
    j: int,
) -> None:
    i = 0
    ten = 0
    nn = 0
    cs = 0
    age = float(P.age_start + j * P.da)
    nvec = np.arange(P.n_parity)
    c_pol = np.asarray(sol.c_pol)
    hR_pol = np.asarray(sol.hR_pol)
    tc = np.asarray(sol.tenure_choice)
    fp = np.asarray(sol.fert_probs)
    tp = getattr(sol, "tenure_probs", None)
    V = np.asarray(getattr(sol, "V", np.empty_like(c_pol)))

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    axes = axes.ravel()
    for zz, z_value in enumerate(z_grid):
        label = f"z={z_value:g}"
        c_line = c_pol[:, ten, i, j, zz, nn, cs]
        tchoice = tc[:, ten, i, j, zz, nn, cs]
        h_line = np.where(
            tchoice <= 0,
            hR_pol[:, ten, i, j, zz, nn, cs],
            np.asarray(P.H_own, dtype=float)[np.maximum(tchoice - 1, 0)],
        )
        fertility_line = fp[:, ten, i, j, zz, :] @ nvec
        if tp is None:
            owner_line = (tchoice > 0).astype(float)
        else:
            owner_line = np.sum(tp[:, ten, i, j, zz, nn, cs, 1:], axis=1)

        valid = V[:, ten, i, j, zz, nn, cs] > -1e9
        c_line = np.where(valid, c_line, np.nan)
        h_line = np.where(valid, h_line, np.nan)
        fertility_line = np.where(valid, fertility_line, np.nan)
        owner_line = np.where(valid, owner_line, np.nan)

        axes[0].plot(b_grid, c_line, lw=1.8, label=label)
        axes[1].plot(b_grid, h_line, lw=1.8, label=label)
        axes[2].plot(b_grid, fertility_line, lw=1.8, label=label)
        axes[3].plot(b_grid, owner_line, lw=1.8, label=label)

    axes[0].set_ylabel("consumption")
    axes[0].set_title("Consumption")
    axes[1].set_ylabel("housing services")
    axes[1].set_title("Housing after tenure choice")
    axes[2].set_ylabel("expected children")
    axes[2].set_xlabel("liquid wealth")
    axes[2].set_title("Fertility choice")
    axes[3].set_ylabel("probability")
    axes[3].set_xlabel("liquid wealth")
    axes[3].set_ylim(-0.05, 1.05)
    axes[3].set_title("Owner-entry policy")
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, ncols=max(1, min(3, len(z_grid))))
    fig.suptitle(f"Childless renter policies at age {age:g}", y=0.995)
    fig.tight_layout()
    fig.savefig(outdir / f"policy_childless_renter_age{int(round(age))}.png", dpi=180)
    plt.close(fig)


def _summary(sol: SimpleNamespace, P: SimpleNamespace) -> dict[str, Any]:
    return {
        "J": P.J,
        "period_years": getattr(P, "period_years", P.da),
        "n_child_stages": P.n_child_stages,
        "markets": P.I,
        "income_process": str(getattr(P, "income_type_transition", "none")),
        "income_states": getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))),
        "income_state_weights": getattr(sol, "type_weights", getattr(P, "z_weights", np.array([1.0]))),
        "income_transition": getattr(sol, "income_transition", getattr(P, "Pi_z", np.eye(1))),
        "income_types": getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))),
        "income_type_weights": getattr(sol, "type_weights", getattr(P, "z_weights", np.array([1.0]))),
        "pti_constraint": bool(getattr(P, "use_pti_constraint", False)),
        "best_max_abs_rel_excess": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "own_rate": float(sol.own_rate),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "tfr": 2.0 * float(getattr(sol, "mean_completed_fertility", np.nan)),
        "mean_completed_fertility": float(getattr(sol, "mean_completed_fertility", np.nan)),
        "tfr_by_income_type": 2.0 * getattr(sol, "mean_fertility_by_income_type", np.array([])),
        "own_rate_by_income_type": getattr(sol, "own_rate_by_income_type", np.array([])),
        "mean_fertility_by_income_type": getattr(sol, "mean_fertility_by_income_type", np.array([])),
        "housing_demand_by_income_type": getattr(sol, "housing_demand_by_income_type", np.array([])),
        "childless_rate": float(getattr(sol, "childless_rate", np.nan)),
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "owner_user_cost": sol.owner_user_cost,
        "owner_asset_price": sol.owner_asset_price,
        "owner_demand_by_size": sol.owner_demand_by_size,
        "rental_demand_by_market": sol.rental_demand_by_market,
        "housing_supply": sol.housing_supply,
        "aggregate_owner_demand": sol.aggregate_owner_demand,
        "aggregate_rental_demand": sol.aggregate_rental_demand,
        "aggregate_housing_demand": sol.aggregate_housing_demand,
        "aggregate_housing_supply": sol.aggregate_housing_supply,
        "aggregate_housing_excess": sol.aggregate_housing_excess,
    }


def _dumps_json(obj: Any) -> str:
    import json

    def convert(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): convert(v) for k, v in value.items()}
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (list, tuple)):
            return [convert(v) for v in value]
        return value

    return json.dumps(convert(obj), indent=2, sort_keys=True)
