"""Diagnostic output for the first-pass intergenerational housing model."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .utils import dumps_json


def write_diagnostics(sol: SimpleNamespace, P: SimpleNamespace, outdir: Path) -> None:
    """Write a small diagnostic packet for a solved model."""

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "summary.json").write_text(dumps_json(_summary(sol, P)))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ages = P.age_start + np.arange(P.J)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ages, sol.own_by_age, lw=2.0)
    ax.axvline(P.fertility_choice_age, color="0.4", ls="--", lw=1.0, label="fertility choice age")
    ax.axvline(P.old_retention_age, color="0.6", ls=":", lw=1.0, label="old-retention age")
    ax.set_xlabel("age")
    ax.set_ylabel("ownership rate")
    ax.set_ylim(0.0, 1.05)
    ax.legend(frameon=False)
    ax.set_title("Ownership by age")
    fig.tight_layout()
    fig.savefig(outdir / "ownership_by_age.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ages, sol.children_by_age, lw=2.0, color="tab:green")
    ax.axvline(P.fertility_choice_age, color="0.4", ls="--", lw=1.0)
    ax.set_xlabel("age")
    ax.set_ylabel("mean completed children")
    ax.set_title("Children by age")
    fig.tight_layout()
    fig.savefig(outdir / "children_by_age.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([0, 1], [float(sol.aggregate_housing_demand), float(sol.aggregate_housing_supply)])
    ax.set_xticks([0, 1], ["demand", "supply"])
    ax.set_ylabel("service units per adult")
    ax.set_title("Aggregate housing-services clearing")
    fig.tight_layout()
    fig.savefig(outdir / "housing_market.png", dpi=180)
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
    renter_labels = [f"R {h:g}" for h in P.renter_h]
    owner_labels = [f"O {h:g}" for h in P.owner_h]
    quantity_labels = renter_labels + owner_labels
    quantity_demands = np.concatenate([sol.rental_demand_by_size, sol.owner_demand_by_size])
    ax.bar(np.arange(len(quantity_demands)), quantity_demands)
    ax.set_xticks(np.arange(len(quantity_demands)), quantity_labels)
    ax.set_ylabel("service units per adult")
    ax.set_title("Housing demand by quantity and tenure")
    fig.tight_layout()
    fig.savefig(outdir / "housing_quantities.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.bar([0], [float(sol.owner_user_cost[0])], width=0.45, label="user cost")
    ax1.set_xticks([0], ["aggregate"])
    ax1.set_ylabel("flow user cost")
    ax2 = ax1.twinx()
    ax2.plot([0], [float(sol.owner_asset_price[0])], marker="s", color="tab:orange", label="asset price")
    ax2.set_ylabel("asset price")
    ax1.set_title("Housing price")
    lines = ax1.get_lines() + ax2.get_lines()
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles + lines, labels + [line.get_label() for line in lines], frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "housing_prices.png", dpi=180)
    plt.close(fig)


def _summary(sol: SimpleNamespace, P: SimpleNamespace) -> dict:
    return {
        "mode": P.mode,
        "converged": bool(getattr(sol, "converged", False)),
        "best_max_abs_rel_excess": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "own_rate": float(sol.own_rate),
        "young_owner_rate": float(sol.young_owner_rate),
        "old_owner_rate": float(sol.old_owner_rate),
        "mean_completed_fertility": float(sol.mean_completed_fertility),
        "childless_rate": float(sol.childless_rate),
        "owner_user_cost": sol.owner_user_cost,
        "owner_asset_price": sol.owner_asset_price,
        "owner_demand_by_size": sol.owner_demand_by_size,
        "rental_demand_by_size": sol.rental_demand_by_size,
        "housing_supply": sol.housing_supply,
        "aggregate_owner_demand": sol.aggregate_owner_demand,
        "aggregate_rental_demand": sol.aggregate_rental_demand,
        "aggregate_housing_demand": sol.aggregate_housing_demand,
        "aggregate_housing_supply": sol.aggregate_housing_supply,
        "aggregate_housing_excess": sol.aggregate_housing_excess,
    }
