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
    ax.set_title("Housing-services clearing")
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


def _summary(sol: SimpleNamespace, P: SimpleNamespace) -> dict[str, Any]:
    return {
        "J": P.J,
        "period_years": getattr(P, "period_years", P.da),
        "n_child_stages": P.n_child_stages,
        "markets": P.I,
        "pti_constraint": bool(getattr(P, "use_pti_constraint", False)),
        "best_max_abs_rel_excess": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "own_rate": float(sol.own_rate),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "mean_completed_fertility": float(getattr(sol, "mean_completed_fertility", np.nan)),
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
