#!/usr/bin/env python3
"""Build one compact decision packet for the sequential-model transition."""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ROOT = ROOT / "output/model/e5f_floor_open_population_transition"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structural-baseline",
        type=Path,
        default=DEFAULT_ROOT / "between_steady_states_static_elastic_long",
    )
    parser.add_argument(
        "--structural-policy",
        type=Path,
        default=DEFAULT_ROOT / "between_steady_states_static_elastic_ltv95_long",
    )
    parser.add_argument(
        "--structural-fixed-stock",
        type=Path,
        default=DEFAULT_ROOT / "between_steady_states_fixed_stock_endpoint_audit",
    )
    parser.add_argument(
        "--historical-retained-supply",
        type=Path,
        default=DEFAULT_ROOT / "historical_accounting_retained_supply",
    )
    parser.add_argument(
        "--historical-fixed-stock",
        type=Path,
        default=DEFAULT_ROOT / "historical_accounting_fixed_stock",
    )
    parser.add_argument(
        "--historical-fitted-supply",
        type=Path,
        default=DEFAULT_ROOT / "historical_accounting_fitted_supply",
    )
    parser.add_argument(
        "--historical-no-decline",
        type=Path,
        default=DEFAULT_ROOT / "historical_accounting_fitted_supply_no_decline",
    )
    parser.add_argument(
        "--feedback-audit",
        type=Path,
        default=ROOT / "output/model/sequential_fertility_price_feedback_audit",
    )
    parser.add_argument(
        "--outside-share-sensitivity",
        type=Path,
        default=DEFAULT_ROOT / "outside_share_sensitivity",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_ROOT / "decision_packet",
    )
    return parser.parse_args()


def load_run(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = path.resolve()
    transition = path / "transition_path.csv"
    summary = path / "summary.json"
    if not transition.exists() or not summary.exists():
        raise FileNotFoundError(f"Incomplete transition packet: {path}")
    return pd.read_csv(transition), json.loads(summary.read_text())


def load_historical(path: Path) -> pd.DataFrame:
    file = path.resolve() / "historical_comparison.csv"
    if not file.exists():
        raise FileNotFoundError(f"Missing historical comparison: {file}")
    return pd.read_csv(file)


def load_summary(path: Path) -> dict[str, Any]:
    file = path.resolve() / "summary.json"
    if not file.exists():
        raise FileNotFoundError(f"Missing run summary: {file}")
    return json.loads(file.read_text())


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def interpolate_row(frame: pd.DataFrame, period: int) -> pd.Series:
    selected = frame.loc[frame["period"] == period]
    if len(selected) != 1:
        raise ValueError(f"Expected one structural row for period {period}")
    return selected.iloc[0]


def make_figure(
    *,
    structural: pd.DataFrame,
    policy: pd.DataFrame,
    structural_summary: dict[str, Any],
    historical_retained: pd.DataFrame,
    historical_fixed: pd.DataFrame,
    historical_fitted: pd.DataFrame,
    historical_no_decline: pd.DataFrame,
    retained_supply_elasticity: float,
    fitted_supply_elasticity: float,
    feedback: pd.DataFrame,
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        "blue": "#24557a",
        "red": "#a34a3a",
        "green": "#3f7254",
        "gold": "#a77b24",
        "gray": "#555555",
    }
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.4), constrained_layout=True)

    observed = historical_fitted.sort_values("calendar_year")
    axes[0, 0].plot(
        observed.calendar_year,
        observed.observed_population_index,
        color="black",
        marker="o",
        lw=1.8,
        label="Census households",
    )
    axes[0, 0].plot(
        observed.calendar_year,
        observed.model_population_index,
        color=colors["blue"],
        marker="s",
        lw=1.8,
        label="Matched by construction",
    )
    axes[0, 0].set(
        title="A. Historical household stock (accounting target)",
        ylabel="Index (2007 = 1)",
    )
    axes[0, 0].legend(frameon=False, fontsize=8)

    axes[0, 1].plot(
        observed.calendar_year,
        observed.observed_real_house_price_index,
        color="black",
        marker="o",
        lw=2.0,
        label="Real house price",
    )
    for frame, label, color, marker in (
        (historical_fixed, "Fixed stock", colors["red"], "^"),
        (
            historical_retained,
            f"Supply elasticity {retained_supply_elasticity:g}",
            colors["green"],
            "s",
        ),
        (
            historical_fitted,
            f"Supply elasticity {fitted_supply_elasticity:g}",
            colors["blue"],
            "D",
        ),
    ):
        frame = frame.sort_values("calendar_year")
        axes[0, 1].plot(
            frame.calendar_year,
            frame.model_asset_price_index,
            color=color,
            marker=marker,
            lw=1.5,
            label=label,
        )
    no_decline = historical_no_decline.sort_values("calendar_year")
    axes[0, 1].plot(
        no_decline.calendar_year,
        no_decline.model_asset_price_index,
        color=colors["blue"],
        lw=1.3,
        ls="--",
        label=f"{fitted_supply_elasticity:g}, no fertility trend",
    )
    axes[0, 1].set(title="B. Historical housing cost", ylabel="Index (2007 = 1)")
    axes[0, 1].legend(frameon=False, fontsize=7)

    endpoint = structural_summary["stationary_open_endpoint"]
    years = structural.years_from_start.to_numpy()
    axes[0, 2].plot(
        years,
        structural.population_index,
        color=colors["blue"],
        lw=2.0,
        label="Transition",
    )
    axes[0, 2].axhline(
        float(endpoint["stationary_population_scale"]),
        color=colors["gray"],
        ls="--",
        lw=1.2,
        label="New open steady state",
    )
    axes[0, 2].set(
        title="C. Structural population transition",
        xlabel="Years after preference shift",
        ylabel="Adult-household mass",
    )
    axes[0, 2].legend(frameon=False, fontsize=8)

    axes[1, 0].plot(
        years,
        structural.asset_price_index,
        color=colors["red"],
        lw=2.0,
        label="Transition",
    )
    axes[1, 0].axhline(
        float(endpoint["price_ratio"]),
        color=colors["gray"],
        ls="--",
        lw=1.2,
        label="New open steady state",
    )
    axes[1, 0].set(
        title="D. Structural housing-price transition",
        xlabel="Years after preference shift",
        ylabel="House-price index",
    )
    axes[1, 0].legend(frameon=False, fontsize=8)

    joined = structural[
        [
            "period",
            "years_from_start",
            "adult_population",
            "owner_rate",
            "dependent_child_owner_rate",
        ]
    ].merge(
        policy[
            [
                "period",
                "adult_population",
                "owner_rate",
                "dependent_child_owner_rate",
            ]
        ],
        on="period",
        suffixes=("_baseline", "_policy"),
        validate="one_to_one",
    )
    joined["owner_effect_pp"] = 100.0 * (
        joined.owner_rate_policy - joined.owner_rate_baseline
    )
    joined["dependent_owner_effect_pp"] = 100.0 * (
        joined.dependent_child_owner_rate_policy
        - joined.dependent_child_owner_rate_baseline
    )
    joined["population_effect_bp"] = 10_000.0 * (
        joined.adult_population_policy / joined.adult_population_baseline - 1.0
    )
    policy_axis = axes[1, 1]
    population_axis = policy_axis.twinx()
    policy_axis.plot(
        joined.years_from_start,
        joined.dependent_owner_effect_pp,
        color=colors["blue"],
        lw=2.0,
        label="Dependent-child ownership (pp)",
    )
    policy_axis.plot(
        joined.years_from_start,
        joined.owner_effect_pp,
        color=colors["green"],
        lw=1.8,
        label="Overall ownership (pp)",
    )
    population_axis.plot(
        joined.years_from_start,
        joined.population_effect_bp,
        color=colors["gold"],
        lw=1.4,
        ls="--",
        label="Population (basis points)",
    )
    policy_axis.set(
        title="E. Dependent-child LTV95 policy",
        xlabel="Years after preference shift",
        ylabel="Ownership effect (percentage points)",
    )
    handles_left, labels_left = policy_axis.get_legend_handles_labels()
    handles_right, labels_right = population_axis.get_legend_handles_labels()
    policy_axis.legend(
        handles_left + handles_right,
        labels_left + labels_right,
        frameon=False,
        fontsize=7,
    )

    highlighted_feedback = {
        1.0: colors["blue"],
        0.25: "#7a5aa6",
        0.10: "#8a5a44",
        0.05: "#d56db3",
        0.02: colors["gray"],
    }
    for factor in sorted(feedback.kappa_factor.unique(), reverse=True):
        rows = feedback.loc[feedback.kappa_factor == factor].sort_values("price_ratio")
        highlight = float(factor) in highlighted_feedback
        axes[1, 2].plot(
            rows.price_ratio,
            rows.reproduction_ratio_B_over_E,
            color=(
                highlighted_feedback[float(factor)]
                if highlight
                else "#c6c6c6"
            ),
            marker="o" if highlight else None,
            ms=3 if highlight else 0,
            lw=1.6 if highlight else 0.8,
            alpha=1.0 if highlight else 0.6,
            label=f"fertility scale x {factor:g}" if highlight else None,
        )
    axes[1, 2].axhline(1.0, color="black", ls="--", lw=1.0)
    axes[1, 2].set(
        title="F. Closed-population existence test",
        xlabel="House-price ratio",
        ylabel="Locally born entrants / required entrants",
        xscale="log",
    )
    axes[1, 2].legend(frameon=False, fontsize=6.5)

    for axis in axes.flat:
        axis.grid(alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "decision_figure.png", dpi=240)
    fig.savefig(outdir / "decision_figure.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    structural, structural_summary = load_run(args.structural_baseline)
    policy, policy_summary = load_run(args.structural_policy)
    _, fixed_summary = load_run(args.structural_fixed_stock)
    historical_retained = load_historical(args.historical_retained_supply)
    historical_fixed = load_historical(args.historical_fixed_stock)
    historical_fitted = load_historical(args.historical_fitted_supply)
    historical_no_decline = load_historical(args.historical_no_decline)
    retained_historical_summary = load_summary(args.historical_retained_supply)
    fitted_historical_summary = load_summary(args.historical_fitted_supply)
    retained_supply_elasticity = float(
        retained_historical_summary["supply_normalization"][
            "transition_housing_supply_elasticity"
        ]
    )
    fitted_supply_elasticity = float(
        fitted_historical_summary["supply_normalization"][
            "transition_housing_supply_elasticity"
        ]
    )
    feedback = pd.read_csv(args.feedback_audit.resolve() / "price_schedules.csv")
    feedback_summary = pd.read_csv(
        args.feedback_audit.resolve() / "sensitivity_summary.csv"
    )
    outside_share_runs: list[tuple[Path, dict[str, Any]]] = []
    for summary_path in sorted(
        args.outside_share_sensitivity.resolve().glob("share_*/summary.json")
    ):
        summary = json.loads(summary_path.read_text())
        endpoint_sensitivity = summary.get("stationary_open_endpoint") or {}
        if endpoint_sensitivity.get("status") != "complete":
            raise RuntimeError(
                f"Incomplete outside-entry sensitivity: {summary_path.parent}"
            )
        outside_share_runs.append((summary_path.parent, summary))
    if not outside_share_runs:
        raise FileNotFoundError(
            "No outside-entry sensitivity packets found under "
            f"{args.outside_share_sensitivity.resolve()}"
        )

    endpoint = structural_summary.get("stationary_open_endpoint") or {}
    if endpoint.get("status") != "complete":
        raise RuntimeError("Structural baseline has no completed open endpoint")
    policy_endpoint = policy_summary.get("stationary_open_endpoint") or {}
    if policy_endpoint.get("status") != "complete":
        raise RuntimeError("Structural policy has no completed open endpoint")
    fixed_endpoint = fixed_summary.get("stationary_open_endpoint") or {}
    if not str(fixed_endpoint.get("status", "")).startswith("no_positive_price_root"):
        raise RuntimeError("Fixed-stock packet does not contain the expected no-root audit")
    feedback_summary["closed_root_found"] = (
        feedback_summary.closed_root_found.astype(str).str.lower() == "true"
    )
    moderate_feedback = feedback_summary.loc[
        feedback_summary.kappa_factor >= 0.25
    ]
    if bool(moderate_feedback.closed_root_found.any()):
        raise RuntimeError(
            "A closed root appeared within the prespecified fourfold slope audit; "
            "rewrite the packet"
        )
    extreme_roots = feedback_summary.loc[feedback_summary.closed_root_found]

    target_2023 = historical_fitted.loc[
        historical_fitted.calendar_year == 2023
    ].iloc[0]
    historical_rows: list[dict[str, Any]] = []
    for label, frame, elasticity, status in (
        ("fixed_stock", historical_fixed, 0.0, "benchmark"),
        (
            "retained_supply",
            historical_retained,
            retained_supply_elasticity,
            "retained_model_value",
        ),
        (
            "endpoint_fitted_supply",
            historical_fitted,
            fitted_supply_elasticity,
            "one_moment_diagnostic",
        ),
    ):
        row = frame.loc[frame.calendar_year == 2023].iloc[0]
        historical_rows.append(
            {
                "case": label,
                "housing_supply_elasticity": elasticity,
                "parameter_status": status,
                "diagnostic_search_lower": 0.0,
                "diagnostic_search_upper": retained_supply_elasticity,
                "near_search_bound": bool(
                    status == "one_moment_diagnostic"
                    and min(
                        elasticity,
                        retained_supply_elasticity - elasticity,
                    )
                    <= 0.05 * retained_supply_elasticity
                ),
                "model_2023_house_price_index": row.model_asset_price_index,
                "observed_2023_real_house_price_index": row.observed_real_house_price_index,
                "gap": row.model_asset_price_index - row.observed_real_house_price_index,
            }
        )
    no_decline_supply_row = historical_no_decline.loc[
        historical_no_decline.calendar_year == 2023
    ].iloc[0]
    historical_rows.append(
        {
            "case": "endpoint_fitted_supply_no_fertility_decline",
            "housing_supply_elasticity": fitted_supply_elasticity,
            "parameter_status": "counterfactual",
            "diagnostic_search_lower": 0.0,
            "diagnostic_search_upper": retained_supply_elasticity,
            "near_search_bound": False,
            "model_2023_house_price_index": (
                no_decline_supply_row.model_asset_price_index
            ),
            "observed_2023_real_house_price_index": (
                no_decline_supply_row.observed_real_house_price_index
            ),
            "gap": (
                no_decline_supply_row.model_asset_price_index
                - no_decline_supply_row.observed_real_house_price_index
            ),
        }
    )
    pd.DataFrame(historical_rows).to_csv(
        outdir / "historical_supply_comparison.csv", index=False
    )

    policy_rows: list[dict[str, Any]] = []
    common_periods = [4, 9, 24, int(structural.period.max())]
    for period in dict.fromkeys(common_periods):
        baseline_row = interpolate_row(structural, period)
        policy_row = interpolate_row(policy, period)
        policy_rows.append(
            {
                "period": period,
                "years_after_shock": baseline_row.years_from_start,
                "population_effect_percent": 100.0
                * (policy_row.adult_population / baseline_row.adult_population - 1.0),
                "house_price_effect_percent": 100.0
                * (policy_row.asset_price / baseline_row.asset_price - 1.0),
                "births_per_adult_effect_percent": 100.0
                * (policy_row.births_per_adult / baseline_row.births_per_adult - 1.0),
                "overall_ownership_effect_pp": 100.0
                * (policy_row.owner_rate - baseline_row.owner_rate),
                "dependent_child_ownership_effect_pp": 100.0
                * (
                    policy_row.dependent_child_owner_rate
                    - baseline_row.dependent_child_owner_rate
                ),
            }
        )
    pd.DataFrame(policy_rows).to_csv(outdir / "policy_effects.csv", index=False)
    pd.DataFrame(
        [
            {
                "case": "baseline",
                "adult_household_scale": endpoint["stationary_population_scale"],
                "house_price_ratio": endpoint["price_ratio"],
                "completed_fertility_statistic": endpoint["tfr_topcoded"],
                "owner_rate": endpoint["own_rate"],
            },
            {
                "case": "dependent_child_ltv95",
                "adult_household_scale": policy_endpoint[
                    "stationary_population_scale"
                ],
                "house_price_ratio": policy_endpoint["price_ratio"],
                "completed_fertility_statistic": policy_endpoint["tfr_topcoded"],
                "owner_rate": policy_endpoint["own_rate"],
            },
        ]
    ).to_csv(outdir / "stationary_comparison.csv", index=False)

    outside_share_rows = [
        {
            "outside_origin_entry_share_at_old_steady_state": float(
                summary["renewal_closure"][
                    "outside_origin_entry_share_at_old_steady_state"
                ]
            ),
            "local_retention_rate": float(
                summary["renewal_closure"]["retention_rho"]
            ),
            "new_steady_state_adult_household_scale": float(
                summary["stationary_open_endpoint"][
                    "stationary_population_scale"
                ]
            ),
            "new_steady_state_house_price_ratio": float(
                summary["stationary_open_endpoint"]["price_ratio"]
            ),
            "new_steady_state_completed_fertility_statistic": float(
                summary["stationary_open_endpoint"]["tfr_topcoded"]
            ),
            "parameter_status": "diagnostic_sensitivity",
        }
        for _, summary in outside_share_runs
    ]
    outside_share_rows.append(
        {
            "outside_origin_entry_share_at_old_steady_state": float(
                structural_summary["renewal_closure"][
                    "outside_origin_entry_share_at_old_steady_state"
                ]
            ),
            "local_retention_rate": float(
                structural_summary["renewal_closure"]["retention_rho"]
            ),
            "new_steady_state_adult_household_scale": float(
                endpoint["stationary_population_scale"]
            ),
            "new_steady_state_house_price_ratio": float(endpoint["price_ratio"]),
            "new_steady_state_completed_fertility_statistic": float(
                endpoint["tfr_topcoded"]
            ),
            "parameter_status": "provisional_baseline",
        }
    )
    outside_share_table = pd.DataFrame(outside_share_rows).sort_values(
        "outside_origin_entry_share_at_old_steady_state"
    )
    outside_share_table.to_csv(outdir / "outside_entry_sensitivity.csv", index=False)

    make_figure(
        structural=structural,
        policy=policy,
        structural_summary=structural_summary,
        historical_retained=historical_retained,
        historical_fixed=historical_fixed,
        historical_fitted=historical_fitted,
        historical_no_decline=historical_no_decline,
        retained_supply_elasticity=retained_supply_elasticity,
        fitted_supply_elasticity=fitted_supply_elasticity,
        feedback=feedback,
        outdir=outdir,
    )

    manifest_paths = {
        "structural_baseline": args.structural_baseline.resolve(),
        "structural_policy": args.structural_policy.resolve(),
        "structural_fixed_stock": args.structural_fixed_stock.resolve(),
        "historical_retained_supply": args.historical_retained_supply.resolve(),
        "historical_fixed_stock": args.historical_fixed_stock.resolve(),
        "historical_fitted_supply": args.historical_fitted_supply.resolve(),
        "historical_no_decline": args.historical_no_decline.resolve(),
        "feedback_audit": args.feedback_audit.resolve(),
    }
    manifest_paths.update(
        {
            f"outside_share_{float(summary['renewal_closure']['outside_origin_entry_share_at_old_steady_state']):g}": path
            for path, summary in outside_share_runs
        }
    )
    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    write_json(
        outdir / "manifest.json",
        {
            "status": "complete",
            "exact_command": command,
            "code": {
                str(path.relative_to(ROOT)): sha256(path)
                for path in (
                    ROOT / "code/model/tools/run_e5f_open_population_transition.py",
                    ROOT / "code/model/tools/audit_e5f_fertility_price_feedback.py",
                    ROOT / "code/model/tools/build_e5f_transition_decision_packet.py",
                    ROOT
                    / "code/data/Spatial_aggregate_withmicrodata/"
                    "build_national_householder_age_path.py",
                )
            },
            "inputs": {
                name: {
                    "path": str(path),
                    "summary_sha256": sha256(path / "summary.json"),
                }
                for name, path in manifest_paths.items()
            },
            "structural_open_endpoint": endpoint,
            "structural_policy_open_endpoint": policy_endpoint,
            "structural_fixed_stock_endpoint": fixed_endpoint,
            "historical_2023": target_2023.to_dict(),
            "feedback_summary": feedback_summary.to_dict(orient="records"),
            "outside_entry_sensitivity": outside_share_table.to_dict(
                orient="records"
            ),
            "policy_case": policy_summary.get("policy_case"),
        },
    )

    fitted_gap = float(
        target_2023.model_asset_price_index
        - target_2023.observed_real_house_price_index
    )
    no_decline_2023 = historical_no_decline.loc[
        historical_no_decline.calendar_year == 2023
    ].iloc[0]
    preference_price_effect = 100.0 * (
        target_2023.model_asset_price_index
        / no_decline_2023.model_asset_price_index
        - 1.0
    )
    policy_endpoint_population_effect = 100.0 * (
        policy_endpoint["stationary_population_scale"]
        / endpoint["stationary_population_scale"]
        - 1.0
    )
    policy_endpoint_ownership_effect = 100.0 * (
        policy_endpoint["own_rate"] - endpoint["own_rate"]
    )
    maximum_feedback = float(
        moderate_feedback.maximum_B_over_E_on_price_grid.max()
    )
    minimum_outside_scale = float(
        outside_share_table.new_steady_state_adult_household_scale.min()
    )
    maximum_outside_scale = float(
        outside_share_table.new_steady_state_adult_household_scale.max()
    )
    minimum_outside_share = float(
        outside_share_table.outside_origin_entry_share_at_old_steady_state.min()
    )
    maximum_outside_share = float(
        outside_share_table.outside_origin_entry_share_at_old_steady_state.max()
    )
    if extreme_roots.empty:
        extreme_feedback_sentence = (
            "Even the deliberately extreme sensitivities down to a 0.02 logit-scale "
            "factor do not restore a closed root."
        )
    else:
        least_extreme_root = extreme_roots.sort_values(
            "kappa_factor", ascending=False
        ).iloc[0]
        extreme_feedback_sentence = (
            "A closed root first reappears only in the deliberately extreme grid "
            f"at a logit-scale factor of {least_extreme_root.kappa_factor:g} "
            f"(a {1.0 / least_extreme_root.kappa_factor:g}-fold smaller noise scale); "
            "this is an identification threshold, not a calibrated alternative."
        )
    readme = f"""# Sequential transition: decision packet

## Bottom line

The minimal quantitative extension is viable as an **open-population transition**. The existing sequential household problem is unchanged: the wrapper adds calendar-time cohort accounting, a delayed birth-to-entry queue, and housing-market clearing at each date. Starting from the model's old open steady state, the fertility-preference shift moves toward a lower open steady state with adult-household mass `{endpoint['stationary_population_scale']:.3f}` and house price `{endpoint['price_ratio']:.3f}` relative to the old equilibrium. This magnitude is not pinned down: moving the outside-origin entrant share from `{minimum_outside_share:.2f}` to `{maximum_outside_share:.2f}` moves the endpoint household scale from `{minimum_outside_scale:.3f}` to `{maximum_outside_scale:.3f}`. The sensitivity is reported in `outside_entry_sensitivity.csv`; the baseline share remains provisional.

The stronger closed-population claim does **not** currently survive the quantitative model. From the estimated fertility logits through a fourfold increase in responsiveness, the largest post-shock locally born renewal ratio is `{maximum_feedback:.3f}`, below replacement at one. {extreme_feedback_sentence} The simple theory should therefore state its existence condition clearly; the calibrated model should use the open closure unless a causal housing-cost-to-fertility moment supports a substantially different response.

## What the six panels establish

1. The 2007--2023 historical bridge exactly matches the published HH-3 household total after mapping it to the model's head ages 18--85 with the national ACS, as well as the ACS four-year householder-age shares. This makes 2007 an observed transition state, not literally a demographic steady state.
2. Fixed housing stock overpredicts the 2023 price increase; the retained long-run supply elasticity underpredicts it. A diagnostic elasticity of `{fitted_supply_elasticity:g}` matches the 2023 endpoint with a gap of `{fitted_gap:.6f}`, but none of the supply cases explains the 2007--2011 price crash. Holding the imposed household path fixed, the fertility-preference decline changes the 2023 model price by only `{preference_price_effect:.3f}` percent relative to no decline.
3. In the clean structural experiment, population and housing prices move from the old open steady state toward the lower open steady state.
4. A permanently fixed inherited housing stock has no usable positive-price terminal root on the audited range. Long-run convergence needs depreciation, vacancy, demolition, or elastic stock adjustment.
5. LTV95 for dependent-child households raises steady-state ownership by `{policy_endpoint_ownership_effect:.3f}` percentage points but steady-state household mass by only `{policy_endpoint_population_effect:.3f}` percent. The tenure/intergenerational mechanism therefore remains central; this policy is not a demographic cure. This is a pure credit-access diagnostic without default pricing or fiscal accounting.
6. Lower housing prices do not restore closed-population replacement under the current sequential calibration or the prespecified fourfold responsiveness check. The extreme rows locate how far the fertility slope would have to move, if a root appears at all.

## Interpretation for the paper

Use the clean steady-state transition as the structural counterfactual and the 2007--2023 age bridge as a separate historical plausibility check. Do not call the historically initialized 2007 economy a steady state. Do not claim a closed lower steady state for the quantitative model. The next model change, if required, is a small housing-stock law; a full redesign of household optimization is not needed.

All household decisions use temporary-equilibrium price expectations, so this packet is not a perfect-foresight welfare exercise. The fitted supply elasticity and the fertility-logit rescalings are explicitly diagnostic, not a recalibration.

Exact build command:

```bash
{command}
```
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print(f"DECISION_PACKET_COMPLETE outdir={outdir}", flush=True)


if __name__ == "__main__":
    main()
