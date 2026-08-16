#!/usr/bin/env python3
"""Build one compact decision packet for the sequential-model transition."""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import shutil
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
        "--option-grid",
        type=Path,
        default=ROOT / "output/model/e5f_option_grid",
    )
    parser.add_argument(
        "--slope-fit-audit",
        type=Path,
        default=ROOT / "output/model/e5f_slope_fit_tradeoff",
    )
    parser.add_argument(
        "--calibration-preflight",
        type=Path,
        default=(
            ROOT
            / "output/model/e5f_corrected_target_preflight_20260816/chain_1_v2"
        ),
    )
    parser.add_argument(
        "--preliminary-canonical-refit",
        type=Path,
        default=(
            ROOT / "output/model/e5f_preliminary_refits_20260816/canonical/chain_1"
        ),
    )
    parser.add_argument(
        "--closed-slope-preliminary-refit",
        type=Path,
        default=(
            ROOT
            / "output/model/e5f_preliminary_refits_20260816/"
            "kappa_both_factor_0p05_100eval/chain_1"
        ),
    )
    parser.add_argument(
        "--closed-slope-feedback-audit",
        type=Path,
        default=(
            ROOT / "output/model/e5f_closed_slope_refit_feedback_audit_20260816"
        ),
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


def calibration_parameter_table(
    winner: dict[str, Any],
    metadata: dict[str, Any],
    *,
    run_label: str,
) -> pd.DataFrame:
    """Report every free and fixed parameter under one transparent contract."""
    theta = {name: float(value) for name, value in winner["theta"].items()}
    rows: list[dict[str, Any]] = []
    active_names: set[str] = set()
    for specification in metadata["active_domain"]:
        name = str(specification["name"])
        active_names.add(name)
        theta_name = "beta" if name == "beta_annual" else name
        estimate = float(theta[theta_name])
        if name == "beta_annual":
            estimate **= 0.25
        lower = float(specification["lower"])
        upper = float(specification["upper"])
        distance = min(estimate - lower, upper - estimate)
        rows.append(
            {
                "run": run_label,
                "parameter": name,
                "estimate": estimate,
                "lower_bound": lower,
                "upper_bound": upper,
                "transform": specification["transform"],
                "distance_to_nearest_bound": distance,
                "distance_as_share_of_range": distance / (upper - lower),
                "near_bound_2pct": distance / (upper - lower) <= 0.02,
                "parameter_status": "estimated",
                "external_restriction": False,
                "external_source": "",
            }
        )
    external_metadata = metadata.get("external_parameter_metadata") or {}
    probe_names = {
        name
        for name in ("kappa_fert", "kappa_fert_continuation")
        if f"probe_fixed_{name}" in metadata
    }
    for name, value in metadata["fixed_parameters"].items():
        if name in active_names:
            continue
        external = external_metadata.get(name) or {}
        is_probe = name in probe_names
        rows.append(
            {
                "run": run_label,
                "parameter": name,
                "estimate": float(value),
                "lower_bound": "",
                "upper_bound": "",
                "transform": "",
                "distance_to_nearest_bound": "",
                "distance_as_share_of_range": "",
                "near_bound_2pct": "",
                "parameter_status": (
                    "fixed_for_closed_slope_probe" if is_probe else "externally_fixed"
                ),
                "external_restriction": not is_probe,
                "external_source": (
                    "diagnostic 0.05 slope factor"
                    if is_probe
                    else external.get("source", "profile convention")
                ),
            }
        )
    return pd.DataFrame(rows)


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
    }
    for factor in sorted(feedback.kappa_factor.unique(), reverse=True):
        if float(factor) < 0.05:
            continue
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
    write_json(outdir / "manifest.json", {"status": "running"})

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
    calibration_preflight = load_summary(args.calibration_preflight)
    calibration_winner = calibration_preflight.get("best_tight") or {}
    calibration_metadata = calibration_preflight.get("metadata") or {}
    if (
        calibration_preflight.get("status") != "complete"
        or not calibration_winner.get("strict_converged", False)
        or not calibration_preflight.get("tight_repeat_check", {}).get(
            "both_strict", False
        )
    ):
        raise RuntimeError("The current-target calibration preflight is not strict")
    preliminary_refit = load_summary(args.preliminary_canonical_refit)
    preliminary_refit_winner = preliminary_refit.get("best_tight") or {}
    if (
        preliminary_refit.get("status") != "complete"
        or int(preliminary_refit.get("n_cases_completed", 0)) != 25
        or not preliminary_refit_winner.get("strict_converged", False)
        or not preliminary_refit.get("tight_repeat_check", {}).get(
            "both_strict", False
        )
    ):
        raise RuntimeError("The preliminary current-target refit is incomplete")
    closed_slope_refit = load_summary(args.closed_slope_preliminary_refit)
    closed_slope_winner = closed_slope_refit.get("best_tight") or {}
    closed_slope_metadata = closed_slope_refit.get("metadata") or {}
    if (
        closed_slope_refit.get("status") != "complete"
        or int(closed_slope_refit.get("n_cases_completed", 0)) != 100
        or not closed_slope_winner.get("strict_converged", False)
        or not closed_slope_refit.get("tight_repeat_check", {}).get(
            "both_strict", False
        )
    ):
        raise RuntimeError("The preliminary closed-slope refit is incomplete")
    for name, expected in (
        ("probe_fixed_kappa_fert", 0.1493420161260494),
        ("probe_fixed_kappa_fert_continuation", 0.08566620151181108),
    ):
        if not np.isclose(
            float(closed_slope_metadata.get(name, np.nan)),
            expected,
            rtol=0.0,
            atol=1e-14,
        ):
            raise RuntimeError(f"Closed-slope probe contract mismatch for {name}")
    closed_slope_feedback = load_summary(args.closed_slope_feedback_audit)
    closed_slope_feedback_rows = closed_slope_feedback.get("rows") or []
    if (
        closed_slope_feedback.get("status") != "complete"
        or [float(value) for value in closed_slope_feedback.get("kappa_factors", [])]
        != [1.0]
        or len(closed_slope_feedback_rows) != 1
    ):
        raise RuntimeError("The refitted closed-slope feedback audit is incomplete")
    if closed_slope_feedback.get("source_sha256") != sha256(
        args.closed_slope_preliminary_refit.resolve() / "summary.json"
    ):
        raise RuntimeError("The refitted-slope feedback audit source hash is stale")
    closed_slope_feedback_row = closed_slope_feedback_rows[0]
    slope_fit_path = args.slope_fit_audit.resolve()
    slope_fit_summary = json.loads((slope_fit_path / "summary.json").read_text())
    slope_fit = pd.read_csv(slope_fit_path / "slope_fit_summary.csv")
    slope_target_fit = pd.read_csv(slope_fit_path / "target_fit_by_factor.csv")
    if slope_fit_summary.get("status") != "complete":
        raise RuntimeError("The slope-fit trade-off audit is incomplete")

    option_grid_runs: dict[tuple[float, float, str], tuple[Path, dict[str, Any]]] = {}
    for summary_path in sorted(args.option_grid.resolve().glob("supply_*/share_*/*/summary.json")):
        summary = json.loads(summary_path.read_text())
        endpoint_grid = summary.get("stationary_open_endpoint") or {}
        gates = summary.get("operator_gates") or {}
        scenario_gates = summary.get("scenario_summary") or {}
        if (
            summary.get("status") != "complete"
            or endpoint_grid.get("status") != "complete"
            or not gates.get("passed", False)
            or float(scenario_gates.get("max_market_residual", float("inf")))
            >= 5e-5
            or float(
                scenario_gates.get("max_abs_mass_accounting_residual", float("inf"))
            )
            >= 1e-12
            or int(scenario_gates.get("max_nonfinite_distribution_count", 1)) != 0
        ):
            raise RuntimeError(f"Option-grid gate failed: {summary_path.parent}")
        elasticity = float(
            summary["supply_normalization"][
                "transition_housing_supply_elasticity"
            ]
        )
        outside_share = float(
            summary["renewal_closure"][
                "outside_origin_entry_share_at_old_steady_state"
            ]
        )
        policy_name = str(summary["policy_case"])
        key = (elasticity, outside_share, policy_name)
        if key in option_grid_runs:
            raise RuntimeError(f"Duplicate option-grid cell: {key}")
        option_grid_runs[key] = (summary_path.parent, summary)
    expected_option_cells = {
        (elasticity, outside_share, policy_name)
        for elasticity in (0.3675, 1.75)
        for outside_share in (0.10, 0.169, 0.40, 0.60)
        for policy_name in ("none", "dependent-child-ltv95")
    }
    if set(option_grid_runs) != expected_option_cells:
        raise RuntimeError(
            "Option grid is incomplete: "
            f"missing={sorted(expected_option_cells - set(option_grid_runs))}, "
            f"extra={sorted(set(option_grid_runs) - expected_option_cells)}"
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

    calibration_target_fit = pd.DataFrame(calibration_winner["target_fit"])
    if len(calibration_target_fit) != int(calibration_metadata["target_count"]):
        raise RuntimeError("Current-target fit table is incomplete")
    calibration_target_fit.to_csv(
        outdir / "provisional_calibration_target_fit.csv", index=False
    )
    calibration_parameter_table(
        calibration_winner,
        calibration_metadata,
        run_label="saved_provisional_fit_under_current_targets",
    ).to_csv(
        outdir / "provisional_calibration_parameters.csv", index=False
    )
    closed_slope_target_fit = pd.DataFrame(closed_slope_winner["target_fit"])
    if len(closed_slope_target_fit) != int(closed_slope_metadata["target_count"]):
        raise RuntimeError("Closed-slope target-fit table is incomplete")
    closed_slope_target_fit.to_csv(
        outdir / "closed_slope_preliminary_target_fit.csv", index=False
    )
    pd.DataFrame([closed_slope_feedback_row]).to_csv(
        outdir / "closed_slope_preliminary_feedback.csv", index=False
    )
    calibration_parameter_table(
        closed_slope_winner,
        closed_slope_metadata,
        run_label="closed_slope_100_evaluation_preliminary_refit",
    ).to_csv(
        outdir / "closed_slope_preliminary_parameters.csv", index=False
    )
    slope_fit.to_csv(outdir / "fertility_slope_fit_tradeoff.csv", index=False)
    slope_target_fit.to_csv(
        outdir / "fertility_slope_target_fit_full.csv", index=False
    )

    option_rows: list[dict[str, Any]] = []
    for elasticity in (0.3675, 1.75):
        for outside_share in (0.10, 0.169, 0.40, 0.60):
            _, baseline_summary = option_grid_runs[
                (elasticity, outside_share, "none")
            ]
            _, policy_summary_grid = option_grid_runs[
                (elasticity, outside_share, "dependent-child-ltv95")
            ]
            baseline_endpoint = baseline_summary["stationary_open_endpoint"]
            policy_endpoint_grid = policy_summary_grid["stationary_open_endpoint"]
            option_rows.append(
                {
                    "housing_supply_elasticity": elasticity,
                    "housing_supply_status": (
                        "one_moment_historical_diagnostic"
                        if elasticity == 0.3675
                        else "retained_model_value"
                    ),
                    "outside_origin_entry_share": outside_share,
                    "outside_share_status": (
                        "provisional_baseline"
                        if outside_share == 0.169
                        else "diagnostic_sensitivity"
                    ),
                    "baseline_adult_household_scale": float(
                        baseline_endpoint["stationary_population_scale"]
                    ),
                    "baseline_house_price_ratio": float(
                        baseline_endpoint["price_ratio"]
                    ),
                    "baseline_completed_fertility_statistic": float(
                        baseline_endpoint["tfr_topcoded"]
                    ),
                    "baseline_owner_rate": float(baseline_endpoint["own_rate"]),
                    "ltv95_adult_household_scale": float(
                        policy_endpoint_grid["stationary_population_scale"]
                    ),
                    "ltv95_house_price_ratio": float(
                        policy_endpoint_grid["price_ratio"]
                    ),
                    "ltv95_completed_fertility_statistic": float(
                        policy_endpoint_grid["tfr_topcoded"]
                    ),
                    "ltv95_owner_rate": float(policy_endpoint_grid["own_rate"]),
                    "ltv95_population_effect_percent": 100.0
                    * (
                        float(policy_endpoint_grid["stationary_population_scale"])
                        / float(baseline_endpoint["stationary_population_scale"])
                        - 1.0
                    ),
                    "ltv95_house_price_effect_percent": 100.0
                    * (
                        float(policy_endpoint_grid["price_ratio"])
                        / float(baseline_endpoint["price_ratio"])
                        - 1.0
                    ),
                    "ltv95_completed_fertility_effect_percent": 100.0
                    * (
                        float(policy_endpoint_grid["tfr_topcoded"])
                        / float(baseline_endpoint["tfr_topcoded"])
                        - 1.0
                    ),
                    "ltv95_ownership_effect_percentage_points": 100.0
                    * (
                        float(policy_endpoint_grid["own_rate"])
                        - float(baseline_endpoint["own_rate"])
                    ),
                }
            )
    option_grid_table = pd.DataFrame(option_rows)
    option_grid_table.to_csv(
        outdir / "closure_supply_policy_grid.csv", index=False
    )
    for suffix in ("png", "pdf"):
        shutil.copy2(
            slope_fit_path / f"slope_fit_tradeoff.{suffix}",
            outdir / f"fertility_slope_fit_tradeoff.{suffix}",
        )

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
        "slope_fit_audit": args.slope_fit_audit.resolve(),
        "calibration_preflight": args.calibration_preflight.resolve(),
        "preliminary_canonical_refit": args.preliminary_canonical_refit.resolve(),
        "closed_slope_preliminary_refit": (
            args.closed_slope_preliminary_refit.resolve()
        ),
        "closed_slope_feedback_audit": args.closed_slope_feedback_audit.resolve(),
    }
    manifest_paths.update(
        {
            f"outside_share_{float(summary['renewal_closure']['outside_origin_entry_share_at_old_steady_state']):g}": path
            for path, summary in outside_share_runs
        }
    )
    manifest_paths.update(
        {
            (
                f"option_grid_xi_{elasticity:g}_outside_{outside_share:g}_"
                f"{policy_name}"
            ): path
            for (elasticity, outside_share, policy_name), (
                path,
                _summary,
            ) in option_grid_runs.items()
        }
    )
    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    write_json(
        outdir / "manifest.json",
        {
            "status": "building",
            "exact_command": command,
            "code": {
                str(path.relative_to(ROOT)): sha256(path)
                for path in (
                    ROOT / "code/model/tools/run_e5f_open_population_transition.py",
                    ROOT / "code/model/tools/audit_e5f_fertility_price_feedback.py",
                    ROOT / "code/model/tools/audit_e5f_slope_fit_tradeoff.py",
                    ROOT / "code/model/tools/build_e5f_transition_decision_packet.py",
                    ROOT
                    / "code/model/intergen_eqscale_seq_optimized/run_e1_chain.py",
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
            "current_target_calibration": {
                "loss": float(calibration_winner["rank_loss"]),
                "market_residual": float(calibration_winner["market_residual"]),
                "target_set": calibration_metadata["target_set"],
                "target_fingerprint": calibration_metadata[
                    "search_target_fingerprint"
                ],
                "target_count": int(calibration_metadata["target_count"]),
                "free_parameter_count": int(
                    calibration_metadata["free_parameter_count"]
                ),
                "preliminary_refit_cases": int(
                    preliminary_refit["n_cases_completed"]
                ),
                "preliminary_refit_loss": float(
                    preliminary_refit_winner["rank_loss"]
                ),
                "closed_slope_preliminary_refit_cases": int(
                    closed_slope_refit["n_cases_completed"]
                ),
                "closed_slope_preliminary_refit_loss": float(
                    closed_slope_winner["rank_loss"]
                ),
            },
            "slope_fit_tradeoff": slope_fit.to_dict(orient="records"),
            "closed_slope_refit_feedback": closed_slope_feedback,
            "closure_supply_policy_grid": option_grid_table.to_dict(
                orient="records"
            ),
            "policy_case": policy_summary.get("policy_case"),
        },
    )

    no_decline_2023 = historical_no_decline.loc[
        historical_no_decline.calendar_year == 2023
    ].iloc[0]
    preference_price_effect = 100.0 * (
        target_2023.model_asset_price_index
        / no_decline_2023.model_asset_price_index
        - 1.0
    )
    maximum_feedback = float(
        moderate_feedback.maximum_B_over_E_on_price_grid.max()
    )
    old_steady_state = structural_summary["old_steady_state"]
    old_psi = float(structural_summary["old_psi_child"])
    new_psi = float(structural_summary["new_psi_child"])
    baseline_outside_share = float(
        structural_summary["renewal_closure"][
            "outside_origin_entry_share_at_old_steady_state"
        ]
    )
    endpoint_stock_decline_percent = 100.0 * (
        1.0
        - float(endpoint["stationary_housing_supply"])
        / float(old_steady_state["housing_supply"])
    )
    current_loss = float(calibration_winner["rank_loss"])
    canonical_refit_loss = float(preliminary_refit_winner["rank_loss"])
    closed_slope_refit_loss = float(closed_slope_winner["rank_loss"])
    slope_root_indicator = (
        slope_fit.closed_root_found.astype(str).str.strip().str.lower() == "true"
    )
    root_rows = slope_fit.loc[slope_root_indicator]
    if root_rows.empty:
        raise RuntimeError("Slope-fit audit found no closed root")
    first_root = root_rows.sort_values("kappa_factor", ascending=False).iloc[0]
    first_root_factor = float(first_root.kappa_factor)
    first_root_fixed_loss = float(first_root.canonical_loss)
    room_target_contribution = float(
        calibration_target_fit.loc[
            calibration_target_fit.moment == "housing_increment_0to1",
            "loss_contribution",
        ].iloc[0]
    )
    largest_current_fit_rows = calibration_target_fit.nlargest(3, "loss_contribution")
    fit_labels = {
        "childless_rate": "childlessness",
        "old_total_wealth_to_annual_income_p90_p50_7684": "old-age wealth dispersion",
        "share_first_births_age30plus": "late first births",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": (
            "the family-size room gradient"
        ),
        "own_family_gap": "the family ownership gap",
    }
    largest_current_fit_text = ", ".join(
        f"{fit_labels.get(str(row.moment), str(row.moment))} "
        f"({float(row.loss_contribution):.1f})"
        for row in largest_current_fit_rows.itertuples()
    )
    endpoint_scale_min = float(option_grid_table.baseline_adult_household_scale.min())
    endpoint_scale_max = float(option_grid_table.baseline_adult_household_scale.max())
    endpoint_price_min = float(option_grid_table.baseline_house_price_ratio.min())
    endpoint_price_max = float(option_grid_table.baseline_house_price_ratio.max())
    policy_population_min = float(
        option_grid_table.ltv95_population_effect_percent.min()
    )
    policy_population_max = float(
        option_grid_table.ltv95_population_effect_percent.max()
    )
    policy_ownership_min = float(
        option_grid_table.ltv95_ownership_effect_percentage_points.min()
    )
    policy_ownership_max = float(
        option_grid_table.ltv95_ownership_effect_percentage_points.max()
    )
    nonsaturated_policy_rows = option_grid_table.loc[
        option_grid_table.baseline_owner_rate < 0.95
    ]
    policy_ownership_nonsaturated_min = float(
        nonsaturated_policy_rows.ltv95_ownership_effect_percentage_points.min()
    )
    if bool(closed_slope_feedback_row["closed_root_found"]):
        closed_slope_root_sentence = (
            "The independently rerun renewal audit still finds a closed root, at "
            f"a house-price ratio of "
            f"{float(closed_slope_feedback_row['closed_root_price_ratio']):.3f}; "
            "the model completed-fertility statistic there is "
            f"{float(closed_slope_feedback_row['closed_root_tfr']):.3f}, and "
            "with the retained static supply curve, the implied household scale "
            f"is only "
            f"{float(closed_slope_feedback_row['closed_root_static_supply_scale']):.3f}."
        )
    else:
        closed_slope_root_sentence = (
            "After reoptimization, the independently rerun renewal audit no longer "
            "finds a closed root on the audited price range."
        )
    shell_continuation = "\\"
    readme = f"""# Sequential transition: evidence and decision menu

## Recommendation

Keep the sequential household model and add the **open-population transition** as a modest quantitative extension. The household problem is unchanged. The wrapper only tracks cohort mass through calendar time, sends births through the existing maturation delay, and clears housing at each date. A birth therefore affects current housing and tenure through the child-room requirement, then affects future adult population after maturation. This is enough to study how a fertility-preference shift changes population, housing costs, tenure, and the incidence of housing policy between two steady states.

Do not make the closed-population feedback the quantitative baseline yet. It is a useful theoretical benchmark, but the fertility response required to obtain its lower steady state is not supported by the present fit; imposing it creates large tradeoffs unless the fertility block is changed and reidentified.

**Advisor version.** I kept the full sequential tenure/fertility model and added only calendar-time population accounting. A permanent fertility-preference shift now generates a transition between open steady states, but the estimated fertility slope is too flat to support the cleaner closed-to-closed story. Across the tested closure and supply compromises, housing credit materially changes ownership and barely changes long-run population; the immediate choice is whether to present that minimal open transition or first add a one-parameter housing-stock adjustment law.

## What the two steady states are

The starting point is a constructed old **open** steady state, not the saved calibration point and not literally the U.S. economy in 2007. The fertility intercept is set to `{old_psi:.4f}`, which gives the model's top-coded completed-fertility statistic `{float(old_steady_state['tfr_topcoded']):.3f}` at the old house price. Local births nevertheless replace only `{float(old_steady_state['reproduction_ratio_B_over_E']):.3f}` of required entrants, so the provisional outside-origin share of `{100.0 * baseline_outside_share:.1f}` percent of old entrants is already needed to keep population constant.

The permanent preference shift lowers the intercept to `{new_psi:.4f}`. The existing adult cohorts and children already in the maturation queue then generate the transition; the same fixed outside flow and local-retention rule determine the new population scale. Thus the experiment is internally coherent, but the preference path and outside-entry share are normalizations, not estimated historical shocks. The separate 2007--2023 exercise is only a dated plausibility check.

This differs from the paper's earlier fixed-quota closure only at the aggregate renewal step. The quota keeps population constant by construction and remains a useful stationary policy benchmark; the fixed outside flow is the smallest change that lets population move while preserving a positive open steady state.

## What was checked

The saved sequential estimate was re-solved under the current 12-moment target contract. It is strict, has 9 free parameters, and has loss `{current_loss:.3f}`. A 25-evaluation local refit returned loss `{canonical_refit_loss:.3f}`: the target revision alone did not reveal a better nearby point. The corrected rooms-after-first-child target contributes only `{room_target_contribution:.3f}` to the loss. The three largest tensions are {largest_current_fit_text}. The rooms target is still under empirical review, so this remains a provisional fit rather than a promoted calibration.

The exact slope audit then asked how much more housing-responsive fertility must become. Holding the other parameters fixed, a closed root first appears in the tested grid when both fertility-logit noise scales are multiplied by `{first_root_factor:g}`---a `{1.0 / first_root_factor:g}`-fold increase in responsiveness. The fit loss rises from `{current_loss:.1f}` to `{first_root_fixed_loss:.1f}`. Allowing the remaining seven parameters to move for 100 local evaluations lowers that diagnostic loss to `{closed_slope_refit_loss:.1f}`, still `{closed_slope_refit_loss / current_loss:.2f}` times the saved fit. {closed_slope_root_sentence} This is a bounded diagnostic refit, not a production calibration; its complete fit and parameter tables are included. A credible closed version would need a causal housing-cost-to-fertility moment and a separate way to fit birth timing and family selection; simply lowering the logit scales is not enough.

Across the 16 open-closure cases---two housing-supply elasticities, four outside-origin entrant shares, and baseline versus LTV95---every numerical gate passes. The implied new steady-state household scale ranges from `{endpoint_scale_min:.3f}` to `{endpoint_scale_max:.3f}`, and the house-price ratio from `{endpoint_price_min:.3f}` to `{endpoint_price_max:.3f}`. Those wide ranges show that the endpoint level is not identified until the outside-entry and housing-supply objects are disciplined.

The policy conclusion is much more stable. LTV95 for dependent-child households raises steady-state ownership by `{policy_ownership_min:.3f}` to `{policy_ownership_max:.3f}` percentage points across the grid, while changing household mass by only `{policy_population_min:.4f}` to `{policy_population_max:.4f}` percent. The smallest ownership response is a saturation case with baseline ownership above 97%; among the nonsaturated cases the range is `{policy_ownership_nonsaturated_min:.3f}` to `{policy_ownership_max:.3f}` points. The intervention remains a tenure and intergenerational-access policy, not a demographic cure.

## Menu of quantitative versions

| Version | What changes | What it buys | What it gives up | Judgment |
|---|---|---|---|---|
| **0. Retained quota benchmark** | Nothing; keep the existing stationary population normalization | Preserves the current policy infrastructure and clean tenure comparisons | Population cannot respond, so it cannot answer the between-steady-states question | **Keep as benchmark/fallback** |
| **A. Open transition, retained supply** | Add cohort accounting and keep the current housing-supply elasticity; use the provisional outside-entry share | Smallest change; transition and policy incidence; preserves the full tenure/intergenerational mechanism | Endpoint population is conditional on outside entry; temporary-equilibrium prices; no welfare | **Use for the next advisor draft** |
| **B. A plus a small housing-stock law** | Let the inherited stock adjust gradually toward the retained long-run supply curve | A more credible slow housing transition and the same long-run endpoint | One adjustment-speed parameter and target; still no perfect-foresight asset pricing | **Best next extension if A is too stark** |
| **C. Closed population** | Endogenize all entry from locally born children and reidentify the fertility slope | Delivers the clean theory's price--fertility feedback | Present slope cannot support it; the preliminary root-capable refit fits much worse | **Theory benchmark, not baseline** |
| **D. 2007--2023 historical bridge** | Initialize with observed head-age shares and choose a diagnostic supply elasticity | Shows demographic momentum and the population/price timing in familiar dates | Population is imposed; the supply elasticity `{fitted_supply_elasticity:g}` is fitted to one endpoint; it misses the 2007--2011 crash | **Appendix plausibility exercise** |
| **E. Perfect-foresight asset prices** | Solve household choices against the whole future price path | Price--rent dynamics and welfare become internally coherent | A materially larger model and computation burden | **Defer** |

Version B can be one equation, $K_{{t+1}}=(1-\\lambda)K_t+\\lambda K^s(P_t)$: $\\lambda=1$ gives Version A and a small $\\lambda$ approaches inherited fixed stock. An explicit vacancy margin is needed when physical stock exceeds occupied demand.

## Reading the figure

- Panels A--B are the historical accounting exercise. The population series is matched by construction. Fixed stock overstates the 2023 price rise; the retained supply elasticity understates it; the one-moment elasticity matches the endpoint but not the housing bust. Holding the household path fixed, the imposed fertility trend changes the 2023 model price by only `{preference_price_effect:.3f}` percent relative to no fertility decline.
- Panels C--D are the clean structural experiment from the old open steady state to a lower one. Housing costs respond as preferences and current choices change, while household population turns only after the inherited child-to-entry pipeline begins to clear. That timing is the link to the paper's intergenerational mechanism. With the provisional closure, the endpoint household scale is `{endpoint['stationary_population_scale']:.3f}` and the house-price ratio is `{endpoint['price_ratio']:.3f}`.
- Panels E--F separate the paper's robust result from the fragile one: credit access materially shifts ownership but barely shifts population, while lower housing costs do not restore closed-population renewal under the saved fertility slopes. Even a fourfold increase in responsiveness reaches only `{maximum_feedback:.3f}` locally born entrants per required entrant.

A permanently fixed inherited stock has no usable terminal root over the audited positive-price range. The elastic endpoint instead has an occupied housing stock `{endpoint_stock_decline_percent:.1f}` percent below the old level. Version A lets that adjustment occur contemporaneously; Version B would spread it through a one-parameter stock law. The 2007 economy is an observed transition state, not a steady state.

All transition choices use temporary-equilibrium price expectations. The packet does not claim a perfect-foresight asset-price path, a price--rent fit, policy welfare, default pricing, or a funded LTV intervention.

## Reproducible evidence

- `decision_figure.pdf`: the six-panel summary.
- `fertility_slope_fit_tradeoff.pdf`: fit versus closed-root tradeoff.
- `provisional_calibration_target_fit.csv` and `provisional_calibration_parameters.csv`: complete saved-fit audit.
- `closed_slope_preliminary_target_fit.csv`, `closed_slope_preliminary_parameters.csv`, and `closed_slope_preliminary_feedback.csv`: complete 100-evaluation slope-refit and root audit.
- `closure_supply_policy_grid.csv`: all 16 endpoint and policy cases.
- `manifest.json`: source hashes, numerical gates, and exact inputs.

Exact bounded closed-slope refit command (from the repository root):

```bash
cd code/model
PYTHONPATH=. E5=1 E5_MATURATION_REPAIR=1 E5F=1 E5_PSI_MIN=0 {shell_continuation}
E2_SEED_RECORD=../../output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json {shell_continuation}
E5_PROBE_FIX_KE=0.1493420161260494 {shell_continuation}
E5_PROBE_FIX_KC=0.08566620151181108 E5_PROBE_SEED_PSI=0.096 {shell_continuation}
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 {shell_continuation}
.venv/bin/python intergen_eqscale_seq_optimized/run_e1_chain.py {shell_continuation}
  --outdir ../../output/model/e5f_preliminary_refits_20260816/kappa_both_factor_0p05_100eval/chain_1 {shell_continuation}
  --seed 2026081632 --max-evals 100 --minutes 45 --initial-step 0.015 {shell_continuation}
  --shrink 0.5 --J 17 --Nb 120 --max-iter-eq 10 --tol-eq 1e-4
```

Exact build command:

```bash
{command}
```
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    completed_manifest = json.loads((outdir / "manifest.json").read_text())
    completed_manifest["status"] = "complete"
    write_json(outdir / "manifest.json", completed_manifest)
    print(f"DECISION_PACKET_COMPLETE outdir={outdir}", flush=True)


if __name__ == "__main__":
    main()
