#!/usr/bin/env python3
"""Match PSID and model decompositions behind the three estate diagnostics.

This is a diagnostic report, not a calibration-target proposal. It re-solves
a strict M3 or M4 winner and decomposes the wealth/income ratios into estate
levels, retirement-income levels, housing and nonhousing components, and
completed-fertility groups. PSID-only marital groups are retained explicitly;
the live model has no marital state.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import numpy as np

from intergen_housing_fertility.calibration import extract_moments
from intergen_housing_fertility.local_panel import jsonable
from intergen_housing_fertility.solver import annual_gross_income_at_state, run_model_cp_dt
from tools.profile_intergen_bequest_reachability import DISPERSION, FAMILY_GAP, MEDIAN
from tools.run_intergen_bequest_exit_chain import arm_contract, common_overrides, load_theta


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_WINNER = (
    ROOT / "output/model/intergen_internal_bequest_recalibration_20260715/report/results.json"
)
DEFAULT_PSID = (
    ROOT
    / "code/data/psid_followup_mar2026/output/intergen_bequest_distribution_diagnostic"
    / "psid_distribution_decomposition.csv"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_bequest_distribution_diagnostic_20260716"

LEVEL_COLUMNS = (
    "income_p10",
    "income_p50",
    "income_p90",
    "total_estate_p50",
    "total_estate_p75",
    "total_estate_p90",
    "nonhousing_p50",
    "nonhousing_p75",
    "nonhousing_p90",
    "housing_p50",
    "housing_p75",
    "housing_p90",
    "liquid_or_debt_p50",
    "liquid_or_debt_p75",
    "liquid_or_debt_p90",
    "gross_housing_p50",
    "gross_housing_p75",
    "gross_housing_p90",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--winner-json", type=Path, default=DEFAULT_WINNER)
    parser.add_argument("--winner-arm", type=str, default="M3")
    parser.add_argument("--arm", choices=("M3", "M4"), default="M3")
    parser.add_argument("--psid-csv", type=Path, default=DEFAULT_PSID)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(ok):
        return math.nan
    x = np.asarray(values[ok], dtype=float)
    w = np.asarray(weights[ok], dtype=float)
    order = np.argsort(x, kind="mergesort")
    x, w = x[order], w[order]
    idx = int(np.searchsorted(np.cumsum(w), probability * float(np.sum(w)), side="left"))
    return float(x[min(idx, x.size - 1)])


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(ok):
        return math.nan
    return float(np.sum(values[ok] * weights[ok]) / np.sum(weights[ok]))


def weighted_covariance(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(ok):
        return math.nan
    xx, yy, ww = x[ok], y[ok], weights[ok]
    mx = float(np.sum(ww * xx) / np.sum(ww))
    my = float(np.sum(ww * yy) / np.sum(ww))
    return float(np.sum(ww * (xx - mx) * (yy - my)) / np.sum(ww))


def weighted_correlation(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    cov = weighted_covariance(x, y, weights)
    vx = weighted_covariance(x, x, weights)
    vy = weighted_covariance(y, y, weights)
    if not all(math.isfinite(v) for v in (cov, vx, vy)) or vx <= 0.0 or vy <= 0.0:
        return math.nan
    return float(cov / math.sqrt(vx * vy))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def build_model_micro(sol: Any, P: Any, price: np.ndarray) -> dict[str, np.ndarray]:
    weights7 = np.asarray(sol.g_beginning_assets_by_current_choice, dtype=float)
    if weights7.ndim != 7:
        raise ValueError("M3 diagnostic requires the seven-dimensional Markov distribution")
    nb, nt, nloc, nj, nz, npar, ncs = weights7.shape
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    if b_grid.size != nb:
        raise ValueError("wealth grid and distribution disagree")
    price_vec = np.asarray(price, dtype=float).reshape(-1)

    b = np.broadcast_to(b_grid[:, None, None, None, None, None, None], weights7.shape)
    housing_by_tenure_location = np.zeros((nt, nloc), dtype=float)
    for ten in range(1, nt):
        housing_by_tenure_location[ten, :] = price_vec * float(P.H_own[ten - 1])
    housing = np.broadcast_to(
        housing_by_tenure_location[None, :, :, None, None, None, None], weights7.shape
    )
    income_state = np.zeros((nloc, nj, nz), dtype=float)
    z_grid = np.asarray(P.z_grid, dtype=float).reshape(-1)
    for i in range(nloc):
        for j in range(nj):
            for zz, z_value in enumerate(z_grid):
                income_state[i, j, zz] = annual_gross_income_at_state(P, i, j, float(z_value))
    income = np.broadcast_to(income_state[None, None, :, :, :, None, None], weights7.shape)
    ages = np.broadcast_to(
        (float(P.age_start) + np.arange(nj) * float(P.da))[None, None, None, :, None, None, None],
        weights7.shape,
    )
    children = np.broadcast_to(
        np.arange(npar, dtype=float)[None, None, None, None, None, :, None], weights7.shape
    )
    positive = np.isfinite(weights7) & (weights7 > 0.0)
    total_estate = b + housing
    # Empirical home equity is net of mortgage debt.  In owner states the
    # model's negative b is collateralized debt, so allocate it against gross
    # housing before comparing components.  Positive owner b remains
    # nonhousing wealth.  Renter debt remains nonhousing debt.  This mapping
    # sums exactly to the model target object b+pH.
    owner = housing > 0.0
    housing_equity = np.where(owner, housing + np.minimum(b, 0.0), 0.0)
    nonhousing = np.where(owner, np.maximum(b, 0.0), b)
    if not np.allclose(housing_equity + nonhousing, total_estate, rtol=0.0, atol=1e-12):
        raise RuntimeError("model component mapping does not sum to b+pH")
    return {
        "weight": weights7[positive],
        "age": ages[positive],
        "children": children[positive],
        "income": income[positive],
        "nonhousing": nonhousing[positive],
        "housing": housing_equity[positive],
        "liquid_or_debt": b[positive],
        "gross_housing": housing[positive],
        "total_estate": total_estate[positive],
    }


def summarize_model_group(
    micro: dict[str, np.ndarray],
    mask: np.ndarray,
    window_name: str,
    group_name: str,
    benchmark_income_median: float,
    window_mass: float,
) -> dict[str, Any]:
    x = {name: values[mask] for name, values in micro.items()}
    w = x["weight"]
    ratio = x["total_estate"] / x["income"]
    quantile = lambda values, probability: weighted_quantile(values, w, probability)
    estate_p90 = quantile(x["total_estate"], 0.90)
    top = x["total_estate"] >= estate_p90
    top_housing_share = float(
        np.sum(w[top] * x["housing"][top]) / np.sum(w[top] * x["total_estate"][top])
    )
    result: dict[str, Any] = {
        "dataset": "MODEL",
        "window": window_name,
        "group": group_name,
        "n_rows": int(np.sum(mask)),
        "n_people": "",
        "mass_or_weight": float(np.sum(w)),
        "share_of_window_weight": float(np.sum(w) / window_mass),
        "benchmark_income_median": benchmark_income_median,
        "income_p10": quantile(x["income"], 0.10),
        "income_p50": quantile(x["income"], 0.50),
        "income_p90": quantile(x["income"], 0.90),
        "total_estate_p50": quantile(x["total_estate"], 0.50),
        "total_estate_p75": quantile(x["total_estate"], 0.75),
        "total_estate_p90": estate_p90,
        "nonhousing_p50": quantile(x["nonhousing"], 0.50),
        "nonhousing_p75": quantile(x["nonhousing"], 0.75),
        "nonhousing_p90": quantile(x["nonhousing"], 0.90),
        "housing_p50": quantile(x["housing"], 0.50),
        "housing_p75": quantile(x["housing"], 0.75),
        "housing_p90": quantile(x["housing"], 0.90),
        "liquid_or_debt_p50": quantile(x["liquid_or_debt"], 0.50),
        "liquid_or_debt_p75": quantile(x["liquid_or_debt"], 0.75),
        "liquid_or_debt_p90": quantile(x["liquid_or_debt"], 0.90),
        "gross_housing_p50": quantile(x["gross_housing"], 0.50),
        "gross_housing_p75": quantile(x["gross_housing"], 0.75),
        "gross_housing_p90": quantile(x["gross_housing"], 0.90),
        "estate_income_ratio_p50": quantile(ratio, 0.50),
        "estate_income_ratio_p75": quantile(ratio, 0.75),
        "estate_income_ratio_p90": quantile(ratio, 0.90),
        "estate_income_ratio_p90_p50": quantile(ratio, 0.90) / quantile(ratio, 0.50),
        "total_estate_income_covariance": weighted_covariance(x["total_estate"], x["income"], w),
        "normalized_estate_income_covariance": weighted_covariance(
            x["total_estate"], x["income"], w
        )
        / benchmark_income_median**2,
        "estate_income_correlation": weighted_correlation(x["total_estate"], x["income"], w),
        "asinh_estate_log_income_correlation": weighted_correlation(
            np.arcsinh(x["total_estate"] / benchmark_income_median),
            np.log(x["income"] / benchmark_income_median),
            w,
        ),
        "mean_housing_share": weighted_mean(x["housing"], w)
        / weighted_mean(x["total_estate"], w),
        "top_estate_decile_housing_share": top_housing_share,
    }
    for column in LEVEL_COLUMNS:
        result[f"{column}_over_benchmark_income"] = float(result[column]) / benchmark_income_median
    return result


def model_rows(sol: Any, P: Any, price: np.ndarray) -> list[dict[str, Any]]:
    micro = build_model_micro(sol, P, price)
    windows = {"ages_65_75": (65.0, 75.0), "ages_76_84": (76.0, 84.0)}
    groups: dict[str, Callable[[np.ndarray], np.ndarray]] = {
        "all": lambda children: np.ones(children.shape, dtype=bool),
        "children_0": lambda children: children == 0,
        "children_1": lambda children: children == 1,
        "children_2plus": lambda children: children >= 2,
    }
    rows: list[dict[str, Any]] = []
    for window_name, (age_lo, age_hi) in windows.items():
        window_mask = (micro["age"] >= age_lo) & (micro["age"] <= age_hi)
        window_mass = float(np.sum(micro["weight"][window_mask]))
        benchmark_income_median = weighted_quantile(
            micro["income"][window_mask], micro["weight"][window_mask], 0.50
        )
        for group_name, group_fn in groups.items():
            mask = window_mask & group_fn(micro["children"])
            rows.append(
                summarize_model_group(
                    micro, mask, window_name, group_name, benchmark_income_median, window_mass
                )
            )
    return rows


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def numeric_row(row: dict[str, Any]) -> dict[str, Any]:
    converted: dict[str, Any] = {}
    text_keys = {"dataset", "window", "group"}
    for key, value in row.items():
        if key in text_keys or value in (None, ""):
            converted[key] = value
            continue
        try:
            converted[key] = float(value)
        except (TypeError, ValueError):
            converted[key] = value
    return converted


def main() -> None:
    args = parse_args()
    theta = load_theta(args.winner_json, args.winner_arm)
    theta["tenure_choice_kappa"] = 0.0
    if args.arm == "M4" and float(theta.get("theta_n", math.nan)) != 0.0:
        raise ValueError("M4 distribution diagnostic requires theta_n=0")
    contract = SimpleNamespace(
        arm=args.arm,
        ltv_terminal=0.0,
        theta1=float(theta["theta1"]),
        seed_theta0=float(theta["theta0"]),
        fixed_theta0=None,
        J=17,
        Nb=120,
        max_iter_eq=40,
        tol_eq=2.5e-5,
    )
    _, _, mechanism = arm_contract(contract)
    overrides = common_overrides(contract, mechanism)
    sol, P, price = run_model_cp_dt({**overrides, **theta}, verbose=not args.quiet)
    residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
    strict = bool(
        getattr(sol, "timings", {}).get("strict_converged", getattr(sol, "converged", False))
        and residual <= float(P.tol_eq)
    )
    if not strict:
        raise RuntimeError(f"{args.arm} diagnostic solve failed strict residual gate: {residual:.6g}")

    moments = extract_moments(sol, P)
    rows_model = model_rows(sol, P, price)
    lookup = {(row["window"], row["group"]): row for row in rows_model}
    all_old = lookup[("ages_76_84", "all")]
    one = lookup[("ages_65_75", "children_1")]
    two_plus = lookup[("ages_65_75", "children_2plus")]
    reproduced = {
        MEDIAN: float(all_old["estate_income_ratio_p50"]),
        DISPERSION: float(all_old["estate_income_ratio_p90_p50"]),
        FAMILY_GAP: float(two_plus["estate_income_ratio_p50"])
        - float(one["estate_income_ratio_p50"]),
    }
    for name, value in reproduced.items():
        if not math.isclose(value, float(moments[name]), rel_tol=0.0, abs_tol=1e-10):
            raise RuntimeError(f"distribution decomposition does not reproduce {name}")

    rows_psid = [numeric_row(row) for row in read_csv(args.psid_csv)]
    combined = rows_psid + rows_model
    headline_keys = {
        ("ages_76_84", "all"),
        ("ages_65_75", "children_1"),
        ("ages_65_75", "children_2plus"),
    }
    headline = [row for row in combined if (row["window"], row["group"]) in headline_keys]
    marital = [row for row in rows_psid if "married" in str(row["group"])]
    psid_lookup = {(row["window"], row["group"]): row for row in rows_psid}
    model_lookup = {(row["window"], row["group"]): row for row in rows_model}
    psid_old = psid_lookup[("ages_76_84", "all")]
    model_old = model_lookup[("ages_76_84", "all")]
    psid_one = psid_lookup[("ages_65_75", "children_1")]
    psid_two = psid_lookup[("ages_65_75", "children_2plus")]
    model_one = model_lookup[("ages_65_75", "children_1")]
    model_two = model_lookup[("ages_65_75", "children_2plus")]
    psid_married = psid_lookup[("ages_76_84", "married")]
    psid_unmarried = psid_lookup[("ages_76_84", "not_married")]

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "model_distribution_decomposition.csv", rows_model)
    write_csv(args.outdir / "matched_distribution_decomposition.csv", combined)
    write_csv(args.outdir / "headline_comparison.csv", headline)
    write_csv(args.outdir / "psid_marital_composition.csv", marital)
    summary = {
        "status": "complete",
        "role": "diagnostic_only_not_calibration_targets",
        "arm": args.arm,
        "winner_json": str(args.winner_json),
        "strict_model_solve": strict,
        "market_residual": residual,
        "price": jsonable(price),
        "target_reproduction": {
            name: {"model_statistic": float(moments[name]), "decomposition": value}
            for name, value in reproduced.items()
        },
        "model_has_marital_state": False,
        "model_rows": rows_model,
    }
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    (args.outdir / "README.md").write_text(
        "\n".join(
            [
                "# Matched estate-distribution diagnostic",
                "",
                "These numbers decompose the existing three estate targets; they are not proposed calibration moments.",
                "",
                f"The strict {args.arm} winner was re-solved at residual `{residual:.3g}`. The decomposition reproduces all three model moments to `1e-10`.",
                "PSID raw levels are real 2022 USD. Model raw levels are normalized model units; columns ending in",
                "`_over_benchmark_income` put both datasets in comparable within-window annual-income units.",
                "The live model has no marital state, so marital decompositions are PSID-only diagnostics.",
                "For the component comparison, negative liquid wealth in owner states is allocated to mortgage debt:",
                "model home equity is `pH + min(b,0)` and model nonhousing wealth is `max(b,0)` for owners;",
                "renter `b` remains nonhousing wealth. Raw model `b` and gross `pH` quantiles are also saved.",
                "",
                "## Headline decomposition",
                "",
                "All level quantiles below are divided by the all-household median income in the same age window.",
                "",
                "| Ages 76--84 | PSID | Model |",
                "|---|---:|---:|",
                f"| Income p10 / p50 / p90 | {float(psid_old['income_p10_over_benchmark_income']):.3f} / 1.000 / {float(psid_old['income_p90_over_benchmark_income']):.3f} | {float(model_old['income_p10_over_benchmark_income']):.3f} / 1.000 / {float(model_old['income_p90_over_benchmark_income']):.3f} |",
                f"| Total estate p50 / p75 / p90 | {float(psid_old['total_estate_p50_over_benchmark_income']):.3f} / {float(psid_old['total_estate_p75_over_benchmark_income']):.3f} / {float(psid_old['total_estate_p90_over_benchmark_income']):.3f} | {float(model_old['total_estate_p50_over_benchmark_income']):.3f} / {float(model_old['total_estate_p75_over_benchmark_income']):.3f} / {float(model_old['total_estate_p90_over_benchmark_income']):.3f} |",
                f"| Nonhousing p50 / p75 / p90 | {float(psid_old['nonhousing_p50_over_benchmark_income']):.3f} / {float(psid_old['nonhousing_p75_over_benchmark_income']):.3f} / {float(psid_old['nonhousing_p90_over_benchmark_income']):.3f} | {float(model_old['nonhousing_p50_over_benchmark_income']):.3f} / {float(model_old['nonhousing_p75_over_benchmark_income']):.3f} / {float(model_old['nonhousing_p90_over_benchmark_income']):.3f} |",
                f"| Housing p50 / p75 / p90 | {float(psid_old['housing_p50_over_benchmark_income']):.3f} / {float(psid_old['housing_p75_over_benchmark_income']):.3f} / {float(psid_old['housing_p90_over_benchmark_income']):.3f} | {float(model_old['housing_p50_over_benchmark_income']):.3f} / {float(model_old['housing_p75_over_benchmark_income']):.3f} / {float(model_old['housing_p90_over_benchmark_income']):.3f} |",
                f"| Estate/income p50 / p75 / p90 | {float(psid_old['estate_income_ratio_p50']):.3f} / {float(psid_old['estate_income_ratio_p75']):.3f} / {float(psid_old['estate_income_ratio_p90']):.3f} | {float(model_old['estate_income_ratio_p50']):.3f} / {float(model_old['estate_income_ratio_p75']):.3f} / {float(model_old['estate_income_ratio_p90']):.3f} |",
                f"| Estate--income correlation | {float(psid_old['estate_income_correlation']):.3f} | {float(model_old['estate_income_correlation']):.3f} |",
                f"| Top-estate-decile housing share | {float(psid_old['top_estate_decile_housing_share']):.3f} | {float(model_old['top_estate_decile_housing_share']):.3f} |",
                "",
                "The PSID raw total-estate p90/p50 is "
                f"`{float(psid_old['total_estate_p90']) / float(psid_old['total_estate_p50']):.3f}` versus "
                f"`{float(model_old['total_estate_p90']) / float(model_old['total_estate_p50']):.3f}` in the model. "
                "Positive estate--income correlation compresses the PSID ratio distribution rather than creating its tail.",
                "",
                "## Family-size and marital decomposition",
                "",
                "| Ages 65--75 estate/income | PSID: 1 child | PSID: 2+ | Model: 1 child | Model: 2+ |",
                "|---|---:|---:|---:|---:|",
                f"| p50 | {float(psid_one['estate_income_ratio_p50']):.3f} | {float(psid_two['estate_income_ratio_p50']):.3f} | {float(model_one['estate_income_ratio_p50']):.3f} | {float(model_two['estate_income_ratio_p50']):.3f} |",
                f"| p75 | {float(psid_one['estate_income_ratio_p75']):.3f} | {float(psid_two['estate_income_ratio_p75']):.3f} | {float(model_one['estate_income_ratio_p75']):.3f} | {float(model_two['estate_income_ratio_p75']):.3f} |",
                f"| p90 | {float(psid_one['estate_income_ratio_p90']):.3f} | {float(psid_two['estate_income_ratio_p90']):.3f} | {float(model_one['estate_income_ratio_p90']):.3f} | {float(model_two['estate_income_ratio_p90']):.3f} |",
                "",
                f"At ages 76--84, married households are {100 * float(psid_married['share_of_window_weight']):.1f}% of PSID weight and have ratio p90/p50 `{float(psid_married['estate_income_ratio_p90_p50']):.3f}`; nonmarried households are {100 * float(psid_unmarried['share_of_window_weight']):.1f}% and have `{float(psid_unmarried['estate_income_ratio_p90_p50']):.3f}`. Both groups remain highly dispersed, so marital composition matters but does not explain the entire tail.",
                "",
                "Files:",
                "- `headline_comparison.csv`: the three target groups, PSID versus model.",
                "- `matched_distribution_decomposition.csv`: all PSID and model groups.",
                "- `psid_marital_composition.csv`: PSID marital and parity-by-marital groups.",
                "- `model_distribution_decomposition.csv`: model-only state decomposition.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
