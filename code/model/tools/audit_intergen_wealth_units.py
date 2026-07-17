#!/usr/bin/env python3
"""Audit wealth, income, rent, and house-price units in the intergen model.

The active one-market model uses 4-year periods and scales income flows to the
period. This script makes the denominator choice explicit so wealth targets are
not accidentally compared across annual and period units.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))


DEFAULT_CACHE = ROOT / "output/model/intergen_current_review/solution_cache.pkl"
DEFAULT_OUTDIR = ROOT / "output/model/intergen_current_review/wealth_units_audit"
DEFAULT_ENTRY_DATA = (
    ROOT
    / "code/data/psid_followup_mar2026/output/entry_wealth_v1/entry_wealth_candidate_targets_v1.csv"
)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sol, P, p_eq = load_solution(args.solution_cache)
    sample_rows = sample_wealth_rows(sol, P, args.entry_data)
    unit_rows = unit_accounting_rows(sol, P, p_eq)
    price_rows = price_accounting_rows(sol, P, p_eq)
    write_csv(args.outdir / "wealth_sample_unit_comparison.csv", sample_rows)
    write_csv(args.outdir / "model_unit_accounting.csv", unit_rows)
    write_csv(args.outdir / "price_rent_accounting.csv", price_rows)
    write_readme(args.outdir, args.solution_cache, args.entry_data, sample_rows, unit_rows, price_rows)
    print(f"Wrote wealth/unit audit to {args.outdir}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solution-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--entry-data", type=Path, default=DEFAULT_ENTRY_DATA)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def load_solution(path: Path) -> tuple[Any, Any, np.ndarray]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    baseline = payload["baseline"]
    return baseline["sol"], baseline["P"], np.asarray(baseline["p_eq"], dtype=float).reshape(-1)


def sample_wealth_rows(sol: Any, P: Any, data_path: Path) -> list[dict[str, Any]]:
    data = read_data_targets(data_path)
    specs = [
        ("young_all_25_30", 25.0, 30.0, False, False),
        ("young_childless_25_35", 25.0, 35.0, True, False),
        ("young_childless_rent_25_35", 25.0, 35.0, True, True),
        ("young_prebirth_rent_25_35", 25.0, 35.0, False, True),
    ]
    rows: list[dict[str, Any]] = []
    for sample, age_lo, age_hi, childless_only, renter_only in specs:
        stats = model_sample_stats(sol, P, age_lo, age_hi, childless_only=childless_only, renter_only=renter_only)
        target = data.get((sample, "liq_nw_to_inc"), {})
        rows.append(
            {
                "sample": sample,
                "psid_liq_nw_to_inc_mean": maybe_float(target.get("mean")),
                "psid_liq_nw_to_inc_median": maybe_float(target.get("p50")),
                "psid_N": target.get("N", ""),
                "model_mass": stats["mass"],
                "model_mean_b": stats["mean_b"],
                "model_median_b": stats["median_b"],
                "model_mean_ratio_period_aftertax_income": stats["mean_ratio_period_aftertax"],
                "model_median_ratio_period_aftertax_income": stats["median_ratio_period_aftertax"],
                "model_mean_ratio_annual_aftertax_income": stats["mean_ratio_annual_aftertax"],
                "model_median_ratio_annual_aftertax_income": stats["median_ratio_annual_aftertax"],
                "model_mean_ratio_annual_gross_income": stats["mean_ratio_annual_gross"],
                "model_median_ratio_annual_gross_income": stats["median_ratio_annual_gross"],
            }
        )
    return rows


def model_sample_stats(
    sol: Any,
    P: Any,
    age_lo: float,
    age_hi: float,
    *,
    childless_only: bool,
    renter_only: bool,
) -> dict[str, float]:
    g = np.asarray(sol.g, dtype=float)
    bg = np.asarray(sol.b_grid, dtype=float)
    ages = float(P.age_start) + np.arange(int(P.J)) * float(P.da)
    z_grid = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
    tau = float(getattr(P, "tau_pay", 0.0))
    vals_b: list[float] = []
    ratio_period: list[float] = []
    ratio_annual_aftertax: list[float] = []
    ratio_annual_gross: list[float] = []
    weights: list[float] = []
    for ib, b in enumerate(bg):
        for j, age in enumerate(ages):
            if age < age_lo or age > age_hi:
                continue
            for zz, z in enumerate(z_grid):
                period_income = float(P.income[0, j]) * (float(z) if j < int(getattr(P, "J_R", P.J)) else 1.0)
                annual_aftertax = period_income / max(period_years, 1e-12)
                annual_gross = annual_aftertax / max(1.0 - tau, 1e-12)
                if childless_only and renter_only:
                    mass = float(np.sum(g[ib, 0, 0, j, zz, 0, :]))
                elif childless_only:
                    mass = float(np.sum(g[ib, :, 0, j, zz, 0, :]))
                elif renter_only:
                    mass = float(np.sum(g[ib, 0, 0, j, zz, :, :]))
                else:
                    mass = float(np.sum(g[ib, :, 0, j, zz, :, :]))
                if mass <= 0:
                    continue
                vals_b.append(float(b))
                ratio_period.append(float(b) / period_income if period_income > 0 else math.nan)
                ratio_annual_aftertax.append(float(b) / annual_aftertax if annual_aftertax > 0 else math.nan)
                ratio_annual_gross.append(float(b) / annual_gross if annual_gross > 0 else math.nan)
                weights.append(mass)
    if not weights:
        return {k: math.nan for k in (
            "mass",
            "mean_b",
            "median_b",
            "mean_ratio_period_aftertax",
            "median_ratio_period_aftertax",
            "mean_ratio_annual_aftertax",
            "median_ratio_annual_aftertax",
            "mean_ratio_annual_gross",
            "median_ratio_annual_gross",
        )}
    w = np.asarray(weights, dtype=float)
    b_arr = np.asarray(vals_b, dtype=float)
    rp = np.asarray(ratio_period, dtype=float)
    rat = np.asarray(ratio_annual_aftertax, dtype=float)
    rg = np.asarray(ratio_annual_gross, dtype=float)
    return {
        "mass": float(np.sum(w)),
        "mean_b": weighted_mean(b_arr, w),
        "median_b": weighted_median(b_arr, w),
        "mean_ratio_period_aftertax": weighted_mean(rp, w),
        "median_ratio_period_aftertax": weighted_median(rp, w),
        "mean_ratio_annual_aftertax": weighted_mean(rat, w),
        "median_ratio_annual_aftertax": weighted_median(rat, w),
        "mean_ratio_annual_gross": weighted_mean(rg, w),
        "median_ratio_annual_gross": weighted_median(rg, w),
    }


def unit_accounting_rows(sol: Any, P: Any, p_eq: np.ndarray) -> list[dict[str, Any]]:
    period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
    tau = float(getattr(P, "tau_pay", 0.0))
    prime_j = int(np.argmax(np.asarray(P.income_age_profile, dtype=float)))
    period_income = float(P.income[0, prime_j])
    annual_aftertax = period_income / period_years
    annual_gross = annual_aftertax / max(1.0 - tau, 1e-12)
    b_entry = float(getattr(P, "b_entry_fixed", math.nan))
    entry_ratio_mean = float(getattr(sol, "entry_wealth_ratio_mean", math.nan))
    entry_grid_mean = float(getattr(sol, "entry_wealth_grid_mean", math.nan))
    return [
        {"object": "period_years", "value": period_years, "interpretation": "years per model decision period"},
        {"object": "scale_flows_to_period", "value": bool(getattr(P, "scale_flows_to_period", False)), "interpretation": "income and flow costs are multiplied by period_years"},
        {"object": "tau_pay", "value": tau, "interpretation": "payroll/income tax wedge applied to working income"},
        {"object": "prime_period_aftertax_income", "value": period_income, "interpretation": "model denominator currently used in many wealth/income stats"},
        {"object": "prime_annual_aftertax_income", "value": annual_aftertax, "interpretation": "period after-tax income divided by period_years"},
        {"object": "prime_annual_gross_income", "value": annual_gross, "interpretation": "annual after-tax income grossed up by 1/(1-tau_pay)"},
        {"object": "entry_wealth_mode", "value": str(getattr(P, "entry_wealth_mode", "scalar")), "interpretation": "entrant wealth closure used by the solver"},
        {"object": "entry_wealth_ratio_mean", "value": entry_ratio_mean, "interpretation": "mean entrant wealth / annual gross income when using external ratio mode"},
        {"object": "entry_wealth_grid_mean", "value": entry_grid_mean, "interpretation": "mean entrant liquid wealth in model goods units after grid scattering"},
        {"object": "b_entry_fixed", "value": b_entry, "interpretation": "legacy scalar entrant stock wealth in model goods units"},
        {"object": "b_entry_over_period_aftertax_income", "value": b_entry / period_income, "interpretation": "entry wealth under current period denominator"},
        {"object": "b_entry_over_annual_aftertax_income", "value": b_entry / annual_aftertax, "interpretation": "entry wealth vs annual after-tax income"},
        {"object": "b_entry_over_annual_gross_income", "value": b_entry / annual_gross, "interpretation": "entry wealth vs annual gross-normalized income"},
        {"object": "reported_young_liquid_wealth_to_income", "value": float(getattr(sol, "young_liquid_wealth_to_income", math.nan)), "interpretation": "current model statistic; period after-tax denominator"},
        {"object": "reported_young_liquid_wealth_to_annual_aftertax_rough", "value": float(getattr(sol, "young_liquid_wealth_to_income", math.nan)) * period_years, "interpretation": "rough rescale if same sample denominator is annual after-tax"},
        {"object": "reported_young_liquid_wealth_to_annual_gross_rough", "value": float(getattr(sol, "young_liquid_wealth_to_income", math.nan)) * period_years * (1.0 - tau), "interpretation": "rough rescale if same sample denominator is annual gross"},
    ]


def price_accounting_rows(sol: Any, P: Any, p_eq: np.ndarray) -> list[dict[str, Any]]:
    period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
    p = float(p_eq[0])
    user_cost_rate = float(getattr(P, "user_cost_rate", math.nan))
    user_cost_period = user_cost_rate * p
    user_cost_annual = user_cost_period / period_years
    rows: list[dict[str, Any]] = [
        {"object": "house_price_per_room_p", "value": p, "interpretation": "asset price of one room/service unit"},
        {"object": "user_cost_rate_4yr", "value": user_cost_rate, "interpretation": "q + delta + tau_H over one 4-year period"},
        {"object": "rent_user_cost_per_room_period", "value": user_cost_period, "interpretation": "per-room rent/user cost over model period"},
        {"object": "rent_user_cost_per_room_annual", "value": user_cost_annual, "interpretation": "per-room annual rent/user cost"},
    ]
    phi = np.asarray(getattr(P, "phi", [0.8]), dtype=float).reshape(-1)
    phi0 = float(phi[0]) if phi.size else 0.8
    for H in np.asarray(P.H_own, dtype=float).reshape(-1):
        price = p * float(H)
        down = (1.0 - phi0) * price
        rent_annual = user_cost_annual * float(H)
        rows.append(
            {
                "object": f"owner_rung_H{H:g}",
                "value": price,
                "interpretation": f"asset price; down payment at phi={phi0:.3f} is {down:.3f}; annual user cost is {rent_annual:.3f}",
            }
        )
    return rows


def read_data_targets(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {(r.get("sample", ""), r.get("moment", "")): r for r in rows}


def write_readme(
    outdir: Path,
    cache_path: Path,
    data_path: Path,
    sample_rows: list[dict[str, Any]],
    unit_rows: list[dict[str, Any]],
    price_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# Wealth And Unit Audit",
        "",
        "Purpose: make denominator and numeraire choices explicit before changing",
        "wealth targets or recalibrating entry wealth.",
        "",
        f"Solution cache: `{relative(cache_path)}`",
        f"PSID entry-wealth summary: `{relative(data_path)}`",
        "",
        "## Main Unit Accounting",
        "",
        "| Object | Value | Interpretation |",
        "|---|---:|---|",
    ]
    for row in unit_rows:
        lines.append(f"| `{row['object']}` | {format_value(row['value'])} | {row['interpretation']} |")
    lines.extend(
        [
            "",
            "## Wealth Targets Under Different Denominators",
            "",
            "| Sample | PSID Mean | PSID Median | Model Mean: Period After-Tax | Model Median: Period After-Tax | Model Mean: Annual Gross | Model Median: Annual Gross |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sample_rows:
        lines.append(
            "| `{sample}` | {dmean} | {dmed} | {mp} | {medp} | {mg} | {medg} |".format(
                sample=row["sample"],
                dmean=format_value(row["psid_liq_nw_to_inc_mean"]),
                dmed=format_value(row["psid_liq_nw_to_inc_median"]),
                mp=format_value(row["model_mean_ratio_period_aftertax_income"]),
                medp=format_value(row["model_median_ratio_period_aftertax_income"]),
                mg=format_value(row["model_mean_ratio_annual_gross_income"]),
                medg=format_value(row["model_median_ratio_annual_gross_income"]),
            )
        )
    lines.extend(
        [
            "",
            "## Price/Rent Scale",
            "",
            "| Object | Value | Interpretation |",
            "|---|---:|---|",
        ]
    )
    for row in price_rows:
        lines.append(f"| `{row['object']}` | {format_value(row['value'])} | {row['interpretation']} |")
    lines.extend(
        [
            "",
            "## Read",
            "",
            "The current model wealth/income statistic uses period after-tax income as",
            "the denominator. Existing PSID targets appear to be annual-income ratios.",
            "Therefore wealth targets should be converted or the model statistic should",
            "be redefined before entry wealth is recalibrated.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")
    manifest = {
        "status": "diagnostic_only_not_calibration",
        "solution_cache": str(cache_path),
        "entry_data": str(data_path),
        "files": [
            "README.md",
            "wealth_sample_unit_comparison.csv",
            "model_unit_accounting.csv",
            "price_rent_accounting.csv",
        ],
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    finite = np.isfinite(values) & np.isfinite(weights)
    if not np.any(finite):
        return math.nan
    return float(np.sum(values[finite] * weights[finite]) / np.sum(weights[finite]))


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    finite = np.isfinite(values) & np.isfinite(weights)
    if not np.any(finite):
        return math.nan
    x = values[finite]
    w = weights[finite]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    return float(x[np.searchsorted(np.cumsum(w), 0.5 * float(np.sum(w)), side="left")])


def maybe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def format_value(value: Any) -> str:
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    try:
        x = float(value)
    except Exception:
        return str(value)
    if not math.isfinite(x):
        return ""
    return f"{x:.3f}"


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except Exception:
        return str(path)


if __name__ == "__main__":
    main()
