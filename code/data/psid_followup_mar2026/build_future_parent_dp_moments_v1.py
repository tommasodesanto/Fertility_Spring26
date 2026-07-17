"""Summarize future-parent down-payment moments for young childless renters.

The companion Stata script constructs the audited ID-level entry sample from
PSID person-years. This script keeps the aggregation logic transparent and
easier to rerun without Stata once the sample CSV exists.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUTDIR = ROOT / "code" / "data" / "psid_followup_mar2026" / "output" / "future_parent_dp_moments_v1"
ENTRY = OUTDIR / "future_parent_dp_entry_sample_v1.csv"
MOMENTS = OUTDIR / "future_parent_dp_moments_v1.csv"
COMPARISON = OUTDIR / "future_parent_dp_comparison_v1.csv"
BINS = OUTDIR / "future_parent_dp_shortfall_bins_v1.csv"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(ENTRY)
    df = add_downpayment_measures(df)

    rows = []
    comparison_rows = []
    bin_rows = []

    for horizon in [35, 40, 45]:
        observed = df[f"observed_by_{horizon}"].eq(1)
        sample = df.loc[observed].copy()
        if sample.empty:
            continue

        group_defs = [
            ("all_entry_childless_renters", np.ones(len(sample), dtype=bool)),
            (f"future_parent_by_{horizon}", sample[f"parent_by_{horizon}"].eq(1).to_numpy()),
            (f"remain_childless_by_{horizon}", sample[f"parent_by_{horizon}"].eq(0).to_numpy()),
        ]
        group_stats = {}
        for group, mask in group_defs:
            stats = summarize(sample.loc[mask], horizon)
            if stats is None:
                continue
            stats.update({"horizon": horizon, "group": group})
            rows.append(stats)
            group_stats[group] = stats

        parent_key = f"future_parent_by_{horizon}"
        childless_key = f"remain_childless_by_{horizon}"
        if parent_key in group_stats and childless_key in group_stats:
            comparison_rows.extend(compare_groups(group_stats[parent_key], group_stats[childless_key], horizon))

        bin_rows.extend(summarize_bins(sample, horizon))

    moments = pd.DataFrame(rows)
    moments.sort_values(["horizon", "group"], inplace=True)
    moments.to_csv(MOMENTS, index=False)

    comparison = pd.DataFrame(comparison_rows)
    comparison.sort_values(["horizon", "metric"], inplace=True)
    comparison.to_csv(COMPARISON, index=False)

    bins = pd.DataFrame(bin_rows)
    bins.sort_values(["horizon", "shortfall20_bin"], inplace=True)
    bins.to_csv(BINS, index=False)

    print(f"Wrote {len(moments):,} rows to {MOMENTS}")
    print(f"Wrote {len(comparison):,} rows to {COMPARISON}")
    print(f"Wrote {len(bins):,} rows to {BINS}")


def add_downpayment_measures(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in [
        "liquid_nw_mean",
        "income_mean",
        "p_target_mean",
        "req_dp20_mean",
        "req_dp035_mean",
        "weight_entry",
    ]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["dp_gap20"] = out["req_dp20_mean"] - out["liquid_nw_mean"]
    out["dp_gap035"] = out["req_dp035_mean"] - out["liquid_nw_mean"]
    out["shortfall20"] = out["dp_gap20"].clip(lower=0.0)
    out["shortfall035"] = out["dp_gap035"].clip(lower=0.0)
    valid_income = out["income_mean"].gt(1000) & out["income_mean"].notna()
    out["shortfall20_to_income"] = np.where(valid_income, out["shortfall20"] / out["income_mean"], np.nan)
    out["shortfall035_to_income"] = np.where(valid_income, out["shortfall035"] / out["income_mean"], np.nan)
    out["liquid_to_reqdp20"] = out["liquid_nw_mean"] / out["req_dp20_mean"]
    out["liquid_to_reqdp035"] = out["liquid_nw_mean"] / out["req_dp035_mean"]
    out["dp_sufficient20"] = out["dp_gap20"].le(0)
    out["dp_sufficient035"] = out["dp_gap035"].le(0)

    ratio = out["shortfall20_to_income"]
    out["shortfall20_bin"] = np.select(
        [
            out["dp_sufficient20"],
            ratio.gt(0) & ratio.le(0.10),
            ratio.gt(0.10) & ratio.le(0.25),
            ratio.gt(0.25),
        ],
        [
            "dp_sufficient",
            "shortfall_le_10pct_income",
            "shortfall_10_25pct_income",
            "shortfall_gt_25pct_income",
        ],
        default="missing_income_or_dp",
    )
    return out


def summarize(df: pd.DataFrame, horizon: int) -> dict[str, float] | None:
    df = df.loc[df["weight_entry"].notna() & df["weight_entry"].gt(0)].copy()
    if df.empty:
        return None
    w = df["weight_entry"].to_numpy(dtype=float)
    parent = df[f"parent_by_{horizon}"].to_numpy(dtype=float)
    liquid = df["liquid_nw_mean"].to_numpy(dtype=float)
    income = df["income_mean"].to_numpy(dtype=float)
    short20 = df["shortfall20_to_income"].to_numpy(dtype=float)
    short035 = df["shortfall035_to_income"].to_numpy(dtype=float)

    positive_short20 = df["shortfall20"].gt(0).to_numpy()
    positive_short035 = df["shortfall035"].gt(0).to_numpy()
    ratio20 = df["shortfall20_to_income"]

    return {
        "n_obs": int(len(df)),
        "weight_sum": float(np.sum(w)),
        f"parent_rate_by_{horizon}": weighted_mean(parent, w),
        "age_entry_mean": weighted_mean(df["age_entry"].to_numpy(dtype=float), w),
        "entry_years_mean": weighted_mean(df["n_entry_years"].to_numpy(dtype=float), w),
        "liquid_nw_mean_k": weighted_mean(liquid, w) / 1000.0,
        "liquid_nw_median_k": weighted_quantile(liquid, w, 0.50) / 1000.0,
        "income_mean_k": weighted_mean(income, w) / 1000.0,
        "income_median_k": weighted_quantile(income, w, 0.50) / 1000.0,
        "target_price_mean_k": weighted_mean(df["p_target_mean"].to_numpy(dtype=float), w) / 1000.0,
        "req_dp20_mean_k": weighted_mean(df["req_dp20_mean"].to_numpy(dtype=float), w) / 1000.0,
        "req_dp035_mean_k": weighted_mean(df["req_dp035_mean"].to_numpy(dtype=float), w) / 1000.0,
        "liquid_to_reqdp20_mean": weighted_mean(df["liquid_to_reqdp20"].to_numpy(dtype=float), w),
        "liquid_to_reqdp20_median": weighted_quantile(df["liquid_to_reqdp20"].to_numpy(dtype=float), w, 0.50),
        "dp_sufficient20_share": weighted_mean(df["dp_sufficient20"].to_numpy(dtype=float), w),
        "dp_sufficient035_share": weighted_mean(df["dp_sufficient035"].to_numpy(dtype=float), w),
        "shortfall20_share": weighted_mean(positive_short20.astype(float), w),
        "shortfall035_share": weighted_mean(positive_short035.astype(float), w),
        "shortfall20_to_income_mean_all": weighted_mean(short20, w),
        "shortfall035_to_income_mean_all": weighted_mean(short035, w),
        "shortfall20_to_income_mean_if_short": weighted_mean(
            df.loc[positive_short20, "shortfall20_to_income"].to_numpy(dtype=float),
            df.loc[positive_short20, "weight_entry"].to_numpy(dtype=float),
        ),
        "shortfall20_le10pct_income_share": weighted_mean((ratio20.gt(0) & ratio20.le(0.10)).to_numpy(dtype=float), w),
        "shortfall20_10to25pct_income_share": weighted_mean((ratio20.gt(0.10) & ratio20.le(0.25)).to_numpy(dtype=float), w),
        "shortfall20_gt25pct_income_share": weighted_mean(ratio20.gt(0.25).to_numpy(dtype=float), w),
    }


def summarize_bins(df: pd.DataFrame, horizon: int) -> list[dict[str, float]]:
    out = []
    total_w = float(df.loc[df["weight_entry"].gt(0), "weight_entry"].sum())
    for bin_name, block in df.groupby("shortfall20_bin", dropna=False):
        block = block.loc[block["weight_entry"].notna() & block["weight_entry"].gt(0)]
        if block.empty:
            continue
        w = block["weight_entry"].to_numpy(dtype=float)
        parent = block[f"parent_by_{horizon}"].to_numpy(dtype=float)
        out.append(
            {
                "horizon": horizon,
                "shortfall20_bin": str(bin_name),
                "n_obs": int(len(block)),
                "weight_sum": float(np.sum(w)),
                "mass_share": float(np.sum(w) / max(total_w, 1e-14)),
                f"parent_rate_by_{horizon}": weighted_mean(parent, w),
                "liquid_nw_mean_k": weighted_mean(block["liquid_nw_mean"].to_numpy(dtype=float), w) / 1000.0,
                "income_mean_k": weighted_mean(block["income_mean"].to_numpy(dtype=float), w) / 1000.0,
                "shortfall20_to_income_mean_all": weighted_mean(
                    block["shortfall20_to_income"].to_numpy(dtype=float), w
                ),
            }
        )
    return out


def compare_groups(parent_stats: dict[str, float], childless_stats: dict[str, float], horizon: int) -> list[dict[str, float]]:
    metrics = [
        "liquid_nw_mean_k",
        "liquid_nw_median_k",
        "income_mean_k",
        "liquid_to_reqdp20_mean",
        "dp_sufficient20_share",
        "shortfall20_share",
        "shortfall20_to_income_mean_all",
        "shortfall20_to_income_mean_if_short",
        "shortfall20_gt25pct_income_share",
    ]
    rows = []
    for metric in metrics:
        rows.append(
            {
                "horizon": horizon,
                "metric": metric,
                "future_parent": parent_stats.get(metric, np.nan),
                "remain_childless": childless_stats.get(metric, np.nan),
                "difference_parent_minus_childless": parent_stats.get(metric, np.nan)
                - childless_stats.get(metric, np.nan),
            }
        )
    return rows


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    return float(np.sum(x[ok] * w[ok]) / np.sum(w[ok]))


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    x = x[ok]
    w = w[ok]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cw = np.cumsum(w)
    return float(x[np.searchsorted(cw, q * cw[-1], side="left")])


if __name__ == "__main__":
    main()
