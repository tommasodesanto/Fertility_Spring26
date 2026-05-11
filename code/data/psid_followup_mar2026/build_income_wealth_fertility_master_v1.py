"""Build PSID income/wealth/fertility bin summaries from the audited sample.

The companion Stata file constructs the sample. This script constructs the
long-format table used by the distributional discipline report.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUTDIR = ROOT / "code" / "data" / "psid_followup_mar2026" / "output" / "income_wealth_fertility_master_v1"
SAMPLE = OUTDIR / "income_wealth_fertility_master_sample_v1.dta"
OUTCSV = OUTDIR / "income_wealth_fertility_master_v1.csv"


WEALTH_VARS = [
    "total_nw_k",
    "liquid_nw_k",
    "broad_nw_k",
    "home_equity_k",
    "income_k",
    "earnings_k",
    "total_nw_to_inc",
    "liquid_nw_to_inc",
    "broad_nw_to_inc",
    "home_equity_to_inc",
]

GROUPS = {
    "all": lambda df: np.ones(len(df), dtype=bool),
    "women": lambda df: df["sex_young"].eq(2),
    "men": lambda df: df["sex_young"].eq(1),
    "renters_2530": lambda df: df["renter_2530"].eq(1),
    "owners_2530": lambda df: df["owner_2530"].eq(1),
}

OUTCOMES = {
    "parent_by_30": lambda df: df["observed_by_30"].eq(1),
    "parent_by_35": lambda df: df["observed_by_35"].eq(1),
    "parent_by_40": lambda df: df["observed_by_40"].eq(1),
    "parent_by_45": lambda df: df["observed_by_45"].eq(1),
    "childless_by_35": lambda df: df["observed_by_35"].eq(1),
    "childless_by_45": lambda df: df["observed_by_45"].eq(1),
    "children_35_40": lambda df: df["children_35_40"].notna(),
    "children_43_50": lambda df: df["children_43_50"].notna(),
    "owner_by_35": lambda df: df["owner_by_35"].notna(),
}


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_stata(SAMPLE, convert_categoricals=False)
    rows = []
    for wvar in WEALTH_VARS:
        for group_name, group_mask_fn in GROUPS.items():
            group_mask = group_mask_fn(df)
            for outcome, outcome_mask_fn in OUTCOMES.items():
                outcome_mask = outcome_mask_fn(df)
                for bin_type, nbins in [("quintile", 5), ("decile", 10)]:
                    keep = (
                        group_mask
                        & outcome_mask
                        & df[wvar].notna()
                        & df[outcome].notna()
                        & df["weight_young"].notna()
                        & df["weight_young"].gt(0)
                    )
                    sample = df.loc[keep, [wvar, outcome, "weight_young"]].copy()
                    if sample.empty:
                        continue

                    x = sample[wvar].to_numpy(dtype=float)
                    wt = sample["weight_young"].to_numpy(dtype=float)
                    p1 = weighted_quantile(x, wt, 0.01)
                    p99 = weighted_quantile(x, wt, 0.99)
                    sample = sample.loc[sample[wvar].between(p1, p99)]
                    if sample.empty:
                        continue

                    x = sample[wvar].to_numpy(dtype=float)
                    wt = sample["weight_young"].to_numpy(dtype=float)
                    cuts = [weighted_quantile(x, wt, q / nbins) for q in range(1, nbins)]
                    sample["bin"] = np.searchsorted(cuts, x, side="right") + 1

                    for b in range(1, nbins + 1):
                        block = sample.loc[sample["bin"].eq(b)]
                        if block.empty:
                            continue
                        xb = block[wvar].to_numpy(dtype=float)
                        yb = block[outcome].to_numpy(dtype=float)
                        wb = block["weight_young"].to_numpy(dtype=float)
                        rows.append(
                            {
                                "wealth_concept": wvar,
                                "group": group_name,
                                "outcome": outcome,
                                "bin_type": bin_type,
                                "bin": b,
                                "n_obs": int(len(block)),
                                "weight_sum": float(np.sum(wb)),
                                "mean_outcome": weighted_mean(yb, wb),
                                "wealth_mean": weighted_mean(xb, wb),
                                "wealth_median": weighted_quantile(xb, wb, 0.50),
                            }
                        )

    out = pd.DataFrame(rows)
    out.sort_values(["wealth_concept", "group", "outcome", "bin_type", "bin"], inplace=True)
    out.to_csv(OUTCSV, index=False)
    print(f"Wrote {len(out):,} rows to {OUTCSV}")


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    return float(np.sum(x * w) / max(float(np.sum(w)), 1e-14))


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
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
