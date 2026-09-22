#!/usr/bin/env python3
"""Direct four-year gross-household-earnings validation from the PSID extract.

This is a diagnostic, not a production target builder.  It uses only the
existing narrow EARNINDRRC extract, keeps complete consecutive annual cells,
and estimates a no-fixed-type persistent AR(1) plus iid-transitory process at
the four-year frequency.  The four offset grids are separate non-overlapping
block definitions.  No interpolation is used for the post-1997 biennial era.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.io.stata import StataReader
from scipy.optimize import least_squares


PERIOD_YEARS = 4
ANNUAL_FIRST = 1984
ANNUAL_LAST = 1997
OFFSETS = (0, 1, 2, 3)
LAGS = (0, 1, 2)
BOOTSTRAP_SEED = 20260921


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def default_extract() -> Path:
    return (
        repo_root()
        / "code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727"
        / "psid_income_md_extract.dta"
    )


def default_output() -> Path:
    return (
        repo_root()
        / "output/model/native_financing_diagnostic_20260919/specification_followup"
        / "earnings_wealth_v1/data_validation"
    )


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    return float(np.sum(x[ok] * w[ok]) / np.sum(w[ok]))


def weighted_cov(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    x = x[ok]
    y = y[ok]
    w = w[ok]
    sw = np.sum(w)
    mx = np.sum(w * x) / sw
    my = np.sum(w * y) / sw
    return float(np.sum(w * (x - mx) * (y - my)) / sw)


def weighted_geometric_mean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        return float("nan")
    return float(np.exp(np.mean(np.log(values))))


def load_extract(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing existing EARNINDRRC extract: {path}")
    df = pd.read_stata(path, convert_categoricals=False)
    required = {"ID", "year", "IW", "AGEREP", "EARNINDRRC"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Extract lacks required columns: {missing}")
    out = df[["ID", "year", "IW", "AGEREP", "EARNINDRRC"]].copy()
    out["id"] = out["ID"].astype(str)
    out["year"] = pd.to_numeric(out["year"], errors="coerce").astype("Int64")
    out["age"] = pd.to_numeric(out["AGEREP"], errors="coerce")
    out["weight"] = pd.to_numeric(out["IW"], errors="coerce")
    out["earnings"] = pd.to_numeric(out["EARNINDRRC"], errors="coerce")
    out = out.dropna(subset=["id", "year", "age", "weight", "earnings"])
    out["year"] = out["year"].astype(int)
    if out.duplicated(["id", "year"]).any():
        raise ValueError("Duplicate person-year observations in the narrow extract")
    out = out[(out["year"] >= ANNUAL_FIRST) & (out["year"] <= ANNUAL_LAST)]
    out = out[(out["weight"] > 0) & (out["earnings"] > 0)]
    if out.empty:
        raise ValueError("No positive annual EARNINDRRC observations remain")
    return out.sort_values(["id", "year"]).reset_index(drop=True)


def build_blocks(annual: pd.DataFrame, offset: int) -> pd.DataFrame:
    """Build complete, non-overlapping four-year cells for one offset grid."""
    starts = range(ANNUAL_FIRST + offset, ANNUAL_LAST - 2, PERIOD_YEARS)
    rows: list[dict[str, Any]] = []
    by_id = {pid: g.set_index("year") for pid, g in annual.groupby("id", sort=False)}
    for pid, g in by_id.items():
        years = set(g.index.astype(int))
        for start in starts:
            block_years = list(range(start, start + PERIOD_YEARS))
            if not set(block_years).issubset(years):
                continue
            cell = g.loc[block_years]
            # ``start`` follows the offset grid, so starts differ by exactly four
            # years.  No adjacent overlapping cells are ever formed here.
            rows.append(
                {
                    "id": pid,
                    "block_start": int(start),
                    "block_end": int(start + PERIOD_YEARS - 1),
                    "age_start": float(cell["age"].iloc[0]),
                    "earnings_total": float(cell["earnings"].sum()),
                    "earnings_mean": float(cell["earnings"].mean()),
                    "block_weight": weighted_geometric_mean(cell["weight"].to_numpy()),
                    "n_annual_obs": int(len(cell)),
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "id", "block_start", "block_end", "age_start", "earnings_total",
                "earnings_mean", "block_weight", "n_annual_obs",
            ]
        )
    blocks = pd.DataFrame(rows)
    blocks["log_earnings_mean"] = np.log(blocks["earnings_mean"])
    blocks["log_earnings_total"] = np.log(blocks["earnings_total"])
    return blocks.sort_values(["id", "block_start"]).reset_index(drop=True)


def weighted_fe_residuals(blocks: pd.DataFrame, frequencies: dict[str, int] | None) -> np.ndarray:
    """Residualize block log income on block-start age and year FE."""
    if blocks.empty:
        return np.array([], dtype=float)
    freq = blocks["id"].map(frequencies).fillna(0.0).to_numpy(float) if frequencies is not None else np.ones(len(blocks))
    base_w = blocks["block_weight"].to_numpy(float)
    row_w = base_w * freq
    y = blocks["log_earnings_mean"].to_numpy(float)
    ages = blocks["age_start"].round().astype(int).astype(str)
    years = blocks["block_start"].astype(int).astype(str)
    X = pd.get_dummies(pd.DataFrame({"age": ages, "year": years}), drop_first=False, dtype=float).to_numpy()
    keep = np.isfinite(y) & np.isfinite(row_w) & (row_w > 0)
    if not np.any(keep):
        return np.full(len(blocks), np.nan)
    sw = np.sqrt(row_w[keep])
    coef, *_ = np.linalg.lstsq(X[keep] * sw[:, None], y[keep] * sw, rcond=None)
    residual = np.full(len(blocks), np.nan)
    residual[keep] = y[keep] - X[keep] @ coef
    return residual


def empirical_moments(blocks: pd.DataFrame, residual: np.ndarray, frequencies: dict[str, int] | None) -> dict[str, Any]:
    if blocks.empty:
        return {"moments": {str(k): float("nan") for k in LAGS}, "counts": {str(k): 0 for k in LAGS}, "persons": {str(k): 0 for k in LAGS}}
    b = blocks.copy()
    b["residual"] = residual
    freq = b["id"].map(frequencies).fillna(0.0).to_numpy(float) if frequencies is not None else np.ones(len(b))
    b["analysis_weight"] = b["block_weight"].to_numpy(float) * freq
    moments = {"0": weighted_cov(b["residual"].to_numpy(), b["residual"].to_numpy(), b["analysis_weight"].to_numpy())}
    counts = {"0": int(np.sum(freq > 0))}
    persons = {"0": int(b.loc[freq > 0, "id"].nunique())}
    for lag in (1, 2):
        left = b[["id", "block_start", "residual", "block_weight"]].copy()
        right = b[["id", "block_start", "residual", "block_weight"]].copy()
        left["join_start"] = left["block_start"] + lag * PERIOD_YEARS
        right = right.rename(columns={"block_start": "join_start", "residual": "residual_lag", "block_weight": "weight_lag"})
        pairs = left.merge(right, on=["id", "join_start"], how="inner")
        if pairs.empty:
            moments[str(lag)] = float("nan")
            counts[str(lag)] = 0
            persons[str(lag)] = 0
            continue
        pair_freq = pairs["id"].map(frequencies).fillna(0.0).to_numpy(float) if frequencies is not None else np.ones(len(pairs))
        pair_w = np.sqrt(pairs["block_weight"].to_numpy(float) * pairs["weight_lag"].to_numpy(float)) * pair_freq
        moments[str(lag)] = weighted_cov(pairs["residual"].to_numpy(), pairs["residual_lag"].to_numpy(), pair_w)
        counts[str(lag)] = int(np.sum(pair_freq > 0))
        persons[str(lag)] = int(pairs.loc[pair_freq > 0, "id"].nunique())
    return {"moments": moments, "counts": counts, "persons": persons}


def fit_process(moments: dict[str, float], counts: dict[str, int]) -> dict[str, Any]:
    """Fit (rho, stationary persistent variance, iid variance), all nonnegative."""
    g0, g1, g2 = (float(moments[str(k)]) for k in LAGS)
    n = np.asarray([counts[str(k)] for k in LAGS], dtype=float)
    informative = np.isfinite([g0, g1, g2]) & (n > 0)
    if int(np.sum(informative)) < 3:
        return {"status": "underidentified_or_missing_moment", "n_informative_moments": int(np.sum(informative))}
    unconstrained_rho = g2 / g1 if np.isfinite(g1) and g1 != 0 else float("nan")
    unconstrained_sp = (g1 * g1 / g2) if np.isfinite(g2) and g2 > 0 else float("nan")
    unconstrained_e = g0 - unconstrained_sp if np.isfinite(unconstrained_sp) else float("nan")
    # Pair-count weighting follows the existing PSID builder. Bounds are explicit:
    # rho in [1e-8, 1-1e-8], both variances in [0, +inf).
    scale = np.sqrt(n / np.sum(n))
    def residual(par: np.ndarray) -> np.ndarray:
        rho, sp, se = par
        pred = np.array([sp + se, sp * rho, sp * rho * rho])
        return scale * (pred - np.array([g0, g1, g2]))
    starts = [
        np.array([
            float(np.clip(unconstrained_rho if np.isfinite(unconstrained_rho) else 0.8, 1e-6, 1 - 1e-6)),
            max(float(unconstrained_sp) if np.isfinite(unconstrained_sp) else max(g0, 1e-6), 1e-8),
            max(float(unconstrained_e) if np.isfinite(unconstrained_e) else 0.0, 1e-8),
        ]),
        np.array([0.8, max(g0 * 0.7, 1e-6), max(g0 * 0.3, 1e-6)]),
        np.array([0.95, max(g0 * 0.9, 1e-6), max(g0 * 0.1, 1e-6)]),
    ]
    fits = [least_squares(residual, x0, bounds=([1e-8, 0.0, 0.0], [1.0 - 1e-8, np.inf, np.inf]), xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=10000) for x0 in starts]
    fit = min(fits, key=lambda x: float(np.sum(x.fun * x.fun)))
    rho, sp, se = map(float, fit.x)
    pred = np.array([sp + se, sp * rho, sp * rho * rho])
    return {
        "status": "ok",
        "rho_4yr": rho,
        "persistent_variance_4yr": sp,
        "iid_variance_4yr": se,
        "iid_sd_4yr": math.sqrt(se),
        "objective_pair_count_weighted": float(np.sum(fit.fun * fit.fun)),
        "unconstrained_rho_4yr": float(unconstrained_rho),
        "unconstrained_persistent_variance_4yr": float(unconstrained_sp),
        "unconstrained_iid_variance_4yr": float(unconstrained_e),
        "unconstrained_iid_negative": bool(np.isfinite(unconstrained_e) and unconstrained_e < 0),
        "iid_at_lower_bound": bool(se <= 1e-10),
        "rho_at_bound": bool(rho <= 1e-7 or rho >= 1.0 - 1e-7),
        "persistent_variance_at_bound": bool(sp <= 1e-10),
        "predicted_covariance": {str(k): float(v) for k, v in zip(LAGS, pred)},
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "parameter_bounds": {"rho_4yr": [1e-8, 1.0 - 1e-8], "persistent_variance_4yr": [0.0, "inf"], "iid_variance_4yr": [0.0, "inf"]},
    }


def one_offset(
    annual: pd.DataFrame,
    offset: int,
    frequencies: dict[str, int] | None = None,
    blocks_override: pd.DataFrame | None = None,
) -> dict[str, Any]:
    blocks = blocks_override.copy() if blocks_override is not None else build_blocks(annual, offset)
    residual = weighted_fe_residuals(blocks, frequencies)
    emp = empirical_moments(blocks, residual, frequencies)
    fit = fit_process(emp["moments"], emp["counts"])
    age_variance_rows = []
    if not blocks.empty:
        age = blocks["age_start"].round().astype(int)
        bins = pd.cut(age, bins=[24, 34, 44, 54, 60], labels=["25-34", "35-44", "45-54", "55-60"])
        freq = blocks["id"].map(frequencies).fillna(0.0).to_numpy(float) if frequencies is not None else np.ones(len(blocks))
        for label, idx in bins.groupby(bins, observed=False).groups.items():
            ii = np.asarray(list(idx), dtype=int)
            w = blocks["block_weight"].to_numpy(float)[ii] * freq[ii]
            x = residual[ii]
            age_variance_rows.append({
                "offset": int(offset),
                "age_bin": str(label),
                "blocks": int(len(ii)),
                "persons": int(blocks.iloc[ii]["id"].nunique()),
                "weighted_mean_residual": weighted_mean(x, w),
                "weighted_residual_log_variance": weighted_cov(x, x, w),
            })
    support = {
        "offset": int(offset),
        "annual_year_min": ANNUAL_FIRST,
        "annual_year_max": ANNUAL_LAST,
        "persons_with_blocks": int(blocks["id"].nunique()) if not blocks.empty else 0,
        "blocks": int(len(blocks)),
        "block_start_min": int(blocks["block_start"].min()) if not blocks.empty else None,
        "block_start_max": int(blocks["block_start"].max()) if not blocks.empty else None,
        "complete_annual_blocks": True,
        "nonoverlap_stride_years": PERIOD_YEARS,
        "pair_counts": emp["counts"],
        "pair_persons": emp["persons"],
    }
    moment_rows = []
    pred = fit.get("predicted_covariance", {})
    for lag in LAGS:
        moment_rows.append({
            "offset": int(offset),
            "lag_4yr": int(lag),
            "empirical_covariance": float(emp["moments"][str(lag)]),
            "fitted_covariance": float(pred.get(str(lag), np.nan)),
            "gap": float(emp["moments"][str(lag)] - pred.get(str(lag), np.nan)) if str(lag) in pred else np.nan,
            "pair_count": int(emp["counts"][str(lag)]),
            "pair_person_count": int(emp["persons"][str(lag)]),
        })
    return {"offset": offset, "support": support, "moments": moment_rows, "fit": fit, "age_variance": age_variance_rows, "blocks": blocks}


def bootstrap(annual: pd.DataFrame, point: list[dict[str, Any]], reps: int, seed: int) -> pd.DataFrame:
    if reps <= 0:
        return pd.DataFrame()
    ids = annual["id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for rep in range(reps):
        draw = rng.choice(ids, size=len(ids), replace=True)
        counts = pd.Series(draw).value_counts().to_dict()
        for rec in point:
            fit = one_offset(
                annual,
                int(rec["offset"]),
                {str(k): int(v) for k, v in counts.items()},
                blocks_override=rec["blocks"],
            )["fit"]
            rows.append({"rep": rep, "offset": int(rec["offset"]), **{k: v for k, v in fit.items() if k in {"rho_4yr", "persistent_variance_4yr", "iid_variance_4yr", "objective_pair_count_weighted", "iid_at_lower_bound", "unconstrained_iid_negative"}}})
    return pd.DataFrame(rows)


def bootstrap_weight_test(point: list[dict[str, Any]]) -> None:
    """Fail closed if omitted persons receive nonzero bootstrap weight."""
    rec = next(r for r in point if r["support"]["pair_counts"]["2"] > 0)
    blocks = rec["blocks"]
    counts_by_person = blocks.groupby("id").size().sort_values(ascending=False)
    selected = str(counts_by_person.index[0])
    frequencies = {selected: 1}
    residual = weighted_fe_residuals(blocks, frequencies)
    omitted = blocks["id"].to_numpy() != selected
    if np.any(np.isfinite(residual[omitted])):
        raise AssertionError("bootstrap omitted-person weight test failed: omitted rows have residuals")
    emp = empirical_moments(blocks, residual, frequencies)
    expected_lag_counts = {}
    expected_lag_persons = {}
    for lag in (1, 2):
        left = blocks[["id", "block_start"]].copy()
        right = blocks[["id", "block_start"]].copy()
        left["join_start"] = left["block_start"] + lag * PERIOD_YEARS
        right = right.rename(columns={"block_start": "join_start"})
        pairs = left.merge(right, on=["id", "join_start"], how="inner")
        selected_pairs = pairs[pairs["id"] == selected]
        expected_lag_counts[str(lag)] = int(len(selected_pairs))
        expected_lag_persons[str(lag)] = int(selected_pairs["id"].nunique())
    if emp["counts"]["0"] != int(counts_by_person.iloc[0]):
        raise AssertionError("bootstrap omitted-person weight test failed: support count mismatch")
    for lag in (1, 2):
        if emp["counts"][str(lag)] != expected_lag_counts[str(lag)] or emp["persons"][str(lag)] != expected_lag_persons[str(lag)]:
            raise AssertionError(f"bootstrap omitted-person weight test failed: lag-{lag} pair support mismatch")


def validate_block_construction(point: list[dict[str, Any]]) -> dict[str, Any]:
    checks = []
    for rec in point:
        blocks = rec["blocks"]
        bad = 0
        for _, group in blocks.groupby("id"):
            starts = np.sort(group["block_start"].to_numpy(dtype=int))
            if starts.size and np.any(np.diff(starts) < PERIOD_YEARS):
                bad += 1
        checks.append({"offset": int(rec["offset"]), "blocks": int(len(blocks)), "persons": int(blocks["id"].nunique()), "overlap_violations": int(bad)})
    return {"passed": all(x["overlap_violations"] == 0 for x in checks), "checks": checks}


def synthetic_observer_recovery(point: list[dict[str, Any]], reps: int = 100) -> dict[str, Any]:
    """Apply the same observer to a known process on the empirical block shape."""
    rec = next(r for r in point if r["support"]["pair_counts"]["2"] > 0)
    blocks = rec["blocks"].copy()
    rho_true = 0.80
    sp_true = 0.50
    se_true = 0.04
    truth = {"rho_4yr": rho_true, "persistent_variance_4yr": sp_true, "iid_variance_4yr": se_true}
    analytic_cov = {"0": sp_true + se_true, "1": sp_true * rho_true, "2": sp_true * rho_true**2}
    fit_rows = []
    cov_rows = []
    for rep in range(int(reps)):
        rng = np.random.default_rng(BOOTSTRAP_SEED + 1 + rep)
        values = np.full(len(blocks), np.nan)
        for _, idx in blocks.sort_values(["id", "block_start"]).groupby("id", sort=False).groups.items():
            ii = np.asarray(list(idx), dtype=int)
            persistent = rng.normal(0.0, math.sqrt(sp_true))
            prev_start = None
            for pos in ii:
                current_start = int(blocks.loc[pos, "block_start"])
                if prev_start is not None:
                    gap_periods = (current_start - prev_start) // PERIOD_YEARS
                    if gap_periods <= 0:
                        raise AssertionError("synthetic block starts are not strictly increasing")
                    rho_gap = rho_true**gap_periods
                    persistent = rho_gap * persistent + rng.normal(0.0, math.sqrt(sp_true * (1.0 - rho_gap**2)))
                values[pos] = persistent + rng.normal(0.0, math.sqrt(se_true))
                prev_start = current_start
        synthetic_blocks = blocks.copy()
        synthetic_blocks["log_earnings_mean"] = values
        residual = weighted_fe_residuals(synthetic_blocks, None)
        emp = empirical_moments(synthetic_blocks, residual, None)
        fit = fit_process(emp["moments"], emp["counts"])
        fit_rows.append({"rep": rep, "rho_4yr": fit.get("rho_4yr"), "persistent_variance_4yr": fit.get("persistent_variance_4yr"), "iid_variance_4yr": fit.get("iid_variance_4yr")})
        cov_rows.extend({"rep": rep, "lag_4yr": int(lag), "observed_covariance": float(emp["moments"][str(lag)])} for lag in LAGS)
    fit_df = pd.DataFrame(fit_rows)
    cov_df = pd.DataFrame(cov_rows)
    parameter_summary = {}
    for name, true in truth.items():
        x = fit_df[name].dropna().to_numpy(float)
        parameter_summary[name] = {"true": true, "reps": int(len(x)), "mean": float(np.mean(x)), "mean_bias": float(np.mean(x) - true), "monte_carlo_se_of_mean": float(np.std(x, ddof=1) / math.sqrt(len(x))), "p025": float(np.quantile(x, .025)), "median": float(np.quantile(x, .5)), "p975": float(np.quantile(x, .975))}
    covariance_summary = {}
    for lag in LAGS:
        x = cov_df.loc[cov_df["lag_4yr"] == lag, "observed_covariance"].dropna().to_numpy(float)
        truth_cov = analytic_cov[str(lag)]
        covariance_summary[str(lag)] = {"analytic_truth": truth_cov, "reps": int(len(x)), "observed_average_covariance_after_FE": float(np.mean(x)), "mean_bias_vs_analytic_truth": float(np.mean(x) - truth_cov), "monte_carlo_se_of_mean": float(np.std(x, ddof=1) / math.sqrt(len(x))), "p025": float(np.quantile(x, .025)), "median": float(np.quantile(x, .5)), "p975": float(np.quantile(x, .975)), "finite_FE_projection_gap": float(np.mean(x) - truth_cov)}
    return {
        "status": "same_observer_monte_carlo_executed",
        "offset_shape": int(rec["offset"]),
        "replications": int(reps),
        "truth": truth,
        "parameter_summary": parameter_summary,
        "covariance_summary": covariance_summary,
        "support": rec["support"],
        "note": "Each synthetic path uses the exact observed block starts for each person. An observed gap of g four-year periods advances persistence with rho^g and innovation variance Vp*(1-rho^(2g)); the reported projection gap is the finite-sample effect of the same age/year FE residualization and missingness pattern.",
    }


def write_outputs(out: Path, annual: pd.DataFrame, point: list[dict[str, Any]], boot: pd.DataFrame, extract: Path, bootstrap_reps: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    support = pd.DataFrame([r["support"] for r in point])
    # Keep nested count fields machine-readable and also provide long tables.
    support.to_json(out / "support.json", orient="records", indent=2)
    support_long = []
    for r in point:
        for lag, count in r["support"]["pair_counts"].items():
            support_long.append({"offset": r["offset"], "lag_4yr": int(lag), "pair_count": count, "pair_person_count": r["support"]["pair_persons"][lag]})
    pd.DataFrame(support_long).to_csv(out / "support_pairs.csv", index=False)
    pd.DataFrame([m for r in point for m in r["moments"]]).to_csv(out / "block_covariances.csv", index=False)
    pd.DataFrame([a for r in point for a in r["age_variance"]]).to_csv(out / "residual_variance_by_age_bin.csv", index=False)
    fit_rows = [{"offset": r["offset"], **r["fit"]} for r in point]
    pd.json_normalize(fit_rows).to_csv(out / "fit.csv", index=False)
    if not boot.empty:
        boot.to_csv(out / "bootstrap_draws.csv", index=False)
        summary = []
        for offset, g in boot.groupby("offset"):
            for p in ("rho_4yr", "persistent_variance_4yr", "iid_variance_4yr"):
                x = pd.to_numeric(g[p], errors="coerce").dropna().to_numpy()
                if x.size == 0:
                    summary.append({"offset": int(offset), "parameter": p, "bootstrap_reps": 0, "p025": np.nan, "median": np.nan, "p975": np.nan, "mean": np.nan, "sd": np.nan})
                    continue
                summary.append({"offset": int(offset), "parameter": p, "bootstrap_reps": int(len(x)), "p025": float(np.quantile(x, .025)), "median": float(np.quantile(x, .5)), "p975": float(np.quantile(x, .975)), "mean": float(np.mean(x)), "sd": float(np.std(x, ddof=1))})
        pd.DataFrame(summary).to_csv(out / "bootstrap_summary.csv", index=False)
    existing = repo_root() / "code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/md_autocovariance_fit.csv"
    candidate = repo_root() / "output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json"
    compare: dict[str, Any] = {"existing_annual_fixed_effect_fit_path": str(existing), "existing_annual_candidate_path": str(candidate)}
    if existing.exists():
        compare["existing_annual_fixed_effect_fit"] = pd.read_csv(existing).to_dict(orient="records")
    if candidate.exists():
        compare["existing_annual_no_fixed_candidate"] = json.loads(candidate.read_text())
    (out / "existing_annual_comparison.json").write_text(json.dumps(compare, indent=2, sort_keys=True))
    block_validation = validate_block_construction(point)
    (out / "block_construction_validation.json").write_text(json.dumps(block_validation, indent=2, sort_keys=True))
    synthetic = synthetic_observer_recovery(point)
    (out / "synthetic_observer_recovery.json").write_text(json.dumps(synthetic, indent=2, sort_keys=True))
    source_labels = StataReader(str(extract)).variable_labels()
    metadata = {
        "status": "diagnostic_only_not_adopted",
        "extract": str(extract),
        "extract_sha256": __import__("hashlib").sha256(extract.read_bytes()).hexdigest(),
        "source_variable_labels": source_labels,
        "timing_contract": {
            "year": source_labels.get("year", ""),
            "EARNINDRRC": source_labels.get("EARNINDRRC", ""),
            "interpretation": "retain survey-year label for year; source labels EARNINDRRC as tax-year earnings. No survey-to-tax-year shift was applied; a common shift would not change within-person covariance but remains a reporting limitation.",
        },
        "earnings_variable": "EARNINDRRC",
        "concept": "RP/spouse combined real gross labor earnings (2022 dollars per source builder)",
        "annual_era_only": [ANNUAL_FIRST, ANNUAL_LAST],
        "post_1997_biennial_years_excluded_from_blocks": [1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019],
        "period_definition": "sum of four observed annual gross earnings; log mean used for covariance because log(sum)=log(mean)+constant",
        "weight_definition": "geometric mean of four annual IW weights for each block; geometric mean of block weights for covariance pairs",
        "residualization": "weighted least squares of log four-year mean on block-start integer age and block-start calendar-year fixed effects",
        "covariance_lags": [0, 1, 2],
        "nonoverlap_stride_years": 4,
        "fit": "stationary persistent AR(1) plus iid transitory at four-year frequency; rho in [1e-8,1-1e-8], both variances >=0",
        "bootstrap_reps": int(bootstrap_reps),
        "bootstrap_seed": int(BOOTSTRAP_SEED),
        "bootstrap_weight_test": {"passed": True, "description": "omitted persons receive frequency zero; selected-person support count matches"},
        "block_construction_validation": block_validation,
        "synthetic_observer_recovery": synthetic,
        "n_input_rows_after_filters": int(len(annual)),
        "n_input_persons_after_filters": int(annual["id"].nunique()),
        "entrant_ages_18_24_represented": bool(((annual["age"] >= 18) & (annual["age"] <= 24)).any()),
        "input_age_min": int(annual["age"].min()),
        "input_age_max": int(annual["age"].max()),
    }
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))
    write_readme(out, metadata, point, boot)


def write_readme(out: Path, metadata: dict[str, Any], point: list[dict[str, Any]], boot: pd.DataFrame) -> None:
    lines = [
        "# Direct four-year gross-household-earnings validation",
        "",
        "Diagnostic only; no target or model specification is adopted.",
        "",
        f"Input: `{metadata['extract']}`; SHA-256 `{metadata['extract_sha256']}`.",
        f"The input has {metadata['n_input_rows_after_filters']:,} positive annual person-years for {metadata['n_input_persons_after_filters']:,} persons after the existing extract filters.",
        "",
        "## Definition",
        "",
        "A block is the arithmetic sum of four observed consecutive annual `EARNINDRRC` values; the logged block mean is used for covariances because dividing every block by four changes only the log intercept. Only annual-era years 1984--1997 are eligible. Post-1997 biennial observations are never interpolated or treated as zeros.",
        "",
        "For each offset `o = 0,1,2,3`, starts are `1984 + o + 4k`; therefore blocks do not overlap within an offset. Block weights are the geometric mean of the four annual `IW` values. Log block means are residualized by weighted least squares on block-start integer age and block-start calendar-year fixed effects. Covariances use geometric-mean block weights and only adjacent non-overlapping blocks at four-year lags 1 and 2.",
        "",
        "The three moments `(gamma_0, gamma_1, gamma_2)` identify the three nonnegative parameters `(persistent variance, rho, iid variance)` absent sampling noise. The optimizer imposes `rho in [1e-8, 1-1e-8]` and both variances `>= 0`; unconstrained implied iid variance and boundary hits are reported rather than clipped silently.",
        "",
        "## Point estimates and support",
        "",
        "| offset | persons | blocks | lag-1 pairs | lag-2 pairs | rho(4yr) | persistent var | iid var | iid boundary | unconstrained iid var |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|",
    ]
    for r in point:
        s = r["support"]
        f = r["fit"]
        lines.append(
            f"| {r['offset']} | {s['persons_with_blocks']} | {s['blocks']} | {s['pair_counts']['1']} | {s['pair_counts']['2']} | {f.get('rho_4yr', float('nan')):.6g} | {f.get('persistent_variance_4yr', float('nan')):.6g} | {f.get('iid_variance_4yr', float('nan')):.6g} | {f.get('iid_at_lower_bound', '')} | {f.get('unconstrained_iid_variance_4yr', float('nan')):.6g} |"
        )
    annual_candidate_path = repo_root() / "output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json"
    if annual_candidate_path.exists():
        annual_candidate = json.loads(annual_candidate_path.read_text())
        annual_params = annual_candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
        lines += [
            "",
            "## Frequency and concept comparison",
            "",
            "The existing no-fixed-type annual candidate reports `rho_annual = "
            f"{annual_params['rho_annual']:.6f}`, persistent variance "
            f"`{annual_params['persistent_variance']:.6f}`, and transitory variance "
            f"`{annual_params['transitory_variance']:.6f}`. The direct estimates above are four-year-frequency covariances from complete annual cells, so their `rho_4yr` and variances are not numerically comparable without an explicit aggregation map. The annual candidate also uses the same EARNINDRRC concept but a different residualization/moment schedule and an endpoint-plus-iid period approximation; this packet does not force annual AR(1) equivalence.",
            "",
            "The existing fixed-effect annual packet's fitted annual `rho` is 0.886345, with fixed-effect variance 0.393053, persistent variance 0.331897, and transitory variance 0.309752. That fixed-effect decomposition is a separate annual observer and is not imposed in this direct four-year no-fixed-type fit.",
        ]
    lines += [
        "",
        "All lag-0/1/2 empirical and fitted covariances, pair counts, and pair-person counts are in `block_covariances.csv`; the full optimizer receipt is in `fit.csv`. Existing annual fixed-effect and no-fixed candidate source values are preserved in `existing_annual_comparison.json` for concept-mismatch review.",
        "",
        "`residual_variance_by_age_bin.csv` reports the weighted residual log variance by block-start age bins 25--34, 35--44, 45--54, and 55--60. The source extract contains no ages 18--24, so it cannot validate the model's entrant-age distribution; that is an external entry restriction rather than evidence of zero entrant risk.",
        "",
        "The Stata variable labels identify `year` as survey year and `EARNINDRRC` as tax-year earnings. The diagnostic retains the survey-year labels and applies no unverified timing shift; a common shift would leave within-person covariance unchanged, but the survey/tax-year distinction limits level and age-profile interpretation.",
        "",
        "`synthetic_observer_recovery.json` applies the same block construction, missingness shape, age/year FE residualization, weights, and covariance estimator to a known four-year AR(1)+iid process. Its finite-sample bias is reported rather than treated as a pass/fail calibration result. `block_construction_validation.json` confirms no overlapping blocks within any offset.",
        "",
        "The deterministic bootstrap-weight test is recorded in `run_metadata.json`; it verifies that omitted persons receive zero frequency rather than the point-estimate default weight.",
        "",
        "## Interpretation limits",
        "",
        "The four-year process is estimated on complete annual cells only, so support is short and selected toward survivors with positive reported labor earnings. It is not an estimate of the post-1997 biennial population. `EARNINDRRC` is gross labor earnings, while the model budget is after-payroll-tax period household resources; the model's tax, pension, transfers, and entry mapping remain separate objects. This direct process must therefore be compared to the existing annual source as a measurement diagnostic, not substituted mechanically.",
        "",
        f"Person bootstrap draws requested: {metadata['bootstrap_reps']}. Bootstrap files are present only when the run is invoked with a positive `--bootstrap-reps`.",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", type=Path, default=default_extract())
    parser.add_argument("--output", type=Path, default=default_output())
    parser.add_argument("--bootstrap-reps", type=int, default=0)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    annual = load_extract(args.extract)
    point = [one_offset(annual, offset) for offset in OFFSETS]
    bootstrap_weight_test(point)
    boot = bootstrap(annual, point, max(0, args.bootstrap_reps), args.bootstrap_seed)
    write_outputs(args.output, annual, point, boot, args.extract, max(0, args.bootstrap_reps))
    for r in point:
        f = r["fit"]
        s = r["support"]
        print(json.dumps({"offset": r["offset"], "persons": s["persons_with_blocks"], "blocks": s["blocks"], "pairs_lag1": s["pair_counts"]["1"], "pairs_lag2": s["pair_counts"]["2"], "rho_4yr": f.get("rho_4yr"), "persistent_variance_4yr": f.get("persistent_variance_4yr"), "iid_variance_4yr": f.get("iid_variance_4yr"), "iid_at_lower_bound": f.get("iid_at_lower_bound"), "unconstrained_iid_variance_4yr": f.get("unconstrained_iid_variance_4yr")}, sort_keys=True))


if __name__ == "__main__":
    main()
