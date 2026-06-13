#!/usr/bin/env python3
"""Proof of concept for quality-adjusted token billing.

The script compares visible output-token use for the winning and losing answers
in public Chatbot Arena preference battles. It uses one common tokenizer for all
models because the object is relative answer length, not exact provider billing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import tiktoken
from datasets import load_dataset
from matplotlib.ticker import PercentFormatter


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / "token_billing_quality"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze whether lower-preferred Arena answers often use more output tokens."
    )
    parser.add_argument(
        "--dataset",
        default="lmarena-ai/arena-human-preference-55k",
        help="Hugging Face dataset name.",
    )
    parser.add_argument("--split", default="train", help="Dataset split to load.")
    parser.add_argument(
        "--encoding",
        default="cl100k_base",
        help="tiktoken encoding used for all prompts and responses.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional row cap for smoke tests.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for tables, figure, and interpretation.",
    )
    return parser.parse_args()


def text_from_arena_cell(value: Any) -> str:
    """Convert Arena JSON-string conversation cells to plain joined text."""
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return ""
        if stripped[0] in "[{":
            try:
                return text_from_arena_cell(json.loads(stripped))
            except json.JSONDecodeError:
                return value
        return value
    if isinstance(value, dict):
        for key in ("content", "text", "value", "message"):
            if key in value:
                return text_from_arena_cell(value[key])
        parts = [text_from_arena_cell(v) for v in value.values()]
        return "\n\n".join(part for part in parts if part)
    if isinstance(value, (list, tuple)):
        parts = [text_from_arena_cell(item) for item in value]
        return "\n\n".join(part for part in parts if part)
    return str(value)


def count_tokens(text: str, encoding: tiktoken.Encoding) -> int:
    if not text:
        return 0
    return len(encoding.encode(text, disallowed_special=()))


def prompt_hash(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def label_winner(df: pd.DataFrame) -> pd.Series:
    conditions = [
        df["winner_model_a"].eq(1),
        df["winner_model_b"].eq(1),
        df["winner_tie"].eq(1),
    ]
    return pd.Series(
        np.select(conditions, ["A", "B", "tie"], default="other"),
        index=df.index,
        name="winner",
    )


def fmt_int(x: float | int) -> str:
    if pd.isna(x):
        return ""
    return f"{int(round(float(x))):,}"


def fmt_float(x: float, digits: int = 3) -> str:
    if pd.isna(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_pct(x: float, digits: int = 1) -> str:
    if pd.isna(x):
        return ""
    return f"{100 * x:.{digits}f}\\%"


def fit_lpm(df_reg: pd.DataFrame, token_var: str) -> Any:
    formula = f"a_win ~ {token_var} + log_prompt_tokens + C(model_a) + C(model_b)"
    duplicate_prompt_rows = len(df_reg) - df_reg["prompt_hash"].nunique()
    if duplicate_prompt_rows > 0 and df_reg["prompt_hash"].nunique() > 1:
        return smf.ols(formula, data=df_reg).fit(
            cov_type="cluster",
            cov_kwds={"groups": df_reg["prompt_hash"]},
        )
    return smf.ols(formula, data=df_reg).fit(cov_type="HC1")


def write_paper_table(stats: dict[str, Any], output_path: Path) -> None:
    rows = [
        ("Non-tie Arena battles", fmt_int(stats["main_n"])),
        (r"$\Pr(T^\ell>T^w)$", fmt_pct(stats["adverse_share"])),
        (r"Mean $T^\ell-T^w$", fmt_int(stats["delta_mean"])),
        (r"Median $T^\ell-T^w$", fmt_int(stats["delta_median"])),
        (r"Median $T^\ell/T^w$", fmt_float(stats["ratio_median"], 2)),
        (r"90th pct. $T^\ell/T^w$", fmt_float(stats["ratio_p90"], 2)),
        ("Tie battles", fmt_int(stats["tie_n"])),
        ("Mean absolute token gap in ties", fmt_int(stats["tie_abs_gap_mean"])),
        (
            r"LPM $\beta$: $\log T^A-\log T^B$",
            f"{fmt_float(stats['lpm_beta'], 3)} ({fmt_float(stats['lpm_se'], 3)})",
        ),
        (
            r"Prompt controls and model FE",
            "Yes",
        ),
        (
            r"Short-prompt $\Pr(T^\ell>T^w)$",
            fmt_pct(stats["short_adverse_share"]),
        ),
        (
            r"Long-prompt $\Pr(T^\ell>T^w)$",
            fmt_pct(stats["long_adverse_share"]),
        ),
    ]

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("\\begin{tabular}{l r}\n")
        handle.write("\\toprule\n")
        handle.write("Statistic & Estimate \\\\\n")
        handle.write("\\midrule\n")
        for label, value in rows:
            handle.write(f"{label} & {value} \\\\\n")
        handle.write("\\bottomrule\n")
        handle.write("\\end{tabular}\n")


def make_figure(main: pd.DataFrame, output_dir: Path, adverse_share: float) -> None:
    delta = main["loser_minus_winner_tokens"].astype(float)
    lo, hi = delta.quantile([0.01, 0.99])
    lo = min(lo, -1.0)
    hi = max(hi, 1.0)
    display_delta = delta.clip(lower=lo, upper=hi)
    bins = np.linspace(lo, hi, 41)
    weights = np.ones(len(display_delta)) / len(display_delta)

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    counts, edges, patches = ax.hist(
        display_delta,
        bins=bins,
        weights=weights,
        edgecolor="white",
        linewidth=0.5,
    )
    for patch, left, right in zip(patches, edges[:-1], edges[1:]):
        midpoint = 0.5 * (left + right)
        patch.set_facecolor("#b44b4b" if midpoint > 0 else "#5f6f84")
    ax.axvline(0, color="#222222", linewidth=1.0)
    ax.axvspan(0, hi, color="#b44b4b", alpha=0.08)
    ax.text(
        0.98,
        0.93,
        rf"$\Pr(T^\ell>T^w)={100 * adverse_share:.1f}\%$",
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#dddddd"},
    )
    ax.set_xlabel("Loser tokens minus winner tokens (winsorized for display)")
    ax.set_ylabel("Share of non-tie battles")
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_title("Token premium for lower-preferred answers")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_dir / "loser_token_premium_hist.pdf")
    fig.savefig(output_dir / "loser_token_premium_hist.png", dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset = load_dataset(args.dataset, split=args.split)
    if args.max_rows is not None:
        dataset = dataset.select(range(min(args.max_rows, len(dataset))))
    df = dataset.to_pandas()

    required_cols = {
        "id",
        "model_a",
        "model_b",
        "prompt",
        "response_a",
        "response_b",
        "winner_model_a",
        "winner_model_b",
        "winner_tie",
    }
    missing_cols = sorted(required_cols - set(df.columns))
    if missing_cols:
        raise ValueError(f"Dataset schema missing required columns: {missing_cols}")

    encoding = tiktoken.get_encoding(args.encoding)
    for source, target in [
        ("prompt", "prompt_text"),
        ("response_a", "response_a_text"),
        ("response_b", "response_b_text"),
    ]:
        df[target] = df[source].map(text_from_arena_cell)

    df["prompt_tokens"] = df["prompt_text"].map(lambda x: count_tokens(x, encoding))
    df["response_a_tokens"] = df["response_a_text"].map(lambda x: count_tokens(x, encoding))
    df["response_b_tokens"] = df["response_b_text"].map(lambda x: count_tokens(x, encoding))
    df["prompt_hash"] = df["prompt_text"].map(prompt_hash)
    df["winner"] = label_winner(df)

    clean = df[
        df["response_a_tokens"].gt(0)
        & df["response_b_tokens"].gt(0)
        & df["prompt_tokens"].gt(0)
        & df["winner"].isin(["A", "B", "tie"])
    ].copy()

    main_sample = clean[clean["winner"].isin(["A", "B"])].copy()
    main_sample["winner_tokens"] = np.where(
        main_sample["winner"].eq("A"),
        main_sample["response_a_tokens"],
        main_sample["response_b_tokens"],
    )
    main_sample["loser_tokens"] = np.where(
        main_sample["winner"].eq("A"),
        main_sample["response_b_tokens"],
        main_sample["response_a_tokens"],
    )
    main_sample["loser_minus_winner_tokens"] = (
        main_sample["loser_tokens"] - main_sample["winner_tokens"]
    )
    main_sample["loser_winner_ratio"] = (
        main_sample["loser_tokens"] / main_sample["winner_tokens"]
    )
    main_sample["adverse_billing"] = main_sample["loser_tokens"].gt(
        main_sample["winner_tokens"]
    )

    df_reg = main_sample.copy()
    df_reg["a_win"] = df_reg["winner"].eq("A").astype(int)
    df_reg["log_token_adv"] = np.log(df_reg["response_a_tokens"]) - np.log(
        df_reg["response_b_tokens"]
    )
    p01, p99 = df_reg["log_token_adv"].quantile([0.01, 0.99])
    df_reg["log_token_adv_winsor"] = df_reg["log_token_adv"].clip(p01, p99)
    df_reg["log_prompt_tokens"] = np.log(df_reg["prompt_tokens"])

    lpm = fit_lpm(df_reg, "log_token_adv")
    lpm_w = fit_lpm(df_reg, "log_token_adv_winsor")

    ties = clean[clean["winner"].eq("tie")].copy()
    ties["tie_abs_token_gap"] = (
        ties["response_a_tokens"] - ties["response_b_tokens"]
    ).abs()
    ties["tie_larger_over_smaller"] = np.maximum(
        ties["response_a_tokens"], ties["response_b_tokens"]
    ) / np.minimum(ties["response_a_tokens"], ties["response_b_tokens"])

    prompt_median = main_sample["prompt_tokens"].median()
    short = main_sample[main_sample["prompt_tokens"].le(prompt_median)]
    long = main_sample[main_sample["prompt_tokens"].gt(prompt_median)]

    stats = {
        "dataset": args.dataset,
        "split": args.split,
        "input_rows": len(df),
        "clean_rows": len(clean),
        "main_n": len(main_sample),
        "tie_n": len(ties),
        "adverse_share": main_sample["adverse_billing"].mean(),
        "delta_mean": main_sample["loser_minus_winner_tokens"].mean(),
        "delta_median": main_sample["loser_minus_winner_tokens"].median(),
        "ratio_mean": main_sample["loser_winner_ratio"].mean(),
        "ratio_median": main_sample["loser_winner_ratio"].median(),
        "ratio_p90": main_sample["loser_winner_ratio"].quantile(0.90),
        "ratio_p95": main_sample["loser_winner_ratio"].quantile(0.95),
        "tie_abs_gap_mean": ties["tie_abs_token_gap"].mean(),
        "tie_abs_gap_median": ties["tie_abs_token_gap"].median(),
        "tie_ratio_median": ties["tie_larger_over_smaller"].median(),
        "prompt_median_tokens": prompt_median,
        "short_adverse_share": short["adverse_billing"].mean(),
        "long_adverse_share": long["adverse_billing"].mean(),
        "lpm_beta": lpm.params["log_token_adv"],
        "lpm_se": lpm.bse["log_token_adv"],
        "lpm_pvalue": lpm.pvalues["log_token_adv"],
        "lpm_r2": lpm.rsquared,
        "lpm_winsor_beta": lpm_w.params["log_token_adv_winsor"],
        "lpm_winsor_se": lpm_w.bse["log_token_adv_winsor"],
        "lpm_winsor_pvalue": lpm_w.pvalues["log_token_adv_winsor"],
        "regression_covariance": lpm.cov_type,
        "duplicate_prompt_rows": len(df_reg) - df_reg["prompt_hash"].nunique(),
        "model_a_count": df_reg["model_a"].nunique(),
        "model_b_count": df_reg["model_b"].nunique(),
    }

    summary_rows = [
        {"statistic": key, "value": value}
        for key, value in stats.items()
        if isinstance(value, (int, float, np.integer, np.floating, str))
    ]
    pd.DataFrame(summary_rows).to_csv(output_dir / "arena_summary_table.csv", index=False)

    regression_rows = pd.DataFrame(
        [
            {
                "specification": "LPM with prompt length and model A/B fixed effects",
                "coefficient": lpm.params["log_token_adv"],
                "std_error": lpm.bse["log_token_adv"],
                "p_value": lpm.pvalues["log_token_adv"],
                "r_squared": lpm.rsquared,
                "n": int(lpm.nobs),
                "covariance": lpm.cov_type,
            },
            {
                "specification": "LPM with winsorized log token advantage",
                "coefficient": lpm_w.params["log_token_adv_winsor"],
                "std_error": lpm_w.bse["log_token_adv_winsor"],
                "p_value": lpm_w.pvalues["log_token_adv_winsor"],
                "r_squared": lpm_w.rsquared,
                "n": int(lpm_w.nobs),
                "covariance": lpm_w.cov_type,
            },
        ]
    )
    regression_rows.to_csv(output_dir / "arena_regression_table.csv", index=False)

    write_paper_table(stats, output_dir / "paper_table.tex")
    make_figure(main_sample, output_dir, stats["adverse_share"])

    with (output_dir / "interpretation.md").open("w", encoding="utf-8") as handle:
        handle.write("# Interpretation\n\n")
        handle.write(
            f"In {stats['main_n']:,} non-tie Arena battles, the lower-preferred "
            f"answer uses more visible output tokens in {100 * stats['adverse_share']:.1f}% "
            "of matched tasks. The mean loser-minus-winner token difference is "
            f"{stats['delta_mean']:.1f}, while the median is {stats['delta_median']:.1f}.\n\n"
        )
        handle.write(
            "The fixed-effect LPM regresses an indicator that response A wins on "
            "A's log token advantage over response B, prompt length, and model A/B "
            f"fixed effects. The coefficient is {stats['lpm_beta']:.3f} "
            f"with standard error {stats['lpm_se']:.3f}. This is descriptive, not "
            "causal: the exercise establishes that visible billable usage and "
            "observed preference are not mechanically aligned in matched tasks.\n"
        )

    print(f"Wrote outputs to {output_dir}")
    print(f"Non-tie battles: {stats['main_n']:,}")
    print(f"Pr(loser uses more tokens): {100 * stats['adverse_share']:.1f}%")
    print(
        "LPM beta on log token advantage: "
        f"{stats['lpm_beta']:.3f} ({stats['lpm_se']:.3f})"
    )


if __name__ == "__main__":
    main()
