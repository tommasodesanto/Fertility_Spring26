#!/usr/bin/env python3
"""Audit whether a steeper fertility-price response restores replacement.

This is an uncalibrated sensitivity around the sequential child-room-floor
estimate.  For each common multiplier on the first-birth and continuation
fertility-logit scales, the script reanchors the preference intercept at the
old stationary housing cost and then traces effective demographic renewal over
housing costs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import audit_closed_reproductive_closure as closure


PRIMARY_RENEWAL_FIELD = "topcode_consistent_B_over_E"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=closure.E5F_FLOOR_SOURCE)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "output/model/sequential_fertility_price_feedback_audit",
    )
    parser.add_argument(
        "--kappa-factors",
        default="1.0,0.75,0.5,0.33,0.25,0.10,0.05,0.02",
        help="Comma-separated common multipliers on both fertility-logit scales.",
    )
    parser.add_argument("--old-fertility-target", type=float, default=2.12)
    parser.add_argument("--low-fertility-target", type=float, default=1.6165)
    parser.add_argument(
        "--old-psi-reference",
        type=float,
        default=0.1062,
        help="Preference intercept used to recover the old equilibrium house price.",
    )
    parser.add_argument("--nb", type=int, default=120)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def bisection(
    evaluate: Callable[[float], float],
    *,
    lower: float,
    upper: float,
    target: float,
    increasing: bool,
    value_tolerance: float,
    max_iter: int = 18,
) -> float:
    f_lower = float(evaluate(lower)) - target
    f_upper = float(evaluate(upper)) - target
    if not increasing:
        f_lower, f_upper = -f_lower, -f_upper
    if f_lower > 0.0 or f_upper < 0.0:
        raise RuntimeError(
            "Root is not bracketed: "
            f"lower_gap={f_lower:.8g}, upper_gap={f_upper:.8g}"
        )
    midpoint = 0.5 * (lower + upper)
    for _ in range(max_iter):
        midpoint = 0.5 * (lower + upper)
        raw_gap = float(evaluate(midpoint)) - target
        if abs(raw_gap) <= value_tolerance:
            return midpoint
        gap = raw_gap if increasing else -raw_gap
        if gap < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return midpoint


def make_plot(price_rows: list[dict[str, Any]], outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    factors = sorted({float(row["kappa_factor"]) for row in price_rows}, reverse=True)
    for factor in factors:
        rows = sorted(
            [row for row in price_rows if float(row["kappa_factor"]) == factor],
            key=lambda row: float(row["price_ratio"]),
        )
        x = np.asarray([float(row["price_ratio"]) for row in rows])
        renewal = np.asarray([float(row[PRIMARY_RENEWAL_FIELD]) for row in rows])
        fertility = np.asarray([float(row["tfr_topcoded"]) for row in rows])
        axes[0].plot(x, renewal, marker="o", lw=1.8, label=f"scale x {factor:g}")
        axes[1].plot(x, fertility, marker="o", lw=1.8, label=f"scale x {factor:g}")
    axes[0].axhline(1.0, color="black", lw=0.9, ls="--")
    axes[0].set(
        title="Does a lower housing cost restore renewal?",
        xlabel="Housing-cost ratio",
        ylabel="Effective locally born entrants / required entrants",
        xscale="log",
    )
    axes[1].set(
        title="Fertility response at the post-shock intercept",
        xlabel="Housing-cost ratio",
        ylabel="Model completed-fertility statistic",
        xscale="log",
    )
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "fertility_price_feedback.png", dpi=220)
    fig.savefig(outdir / "fertility_price_feedback.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    factors = [float(item) for item in args.kappa_factors.split(",") if item.strip()]
    if not factors or any(not 0.0 < value <= 1.0 for value in factors):
        raise ValueError("Kappa factors must lie in (0,1]")
    if args.smoke:
        factors = factors[:2]
    price_ratios = (
        (0.01, 0.10, 0.50, 1.00)
        if args.smoke
        else (0.0002, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00)
    )

    started = time.perf_counter()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve()
    winner, metadata = closure.load_winner(source)
    theta = closure.theta_from_winner(winner)
    chain = closure.load_chain(profile="e5f-floor")
    base = closure.make_overrides(chain, theta, nb=args.nb, profile="e5f-floor")
    old_price_overrides = dict(base)
    old_price_overrides["psi_child"] = float(args.old_psi_reference)
    baseline_solution, baseline_parameters, baseline_price, baseline_seconds = (
        closure.solve_ge(chain, old_price_overrides)
    )
    old_price = float(baseline_price[0])
    baseline = closure.readout(
        chain,
        baseline_solution,
        baseline_parameters,
        baseline_price,
        label="estimated_sequential_baseline",
        price_ratio=1.0,
        psi_child=float(args.old_psi_reference),
        elapsed=baseline_seconds,
    )
    original_kappa_first = float(theta["kappa_fert"])
    original_kappa_continuation = float(theta["kappa_fert_continuation"])

    summaries: list[dict[str, Any]] = []
    price_rows: list[dict[str, Any]] = []
    total_solves = 1
    for factor_index, factor in enumerate(factors, start=1):
        factor_started = time.perf_counter()
        scaled_base = dict(base)
        scaled_base["kappa_fert"] = original_kappa_first * factor
        scaled_base["kappa_fert_continuation"] = (
            original_kappa_continuation * factor
        )
        cache: dict[tuple[float, float], dict[str, Any]] = {}

        def evaluate(psi: float, asset_price: float) -> dict[str, Any]:
            nonlocal total_solves
            key = (round(float(psi), 10), round(float(asset_price), 10))
            if key not in cache:
                solution, parameters, price, elapsed = closure.solve_pe(
                    chain,
                    scaled_base,
                    asset_price=float(asset_price),
                    psi_child=float(psi),
                )
                cache[key] = closure.readout(
                    chain,
                    solution,
                    parameters,
                    price,
                    label="fertility_price_slope_sensitivity",
                    price_ratio=float(asset_price / old_price),
                    psi_child=float(psi),
                    elapsed=elapsed,
                )
                total_solves += 1
            return cache[key]

        psi_old = bisection(
            lambda psi: float(evaluate(psi, old_price)["tfr_topcoded"]),
            lower=-5.0,
            upper=5.0,
            target=float(args.old_fertility_target),
            increasing=True,
            value_tolerance=2e-4,
        )
        old_anchor = evaluate(psi_old, old_price)
        psi_replacement = bisection(
            lambda psi: float(evaluate(psi, old_price)[PRIMARY_RENEWAL_FIELD]),
            lower=-5.0,
            upper=5.0,
            target=1.0,
            increasing=True,
            value_tolerance=2e-5,
        )
        replacement_anchor = evaluate(psi_replacement, old_price)
        psi_low = bisection(
            lambda psi: float(evaluate(psi, old_price)["tfr_topcoded"]),
            lower=-5.0,
            upper=psi_old,
            target=float(args.low_fertility_target),
            increasing=True,
            value_tolerance=2e-4,
        )
        low_anchor = evaluate(psi_low, old_price)

        factor_price_rows: list[dict[str, Any]] = []
        for ratio in price_ratios:
            row = dict(evaluate(psi_low, old_price * ratio))
            row.update(
                kappa_factor=factor,
                kappa_fert=scaled_base["kappa_fert"],
                kappa_fert_continuation=scaled_base["kappa_fert_continuation"],
                anchor="post_shock_fertility_at_old_price",
            )
            factor_price_rows.append(row)
            price_rows.append(row)
            write_csv(outdir / "price_schedules.csv", price_rows)

        ordered = sorted(factor_price_rows, key=lambda row: float(row["asset_price"]))
        root = None
        root_bracket = None
        for left, right in zip(ordered[:-1], ordered[1:]):
            left_gap = float(left[PRIMARY_RENEWAL_FIELD]) - 1.0
            right_gap = float(right[PRIMARY_RENEWAL_FIELD]) - 1.0
            if left_gap == 0.0 or left_gap * right_gap <= 0.0:
                root_bracket = (float(left["asset_price"]), float(right["asset_price"]))
                break
        if root_bracket is not None:
            root_price = bisection(
                lambda price: float(evaluate(psi_low, price)[PRIMARY_RENEWAL_FIELD]),
                lower=root_bracket[0],
                upper=root_bracket[1],
                target=1.0,
                increasing=False,
                value_tolerance=2e-5,
            )
            root = dict(evaluate(psi_low, root_price))

        summary = {
            "kappa_factor": factor,
            "kappa_fert": scaled_base["kappa_fert"],
            "kappa_fert_continuation": scaled_base["kappa_fert_continuation"],
            "psi_matching_old_fertility": psi_old,
            "old_anchor_tfr": old_anchor["tfr_topcoded"],
            "old_anchor_B_over_E": old_anchor[PRIMARY_RENEWAL_FIELD],
            "old_anchor_raw_state_B_over_E": old_anchor[
                "reproduction_ratio_B_over_E"
            ],
            "old_anchor_childlessness": old_anchor["childless_rate"],
            "psi_imposing_replacement": psi_replacement,
            "replacement_anchor_tfr": replacement_anchor["tfr_topcoded"],
            "replacement_anchor_childlessness": replacement_anchor["childless_rate"],
            "psi_matching_low_fertility": psi_low,
            "low_anchor_tfr": low_anchor["tfr_topcoded"],
            "low_anchor_B_over_E": low_anchor[PRIMARY_RENEWAL_FIELD],
            "low_anchor_raw_state_B_over_E": low_anchor[
                "reproduction_ratio_B_over_E"
            ],
            "low_anchor_childlessness": low_anchor["childless_rate"],
            "maximum_B_over_E_on_price_grid": max(
                float(row[PRIMARY_RENEWAL_FIELD])
                for row in factor_price_rows
            ),
            "maximum_raw_state_B_over_E_on_price_grid": max(
                float(row["reproduction_ratio_B_over_E"])
                for row in factor_price_rows
            ),
            "closed_root_found": root is not None,
            "closed_root_price_ratio": (
                float(root["price_ratio"]) if root is not None else math.nan
            ),
            "closed_root_tfr": (
                float(root["tfr_topcoded"]) if root is not None else math.nan
            ),
            "closed_root_static_supply_scale": (
                float(root["scale_mapping_Hs_over_D"])
                if root is not None
                else math.nan
            ),
            "closed_root_fixed_old_stock_scale": (
                float(baseline["housing_supply"])
                / float(root["housing_demand_per_adult"])
                if root is not None
                else math.nan
            ),
            "factor_solve_count": len(cache),
            "elapsed_seconds": time.perf_counter() - factor_started,
        }
        summaries.append(summary)
        write_csv(outdir / "sensitivity_summary.csv", summaries)
        write_json(
            outdir / "latest_completed_factor.json",
            {
                "status": "running",
                "completed": factor_index,
                "requested": len(factors),
                "latest": summary,
            },
        )
        print(
            f"FACTOR {factor_index}/{len(factors)} x={factor:g} "
            f"old_B/E={summary['old_anchor_B_over_E']:.4f} "
            f"low_max_B/E={summary['maximum_B_over_E_on_price_grid']:.4f} "
            f"root={summary['closed_root_found']}",
            flush=True,
        )

    make_plot(price_rows, outdir)
    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    result = {
        "status": "complete",
        "interpretation": (
            "Uncalibrated slope sensitivity. Each row changes both fertility-logit "
            "scales and reanchors psi_child at the old house price. It is not a "
            "promotable calibration without a causal housing-cost-to-fertility moment."
        ),
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_metadata": metadata,
        "old_asset_price": old_price,
        "old_psi_reference": args.old_psi_reference,
        "old_fertility_target": args.old_fertility_target,
        "low_fertility_target": args.low_fertility_target,
        "primary_renewal_accounting": PRIMARY_RENEWAL_FIELD,
        "baseline": baseline,
        "original_kappa_fert": original_kappa_first,
        "original_kappa_fert_continuation": original_kappa_continuation,
        "kappa_factors": factors,
        "price_ratios": price_ratios,
        "total_model_solves": total_solves,
        "elapsed_seconds": time.perf_counter() - started,
        "rows": summaries,
        "exact_command": command,
    }
    write_json(outdir / "summary.json", result)
    (outdir / "README.md").write_text(
        "# Fertility-price feedback sensitivity\n\n"
        "This packet asks whether a lower housing price can restore demographic "
        "renewal after the fertility preference shift when the fertility choice "
        "is made more price-responsive. Both fertility-logit scales are changed "
        "together and the preference intercept is reanchored at the old price. "
        "Therefore this is an identification diagnostic, not a recalibration.\n\n"
        "The decisive columns are `maximum_B_over_E_on_price_grid` and "
        "`closed_root_found` in `sensitivity_summary.csv`. A different slope is "
        "not admissible for the paper until it is disciplined by causal "
        "housing-cost variation. Renewal uses the same 3+ top-bin child units "
        "as the completed-fertility statistic; raw three-child-state flows are "
        "reported alongside it.\n\n"
        f"Exact command:\n\n```bash\n{command}\n```\n",
        encoding="utf-8",
    )
    print("FEEDBACK_AUDIT_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
