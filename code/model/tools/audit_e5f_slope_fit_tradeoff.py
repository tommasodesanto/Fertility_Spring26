#!/usr/bin/env python3
"""Map the fit cost of making sequential fertility more price-responsive.

This is a fixed-parameter diagnostic, not a recalibration.  At each common
factor on the two fertility-logit noise scales, it holds every other saved
parameter fixed and chooses the nonnegative fertility intercept that comes as
close as possible to the current completed-fertility target.  It then reports
the complete current target-fit table and joins the previously audited closed-
population existence result.
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
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import audit_closed_reproductive_closure as closure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=closure.E5F_FLOOR_SOURCE)
    parser.add_argument(
        "--existence-audit",
        type=Path,
        default=ROOT / "output/model/sequential_fertility_price_feedback_audit",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "output/model/e5f_slope_fit_tradeoff",
    )
    parser.add_argument("--factors", default="1.0,0.5,0.25,0.10,0.05")
    parser.add_argument("--psi-upper", type=float, default=0.5)
    parser.add_argument("--tfr-tolerance", type=float, default=2e-4)
    parser.add_argument("--max-bisection-iterations", type=int, default=14)
    parser.add_argument("--nb", type=int, default=120)
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
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
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    factors = [float(item) for item in args.factors.split(",") if item.strip()]
    if not factors or any(not 0.0 < factor <= 1.0 for factor in factors):
        raise ValueError("Every factor must lie in (0, 1]")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve()
    winner, source_metadata = closure.load_winner(source)
    theta = closure.theta_from_winner(winner)
    chain = closure.load_chain(profile="e5f-floor")
    base = closure.make_overrides(
        chain, theta, nb=int(args.nb), profile="e5f-floor"
    )

    from intergen_eqscale_seq_optimized.e5_profile import (
        E5_PROFILE_NAME,
        E5_TARGET_PROVENANCE,
        e5_target_system,
    )

    system = e5_target_system()
    targets = system.targets_dict()
    weights = system.weights_dict()
    target_tfr = float(targets["tfr"])
    original_first = float(theta["kappa_fert"])
    original_continuation = float(theta["kappa_fert_continuation"])

    existence_rows = {
        float(row["kappa_factor"]): row
        for row in csv.DictReader(
            (args.existence_audit.resolve() / "sensitivity_summary.csv").open()
        )
    }
    missing = sorted(set(factors) - set(existence_rows))
    if missing:
        raise ValueError(f"Existence audit lacks factor(s): {missing}")
    existence_summary_path = args.existence_audit.resolve() / "summary.json"
    existence_summary = json.loads(existence_summary_path.read_text())
    if (
        existence_summary.get("status") != "complete"
        or existence_summary.get("source_sha256")
        != hashlib.sha256(source.read_bytes()).hexdigest()
    ):
        raise RuntimeError("The closed-root existence audit is incomplete or stale")

    started = time.perf_counter()
    factor_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    total_solves = 0
    baseline_gate: dict[str, Any] | None = None
    write_json(
        outdir / "summary.json",
        {
            "status": "running",
            "source": str(source),
            "requested_factors": factors,
        },
    )

    for factor_index, factor in enumerate(factors, start=1):
        factor_started = time.perf_counter()
        scaled = dict(base)
        scaled["kappa_fert"] = original_first * factor
        scaled["kappa_fert_continuation"] = original_continuation * factor
        cache: dict[float, tuple[Any, Any, np.ndarray, dict[str, float], float]] = {}

        def evaluate(psi: float) -> tuple[Any, Any, np.ndarray, dict[str, float], float]:
            nonlocal total_solves
            key = round(float(psi), 12)
            if key not in cache:
                overrides = dict(scaled)
                overrides["psi_child"] = float(psi)
                solution, parameters, price, elapsed = closure.solve_ge(
                    chain, overrides
                )
                timings = dict(getattr(solution, "timings", {}))
                residual = float(
                    getattr(solution, "best_max_abs_rel_excess", math.inf)
                )
                if (
                    not bool(
                        timings.get(
                            "strict_converged",
                            getattr(solution, "converged", False),
                        )
                    )
                    or residual > 2.5e-5
                ):
                    raise RuntimeError(
                        "Strict slope-fit solve failed: "
                        f"factor={factor:g}, psi={psi:.8f}, residual={residual:.6g}"
                    )
                moments = {
                    name: float(value)
                    for name, value in chain.extract_moments(
                        solution, parameters
                    ).items()
                    if np.isscalar(value)
                }
                cache[key] = (solution, parameters, price, moments, elapsed)
                total_solves += 1
                print(
                    "SLOPE_FIT_SOLVE "
                    f"factor={factor:g} psi={psi:.8f} "
                    f"tfr={moments['tfr']:.6f}",
                    flush=True,
                )
            return cache[key]

        at_lower = evaluate(0.0)
        lower_tfr = float(at_lower[3]["tfr"])
        if math.isclose(factor, 1.0, rel_tol=0.0, abs_tol=1e-14):
            source_moments = winner.get("moments") or {}
            differences = {
                name: float(at_lower[3][name]) - float(source_moments[name])
                for name in targets
            }
            baseline_gate = {
                "maximum_absolute_moment_difference": max(
                    abs(value) for value in differences.values()
                ),
                "price_difference": float(at_lower[2][0])
                - float(winner["price"]),
                "differences": differences,
            }
            if (
                baseline_gate["maximum_absolute_moment_difference"] > 2e-8
                or abs(baseline_gate["price_difference"]) > 2e-6
            ):
                raise RuntimeError(f"Source nesting gate failed: {baseline_gate}")

        intercept_status = "target_matched"
        if lower_tfr >= target_tfr:
            chosen_psi = 0.0
            chosen = at_lower
            intercept_status = "nonnegative_boundary_above_target"
        else:
            at_upper = evaluate(float(args.psi_upper))
            if float(at_upper[3]["tfr"]) < target_tfr:
                raise RuntimeError(
                    f"TFR target is not bracketed for factor {factor:g}"
                )
            lower, upper = 0.0, float(args.psi_upper)
            chosen_psi = 0.5 * (lower + upper)
            chosen = evaluate(chosen_psi)
            for _ in range(int(args.max_bisection_iterations)):
                gap = float(chosen[3]["tfr"]) - target_tfr
                if abs(gap) <= float(args.tfr_tolerance):
                    break
                if gap < 0.0:
                    lower = chosen_psi
                else:
                    upper = chosen_psi
                chosen_psi = 0.5 * (lower + upper)
                chosen = evaluate(chosen_psi)

        solution, parameters, price, moments, _ = chosen
        rows_for_factor: list[dict[str, Any]] = []
        canonical_loss = 0.0
        for name, target in targets.items():
            gap = float(moments[name]) - float(target)
            contribution = float(weights[name]) * gap**2
            canonical_loss += contribution
            row = {
                "kappa_factor": factor,
                "moment": name,
                "target": float(target),
                "model": float(moments[name]),
                "gap": gap,
                "weight": float(weights[name]),
                "loss_contribution": contribution,
            }
            rows_for_factor.append(row)
            fit_rows.append(row)

        existence = existence_rows[factor]
        closed_root = str(existence["closed_root_found"]).strip().lower() == "true"
        top_rows = sorted(
            rows_for_factor,
            key=lambda row: float(row["loss_contribution"]),
            reverse=True,
        )[:3]
        factor_row = {
            "kappa_factor": factor,
            "kappa_fert": original_first * factor,
            "kappa_fert_continuation": original_continuation * factor,
            "psi_child": chosen_psi,
            "psi_status": intercept_status,
            "tfr": float(moments["tfr"]),
            "tfr_gap": float(moments["tfr"]) - target_tfr,
            "childless_rate": float(moments["childless_rate"]),
            "mean_age_first_birth": float(moments["mean_age_first_birth"]),
            "share_first_births_age30plus": float(
                moments["share_first_births_age30plus"]
            ),
            "housing_increment_0to1": float(moments["housing_increment_0to1"]),
            "family_size_room_gradient": float(
                moments["prime30_55_parent_3plus_minus_1to2_mean_rooms"]
            ),
            "own_family_gap": float(moments["own_family_gap"]),
            "own_rate": float(moments["own_rate"]),
            "asset_price": float(price[0]),
            "canonical_loss": canonical_loss,
            "largest_loss_moment_1": top_rows[0]["moment"],
            "largest_loss_contribution_1": top_rows[0]["loss_contribution"],
            "largest_loss_moment_2": top_rows[1]["moment"],
            "largest_loss_contribution_2": top_rows[1]["loss_contribution"],
            "largest_loss_moment_3": top_rows[2]["moment"],
            "largest_loss_contribution_3": top_rows[2]["loss_contribution"],
            "maximum_postshock_B_over_E": float(
                existence["maximum_B_over_E_on_price_grid"]
            ),
            "closed_root_found": closed_root,
            "closed_root_price_ratio": (
                float(existence["closed_root_price_ratio"])
                if closed_root
                else math.nan
            ),
            "factor_solve_count": len(cache),
            "factor_elapsed_seconds": time.perf_counter() - factor_started,
        }
        factor_rows.append(factor_row)
        write_csv(outdir / "slope_fit_summary.csv", factor_rows)
        write_csv(outdir / "target_fit_by_factor.csv", fit_rows)
        write_json(
            outdir / "latest_completed_factor.json",
            {"completed": factor_index, "total": len(factors), "latest": factor_row},
        )

    if baseline_gate is None:
        raise RuntimeError("The factor grid must include the estimated scale 1")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = sorted(factor_rows, key=lambda row: float(row["kappa_factor"]))
    x = np.asarray([float(row["kappa_factor"]) for row in ordered])
    loss = np.asarray([float(row["canonical_loss"]) for row in ordered])
    late = np.asarray([float(row["share_first_births_age30plus"]) for row in ordered])
    childless = np.asarray([float(row["childless_rate"]) for row in ordered])
    renewal = np.asarray([float(row["maximum_postshock_B_over_E"]) for row in ordered])
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)
    axes[0].plot(x, loss, marker="o", lw=2)
    axes[0].set(xscale="log", yscale="log", title="Full calibration loss", xlabel="Fertility-logit scale factor")
    axes[1].plot(x, childless, marker="o", lw=2, label="Childlessness")
    axes[1].plot(x, late, marker="s", lw=2, label="First births at 30+")
    axes[1].axhline(targets["childless_rate"], color="C0", ls="--", lw=1)
    axes[1].axhline(targets["share_first_births_age30plus"], color="C1", ls="--", lw=1)
    axes[1].set(xscale="log", title="Fertility composition", xlabel="Fertility-logit scale factor")
    axes[1].legend(frameon=False)
    axes[2].plot(x, renewal, marker="o", lw=2)
    axes[2].axhline(1.0, color="black", ls="--", lw=1)
    axes[2].set(xscale="log", title="Best renewal over audited prices", xlabel="Fertility-logit scale factor")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "slope_fit_tradeoff.png", dpi=220)
    fig.savefig(outdir / "slope_fit_tradeoff.pdf")
    plt.close(fig)

    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    write_json(
        outdir / "summary.json",
        {
            "status": "complete",
            "interpretation": "fixed-other-parameter slope and intercept diagnostic; not a recalibration",
            "source": str(source),
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "source_target_set": source_metadata.get("target_set"),
            "existence_audit": str(args.existence_audit.resolve()),
            "existence_audit_summary_sha256": hashlib.sha256(
                existence_summary_path.read_bytes()
            ).hexdigest(),
            "current_target_profile": E5_PROFILE_NAME,
            "current_target_fingerprint": system.fingerprint,
            "current_target_provenance": E5_TARGET_PROVENANCE,
            "target_count": system.count,
            "estimated_source_parameter_count": 9,
            "rooms_target_status": "empirical_hold",
            "baseline_nesting_gate": baseline_gate,
            "factors": factor_rows,
            "total_model_solves": total_solves,
            "elapsed_seconds": time.perf_counter() - started,
            "exact_command": command,
        },
    )
    print(f"SLOPE_FIT_TRADEOFF_COMPLETE outdir={outdir}", flush=True)


if __name__ == "__main__":
    main()
