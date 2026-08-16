#!/usr/bin/env python3
"""Audit the proposed closed reproductive closure in a sequential model.

This is a diagnostic, not a calibration or a transition solver.  It evaluates
normalized stationary fixed-price solutions and records

    E(r) = required young-adult entrant households,
    B(r) = mature city-born entrant households,
    R(r) = B(r) - E(r), and
    F(r) = H^S(r) / D(r).

The historical default is the retained shared-clock E5b point.  The
``e5f-floor`` profile instead routes the exact child-room-floor architecture
shown in the August 7 Raquel slides.  A constructed preference normalization
is retained only for the historical full audit; ``--schedule-only`` stops
after the actual calibrated reproductive and housing schedules.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
DEFAULT_SOURCE = (
    ROOT
    / "output/model/closed_reproductive_closure_audit_20260811/source/"
    "incumbent_current_contract_preflight_summary.json"
)
FALLBACK_SOURCE = ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
REPAIRED_SOURCE = (
    ROOT
    / "output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json"
)
E5F_FLOOR_SOURCE = (
    ROOT
    / "output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/closed_reproductive_closure_audit_20260811"
PRICE_RATIOS = (0.01, 0.025, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 1.00, 1.10, 1.25, 1.50, 2.00)
E5F_PRICE_RATIOS = (0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 1.00, 1.10, 1.25, 1.50, 2.00, 3.00)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def configure_environment(*, profile: str) -> None:
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E3_TFR_TOP_BIN_WEIGHT"] = "3.602359422009"
    for name in ("E5_MATURATION_REPAIR", "E5F", "E6A", "E6B", "E6C"):
        os.environ.pop(name, None)
    if profile in {"repaired-e5", "e5f-floor"}:
        os.environ["E5_MATURATION_REPAIR"] = "1"
    if profile == "e5f-floor":
        os.environ["E5F"] = "1"


def load_chain(*, profile: str = "maintained-e5b") -> Any:
    configure_environment(profile=profile)
    if str(MODEL_ROOT) not in sys.path:
        sys.path.insert(0, str(MODEL_ROOT))
    from intergen_eqscale_seq_optimized import run_e1_chain

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    return chain


def load_winner(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data.get("best_tight"), dict):
        winner = data["best_tight"]
        metadata = data.get("metadata", {})
    elif isinstance(data.get("winners", {}).get("E1"), dict):
        winner = data["winners"]["E1"]
        metadata = data.get("winner_metadata", {})
    else:
        raise ValueError(f"Cannot locate a winner record in {path}")
    if not bool(winner.get("strict_converged", False)) or winner.get("status") != "ok":
        raise RuntimeError(f"Source winner is not a strict accepted solution: {path}")
    return winner, metadata


def theta_from_winner(winner: dict[str, Any]) -> dict[str, float]:
    return {str(key): float(value) for key, value in winner["theta"].items()}


def make_overrides(
    chain: Any,
    theta: dict[str, float],
    *,
    nb: int,
    profile: str = "maintained-e5b",
) -> dict[str, Any]:
    overrides = chain.common_overrides(
        argparse.Namespace(J=17, Nb=int(nb), max_iter_eq=40, tol_eq=2.5e-5)
    )
    overrides.update(theta)
    if profile == "e5f-floor":
        overrides["child_state_mode"] = "independent_count"
        if not bool(overrides.get("child_room_floor", False)):
            raise RuntimeError("E5F routing failed to activate the child-room floor")
    return overrides


def solve_ge(chain: Any, base: dict[str, Any]) -> tuple[Any, Any, np.ndarray, float]:
    started = time.perf_counter()
    sol, parameters, price = chain.run_model_cp_dt(dict(base), verbose=False)
    return sol, parameters, np.asarray(price, dtype=float).reshape(-1), time.perf_counter() - started


def solve_pe(
    chain: Any,
    base: dict[str, Any],
    *,
    asset_price: float,
    psi_child: float | None = None,
) -> tuple[Any, Any, np.ndarray, float]:
    overrides = dict(base)
    overrides.update(solve_mode="pe", p_fixed=np.array([float(asset_price)]))
    if psi_child is not None:
        overrides["psi_child"] = float(psi_child)
    started = time.perf_counter()
    sol, parameters, price = chain.run_model_cp_dt(overrides, verbose=False)
    return sol, parameters, np.asarray(price, dtype=float).reshape(-1), time.perf_counter() - started


def solve_fiscally_balanced(
    chain: Any,
    base: dict[str, Any],
    *,
    asset_price: float | None,
) -> tuple[Any, Any, np.ndarray, float, float, int]:
    """Solve the rebated 1-percent property-tax baseline at a fixed price or in GE."""
    from run_intergen_funded_property_tax_test import case_overrides

    cache: dict[float, tuple[Any, Any, np.ndarray]] = {}
    started = time.perf_counter()

    def solve_at(transfer: float) -> tuple[Any, Any, np.ndarray]:
        key = float(round(transfer, 12))
        if key not in cache:
            overrides = case_overrides(base, 0.01, False, float(transfer), False)
            if asset_price is not None:
                overrides.update(solve_mode="pe", p_fixed=np.array([float(asset_price)]))
            cache[key] = chain.run_model_cp_dt(overrides, verbose=False)
        return cache[key]

    def residual(transfer: float) -> float:
        return float(solve_at(transfer)[0].property_tax_budget_residual)

    lower, upper = 0.0, 0.5
    f_lower, f_upper = residual(lower), residual(upper)
    while f_upper > 0.0 and upper < 8.0:
        upper *= 2.0
        f_upper = residual(upper)
    if f_lower < 0.0 or f_upper > 0.0:
        raise RuntimeError(
            "Could not bracket the rebated-baseline transfer: "
            f"F(0)={f_lower:.6g}, F({upper:g})={f_upper:.6g}"
        )

    transfer = 0.5 * (lower + upper)
    for _ in range(30):
        denominator = f_upper - f_lower
        secant = lower - f_lower * (upper - lower) / denominator if denominator != 0.0 else math.nan
        span = upper - lower
        if (
            not math.isfinite(secant)
            or secant <= lower + 0.05 * span
            or secant >= upper - 0.05 * span
        ):
            secant = 0.5 * (lower + upper)
        transfer = float(secant)
        f_transfer = residual(transfer)
        if abs(f_transfer) <= 2.5e-5:
            break
        if f_transfer > 0.0:
            lower, f_lower = transfer, f_transfer
        else:
            upper, f_upper = transfer, f_transfer
    solution, parameters, price = solve_at(transfer)
    return (
        solution,
        parameters,
        np.asarray(price, dtype=float).reshape(-1),
        time.perf_counter() - started,
        transfer,
        len(cache),
    )


def shared_clock_mature_flows_by_parity(sol: Any, parameters: Any) -> np.ndarray:
    """Reconstruct the code's shared-clock mature flow, split by parity."""
    if str(parameters.child_state_mode).strip().lower() != "shared_clock":
        return np.full(int(parameters.n_parity), np.nan)
    g = np.asarray(sol.g, dtype=float)
    if g.ndim != 7:
        raise RuntimeError(f"Unexpected maintained distribution shape: {g.shape}")
    K = int(parameters.n_child_stages)
    csm1, csm2 = K + 1, K + 2
    ecf = float(parameters.entrant_conversion_factor)
    flows = np.zeros(int(parameters.n_parity), dtype=float)
    for j in range(int(parameters.J) - 1):
        survival = (
            float(parameters.survival_probs[j])
            if bool(getattr(parameters, "use_age_survival", False))
            else 1.0
        )
        for parity in range(1, int(parameters.n_parity)):
            target = csm1 if parity == 1 else csm2
            probability = float(parameters.Pi_child[K, target, parity])
            mass = float(np.sum(g[:, :, :, j, :, parity, K]))
            flows[parity] += ecf * parity * probability * survival * mass
    actual = float(sol.entrants_mature_total)
    if not math.isclose(float(np.sum(flows)), actual, rel_tol=0.0, abs_tol=1e-10):
        raise RuntimeError(
            "Parity-flow reconstruction failed: "
            f"reconstructed={np.sum(flows):.12g}, actual={actual:.12g}"
        )
    return flows


def topcode_consistent_renewal_accounting(
    sol: Any,
    parameters: Any,
) -> dict[str, Any]:
    """Put measured fertility and demographic renewal in the same child units.

    The sequential model's last parity state represents ``3+`` births.  Its
    calibration moment therefore weights that state by the observed mean in
    the top-coded bin, while the transition kernel can only carry three
    explicit children.  For population renewal, we assign the extra top-bin
    children the same maturation yield as explicit children.  This leaves the
    household problem unchanged and makes the demographic accounting
    consistent with the reported completed-fertility statistic.
    """
    raw_births = float(getattr(sol, "total_births_kfe", math.nan))
    raw_mature = float(getattr(sol, "entrants_mature_total", math.nan))
    conversion = float(getattr(parameters, "entrant_conversion_factor", math.nan))
    if not (
        math.isfinite(raw_births)
        and raw_births > 0.0
        and math.isfinite(raw_mature)
        and raw_mature >= 0.0
        and math.isfinite(conversion)
        and conversion > 0.0
    ):
        raise RuntimeError(
            "Invalid raw birth or maturation flow for renewal accounting: "
            f"births={raw_births}, mature={raw_mature}, conversion={conversion}"
        )

    top_state = int(parameters.n_parity) - 1
    top_weight = float(getattr(parameters, "tfr_top_bin_weight", top_state))
    fertility_units = str(
        getattr(parameters, "fertility_units", "parity2x")
    ).strip().lower()
    extra_per_top_family = (
        max(top_weight - float(top_state), 0.0)
        if fertility_units == "literal_topcode"
        else 0.0
    )
    maturation_yield = raw_mature / (conversion * raw_births)
    if not 0.0 <= maturation_yield <= 1.0 + 1e-10:
        raise RuntimeError(
            "Invalid explicit-child maturation yield: "
            f"{maturation_yield:.12g}"
        )

    top_entry_birth_flow = 0.0
    top_flow_share = 0.0
    method = "raw_explicit_children"
    adjusted_births = raw_births
    adjusted_mature = raw_mature

    if extra_per_top_family > 0.0:
        if bool(getattr(parameters, "sequential_births", False)):
            if top_state != 3:
                raise RuntimeError(
                    "Sequential top-code accounting currently requires parity states 0/1/2/3+."
                )
            third_births_by_age = np.asarray(
                getattr(parameters, "_third_births_by_age", np.array([])),
                dtype=float,
            ).reshape(-1)
            if third_births_by_age.size == 0 or not np.all(
                np.isfinite(third_births_by_age)
            ):
                raise RuntimeError(
                    "Sequential top-bin entry flow is unavailable after the KFE solve."
                )
            top_entry_birth_flow = float(np.sum(third_births_by_age))
            adjusted_births = (
                raw_births + extra_per_top_family * top_entry_birth_flow
            )
            adjusted_mature = conversion * maturation_yield * adjusted_births
            top_flow_share = top_entry_birth_flow / raw_births
            method = "sequential_top_bin_birth_flow_common_maturation_yield"
        else:
            parity_flows = shared_clock_mature_flows_by_parity(sol, parameters)
            if not np.all(np.isfinite(parity_flows)):
                raise RuntimeError(
                    "Top-code-consistent mature flow could not be reconstructed."
                )
            adjusted_mature = float(
                np.sum(parity_flows[:-1])
                + parity_flows[-1] * top_weight / float(top_state)
            )
            top_flow_share = float(
                parity_flows[-1] / max(np.sum(parity_flows), 1e-15)
            )
            method = "shared_clock_mature_flow_by_parity"

    return {
        "method": method,
        "raw_birth_children": raw_births,
        "top_bin_entry_birth_flow": top_entry_birth_flow,
        "extra_children_per_top_bin_family": extra_per_top_family,
        "topcode_adjusted_birth_children": adjusted_births,
        "raw_mature_entrant_households": raw_mature,
        "topcode_adjusted_mature_entrant_households": adjusted_mature,
        "explicit_child_maturation_yield": maturation_yield,
        "top_bin_flow_share": top_flow_share,
    }


def readout(
    chain: Any,
    sol: Any,
    parameters: Any,
    price: np.ndarray,
    *,
    label: str,
    price_ratio: float,
    psi_child: float,
    elapsed: float,
) -> dict[str, Any]:
    moments = chain.extract_moments(sol, parameters)
    total_mass = float(sol.total_mass)
    entry = float(sol.entry_rate) / total_mass
    mature = float(sol.entrants_mature_total) / total_mass
    births = float(sol.total_births_kfe) / total_mass
    demand = float(np.sum(np.asarray(sol.housing_demand, dtype=float))) / total_mass
    supply = float(np.sum(np.asarray(sol.housing_supply, dtype=float)))
    renewal = topcode_consistent_renewal_accounting(sol, parameters)
    adjusted_births = float(renewal["topcode_adjusted_birth_children"]) / total_mass
    adjusted_mature = float(
        renewal["topcode_adjusted_mature_entrant_households"]
    ) / total_mass
    top_entry_birth_flow = float(renewal["top_bin_entry_birth_flow"]) / total_mass
    top_flow_share = float(renewal["top_bin_flow_share"])
    return {
        "label": label,
        "price_ratio": float(price_ratio),
        "asset_price": float(price[0]),
        "housing_user_cost": float(parameters.user_cost_rate * price[0]),
        "psi_child": float(psi_child),
        "total_mass": total_mass,
        "entry_households_per_adult": entry,
        "birth_children_per_adult": births,
        "top_bin_entry_birth_flow_per_adult": top_entry_birth_flow,
        "topcode_adjusted_birth_children_per_adult": adjusted_births,
        "mature_entrant_households_per_adult": mature,
        "topcode_adjusted_mature_entrant_households_per_adult": adjusted_mature,
        "reproduction_ratio_B_over_E": mature / entry,
        "raw_state_B_over_E": mature / entry,
        "reproduction_residual_B_minus_E": mature - entry,
        "topcode_adjusted_B_over_E_diagnostic": adjusted_mature / entry,
        "topcode_consistent_B_over_E": adjusted_mature / entry,
        "renewal_child_accounting_method": str(renewal["method"]),
        "extra_children_per_top_bin_family": float(
            renewal["extra_children_per_top_bin_family"]
        ),
        "threeplus_share_of_mature_flow": top_flow_share,
        "maturation_exit_yield": mature / max(float(parameters.entrant_conversion_factor) * births, 1e-15),
        "entrant_conversion_factor": float(parameters.entrant_conversion_factor),
        "housing_demand_per_adult": demand,
        "housing_supply": supply,
        "scale_mapping_Hs_over_D": supply / demand,
        "tfr_topcoded": float(moments["tfr"]),
        "mean_completed_fertility_raw": float(moments["mean_completed_fertility"]),
        "childless_rate": float(moments["childless_rate"]),
        "own_rate": float(moments["own_rate"]),
        "solve_seconds": float(elapsed),
    }


def verify_source_solution(chain: Any, sol: Any, parameters: Any, price: np.ndarray, winner: dict[str, Any]) -> None:
    if abs(float(price[0]) - float(winner["price"])) > 2e-6:
        raise RuntimeError(f"Source price gate failed: solved={price[0]}, source={winner['price']}")
    moments = chain.extract_moments(sol, parameters)
    failures = []
    for row in winner.get("target_fit", []):
        name = str(row["moment"])
        difference = abs(float(moments[name]) - float(row["model"]))
        if difference > 1e-6:
            failures.append((name, difference))
    if failures:
        raise RuntimeError(f"Source target-moment gate failed: {failures}")


def checkpoint(outdir: Path, phase: str, rows: list[dict[str, Any]]) -> None:
    write_json(
        outdir / "latest_completed_case.json",
        {"phase": phase, "completed_cases": len(rows), "latest": rows[-1]},
    )


def monotone(values: list[float], *, increasing: bool, tolerance: float = 1e-10) -> bool:
    differences = np.diff(np.asarray(values, dtype=float))
    return bool(np.all(differences >= -tolerance) if increasing else np.all(differences <= tolerance))


def log_slope(rows: list[dict[str, Any]], field: str, low_ratio: float = 0.8, high_ratio: float = 1.25) -> float:
    low = min(rows, key=lambda row: abs(float(row["price_ratio"]) - low_ratio))
    high = min(rows, key=lambda row: abs(float(row["price_ratio"]) - high_ratio))
    return float(
        math.log(float(high[field]) / float(low[field]))
        / math.log(float(high["housing_user_cost"]) / float(low["housing_user_cost"]))
    )


def bisection(
    evaluator: Callable[[float], float],
    *,
    lower: float,
    upper: float,
    target: float,
    increasing: bool,
    tolerance: float = 2e-5,
    max_iterations: int = 18,
    label: str,
) -> float:
    flo, fhi = evaluator(lower) - target, evaluator(upper) - target
    if increasing:
        bracketed = flo <= 0.0 <= fhi
    else:
        bracketed = flo >= 0.0 >= fhi
    if not bracketed:
        raise RuntimeError(
            f"{label} root is not bracketed: lower residual={flo:.6g}, upper residual={fhi:.6g}"
        )
    lo, hi = float(lower), float(upper)
    for iteration in range(1, max_iterations + 1):
        mid = 0.5 * (lo + hi)
        residual = evaluator(mid) - target
        print(f"{label.upper()} iteration={iteration} x={mid:.12g} residual={residual:.6g}", flush=True)
        if abs(residual) <= tolerance:
            return mid
        if increasing:
            if residual < 0.0:
                lo = mid
            else:
                hi = mid
        else:
            if residual > 0.0:
                lo = mid
            else:
                hi = mid
    return 0.5 * (lo + hi)


def parameter_rows(theta: dict[str, float], metadata: dict[str, Any]) -> list[dict[str, Any]]:
    domain = {str(row["name"]): row for row in metadata.get("active_domain", [])}
    rows: list[dict[str, Any]] = []
    for name, spec in domain.items():
        theta_name = "beta" if name == "beta_annual" else name
        estimate = float(theta[theta_name]) ** 0.25 if name == "beta_annual" else float(theta[theta_name])
        lower, upper = float(spec["lower"]), float(spec["upper"])
        distance = min(estimate - lower, upper - estimate)
        span = upper - lower
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower_bound": lower,
                "upper_bound": upper,
                "transform": spec["transform"],
                "distance_as_share_of_range": distance / span,
                "near_bound_2pct": distance / span <= 0.02,
                "status": "estimated",
            }
        )
    external = metadata.get("external_parameter_metadata", {})
    for name, spec in external.items():
        rows.append(
            {
                "parameter": name,
                "estimate": float(theta[name]),
                "lower_bound": "",
                "upper_bound": "",
                "transform": "",
                "distance_as_share_of_range": "",
                "near_bound_2pct": "",
                "status": f"externally fixed: {spec.get('source', '')}",
            }
        )
    return rows


def run_architecture_only(args: argparse.Namespace) -> None:
    chain = load_chain(profile="repaired-e5")
    main_winner, _ = load_winner(args.source)
    repaired_winner, _ = load_winner(REPAIRED_SOURCE)
    records = []
    for label, winner in (
        ("same_retained_theta_independent_maturation", main_winner),
        ("reestimated_Aug5_independent_maturation", repaired_winner),
    ):
        base = make_overrides(
            chain,
            theta_from_winner(winner),
            nb=args.nb,
            profile="repaired-e5",
        )
        sol, parameters, price, elapsed = solve_ge(chain, base)
        records.append(
            readout(
                chain,
                sol,
                parameters,
                price,
                label=label,
                price_ratio=math.nan,
                psi_child=float(base["psi_child"]),
                elapsed=elapsed,
            )
        )
        print(f"ARCHITECTURE {label} B/E={records[-1]['reproduction_ratio_B_over_E']:.9f}", flush=True)
    write_json(args.architecture_json, records)


def make_schedule_plot(outdir: Path, rows: list[dict[str, Any]], *, title: str) -> None:
    ratios = np.asarray([row["price_ratio"] for row in rows], dtype=float)
    raw_reproduction = np.asarray(
        [row["reproduction_ratio_B_over_E"] for row in rows], dtype=float
    )
    reproduction = np.asarray(
        [row["topcode_consistent_B_over_E"] for row in rows], dtype=float
    )
    scale = np.asarray([row["scale_mapping_Hs_over_D"] for row in rows], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    axes[0].plot(
        ratios,
        reproduction,
        marker="o",
        lw=2,
        color="#21476b",
        label="3+ bin in measured child units",
    )
    axes[0].plot(
        ratios,
        raw_reproduction,
        marker="s",
        lw=1.4,
        ls="--",
        color="#7f8c8d",
        label="Three explicit child states",
    )
    axes[0].axhline(1.0, color="black", lw=1, ls=":")
    axes[0].set(
        xscale="log",
        xlabel="Housing-cost ratio",
        ylabel=r"Effective reproduction $B/E$",
        title="Fertility and replacement",
    )
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(ratios, scale, marker="o", lw=2, color="#b14f32")
    axes[1].set(
        xscale="log",
        yscale="log",
        xlabel="Housing-cost ratio",
        ylabel=r"Stationary scale $H^S/D$",
        title="Housing clearing conditional on cost",
    )
    axes[1].grid(alpha=0.25)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outdir / "calibrated_schedule.png", dpi=220)
    fig.savefig(outdir / "calibrated_schedule.pdf")
    plt.close(fig)


def make_plots(
    outdir: Path,
    price_rows: list[dict[str, Any]],
    normalized_rows: list[dict[str, Any]],
    shock_rows: list[dict[str, Any]],
    architecture_rows: list[dict[str, Any]],
    constructed: dict[str, Any],
) -> None:
    ratios = np.asarray([row["price_ratio"] for row in price_rows], dtype=float)
    raw = np.asarray([row["reproduction_ratio_B_over_E"] for row in price_rows], dtype=float)
    adjusted = np.asarray([row["topcode_adjusted_B_over_E_diagnostic"] for row in price_rows], dtype=float)
    scale = np.asarray([row["scale_mapping_Hs_over_D"] for row in price_rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    axes[0].plot(ratios, raw, marker="o", lw=2, label="Model KFE flow")
    axes[0].plot(ratios, adjusted, marker="s", lw=1.6, ls="--", label="3+ top-code diagnostic")
    axes[0].axhline(1.0, color="black", lw=1, ls=":")
    axes[0].set(xscale="log", xlabel="Asset price / maintained price", ylabel=r"Effective reproduction $B/E$", title="No reproductive root at maintained preferences")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(ratios, scale, marker="o", lw=2, color="#a34a28")
    axes[1].axhline(1.0, color="black", lw=1, ls=":")
    axes[1].set(xscale="log", yscale="log", xlabel="Asset price / maintained price", ylabel=r"Scale map $H^S/D$", title="Housing block assigns scale conditional on price")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "maintained_closure_grid.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    for rows, label, style in (
        (normalized_rows, "Replacement-normalized preference", "-"),
        (shock_rows, "Half-capacity taste deterioration", "--"),
    ):
        x = np.asarray([row["price_ratio"] for row in rows], dtype=float)
        y = np.asarray([row["reproduction_ratio_B_over_E"] for row in rows], dtype=float)
        axes[0].plot(x, y, marker="o", lw=1.8, ls=style, label=label)
    axes[0].axhline(1.0, color="black", lw=1, ls=":")
    axes[0].axvline(float(constructed["new_root_price_ratio"]), color="#a34a28", lw=1, ls=":")
    axes[0].set(xscale="log", xlabel="Asset price / maintained price", ylabel=r"Effective reproduction $B/E$", title="Constructed closure: only a very small shock is recoverable")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    components = [
        float(constructed["price_effect_log_scale"]),
        float(constructed["direct_composition_effect_log_scale"]),
        float(constructed["total_log_scale_change"]),
    ]
    colors = ["#4472c4", "#70ad47", "#a34a28"]
    axes[1].bar(["Price\nmovement", "Direct taste /\ncomposition", "Total"], components, color=colors)
    axes[1].axhline(0.0, color="black", lw=1)
    axes[1].set(ylabel="Log change in population scale", title="Population-scale decomposition")
    axes[1].grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "constructed_normalized_closure.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    labels = [
        "Maintained\nshared clock",
        "Same theta,\nindependent",
        "Re-estimated\nindependent",
    ]
    values = [float(row["reproduction_ratio_B_over_E"]) for row in architecture_rows]
    ax.bar(labels, values, color=["#4472c4", "#ed7d31", "#70ad47"])
    ax.axhline(1.0, color="black", lw=1, ls=":")
    for index, value in enumerate(values):
        ax.text(index, value + 0.015, f"{value:.3f}", ha="center", va="bottom")
    ax.set(ylabel=r"Effective reproduction $B/E$", title="Maturation architecture is quantitatively material")
    ax.set_ylim(0.0, max(values) * 1.18)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "architecture_sensitivity.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("maintained-e5b", "e5f-floor"),
        default="maintained-e5b",
    )
    parser.add_argument("--source", type=Path)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--nb", type=int, default=120)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--fiscal-convention",
        choices=("calibration-unrebated", "rebated-1pct"),
        default="calibration-unrebated",
    )
    parser.add_argument(
        "--schedule-only",
        action="store_true",
        help="Stop after the calibrated price schedule and any bracketed closed root.",
    )
    parser.add_argument("--skip-architecture", action="store_true")
    parser.add_argument("--architecture-only", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--architecture-json", type=Path, default=DEFAULT_OUTDIR / "architecture_helper.json", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.source is None:
        if args.profile == "e5f-floor":
            args.source = E5F_FLOOR_SOURCE
        else:
            args.source = DEFAULT_SOURCE if DEFAULT_SOURCE.exists() else FALLBACK_SOURCE
    args.source = args.source.resolve()
    args.outdir = args.outdir.resolve()
    args.architecture_json = args.architecture_json.resolve()

    if args.architecture_only:
        run_architecture_only(args)
        return

    args.outdir.mkdir(parents=True, exist_ok=True)
    winner, metadata = load_winner(args.source)
    theta = theta_from_winner(winner)
    chain = load_chain(profile=args.profile)
    base = make_overrides(chain, theta, nb=args.nb, profile=args.profile)

    print(f"BASELINE_START source={args.source}", flush=True)
    if args.fiscal_convention == "rebated-1pct":
        sol0, parameters0, price0, elapsed0, transfer0, fiscal_evaluations0 = solve_fiscally_balanced(
            chain,
            base,
            asset_price=None,
        )
    else:
        sol0, parameters0, price0, elapsed0 = solve_ge(chain, base)
        transfer0, fiscal_evaluations0 = 0.0, 1
        verify_source_solution(chain, sol0, parameters0, price0, winner)
    baseline = readout(
        chain,
        sol0,
        parameters0,
        price0,
        label=f"{args.profile}_GE",
        price_ratio=1.0,
        psi_child=float(theta["psi_child"]),
        elapsed=elapsed0,
    )
    baseline.update(
        fiscal_convention=args.fiscal_convention,
        lump_sum_transfer_period_units=float(transfer0),
        fiscal_root_evaluations=int(fiscal_evaluations0),
        fiscal_residual=float(getattr(sol0, "property_tax_budget_residual", math.nan)),
    )
    write_csv(args.outdir / "baseline_accounting.csv", [baseline])
    print(
        "BASELINE_DONE "
        f"raw_B/E={baseline['reproduction_ratio_B_over_E']:.9f} "
        f"topcode_B/E={baseline['topcode_consistent_B_over_E']:.9f}",
        flush=True,
    )

    if args.smoke:
        ratios = (0.10, 1.0, 2.0) if args.profile == "e5f-floor" else (0.35, 1.0, 1.50)
    elif args.profile == "e5f-floor":
        ratios = E5F_PRICE_RATIOS
    else:
        ratios = PRICE_RATIOS
    price_rows: list[dict[str, Any]] = []
    for index, ratio in enumerate(ratios, start=1):
        if args.fiscal_convention == "rebated-1pct":
            sol, parameters, price, elapsed, transfer, fiscal_evaluations = solve_fiscally_balanced(
                chain,
                base,
                asset_price=float(price0[0] * ratio),
            )
        else:
            sol, parameters, price, elapsed = solve_pe(
                chain, base, asset_price=float(price0[0] * ratio), psi_child=float(theta["psi_child"])
            )
            transfer, fiscal_evaluations = 0.0, 1
        row = readout(
            chain,
            sol,
            parameters,
            price,
            label=f"{args.profile}_calibrated_price_grid",
            price_ratio=float(ratio),
            psi_child=float(theta["psi_child"]),
            elapsed=elapsed,
        )
        row.update(
            fiscal_convention=args.fiscal_convention,
            lump_sum_transfer_period_units=float(transfer),
            fiscal_root_evaluations=int(fiscal_evaluations),
            fiscal_residual=float(getattr(sol, "property_tax_budget_residual", math.nan)),
        )
        price_rows.append(row)
        write_csv(args.outdir / "price_grid.csv", price_rows)
        checkpoint(args.outdir, "maintained_price_grid", price_rows)
        print(
            f"PRICE_CASE {index}/{len(ratios)} ratio={ratio:.4g} "
            f"raw_B/E={row['reproduction_ratio_B_over_E']:.9f} "
            f"topcode_B/E={row['topcode_consistent_B_over_E']:.9f} "
            f"F={row['scale_mapping_Hs_over_D']:.9g}",
            flush=True,
        )

    if args.smoke:
        if args.profile == "e5f-floor" and args.fiscal_convention == "rebated-1pct":
            quota_path = (
                ROOT
                / "output/model/eqscale_seq_e5_policy_quota_closure_20260807/raw_runs/floor/quota/benchmark_quota_objects.json"
            )
            quota = json.loads(quota_path.read_text(encoding="utf-8"))
            for field, expected in (
                ("entry_households_per_adult", float(quota["baseline_entry_flow"])),
                ("mature_entrant_households_per_adult", float(quota["baseline_mature_cityborn_flow"])),
            ):
                if abs(float(baseline[field]) - expected) > 1e-7:
                    raise RuntimeError(
                        f"E5F baseline nesting failed for {field}: "
                        f"solved={baseline[field]:.12g}, expected={expected:.12g}"
                    )
        if args.fiscal_convention == "rebated-1pct":
            fiscal_residuals = [abs(float(row["fiscal_residual"])) for row in [baseline, *price_rows]]
            if max(fiscal_residuals) > 2.5e-5:
                raise RuntimeError(
                    f"Rebated-baseline fiscal gate failed: max residual={max(fiscal_residuals):.6g}"
                )
        write_json(
            args.outdir / "smoke_summary.json",
            {
                "status": "passed",
                "profile": args.profile,
                "fiscal_convention": args.fiscal_convention,
                "source": args.source,
                "baseline": baseline,
                "price_grid": price_rows,
            },
        )
        make_schedule_plot(args.outdir, price_rows, title=f"{args.profile}: smoke schedule")
        print("SMOKE_COMPLETE", flush=True)
        return

    if args.schedule_only:
        primary_ratio_field = (
            "topcode_consistent_B_over_E"
            if args.profile == "e5f-floor"
            else "reproduction_ratio_B_over_E"
        )
        reproduction = [float(row[primary_ratio_field]) for row in price_rows]
        raw_reproduction = [
            float(row["reproduction_ratio_B_over_E"]) for row in price_rows
        ]
        if not all(math.isfinite(value) and value > 0.0 for value in reproduction):
            raise RuntimeError("Nonfinite or nonpositive reproduction on the calibrated schedule")
        ordered = sorted(price_rows, key=lambda row: float(row["asset_price"]))
        root_bracket: tuple[dict[str, Any], dict[str, Any]] | None = None
        for left, right in zip(ordered[:-1], ordered[1:]):
            left_residual = float(left[primary_ratio_field]) - 1.0
            right_residual = float(right[primary_ratio_field]) - 1.0
            if left_residual == 0.0 or left_residual * right_residual <= 0.0:
                root_bracket = (left, right)
                break

        root_row: dict[str, Any] | None = None
        if root_bracket is not None:
            row_cache = {round(float(row["asset_price"]), 12): row for row in price_rows}

            def evaluate_root(asset_price: float) -> float:
                key = round(float(asset_price), 12)
                if key not in row_cache:
                    if args.fiscal_convention == "rebated-1pct":
                        sol, parameters, price, elapsed, transfer, fiscal_evaluations = solve_fiscally_balanced(
                            chain,
                            base,
                            asset_price=float(asset_price),
                        )
                    else:
                        sol, parameters, price, elapsed = solve_pe(
                            chain,
                            base,
                            asset_price=float(asset_price),
                            psi_child=float(theta["psi_child"]),
                        )
                        transfer, fiscal_evaluations = 0.0, 1
                    row_cache[key] = readout(
                        chain,
                        sol,
                        parameters,
                        price,
                        label=f"{args.profile}_closed_root_search",
                        price_ratio=float(asset_price / price0[0]),
                        psi_child=float(theta["psi_child"]),
                        elapsed=elapsed,
                    )
                    row_cache[key].update(
                        fiscal_convention=args.fiscal_convention,
                        lump_sum_transfer_period_units=float(transfer),
                        fiscal_root_evaluations=int(fiscal_evaluations),
                        fiscal_residual=float(
                            getattr(sol, "property_tax_budget_residual", math.nan)
                        ),
                    )
                    write_csv(args.outdir / "root_search.csv", list(row_cache.values()))
                return float(row_cache[key][primary_ratio_field])

            left, right = root_bracket
            increasing = float(right[primary_ratio_field]) > float(
                left[primary_ratio_field]
            )
            root_price = bisection(
                evaluate_root,
                lower=float(left["asset_price"]),
                upper=float(right["asset_price"]),
                target=1.0,
                increasing=increasing,
                tolerance=5e-7,
                label="calibrated_closed_price",
            )
            evaluate_root(root_price)
            root_row = dict(row_cache[round(float(root_price), 12)])
            root_row["static_supply_population_scale"] = float(root_row["scale_mapping_Hs_over_D"])
            root_row["fixed_baseline_stock_population_scale"] = float(
                baseline["housing_supply"] / root_row["housing_demand_per_adult"]
            )
            write_csv(args.outdir / "closed_root.csv", [root_row])

        summary = {
            "status": "complete",
            "profile": args.profile,
            "fiscal_convention": args.fiscal_convention,
            "source": str(args.source),
            "baseline": baseline,
            "primary_renewal_accounting": primary_ratio_field,
            "price_ratios": [float(row["price_ratio"]) for row in price_rows],
            "B_over_E_min": min(reproduction),
            "B_over_E_max": max(reproduction),
            "B_over_E_decreases_with_price": monotone(reproduction, increasing=False),
            "raw_state_B_over_E_min": min(raw_reproduction),
            "raw_state_B_over_E_max": max(raw_reproduction),
            "closed_root_found": root_row is not None,
            "closed_root": root_row,
            "baseline_full_retention_outside_share": 1.0
            - float(baseline[primary_ratio_field]),
            "retention_that_nests_16p9_outside_share": (1.0 - 0.169)
            / float(baseline[primary_ratio_field]),
            "local_elasticity_B_over_E_to_user_cost": log_slope(
                price_rows, primary_ratio_field
            ),
            "local_elasticity_scale_map_to_user_cost": log_slope(
                price_rows, "scale_mapping_Hs_over_D"
            ),
            "maximum_absolute_fiscal_residual": max(
                abs(float(row["fiscal_residual"])) for row in [baseline, *price_rows]
            ),
            "runtime_contract": {"J": 17, "Nb": args.nb, "max_iter_eq": 40, "tol_eq": 2.5e-5},
        }
        write_json(args.outdir / "schedule_summary.json", summary)
        plot_title = (
            "Stationary schedules in the sequential model"
            if args.profile == "e5f-floor"
            else "Stationary fertility and housing schedules"
        )
        make_schedule_plot(args.outdir, price_rows, title=plot_title)
        print("SCHEDULE_COMPLETE", flush=True)
        return

    cache: dict[tuple[float, float], dict[str, Any]] = {}

    def evaluate(psi: float, asset_price: float, label: str) -> dict[str, Any]:
        key = (round(float(psi), 12), round(float(asset_price), 12))
        if key not in cache:
            sol, parameters, price, elapsed = solve_pe(
                chain, base, asset_price=float(asset_price), psi_child=float(psi)
            )
            cache[key] = readout(
                chain,
                sol,
                parameters,
                price,
                label=label,
                price_ratio=float(asset_price / price0[0]),
                psi_child=float(psi),
                elapsed=elapsed,
            )
        return cache[key]

    psi_rep = bisection(
        lambda psi: float(evaluate(psi, float(price0[0]), "psi_root_search")["reproduction_ratio_B_over_E"]),
        lower=-3.0,
        upper=3.0,
        target=1.0,
        increasing=True,
        label="psi_replacement",
    )
    normalized_rows: list[dict[str, Any]] = []
    for index, ratio in enumerate(PRICE_RATIOS, start=1):
        row = dict(evaluate(psi_rep, float(price0[0] * ratio), "replacement_normalized_grid"))
        normalized_rows.append(row)
        write_csv(args.outdir / "replacement_normalized_price_grid.csv", normalized_rows)
        checkpoint(args.outdir, "replacement_normalized_price_grid", normalized_rows)
        print(
            f"NORMALIZED_CASE {index}/{len(PRICE_RATIOS)} ratio={ratio:.4g} "
            f"B/E={row['reproduction_ratio_B_over_E']:.9f}",
            flush=True,
        )

    normalized_at_old_price = min(normalized_rows, key=lambda row: abs(float(row["price_ratio"]) - 1.0))
    normalized_at_floor = min(normalized_rows, key=lambda row: float(row["price_ratio"]))
    feedback_capacity = float(normalized_at_floor["reproduction_ratio_B_over_E"]) - float(
        normalized_at_old_price["reproduction_ratio_B_over_E"]
    )
    if feedback_capacity <= 0.0:
        raise RuntimeError(f"Replacement-normalized price feedback is not positive: {feedback_capacity}")
    constructed_shortfall = 0.5 * feedback_capacity
    target_ratio_at_old_price = 1.0 - constructed_shortfall
    psi_shock = bisection(
        lambda psi: float(evaluate(psi, float(price0[0]), "shock_psi_search")["reproduction_ratio_B_over_E"]),
        lower=-3.0,
        upper=psi_rep,
        target=target_ratio_at_old_price,
        increasing=True,
        label="psi_half_capacity_shock",
    )

    shock_rows: list[dict[str, Any]] = []
    for index, ratio in enumerate(PRICE_RATIOS, start=1):
        row = dict(evaluate(psi_shock, float(price0[0] * ratio), "half_capacity_shock_grid"))
        shock_rows.append(row)
        write_csv(args.outdir / "half_capacity_shock_price_grid.csv", shock_rows)
        checkpoint(args.outdir, "half_capacity_shock_price_grid", shock_rows)
        print(
            f"SHOCK_CASE {index}/{len(PRICE_RATIOS)} ratio={ratio:.4g} "
            f"B/E={row['reproduction_ratio_B_over_E']:.9f}",
            flush=True,
        )

    p_floor = float(price0[0] * min(PRICE_RATIOS))
    new_root_price = bisection(
        lambda asset_price: float(evaluate(psi_shock, asset_price, "shock_price_root_search")["reproduction_ratio_B_over_E"]),
        lower=p_floor,
        upper=float(price0[0]),
        target=1.0,
        increasing=False,
        label="new_price_root",
    )
    old_endpoint = evaluate(psi_rep, float(price0[0]), "old_constructed_endpoint")
    price_only = evaluate(psi_rep, new_root_price, "price_only_counterfactual")
    new_endpoint = evaluate(psi_shock, new_root_price, "new_constructed_endpoint")
    old_scale = float(old_endpoint["scale_mapping_Hs_over_D"])
    price_scale = float(price_only["scale_mapping_Hs_over_D"])
    new_scale = float(new_endpoint["scale_mapping_Hs_over_D"])
    constructed = {
        "interpretation": (
            "Diagnostic only: old psi is chosen to impose replacement at the maintained price; "
            "the taste shock uses one half of the maximum reproductive feedback observed between "
            "one percent and one hundred percent of that price."
        ),
        "tested_price_floor_ratio": min(PRICE_RATIOS),
        "replacement_normalized_psi_child": psi_rep,
        "maintained_calibrated_psi_child": float(theta["psi_child"]),
        "psi_normalization_change": psi_rep - float(theta["psi_child"]),
        "replacement_normalized_tfr": float(old_endpoint["tfr_topcoded"]),
        "replacement_normalized_childless_rate": float(old_endpoint["childless_rate"]),
        "max_reproduction_feedback_over_tested_range": feedback_capacity,
        "constructed_reproduction_shortfall_at_old_price": constructed_shortfall,
        "postshock_psi_child": psi_shock,
        "new_root_asset_price": new_root_price,
        "new_root_price_ratio": new_root_price / float(price0[0]),
        "old_scale_mapping": old_scale,
        "price_only_scale_mapping": price_scale,
        "new_scale_mapping": new_scale,
        "price_effect_log_scale": math.log(price_scale / old_scale),
        "direct_composition_effect_log_scale": math.log(new_scale / price_scale),
        "total_log_scale_change": math.log(new_scale / old_scale),
        "new_to_old_population_scale_ratio": new_scale / old_scale,
        "old_endpoint": old_endpoint,
        "price_only_counterfactual": price_only,
        "new_endpoint": new_endpoint,
    }
    write_json(args.outdir / "constructed_normalized_closure.json", constructed)

    architecture_rows = [dict(baseline)]
    architecture_rows[0]["label"] = "maintained_shared_clock"
    if not args.skip_architecture:
        helper_json = args.outdir / "architecture_helper.json"
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--architecture-only",
            "--source",
            str(args.source),
            "--architecture-json",
            str(helper_json),
            "--nb",
            str(args.nb),
        ]
        subprocess.run(command, check=True, cwd=ROOT)
        architecture_rows.extend(json.loads(helper_json.read_text(encoding="utf-8")))
        helper_json.unlink()
    write_csv(args.outdir / "architecture_sensitivity.csv", architecture_rows)

    target_rows = [dict(row) for row in winner.get("target_fit", [])]
    if target_rows:
        write_csv(args.outdir / "source_target_fit_full.csv", target_rows)
    parameter_table = parameter_rows(theta, metadata)
    if parameter_table:
        write_csv(args.outdir / "source_parameter_table_full.csv", parameter_table)

    raw_values = [float(row["reproduction_ratio_B_over_E"]) for row in price_rows]
    adjusted_values = [float(row["topcode_adjusted_B_over_E_diagnostic"]) for row in price_rows]
    f_values = [float(row["scale_mapping_Hs_over_D"]) for row in price_rows]
    summary = {
        "status": "complete",
        "source": str(args.source),
        "source_loss": float(winner["rank_loss"]),
        "source_price": float(winner["price"]),
        "source_strict_converged": bool(winner["strict_converged"]),
        "baseline": baseline,
        "maintained_grid": {
            "price_ratios": list(PRICE_RATIOS),
            "raw_B_over_E_min": min(raw_values),
            "raw_B_over_E_max": max(raw_values),
            "topcode_adjusted_B_over_E_max": max(adjusted_values),
            "raw_reproduction_decreases_with_price": monotone(raw_values, increasing=False),
            "scale_mapping_increases_with_price": monotone(f_values, increasing=True),
            "raw_reproduction_root_found": min(raw_values) <= 1.0 <= max(raw_values),
            "topcode_adjusted_root_found": min(adjusted_values) <= 1.0 <= max(adjusted_values),
            "local_elasticity_B_over_E_to_user_cost": log_slope(price_rows, "reproduction_ratio_B_over_E"),
            "local_elasticity_scale_map_to_user_cost": log_slope(price_rows, "scale_mapping_Hs_over_D"),
        },
        "constructed_normalization": constructed,
        "architecture_sensitivity": architecture_rows,
        "runtime_contract": {"J": 17, "Nb": args.nb, "max_iter_eq": 40, "tol_eq": 2.5e-5},
    }
    write_json(args.outdir / "summary.json", summary)
    make_plots(args.outdir, price_rows, normalized_rows, shock_rows, architecture_rows, constructed)
    print("AUDIT_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
