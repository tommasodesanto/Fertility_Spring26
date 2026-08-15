#!/usr/bin/env python3
"""Calendar-time population transition for the circulated one-shot model.

This is deliberately a temporary-equilibrium first pass.  At each date the
household Bellman problem is solved as if the current price and primitives were
permanent, while the unnormalised cohort distribution is carried forward one
calendar period with the maintained one-shot transition kernels.  It is not a
perfect-foresight asset-pricing solution.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import shlex
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from intergen_housing_fertility_optimized import solver as model
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_OUTPUT = ROOT / "output/model/dynamic_population_transition/current_model"
DEFAULT_INITIAL_STATE_PACKET = ROOT / "output/model/current_one_shot_transition_initial_state"
SOLVER_PACKAGE = "intergen_housing_fertility_optimized"

# Exact selected-model values independently reproduced from the Jul-23
# source.  These are source-contract gates, not calibration targets.
EXPECTED_ENTRY_E = 0.061733455656451254
EXPECTED_MATURE_B = 0.060701144771170716
EXPECTED_B_OVER_E = EXPECTED_MATURE_B / EXPECTED_ENTRY_E
FLOW_GATE_ATOL = 5e-10
NESTING_GATE_L1 = 5e-10
MASS_GATE_ATOL = 2e-8
FEASIBILITY_PROJECTION_GATE = 1e-6
DISTRIBUTION_NEGATIVE_ATOL = 1e-14


@dataclass
class TransitionMaps:
    lmm_idx: np.ndarray
    lmm_wt: np.ndarray
    tmx_idx: np.ndarray
    tmx_wt: np.ndarray


@dataclass
class PolicyBundle:
    V: np.ndarray
    c_pol: np.ndarray
    hR_pol: np.ndarray
    bp_pol: np.ndarray
    tenure_choice: np.ndarray
    tenure_probs: np.ndarray | None
    loc_probs: np.ndarray
    fert_probs: np.ndarray
    fert_value: np.ndarray
    price: np.ndarray
    maps: TransitionMaps


@dataclass
class PeriodEvaluation:
    policy: PolicyBundle
    g_pre: np.ndarray
    g_post_fertility: np.ndarray
    g_current: np.ndarray
    births: float
    births_by_loc: np.ndarray
    demand_by_loc: np.ndarray
    supply_by_loc: np.ndarray
    relative_market_residual: float
    feasibility_projection_mass: float


class SolveCounter:
    def __init__(self) -> None:
        self.bellman = 0
        self.stationary = 0

    @property
    def total(self) -> int:
        return self.bellman + self.stationary


@dataclass(frozen=True)
class HousingSupplyRule:
    """Date-0-normalized aggregate housing supply for the transition."""

    mode: str
    initial_price: float
    initial_stock: float
    elasticity: float

    def quantity(self, price: np.ndarray) -> np.ndarray:
        value = float(np.asarray(price, dtype=float).reshape(-1)[0])
        if self.mode == "fixed-stock":
            return np.array([self.initial_stock])
        if self.mode == "static-elastic":
            return np.array(
                [self.initial_stock * (value / self.initial_price) ** self.elasticity]
            )
        raise ValueError(f"Unsupported housing-supply mode: {self.mode}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--periods", type=int, default=30)
    parser.add_argument(
        "--initialization",
        choices=("observed-age", "stationary"),
        default="stationary",
    )
    parser.add_argument(
        "--initial-state-packet",
        type=Path,
        default=DEFAULT_INITIAL_STATE_PACKET,
        help="Directory or NPZ containing g0_observed_age_pre_fertility and sibling metadata.json.",
    )
    parser.add_argument(
        "--closure",
        choices=("closed", "nesting", "fixed-entry", "all"),
        default="all",
        help="Adult-entry accounting rule. `all` runs the three rules from the same initial state.",
    )
    parser.add_argument(
        "--housing-supply",
        choices=("fixed-stock", "static-elastic"),
        default="fixed-stock",
        help="Fixed stock is the main transition; static-elastic is a diagnostic.",
    )
    parser.add_argument("--renewal-retention", type=float, default=1.0)
    parser.add_argument(
        "--entrant-conversion-factor",
        type=float,
        default=1.0,
        help="Mature-child units converted into entrant-household units; maintained default is one.",
    )
    parser.add_argument("--market-tol", type=float, default=5e-5)
    parser.add_argument("--market-max-iter", type=int, default=18)
    parser.add_argument("--root-audit", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run the exact J=17, Nb=120 model for at most four calendar periods and skip the root audit.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def distribution_health(arrays: dict[str, np.ndarray]) -> dict[str, Any]:
    """Return finite/nonnegative diagnostics for calendar-time mass arrays."""

    nonfinite_by_array: dict[str, int] = {}
    minimum_by_array: dict[str, float | None] = {}
    for name, raw in arrays.items():
        values = np.asarray(raw, dtype=float)
        finite = np.isfinite(values)
        nonfinite_by_array[name] = int(values.size - np.count_nonzero(finite))
        minimum_by_array[name] = (
            float(np.min(values[finite])) if np.any(finite) else None
        )
    finite_minima = [value for value in minimum_by_array.values() if value is not None]
    return {
        "nonfinite_distribution_count": int(sum(nonfinite_by_array.values())),
        "min_distribution_mass": min(finite_minima) if finite_minima else None,
        "nonfinite_by_array": nonfinite_by_array,
        "minimum_by_array": minimum_by_array,
    }


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_selected_source(path: Path) -> tuple[dict[str, Any], dict[str, Any], str]:
    raw = path.read_bytes()
    payload = json.loads(raw)
    if payload.get("profile") != NEW_MOMENT_PROFILE_NAME:
        raise ValueError(
            f"Expected profile {NEW_MOMENT_PROFILE_NAME!r}; found {payload.get('profile')!r}."
        )
    selected = payload.get("selected")
    if not isinstance(selected, dict) or selected.get("status") != "ok":
        raise ValueError(f"No successful selected result in {path}")
    theta = selected.get("theta")
    if not isinstance(theta, dict) or not theta:
        raise ValueError(f"Selected result has no theta in {path}")
    return payload, selected, hashlib.sha256(raw).hexdigest()


def build_transition_maps(
    price: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> TransitionMaps:
    price = np.asarray(price, dtype=float).reshape(-1)
    nb = len(b_grid)
    nt = 1 + int(P.n_house)
    hc = np.zeros((P.I, nt))
    he = np.zeros((P.I, nt))
    for i in range(P.I):
        for tenure in range(1, nt):
            size = float(P.H_own[tenure - 1])
            hc[i, tenure] = price[i] * size
            he[i, tenure] = (1.0 - float(P.psi)) * price[i] * size
    lmm_idx = np.zeros((P.I, nt, nb), dtype=np.int64)
    lmm_wt = np.zeros((P.I, nt, nb))
    for i in range(P.I):
        for tenure in range(nt):
            wealth = np.clip(b_grid + he[i, tenure], b_grid[0], b_grid[-1])
            lmm_idx[i, tenure], lmm_wt[i, tenure] = model.interp_indices(b_grid, wealth)
    tmx_idx, tmx_wt = model.build_forward_tenure_transition_maps(
        P,
        b_grid,
        hc,
        he,
        shared.phi_choice,
        shared.birth_dp,
        shared.birth_entry_grant,
    )
    return TransitionMaps(lmm_idx, lmm_wt, tmx_idx, tmx_wt)


def policy_from_solution(
    solution: SimpleNamespace,
    price: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> PolicyBundle:
    return PolicyBundle(
        solution.V,
        solution.c_pol,
        solution.hR_pol,
        solution.bp_pol,
        solution.tenure_choice,
        solution.tenure_probs,
        solution.loc_probs,
        solution.fert_probs,
        solution.fert_value,
        np.asarray(price, dtype=float).reshape(-1).copy(),
        build_transition_maps(price, P, b_grid, shared),
    )


def solve_policy(
    price: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    counter: SolveCounter,
) -> PolicyBundle:
    price = np.asarray(price, dtype=float).reshape(-1)
    objects = model.solve_bellman_full_markov_income(
        P.user_cost_rate * price, price, P, b_grid, shared
    )
    counter.bellman += 1
    V, c_pol, hR_pol, bp_pol, tenure_choice, tenure_probs, loc_probs, fert_probs, fert_value, _ = objects
    return PolicyBundle(
        V,
        c_pol,
        hR_pol,
        bp_pol,
        tenure_choice,
        tenure_probs,
        loc_probs,
        fert_probs,
        fert_value,
        price.copy(),
        build_transition_maps(price, P, b_grid, shared),
    )


def entrant_cohort(
    entry_by_loc: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
) -> np.ndarray:
    z_grid, z_weights, _ = model.income_transition_values(P)
    nb = len(b_grid)
    nt = 1 + int(P.n_house)
    out = np.zeros((nb, nt, P.I, len(z_grid), P.n_parity, P.n_child_states))
    for i in range(P.I):
        for zz, z_value in enumerate(z_grid):
            idx, weights = model.entry_wealth_grid_weights(
                b_grid, P, i=i, j=0, z_value=float(z_value)
            )
            for node, weight in zip(idx, weights):
                out[int(node), 0, i, zz, 0, 0] += (
                    float(entry_by_loc[i]) * float(z_weights[zz]) * float(weight)
                )
    return out


def gate_pre_fertility_distribution(
    g_pre: np.ndarray,
    policy: PolicyBundle,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> tuple[np.ndarray, float]:
    out = np.asarray(g_pre, dtype=float).copy()
    projected_mass = 0.0
    for j in range(P.J):
        # A surprise price change can leave an inherited debt node just beyond
        # the new Bellman feasibility frontier.  Use the core model's own
        # upward-wealth frontier projection for those inherited nodes, then
        # retain the strict dead-mass gate as a backstop.  The projected mass
        # is reported and hard-gated in the transition output.
        projected_mass += model._censor_entry_dead_mass(
            out[:, :, :, j, :, :, :], policy.V[:, :, :, j, :, :, :]
        )
        model._gate_dead_mass_at_age(
            out[:, :, :, j, :, :, :],
            policy.V[:, :, :, j, :, :, :],
            j,
            f"calendar_period_age_{P.age_start + j * P.da:g}",
            P.user_cost_rate * policy.price,
            policy.price,
            P,
            b_grid,
            shared,
            markov_income=True,
        )
    return out, projected_mass


def apply_fertility(
    g_pre: np.ndarray,
    fert_probs: np.ndarray,
    P: SimpleNamespace,
) -> tuple[np.ndarray, float, np.ndarray]:
    out = g_pre.copy()
    births = 0.0
    births_by_loc = np.zeros(P.I)
    parity_values = np.arange(P.n_parity).reshape(1, 1, 1, P.n_parity)
    for j in range(P.J):
        if not (P.A_f_start <= j + 1 <= P.A_f_end):
            continue
        for zz in range(g_pre.shape[4]):
            childless = g_pre[:, :, :, j, zz, 0, 0]
            probabilities = fert_probs[:, :, :, j, zz, :]
            split = childless[:, :, :, None] * probabilities
            births_cell = np.sum(split * parity_values, axis=3)
            births += float(np.sum(births_cell))
            births_by_loc += np.sum(births_cell, axis=(0, 1))
            out[:, :, :, j, zz, 0, 0] = 0.0
            out[:, :, :, j, zz, 0, 0] += split[:, :, :, 0]
            out[:, :, :, j, zz, 1:, 1] += split[:, :, :, 1:]
    return out, births, births_by_loc


def housing_demand_by_location(
    g_current: np.ndarray,
    hR_pol: np.ndarray,
    P: SimpleNamespace,
) -> np.ndarray:
    demand = np.zeros(P.I)
    nt = 1 + int(P.n_house)
    for i in range(P.I):
        demand[i] += float(
            np.sum(g_current[:, 0, i, :, :, :, :] * hR_pol[:, 0, i, :, :, :, :])
        )
        for tenure in range(1, nt):
            demand[i] += float(np.sum(g_current[:, tenure, i, :, :, :, :])) * float(
                P.H_own[tenure - 1]
            )
    return demand


def evaluate_period(
    price: np.ndarray,
    g_pre: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    counter: SolveCounter,
    supply_rule: HousingSupplyRule | None = None,
    supplied_policy: PolicyBundle | None = None,
) -> PeriodEvaluation:
    policy = supplied_policy or solve_policy(price, P, b_grid, shared, counter)
    gated, projected_mass = gate_pre_fertility_distribution(g_pre, policy, P, b_grid, shared)
    g_post, births, births_by_loc = apply_fertility(gated, policy.fert_probs, P)
    g_current = model.realize_current_cross_section(
        g_post,
        policy.loc_probs,
        policy.tenure_choice,
        policy.tenure_probs,
        policy.maps.lmm_idx,
        policy.maps.lmm_wt,
        policy.maps.tmx_idx,
        policy.maps.tmx_wt,
        use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
    )
    demand = housing_demand_by_location(g_current, policy.hR_pol, P)
    supply = (
        supply_rule.quantity(policy.price)
        if supply_rule is not None
        else P.H0 * ((P.user_cost_rate * policy.price) / P.r_bar) ** P.xi_supply
    )
    residual = float(np.max(np.abs(demand - supply) / np.maximum(supply, 1e-12)))
    return PeriodEvaluation(
        policy,
        gated,
        g_post,
        g_current,
        births,
        births_by_loc,
        demand,
        supply,
        residual,
        projected_mass,
    )


def clear_scalar_housing_market(
    g_pre: np.ndarray,
    price_guess: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    counter: SolveCounter,
    market_tol: float,
    max_iter: int,
    supply_rule: HousingSupplyRule,
    initial_policy: PolicyBundle | None = None,
) -> PeriodEvaluation:
    if P.I != 1:
        raise NotImplementedError("The calendar transition currently supports the circulated one-market model only.")

    def full(price: float, supplied: PolicyBundle | None = None) -> PeriodEvaluation:
        return evaluate_period(
            np.array([price]),
            g_pre,
            P,
            b_grid,
            shared,
            counter,
            supply_rule=supply_rule,
            supplied_policy=supplied,
        )

    def signed(evaluation: PeriodEvaluation) -> float:
        return float(evaluation.demand_by_loc[0] - evaluation.supply_by_loc[0])

    first = full(float(price_guess), initial_policy)
    if first.relative_market_residual <= market_tol:
        return first
    first_excess = signed(first)
    del first

    p_min = max(float(getattr(P, "p_min", 1e-4)), 1e-8)
    p_max = float(getattr(P, "p_max", 100.0))
    lo = hi = float(np.clip(price_guess, p_min, p_max))
    ex_lo = ex_hi = first_excess
    expansion = 1.18
    for _ in range(18):
        if ex_lo * ex_hi <= 0.0:
            break
        if first_excess > 0.0:
            hi = min(p_max, hi * expansion)
            candidate = full(hi)
            ex_hi = signed(candidate)
        else:
            lo = max(p_min, lo / expansion)
            candidate = full(lo)
            ex_lo = signed(candidate)
        if candidate.relative_market_residual <= market_tol:
            return candidate
        del candidate
    if ex_lo * ex_hi > 0.0:
        raise RuntimeError(
            f"Could not bracket the contemporaneous housing-market root from p={price_guess:.8g}; "
            f"excess endpoints=({ex_lo:.6g},{ex_hi:.6g})."
        )

    midpoint = 0.5 * (lo + hi)
    for _ in range(max_iter):
        midpoint = 0.5 * (lo + hi)
        candidate = full(midpoint)
        excess = signed(candidate)
        if candidate.relative_market_residual <= market_tol:
            return candidate
        del candidate
        if ex_lo * excess <= 0.0:
            hi, ex_hi = midpoint, excess
        else:
            lo, ex_lo = midpoint, excess
    final = full(midpoint)
    if final.relative_market_residual > market_tol:
        raise RuntimeError(
            f"Housing market did not clear: residual={final.relative_market_residual:.3e}."
        )
    return final


def reconstruct_stationary_pre_fertility(
    solution: SimpleNamespace,
    policy: PolicyBundle,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> tuple[np.ndarray, dict[str, float]]:
    saved_post = np.asarray(solution.g_beginning_distribution, dtype=float)
    g_pre = np.zeros_like(saved_post)
    g_pre[:, :, :, 0, :, :, :] = entrant_cohort(
        np.asarray(solution.entry_by_loc, dtype=float), P, b_grid
    )
    g_pre, projected_mass = gate_pre_fertility_distribution(g_pre, policy, P, b_grid, shared)
    _, _, Pi_z = model.income_transition_values(P)
    stochastic = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    for j in range(1, P.J):
        survival = float(P.survival_probs[j - 1]) if bool(P.use_age_survival) else 1.0
        g_pre[:, :, :, j, :, :, :] = model.advance_cohort_one_period_markov_income(
            survival * saved_post[:, :, :, j - 1, :, :, :],
            j - 1,
            policy.loc_probs,
            policy.tenure_choice,
            policy.tenure_probs,
            policy.bp_pol,
            P,
            b_grid,
            shared,
            policy.maps.lmm_idx,
            policy.maps.lmm_wt,
            policy.maps.tmx_idx,
            policy.maps.tmx_wt,
            stochastic,
            P.Pi_child if stochastic else None,
            Pi_z,
        )
    reconstructed_post, births, _ = apply_fertility(g_pre, policy.fert_probs, P)
    diagnostics = {
        "stationary_post_fertility_nesting_l1": float(np.sum(np.abs(reconstructed_post - saved_post))),
        "stationary_post_fertility_nesting_max_abs": float(np.max(np.abs(reconstructed_post - saved_post))),
        "reconstructed_baseline_births": float(births),
        "stationary_feasibility_projection_mass": float(projected_mass),
    }
    return g_pre, diagnostics


def load_initial_state_packet(
    packet: Path,
    expected_shape: tuple[int, ...],
    stationary_g_pre: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load the independently built 2023 pre-fertility initial state."""

    packet = packet.resolve()
    if packet.is_dir():
        candidates = sorted(packet.glob("*.npz"))
        observed_candidates: list[Path] = []
        stationary_candidates: list[Path] = []
        for candidate in candidates:
            with np.load(candidate, allow_pickle=False) as candidate_payload:
                if "g0_observed_age_pre_fertility" in candidate_payload.files:
                    observed_candidates.append(candidate)
                if "g0_stationary_pre_fertility" in candidate_payload.files:
                    stationary_candidates.append(candidate)
        if len(observed_candidates) != 1:
            raise ValueError(
                "Initial-state directory must contain exactly one NPZ with "
                f"g0_observed_age_pre_fertility; found {observed_candidates}."
            )
        npz_path = observed_candidates[0]
        metadata_candidates = sorted(packet.glob("*.json"))
        preferred_metadata = packet / "metadata.json"
        if preferred_metadata.is_file():
            metadata_path = preferred_metadata
        elif len([p for p in metadata_candidates if "observed" in p.name.lower()]) == 1:
            metadata_path = [p for p in metadata_candidates if "observed" in p.name.lower()][0]
        elif len(metadata_candidates) == 1:
            metadata_path = metadata_candidates[0]
        else:
            raise ValueError(
                "Initial-state directory needs metadata.json or exactly one JSON metadata file; "
                f"found {metadata_candidates}."
            )
    else:
        npz_path = packet
        metadata_path = packet.parent / "metadata.json"
        stationary_candidates = sorted(packet.parent.glob("*stationary*.npz"))
    if not npz_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(
            f"Initial-state packet requires {npz_path} and sibling {metadata_path}."
        )
    npz_raw = npz_path.read_bytes()
    metadata_raw = metadata_path.read_bytes()
    metadata = json.loads(metadata_raw)
    with np.load(npz_path, allow_pickle=False) as payload:
        key = "g0_observed_age_pre_fertility"
        if key not in payload.files:
            raise KeyError(f"Initial-state NPZ is missing {key!r}; keys={payload.files}")
        observed = np.asarray(payload[key], dtype=float).copy()
        packet_stationary = (
            np.asarray(payload["g0_stationary_pre_fertility"], dtype=float)
            if "g0_stationary_pre_fertility" in payload.files
            else None
        )
        packet_p0 = float(np.asarray(payload["date0_equilibrium_price"]).reshape(-1)[0])
        required_scalars = (
            "date0_housing_stock_K0",
            "closure_reference_E0",
            "closure_reference_B0",
            "closure_reference_M0",
        )
        missing_scalars = [name for name in required_scalars if name not in payload.files]
        if missing_scalars:
            raise KeyError(f"Initial-state NPZ is missing closure scalars: {missing_scalars}")
        packet_K0 = float(np.asarray(payload["date0_housing_stock_K0"]).reshape(-1)[0])
        closure_E0 = float(np.asarray(payload["closure_reference_E0"]).reshape(-1)[0])
        closure_B0 = float(np.asarray(payload["closure_reference_B0"]).reshape(-1)[0])
        closure_M0 = float(np.asarray(payload["closure_reference_M0"]).reshape(-1)[0])
    if packet_stationary is None and len(stationary_candidates) == 1:
        with np.load(stationary_candidates[0], allow_pickle=False) as stationary_payload:
            if "g0_stationary_pre_fertility" in stationary_payload.files:
                packet_stationary = np.asarray(
                    stationary_payload["g0_stationary_pre_fertility"], dtype=float
                )
    if observed.shape != expected_shape:
        raise ValueError(
            f"Observed initial-state shape {observed.shape} does not match {expected_shape}."
        )
    if np.any(~np.isfinite(observed)) or np.any(observed < 0.0):
        raise ValueError("Observed initial state contains nonfinite or negative mass.")
    initial_mass = float(np.sum(observed))
    if abs(initial_mass - float(np.sum(stationary_g_pre))) > 1e-8:
        raise ValueError(
            "Observed initial state must retain the stationary total-mass normalization: "
            f"observed={initial_mass:.12g}, stationary={np.sum(stationary_g_pre):.12g}."
        )
    if packet_stationary is None:
        raise ValueError("Initial-state packet must include g0_stationary_pre_fertility.")
    if packet_stationary.shape != expected_shape:
        raise ValueError(
            f"Packet stationary shape {packet_stationary.shape} does not match {expected_shape}."
        )
    stationary_packet_l1 = float(np.sum(np.abs(packet_stationary - stationary_g_pre)))
    if stationary_packet_l1 > NESTING_GATE_L1:
        raise ValueError(
            "Initial-state packet does not nest the driver's stationary pre-fertility state: "
            f"L1={stationary_packet_l1:.3e}."
        )
    acs_metadata = metadata.get("observed_age_initializer", {}).get("acs", {})
    sample_filter = str(acs_metadata.get("sample_filter", ""))
    if (
        int(acs_metadata.get("year", -1)) != 2023
        or str(acs_metadata.get("weight", "")).upper() != "HHWT"
        or "PERNUM=1" not in sample_filter
        or "RELATE=1" not in sample_filter
        or "AGE=18..85" not in sample_filter
    ):
        raise ValueError(
            "Initial-state metadata fails the 2023 household-head contract: "
            f"year={acs_metadata.get('year')}, weight={acs_metadata.get('weight')}, "
            f"sample_filter={sample_filter!r}."
        )
    details = {
        "packet": str(packet),
        "npz": str(npz_path),
        "npz_sha256": hashlib.sha256(npz_raw).hexdigest(),
        "metadata": str(metadata_path),
        "metadata_sha256": hashlib.sha256(metadata_raw).hexdigest(),
        "packet_metadata": metadata,
        "array_key": "g0_observed_age_pre_fertility",
        "array_shape": list(observed.shape),
        "initial_mass": initial_mass,
        "packet_stationary_nesting_l1": stationary_packet_l1,
        "timing": "pre_fertility",
        "date0_equilibrium_price": packet_p0,
        "date0_observed_age_housing_demand": packet_K0,
        "closure_reference_E0": closure_E0,
        "closure_reference_B0": closure_B0,
        "closure_reference_M0": closure_M0,
        "empirical_contract_gate": {
            "passed": True,
            "year": 2023,
            "weight": "HHWT",
            "head_definition": "PERNUM=1 and RELATE=1",
            "ages": "18--85 in closed four-year bins",
        },
    }
    return observed, details


def normalize_date0_housing_supply(
    g_pre: np.ndarray,
    baseline_policy: PolicyBundle,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    mode: str,
) -> tuple[HousingSupplyRule, dict[str, Any]]:
    evaluation = evaluate_period(
        baseline_policy.price,
        g_pre,
        P,
        b_grid,
        shared,
        SolveCounter(),
        supply_rule=None,
        supplied_policy=baseline_policy,
    )
    p0 = float(baseline_policy.price[0])
    demand0 = float(evaluation.demand_by_loc[0])
    original_supply0 = float(
        P.H0[0] * ((P.user_cost_rate * p0) / P.r_bar[0]) ** P.xi_supply[0]
    )
    rule = HousingSupplyRule(mode, p0, demand0, float(P.xi_supply[0]))
    normalized_supply0 = float(rule.quantity(np.array([p0]))[0])
    residual = abs(demand0 - normalized_supply0) / max(normalized_supply0, 1e-15)
    diagnostics = {
        "status": "empirical_initialization_normalization_not_recalibration",
        "housing_supply_mode": mode,
        "retained_date0_asset_price": p0,
        "date0_housing_demand": demand0,
        "date0_normalized_housing_stock": demand0,
        "original_stationary_supply_at_date0_price": original_supply0,
        "static_supply_intercept_scale": demand0 / max(original_supply0, 1e-15),
        "date0_market_residual_after_normalization": residual,
        "rule": (
            "fixed K0 at date-0 demand"
            if mode == "fixed-stock"
            else "rescale the static supply intercept so supply at p0 equals date-0 demand"
        ),
    }
    if residual > 1e-12:
        raise RuntimeError(f"Date-0 supply normalization failed: residual={residual:.3e}")
    return rule, diagnostics


def advance_calendar_distribution(
    evaluation: PeriodEvaluation,
    entry_by_loc_next: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    policy = evaluation.policy
    g_post = evaluation.g_post_fertility
    next_pre = np.zeros_like(g_post)
    mature_by_loc = np.zeros(P.I)
    deaths = 0.0
    _, _, Pi_z = model.income_transition_values(P)
    stochastic = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    K = int(P.n_child_stages)
    for j in range(P.J - 1):
        survival = float(P.survival_probs[j]) if bool(P.use_age_survival) else 1.0
        cohort = g_post[:, :, :, j, :, :, :]
        deaths += (1.0 - survival) * float(np.sum(cohort))
        survivors = survival * cohort
        next_pre[:, :, :, j + 1, :, :, :] = model.advance_cohort_one_period_markov_income(
            survivors,
            j,
            policy.loc_probs,
            policy.tenure_choice,
            policy.tenure_probs,
            policy.bp_pol,
            P,
            b_grid,
            shared,
            policy.maps.lmm_idx,
            policy.maps.lmm_wt,
            policy.maps.tmx_idx,
            policy.maps.tmx_wt,
            stochastic,
            P.Pi_child if stochastic else None,
            Pi_z,
        )
        # Reuse the exact transition kernel on the last active child stage.
        # The resulting mass in an absorbing matured state is this period's
        # newly matured child flow; multiplying by parity matches the core KFE.
        active_last = np.zeros_like(survivors)
        active_last[:, :, :, :, 1:, K] = survivors[:, :, :, :, 1:, K]
        matured_component = model.advance_cohort_one_period_markov_income(
            active_last,
            j,
            policy.loc_probs,
            policy.tenure_choice,
            policy.tenure_probs,
            policy.bp_pol,
            P,
            b_grid,
            shared,
            policy.maps.lmm_idx,
            policy.maps.lmm_wt,
            policy.maps.tmx_idx,
            policy.maps.tmx_wt,
            stochastic,
            P.Pi_child if stochastic else None,
            Pi_z,
        )
        for parity in range(1, P.n_parity):
            matured_state = K + 1 if parity == 1 else K + 2
            mature_by_loc += parity * np.sum(
                matured_component[:, :, :, :, parity, matured_state], axis=(0, 1, 3)
            )
    deaths += float(np.sum(g_post[:, :, :, P.J - 1, :, :, :]))
    next_pre[:, :, :, 0, :, :, :] = entrant_cohort(entry_by_loc_next, P, b_grid)
    expected_mass = float(np.sum(g_post)) - deaths + float(np.sum(entry_by_loc_next))
    mass_residual = float(np.sum(next_pre)) - expected_mass
    return next_pre, mature_by_loc, deaths, mass_residual


def baseline_operator_gates(
    solution: SimpleNamespace,
    policy: PolicyBundle,
    g_pre: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> dict[str, Any]:
    evaluation = evaluate_period(
        policy.price, g_pre, P, b_grid, shared, SolveCounter(), supplied_policy=policy
    )
    E = float(solution.entry_rate)
    B = float(solution.entrants_mature_total)
    outside = E - B
    outside_shares = np.asarray(P.entry_shares, dtype=float)
    provisional_next, mature_by_loc, _, mass_residual = advance_calendar_distribution(
        evaluation,
        np.zeros(P.I),
        P,
        b_grid,
        shared,
    )
    # Replace the deliberately zero age-0 cohort with the baseline-renewal flow.
    entry_by_loc = outside * outside_shares + mature_by_loc
    provisional_next[:, :, :, 0, :, :, :] = entrant_cohort(entry_by_loc, P, b_grid)
    nesting_l1 = float(np.sum(np.abs(provisional_next - g_pre)))
    mature_reconstruction = float(np.sum(mature_by_loc))
    gates = {
        "baseline_entry_flow_E": E,
        "baseline_mature_entrant_flow_B": B,
        "baseline_reproduction_B_over_E": B / E,
        "reconstructed_mature_entrant_flow_B": mature_reconstruction,
        "mature_flow_reconstruction_abs_error": abs(mature_reconstruction - B),
        "one_step_constant_path_nesting_l1": nesting_l1,
        "one_step_zero_entry_mass_residual": mass_residual,
        "expected_entry_flow_E": EXPECTED_ENTRY_E,
        "expected_mature_flow_B": EXPECTED_MATURE_B,
        "expected_reproduction_B_over_E": EXPECTED_B_OVER_E,
    }
    gates["passed"] = bool(
        abs(E - EXPECTED_ENTRY_E) <= FLOW_GATE_ATOL
        and abs(B - EXPECTED_MATURE_B) <= FLOW_GATE_ATOL
        and abs(B / E - EXPECTED_B_OVER_E) <= FLOW_GATE_ATOL
        and abs(mature_reconstruction - B) <= FLOW_GATE_ATOL
        and nesting_l1 <= NESTING_GATE_L1
    )
    return gates


def distribution_rows(
    scenario: str,
    period: int,
    g_post: np.ndarray,
    P: SimpleNamespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    total = float(np.sum(g_post))
    age_rows: list[dict[str, Any]] = []
    for j in range(P.J):
        mass = float(np.sum(g_post[:, :, :, j, :, :, :]))
        age_rows.append(
            {
                "scenario": scenario,
                "period": period,
                "years_from_start": period * float(P.period_years),
                "age": float(P.age_start + j * P.da),
                "adult_mass": mass,
                "adult_share": mass / max(total, 1e-15),
            }
        )
    K = int(P.n_child_stages)
    labels = {0: "no_children_at_home"}
    labels.update({k: f"active_child_stage_{k}" for k in range(1, K + 1)})
    labels[K + 1] = "matured_one_child_state"
    labels[K + 2] = "matured_two_plus_state"
    child_rows: list[dict[str, Any]] = []
    for parity in range(P.n_parity):
        for child_state in range(P.n_child_states):
            mass = float(np.sum(g_post[:, :, :, :, :, parity, child_state]))
            child_rows.append(
                {
                    "scenario": scenario,
                    "period": period,
                    "years_from_start": period * float(P.period_years),
                    "parity_state": parity,
                    "child_state": child_state,
                    "child_state_label": labels.get(child_state, f"state_{child_state}"),
                    "household_mass": mass,
                    "household_share": mass / max(total, 1e-15),
                    "child_units": parity * mass,
                }
            )
    return age_rows, child_rows


def run_scenario(
    label: str,
    closure: str,
    initial_g_pre: np.ndarray,
    baseline_policy: PolicyBundle,
    baseline_price: float,
    baseline_E: float,
    baseline_B: float,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    periods: int,
    retention: float,
    conversion: float,
    supply_rule: HousingSupplyRule,
    market_tol: float,
    market_max_iter: int,
    outdir: Path,
    counter: SolveCounter,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    started = time.perf_counter()
    P = copy.deepcopy(base_parameters)
    shared = model.precompute_shared(P, b_grid)
    g_pre = initial_g_pre.copy()
    outside_shares = np.asarray(P.entry_shares, dtype=float).reshape(-1)
    outside_shares = outside_shares / np.sum(outside_shares)
    if closure == "closed":
        outside_flow = 0.0
    elif closure == "nesting":
        outside_flow = baseline_E - retention * conversion * baseline_B
        if outside_flow < -FLOW_GATE_ATOL:
            raise ValueError(
                "Nesting closure implies a negative outside flow: "
                f"E={baseline_E:.8g}, rho*c*B={retention * conversion * baseline_B:.8g}."
            )
        outside_flow = max(outside_flow, 0.0)
    elif closure == "fixed-entry":
        outside_flow = None
    else:
        raise ValueError(f"Unsupported closure: {closure}")
    path_rows: list[dict[str, Any]] = []
    age_rows: list[dict[str, Any]] = []
    child_rows: list[dict[str, Any]] = []
    price_guess = baseline_price
    initial_policy = baseline_policy
    scenario_start_solves = counter.total
    max_mass_residual = 0.0
    max_market_residual = 0.0
    max_feasibility_projection = 0.0
    min_distribution_mass = math.inf
    max_nonfinite_distribution_count = 0
    for period in range(periods):
        evaluation = clear_scalar_housing_market(
            g_pre,
            price_guess,
            P,
            b_grid,
            shared,
            counter,
            market_tol,
            market_max_iter,
            supply_rule,
            initial_policy=initial_policy,
        )
        initial_policy = evaluation.policy
        price_guess = float(evaluation.policy.price[0])
        entry_flow = float(np.sum(evaluation.g_pre[:, :, :, 0, :, :, :]))

        # Obtain all incumbent transitions and the mature flow first; age-0 is
        # then filled by the renewal valve M + rho*c*B_t.
        empty_next, mature_children_by_loc, deaths, _ = advance_calendar_distribution(
            evaluation, np.zeros(P.I), P, b_grid, shared
        )
        effective_mature_by_loc = conversion * mature_children_by_loc
        retained_by_loc = retention * effective_mature_by_loc
        if closure == "fixed-entry":
            entrants_next_by_loc = baseline_E * outside_shares
        else:
            entrants_next_by_loc = outside_flow * outside_shares + retained_by_loc
        empty_next[:, :, :, 0, :, :, :] = entrant_cohort(entrants_next_by_loc, P, b_grid)
        next_pre = empty_next
        health = distribution_health(
            {
                "pre_fertility": evaluation.g_pre,
                "post_fertility": evaluation.g_post_fertility,
                "current_choice": evaluation.g_current,
                "next_pre_fertility": next_pre,
            }
        )
        period_min_mass = health["min_distribution_mass"]
        period_nonfinite_count = int(health["nonfinite_distribution_count"])
        if period_min_mass is not None:
            min_distribution_mass = min(min_distribution_mass, float(period_min_mass))
        max_nonfinite_distribution_count = max(
            max_nonfinite_distribution_count, period_nonfinite_count
        )
        if (
            period_nonfinite_count > 0
            or period_min_mass is None
            or float(period_min_mass) < -DISTRIBUTION_NEGATIVE_ATOL
        ):
            raise RuntimeError(
                f"{label} t={period}: distribution-health gate failed: {health}"
            )
        expected_next_mass = (
            float(np.sum(evaluation.g_post_fertility)) - deaths + float(np.sum(entrants_next_by_loc))
        )
        mass_residual = float(np.sum(next_pre)) - expected_next_mass
        max_mass_residual = max(max_mass_residual, abs(mass_residual))
        max_market_residual = max(max_market_residual, evaluation.relative_market_residual)
        max_feasibility_projection = max(
            max_feasibility_projection, evaluation.feasibility_projection_mass
        )

        adult_population = float(np.sum(evaluation.g_post_fertility))
        current_mass = float(np.sum(evaluation.g_current))
        owner_rate = float(np.sum(evaluation.g_current[:, 1:, :, :, :, :, :])) / max(current_mass, 1e-15)
        ages = P.age_start + np.arange(P.J) * P.da
        age_mass = np.sum(evaluation.g_post_fertility, axis=(0, 1, 2, 4, 5, 6))
        mean_age = float(np.sum(ages * age_mass) / max(np.sum(age_mass), 1e-15))
        mature_children = float(np.sum(mature_children_by_loc))
        effective_B = conversion * mature_children
        entrants_next = float(np.sum(entrants_next_by_loc))
        demand = float(np.sum(evaluation.demand_by_loc))
        supply = float(np.sum(evaluation.supply_by_loc))
        row = {
            "scenario": label,
            "period": period,
            "years_from_start": period * float(P.period_years),
            "asset_price": price_guess,
            "asset_price_index": price_guess / baseline_price,
            "housing_user_cost": float(P.user_cost_rate * price_guess),
            "adult_population": adult_population,
            "population_index": adult_population / max(float(np.sum(initial_g_pre)), 1e-15),
            "mean_adult_age": mean_age,
            "entry_flow_E": entry_flow,
            "birth_children": float(evaluation.births),
            "births_over_entry": float(evaluation.births) / max(entry_flow, 1e-15),
            "mature_children_raw": mature_children,
            "effective_mature_entrant_flow_B": effective_B,
            "mature_to_current_entry_flow_ratio_diagnostic": effective_B / max(entry_flow, 1e-15),
            "closure": closure,
            "housing_supply_mode": supply_rule.mode,
            "outside_entry_flow_M": outside_flow,
            "retained_mature_entrants": float(np.sum(retained_by_loc)),
            "entrant_flow_next": entrants_next,
            "adult_deaths": deaths,
            "housing_demand": demand,
            "housing_demand_per_adult": demand / max(adult_population, 1e-15),
            "housing_supply": supply,
            "relative_market_residual": evaluation.relative_market_residual,
            "owner_rate": owner_rate,
            "mass_accounting_residual": mass_residual,
            "feasibility_frontier_projection_mass": evaluation.feasibility_projection_mass,
            "min_distribution_mass": period_min_mass,
            "nonfinite_distribution_count": period_nonfinite_count,
        }
        path_rows.append(row)
        period_age, period_child = distribution_rows(
            label, period, evaluation.g_post_fertility, P
        )
        age_rows.extend(period_age)
        child_rows.extend(period_child)
        write_json_atomic(
            outdir / "latest_completed_period.json",
            {
                "status": "running",
                "scenario": label,
                "completed_period": period,
                "periods_requested": periods,
                "model_solve_count": counter.total,
                "latest": row,
            },
        )
        print(
            f"{label} t={period:02d} p={price_guess:.6f} pop={adult_population:.6f} "
            f"births={evaluation.births:.6f} mature={effective_B:.6f} "
            f"market={evaluation.relative_market_residual:.2e}",
            flush=True,
        )
        g_pre = next_pre

    scenario_summary = {
        "scenario": label,
        "closure": closure,
        "housing_supply_mode": supply_rule.mode,
        "outside_entry_flow_M": outside_flow,
        "renewal_retention": retention,
        "entrant_conversion_factor": conversion,
        "periods": periods,
        "model_solve_count": counter.total - scenario_start_solves,
        "elapsed_seconds": time.perf_counter() - started,
        "max_market_residual": max_market_residual,
        "max_abs_mass_accounting_residual": max_mass_residual,
        "max_feasibility_frontier_projection_mass": max_feasibility_projection,
        "min_distribution_mass": (
            min_distribution_mass if math.isfinite(min_distribution_mass) else None
        ),
        "max_nonfinite_distribution_count": max_nonfinite_distribution_count,
    }
    if max_mass_residual > MASS_GATE_ATOL:
        raise RuntimeError(
            f"{label}: mass-accounting gate failed, max residual={max_mass_residual:.3e}."
        )
    if max_feasibility_projection > FEASIBILITY_PROJECTION_GATE:
        raise RuntimeError(
            f"{label}: feasibility projection gate failed, max mass={max_feasibility_projection:.3e}."
        )
    return path_rows, age_rows, child_rows, scenario_summary


def reproduction_readout(
    price: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    counter: SolveCounter,
) -> tuple[float, dict[str, Any]]:
    solution = model.solve_markov_income_at_prices(
        np.array([price]), P, b_grid, verbose=False, fast_stats=True, SD=shared
    )
    counter.stationary += 1
    E = float(solution.entry_rate)
    B = float(solution.entrants_mature_total)
    ratio = B / max(E, 1e-15)
    row = {
        "asset_price": price,
        "housing_user_cost": float(P.user_cost_rate * price),
        "entry_flow_E": E,
        "mature_entrant_flow_B": B,
        "reproduction_B_over_E": ratio,
        "tfr": 2.0 * float(solution.mean_parity),
        "adult_population": float(solution.total_mass),
        "housing_demand_per_adult": float(np.sum(solution.housing_demand)),
        "housing_supply": float(np.sum(solution.housing_supply)),
    }
    return ratio - 1.0, row


def fixed_price_root_audit(
    label: str,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    baseline_price: float,
    observed_date0_stock: float,
    stationary_baseline_supply: float,
    counter: SolveCounter,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    shared = model.precompute_shared(P, b_grid)
    rows: list[dict[str, Any]] = []
    cache: dict[float, tuple[float, dict[str, Any]]] = {}

    def evaluate(price: float) -> tuple[float, dict[str, Any]]:
        key = round(float(price), 12)
        if key not in cache:
            cache[key] = reproduction_readout(price, P, b_grid, shared, counter)
            cache[key][1]["scenario"] = label
            cache[key][1]["price_ratio_to_current"] = price / baseline_price
        return cache[key]

    lo = max(float(getattr(P, "p_min", 1e-4)), 0.20 * baseline_price)
    hi = min(float(getattr(P, "p_max", 100.0)), 1.25 * baseline_price)
    f_lo, row_lo = evaluate(lo)
    f_hi, row_hi = evaluate(hi)
    rows.extend((dict(row_lo), dict(row_hi)))
    if f_lo * f_hi > 0.0:
        return rows, {
            "scenario": label,
            "status": "not_bracketed",
            "lower_price": lo,
            "upper_price": hi,
            "lower_B_over_E": f_lo + 1.0,
            "upper_B_over_E": f_hi + 1.0,
        }
    midpoint = 0.5 * (lo + hi)
    root_row: dict[str, Any] = {}
    for _ in range(24):
        midpoint = 0.5 * (lo + hi)
        f_mid, root_row = evaluate(midpoint)
        if abs(f_mid) <= 5e-9:
            break
        if f_lo * f_mid <= 0.0:
            hi, f_hi = midpoint, f_mid
        else:
            lo, f_lo = midpoint, f_mid
    rows.extend(dict(value[1]) for _, value in sorted(cache.items()) if value[1] not in (row_lo, row_hi))
    demand = float(root_row["housing_demand_per_adult"])
    static_supply = float(root_row["housing_supply"])
    result = {
        "scenario": label,
        "status": "root_found",
        "closed_reproduction_root_price": float(root_row["asset_price"]),
        "root_price_ratio_to_current": float(root_row["asset_price"]) / baseline_price,
        "root_reproduction_B_over_E": float(root_row["reproduction_B_over_E"]),
        "root_tfr": float(root_row["tfr"]),
        "root_housing_demand_per_adult": demand,
        "diagnostic_static_supply_population_scale": static_supply / max(demand, 1e-15),
        "fixed_observed_date0_stock_population_scale": observed_date0_stock / max(demand, 1e-15),
        "fixed_stationary_calibrated_stock_population_scale_diagnostic": (
            stationary_baseline_supply / max(demand, 1e-15)
        ),
        "static_supply_interpretation": (
            "Diagnostic endpoint under the maintained static supply curve; extreme and not a forecast."
        ),
        "fixed_stock_interpretation": (
            "Endpoint scale holding the observed-state date-0 normalized housing stock fixed."
        ),
    }
    return rows, result


def make_plots(
    path_rows: list[dict[str, Any]],
    age_rows: list[dict[str, Any]],
    child_rows: list[dict[str, Any]],
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scenarios = list(dict.fromkeys(row["scenario"] for row in path_rows))
    human_label = {
        scenario: (
            "Closed"
            if scenario.startswith("endogenous_closed")
            else "Open"
            if scenario.startswith("endogenous_open")
            else "Fixed entrants"
            if scenario.startswith("fixed_entry_benchmark")
            else scenario
        )
        for scenario in scenarios
    }
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for scenario in scenarios:
        rows = [row for row in path_rows if row["scenario"] == scenario]
        x = [row["years_from_start"] for row in rows]
        label = human_label[scenario]
        axes[0, 0].plot(x, [row["population_index"] for row in rows], marker="o", label=label)
        axes[0, 1].plot(x, [row["asset_price_index"] for row in rows], marker="o", label=label)
        axes[1, 0].plot(x, [row["entrant_flow_next"] for row in rows], marker="o", label=f"{label}: entrants")
        axes[1, 0].plot(x, [row["adult_deaths"] for row in rows], linestyle="--", label=f"{label}: deaths")
        axes[1, 1].plot(x, [row["housing_demand_per_adult"] for row in rows], marker="o", label=label)
    axes[0, 0].set(title="Adult population", ylabel="Index (initial = 1)")
    axes[0, 1].set(title="House price", ylabel="Index (date 0 = 1)")
    axes[1, 0].set(title="Adult renewal flows", ylabel="Household mass", xlabel="Years")
    axes[1, 1].set(title="Housing intensity", ylabel="Demand per adult", xlabel="Years")
    panel_values = (
        [float(row["population_index"]) for row in path_rows],
        [float(row["asset_price_index"]) for row in path_rows],
        [float(row[field]) for row in path_rows for field in ("entrant_flow_next", "adult_deaths")],
        [float(row["housing_demand_per_adult"]) for row in path_rows],
    )
    for ax, raw_values in zip(axes.flat, panel_values):
        values = np.asarray(raw_values, dtype=float)
        lower, upper = float(np.min(values)), float(np.max(values))
        span = upper - lower
        pad = 0.06 * span if span > 1e-8 else max(0.005 * max(abs(lower), 1.0), 1e-3)
        ax.set_ylim(lower - pad, upper + pad)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(outdir / "transition_paths.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(scenarios), figsize=(5 * len(scenarios), 4), squeeze=False, constrained_layout=True)
    for ax, scenario in zip(axes[0], scenarios):
        rows = [row for row in age_rows if row["scenario"] == scenario]
        periods = sorted({int(row["period"]) for row in rows})
        ages = sorted({float(row["age"]) for row in rows})
        matrix = np.array(
            [[next(row["adult_share"] for row in rows if row["period"] == t and row["age"] == age) for age in ages] for t in periods]
        )
        image = ax.imshow(matrix, aspect="auto", origin="lower", cmap="viridis")
        ax.set(title=f"Age distribution: {human_label[scenario]}", xlabel="Age index", ylabel="Calendar period")
        fig.colorbar(image, ax=ax, label="Adult share")
    fig.savefig(outdir / "adult_age_distribution.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    for scenario in scenarios:
        rows = [
            row for row in child_rows
            if row["scenario"] == scenario and row["child_state_label"].startswith("active_child_stage")
        ]
        by_period: dict[int, float] = {}
        for row in rows:
            by_period[int(row["period"])] = by_period.get(int(row["period"]), 0.0) + float(row["child_units"])
        ax.plot(sorted(by_period), [by_period[t] for t in sorted(by_period)], marker="o", label=human_label[scenario])
    ax.set(title="Children in the active pipeline", xlabel="Calendar period", ylabel="Child units")
    ax.legend(frameon=False)
    fig.savefig(outdir / "child_pipeline.png", dpi=180)
    plt.close(fig)


def write_run_readme(outdir: Path, summary: dict[str, Any]) -> None:
    command = " ".join(shlex.quote(arg) for arg in [sys.executable, *sys.argv])
    scenario_lines = "\n".join(
        f"- `{item['scenario']}`: closure `{item['closure']}`, supply `{item['housing_supply_mode']}`, "
        f"{item['periods']} periods, {item['model_solve_count']} calendar Bellman solves, "
        f"{item['elapsed_seconds']:.1f} seconds."
        for item in summary["scenario_summaries"]
    )
    root_lines = ""
    for root in summary.get("fixed_price_reproduction_root_audit", []):
        if root.get("status") == "root_found":
            root_lines += (
                f"\n- Reproductive root: price `{root['closed_reproduction_root_price']:.10f}` "
                f"(`{root['root_price_ratio_to_current']:.6f}` of date-0), "
                f"B/E `{root['root_reproduction_B_over_E']:.9f}`, TFR `{root['root_tfr']:.5f}`.\n"
                f"- Population scale at that endpoint: `{root['fixed_observed_date0_stock_population_scale']:.5f}` "
                "holding observed date-0 stock fixed; the static-supply value is diagnostic, not a forecast.\n"
            )
    text = f"""# Dynamic population transition packet

Status: **{summary['status']}**. This packet carries the selected one-shot household model through
calendar time without renormalizing population. It is a temporary-equilibrium/static-expectations
calculation, not a perfect-foresight asset-price transition.

## Exact command

```bash
{command}
```

## Initialization and scenarios

- Initial state: `{summary['initialization']}`. The observed initializer uses the independently built
  2023 ACS MMS household-head packet (HHWT; PERNUM=1 and RELATE=1; ages 18--85).
- Date-0 housing stock is normalized to observed-state demand at the retained price; this is an
  initialization convention, not a recalibration.
- Current fertility preferences are held fixed. No dated fertility-preference shock is imposed.
{scenario_lines}

## Acceptance gates and runtime

- Stationary flow diagnostics: E `{summary['baseline_entry_flow_E']:.14f}`, B
  `{summary['baseline_mature_entrant_flow_B']:.14f}`, B/E
  `{summary['baseline_reproduction_B_over_E']:.12f}`.
- Stationary reconstruction L1: `{summary['stationary_post_fertility_nesting_l1']:.3e}`;
  one-step constant-path nesting L1: `{summary['one_step_constant_path_nesting_l1']:.3e}`.
- Maximum market residual: `{summary['max_market_residual']:.3e}` (required no larger than
  `{summary['gate_tolerances']['market_relative']:.3e}`).
- Maximum mass-accounting residual: `{summary['max_abs_mass_accounting_residual']:.3e}`.
- Total model solve count, including the initial exact GE solve: `{summary['model_solve_count']}`;
  total elapsed time: `{summary['elapsed_seconds']:.1f}` seconds.
- `latest_completed_period.json` is replaced after every period. A production launcher should stop
  at its declared wall-clock budget; absence of a heartbeat for 30 minutes is unhealthy.
{root_lines}
## Interpretation limits

- Main results use the date-0 fixed housing stock. The maintained static-elastic supply curve is only
  a separately labeled diagnostic.
- Off stationarity, mature children divided by the current age-18--21 household-head stock is only a
  flow ratio, not an estimate of a demographic reproduction rate. Stationary B/E is reported only in
  the stationary/root audit.
- The selected calibration is provisional under the live empirical hold.
- The shared child clock is the one already in the one-shot model; this wrapper does not invent
  separate child birth vintages.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.periods < 1:
        raise ValueError("--periods must be positive")
    if args.smoke:
        args.periods = min(args.periods, 4)
        args.root_audit = False
    if args.renewal_retention < 0.0 or args.entrant_conversion_factor < 0.0:
        raise ValueError("Retention and conversion factors must be nonnegative.")

    started = time.perf_counter()
    source = args.source.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    payload, selected, source_sha = load_selected_source(source)
    theta = {str(key): float(value) for key, value in selected["theta"].items()}
    overrides = new_moment_overrides(tight=True, optimized=True)
    overrides.update(theta)

    print("Solving the exact selected one-shot model point...", flush=True)
    baseline_solution, base_parameters, price = model.run_model_cp_dt(overrides, verbose=False)
    if int(base_parameters.J) != 17 or int(base_parameters.Nb) != 120:
        raise RuntimeError(
            f"Selected-model dimension gate failed: J={base_parameters.J}, Nb={base_parameters.Nb}."
        )
    b_grid = np.asarray(baseline_solution.b_grid, dtype=float)
    baseline_shared = model.precompute_shared(base_parameters, b_grid)
    baseline_policy = policy_from_solution(
        baseline_solution, price, base_parameters, b_grid, baseline_shared
    )
    initial_g_pre, nesting = reconstruct_stationary_pre_fertility(
        baseline_solution,
        baseline_policy,
        base_parameters,
        b_grid,
        baseline_shared,
    )
    gates = baseline_operator_gates(
        baseline_solution,
        baseline_policy,
        initial_g_pre,
        base_parameters,
        b_grid,
        baseline_shared,
    )
    gates.update(nesting)
    gates["passed"] = bool(
        gates["passed"]
        and gates["stationary_post_fertility_nesting_l1"] <= NESTING_GATE_L1
        and abs(gates["reconstructed_baseline_births"] - float(baseline_solution.total_births_kfe))
        <= FLOW_GATE_ATOL
    )
    if not gates["passed"]:
        write_json_atomic(
            outdir / "summary.json",
            {
                "status": "failed_baseline_contract",
                "source": source,
                "source_sha256": source_sha,
                "baseline_gates": gates,
            },
        )
        raise RuntimeError(f"Selected-model baseline contract failed: {gates}")

    stationary_g_pre = initial_g_pre
    if args.initialization == "observed-age":
        initial_g_pre, initialization_details = load_initial_state_packet(
            args.initial_state_packet,
            tuple(stationary_g_pre.shape),
            stationary_g_pre,
        )
    else:
        initialization_details = {
            "source": "stationary selected-model distribution",
            "initial_mass": float(np.sum(initial_g_pre)),
            "warning": "This state has no demographic momentum and is retained only as a nesting benchmark.",
        }
    supply_rule, date0_supply_normalization = normalize_date0_housing_supply(
        initial_g_pre,
        baseline_policy,
        base_parameters,
        b_grid,
        baseline_shared,
        args.housing_supply,
    )
    if args.initialization == "observed-age":
        packet_price_error = abs(
            float(initialization_details["date0_equilibrium_price"])
            - float(np.asarray(price).reshape(-1)[0])
        )
        packet_stock_error = abs(
            float(initialization_details["date0_observed_age_housing_demand"])
            - float(supply_rule.initial_stock)
        )
        packet_E0_error = abs(
            float(initialization_details["closure_reference_E0"])
            - float(baseline_solution.entry_rate)
        )
        packet_B0_error = abs(
            float(initialization_details["closure_reference_B0"])
            - float(baseline_solution.entrants_mature_total)
        )
        packet_M0_error = abs(
            float(initialization_details["closure_reference_M0"])
            - (
                float(baseline_solution.entry_rate)
                - float(baseline_solution.entrants_mature_total)
            )
        )
        date0_supply_normalization["packet_price_abs_error"] = packet_price_error
        date0_supply_normalization["packet_stock_abs_error"] = packet_stock_error
        date0_supply_normalization["packet_E0_abs_error"] = packet_E0_error
        date0_supply_normalization["packet_B0_abs_error"] = packet_B0_error
        date0_supply_normalization["packet_M0_abs_error"] = packet_M0_error
        date0_supply_normalization["packet_contract_passed"] = bool(
            max(
                packet_price_error,
                packet_stock_error,
                packet_E0_error,
                packet_B0_error,
                packet_M0_error,
            )
            <= 1e-10
        )
        if not date0_supply_normalization["packet_contract_passed"]:
            raise RuntimeError(
                "Observed-state packet date-0 normalization does not match the driver: "
                f"price error={packet_price_error:.3e}, stock error={packet_stock_error:.3e}."
            )
    if date0_supply_normalization["date0_market_residual_after_normalization"] > 1e-12:
        raise RuntimeError("Observed-state date-0 housing normalization gate failed.")

    closures = (
        ["closed", "nesting", "fixed-entry"]
        if args.closure == "all"
        else [args.closure]
    )
    closure_labels = {
        "closed": "endogenous_closed",
        "nesting": "endogenous_open",
        "fixed-entry": "fixed_entry_benchmark",
    }
    scenarios = [
        (f"{closure_labels[closure]}_{args.housing_supply}", closure)
        for closure in closures
    ]

    counter = SolveCounter()
    all_paths: list[dict[str, Any]] = []
    all_ages: list[dict[str, Any]] = []
    all_children: list[dict[str, Any]] = []
    scenario_summaries: list[dict[str, Any]] = []
    for label, closure in scenarios:
        paths, ages, children, scenario_summary = run_scenario(
            label,
            closure,
            initial_g_pre,
            baseline_policy,
            float(np.asarray(price).reshape(-1)[0]),
            float(baseline_solution.entry_rate),
            float(baseline_solution.entrants_mature_total),
            base_parameters,
            b_grid,
            int(args.periods),
            float(args.renewal_retention),
            float(args.entrant_conversion_factor),
            supply_rule,
            float(args.market_tol),
            int(args.market_max_iter),
            outdir,
            counter,
        )
        all_paths.extend(paths)
        all_ages.extend(ages)
        all_children.extend(children)
        scenario_summaries.append(scenario_summary)

    period0_rows = [row for row in all_paths if int(row["period"]) == 0]
    period0_compare_fields = (
        "asset_price",
        "housing_user_cost",
        "adult_population",
        "mean_adult_age",
        "entry_flow_E",
        "birth_children",
        "mature_children_raw",
        "effective_mature_entrant_flow_B",
        "mature_to_current_entry_flow_ratio_diagnostic",
        "housing_demand",
        "housing_supply",
        "owner_rate",
    )
    period0_max_differences = {
        field: max(float(row[field]) for row in period0_rows)
        - min(float(row[field]) for row in period0_rows)
        for field in period0_compare_fields
    }
    period0_common_state_gate = {
        "passed": bool(max(period0_max_differences.values(), default=0.0) <= 1e-12),
        "tolerance": 1e-12,
        "max_differences": period0_max_differences,
    }
    if not period0_common_state_gate["passed"]:
        raise RuntimeError(
            f"Scenario date-0 common-state gate failed: {period0_common_state_gate}"
        )

    root_rows: list[dict[str, Any]] = []
    root_results: list[dict[str, Any]] = []
    if args.root_audit:
        rows, result = fixed_price_root_audit(
            "selected_one_shot_model",
            copy.deepcopy(base_parameters),
            b_grid,
            float(np.asarray(price).reshape(-1)[0]),
            float(supply_rule.initial_stock),
            float(np.sum(baseline_solution.housing_supply)),
            counter,
        )
        root_rows.extend(rows)
        root_results.append(result)

    write_csv(outdir / "transition_path.csv", all_paths)
    write_csv(outdir / "adult_age_distribution.csv", all_ages)
    write_csv(outdir / "child_pipeline.csv", all_children)
    if root_rows:
        write_csv(outdir / "reproduction_schedule.csv", root_rows)
    make_plots(all_paths, all_ages, all_children, outdir)

    elapsed = time.perf_counter() - started
    max_market_residual = max(s["max_market_residual"] for s in scenario_summaries)
    if max_market_residual > float(args.market_tol):
        raise RuntimeError(
            f"Production market-clearing gate failed: {max_market_residual:.3e} > {args.market_tol:.3e}."
        )
    driver_path = Path(__file__).resolve()
    solver_path = Path(model.__file__).resolve()
    summary = {
        "status": "complete",
        "method": "deterministic cohort-mass transition with temporary-equilibrium/static expectations",
        "not_perfect_foresight": True,
        "solver_package": SOLVER_PACKAGE,
        "source": source,
        "source_sha256": source_sha,
        "driver": str(driver_path),
        "driver_sha256": hashlib.sha256(driver_path.read_bytes()).hexdigest(),
        "solver_source": str(solver_path),
        "solver_source_sha256": hashlib.sha256(solver_path.read_bytes()).hexdigest(),
        "source_profile": payload.get("profile"),
        "source_selected_price": float(selected["price"]),
        "source_selected_loss": float(selected["rank_loss"]),
        "theta": theta,
        "J": int(base_parameters.J),
        "Nb": int(base_parameters.Nb),
        "period_years": float(base_parameters.period_years),
        "periods": int(args.periods),
        "initialization": args.initialization,
        "initial_state": initialization_details,
        "date0_housing_normalization": date0_supply_normalization,
        "date0_common_state_across_closures_gate": period0_common_state_gate,
        "housing_supply_mode": args.housing_supply,
        "scenarios": [label for label, _ in scenarios],
        "baseline_entry_flow_E": gates["baseline_entry_flow_E"],
        "baseline_mature_entrant_flow_B": gates["baseline_mature_entrant_flow_B"],
        "baseline_reproduction_B_over_E": gates["baseline_reproduction_B_over_E"],
        "stationary_post_fertility_nesting_l1": gates["stationary_post_fertility_nesting_l1"],
        "one_step_constant_path_nesting_l1": gates["one_step_constant_path_nesting_l1"],
        "baseline_gates": gates,
        "gate_tolerances": {
            "flow_absolute": FLOW_GATE_ATOL,
            "distribution_l1": NESTING_GATE_L1,
            "mass_accounting_absolute": MASS_GATE_ATOL,
            "feasibility_frontier_projection_mass": FEASIBILITY_PROJECTION_GATE,
            "distribution_negative_absolute": DISTRIBUTION_NEGATIVE_ATOL,
            "market_relative": float(args.market_tol),
        },
        "renewal_closures": {
            "endogenous_closed": "E_(t+1) = retention * conversion * mature_children_t; M=0",
            "endogenous_open": "E_(t+1) = M + retention * conversion * mature_children_t; M=E0-rho*c*B0",
            "fixed_entry_benchmark": "E_(t+1) = E0; accounting benchmark isolating endogenous renewal",
            "renewal_retention": float(args.renewal_retention),
            "entrant_conversion_factor": float(args.entrant_conversion_factor),
        },
        "expectations_and_price_assumptions": {
            "household_expectations": "current price and current primitives treated as permanent at each date",
            "asset_price_rule": "r_t = (q + delta + tau_H) p_t; no transition capital-gain term",
            "housing_market": (
                "contemporaneous demand equals fixed date-0 stock K0"
                if args.housing_supply == "fixed-stock"
                else "contemporaneous demand equals date-0-normalized static-elastic supply (diagnostic)"
            ),
            "population_normalization": "none after initialization",
        },
        "scenario_summaries": scenario_summaries,
        "fixed_price_reproduction_root_audit": root_results,
        "model_solve_count": counter.total + 1,
        "model_solve_count_breakdown": {
            "calendar_bellman_solves": counter.bellman,
            "stationary_root_audit_solves": counter.stationary,
            "initial_exact_source_ge_solve": 1,
        },
        "elapsed_seconds": elapsed,
        "max_market_residual": max_market_residual,
        "max_abs_mass_accounting_residual": max(
            s["max_abs_mass_accounting_residual"] for s in scenario_summaries
        ),
        "max_feasibility_frontier_projection_mass": max(
            s["max_feasibility_frontier_projection_mass"] for s in scenario_summaries
        ),
        "min_distribution_mass": min(
            float(s["min_distribution_mass"])
            for s in scenario_summaries
            if s["min_distribution_mass"] is not None
        ),
        "max_nonfinite_distribution_count": max(
            s["max_nonfinite_distribution_count"] for s in scenario_summaries
        ),
        "artifacts": [
            "transition_path.csv",
            "adult_age_distribution.csv",
            "child_pipeline.csv",
            "transition_paths.png",
            "adult_age_distribution.png",
            "child_pipeline.png",
            "README.md",
        ] + (["reproduction_schedule.csv"] if root_rows else []),
        "caveats": [
            "The selected calibration remains diagnostic/provisional under the live empirical hold.",
            "Current fertility preferences are held fixed; the exercise does not invent a dated preference shock.",
            "The transition is not a perfect-foresight solution for a sequence of asset prices.",
        ],
    }
    write_json_atomic(outdir / "summary.json", summary)
    write_run_readme(outdir, summary)
    write_json_atomic(
        outdir / "latest_completed_period.json",
        {
            "status": "complete",
            "periods": int(args.periods),
            "scenarios": [label for label, _ in scenarios],
            "model_solve_count": counter.total + 1,
            "elapsed_seconds": elapsed,
            "summary": str(outdir / "summary.json"),
        },
    )
    print(
        f"complete: periods={args.periods}, scenarios={len(scenarios)}, "
        f"model_solves={counter.total + 1}, elapsed={elapsed:.1f}s, output={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
