#!/usr/bin/env python3
"""Fixed-price numerical asset-grid diagnosis for the selected soft reference.

This driver changes only the wealth grid. It does not recalibrate parameters,
solve a price root, or change the economic specification. It performs four
sequential one-core native solves and writes a checkpoint after each case.
"""
from __future__ import annotations

import os

# Set numerical thread limits before importing NumPy or the model.
for _name in (
    "NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import copy
import hashlib
import json
import resource
import signal
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1"
SAVED_CASE = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz"
TIME_BUDGET_SECONDS = 600.0
CASE_BUDGET_SECONDS = 180.0
MEMORY_LIMIT_BYTES = 24 * 1024**3
SUPPORT_CUTOFF = 1e-12
TAIL_CUTOFF = 33.66055
OCCUPIED_LOWER = -6.4


def _json_dump(path: Path, value: object) -> None:
    def safe(x: object) -> object:
        if isinstance(x, dict):
            return {str(k): safe(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [safe(v) for v in x]
        if isinstance(x, np.ndarray):
            return safe(x.tolist())
        if isinstance(x, np.generic):
            return safe(x.item())
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def _memory_limit() -> dict[str, object]:
    """Apply a portable address-space cap only if it does not lower current use."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        # /proc is not available on every supported platform. If unavailable,
        # leave the limit alone and record that the platform did not expose it.
        status = Path("/proc/self/status")
        current = None
        if status.exists():
            for line in status.read_text().splitlines():
                if line.startswith("VmSize:"):
                    current = int(line.split()[1]) * 1024
                    break
        if current is None or current >= MEMORY_LIMIT_BYTES:
            return {"applied": False, "reason": "current virtual size unavailable or already at/above 24 GiB"}
        target_hard = MEMORY_LIMIT_BYTES if hard == resource.RLIM_INFINITY else min(hard, MEMORY_LIMIT_BYTES)
        target_soft = MEMORY_LIMIT_BYTES if soft == resource.RLIM_INFINITY else min(target_hard, MEMORY_LIMIT_BYTES)
        if target_soft < current:
            return {"applied": False, "reason": "24 GiB cap is below current virtual size"}
        resource.setrlimit(resource.RLIMIT_AS, (target_soft, target_hard))
        return {"applied": True, "bytes": int(target_soft), "current_vms_bytes": int(current)}
    except (AttributeError, OSError, ValueError) as exc:
        return {"applied": False, "reason": f"resource limit unavailable: {type(exc).__name__}: {exc}"}


def _grid_variants(old: np.ndarray) -> dict[str, np.ndarray]:
    end_extra = np.array([3500.0, 4000.0, 5000.0, 6000.0])
    tail_midpoints = (old[:-1][old[1:] > TAIL_CUTOFF] + old[1:][old[1:] > TAIL_CUTOFF]) / 2.0
    core_mask = (old[:-1] < TAIL_CUTOFF) & (old[1:] > OCCUPIED_LOWER)
    core_midpoints = (old[:-1][core_mask] + old[1:][core_mask]) / 2.0
    return {
        "baseline120": old.copy(),
        "endpoint_extension": np.unique(np.concatenate([old, end_extra])),
        "upper_tail_refinement": np.unique(np.concatenate([old, tail_midpoints])),
        "occupied_region_refinement": np.unique(np.concatenate([old, core_midpoints])),
    }


def _followup_grid_variants(old: np.ndarray) -> tuple[dict[str, np.ndarray], float]:
    exact_core = float(old[int(np.searchsorted(old, TAIL_CUTOFF, side="left"))])
    strict_tail = old[:-1] >= exact_core
    tail_midpoints = (old[:-1][strict_tail] + old[1:][strict_tail]) / 2.0
    occupied = (old[:-1] < exact_core) & (old[1:] > OCCUPIED_LOWER)
    left, right = old[:-1][occupied], old[1:][occupied]
    quarter_points = np.concatenate([left + (right - left) / 4.0,
                                     left + (right - left) / 2.0,
                                     left + 3.0 * (right - left) / 4.0])
    return {
        "strict_upper_tail_refinement": np.unique(np.concatenate([old, tail_midpoints])),
        "occupied_region_quarter_refinement": np.unique(np.concatenate([old, quarter_points])),
    }, exact_core


def _input_fingerprint(P: SimpleNamespace) -> str:
    grid_fields = {"Nb", "b_min", "b_max", "earnings_transaction_grid",
                   "fixed_reference_entry_grid", "fixed_reference_entry_conditional"}
    payload = {key: repr(value) for key, value in sorted(vars(P).items())
               if key not in grid_fields and not key.startswith("_")}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _followup_input_fingerprint(P: SimpleNamespace) -> str:
    """Fingerprint economics while omitting known per-process evidence scratch path."""
    excluded = {"Nb", "b_min", "b_max", "earnings_transaction_grid",
                "fixed_reference_entry_grid", "fixed_reference_entry_conditional",
                "native_inherited_distribution_evidence_dir"}
    payload = {key: repr(value) for key, value in sorted(vars(P).items())
               if key not in excluded and not key.startswith("_")}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _install_grid(P: SimpleNamespace, old: np.ndarray, new: np.ndarray) -> None:
    positions = np.searchsorted(new, old)
    if np.any(positions >= len(new)) or not np.array_equal(new[positions], old):
        raise AssertionError("every original wealth-grid node must remain exactly present")
    original_conditional = np.asarray(P.fixed_reference_entry_conditional, dtype=float)
    expanded = np.zeros((len(new), original_conditional.shape[1]), dtype=float)
    expanded[positions, :] = original_conditional
    P.Nb = int(len(new))
    P.b_min = float(new[0])
    P.b_max = float(new[-1])
    P.earnings_transaction_grid = new.copy()
    P.fixed_reference_entry_grid = new.copy()
    P.fixed_reference_entry_conditional = expanded
    P.native_explicit_transaction_grid = True

    np.testing.assert_array_equal(new[positions], old)
    np.testing.assert_allclose(expanded.sum(axis=0), original_conditional.sum(axis=0), rtol=0, atol=2e-14)
    for power in (1, 2):
        old_moment = np.sum(original_conditional * old[:, None] ** power, axis=0)
        new_moment = np.sum(expanded * new[:, None] ** power, axis=0)
        np.testing.assert_allclose(new_moment, old_moment, rtol=0, atol=2e-11)
    # No entrant mass may be relocated, clipped, or censored by a grid change.
    if not np.array_equal(expanded[positions], original_conditional) or np.any(expanded.sum(axis=0) != original_conditional.sum(axis=0)):
        raise AssertionError("entrant mass changed while expanding the grid")


def _solution_arrays(sol: SimpleNamespace) -> dict[str, np.ndarray]:
    names = (
        "V", "c_pol", "bp_pol", "hR_pol", "tenure_probs", "fert_probs",
        "bp_pol_stay", "c_pol_stay", "g_beginning_distribution", "g",
        "g_stay_distribution", "b_grid",
    )
    values = {}
    for name in names:
        value = getattr(sol, name, None)
        if value is not None:
            values[name] = np.asarray(value)
    return values


def _save_solution(sol: SimpleNamespace, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Uncompressed NPZ avoids spending a large share of the short budget on I/O.
    np.savez(path, **_solution_arrays(sol))


def _exact_baseline_check(sol: SimpleNamespace, saved: SimpleNamespace) -> dict[str, object]:
    checked = (
        "V", "c_pol", "bp_pol", "hR_pol", "tenure_probs", "fert_probs",
        "bp_pol_stay", "c_pol_stay", "g", "g_stay_distribution",
        "g_beginning_distribution",
    )
    result = {}
    for name in checked:
        native = getattr(sol, name, None)
        reference = getattr(saved, name, None)
        if native is None or reference is None:
            result[name] = {"available": False}
            continue
        if np.asarray(native).shape != np.asarray(reference).shape:
            result[name] = {"available": True, "exact": False, "reason": "shape mismatch",
                            "native_shape": list(np.asarray(native).shape),
                            "saved_shape": list(np.asarray(reference).shape)}
            continue
        equal = bool(np.array_equal(native, reference, equal_nan=True))
        result[name] = {"available": True, "exact": equal}
        if not equal:
            result[name]["max_abs_difference"] = float(np.nanmax(np.abs(np.asarray(native) - np.asarray(reference))))
    return result


def _grid_node_map(old: np.ndarray, new: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(new, old)
    if np.any(idx >= new.size) or not np.array_equal(new[idx], old):
        raise AssertionError("candidate grid lost an inherited node")
    return idx


def _policy_differences(base: SimpleNamespace, candidate: SimpleNamespace,
                        old_grid: np.ndarray, new_grid: np.ndarray,
                        P: SimpleNamespace) -> dict[str, object]:
    """Compare tenure-probability expected policies at inherited beginning states."""
    mass = np.asarray(base.g_beginning_distribution, dtype=float)
    probs = np.asarray(base.tenure_probs, dtype=float)
    total = float(mass.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("baseline beginning distribution has no positive mass")
    result: dict[str, object] = {
        "weighting": "baseline g_beginning_distribution; policies average over saved tenure_probs at each inherited state",
        "population_mass": total,
        "strictly_positive_state_cells": int(np.count_nonzero(mass > 0.0)),
        "state_cells_at_or_above_1e-12_population": int(np.count_nonzero(mass >= SUPPORT_CUTOFF * total)),
        "transaction_map": "x=b when tenure is unchanged; otherwise x=b+(1-psi)*price*H_old-price*H_new, clipped as in native forward mapping",
        "grid_node_comparison": "candidate expected policies are evaluated on every original beginning-state wealth node",
        "policies": {},
    }
    # The pinned household builder defines heq=(1-psi)*price*H and hcost=price*H;
    # distribution.build_forward_tenure_transition_maps uses x=b+heq_old-hcost_new.
    price = float(P.reference_price)
    psi = float(P.psi)
    h = np.concatenate(([0.0], np.asarray(P.H_own, dtype=float).reshape(-1)))
    sale = (1.0 - psi) * price * h
    purchase = price * h
    nt = mass.shape[1]
    locs = mass.shape[2]
    if probs.shape != mass.shape + (nt,):
        raise ValueError(f"tenure probabilities do not align with beginning distribution: {probs.shape}")

    def interp_axis0(values: np.ndarray, grid: np.ndarray, points: np.ndarray) -> np.ndarray:
        clipped = np.clip(points, grid[0], grid[-1])
        idx = np.clip(np.searchsorted(grid, clipped, side="right") - 1, 0, grid.size - 2)
        wt = (clipped - grid[idx]) / (grid[idx + 1] - grid[idx])
        return (1.0 - wt.reshape((-1,) + (1,) * (values.ndim - 1))) * values[idx] + wt.reshape((-1,) + (1,) * (values.ndim - 1)) * values[idx + 1]

    common = _grid_node_map(old_grid, new_grid)
    expected: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name in ("consumption", "next_wealth", "housing_services"):
        expected[name] = (np.zeros_like(mass, dtype=float), np.zeros_like(mass, dtype=float))
    for loc in range(locs):
        for old_tenure in range(nt):
            for new_tenure in range(nt):
                candidate_probs = np.asarray(candidate.tenure_probs, dtype=float)[common, old_tenure, loc, ..., new_tenure]
                if (not np.any(probs[:, old_tenure, loc, ..., new_tenure] > 0.0)
                        and not np.any(candidate_probs > 0.0)):
                    continue
                x = old_grid.copy() if old_tenure == new_tenure else old_grid + sale[old_tenure] - purchase[new_tenure]
                mapped = np.clip(x, old_grid[0], old_grid[-1])
                mapped_candidate = np.clip(x, new_grid[0], new_grid[-1])
                for sol, grid, slot in ((base, old_grid, 0), (candidate, new_grid, 1)):
                    probabilities = np.asarray(sol.tenure_probs, dtype=float)
                    if slot == 1:
                        probabilities = probabilities[common]
                    choice_probability = probabilities[:, old_tenure, loc, ..., new_tenure]
                    c_source = np.asarray(sol.c_pol, dtype=float)[:, new_tenure, loc, ...]
                    b_source = np.asarray(sol.bp_pol, dtype=float)[:, new_tenure, loc, ...]
                    if old_tenure == new_tenure and old_tenure > 0:
                        stay_c = getattr(sol, "c_pol_stay", None)
                        stay_b = getattr(sol, "bp_pol_stay", None)
                        if stay_c is not None:
                            c_source = np.asarray(stay_c, dtype=float)[:, new_tenure, loc, ...]
                        if stay_b is not None:
                            b_source = np.asarray(stay_b, dtype=float)[:, new_tenure, loc, ...]
                    source_grid = grid
                    points = mapped if slot == 0 else mapped_candidate
                    c_at_x = interp_axis0(c_source, source_grid, points)
                    b_at_x = interp_axis0(b_source, source_grid, points)
                    if new_tenure == 0:
                        h_source = np.asarray(sol.hR_pol, dtype=float)[:, 0, loc, ...]
                        h_at_x = interp_axis0(h_source, source_grid, points)
                    else:
                        h_at_x = np.full_like(c_at_x, float(P.H_own[new_tenure - 1]))
                    for arr, value in ((expected["consumption"][slot], c_at_x),
                                       (expected["next_wealth"][slot], b_at_x),
                                       (expected["housing_services"][slot], h_at_x)):
                        arr[:, old_tenure, loc, ...] += choice_probability * value

    for name, (base_expected, candidate_expected) in expected.items():
        diff = np.abs(candidate_expected - base_expected)
        records: dict[str, object] = {}
        for label, selected in (("strictly_positive", mass > 0.0),
                                ("mass_ge_1e-12_population", mass >= SUPPORT_CUTOFF * total)):
            weights = mass[selected]
            values = diff[selected]
            if weights.size == 0 or float(weights.sum()) <= 0:
                records[label] = {"cells": int(weights.size), "mass": float(weights.sum()),
                                  "weighted_mean_abs": None, "max_abs": None, "max_state": None}
                continue
            flat = int(np.argmax(values))
            state_idx = np.unravel_index(np.flatnonzero(selected)[flat], mass.shape)
            records[label] = {
                "cells": int(weights.size), "mass": float(weights.sum()),
                "weighted_mean_abs": float(np.average(values, weights=weights)),
                "max_abs": float(values[flat]),
                "max_state": {"wealth": float(old_grid[state_idx[0]]),
                              "tenure_index": int(state_idx[1]), "location_index": int(state_idx[2]),
                              "age_income_children_childstate_indices": [int(x) for x in state_idx[3:]]},
            }
        result["policies"][name] = records
    return result


def _distribution_comparison(base: SimpleNamespace, candidate: SimpleNamespace,
                             old_grid: np.ndarray, new_grid: np.ndarray) -> dict[str, object]:
    old = np.asarray(base.g_cross_sectional_wealth_distribution, dtype=float).sum(axis=tuple(range(1, np.asarray(base.g_cross_sectional_wealth_distribution).ndim)))
    new = np.asarray(candidate.g_cross_sectional_wealth_distribution, dtype=float).sum(axis=tuple(range(1, np.asarray(candidate.g_cross_sectional_wealth_distribution).ndim)))
    old_mass = float(old.sum())
    new_mass = float(new.sum())
    if old_mass <= 0 or new_mass <= 0:
        raise ValueError("wealth distribution has no mass")
    old = old / old_mass
    new = new / new_mass
    union = np.unique(np.concatenate([old_grid, new_grid]))
    old_on_union = np.zeros(union.size, dtype=float)
    new_on_union = np.zeros(union.size, dtype=float)
    old_on_union[np.searchsorted(union, old_grid)] = old
    new_on_union[np.searchsorted(union, new_grid)] = new
    # For discrete wealth distributions, the CDF is constant between adjacent
    # support points; integrate its absolute gap exactly on the union grid.
    cdf_gap = np.cumsum(old_on_union - new_on_union)
    w1 = float(np.sum(np.abs(cdf_gap[:-1]) * np.diff(union)))
    return {"baseline_mass": old_mass, "candidate_mass": new_mass,
            "baseline_normalized_mass": float(old.sum()),
            "candidate_normalized_mass": float(new.sum()), "wealth_cdf_wasserstein1": w1}


def _aggregate_comparison(base: SimpleNamespace, candidate: SimpleNamespace,
                          P: SimpleNamespace) -> dict[str, object]:
    from model_policy_tools import aggregate_solution

    old = aggregate_solution(base, houses=P.H_own, age_start=P.age_start,
                             period_years=P.period_years)
    new = aggregate_solution(candidate, houses=P.H_own, age_start=P.age_start,
                             period_years=P.period_years)
    old_overall = old["overall"]
    new_overall = new["overall"]
    differences = {}
    for key, new_value in new_overall.items():
        old_value = old_overall.get(key)
        if isinstance(new_value, (int, float, np.number)) and isinstance(old_value, (int, float, np.number)):
            old_float, new_float = float(old_value), float(new_value)
            differences[key] = {
                "baseline": old_float, "candidate": new_float,
                "absolute_change": new_float - old_float,
                "relative_change": ((new_float - old_float) / abs(old_float)) if old_float != 0.0 else None,
            }
    return {"units": old.get("units", new.get("units")),
            "baseline_population_mass": old.get("population_mass"),
            "candidate_population_mass": new.get("population_mass"),
            "baseline_overall": old_overall,
            "candidate_overall": new_overall,
            "overall_differences": differences}


def _renter_slice_difference(base: SimpleNamespace, candidate: SimpleNamespace,
                             old_grid: np.ndarray, new_grid: np.ndarray) -> dict[str, object]:
    common = _grid_node_map(old_grid, new_grid)
    # Age 30 (index 3), income state five (index 4), renter, childless.
    selector = (0, 0, 3, 4, 0, 0)
    result = {"definition": "raw renter conditional policy; age 30, income state 5, childless",
              "max_abs_change": {}}
    for name in ("c_pol", "bp_pol", "hR_pol"):
        old_values = np.asarray(getattr(base, name), dtype=float)[(slice(None),) + selector]
        new_values = np.asarray(getattr(candidate, name), dtype=float)[(common,) + selector]
        delta = np.abs(new_values - old_values)
        idx = int(np.argmax(delta))
        result["max_abs_change"][name] = {"value": float(delta[idx]), "wealth": float(old_grid[idx]),
                                          "baseline": float(old_values[idx]), "candidate": float(new_values[idx])}
    return result


def _load_saved_npz(path: Path, price: float) -> SimpleNamespace:
    with np.load(path, allow_pickle=False) as archive:
        values = {key: archive[key].copy() for key in archive.files}
    # The pinned distribution builder assigns both diagnostics from the same
    # native `g.copy()`. v1's minimal archive omitted the alias field.
    if "g_cross_sectional_wealth_distribution" not in values:
        values["g_cross_sectional_wealth_distribution"] = values["g_beginning_distribution"].copy()
    values.update(price=float(price), timing="transaction_inside", case_id="v1_occupied_region_refinement_214")
    return SimpleNamespace(**values)


def _exact_archive_replay(path: Path, saved: SimpleNamespace) -> dict[str, object]:
    checked = ("V", "c_pol", "bp_pol", "hR_pol", "tenure_probs", "fert_probs",
               "bp_pol_stay", "c_pol_stay", "g", "g_stay_distribution",
               "g_beginning_distribution")
    result: dict[str, object] = {}
    with np.load(path, allow_pickle=False) as archive:
        for name in checked:
            if name not in archive.files or getattr(saved, name, None) is None:
                result[name] = {"available": False}
                continue
            candidate = archive[name]
            reference = np.asarray(getattr(saved, name))
            exact = bool(candidate.shape == reference.shape and np.array_equal(candidate, reference, equal_nan=True))
            result[name] = {"available": True, "exact": exact,
                            "shape": list(candidate.shape)}
    return result


def _run_followup(saved: SimpleNamespace, P0: SimpleNamespace, old_grid: np.ndarray,
                  solver: SimpleNamespace, solve_at_price) -> int:
    followup_out = OUT / "followup_run2"
    prior_attempt_out = OUT / "followup"
    previous_summary_path = OUT / "summary.json"
    previous_summary = json.loads(previous_summary_path.read_text())
    old_summary = previous_summary.get("cases", {}).get("occupied_region_refinement", {})
    baseline_summary = previous_summary.get("cases", {}).get("baseline120", {})
    cached_path = OUT / "occupied_region_refinement" / "solution_arrays.npz"
    cached_baseline_path = OUT / "baseline120" / "solution_arrays.npz"
    if previous_summary.get("status") != "complete" or old_summary.get("grid_nodes") != 214:
        raise RuntimeError("the completed v1 214-node occupied-region case is unavailable or has changed")
    if not cached_path.exists():
        raise FileNotFoundError(f"cached v1 214-node arrays are missing: {cached_path}")
    if not cached_baseline_path.exists():
        raise FileNotFoundError(f"cached v1 baseline arrays are missing: {cached_baseline_path}")
    if baseline_summary.get("all_available_saved_arrays_exact") is not True:
        raise RuntimeError("the v1 baseline does not carry a prior exact 11-array replay receipt")
    prior_receipt = baseline_summary.get("saved_array_replay", {})
    if len(prior_receipt) != 11 or not all(item.get("available") and item.get("exact") for item in prior_receipt.values()):
        raise RuntimeError("the v1 baseline exact replay receipt is incomplete or nonexact")
    canonical_digest = hashlib.sha256(SAVED_CASE.read_bytes()).hexdigest()
    if canonical_digest != previous_summary.get("saved_arrays_sha256"):
        raise RuntimeError("canonical saved soft baseline hash differs from the v1 replay receipt")
    replay = _exact_archive_replay(cached_baseline_path, saved)
    if len(replay) != 11 or not all(item.get("available") and item.get("exact") for item in replay.values()):
        raise RuntimeError(f"cached v1 baseline120 arrays no longer exactly match canonical saved arrays: {replay}")
    if float(P0.reference_price) != float(previous_summary.get("price")):
        raise RuntimeError("current selected reference price differs from the completed v1 diagnosis")
    if not np.allclose(np.asarray(P0.H0, dtype=float).reshape(-1),
                       np.asarray(previous_summary.get("H0"), dtype=float).reshape(-1), rtol=0, atol=1e-12):
        raise RuntimeError("current derived H0 differs from the completed v1 diagnosis")
    cached_hash = hashlib.sha256(cached_path.read_bytes()).hexdigest()
    cached214 = _load_saved_npz(cached_path, float(P0.reference_price))
    variants, exact_core = _followup_grid_variants(old_grid)
    expected_old214 = _grid_variants(old_grid)["occupied_region_refinement"]
    if not np.array_equal(np.asarray(cached214.b_grid), expected_old214):
        raise RuntimeError("cached v1 214-node wealth grid differs from the recorded midpoint refinement")
    if not np.array_equal(np.asarray(cached214.g_cross_sectional_wealth_distribution),
                          np.asarray(cached214.g_beginning_distribution)):
        raise RuntimeError("cached 214-case wealth-distribution alias is not its beginning distribution")
    if not np.array_equal(np.asarray(cached214.b_grid)[_grid_node_map(old_grid, cached214.b_grid)], old_grid):
        raise RuntimeError("cached 214 grid lost an original wealth node")
    expected_input_fingerprint = _followup_input_fingerprint(P0)

    started = time.monotonic()
    memory_record = _memory_limit()
    followup_out.mkdir(parents=True, exist_ok=True)
    metadata: dict[str, object] = {
        "status": "running",
        "classification": "follow-up fixed-price grid diagnosis; no recalibration or general-equilibrium root",
        "reference_case": "authenticated selected soft checkpoint chain16/case0046",
        "selection_source_key": "chain_16/search/best_so_far.json",
        "price": float(P0.reference_price),
        "H0": np.asarray(P0.H0, dtype=float).reshape(-1).tolist(),
        "input_fingerprint_excluding_grid_and_exact_evidence_scratch_path": expected_input_fingerprint,
        "excluded_process_local_field": "native_inherited_distribution_evidence_dir",
        "baseline_authentication": {
            "canonical_saved_baseline_sha256": canonical_digest,
            "prior_v1_11_array_native_exact_replay_receipt": prior_receipt,
            "cached_v1_baseline120_against_canonical_arrays": replay,
            "current_reference_price_matches_v1": True,
            "current_H0_matches_v1_with_absolute_tolerance_1e-12": True,
            "current_loader_reverifies_selected_sources_and_target_contract": True,
        },
        "v1_214_saved_arrays": str(cached_path.relative_to(ROOT)),
        "v1_214_saved_arrays_sha256": cached_hash,
        "v1_214_wealth_distribution_alias": {
            "stored_in_minimal_npz": False,
            "reconstructed_from": "g_beginning_distribution.copy()",
            "pinned_native_source": "output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source/small_credit_lab/engine/distribution.py:1294-1295",
            "basis": "native solver assigns g_cross_sectional_wealth_distribution = g.copy() immediately after capturing g_beginning_distribution = g.copy()",
        },
        "preserved_prior_failed_attempt_log": str(prior_attempt_out.relative_to(ROOT)),
        "v1_boundary_note": "The preserved v1 tail refinement used old[1:] > 33.66055 and therefore included the interval ending at 33.6605536290589. This follow-up strict tail grid starts at the original node old[searchsorted(old, 33.66055)] = exact core cutoff.",
        "exact_core_cutoff": exact_core,
        "occupied_interval_rule": "original intervals overlapping [-6.4, exact core cutoff], with quarter, midpoint, and three-quarter nodes",
        "resource_limit": memory_record,
        "wall_budget_seconds": TIME_BUDGET_SECONDS,
        "case_budget_seconds": CASE_BUDGET_SECONDS,
        "cached_214_distribution_alias_verified": True,
        "cases": {},
    }
    _json_dump(followup_out / "summary.json", metadata)
    cases: dict[str, SimpleNamespace] = {}
    for case_name, grid in variants.items():
        if time.monotonic() - started >= TIME_BUDGET_SECONDS:
            raise TimeoutError("10-minute follow-up budget reached before next case")
        P = copy.deepcopy(P0)
        _install_grid(P, old_grid, grid)
        if _followup_input_fingerprint(P) != expected_input_fingerprint:
            raise RuntimeError(f"non-grid input changed in follow-up case {case_name}")
        case_start = time.monotonic()
        remaining = max(0.1, TIME_BUDGET_SECONDS - (case_start - started))
        _set_deadline(min(CASE_BUDGET_SECONDS, remaining))
        try:
            sol = solve_at_price(P, grid, solver, float(P0.reference_price))
        finally:
            _clear_deadline()
        elapsed = time.monotonic() - case_start
        if elapsed > CASE_BUDGET_SECONDS:
            raise TimeoutError(f"{case_name} solve exceeded its 180-second limit")
        if float(getattr(sol, "entry_censored_mass", 0.0)) > 0.0 or float(getattr(sol, "entry_censored_share", 0.0)) > 0.0:
            raise RuntimeError(f"{case_name} has entrant censorship")
        cases[case_name] = sol
        case_dir = followup_out / case_name
        remaining = max(0.1, TIME_BUDGET_SECONDS - (time.monotonic() - started))
        _set_deadline(remaining)
        try:
            _save_solution(sol, case_dir / "solution_arrays.npz")
            from small_credit_lab.engine import diagnostics
            diagnostics.write_diagnostics(sol, P, case_dir / "diagnostics")
        finally:
            _clear_deadline()

        baselines = {"canonical_baseline120": (saved, old_grid),
                     "v1_occupied_refinement214": (cached214, np.asarray(cached214.b_grid, dtype=float))}
        comparisons = {}
        for baseline_name, (baseline, baseline_grid) in baselines.items():
            if case_name == "strict_upper_tail_refinement" and baseline_name != "canonical_baseline120":
                continue
            comparisons[baseline_name] = {
                "expected_policy_differences": _policy_differences(baseline, sol, baseline_grid, grid, P0),
                "selected_renter_slice": _renter_slice_difference(baseline, sol, baseline_grid, grid),
                "pooled_wealth_distribution": _distribution_comparison(baseline, sol, baseline_grid, grid),
                "whole_population_aggregates": _aggregate_comparison(baseline, sol, P0),
            }
        record = {
            "status": "complete", "elapsed_solve_seconds": elapsed,
            "grid_nodes": int(grid.size), "grid_min": float(grid[0]), "grid_max": float(grid[-1]),
            "total_population_mass": float(getattr(sol, "total_mass", np.nan)),
            "g_beginning_mass": float(np.asarray(sol.g_beginning_distribution).sum()),
            "upper_endpoint_mass": float(np.asarray(sol.g_beginning_distribution)[-1].sum()),
            "entry_censored_mass": float(getattr(sol, "entry_censored_mass", 0.0)),
            "comparisons": comparisons,
        }
        metadata["cases"][case_name] = record
        metadata["elapsed_seconds"] = time.monotonic() - started
        _json_dump(followup_out / "summary.json", metadata)
        _json_dump(followup_out / f"{case_name}.json", record)
    metadata["status"] = "complete"
    metadata["elapsed_seconds"] = time.monotonic() - started
    _json_dump(followup_out / "summary.json", metadata)
    return 0


def _supplemental_plot(base: SimpleNamespace, cases: dict[str, SimpleNamespace],
                       old_grid: np.ndarray, outpath: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    j, z, loc, ten, n, cs = 3, 4, 0, 0, 0, 0
    left_x = (-8.0, 80.0)
    base_y = np.asarray(base.c_pol[:, ten, loc, j, z, n, cs], dtype=float)
    visible_y = [base_y[(old_grid >= left_x[0]) & (old_grid <= left_x[1])]]
    axes[0].plot(old_grid, base_y, "o-", ms=2.5, lw=1, label="baseline120")
    for name, sol in cases.items():
        if name == "baseline120":
            continue
        grid = np.asarray(sol.b_grid)
        series = np.asarray(sol.c_pol[:, ten, loc, j, z, n, cs], dtype=float)
        visible_y.append(series[(grid >= left_x[0]) & (grid <= left_x[1])])
        axes[0].plot(grid, series, ".-", ms=2, lw=0.9, label=name)
    visible = np.concatenate(visible_y)
    visible = visible[np.isfinite(visible)]
    if visible.size:
        low, high = float(np.min(visible)), float(np.max(visible))
        pad = max(0.06 * (high - low), 0.02 * max(abs(low), abs(high), 1.0))
        axes[0].set_ylim(low - pad, high + pad)
    axes[0].set_xlim(*left_x)
    axes[0].set_xlabel("Financial assets b (model units)")
    axes[0].set_ylabel("Consumption per model period")
    axes[0].set_title("Age 30, income state 5, childless renter")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(old_grid, base.c_pol[:, ten, loc, j, z, n, cs], "o-", ms=2, lw=1, label="baseline120")
    for name, sol in cases.items():
        if name == "baseline120":
            continue
        grid = np.asarray(sol.b_grid)
        axes[1].plot(grid, sol.c_pol[:, ten, loc, j, z, n, cs], ".-", ms=2, lw=0.9, label=name)
    axes[1].set_xlim(0, 6000)
    axes[1].set_xlabel("Financial assets b (model units)")
    axes[1].set_ylabel("Consumption per model period")
    axes[1].set_title("Upper tail")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def _set_deadline(seconds: float) -> None:
    if hasattr(signal, "setitimer"):
        def expired(_signum: int, _frame: object) -> None:
            raise TimeoutError("native work exceeded its wall-clock limit")
        signal.signal(signal.SIGALRM, expired)
        signal.setitimer(signal.ITIMER_REAL, max(0.1, seconds))


def _clear_deadline() -> None:
    if hasattr(signal, "setitimer"):
        signal.setitimer(signal.ITIMER_REAL, 0.0)


def main() -> int:
    args = sys.argv[1:]
    if any(arg not in {"--followup", "-h", "--help"} for arg in args):
        raise SystemExit("usage: diagnose_playground_asset_grid.py [--followup]")
    if "-h" in args or "--help" in args:
        print(__doc__)
        print("Default runs v1 four-case grid diagnosis; --followup runs the two-case refinement against saved baselines.")
        return 0
    started = time.monotonic()
    resource_record = _memory_limit()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "case_summaries").mkdir(exist_ok=True)

    # Import the hash-checked playground only after thread controls and limits.
    sys.path.insert(0, str(ROOT / "code/model/tools"))
    from model_playground import load_reference_model, load_saved_solution, solve_at_price

    saved = load_saved_solution("soft")
    P0, old_grid, solver = load_reference_model()
    old_grid = np.asarray(old_grid, dtype=float)
    if not np.array_equal(old_grid, np.asarray(saved.b_grid, dtype=float)):
        raise RuntimeError("native builder grid differs from the saved selected soft reference")
    if "--followup" in args:
        return _run_followup(saved, P0, old_grid, solver, solve_at_price)
    grids = _grid_variants(old_grid)
    entry_before = np.asarray(P0.fixed_reference_entry_conditional, dtype=float).copy()
    economic_fingerprint = hashlib.sha256(json.dumps({
        key: repr(value) for key, value in sorted(vars(P0).items())
        if key not in {"Nb", "b_min", "b_max", "earnings_transaction_grid", "fixed_reference_entry_grid", "fixed_reference_entry_conditional"}
        and not key.startswith("_")
    }, sort_keys=True, default=str).encode()).hexdigest()
    meta: dict[str, object] = {
        "status": "running", "classification": "fixed-price numerical grid comparison; no recalibration or general-equilibrium root",
        "price": float(P0.reference_price), "H0": np.asarray(P0.H0, dtype=float).reshape(-1).tolist(),
        "reference_case": "soft financing, original transaction-inside timing, authenticated selected soft checkpoint chain16/case0046",
        "source": "model_playground.load_reference_model; authenticated frozen native engine",
        "saved_arrays_sha256": hashlib.sha256(SAVED_CASE.read_bytes()).hexdigest(),
        "economic_inputs_fingerprint_excluding_grid": economic_fingerprint,
        "resource_limit": resource_record, "wall_budget_seconds": TIME_BUDGET_SECONDS,
        "case_budget_seconds": CASE_BUDGET_SECONDS, "cases": {},
        "grid_definitions": {k: {"nodes": int(v.size), "minimum": float(v[0]), "maximum": float(v[-1])} for k, v in grids.items()},
        "entry_grid_contract": {"original_nodes": int(old_grid.size), "entry_conditional_shape": list(entry_before.shape),
                                "entry_row_mass_preserved": True, "first_and_second_moments_preserved": True,
                                "no_relocation_or_censoring": True},
    }
    _json_dump(OUT / "summary.json", meta)
    solutions: dict[str, SimpleNamespace] = {}
    try:
        for case_name, grid in grids.items():
            elapsed_total = time.monotonic() - started
            if elapsed_total >= TIME_BUDGET_SECONDS:
                raise TimeoutError("10-minute total wall budget reached before next case")
            P = copy.deepcopy(P0)
            _install_grid(P, old_grid, grid)
            if economic_fingerprint != hashlib.sha256(json.dumps({
                key: repr(value) for key, value in sorted(vars(P).items())
                if key not in {"Nb", "b_min", "b_max", "earnings_transaction_grid", "fixed_reference_entry_grid", "fixed_reference_entry_conditional"}
                and not key.startswith("_")
            }, sort_keys=True, default=str).encode()).hexdigest():
                raise RuntimeError(f"non-grid input changed in {case_name}")
            case_start = time.monotonic()
            remaining_total = max(0.1, TIME_BUDGET_SECONDS - (case_start - started))
            _set_deadline(min(CASE_BUDGET_SECONDS, remaining_total))
            try:
                sol = solve_at_price(P, grid, solver, float(P0.reference_price))
            finally:
                _clear_deadline()
            elapsed = time.monotonic() - case_start
            if elapsed > CASE_BUDGET_SECONDS:
                raise TimeoutError(f"{case_name} solve took {elapsed:.1f}s; per-case limit is 180s")
            solutions[case_name] = sol
            case_dir = OUT / case_name
            remaining_total = max(0.1, TIME_BUDGET_SECONDS - (time.monotonic() - started))
            _set_deadline(remaining_total)
            try:
                _save_solution(sol, case_dir / "solution_arrays.npz")
                from small_credit_lab.engine import diagnostics
                diagnostics.write_diagnostics(sol, P, case_dir / "diagnostics")
            finally:
                _clear_deadline()

            record: dict[str, object] = {
                "status": "complete", "elapsed_seconds": elapsed,
                "grid_nodes": int(grid.size), "grid_min": float(grid[0]), "grid_max": float(grid[-1]),
                "price": float(P0.reference_price), "H0": np.asarray(P0.H0, dtype=float).reshape(-1).tolist(),
                "total_population_mass": float(getattr(sol, "total_mass", np.nan)),
                "g_beginning_mass": float(np.asarray(sol.g_beginning_distribution).sum()),
                "upper_endpoint_distribution_mass": float(np.asarray(sol.g_beginning_distribution)[-1].sum()),
                "entry_censored_mass": float(getattr(sol, "entry_censored_mass", 0.0)),
                "entry_censored_share": float(getattr(sol, "entry_censored_share", 0.0)),
            }
            if record["entry_censored_mass"] > 0.0 or record["entry_censored_share"] > 0.0:
                raise RuntimeError(f"{case_name} has nonzero native entrant censorship")
            if case_name == "baseline120":
                record["population_aggregates"] = _aggregate_comparison(sol, sol, P)
                record["saved_array_replay"] = _exact_baseline_check(sol, saved)
                required_exact = all(
                    item.get("exact", False)
                    for item in record["saved_array_replay"].values()
                    if item.get("available")
                )
                record["all_available_saved_arrays_exact"] = bool(required_exact)
                meta["cases"][case_name] = record
                _json_dump(OUT / "case_summaries" / f"{case_name}.json", record)
                if not required_exact:
                    raise RuntimeError("native baseline does not exactly replay every available saved reference array; candidate comparisons stopped")
            else:
                record["inherited_state_policy_differences"] = _policy_differences(
                    solutions["baseline120"], sol, old_grid, grid, P0
                )
                record["selected_renter_slice"] = _renter_slice_difference(
                    solutions["baseline120"], sol, old_grid, grid
                )
                record["cross_sectional_wealth_distribution"] = _distribution_comparison(
                    solutions["baseline120"], sol, old_grid, grid
                )
                record["population_aggregates"] = _aggregate_comparison(
                    solutions["baseline120"], sol, P0
                )
            meta["cases"][case_name] = record
            meta["elapsed_seconds"] = time.monotonic() - started
            _json_dump(OUT / "summary.json", meta)
            _json_dump(OUT / "case_summaries" / f"{case_name}.json", record)
            if time.monotonic() - started >= TIME_BUDGET_SECONDS:
                raise TimeoutError("10-minute total wall budget reached after completed case")
        _supplemental_plot(solutions["baseline120"], solutions, old_grid, OUT / "supplemental_grid_policy_comparison.png")
        meta["status"] = "complete"
    except Exception as exc:
        meta["status"] = "failed"
        meta["failure"] = f"{type(exc).__name__}: {exc}"
        meta["elapsed_seconds"] = time.monotonic() - started
        _json_dump(OUT / "summary.json", meta)
        raise
    meta["elapsed_seconds"] = time.monotonic() - started
    _json_dump(OUT / "summary.json", meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
