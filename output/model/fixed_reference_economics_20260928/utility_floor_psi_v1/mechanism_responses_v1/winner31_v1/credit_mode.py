"""Isolated mode and support gates for the current floor candidate experiments.

This module does not edit the engine. The `lifetime_repayment_only` arm uses
the existing native solvency mode with the DUE stayer rule disabled.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
import numpy as np


PRICE_FACTORS = (0.99, 1.0, 1.01)
BOUNDARY_ATOL = 1e-10
OCCUPIED_MASS_TOL = 0.0


def configure_credit_arm(P, mode: str, validate_native_mode, bind_engine_credit):
    """Return a copy with only the prescribed credit-rule switch applied."""
    Q = copy.deepcopy(P)
    if mode == "reference":
        # Bind the actual solved-candidate mode before checking its credit fields.
        bind_engine_credit(Q, "corrected", 0.0)
        if bool(getattr(Q, "native_solvency_credit", False)) or not bool(getattr(Q, "native_due_stayer_credit", False)):
            raise ValueError("Reference arm differs from the candidate DUE contract")
        if float(Q.unsecured_credit_limit) != 0.0:
            raise ValueError("Reference arm failed to bind the candidate's zero unsecured-credit limit")
        return Q
    if mode != "lifetime_repayment_only":
        raise ValueError(f"Unknown credit regime: {mode}")
    Q.native_due_stayer_credit = False
    Q.native_solvency_credit = True
    Q.unsecured_credit_limit = None
    if bool(getattr(Q, "native_due_stayer_credit", True)):
        raise ValueError("DUE stayer rule remains enabled in the natural-credit arm")
    if getattr(Q, "unsecured_credit_limit", 0.0) is not None:
        raise ValueError("A fixed unsecured-credit limit remains in the natural-credit arm")
    if not validate_native_mode(Q):
        raise ValueError("Native solvency mode did not validate")
    return Q


@contextmanager
def trace_native_support(household_module):
    """Observe the solver's native support-floor calls without changing results."""
    original_savings = household_module._savings_stage
    original_support = household_module.native_solvency_support_floor
    current = {}
    calls = []

    def savings(*args, **kwargs):
        if len(args) < 10:
            raise RuntimeError("Native saving-stage signature changed")
        P, b_grid, _, _, j, z_value = args[1], args[2], args[3], args[4], args[6], args[7]
        current.update(P=P, b_grid=np.asarray(b_grid), age=int(j), z_value=float(z_value), tenure_call=0)
        try:
            result = original_savings(*args, **kwargs)
        finally:
            expected = 1 + int(P.n_house)
            observed = current.get("tenure_call", 0)
            current.clear()
        if observed != expected:
            raise RuntimeError(f"Expected {expected} native support calls, saw {observed}")
        return result

    def support(values, grid):
        floor, dead = original_support(values, grid)
        if not current:
            raise RuntimeError("Native support floor called outside observed saving stage")
        k = int(current["tenure_call"])
        if k > int(current["P"].n_house):
            raise RuntimeError("Unexpected native support-floor branch")
        z = np.asarray(current["P"].z_grid, dtype=float)
        zz = int(np.argmin(np.abs(z - current["z_value"])))
        if abs(float(z[zz]) - current["z_value"]) > 1e-12:
            raise RuntimeError("Saving-stage income is absent from the candidate Markov grid")
        calls.append({"age": current["age"], "income_index": zz,
                      "tenure_index": k, "floor": np.asarray(floor, dtype=float).copy()})
        current["tenure_call"] = k + 1
        return floor, dead

    household_module._savings_stage = savings
    household_module.native_solvency_support_floor = support
    try:
        yield calls
    finally:
        household_module._savings_stage = original_savings
        household_module.native_solvency_support_floor = original_support


def lower_grid_diagnostic(solution, b_grid, support_calls=None, P=None, realized_distribution=None):
    """Measure lower-endpoint choices on positive-mass beginning states.

    The policy occupancy weights are the calendar evaluation's realized `g_current`; the
    policies are its `bp_pol`. Any positive occupied mass selecting b_min fails.
    Missing/misaligned arrays are unresolved and fail closed.
    """
    bg = np.asarray(b_grid, dtype=float)
    g = getattr(solution, "g_beginning_distribution", None)
    bp = getattr(solution, "bp_pol", None)
    if g is None or bp is None:
        return {"status": "unresolved", "reason": "missing beginning distribution or saving policy"}
    g = np.asarray(g, dtype=float)
    bp = np.asarray(bp, dtype=float)
    if g.shape != bp.shape or g.shape[0] != bg.size or not np.isfinite(g).all() or not np.isfinite(bp).all():
        return {"status": "unresolved", "reason": "policy/distribution axes do not match the wealth grid"}
    if np.any(g < 0.0) or float(g.sum()) <= 0.0:
        return {"status": "unresolved", "reason": "invalid occupied-state weights"}
    if realized_distribution is None:
        return {"status":"unresolved","reason":"realized post-tenure saving-policy occupancy was not supplied"}
    realized=np.asarray(realized_distribution,dtype=float)
    if realized.shape!=bp.shape or not np.isfinite(realized).all() or (realized<0).any() or realized.sum()<=0:
        return {"status":"unresolved","reason":"realized post-tenure saving-policy occupancy is invalid"}
    occupied = realized > 0.0
    tol = max(BOUNDARY_ATOL, 1e-12 * float(bg[-1] - bg[0]))
    at_floor = np.abs(bp - float(bg[0])) <= tol
    mass = float(realized[occupied & at_floor].sum())
    share = mass / float(realized.sum())
    support_floor_mass = 0.0
    support_floor_hits = 0
    if support_calls is None or not support_calls:
        return {"status": "unresolved", "reason": "native support-floor calls were not observed",
                "occupied_policy_boundary_share": share}
    if P is None:
        return {"status": "unresolved", "reason": "model parameters needed for complete support-call coverage were not supplied",
                "occupied_policy_boundary_share": share}
    expected = {(j, zz, k) for j in range(int(P.J)) for zz in range(int(P.Nz))
                for k in range(1 + int(P.n_house))}
    observed = {(int(c["age"]), int(c["income_index"]), int(c["tenure_index"])) for c in support_calls}
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        return {"status": "unresolved", "reason": "native support trace did not cover every age/income/destination-tenure branch",
                "occupied_policy_boundary_share": share,
                "support_trace_coverage": {"expected_branches": len(expected), "observed_calls": len(support_calls),
                                           "unique_branches": len(observed), "duplicate_calls": len(support_calls)-len(observed), "missing": [list(x) for x in missing],
                                           "extra": [list(x) for x in extra]}}
    for call in support_calls:
        j, zz, tenure = call["age"], call["income_index"], call["tenure_index"]
        floor = np.asarray(call["floor"], dtype=float)
        # tenure in this trace identifies the destination renter/owner branch.
        # Beginning-distribution tenure is inherited and must be integrated out.
        family_mass = g[:, :, 0, j, zz, :, :].sum(axis=(0, 1)).reshape(-1, order="F")
        if floor.shape != family_mass.shape:
            return {"status": "unresolved", "reason": "support floor and occupied family-state axes differ"}
        hit = floor <= float(bg[0]) + tol
        support_floor_hits += int(np.count_nonzero(hit & (family_mass > 0.0)))
        support_floor_mass = max(support_floor_mass, float(family_mass[hit].sum()) / float(g.sum()))
    status = ("occupied_support_pass_unoccupied_alternatives_unverified"
              if share <= OCCUPIED_MASS_TOL and support_floor_mass <= OCCUPIED_MASS_TOL
              else "support_uncertified")
    return {
        "status": status,
        "guard": "positive-mass beginning states with bp_pol within tol of b_grid[0]",
        "distribution_source": "evaluation.g_current (realized post-tenure occupancy)",
        "support_family_mass_source":"solution.g_beginning_distribution (post-fertility inherited tenure aggregated)",
        "policy_source": "solution.bp_pol",
        "threshold": {"absolute_atol": BOUNDARY_ATOL, "scaled_atol": 1e-12,
                      "applied_atol": tol, "occupied_mass_share_fail_if_gt": OCCUPIED_MASS_TOL},
        "b_min": float(bg[0]),
        "occupied_boundary_mass": mass,
        "occupied_boundary_share": share,
        "occupied_state_count": int(occupied.sum()),
        "support_floor_source": "native_solvency_support_floor(Vc, b_grid) inside engine.household._savings_stage",
        "support_floor_relevance": "for each destination branch, tested against positive-mass family states after aggregating inherited tenure; family states use Fortran order",
        "support_trace_coverage": {"expected_branches": len(expected), "observed_calls": len(support_calls),
                                   "unique_branches": len(observed), "duplicate_calls": len(support_calls)-len(observed), "complete": True},
        "support_trace_tenure_meaning": "destination renter/owner branch, not inherited beginning-distribution tenure",
        "support_floor_boundary_cells": support_floor_hits,
        "occupied_support_floor_boundary_share_max": support_floor_mass,
        "limitations": ["This guard covers occupied family states and cannot rule out effects from unoccupied continuation alternatives.",
                        "The native solvency cutoff remains an approximation; a passing occupied-state guard does not certify the lifetime-repayment-only regime over the full unbounded support."],
    }
