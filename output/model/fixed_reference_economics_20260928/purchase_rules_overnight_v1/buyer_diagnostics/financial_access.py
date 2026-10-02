"""Pointwise financial access for matched first-birth renter states; no solve.

Call from an authenticated native postcheck with its actual P, SD, b_grid,
prices, pre-fertility mass and first-birth attempt probabilities. The caller
must load the matching isolated engine package as ``refactor_lab``.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np


def matched_first_birth_access(P, SD, b_grid, prices, g_pre, fert_probs, *, rule,
                               observed_tenure_probs=None, observed_phi=None):
    """Compare 80% and 100% owner *financial* access on identical birth states.

    Returns aggregate flow shares and arrays of candidate-owner feasibility.
    The function deliberately does not infer whether a household wants to buy.
    It also does not evaluate continuation values or numerical interpolation of
    the owner value function, which can rule out financially feasible options.
    """
    from refactor_lab.engine.household import (
        _build_housing_stage_ctx, child_earnings_multiplier,
        child_earnings_penalty_active, children_at_home_count,
        native_due_death_floor,
    )
    from refactor_lab.engine.parameters import get_fecundity_by_age
    from refactor_lab.engine.shared import income_at_state

    if rule not in {"hard", "quarter"}:
        raise ValueError("rule must be hard or quarter")
    if not bool(getattr(P, "native_purchase_income", False)):
        raise ValueError("Expected native purchase-income accounting")
    if bool(getattr(P, "use_pti_constraint", False)):
        raise ValueError("PTI rule needs its separate exact screen")
    if np.any(SD.birth_dp) or np.any(SD.birth_entry_grant):
        raise ValueError("Buyer waivers or grants require their own transaction accounting")
    if rule == "quarter" and abs(float(getattr(P, "experimental_purchase_saving_fraction", 0.25)) - 0.25) > 1e-12:
        raise ValueError("Expected quarter-saving fraction 0.25")
    if rule == "hard" and abs(float(getattr(P, "experimental_purchase_saving_fraction", 1.0)) - 1.0) > 1e-12:
        raise ValueError("Expected no quarter-saving adjustment in hard rule")
    if g_pre.shape != (len(b_grid), 1 + P.n_house, P.I, P.J,
                       len(np.asarray(P.z_grid)), P.n_parity, P.n_child_states):
        # If the native income-state field has a different name, the caller
        # should pass an exact dimension check rather than silently reshape.
        raise ValueError("Pre-fertility distribution dimension mismatch")
    if fert_probs.shape != g_pre.shape[:-2] + (P.n_parity,):
        raise ValueError("Fertility probability dimension mismatch")
    fec = get_fecundity_by_age(P)
    if float(np.sum(g_pre[..., 0, 1:])) > 1e-10:
        raise ValueError("Occupied childless readiness states require explicit branch weights")
    if P.n_parity < 2 or P.n_child_states < 2:
        raise ValueError("First-birth destination state is missing")
    bg = np.asarray(b_grid, float)
    q = np.asarray(prices, float)
    if q.shape != (P.I,) or not np.isfinite(q).all():
        raise ValueError("Price shape or values invalid")
    weights = g_pre[:, 0, :, :, :, 0, 0] * fert_probs[:, 0, :, :, :, 1] * fec[None, None, :, None]
    all_birth_weights = (g_pre[..., 0, 0] * fert_probs[..., 1]
                         * fec[None, None, None, :, None])
    all_first_birth_mass = float(np.sum(all_birth_weights))
    if np.any(weights < -1e-12) or not np.isfinite(weights).all():
        raise ValueError("Invalid first-birth flow")
    masks = {}
    margin_masks = {}
    details = {}
    for phi in (0.8, 1.0):
        trial_SD = copy.copy(SD)
        choice = np.ones_like(SD.phi_choice)
        choice[:, 1:, :, :] = phi
        trial_SD.phi_choice = choice
        ctx = _build_housing_stage_ctx(P, bg, trial_SD, q, True, False)
        possible = np.zeros(weights.shape + (P.n_house,), dtype=bool)
        cmin_possible = np.zeros_like(possible)
        reasons = {"screen": 0.0, "transaction_support": 0.0,
                   "housing_floor": 0.0, "saving_budget": 0.0}
        family_col = 1 + P.n_parity  # (children-ever-born=1, at-home=1), Fortran order
        cb = float(SD.cb_flat[0, family_col])
        hb = float(SD.hb_flat[0, family_col])
        grant_floor = float(SD.gb_flat[0, family_col])
        for i in range(P.I):
            for j in range(P.J):
                for zz, z in enumerate(np.asarray(P.z_grid)):
                    y = float(income_at_state(P, i, j, float(z)))
                    yadj = 0.0
                    if child_earnings_penalty_active(P) and j < int(getattr(P, "J_R", P.J)):
                        m = children_at_home_count(1, 1, P)
                        yadj = float(P.income[i, j]) * float(z) * (child_earnings_multiplier(P, j, m) - 1.0)
                    income = y + yadj
                    for dest in range(1, P.n_house + 1):
                        Q = float(ctx.hcost[i, dest])
                        x = bg - Q  # renter origin: A=b, no sale
                        bf = -phi * Q
                        dp = (1 - phi) * Q
                        if rule == "hard":
                            screen = (bg >= dp) & (x >= max(bf, bg[0]))
                            floor = np.full_like(bg, bf)
                        else:
                            # The original soft purchase-income screen remains
                            # in the quarter variant's tenure-choice kernel.
                            screen = (bg >= dp - y / P.R_gross) & (
                                x >= max(bf - y / P.R_gross, bg[0]))
                            floor = bf + 3.0 * np.maximum(0.0, bf - x)
                        support = (x >= bg[0]) & (x <= bg[-1])
                        housing = (not ctx.strict_owner_hbar_feasibility or
                                   ctx.hsrv[i, dest] > ctx.owner_h_bar_scale * hb)
                        death_floor = native_due_death_floor(P, j, q[i], P.H_own[dest - 1])
                        floor = np.maximum.reduce((floor, np.full_like(bg, bf),
                                                   np.full_like(bg, death_floor),
                                                   np.full_like(bg, bg[0])))
                        Rv = P.R_gross * x + income
                        Rvt = P.R_gross * np.maximum(x, 0.0) + income
                        transfer = (np.clip(grant_floor - Rvt, 0.0, grant_floor)
                                    if grant_floor > 0 else np.zeros_like(Rvt))
                        resources = Rv + transfer - ctx.ocst[i, dest]
                        # The optimizer uses hi=resources-cb-1e-6, but when
                        # hi<lo it evaluates lo itself. The owner's utility
                        # kernel rejects residual consumption <=1e-10.
                        budget = (resources - cb - floor) > 1e-10
                        idx = (slice(None), i, j, zz, dest - 1)
                        possible[idx] = screen & support & housing & budget
                        cmin_possible[idx] = screen & support & housing & (
                            resources - cb - P.c_min > floor)
                        active = weights[:, i, j, zz]
                        reasons["screen"] += float(np.sum(active * ~screen))
                        reasons["transaction_support"] += float(np.sum(active * screen * ~support))
                        if not housing:
                            reasons["housing_floor"] += float(np.sum(active * screen * support))
                        reasons["saving_budget"] += float(np.sum(active * screen * support * housing * ~budget))
        masks[phi] = possible
        margin_masks[phi] = cmin_possible
        any_possible = possible.any(axis=-1)
        any_cmin = cmin_possible.any(axis=-1)
        total = float(weights.sum())
        details[str(phi)] = {
            "eligible_share": float(np.sum(weights * any_possible) / total) if total else None,
            "origin_renter_financially_ineligible_share_among_all_first_births": (
                float(np.sum(weights * ~any_possible) / all_first_birth_mass)
                if all_first_birth_mass else None),
            "eligible_share_with_cmin_margin": float(np.sum(weights * any_cmin) / total) if total else None,
            "candidate_failure_weighted_counts_overlapping": reasons,
        }
    total = float(weights.sum())
    if np.any(masks[0.8] & ~masks[1.0]):
        raise RuntimeError("Raising financed share removed a pointwise feasible owner product")
    gain = (~masks[0.8].any(axis=-1)) & masks[1.0].any(axis=-1)
    observed_check = None
    if observed_tenure_probs is not None:
        if observed_phi not in (0.8, 1.0):
            raise ValueError("Observed tenure audit requires observed_phi 0.8 or 1.0")
        if observed_tenure_probs.shape != g_pre.shape + (1 + P.n_house,):
            raise ValueError("Observed tenure probability dimensions mismatch")
        owner_prob = observed_tenure_probs[:, 0, :, :, :, 1, 1, 1:]
        chosen = weights[..., None] * owner_prob
        chosen_mass = float(np.sum(chosen))
        unsupported_mass = float(np.sum(chosen * ~masks[observed_phi]))
        observed_check = {
            "observed_phi": observed_phi,
            "birth_branch_owner_choice_mass": chosen_mass,
            "chosen_owner_mass_outside_pointwise_financial_map": unsupported_mass,
            "share_outside_map": unsupported_mass / chosen_mass if chosen_mass else None,
        }
    return {
        "rule": rule,
        "first_birth_origin_renter_flow_mass": total,
        "all_origin_first_birth_flow_mass": all_first_birth_mass,
        "renter_origin_share_of_first_births": total / all_first_birth_mass if all_first_birth_mass else None,
        "phi": details,
        "share_ineligible_at_80_eligible_at_100": float(np.sum(weights * gain) / total) if total else None,
        "share_ineligible_at_80_eligible_at_100_among_all_first_births": (
            float(np.sum(weights * gain) / all_first_birth_mass) if all_first_birth_mass else None),
        "financial_access_only": True,
        "continuation_value_and_policy_choice_not_checked": True,
        "solver_consumption_surplus_threshold": 1e-10,
        "cmin_margin_reported_separately": True,
        "observed_choice_audit": observed_check,
        "weights": weights,
        "owner_feasible_at_80": masks[0.8],
        "owner_feasible_at_100": masks[1.0],
    }


def write_compact_summary(path, result):
    """Save the scalar receipt; arrays stay with the authenticated native stage."""
    omitted = {"weights", "owner_feasible_at_80", "owner_feasible_at_100"}
    record = {k: v for k, v in result.items() if k not in omitted}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
