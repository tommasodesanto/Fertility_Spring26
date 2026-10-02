"""Household stage of the extracted solver (bodies byte-identical; see split_receipt.json)."""
from __future__ import annotations
import copy
import math
import time
from types import SimpleNamespace
from typing import Any
import numpy as np
from . import joint_nested
from .adult_entry import adjusted_births, potential_entry_households
from .warm_price import search_warm_price
from .child_preferences import apply_child_preferences
from .parameters import (bequest_utility_net_active, child_earnings_multiplier, child_earnings_penalty_active, children_at_home_count, estate_housing_value, estate_receiver_active, estate_transfer_at_age, get_fecundity_by_age, independent_child_maturation_active, mortgage_stay_floor_active, parent_age_maturation_active, readiness_childless_states, readiness_cumulative_probability, readiness_gate_active, readiness_settled_state, readiness_transition_hazard, rental_wedge_active, unsecured_debt_floor)
from .kernels import (NUMBA_AVAILABLE, full_owner_block_kernel, full_renter_block_kernel, location_logit_kernel, scatter_cols_kernel, scatter_cols_sameidx_kernel, scatter_vec_kernel, tenure_choice_kernel, tenure_logit_kernel)
from .utils import (decode_flat_family_state, flat_nc, interp_indices, interp_on_grid, interp_vector, logsumexp, make_grid, make_value_interp, scatter_redistribute, scatter_redistribute_cols, scatter_redistribute_cols_sameidx, unflat_nc, weighted_median_from_cells, weighted_quantile)
from .shared import (DEAD_VALUE_CUTOFF, get_completed_fertility, income_at_state, income_transition_values)


def debt_rule_at_age(P: SimpleNamespace, current_unsecured: Any, j: int) -> np.ndarray:
    """Evaluate the unsecured component of the next-period debt floor."""

    return unsecured_debt_floor(
        current_unsecured,
        float(np.asarray(P.debt_taper_weights)[int(j) + 1]),
        float(np.asarray(P.debt_caps)[int(j) + 1]),
    )


def fixed_unsecured_credit_active(P: SimpleNamespace) -> bool:
    return getattr(P, "unsecured_credit_limit", None) is not None


def renter_borrowing_floor(P: SimpleNamespace, b: Any, j: int) -> np.ndarray:
    """Renter floor; scalar credit is separate from the legacy rollover rule."""

    credit = getattr(P, "unsecured_credit_limit", None)
    if credit is None:
        return debt_rule_at_age(P, b, j)
    floor = -float(credit)
    death_possible = int(j) == int(P.J) - 1 or (
        bool(getattr(P, "use_age_survival", False))
        and float(np.asarray(P.survival_probs)[int(j)]) < 1.0
    )
    # Estates remain non-negative: a death branch makes the effective floor
    # max(-D, 0), independently of the age-taper arrays.
    if death_possible:
        floor = max(floor, 0.0)
    return np.zeros_like(np.asarray(b, dtype=float)) + floor


def native_due_owner_floor(b, collateral_floor, *, death_floor=-np.inf):
    """DUE stayer rule in asset units: service interest, do not grow excess debt.

    The separate net-estate floor preserves no negative estates at possible
    death; it can require repayment following a sufficiently large price fall.
    """
    return np.maximum(np.minimum(np.asarray(b, dtype=float),
                                 np.asarray(collateral_floor, dtype=float)), death_floor)


def native_due_death_floor(P, j, price, house):
    death_possible = (j == int(P.J) - 1 or
        (bool(getattr(P, "use_age_survival", False)) and float(P.survival_probs[j]) < 1.0))
    return -(1.0 - float(P.psi)) * float(price) * float(house) if death_possible else -np.inf


def owner_borrowing_floor(
    P: SimpleNamespace,
    b: Any,
    collateral_floor: Any,
    j: int,
    *,
    stay_on: bool = False,
    stay_orig: bool = False,
    amort: float = 0.0,
) -> np.ndarray:
    """Owner floor after separating secured from unsecured debt.

    Prices and therefore ``collateral_floor`` are those of the current solver
    iterate.  This convention is inert in stationary equilibrium.
    With ``stay_on``, the stayer rule replaces the taper/line rollover: debt
    may not rise (``stay_orig``), and with ``amort`` it must fall by at least
    that share; with no debt the collateral floor applies under
    ``stay_orig``.  Defaults reproduce the legacy floor bit for bit.
    """

    b_arr = np.asarray(b, dtype=float)
    bf_arr = effective_owner_collateral_floor(P, collateral_floor, j)
    if bool(getattr(P, "native_purchase_income", False)):
        if stay_on:
            if bool(getattr(P, "native_due_stayer_credit", False)):
                return native_due_owner_floor(b_arr, bf_arr)
            raise ValueError("native purchase-income floor does not combine stayer rules")
        return np.broadcast_to(bf_arr, np.broadcast_shapes(b_arr.shape, bf_arr.shape)).copy()
    if not stay_on:
        current_unsecured = b_arr - bf_arr
        return bf_arr + debt_rule_at_age(P, current_unsecured, j)
    standard = bf_arr + debt_rule_at_age(P, b_arr - bf_arr, j)
    amort_floor = b_arr * (1.0 - float(amort))
    if stay_orig:
        return np.where(b_arr < 0.0, np.maximum(amort_floor, b_arr), bf_arr)
    return np.where(b_arr < 0.0, np.maximum(amort_floor, standard), standard)


def effective_owner_collateral_floor(
    P: SimpleNamespace,
    collateral_floor: Any,
    j: int,
) -> np.ndarray:
    """Apply the opt-in next-period post-retirement LTV schedule.

    ``collateral_floor`` uses the ordinary financed-share convention
    ``-phi * p * H``.  Bellman choices at age index ``j`` select next-period
    liquid wealth, so the schedule is evaluated at ``j + 1``.  With the
    experiment switched off, the multiplier is one and this helper is inert.
    """

    base = np.asarray(collateral_floor, dtype=float)
    multipliers = np.asarray(getattr(P, "owner_ltv_multipliers", np.ones(int(P.J) + 1)), dtype=float)
    next_idx = min(max(int(j) + 1, 0), multipliers.size - 1)
    return float(multipliers[next_idx]) * base


def birth_destination_child_state(P: SimpleNamespace, current_child_state: int) -> int:
    """Child-state destination after a successful upward birth attempt.

    The historical shared-clock architecture resets every birth to its single
    active child state.  Only the independent-count repair increments the
    number of children currently at home.
    """
    if independent_child_maturation_active(P):
        return int(current_child_state) + 1
    return 1


def _build_housing_stage_ctx(
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    p_hat: np.ndarray,
    use_full_kernel: bool,
    exhaustive_saving: bool,
) -> SimpleNamespace:
    """Loop-invariant arrays for the housing/saving + tenure/location stages.

    Packs the objects the per-age housing block needs so the same code can run
    once under the standard continuation ``Vc`` and — at fertile ages in
    parent-age mode — a second time under the newborn-exempt continuation
    ``Vc_ex``.  Bodies below are verbatim moves out of
    ``solve_bellman_full_markov_income`` so the standard path stays bitwise
    identical.
    """
    Nb = len(b_grid)
    nt = 1 + P.n_house
    I = P.I
    npar = P.n_parity
    ncs = P.n_child_states
    owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    owner_service_premium = max(float(getattr(P, "chi", 1.0)), 1e-8)
    strict_owner_hbar_feasibility = int(
        bool(getattr(P, "child_room_floor", False))
        and (
            float(getattr(P, "hbar_child_rooms", 0.0)) > 0.0
            or float(getattr(P, "hbar_first_child_jump", 0.0)) > 0.0
        )
    )
    owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
    owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
    owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))

    phi_choice = SD.phi_choice
    hcost = np.zeros((I, nt))
    heq = np.zeros((I, nt))
    dp_arr = np.zeros((I, nt, npar, ncs))
    bmo = np.zeros((I, nt, npar, ncs))
    hsrv = np.zeros((I, nt))
    ocst = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hcost[i, ten] = p_hat[i] * hs
            heq[i, ten] = (1 - P.psi) * p_hat[i] * hs
            hsrv[i, ten] = hs
            extra_size_cost = owner_size_cost * p_hat[i] * max(hs - owner_size_cost_ref, 0.0) ** owner_size_cost_power
            ocst[i, ten] = (P.delta + P.tau_H) * p_hat[i] * hs + extra_size_cost
            for nn in range(npar):
                for cs in range(ncs):
                    phi_ncs = phi_choice[i, ten, nn, cs]
                    dp_arr[i, ten, nn, cs] = (1 - phi_ncs) * hcost[i, ten]
                    bmo[i, ten, nn, cs] = -phi_ncs * hcost[i, ten]

    loc_shift = np.zeros((I, I))
    for io in range(I):
        for id_ in range(I):
            move_cost = P.mu_stay if id_ == io else P.mu_move
            loc_shift[io, id_] = P.E_loc[id_] - move_cost

    iidx = np.zeros((Nb, I, nt), dtype=np.int64)
    iwt = np.zeros((Nb, I, nt))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + heq[io, to], b_grid[0], b_grid[-1])
            idx, wt = interp_indices(b_grid, ba)
            iidx[:, io, to] = idx
            iwt[:, io, to] = wt

    gs_tol = 1e-3
    gs_alpha1 = (3 - math.sqrt(5)) / 2
    gs_alpha2 = (math.sqrt(5) - 1) / 2

    return SimpleNamespace(
        hcost=hcost,
        current_prices=np.asarray(p_hat, dtype=float).copy(),
        heq=heq,
        dp_arr=dp_arr,
        bmo=bmo,
        hsrv=hsrv,
        ocst=ocst,
        loc_shift=loc_shift,
        iidx=iidx,
        iwt=iwt,
        cb_v=np.ascontiguousarray(SD.cb_flat.reshape(-1)),
        hb_v=np.ascontiguousarray(SD.hb_flat.reshape(-1)),
        psi_v_flat=np.ascontiguousarray(SD.psi_flat.reshape(-1)),
        gb_v=np.ascontiguousarray(SD.gb_flat.reshape(-1)),
        alpha_v=np.ascontiguousarray(SD.alpha_flat.reshape(-1)),
        esc_v=np.ascontiguousarray(SD.escale_flat.reshape(-1)),
        b=b_grid.reshape(-1, 1),
        gs_tol=gs_tol,
        gs_alpha1=gs_alpha1,
        gs_alpha2=gs_alpha2,
        use_full_kernel=bool(use_full_kernel),
        exhaustive_saving=bool(exhaustive_saving),
        strict_owner_hbar_feasibility=strict_owner_hbar_feasibility,
        owner_h_bar_scale=owner_h_bar_scale,
        owner_service_premium=owner_service_premium,
    )


def _savings_stage(
    Vc_arr: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    ctx: SimpleNamespace,
    r_hat: np.ndarray,
    j: int,
    z_value: float,
    s_next: float,
    D_next: float,
    renter_floor: np.ndarray,
    stay_floor: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Housing/saving stage for a given continuation array.

    Runs the renter + owner savings optimizations with ``Vc_arr`` as the
    continuation (standard ``Vc`` or newborn-exempt ``Vc_ex``) and returns
    ``(Vd, cd, hd, bd)``.  Policies are re-optimized, so calling this with
    ``Vc_ex`` yields the exact exempt values, not an envelope shift.
    With ``stay_floor``, owner cells use the selected stayer rule: native DUE
    grandfathering with a separate death-solvency floor, or the legacy no
    cash-out/amortization experiment. Renter cells are unchanged.
    """
    if bool(getattr(P, "native_due_stayer_credit", False)) and not ctx.use_full_kernel:
        raise NotImplementedError("Native DUE saving requires the full owner kernel with explicit death solvency")
    natural_credit = bool(getattr(P, "native_solvency_credit", False))
    fixed_credit = fixed_unsecured_credit_active(P)
    if fixed_credit and natural_credit:
        raise ValueError("unsecured_credit_limit cannot be combined with native_solvency_credit")
    Nb = len(b_grid)
    I = P.I
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    beta = P.beta
    Rg = P.R_gross
    sigma = P.sigma
    alpha = P.alpha_cons
    oms = 1.0 - sigma
    b = ctx.b
    interp_method = str(getattr(P, "interp_method", "linear")).lower()
    use_value_kernel = NUMBA_AVAILABLE and interp_method == "linear"
    stay_orig_on = bool(getattr(P, "mortgage_origination_only", False))
    amort_rate = float(getattr(P, "mortgage_amortization", 0.0)) if stay_floor else 0.0
    due_stay = bool(stay_floor) and bool(getattr(P, "native_due_stayer_credit", False))
    purchase_saving_fraction = float(getattr(P, "experimental_purchase_saving_fraction", 1.0))
    if not 0.0 < purchase_saving_fraction <= 1.0:
        raise ValueError("purchase saving fraction must be in (0, 1]")
    stay_flag = int(bool(stay_floor) and not due_stay)
    stay_orig_flag = int(bool(stay_floor) and stay_orig_on)
    Vd = np.zeros((Nb, ctx.hcost.shape[1], I, npar, ncs))
    cd = np.zeros_like(Vd)
    hd = np.zeros_like(Vd)
    bd = np.zeros_like(Vd)
    Vo_nc = np.zeros((Nb, nc))
    bp_nc = np.zeros((Nb, nc))

    for i in range(I):
        yj = income_at_state(P, i, j, float(z_value))
        ri = r_hat[i]
        Rv = Rg * b + yj
        Rv_test = Rg * np.maximum(b, 0.0) + yj
        # Children-at-home earnings adjustment per family cell: working-age
        # earnings fall by the cell's penalty while fiscal transfers do not.
        # All zeros when the switch is off, so resources match bit for bit.
        pen_on = child_earnings_penalty_active(P)
        if pen_on and int(j) < int(getattr(P, "J_R", P.J)):
            base_earn = float(P.income[i, j]) * float(z_value)
            yadj_v = np.empty(nc)
            for _c in range(nc):
                _nn, _cs = decode_flat_family_state(_c, npar)
                _m = children_at_home_count(_nn, _cs, P)
                yadj_v[_c] = base_earn * (child_earnings_multiplier(P, int(j), _m) - 1.0)
            pen_flag = 1
        else:
            yadj_v = np.zeros(nc)
            pen_flag = 0
        hRmax = P.hR_max
        Vcr = flat_nc(Vc_arr[:, 0, i, :, :], Nb, nc)
        Rv1d_full = np.ascontiguousarray(Rv[:, 0])
        Rvt1d_full = np.ascontiguousarray(Rv_test[:, 0])
        wedge_on = rental_wedge_active(P)
        wedge_w0 = float(getattr(P, "rental_wedge_intercept", 0.0))
        wedge_w1 = float(getattr(P, "rental_wedge_slope", 0.0))
        wedge_hk = float(getattr(P, "rental_wedge_knee", 6.0))
        if wedge_on and ctx.exhaustive_saving:
            raise NotImplementedError("Rental wedge requires the golden-section renter block")
        if ctx.use_full_kernel:
            bp_prev_r = np.zeros((Nb, nc))
            has_prev_r = 0
            natural_floor, natural_dead = (native_solvency_support_floor(Vcr, b_grid)
                if natural_credit else (None, None))
            Vo_nc, bp_nc, co_nc, ho_nc = full_renter_block_kernel(
                Rv1d_full, Rvt1d_full, Vcr, bp_prev_r, has_prev_r, b_grid,
                ctx.cb_v, ctx.hb_v, ctx.psi_v_flat, ctx.gb_v, ctx.alpha_v, ctx.esc_v,
                ri, hRmax, P.c_min, P.c_bar_0, P.h_bar_0,
                alpha, oms, beta, s_next, D_next, ctx.gs_alpha1, ctx.gs_alpha2, ctx.gs_tol,
                int(ctx.exhaustive_saving),
                np.ascontiguousarray(yadj_v), pen_flag,
                int(wedge_on), wedge_w0, wedge_w1, wedge_hk,
                bool(getattr(P, "native_exact_allocation_output", False)),
                natural_floor,
                float(renter_floor[0]) if fixed_credit else -np.inf,
            )
            if natural_credit:
                Vo_nc[:, natural_dead] = -1e10
        else:
            Kr = (alpha**alpha * ((1 - alpha) / ri) ** (1 - alpha)) ** oms
            d_nc = SD.cb_flat + ri * SD.hb_flat
            if wedge_on:
                wedge_w0 = float(getattr(P, "rental_wedge_intercept", 0.0))
                wedge_w1 = float(getattr(P, "rental_wedge_slope", 0.0))
                wedge_hk = float(getattr(P, "rental_wedge_knee", 6.0))
            if pen_flag:
                Rv_nc = Rv + yadj_v.reshape(1, -1)
                Rv_test_nc = Rv_test + yadj_v.reshape(1, -1)
            else:
                Rv_nc = Rv
                Rv_test_nc = Rv_test
            Rv_eff_nc = Rv_nc + np.clip(SD.gb_flat - Rv_test_nc, 0.0, SD.gb_flat)
            cap_nc = ri * (hRmax - SD.hb_flat) / (1 - alpha)
            for c in range(nc):
                Vbar = Vcr[:, c]
                dc = d_nc[0, c]
                pc = SD.psi_flat[0, c]
                cc = cap_nc[0, c]
                cb_c = SD.cb_flat[0, c]
                hb_c = SD.hb_flat[0, c]
                ht_cap_c = max(hRmax - hb_c, 1e-10)
                lo = renter_floor.copy()
                if wedge_on:
                    wedge_cc = cb_c + hb_c * (ri + wedge_w0) + wedge_w1 * hb_c * max(hb_c - wedge_hk, 0.0)
                    hi = np.maximum(Rv_eff_nc[:, c] - wedge_cc - 1e-6, lo)
                    bp, val = golden_renter_wedge(
                        lo, hi, Rv_eff_nc[:, c], Vbar, b_grid, cb_c, hb_c, pc,
                        ri, hRmax, wedge_w0, wedge_w1, wedge_hk, alpha, oms, beta,
                        ctx.gs_alpha1, ctx.gs_alpha2, ctx.gs_tol, interp_method, 1.0,
                    )
                else:
                    hi = np.maximum(Rv_eff_nc[:, c] - dc - 1e-6, lo)
                    bp, val = golden_renter(
                        lo, hi, Rv_eff_nc[:, c], Vbar, b_grid, dc, pc, cc, cb_c, hb_c,
                        ri, hRmax, ht_cap_c, Kr, alpha, oms, beta,
                        ctx.gs_alpha1, ctx.gs_alpha2, ctx.gs_tol, interp_method, 1.0,
                    )
                bp_nc[:, c] = bp
                Vo_nc[:, c] = val
            if wedge_on:
                _, ct_nc, ht_nc = renter_wedge_flow_py(
                    Rv_eff_nc - bp_nc, SD.cb_flat, SD.hb_flat, ri,
                    wedge_w0, wedge_w1, wedge_hk, hRmax, alpha, oms, 1.0,
                )
                bad = (Rv_eff_nc - bp_nc - SD.cb_flat - (
                    SD.hb_flat * (ri + wedge_w0)
                    + wedge_w1 * SD.hb_flat * np.maximum(SD.hb_flat - wedge_hk, 0.0)
                )) <= 1e-10
                co_nc = SD.cb_flat + np.maximum(ct_nc, P.c_min)
                ho_nc = SD.hb_flat + np.maximum(ht_nc, 0.01)
                co_nc[bad] = P.c_bar_0 + P.c_min
                ho_nc[bad] = P.h_bar_0 + 0.01
            else:
                surplus_nc = Rv_eff_nc - d_nc - bp_nc
                ct_nc = alpha * np.maximum(surplus_nc, 1e-10)
                ht_nc = (1 - alpha) / ri * np.maximum(surplus_nc, 1e-10)
                cm = (SD.hb_flat + ht_nc) > hRmax
                if np.any(cm):
                    ct_cap = np.maximum(Rv_eff_nc - SD.cb_flat - ri * hRmax - bp_nc, 1e-10)
                    hcap = np.tile(np.maximum(hRmax - SD.hb_flat, 1e-10), (Nb, 1))
                    ct_nc[cm] = ct_cap[cm]
                    ht_nc[cm] = hcap[cm]
                co_nc = SD.cb_flat + np.maximum(ct_nc, P.c_min)
                ho_nc = SD.hb_flat + np.maximum(ht_nc, 0.01)
                bad = surplus_nc <= 1e-10
                co_nc[bad] = P.c_bar_0 + P.c_min
                ho_nc[bad] = P.h_bar_0 + 0.01
        Vd[:, 0, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
        bd[:, 0, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
        cd[:, 0, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
        hd[:, 0, i, :, :] = unflat_nc(ho_nc, Nb, npar, ncs)

        for ten in range(1, ctx.hcost.shape[1]):
            oc = ctx.ocst[i, ten]
            hsv = ctx.hsrv[i, ten]
            Vco = flat_nc(Vc_arr[:, ten, i, :, :], Nb, nc)
            if ctx.use_full_kernel:
                bf_v_base = ctx.bmo[i, ten, :, :].reshape(-1, order="F")
                bf_v = np.ascontiguousarray(
                    effective_owner_collateral_floor(P, bf_v_base, j)
                )
                if natural_credit:
                    bf_v, natural_dead = native_solvency_support_floor(Vco, b_grid)
                bp_prev_o = np.zeros((Nb, nc))
                has_prev_o = 0
                Vo_nc, bp_nc, co_nc = full_owner_block_kernel(
                    Rv1d_full, Rvt1d_full, Vco, bp_prev_o, has_prev_o, b_grid,
                    ctx.cb_v, ctx.hb_v, ctx.psi_v_flat, ctx.gb_v, ctx.alpha_v, ctx.esc_v, bf_v,
                    oc, hsv, ctx.owner_h_bar_scale, ctx.owner_service_premium, P.c_min,
                    alpha, oms, beta,
                    0.0 if getattr(P, "native_purchase_income", False) else s_next,
                    0.0 if getattr(P, "native_purchase_income", False) else D_next,
                    ctx.gs_alpha1, ctx.gs_alpha2, ctx.gs_tol,
                    ctx.strict_owner_hbar_feasibility, int(ctx.exhaustive_saving),
                    np.ascontiguousarray(yadj_v), pen_flag,
                    stay_flag, stay_orig_flag, amort_rate,
                    bool(getattr(P, "native_exact_allocation_output", False)),
                    due_stay,
                    native_due_death_floor(P, j, ctx.current_prices[i], P.H_own[ten - 1]) if due_stay else -np.inf,
                    1.0 if stay_floor else purchase_saving_fraction,
                )
                if natural_credit:
                    Vo_nc[:, natural_dead] = -1e10
            else:
                for c in range(nc):
                    Vbar = Vco[:, c]
                    cb_c = SD.cb_flat[0, c]
                    pc = SD.psi_flat[0, c]
                    owner_residual_h = hsv - ctx.owner_h_bar_scale * SD.hb_flat[0, c]
                    nn_c, cs_c = decode_flat_family_state(c, npar)
                    bf_c = ctx.bmo[i, ten, nn_c, cs_c]
                    lo = np.maximum(
                        owner_borrowing_floor(
                            P, b_grid, bf_c, j,
                            stay_on=stay_floor, stay_orig=stay_orig_on, amort=amort_rate,
                        ),
                        b_grid[0],
                    )
                    if not stay_floor and purchase_saving_fraction < 1.0:
                        quarter_floor = bf_c + (1.0 - purchase_saving_fraction) / purchase_saving_fraction * np.maximum(0.0, bf_c - b_grid)
                        lo = np.maximum(lo, quarter_floor)
                    hi = np.maximum(Rv_eff_nc[:, c] - oc - cb_c - 1e-6, lo)
                    if ctx.strict_owner_hbar_feasibility and owner_residual_h <= 0.0:
                        bp = lo.copy()
                        val = np.full(Nb, -1e10)
                    else:
                        Ko_c = (
                            ctx.owner_service_premium * max(owner_residual_h, 1e-10)
                        ) ** ((1 - alpha) * oms)
                        bp, val = golden_owner(
                            lo, hi, Rv_eff_nc[:, c], Vbar, b_grid, oc, cb_c, pc,
                            Ko_c, alpha, oms, beta, ctx.gs_alpha1, ctx.gs_alpha2, ctx.gs_tol, interp_method, 1.0,
                        )
                    bp_nc[:, c] = bp
                    Vo_nc[:, c] = val
                co_nc = SD.cb_flat + np.maximum(Rv_eff_nc - oc - SD.cb_flat - bp_nc, P.c_min)
            if ctx.exhaustive_saving and not bool(getattr(P, "native_exact_allocation_output", False)):
                if pen_flag:
                    res_base = Rv + yadj_v.reshape(1, -1)
                    res_test = Rv_test + yadj_v.reshape(1, -1)
                else:
                    res_base = Rv
                    res_test = Rv_test
                resources = res_base + np.clip(SD.gb_flat - res_test, 0.0, SD.gb_flat)
                co_nc = joint_nested.owner_consumption_from_solution(
                    resources, oc, bp_nc, SD.cb_flat, Vo_nc, co_nc)
            Vd[:, ten, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
            bd[:, ten, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
            cd[:, ten, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
    return Vd, cd, hd, bd


def _tenure_location_stage(
    Vd: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    ctx: SimpleNamespace,
    dp_choice: np.ndarray,
    Vd_stay: np.ndarray | None = None,
    bmo_purchase: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
    """Tenure choice + location logit for savings values ``Vd``.

    Returns ``(VH, tcj, prj_or_None, VI, lpj)`` where ``prj_or_None`` is the
    tenure-probability block when the tenure logit is active and ``None``
    otherwise (the caller then leaves ``tenure_probs`` untouched, as before).
    ``Vd_stay`` supplies the stayer (to == tn) values; ``None`` reads the
    standard ``Vd`` there, exactly as before.
    """
    Nb = len(b_grid)
    I = P.I
    npar = P.n_parity
    ncs = P.n_child_states
    nt = ctx.hcost.shape[1]
    purchase_floor = ctx.bmo if bmo_purchase is None else bmo_purchase
    transaction_support = bool(getattr(P, "native_purchase_income", False))
    require_owner_sale_solvency = fixed_unsecured_credit_active(P)
    birth_entry_grant = SD.birth_entry_grant
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    if Vd_stay is None:
        Vd_stay = Vd

    if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj, prj = tenure_logit_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, tenure_choice_kappa, Vd_stay, transaction_support, require_owner_sale_solvency
        )
        prj_full: np.ndarray | None = prj
    elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj = tenure_choice_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, Vd_stay, False, transaction_support, require_owner_sale_solvency
        )
        prj_full = None
    else:
        VH = np.zeros((Nb, nt, I, npar, ncs))
        tcj = np.zeros((Nb, nt, I, npar, ncs), dtype=np.int16)
        prj_full = np.zeros((Nb, nt, I, npar, ncs, nt)) if use_tenure_logit else None
        for id_ in range(I):
            for to in range(nt):
                sp = ctx.heq[id_, to] if to > 0 else 0.0
                Vopt = np.zeros((Nb, npar, ncs, nt))
                if to == 0:
                    Vopt[:, :, :, 0] = Vd[:, 0, id_, :, :]
                else:
                    ba = np.clip(b_grid + sp, b_grid[0], b_grid[-1])
                    Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)
                    if require_owner_sale_solvency:
                        Vopt[b_grid + sp < 0.0, :, :, 0] = -1e10
                for tn in range(1, nt):
                    hc = ctx.hcost[id_, tn]
                    Vow = Vd_stay[:, tn, id_, :, :] if to == tn else Vd[:, tn, id_, :, :]
                    if to == tn:
                        Vopt[:, :, :, tn] = Vow
                    elif to == 0:
                        bab = b_grid - hc
                        Vb = interp_on_grid(b_grid, Vow, bab)
                        for nn in range(npar):
                            for cs in range(ncs):
                                dpn = dp_choice[id_, tn, nn, cs]
                                bmn = ctx.bmo[id_, tn, nn, cs]
                                if SD.birth_dp[nn, cs, to, tn]:
                                    bag = np.maximum(bab, bmn)
                                    Vb[:, nn, cs] = interp_vector(b_grid, Vow[:, nn, cs], bag)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    Vg = interp_vector(b_grid, Vow[:, nn, cs], babg)
                                    inf_m = ((b_grid + gfix) < dpn) | (babg < bmn)
                                    Vg[inf_m] = -1e10
                                    Vb[:, nn, cs] = Vg
                                else:
                                    inf_m = (b_grid < dpn) | (bab < bmn)
                                    Vb[inf_m, nn, cs] = -1e10
                        Vopt[:, :, :, tn] = Vb
                    else:
                        bar = b_grid + sp - hc
                        Vrs = interp_on_grid(b_grid, Vow, bar)
                        for nn in range(npar):
                            for cs in range(ncs):
                                dpn = dp_choice[id_, tn, nn, cs]
                                bmn = ctx.bmo[id_, tn, nn, cs]
                                dpc = dpn - sp
                                inf_m = (b_grid < dpc) | (bar < bmn)
                                Vrs[inf_m, nn, cs] = -1e10
                        Vopt[:, :, :, tn] = Vrs
                tc = np.argmax(Vopt, axis=3)
                if use_tenure_logit:
                    ls, pr = logsumexp(Vopt / tenure_choice_kappa, axis=3)
                    pr[np.max(Vopt, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                    VH[:, to, id_, :, :] = tenure_choice_kappa * ls
                    assert prj_full is not None
                    prj_full[:, to, id_, :, :, :] = pr.astype(np.float32)
                else:
                    VH[:, to, id_, :, :] = np.max(Vopt, axis=3)
                tcj[:, to, id_, :, :] = tc

    kl = P.kappa_loc
    if NUMBA_AVAILABLE and bool(getattr(P, "use_loc_kernel", True)):
        VI, lpj = location_logit_kernel(VH, ctx.iidx, ctx.iwt, ctx.loc_shift, kl)
    else:
        VI = np.zeros((Nb, nt, I, npar, ncs))
        lpj = np.zeros((Nb, nt, I, I, npar, ncs))
        for io in range(I):
            for to in range(nt):
                Va = np.zeros((Nb, I, npar, ncs))
                Va[:, io, :, :] = VH[:, to, io, :, :]
                idx = ctx.iidx[:, io, to]
                wt = ctx.iwt[:, io, to]
                for id_ in range(I):
                    if id_ == io:
                        continue
                    Vdst = VH[:, 0, id_, :, :]
                    Va[:, id_, :, :] = (1 - wt)[:, None, None] * Vdst[idx, :, :] + wt[:, None, None] * Vdst[idx + 1, :, :]
                la = Va.copy()
                for id_ in range(I):
                    la[:, id_, :, :] += ctx.loc_shift[io, id_]
                la = la / kl
                ls, pr = logsumexp(la, axis=1)
                dead_loc = np.max(Va, axis=1) <= DEAD_VALUE_CUTOFF
                pr = np.where(dead_loc[:, None, :, :], 0.0, pr)
                VI[:, to, io, :, :] = kl * ls
                lpj[:, to, io, :, :, :] = pr
    return VH, tcj, prj_full, VI, lpj


def native_solvency_support_floor(values: np.ndarray, grid: np.ndarray):
    """Conservative feasible-node boundary; native value cutoff is approximate.

    This does not distinguish very negative finite utility from infeasibility.
    Reject holes instead of interpolating across infeasible continuation nodes.
    """
    v = np.asarray(values)
    bg = np.asarray(grid)
    if (v.ndim != 2 or v.shape[0] != bg.size or bg.size < 2
            or not np.isfinite(bg).all() or np.any(np.diff(bg) <= 0)
            or not np.isfinite(v).all()):
        raise ValueError("Invalid natural-credit continuation or grid")
    feasible = v > DEAD_VALUE_CUTOFF
    if np.any(feasible[:-1] & ~feasible[1:]):
        raise ValueError("Natural-credit feasible support must be an upper interval")
    dead = ~feasible.any(axis=0)
    floor = bg[np.argmax(feasible, axis=0)].astype(float)
    floor[dead] = bg[-1]
    return np.ascontiguousarray(floor), dead


def native_solvency_death_mask(grid, current_house_value, selling_cost):
    """Post-saving estate at the current decision price, net of liquidation."""
    return np.asarray(grid) + (1.0 - float(selling_cost)) * float(current_house_value) < 0.0


def native_solvency_continuation(next_values, probabilities, death_values, survival):
    """Expected continuation with infeasibility on every reachable branch.

    next_values has income first; caller chooses stationary or dated values.
    Death is valued at today's price, independently of next-date prices.
    """
    probabilities = np.asarray(probabilities)
    if (not 0 <= survival <= 1 or len(next_values) != len(probabilities)
            or not np.isfinite(probabilities).all() or np.any(probabilities < 0)
            or not np.isclose(probabilities.sum(), 1.0)):
        raise ValueError("Invalid natural-credit transition probabilities")
    values = np.zeros_like(death_values)
    bad = np.zeros_like(death_values, dtype=bool)
    for weight, continuation in zip(probabilities, next_values):
        if weight > 0.0:
            values += weight * continuation
            bad |= continuation <= DEAD_VALUE_CUTOFF
    values = survival * values + (1.0 - survival) * death_values
    bad = ((survival > 0.0) & bad) | ((survival < 1.0) & (death_values <= DEAD_VALUE_CUTOFF))
    values[bad] = -1e10
    return values


def validate_native_solvency_mode(P):
    """Default-off experimental mode, restricted to the reviewed household contract."""
    enabled = bool(getattr(P, "native_solvency_credit", False))
    if not enabled:
        return False
    required = ("native_purchase_income", "native_fixed_reference_entry",
                "native_explicit_transaction_grid", "native_exact_allocation_output",
                "exhaustive_saving_control")
    if (not all(bool(getattr(P, name, False)) for name in required)
            or parent_age_maturation_active(P) or mortgage_stay_floor_active(P)
            or bequest_utility_net_active(P) or estate_receiver_active(P)):
        raise ValueError("Natural credit requires reviewed native household mode without additional estate/stayer/maturation mechanisms")
    return True


def solve_bellman_full_markov_income(
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    continuation_V: np.ndarray | None = None,
):
    t0 = time.perf_counter()
    natural_credit = validate_native_solvency_mode(P)
    due_credit = bool(getattr(P, "native_due_stayer_credit", False))
    if due_credit and (natural_credit or not bool(getattr(P, "native_purchase_income", False))
                       or mortgage_stay_floor_active(P) or int(P.I) != 1
                       or not bool(getattr(P, "use_postdecision_current_distribution", True))):
        raise ValueError("DUE stayer mode requires one-market native accounting, no legacy stayer or natural-credit mode")
    purchase_income = bool(getattr(P, "native_purchase_income", False))
    if purchase_income and (
            bool(getattr(P, "joint_nested_choice", False))
            or bool(getattr(P, "use_pti_constraint", False))
            or float(P.lambda_d) != 0.0
            or not bool(getattr(P, "use_tenure_kernel", True))
            or not bool(getattr(P, "use_full_kernel", True)) or not NUMBA_AVAILABLE
            or str(getattr(P, "interp_method", "linear")) != "linear"
            or bool(getattr(P, "mortgage_origination_only", False))
            or float(getattr(P, "mortgage_amortization", 0.0)) != 0.0
            or np.any(np.asarray(getattr(P, "owner_ltv_multipliers", [1.0])) != 1.0)
            or child_earnings_penalty_active(P) or rental_wedge_active(P)
            or np.any(SD.birth_dp) or np.any(SD.birth_entry_grant)):
        raise ValueError("unsupported mechanism combined with native purchase-income timing")
    fec = get_fecundity_by_age(P)
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nh = P.n_house
    nt = 1 + nh
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    z_grid, z_weights, Pi_z = income_transition_values(P)
    Nz = len(z_grid)
    beta = P.beta
    Rg = P.R_gross
    sigma = P.sigma
    alpha = P.alpha_cons
    oms = 1.0 - sigma
    owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    owner_service_premium = max(float(getattr(P, "chi", 1.0)), 1e-8)
    strict_owner_hbar_feasibility = int(
        bool(getattr(P, "child_room_floor", False))
        and (
            float(getattr(P, "hbar_child_rooms", 0.0)) > 0.0
            or float(getattr(P, "hbar_first_child_jump", 0.0)) > 0.0
        )
    )
    owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
    owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
    owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    b = b_grid.reshape(-1, 1)
    interp_method = str(getattr(P, "interp_method", "linear")).lower()
    # Non-linear value interpolation only exists in the Python eval path, so a
    # non-linear interp method routes off the compiled full-Bellman kernel.
    use_full_kernel = NUMBA_AVAILABLE and bool(getattr(P, "use_full_kernel", True)) and interp_method == "linear"
    if str(getattr(P, "preference_spec", "stone_geary")) == "eqscale" and not use_full_kernel:
        raise NotImplementedError("eqscale preferences: full-kernel markov path only")
    cb_v = np.ascontiguousarray(SD.cb_flat.reshape(-1))
    hb_v = np.ascontiguousarray(SD.hb_flat.reshape(-1))
    psi_v_flat = np.ascontiguousarray(SD.psi_flat.reshape(-1))
    gb_v = np.ascontiguousarray(SD.gb_flat.reshape(-1))
    alpha_v = np.ascontiguousarray(SD.alpha_flat.reshape(-1))
    esc_v = np.ascontiguousarray(SD.escale_flat.reshape(-1))

    V = np.zeros((Nb, nt, I, J, Nz, npar, ncs))
    joint_active = bool(getattr(P, "joint_nested_choice", False))
    exhaustive_saving = joint_active or bool(getattr(P, "exhaustive_saving_control", False))
    if getattr(P, "two_shock_choice", False) and not joint_active:
        raise ValueError("Two-shock experiment requires joint choice mass accounting")
    if getattr(P, "fertility_nest_choice", False) and not joint_active:
        raise ValueError("Fertility-nest experiment requires joint choice mass accounting")
    if exhaustive_saving and not use_full_kernel:
        raise ValueError("Joint nested choice requires exhaustive compiled saving kernels")
    joint = joint_nested.allocate(V.shape, P) if joint_active else None
    if continuation_V is not None:
        continuation_V = np.asarray(continuation_V, dtype=float)
        if continuation_V.shape != V.shape:
            raise ValueError(
                "Calendar-time continuation value has shape "
                f"{continuation_V.shape}; expected {V.shape}."
            )
        if not np.all(np.isfinite(continuation_V)):
            raise ValueError("Calendar-time continuation value must be finite.")
    c_pol = np.zeros_like(V)
    hR_pol = np.zeros_like(V)
    bp_pol = np.ones_like(V)
    tenure_choice = np.zeros((Nb, nt, I, J, Nz, npar, ncs), dtype=np.int16)
    tenure_probs = (
        np.zeros((Nb, nt, I, J, Nz, npar, ncs, nt), dtype=float if joint_active else np.float32)
        if use_tenure_logit
        else None
    )
    loc_probs = np.zeros((Nb, nt, I, I, J, Nz, npar, ncs))
    fert_probs = np.zeros((Nb, nt, I, J, Nz, npar))
    # alt x attempting-parity slot: slot nn-1 holds the {stop, try}
    # probabilities of the upward attempt from parity nn (slot 0 = second
    # birth, slot 1 = third birth under n_parity=4).  The repaired child-count
    # mode adds a final at-home-count axis; historical mode keeps this shape.
    # Shape is unchanged for npar in {3, 4}; retained on P to avoid changing
    # the established Bellman return contract used by the GE loop.
    if independent_child_maturation_active(P):
        fert2_probs = np.zeros((Nb, nt, I, J, Nz, 2, max(npar - 2, 2), ncs))
    else:
        fert2_probs = np.zeros((Nb, nt, I, J, Nz, 2, max(npar - 2, 2)))
    fert_value = np.zeros((Nb, nt, I, J, Nz))

    ctx = _build_housing_stage_ctx(P, b_grid, SD, p_hat, use_full_kernel, exhaustive_saving)

    stay_active = mortgage_stay_floor_active(P) or due_credit
    if stay_active and joint_active:
        raise NotImplementedError("stayer mortgage floors: sequential Bellman path only")
    bp_pol_stay: np.ndarray | None = np.ones_like(bp_pol) if stay_active else None
    c_pol_stay = np.zeros_like(c_pol) if due_credit else None

    Vbq = np.zeros((Nb, nt, I, npar, ncs))
    for i in range(I):
        for ten in range(nt):
            # Gross by default (bitwise nesting); net of the selling cost when
            # the utility-side switch holds (alone or via the estate transfer).
            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
            if ten > 0 and bequest_utility_net_active(P):
                hv = estate_housing_value(P, float(p_hat[i]), float(P.H_own[ten - 1]), for_accounting=False)
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)
                    if natural_credit:
                        gross_house_value = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
                        dead_estate = native_solvency_death_mask(b_grid, gross_house_value, P.psi)
                        Vbq[dead_estate, ten, i, nn, cs] = -1e10

    for j in range(J - 1, -1, -1):
        in_fert = (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end)
        s_next = float(P.debt_taper_weights[j + 1])
        D_next = float(P.debt_caps[j + 1])
        renter_floor = np.maximum(renter_borrowing_floor(P, b_grid, j), b_grid[0])
        for zz, z_value in enumerate(z_grid):
            if j == J - 1:
                Vnr = Vbq
            else:
                Vnr = np.zeros((Nb, nt, I, npar, ncs))
                next_values = V if continuation_V is None else continuation_V
                for znext in range(Nz):
                    transition_weight = Pi_z[zz, znext]
                    if transition_weight > 0.0:
                        Vnr += transition_weight * next_values[
                            :, :, :, j + 1, znext, :, :
                        ]
                if bool(getattr(P, "use_age_survival", False)):
                    survival = float(P.survival_probs[j])
                    Vnr = survival * Vnr + (1.0 - survival) * Vbq
            if natural_credit and j < J - 1:
                survival = float(P.survival_probs[j]) if bool(getattr(P, "use_age_survival", False)) else 1.0
                # Income is the leading axis for the strict support operator.
                dated_next = np.moveaxis(next_values[:, :, :, j + 1, :, :, :], 3, 0)
                Vnr = native_solvency_continuation(dated_next, Pi_z[zz], Vbq, survival)
            Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs, age_index=j)
            if natural_credit:
                child_bad = apply_child_aging((Vnr <= DEAD_VALUE_CUTOFF).astype(float),
                    P, Nb, nt, I, npar, ncs, age_index=j) > 0.0
                Vc[child_bad] = -1e10
            # Parent-age newborn exemption (m-d): continuation with the
            # birth-period child safe.  At fertile ages the housing/saving +
            # tenure/location stages are re-solved under Vc_ex and the success
            # branch reads VI_ex at birth destinations; the wait branch and
            # all non-destination uses keep VI.  Constant mode skips this
            # (Vc_ex is None) bit for bit.
            Vc_ex = (
                apply_child_aging_exempt(Vnr, P, Nb, nt, I, npar, ncs, j)
                if parent_age_maturation_active(P)
                and independent_child_maturation_active(P)
                else None
            )
            Vd, cd, hd, bd = _savings_stage(
                Vc, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                s_next, D_next, renter_floor,
            )
            if stay_active:
                assert bp_pol_stay is not None
                Vd_s, cd_s, _, bd_s = _savings_stage(
                    Vc, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                    s_next, D_next, renter_floor, stay_floor=True,
                )
                bp_pol_stay[:, :, :, j, zz, :, :] = bd_s
                if c_pol_stay is not None:
                    c_pol_stay[:, :, :, j, zz, :, :] = cd_s
            else:
                Vd_s = Vd

            c_pol[:, :, :, j, zz, :, :] = cd
            hR_pol[:, :, :, j, zz, :, :] = hd
            bp_pol[:, :, :, j, zz, :, :] = bd

            dp_choice = ctx.dp_arr
            bmo_purchase = None
            if purchase_income:
                income_for_purchase = np.array([
                    income_at_state(P, i, j, float(z_value)) for i in range(I)
                ], dtype=float).reshape(I, 1, 1, 1) / Rg
                dp_choice = ctx.dp_arr - income_for_purchase
                bmo_purchase = np.maximum(ctx.bmo - income_for_purchase, b_grid[0])
            if natural_credit:
                dp_choice = np.full_like(ctx.dp_arr, -np.inf)
                bmo_purchase = np.full_like(ctx.bmo, b_grid[0])
            if bool(getattr(P, "use_pti_constraint", False)):
                income_j = np.array([income_at_state(P, i, j, float(z_value)) for i in range(I)], dtype=float)
                dp_choice = pti_adjusted_downpayment(ctx.dp_arr, ctx.hcost, income_j, P, b_grid)

            if joint_active:
                joint_result = joint_nested.bellman_block(
                    Vd, (b_grid, ctx.heq, ctx.hcost, dp_choice, ctx.bmo, SD.birth_dp, SD.birth_entry_grant),
                    P, j, fec, tenure_choice_kernel,
                )
                value, joint_prob, product, wait_prob = joint_result[:4]
                if getattr(P, "fertility_nest_choice", False):
                    joint.failure_probabilities[:, :, :, j, zz] = joint_result[4]
                V[:, :, :, j, zz] = value
                joint.probabilities[:, :, :, j, zz] = joint_prob
                joint.products[:, :, :, j, zz] = product
                joint.wait_probabilities[:, :, :, j, zz] = wait_prob
                # This wait-menu kernel is a fallback only. Every forward call
                # replaces it with exact joint-selected mass for its own pool.
                for product_index in range(nt):
                    tenure_probs[:, :, :, j, zz, :, :, product_index] = np.sum(
                        wait_prob * (product == product_index), axis=-1
                    )
                tenure_choice[:, :, :, j, zz] = (np.argmax(wait_prob, axis=-1)
                    if getattr(P, "fertility_nest_choice", False) else product[..., 0])
                loc_probs[:, :, 0, 0, j, zz] = (value[:, :, 0] > DEAD_VALUE_CUTOFF)
                fert_probs[:, :, :, j, zz, :2] = joint_nested.action_marginals(joint_prob[..., 0, 0, :, :])
                fert_value[:, :, :, j, zz] = value[..., 0, 0]
                for nn in range(1, npar - 1):
                    for cs in range(nn + 1):
                        fert2_probs[:, :, :, j, zz, :, nn - 1, cs] = joint_nested.action_marginals(joint_prob[..., nn, cs, :, :])
                continue

            VH, tcj, prj_full, VI, lpj = _tenure_location_stage(
                Vd, P, b_grid, SD, ctx, dp_choice, Vd_s, bmo_purchase,
            )
            if prj_full is not None:
                tenure_probs[:, :, :, j, zz, :, :, :] = prj_full
            tenure_choice[:, :, :, j, zz, :, :] = tcj
            loc_probs[:, :, :, :, j, zz, :, :] = lpj

            # Exact newborn-exempt success values: re-solve the housing/saving
            # + tenure/location stages under Vc_ex and read the success branch
            # off VI_ex at birth-destination states.  Only at fertile ages in
            # parent-age mode, so constant mode stays bitwise identical.
            VI_ex = None
            if (
                in_fert
                and Vc_ex is not None
                and bool(getattr(P, "sequential_births", False))
            ):
                Vd_ex, _, _, _ = _savings_stage(
                    Vc_ex, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                    s_next, D_next, renter_floor,
                )
                if stay_active:
                    Vd_ex_s, _, _, _ = _savings_stage(
                        Vc_ex, P, b_grid, SD, ctx, r_hat, j, float(z_value),
                        s_next, D_next, renter_floor, stay_floor=True,
                    )
                else:
                    Vd_ex_s = Vd_ex
                _, _, _, VI_ex, _ = _tenure_location_stage(
                    Vd_ex, P, b_grid, SD, ctx, dp_choice, Vd_ex_s, bmo_purchase,
                )

            if in_fert:
                pi_j = float(fec[j])
                if bool(getattr(P, "sequential_births", False)):
                    Vfa = np.empty((Nb, nt, I, 2))
                    settled_cs = readiness_settled_state(P)
                    Vfa[:, :, :, 0] = VI[:, :, :, 0, settled_cs]
                    # Exact exempt success value: re-optimized housing/saving +
                    # tenure/location value with the newborn safe next period.
                    if VI_ex is not None:
                        first_dest = VI_ex[:, :, :, 1, 1]
                    else:
                        first_dest = VI[:, :, :, 1, 1]
                    Vfa[:, :, :, 1] = (
                        pi_j
                        * (
                            first_dest
                            - float(P.first_birth_fixed_cost)
                        )
                        + (1.0 - pi_j) * VI[:, :, :, 0, settled_cs]
                    )
                    lf = Vfa / P.kappa_fert
                    ls, pr = logsumexp(lf, axis=3)
                    pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                    fert_probs[:, :, :, j, zz, :2] = pr
                    fert_value[:, :, :, j, zz] = P.kappa_fert * ls
                    if readiness_gate_active(P):
                        # Unsettled households cannot attempt a first birth.
                        # Settled households retain the existing entry logit.
                        V[:, :, :, j, zz, 0, 0] = VI[:, :, :, 0, 0]
                        V[:, :, :, j, zz, 0, 1] = fert_value[:, :, :, j, zz]
                    else:
                        V[:, :, :, j, zz, 0, 0] = fert_value[:, :, :, j, zz]
                    V[:, :, :, j, zz, 1:, :] = VI[:, :, :, 1:, :]
                    childless_copy_start = 2 if readiness_gate_active(P) else 1
                    V[:, :, :, j, zz, 0, childless_copy_start:] = VI[
                        :, :, :, 0, childless_copy_start:
                    ]
                    # Entry (childless wait/try) keeps kappa_fert; upward attempts
                    # at every parity use the continuation scale when set — margin-specific Gumbel scales on the same sequential choice tree.
                    kf_cont_raw = getattr(P, "kappa_fert_continuation", None)
                    kf_cont = float(P.kappa_fert) if kf_cont_raw is None else float(kf_cont_raw)
                    # Parity nn may try for birth nn+1.  Under the historical
                    # shared clock this is only child state 1.  Under the
                    # repaired specification, cs is the current at-home count
                    # and a birth maps (nn, cs) to (nn+1, cs+1).
                    for nn in range(1, npar - 1):
                        child_states = range(0, nn + 1) if independent_child_maturation_active(P) else (1,)
                        for cs in child_states:
                            destination_cs = birth_destination_child_state(P, cs)
                            V2 = np.empty((Nb, nt, I, 2))
                            V2[:, :, :, 0] = VI[:, :, :, nn, cs]
                            if VI_ex is not None:
                                cont_dest = VI_ex[:, :, :, nn + 1, destination_cs]
                            else:
                                cont_dest = VI[:, :, :, nn + 1, destination_cs]
                            V2[:, :, :, 1] = (
                                pi_j * cont_dest
                                + (1.0 - pi_j) * VI[:, :, :, nn, cs]
                            )
                            l2, p2 = logsumexp(V2 / kf_cont, axis=3)
                            p2[np.max(V2, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                            if independent_child_maturation_active(P):
                                fert2_probs[:, :, :, j, zz, :, nn - 1, cs] = p2
                            else:
                                fert2_probs[:, :, :, j, zz, :, nn - 1] = p2
                            V[:, :, :, j, zz, nn, cs] = kf_cont * l2
                else:
                    Vfa = np.zeros((Nb, nt, I, npar))
                    Vfa[:, :, :, 0] = VI[:, :, :, 0, 0]
                    if pi_j < 1.0:
                        for nn in range(1, npar):
                            Vfa[:, :, :, nn] = pi_j * VI[:, :, :, nn, 1] + (1.0 - pi_j) * VI[:, :, :, 0, 0]
                    else:
                        for nn in range(1, npar):
                            Vfa[:, :, :, nn] = VI[:, :, :, nn, 1]
                    lf = Vfa / P.kappa_fert
                    ls, pr = logsumexp(lf, axis=3)
                    pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                    fert_probs[:, :, :, j, zz, :] = pr
                    fert_value[:, :, :, j, zz] = P.kappa_fert * ls
                    V[:, :, :, j, zz, 0, 0] = fert_value[:, :, :, j, zz]
                    V[:, :, :, j, zz, 1:, :] = VI[:, :, :, 1:, :]
                    V[:, :, :, j, zz, 0, 1:] = VI[:, :, :, 0, 1:]
            else:
                V[:, :, :, j, zz, :, :] = VI

    P._fert2_probs = fert2_probs
    P._joint_choice = joint
    P._bp_pol_stay = bp_pol_stay
    P._c_pol_stay = c_pol_stay
    return (
        V,
        c_pol,
        hR_pol,
        bp_pol,
        tenure_choice,
        tenure_probs,
        loc_probs,
        fert_probs,
        fert_value,
        {"bellman": time.perf_counter() - t0},
    )


def renter_wedge_flow_py(S_raw, cbc, hbc, ri, w0, w1, hk, hRmax, alpha, oms, es=1.0):
    """Numpy mirror of kernels.renter_wedge_flow (see its docstring for the
    economics). All inputs broadcastable arrays or scalars; returns
    ``(u_flow, ct, ht)`` with infeasible cells marked by u_flow = -1e10 and
    (ct, ht) = 0. Cross-tested against the numba version."""
    S_raw = np.asarray(S_raw, dtype=float)
    cbc = np.asarray(cbc, dtype=float)
    hbc = np.asarray(hbc, dtype=float)
    ri1 = float(ri) + float(w0)
    w1 = float(w1)
    hk = float(hk)
    S_raw, cbc, hbc = np.broadcast_arrays(S_raw, cbc, hbc)
    Cc = hbc * ri1 + w1 * hbc * np.maximum(hbc - hk, 0.0)
    S = S_raw - cbc - Cc
    feasible = S > 1e-10
    u = np.full_like(S, -1e10)
    ct = np.zeros_like(S)
    ht = np.zeros_like(S)
    if not np.any(feasible):
        return u, ct, ht
    Sf = S[feasible]
    hbc_f = hbc[feasible]
    cbc_f = cbc[feasible]
    S_raw_f = S_raw[feasible]
    Kr1 = (alpha**alpha * ((1.0 - alpha) / ri1) ** (1.0 - alpha)) ** oms
    ht_a = (1.0 - alpha) * Sf / ri1
    h_a = hbc_f + ht_a
    ua = Kr1 * np.maximum(Sf, 1e-10) ** oms / oms
    if es != 1.0:
        ua = es * ua
    valid_a = h_a <= hk
    u_f = np.where(valid_a, ua, -1e10)
    ct_f = np.where(valid_a, alpha * Sf, 0.0)
    ht_f = np.where(valid_a, ht_a, 0.0)
    # Kink bunching candidate.
    ht_k = hk - hbc_f
    valid_k = (hbc_f < hk) & (hk <= float(hRmax)) & (Sf >= ri1 * ht_k)
    ct_k = np.maximum(Sf - ri1 * ht_k, 1e-10)
    u_k = (ct_k**alpha * np.maximum(ht_k, 1e-10) ** (1.0 - alpha)) ** oms / oms
    if es != 1.0:
        u_k = es * u_k
    take_k = valid_k & (u_k > u_f)
    u_f = np.where(take_k, u_k, u_f)
    ct_f = np.where(take_k, ct_k, ct_f)
    ht_f = np.where(take_k, ht_k, ht_f)
    # Above-knee quadratic candidate.
    S_eff = Sf + np.where(hbc_f < hk, w1 * hbc_f * (hk - hbc_f), 0.0)
    A2 = ri1 + w1 * (2.0 * hbc_f - hk)
    quad_b = w1 * (1.0 + alpha)
    if quad_b > 0.0:
        disc = A2 * A2 + 4.0 * quad_b * (1.0 - alpha) * S_eff
        ht_b = (-A2 + np.sqrt(np.maximum(disc, 0.0))) / (2.0 * quad_b)
    else:
        ht_b = (1.0 - alpha) * S_eff / A2
    h_b = hbc_f + ht_b
    valid_b = h_b > hk
    ct_b = np.maximum(S_eff - A2 * ht_b - w1 * ht_b**2, 1e-10)
    u_b = (ct_b**alpha * np.maximum(ht_b, 1e-10) ** (1.0 - alpha)) ** oms / oms
    if es != 1.0:
        u_b = es * u_b
    take_b = valid_b & (u_b > u_f)
    u_f = np.where(take_b, u_b, u_f)
    ct_f = np.where(take_b, ct_b, ct_f)
    ht_f = np.where(take_b, ht_b, ht_f)
    # Cap constraint.
    h_best = hbc_f + ht_f
    over_cap = h_best > float(hRmax)
    if np.any(over_cap):
        Ccap = float(hRmax) * ri1 + w1 * float(hRmax) * max(float(hRmax) - hk, 0.0)
        ct_cap = np.maximum(S_raw_f - cbc_f - Ccap, 1e-10)
        ht_use = np.maximum(float(hRmax) - hbc_f, 1e-10)
        u_cap = (ct_cap**alpha * ht_use ** (1.0 - alpha)) ** oms / oms
        if es != 1.0:
            u_cap = es * u_cap
        u_f = np.where(over_cap, u_cap, u_f)
        ct_f = np.where(over_cap, ct_cap, ct_f)
        ht_f = np.where(over_cap, ht_use, ht_f)
    u[feasible] = u_f
    ct[feasible] = ct_f
    ht[feasible] = ht_f
    return u, ct, ht


def eval_renter_wedge(bp, Rv, Vbar, b_grid, cb_c, hb_c, pc, ri, hRmax, w0, w1, hk,
                      alpha, oms, beta, vinterp=None, es=1.0):
    """Renter objective at savings bp under the rental wedge (Python path)."""
    if vinterp is None:
        vinterp = make_value_interp(b_grid, Vbar, "linear")
    u_flow, _, _ = renter_wedge_flow_py(Rv - bp, cb_c, hb_c, ri, w0, w1, hk, hRmax, alpha, oms, es)
    f = u_flow + pc + beta * vinterp(bp)
    bad = u_flow <= -1e9
    if np.any(bad):
        f = np.where(bad, -1e10, f)
    return f


def golden_renter_wedge(lo, hi, Rv, Vbar, b_grid, cb_c, hb_c, pc, ri, hRmax, w0, w1, hk,
                        alpha, oms, beta, a1, a2, tol, method="linear", es=1.0):
    """Golden-section renter savings search under the rental wedge."""
    vinterp = make_value_interp(b_grid, Vbar, method)
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_renter_wedge(x1, Rv, Vbar, b_grid, cb_c, hb_c, pc, ri, hRmax, w0, w1, hk, alpha, oms, beta, vinterp, es)
    f2 = eval_renter_wedge(x2, Rv, Vbar, b_grid, cb_c, hb_c, pc, ri, hRmax, w0, w1, hk, alpha, oms, beta, vinterp, es)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.clip(np.where(bt, x2 + d, x1 - d), lo, hi)
        fe = eval_renter_wedge(xe, Rv, Vbar, b_grid, cb_c, hb_c, pc, ri, hRmax, w0, w1, hk, alpha, oms, beta, vinterp, es)
        x1n = np.where(bt, x2, xe)
        f1n = np.where(bt, f2, fe)
        x2n = np.where(bt, xe, x1)
        f2n = np.where(bt, fe, f1)
        d = d * a2
        x1, x2, f1, f2 = x1n, x2n, f1n, f2n
    bt = f2 >= f1
    return np.where(bt, x2, x1), np.maximum(f1, f2)


def golden_renter(lo, hi, Rv, Vbar, b_grid, dc, pc, cc, cb_c, hb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, a1, a2, tol, method="linear", es=1.0):
    vinterp = make_value_interp(b_grid, Vbar, method)
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_renter(x1, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, vinterp, es)
    f2 = eval_renter(x2, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, vinterp, es)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.clip(np.where(bt, x2 + d, x1 - d), lo, hi)
        fe = eval_renter(xe, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, vinterp, es)
        x1n = np.where(bt, x2, xe)
        f1n = np.where(bt, f2, fe)
        x2n = np.where(bt, xe, x1)
        f2n = np.where(bt, fe, f1)
        d = d * a2
        x1, x2, f1, f2 = x1n, x2n, f1n, f2n
    bt = f2 >= f1
    return np.where(bt, x2, x1), np.maximum(f1, f2)


def eval_renter(bp, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, vinterp=None, es=1.0):
    if vinterp is None:
        vinterp = make_value_interp(b_grid, Vbar, "linear")
    surplus = Rv - dc - bp
    ss = np.maximum(surplus, 1e-10)
    if es != 1.0:
        f = es * Kr * ss ** oms / oms + pc + beta * vinterp(bp)
    else:
        f = Kr * ss ** oms / oms + pc + beta * vinterp(bp)
    cm = surplus > cc
    if np.any(cm):
        ct = np.maximum(Rv[cm] - cb_c - ri * hRmax - bp[cm], 1e-10)
        if es != 1.0:
            f[cm] = es * (ct**alpha * ht_cap_c ** (1 - alpha)) ** oms / oms + pc + beta * vinterp(bp[cm])
        else:
            f[cm] = (ct**alpha * ht_cap_c ** (1 - alpha)) ** oms / oms + pc + beta * vinterp(bp[cm])
    f[surplus <= 1e-10] = -1e10
    return f


def golden_owner(lo, hi, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, a1, a2, tol, method="linear", es=1.0):
    vinterp = make_value_interp(b_grid, Vbar, method)
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_owner(x1, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, vinterp, es)
    f2 = eval_owner(x2, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, vinterp, es)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.clip(np.where(bt, x2 + d, x1 - d), lo, hi)
        fe = eval_owner(xe, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, vinterp, es)
        x1n = np.where(bt, x2, xe)
        f1n = np.where(bt, f2, fe)
        x2n = np.where(bt, xe, x1)
        f2n = np.where(bt, fe, f1)
        d = d * a2
        x1, x2, f1, f2 = x1n, x2n, f1n, f2n
    bt = f2 >= f1
    return np.where(bt, x2, x1), np.maximum(f1, f2)


def eval_owner(bp, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, vinterp=None, es=1.0):
    if vinterp is None:
        vinterp = make_value_interp(b_grid, Vbar, "linear")
    ct_raw = Rv - oc - cb_c - bp
    ct = np.maximum(ct_raw, 1e-10)
    if es != 1.0:
        f = es * Ko_c * ct ** (alpha * oms) / oms + pc + beta * vinterp(bp)
    else:
        f = Ko_c * ct ** (alpha * oms) / oms + pc + beta * vinterp(bp)
    f[ct_raw <= 1e-10] = -1e10
    return f


def build_forward_tenure_transition_maps(
    P: SimpleNamespace,
    b_grid: np.ndarray,
    hc: np.ndarray,
    he: np.ndarray,
    phi_choice: np.ndarray,
    birth_dp: np.ndarray,
    birth_entry_grant: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map pre-transaction wealth into the conditional tenure branch.

    The map mirrors the Bellman tenure-choice accounting. In particular, a
    fixed birth-linked entry grant is added only when a renter buys, before
    the owner borrowing floor is applied. With a zero grant this is exactly
    the historical transaction map.
    """

    if bool(getattr(P, "native_purchase_income", False)):
        if np.any(birth_dp) or np.any(birth_entry_grant):
            raise ValueError("native purchase-income maps do not combine grants or waivers")
        bg = np.asarray(b_grid, dtype=float)
        I, nt = hc.shape
        npar, ncs = phi_choice.shape[2:]
        idx = np.zeros((I, nt, nt, npar, ncs, len(bg)), dtype=np.int64)
        wt = np.zeros_like(idx, dtype=float)
        for i in range(I):
            for old in range(nt):
                for new in range(nt):
                    # No collateral clipping: income enters spending once, later.
                    x = bg if old == new else bg + he[i, old] - hc[i, new]
                    ii, ww = interp_indices(bg, np.clip(x, bg[0], bg[-1]))
                    idx[i, old, new, :, :, :] = ii
                    wt[i, old, new, :, :, :] = ww
        return idx, wt
    bg = np.asarray(b_grid, dtype=float)
    bmin = float(bg[0])
    bmax = float(bg[-1])
    I, nt = hc.shape
    npar = phi_choice.shape[2]
    ncs = phi_choice.shape[3]
    Nb = bg.size
    tmx_idx = np.zeros((I, nt, nt, npar, ncs, Nb), dtype=np.int64)
    tmx_wt = np.zeros((I, nt, nt, npar, ncs, Nb))
    use_grant = bool(getattr(P, "propagate_birth_entry_grant", True))

    for nn in range(npar):
        for cs in range(ncs):
            for id_ in range(I):
                for to in range(nt):
                    sale_proceeds = he[id_, to]
                    for tn in range(nt):
                        financed_share = phi_choice[id_, tn, nn, cs] if tn > 0 else 1.0
                        if tn == to:
                            branch_wealth = bg.copy()
                        elif to == 0 and tn > 0:
                            purchase_cost = hc[id_, tn]
                            # Bellman precedence: the down-payment waiver branch
                            # wins when both policy hooks are configured.
                            grant = (
                                birth_entry_grant[id_, tn, nn, cs]
                                if use_grant and not birth_dp[nn, cs, to, tn]
                                else 0.0
                            )
                            branch_wealth = np.clip(
                                np.maximum(bg - purchase_cost + grant, -financed_share * purchase_cost),
                                bmin,
                                bmax,
                            )
                        elif to > 0 and tn == 0:
                            # Post-liquidation wealth is a renter position, so
                            # any residual unsecured debt must carry through.
                            branch_wealth = np.clip(bg + sale_proceeds, bmin, bmax)
                        else:
                            purchase_cost = hc[id_, tn]
                            branch_wealth = np.clip(
                                np.maximum(
                                    bg + sale_proceeds - purchase_cost,
                                    -financed_share * purchase_cost,
                                ),
                                bmin,
                                bmax,
                            )
                        tmx_idx[id_, to, tn, nn, cs, :], tmx_wt[id_, to, tn, nn, cs, :] = interp_indices(
                            bg, branch_wealth
                        )
    return tmx_idx, tmx_wt


def _child_Pa_for_age(P: SimpleNamespace, Pia, age_index) -> np.ndarray:
    """Standard child transition for age ``j`` under parent_age, else ``Pia``.

    In constant mode this returns ``Pia`` untouched, so every existing code
    path is bitwise identical.  In parent_age mode it returns the
    age-specific standard matrix ``P.Pi_child_by_age[j]`` (draw on all ``m``
    children at home).  A ``None`` age falls back to ``Pia``.
    """
    if (
        age_index is not None
        and parent_age_maturation_active(P)
        and independent_child_maturation_active(P)
    ):
        by_age = getattr(P, "Pi_child_by_age", None)
        if by_age is not None:
            return np.asarray(by_age[int(age_index)])
    return Pia


def _child_Pa_exempt_for_age(P: SimpleNamespace, age_index):
    """Newborn-exempt child transition for age ``j`` (parent_age only).

    Returns ``P.Pi_child_exempt_by_age[j]`` (draw on ``m - d`` with the
    birth-period child safe), or ``None`` in constant mode.
    """
    if (
        age_index is not None
        and parent_age_maturation_active(P)
        and independent_child_maturation_active(P)
    ):
        ex_age = getattr(P, "Pi_child_exempt_by_age", None)
        if ex_age is not None:
            return np.asarray(ex_age[int(age_index)])
    return None


def _blended_child_Pi_for_cell(
    P: SimpleNamespace,
    Pia,
    age_index,
    newborn_frac: float,
) -> np.ndarray:
    """Blend standard and exempt rows by the newborn share ``f`` of a cell.

    Implements the m-d exemption without a one-bit state expansion: the
    fraction ``f`` of post-birth mass that just arrived faces the exempt
    row, the remainder faces the standard row.  The blend is exact for
    cell totals and entrant flows (both linear in mass); only the
    within-cell wealth split is pooled.  With ``f = 0`` (or constant mode)
    this equals the standard matrix bit for bit.
    """
    Pi_std = np.asarray(_child_Pa_for_age(P, Pia, age_index))
    if not parent_age_maturation_active(P):
        return Pi_std
    f = float(np.clip(float(newborn_frac), 0.0, 1.0))
    if f <= 0.0:
        return Pi_std
    Pi_ex = _child_Pa_exempt_for_age(P, age_index)
    if Pi_ex is None:
        return Pi_std
    return (1.0 - f) * Pi_std + f * np.asarray(Pi_ex)


def apply_child_aging(Vn, P, Nb, nt, I, npar, ncs, age_index=None):
    Vc = np.zeros((Nb, nt, I, npar, ncs))
    K = P.n_child_stages
    if P.use_stochastic_aging and hasattr(P, "Pi_child"):
        Pa = _child_Pa_for_age(P, P.Pi_child, age_index)
        for nn in range(npar):
            if readiness_gate_active(P) and nn == 0 and age_index is not None:
                current_age = float(P.age_start) + float(age_index) * float(P.da)
                next_age = current_age + float(P.da)
                hazard = readiness_transition_hazard(P, current_age, next_age)
                Vc[:, :, :, 0, 0] = (
                    (1.0 - hazard) * Vn[:, :, :, 0, 0]
                    + hazard * Vn[:, :, :, 0, 1]
                )
                Vc[:, :, :, 0, 1] = Vn[:, :, :, 0, 1]
                if ncs > 2:
                    Vc[:, :, :, 0, 2:] = Vn[:, :, :, 0, 2:]
                continue
            Pi = Pa[:, :, nn]
            Vnn = np.reshape(Vn[:, :, :, nn, :], (-1, ncs), order="F")
            Vc[:, :, :, nn, :] = np.reshape(Vnn @ Pi.T, (Nb, nt, I, ncs), order="F")
    else:
        csm1 = K + 1
        csm2 = K + 2
        for nn in range(npar):
            for cs in range(ncs):
                if cs == 0:
                    csn = 0
                elif cs >= csm1:
                    csn = cs
                elif cs < K:
                    csn = cs + 1
                else:
                    csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                Vc[:, :, :, nn, cs] = Vn[:, :, :, nn, csn]
    return Vc


def apply_child_aging_exempt(Vn, P, Nb, nt, I, npar, ncs, age_index):
    """Newborn-exempt continuation for parent_age Bellman birth values.

    Applies the exempt matrix (draw on ``m - d``) to every cell.  Callers
    use the resulting ``Vc_ex`` only for destination states reached by a
    current-period birth; stayer states keep the standard ``Vc``.  This is
    the m-d implementation (no one-bit state expansion: the flag would
    double the child state and every downstream policy array).  In constant
    mode there is no exempt matrix and this falls back to ``apply_child_aging``.
    """
    Pa_ex = _child_Pa_exempt_for_age(P, age_index)
    if Pa_ex is None:
        return apply_child_aging(Vn, P, Nb, nt, I, npar, ncs, age_index=age_index)
    Vc = np.zeros((Nb, nt, I, npar, ncs))
    for nn in range(npar):
        Pi = np.asarray(Pa_ex)[:, :, nn]
        Vnn = np.reshape(Vn[:, :, :, nn, :], (-1, ncs), order="F")
        Vc[:, :, :, nn, :] = np.reshape(Vnn @ Pi.T, (Nb, nt, I, ncs), order="F")
    return Vc


def bequest_utility_vec(b, nk, P):
    b_gross = np.maximum(b, 0.0)
    estate_tax_rate = min(max(float(getattr(P, "estate_tax_rate", 0.0)), 0.0), 0.999)
    estate_tax_exemption = max(float(getattr(P, "estate_tax_exemption", 0.0)), 0.0)
    taxable = np.maximum(b_gross - estate_tax_exemption, 0.0)
    b = np.maximum(b_gross - estate_tax_rate * taxable, 0.0)
    spec = str(getattr(P, "bequest_spec", "linear_child_scale")).strip().lower()
    utility_wealth = b
    if spec in {"linear_child_scale", "legacy", "current"}:
        scale = P.theta0 * max(1 + P.theta_n * nk, 0)
    elif spec in {"parent_gated_luxury", "parent_gated"}:
        scale = P.theta0 if int(nk) >= 1 else 0.0
    elif spec in {"equal_division_luxury", "equal_division"}:
        n_children = int(nk)
        scale = P.theta0 * n_children if n_children >= 1 else 0.0
        utility_wealth = b / n_children if n_children >= 1 else np.zeros_like(b)
    else:
        raise ValueError(f"unknown bequest_spec: {spec}")
    normalize = bool(getattr(P, "normalize_bequest_utility", False))
    if spec in {
        "parent_gated_luxury",
        "parent_gated",
        "equal_division_luxury",
        "equal_division",
    } and not normalize:
        raise ValueError("parent-gated and equal-division bequests require zero-estate normalization")
    if abs(P.sigma - 1) < 1e-6:
        utility = np.log(P.theta1 + utility_wealth)
        if normalize:
            utility = utility - np.log(P.theta1)
        return scale * utility
    utility = (P.theta1 + utility_wealth) ** (1 - P.sigma) / (1 - P.sigma)
    if normalize:
        utility = utility - P.theta1 ** (1 - P.sigma) / (1 - P.sigma)
    return scale * utility


def pti_adjusted_downpayment(dp_arr, hcost, income, P, b_grid):
    """Optional underwriting screen based on actual transaction debt.

    The core collateral constraint is `dp_arr`: cash available before purchase
    must cover `(1 - phi) * pH`, and branch liquid wealth may not fall below
    `-phi * pH`. If PTI is enabled, use actual transaction debt
    `D=max(pH-cash, 0)`, not the maximum allowed LTV debt, so extra cash can
    relax the payment screen.
    """
    out = np.array(dp_arr, dtype=float, copy=True)
    pti_limit = max(float(getattr(P, "pti_limit", 1.0)), 0.0)
    q = max(float(getattr(P, "q", 0.0)), 0.0)
    tau_h = max(float(getattr(P, "tau_H", 0.0)), 0.0)
    incomes = np.asarray(income, dtype=float).reshape(-1)
    big_dp = max(float(b_grid[-1]) + 10.0 * np.max(np.maximum(hcost, 1.0)), 1e8)
    for i in range(out.shape[0]):
        y = float(incomes[i]) if i < incomes.size and np.isfinite(incomes[i]) else 0.0
        for ten in range(1, out.shape[1]):
            house_cost = float(hcost[i, ten])
            tax_payment = tau_h * house_cost
            allowed_debt_payment = pti_limit * y - tax_payment
            if allowed_debt_payment < 0.0:
                out[i, ten, :, :] = big_dp
            elif q > 1e-12:
                max_debt_pti = allowed_debt_payment / q
                min_cash_pti = house_cost - max_debt_pti
                out[i, ten, :, :] = np.maximum(out[i, ten, :, :], min_cash_pti)
    return out


def current_child_bin_dt(nn, cs, dep_last, high_cutoff=2, child_state_mode="shared_clock"):
    """Family-size room bin: 2 = no dependent child present, 3 = small
    family, 4 = large family.  ``high_cutoff`` is the parity at which the
    large-family bin starts (default 2 = legacy parity bins; 3 under the
    literal-parity convention, where bin 3 is "1-2 children" and bin 4 is
    "3+")."""
    if str(child_state_mode).strip().lower() == "independent_count":
        current_n = max(int(cs), 0) if int(cs) <= int(nn) else 0
        if current_n <= 0:
            return 2
        return 3 if current_n < high_cutoff else 4
    if cs == 0 or cs > dep_last:
        return 2
    current_n = max(nn, 0)
    if current_n <= 0:
        return 2
    if current_n < high_cutoff:
        return 3
    return 4
