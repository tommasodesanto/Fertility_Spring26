# Selected frozen implementation excerpts

These are read-only excerpts from the Torch source root used by the authenticated block0506 reference. The complete parent-file SHA256 identities are recorded in calculation_receipt.json. Ellipses between sections indicate omitted implementation, not a complete solver. Optional branches appear in source; use the saved parameter table and the prompt to determine which are active. No new model import, solve or source modification was performed to assemble this file.


## code/model/intergen_eqscale_seq_optimized/child_preferences.py (lines 1-58)

```python
"""Explicit optional child-benefit curvature and compensated housing shares."""
from __future__ import annotations

import math
import numpy as np


def apply_child_preferences(P, alpha, benefit, material_multiplier):
    """Apply the declared specification before native family-type compression.

    P.psi_child stores b in b*m**(1-kappa), where m is children at home.
    With compensated shares, material utility is CRRA of A(m)*Q/e(m),
    Q=c**alpha(m)*s**(1-alpha(m)). A=K(alpha0,r*)/K(alpha(m),r*) uses
    a fixed reference rent, never the current equilibrium rent.
    Absent both options, no array or floating-point operation is changed.
    """
    curvature = float(getattr(P, "child_benefit_curvature", 0.0))
    compensated = bool(getattr(P, "compensated_child_housing_shares", False))
    if not math.isfinite(curvature) or not 0 <= curvature < 1:
        raise ValueError("Child-benefit curvature must be finite and in [0,1)")
    if curvature == 0 and not compensated:
        return
    if str(getattr(P, "child_state_mode", "")) != "independent_count":
        raise ValueError("Native child preferences require independent children-at-home states")
    if getattr(P, "utility_comparison_arm", None) is not None:
        raise ValueError("Choose native child preferences or the legacy comparison adapter, not both")
    if curvature != 0:
        for n in range(int(P.n_parity)):
            for m in range(1, min(n + 1, int(P.n_child_states))):
                benefit[n, m] = P.psi_child * float(m) ** (1.0 - curvature)
    if not compensated:
        return
    if (str(getattr(P, "preference_spec", "")) != "eqscale"
            or str(getattr(P, "eqscale_form", "")) != "power"
            or bool(getattr(P, "child_room_floor", False))
            or float(getattr(P, "hbar_first_child_jump", 0.0)) != 0
            or float(getattr(P, "hbar_child_rooms", 0.0)) != 0
            or float(P.delta_alpha) != 0):
        raise ValueError("Compensated first-child shares require power equivalence scale, no floor and no later loading")
    base_alpha = float(P.alpha_cons)
    loading = float(P.delta_alpha_jump)
    reference_rent = float(getattr(P, "utility_reference_rent", math.nan))
    if (not math.isfinite(reference_rent) or reference_rent <= 0
            or not math.isfinite(base_alpha) or not math.isfinite(loading)
            or loading < 0 or not .05 <= base_alpha - loading <= base_alpha <= .95):
        raise ValueError("Explicit positive reference rent and unclipped interior housing shares required")
    log_k = alpha * np.log(alpha) + (1 - alpha) * np.log((1 - alpha) / reference_rent)
    log_k0 = (base_alpha * math.log(base_alpha)
              + (1 - base_alpha) * math.log((1 - base_alpha) / reference_rent))
    factors = np.where(alpha == base_alpha, 1.0, np.exp(log_k0 - log_k))
    material_multiplier *= factors ** (1.0 - float(P.sigma))
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 137-195)

```python
def renter_borrowing_floor(P: SimpleNamespace, b: Any, j: int) -> np.ndarray:
    """Renter floor; all renter debt is unsecured."""

    return debt_rule_at_age(P, b, j)


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

```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 2549-2618)

```python
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            if independent_child_maturation_active(P):
                nk = cs if cs <= nn else 0
                kp = nk > 0
            else:
                nk = nn
                kp = (cs >= 1) and (cs < csm1)
            if readiness_gate_active(P) and nn == 0 and cs == 1:
                # E6c reuses the otherwise invalid (childless, cs=1) cell for
                # the settled state. It has childless preferences, not the
                # child-at-home consumption and housing adjustments.
                kp = False
            if kp:
                c_bar[nn, cs] = P.c_bar_0 + P.c_bar_n * nk
                if str(P.child_housing_spec).lower() == "linear_only":
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_n * nk
                else:
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * nk
                psi_v[nn, cs] = P.psi_child * nk
                g_bar[nn, cs] = g0 + gn * nk
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    if child_room_floor_active:
                        h_bar[nn, cs] = (
                            float(P.hbar_first_child_jump)
                            + float(P.hbar_child_rooms) * nk
                        )
                        alpha_bar[nn, cs] = P.alpha_cons
                    else:
                        h_bar[nn, cs] = 0.0
                        alpha_bar[nn, cs] = np.clip(
                            P.alpha_cons - (P.delta_alpha_jump + P.delta_alpha * nk), 0.05, 0.95
                        )
                    if eqscale_form == "power":
                        # Imposed Scholz-Seshadri-Khitatrakun (2006, JPE 114(4), p.619;
                        # Citro-Michael 1995) scale relative to a childless couple:
                        #   e(n) = ((2 + 0.7 n)/2)**0.7.
                        # Flow utility is multiplied by escale, u = escale * x**(1-sigma)/(1-sigma),
                        # while per-equivalent CRRA utility is u(x/e) = e**(sigma-1) * x**(1-sigma)/(1-sigma),
                        # so the multiplier is e**(sigma-1); at the baseline sigma = 2 the
                        # multiplier equals the scale itself. n is the literal parity state
                        # (under L4, nn = 3 is the top-coded 3+ bin, scaled at n = 3).
                        escale[nn, cs] = (((2.0 + 0.7 * nk) / 2.0) ** 0.7) ** (float(P.sigma) - 1.0)
                    elif eqscale_form == "sqrt":
                        # Declared robustness alternative: square-root household-size scale,
                        # e(n) = sqrt((2 + n)/2) relative to a childless couple.
                        escale[nn, cs] = (((2.0 + nk) / 2.0) ** 0.5) ** (float(P.sigma) - 1.0)
                    else:
                        escale[nn, cs] = 1.0 + P.gamma_e * nk
            else:
                c_bar[nn, cs] = P.c_bar_0
                h_bar[nn, cs] = P.h_bar_0
                g_bar[nn, cs] = g0
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    h_bar[nn, cs] = 0.0

    apply_child_preferences(P, alpha_bar, psi_v, escale)
    triples = np.column_stack(
        [
            c_bar.reshape(-1, order="F"),
            h_bar.reshape(-1, order="F"),
            psi_v.reshape(-1, order="F"),
        ]
    )
    unique_triples, type_map = np.unique(triples, axis=0, return_inverse=True)
    birth_dp = np.zeros((P.n_parity, P.n_child_states, 1 + P.n_house, 1 + P.n_house), dtype=bool)
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 3013-3130)

```python
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
    birth_entry_grant = SD.birth_entry_grant
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    if Vd_stay is None:
        Vd_stay = Vd

    if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj, prj = tenure_logit_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, tenure_choice_kappa, Vd_stay, transaction_support
        )
        prj_full: np.ndarray | None = prj
    elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
        VH, tcj = tenure_choice_kernel(
            Vd, b_grid, ctx.heq, ctx.hcost, dp_choice, purchase_floor, SD.birth_dp, birth_entry_grant, Vd_stay, False, transaction_support
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
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 3348-3555)

```python
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
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 4280-4297)

```python
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
```

## code/model/intergen_eqscale_seq_optimized/solver.py (lines 4322-4335)

```python
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
```

## code/model/intergen_eqscale_seq_optimized/parameters.py (lines 885-925)

```python
        jr = int(getattr(P, "J_R", P.J))
        work_mean = float(np.mean(income_profile[: max(jr, 1)]))
        if work_mean > 0:
            income_profile = income_profile / work_mean
    return income_profile


def get_fecundity_by_age(P: SimpleNamespace) -> np.ndarray:
    """Per-period conception probability by age index j (length J).

    omega1 == 0 -> all ones (production behavior, including beyond the
    terminal age): this exact rule is the bitwise-nesting guarantee.
    """
    J = int(P.J)
    w1 = float(getattr(P, "fecundity_omega1", 0.0))
    if w1 == 0.0:
        return np.ones(J, dtype=float)
    w2 = float(getattr(P, "fecundity_omega2", 0.0))
    terminal = float(getattr(P, "fecundity_terminal_age", 45.0))
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)
    pi = 1.0 - w1 * np.exp(w2 * (ages - float(P.age_start)))
    pi = np.clip(pi, 0.0, 1.0)
    terminal_decay = float(getattr(P, "fecundity_terminal_decay", 0.0))
    if terminal_decay < 0.0 or not np.isfinite(terminal_decay):
        raise ValueError("fecundity_terminal_decay must be finite and nonnegative.")
    if terminal_decay > 0.0:
        tail_start = float(getattr(P, "fecundity_tail_start_age", 40.0))
        if not np.isfinite(tail_start):
            raise ValueError("fecundity_tail_start_age must be finite.")
        pi *= np.exp(-terminal_decay * np.maximum(ages - tail_start, 0.0))
    pi[ages >= terminal] = 0.0
    return pi


def fecundity_active(P: SimpleNamespace) -> bool:
    return float(getattr(P, "fecundity_omega1", 0.0)) != 0.0


def readiness_gate_active(P: SimpleNamespace) -> bool:
    """Whether the default-off E6c childless readiness state is active."""
    return bool(getattr(P, "readiness_gate_enabled", False))
```

## code/model/tools/run_e5f_perfect_foresight_transition.py (lines 358-386)

```python
def rents_from_asset_prices(
    prices: Sequence[float], terminal_price: float, P: SimpleNamespace
) -> np.ndarray:
    """Apply the one-period owner/renter no-arbitrage identity.

    Owners earn next period's asset price and pay depreciation and property
    tax.  Hence r_t + p_{t+1} = (R + delta + tau_H) p_t.  At a constant
    price this reduces exactly to the stationary user-cost identity used by
    the existing model.
    """
    current = np.asarray(prices, dtype=float).reshape(-1)
    if current.size < 1 or np.any(~np.isfinite(current)) or np.any(current <= 0.0):
        raise ValueError("Asset prices must be finite and strictly positive.")
    terminal = float(terminal_price)
    if not math.isfinite(terminal) or terminal <= 0.0:
        raise ValueError("The terminal asset price must be finite and positive.")
    next_prices = np.r_[current[1:], terminal]
    carrying_factor = float(P.R_gross) + float(P.delta) + float(P.tau_H)
    user_cost = float(getattr(P, "user_cost_rate", carrying_factor - 1.0))
    if not math.isclose(user_cost, carrying_factor - 1.0, rel_tol=0.0, abs_tol=2e-14):
        raise ValueError("Stationary user cost disagrees with interest, depreciation and tax")
    # The equivalent expression avoids subtracting two price-sized terms.
    # Constant paths reproduce the stationary Bellman's rent bit for bit.
    rents = user_cost * current + (current - next_prices)
    if np.any(~np.isfinite(rents)) or np.any(rents <= 0.0):
        raise ValueError(
            "The candidate asset-price path implies a nonpositive renter price."
        )
    return rents
```

## code/model/tools/run_e5f_open_population_transition.py (lines 767-805)

```python
def calendar_topcode_birth_accounting(
    g_pre: np.ndarray,
    g_post: np.ndarray,
    explicit_births: float,
    P: SimpleNamespace,
) -> dict[str, float]:
    """Translate the explicit 0/1/2/3+ state into measured child units.

    Entry into the last parity state identifies families reaching the 3+ bin.
    The additional children represented by that bin are used only for aggregate
    population renewal; household choices continue to use the existing state.
    """
    top_state = int(P.n_parity) - 1
    if top_state != 3 or str(getattr(P, "fertility_units", "")).lower() != "literal_topcode":
        return {
            "explicit_birth_children": float(explicit_births),
            "top_bin_entry_birth_flow": 0.0,
            "topcode_adjusted_birth_children": float(explicit_births),
        }
    top_weight = float(getattr(P, "tfr_top_bin_weight", top_state))
    top_before = float(np.sum(g_pre[:, :, :, :, :, top_state, :]))
    top_after = float(np.sum(g_post[:, :, :, :, :, top_state, :]))
    top_entry = top_after - top_before
    if top_entry < -2e-12:
        raise RuntimeError(
            f"Top-bin mass fell during the fertility stage: {top_entry:.3e}"
        )
    top_entry = max(top_entry, 0.0)
    adjusted = float(explicit_births) + (top_weight - top_state) * top_entry
    if adjusted + 1e-14 < float(explicit_births):
        raise RuntimeError("Top-code adjustment reduced the birth flow")
    return {
        "explicit_birth_children": float(explicit_births),
        "top_bin_entry_birth_flow": top_entry,
        "topcode_adjusted_birth_children": adjusted,
    }


def children_at_home_units(distribution: np.ndarray, P: SimpleNamespace) -> float:
```

# Transition closure evidence

Source: output/model/fixed_reference_transition_20260928/preparation_v1/transition_readiness.md. This is the preparation team's current closure assessment, not a completed transition result. The selected sections below are reproduced verbatim.

## Closures and their status

| Object | Classification | Contract for this preparation |
|---|---|---|
| Preferences | Estimated or calibration-normalized, then fixed | All saved preferences; no post-shock normalization. Historical preference shocks are a separate, pending specification. |
| Earnings | Externally estimated, fixed | Approved B15 15-state Markov process and saved age profile; use completed measurement audit, not constructor defaults. |
| Initial population | Empirically normalized reference | Exact saved pre-choice distribution, mass one. No age reweighting, reset or rescaling along a path. |
| Adult entry | Author-fixed reduced-form normalization | Half of each birth vintage enters after 16 years and half after 20; adjusted children map to households at 1/2.1. Preserve both raw and adjusted queues. This is not a new estimate of child survival. |
| Entry wealth/income | Empirically normalized; coupling approximation retained | Full saved conditional entrant distribution and grid. No zero-wealth fallback or frontier censoring of inherited population. |
| Survival | Externally fixed | Exact saved age survival and terminal death; no mortality change. |
| Geography and outside entry | National target scope; closed-path assumption explicit | One pooled market; no spatial migration in I=1. No new immigration, retention parameter, age bridge or old quota defaults. National calibration does not itself validate a population forecast. |
| PAYGO pension | Empirical baseline normalization; fixed tax closure | Hold saved payroll tax at 8.028%; solve each date's equal pension from actual worker/retiree masses. Baseline benefit 0.918 is a starting value, not a fixed transition benefit. |
| Property tax | Externally fixed | Saved annual rate 1.060%, period rate 4.239%, zero household rebate. Do not import old equal rebates. |
| Estates/entry funding | Outstanding substantive settlement; provisional reference ledger retained | Net positive estates fund actual next-period positive entrant assets; residual sink; funding shortage and negative estates fail. Donor utility stays unchanged. Lender counterparties, recipient law and physical settlement remain unresolved. |
| Housing, credit experiment | Estimated intercept, fixed external elasticity | Retain saved absolute supply curve (elasticity 0.630). Do not rescale it by population. |
| Housing, fixed-stock experiment | Author-fixed counterfactual | H equals actual reference supply at the reference equilibrium price, not H0. Prices/rents and individual housing/tenure choices may change. Constant gross stock entails replacement of depreciation; it is not zero gross construction. |
| Credit experiment | Author-adopted comparison; implementation outstanding | Replace artificial purchaser and incumbent debt restrictions by lifetime no-default solvency; keep interest, repayment, prices and death settlement. No arbitrary debt-floor substitute. |
| Terminal population | Endogenous equilibrium object; outstanding | Solve renewal, PAYGO and housing jointly at fixed preferences. A normalized stationary price root is not enough. |

The reference has a small measured renewal discrepancy: entry exceeds potential
birth-derived entry by \(4.889\times10^{-8}\) per reference household. Seed
prehistory from actual saved entry, then let actual births enter the queue.
Report resulting no-shock drift rather than adjusting child benefit or queues.

## Natural solvency and the new steady state

A natural borrowing limit is the most debt a household can repay under every
modeled future event with positive probability. Construct feasible sets backward
at each existing decision node, preserving the order of fertility realizations
and subsequent housing choices. Feasibility must be a separate Boolean object;
large negative utility is not itself proof of infeasibility.

For each chosen saving/tenure branch, require the resulting successor state to
be feasible for every reachable income and family-state outcome. Whenever death
has positive probability, require the inherited timing's net liquidation estate
\(b'+(1-\psi_{sell})q_t h\geq0\). The terminal condition is the same repayment
condition. With positive death risk and no default or life insurance, this
condition can itself exclude unsecured debt even when future earnings are
positive; removing an artificial limit does not authorize unpaid death-state
liabilities. The current model values this estate at today's price after saving;
changing its timing would be an additional economic change.

The current `native_solvency_credit` uses value cutoffs, first feasible grid
nodes and the grid's lower end; it is a prototype, not yet a certified natural
limit. It also conflicts with saved `native_due_stayer_credit=True`. Replacing
the incumbent-owner restriction is part of the authorized credit change, but
must be explicit. Setting financed share to one only removes a down payment;
it does not remove the collateral-based debt limit. Verify the new feasible
sets against terminal/worst-income cases and an expanded/refined debt grid,
without altering the reference entry distribution or relaxing occupied gates.

For a closed positive stationary population, let \(B(q,p)\), \(E(q,p)\) and
\(d(q,p)\) denote adjusted births, new household entry and housing demand per
normalized household, with pension \(p\). At fixed preferences the endpoint
must satisfy
\[
B(q,p)/(2.1E(q,p))=1,\qquad
\tau_{pay}Y_W(q,p)=p N_R(q,p),\qquad
N d(q,p)=H^S(q).
\]
Here \(Y_W\) is gross worker earnings and \(N_R\) is retiree exposure in the
normalized distribution. After solving renewal and PAYGO, housing determines
the level \(N=H^S(q)/d(q,p)\); for fixed stock replace \(H^S(q)\) by \(\bar H\).
One-step native distribution and queue reproduction must then pass. If no
positive root is found, report that failure and the search range; do not force
replacement fertility, rescale entrants or claim nonexistence from a timeout.

## Supply scaling: a candidate endpoint, not a transition

For a pure 10% increase in the housing-supply intercept, a candidate endpoint
keeps prices, pension, household policies and the normalized distribution
unchanged and multiplies population, both entry queues, births, deaths, housing
demand and fiscal flows by 1.1. This follows from the level-linear forward
population operator and the per-household earnings and entry laws; both sides
of PAYGO and the provisional estate-funding ledger scale by the same factor.
The absolute supply schedule also scales by 1.1. The reference's small renewal
residual scales in levels and is unchanged in relative terms.

This algebra applies to the present single-market, exogenous-earnings closure
with no outside entry and no fixed aggregate transfer. A fixed immigration
flow, aggregate fiscal grant, population-dependent wages or amenities, a
non-scaling estate allocation, or normalization of population at each date
would break it. The reviewed estate ledger is homogeneous in its economic
flows; its absolute numerical tolerance is not an economic transfer. Native
one-step and scaling checks are still required, and no uniqueness, stability
or transition claim follows. The economic-analysis chat owns this supply
comparison; no additional supply solve was launched here.

## What each calculation can establish
