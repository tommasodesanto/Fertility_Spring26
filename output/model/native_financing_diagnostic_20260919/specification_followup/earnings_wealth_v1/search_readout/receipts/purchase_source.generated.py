def solve_bellman_full_markov_income(
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    continuation_V: np.ndarray | None = None,
):
    t0 = time.perf_counter()
    if (bool(getattr(P, "joint_nested_choice", False))
            or bool(getattr(P, "use_pti_constraint", False))
            or float(P.lambda_d) != 0.0
            or not bool(getattr(P, "use_tenure_kernel", True))
            or not bool(getattr(P, "use_full_kernel", True))
            or not NUMBA_AVAILABLE
            or str(getattr(P, "interp_method", "linear")) != "linear"
            or np.any(SD.birth_dp) or np.any(SD.birth_entry_grant)):
        raise ValueError("unsupported mechanism combined with purchase-income diagnostic")
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

    phi_choice = SD.phi_choice
    birth_entry_grant = SD.birth_entry_grant
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

    Vbq = np.zeros((Nb, nt, I, npar, ncs))
    for i in range(I):
        for ten in range(nt):
            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)

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
            Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs, age_index=j)
            Vd = np.zeros((Nb, nt, I, npar, ncs))
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
                hRmax = P.hR_max
                Vcr = flat_nc(Vc[:, 0, i, :, :], Nb, nc)
                Rv1d_full = np.ascontiguousarray(Rv[:, 0])
                Rvt1d_full = np.ascontiguousarray(Rv_test[:, 0])
                if use_full_kernel:
                    bp_prev_r = np.zeros((Nb, nc))
                    has_prev_r = 0
                    Vo_nc, bp_nc, co_nc, ho_nc = full_renter_block_kernel(
                        Rv1d_full, Rvt1d_full, Vcr, bp_prev_r, has_prev_r, b_grid,
                        cb_v, hb_v, psi_v_flat, gb_v, alpha_v, esc_v,
                        ri, hRmax, P.c_min, P.c_bar_0, P.h_bar_0,
                        alpha, oms, beta, s_next, D_next, gs_alpha1, gs_alpha2, gs_tol,
                        int(exhaustive_saving),
                    )
                else:
                    Kr = (alpha**alpha * ((1 - alpha) / ri) ** (1 - alpha)) ** oms
                    d_nc = SD.cb_flat + ri * SD.hb_flat
                    Rv_eff_nc = Rv + np.clip(SD.gb_flat - Rv_test, 0.0, SD.gb_flat)
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
                        hi = np.maximum(Rv_eff_nc[:, c] - dc - 1e-6, lo)
                        bp, val = golden_renter(
                            lo, hi, Rv_eff_nc[:, c], Vbar, b_grid, dc, pc, cc, cb_c, hb_c,
                            ri, hRmax, ht_cap_c, Kr, alpha, oms, beta,
                            gs_alpha1, gs_alpha2, gs_tol, interp_method, 1.0,
                        )
                        bp_nc[:, c] = bp
                        Vo_nc[:, c] = val
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

                for ten in range(1, nt):
                    oc = ocst[i, ten]
                    hsv = hsrv[i, ten]
                    Vco = flat_nc(Vc[:, ten, i, :, :], Nb, nc)
                    if use_full_kernel:
                        bf_v_base = bmo[i, ten, :, :].reshape(-1, order="F")
                        bf_v = np.ascontiguousarray(
                            effective_owner_collateral_floor(P, bf_v_base, j)
                        )
                        bp_prev_o = np.zeros((Nb, nc))
                        has_prev_o = 0
                        Vo_nc, bp_nc, co_nc = full_owner_block_kernel(
                            Rv1d_full, Rvt1d_full, Vco, bp_prev_o, has_prev_o, b_grid,
                            cb_v, hb_v, psi_v_flat, gb_v, alpha_v, esc_v, bf_v,
                            oc, hsv, owner_h_bar_scale, owner_service_premium, P.c_min,
                            alpha, oms, beta, 0.0, 0.0, gs_alpha1, gs_alpha2, gs_tol,
                            strict_owner_hbar_feasibility, int(exhaustive_saving),
                        )
                    else:
                        for c in range(nc):
                            Vbar = Vco[:, c]
                            cb_c = SD.cb_flat[0, c]
                            pc = SD.psi_flat[0, c]
                            owner_residual_h = hsv - owner_h_bar_scale * SD.hb_flat[0, c]
                            nn_c, cs_c = decode_flat_family_state(c, npar)
                            bf_c = bmo[i, ten, nn_c, cs_c]
                            lo = np.maximum(owner_borrowing_floor(P, b_grid, bf_c, j), b_grid[0])
                            hi = np.maximum(Rv_eff_nc[:, c] - oc - cb_c - 1e-6, lo)
                            if strict_owner_hbar_feasibility and owner_residual_h <= 0.0:
                                bp = lo.copy()
                                val = np.full(Nb, -1e10)
                            else:
                                Ko_c = (
                                    owner_service_premium * max(owner_residual_h, 1e-10)
                                ) ** ((1 - alpha) * oms)
                                bp, val = golden_owner(
                                    lo, hi, Rv_eff_nc[:, c], Vbar, b_grid, oc, cb_c, pc,
                                    Ko_c, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol, interp_method, 1.0,
                                )
                            bp_nc[:, c] = bp
                            Vo_nc[:, c] = val
                        co_nc = SD.cb_flat + np.maximum(Rv_eff_nc - oc - SD.cb_flat - bp_nc, P.c_min)
                    if exhaustive_saving:
                        resources = Rv + np.clip(SD.gb_flat - Rv_test, 0.0, SD.gb_flat)
                        co_nc = joint_nested.owner_consumption_from_solution(
                            resources, oc, bp_nc, SD.cb_flat, Vo_nc, co_nc)
                    Vd[:, ten, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                    bd[:, ten, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
                    cd[:, ten, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)

            c_pol[:, :, :, j, zz, :, :] = cd
            hR_pol[:, :, :, j, zz, :, :] = hd
            bp_pol[:, :, :, j, zz, :, :] = bd

            income_for_purchase = np.array([
                income_at_state(P, i, j, float(z_value)) for i in range(I)
            ], dtype=float).reshape(I, 1, 1, 1) / Rg
            dp_choice = dp_arr - income_for_purchase
            # Require actual post-transaction wealth to lie on the solved grid.
            bmo_purchase = np.maximum(bmo - income_for_purchase, b_grid[0])
            if bool(getattr(P, "use_pti_constraint", False)):
                income_j = np.array([income_at_state(P, i, j, float(z_value)) for i in range(I)], dtype=float)
                dp_choice = pti_adjusted_downpayment(dp_arr, hcost, income_j, P, b_grid)

            if joint_active:
                joint_result = joint_nested.bellman_block(
                    Vd, (b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant),
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

            if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
                VH, tcj, prj = tenure_logit_kernel(
                    Vd, b_grid, heq, hcost, dp_choice, bmo_purchase, SD.birth_dp, birth_entry_grant, tenure_choice_kappa
                )
                tenure_probs[:, :, :, j, zz, :, :, :] = prj
            elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
                VH, tcj = tenure_choice_kernel(
                    Vd, b_grid, heq, hcost, dp_choice, bmo_purchase, SD.birth_dp, birth_entry_grant
                )
            else:
                VH = np.zeros((Nb, nt, I, npar, ncs))
                tcj = np.zeros((Nb, nt, I, npar, ncs), dtype=np.int16)
                for id_ in range(I):
                    for to in range(nt):
                        sp = heq[id_, to] if to > 0 else 0.0
                        Vopt = np.zeros((Nb, npar, ncs, nt))
                        if to == 0:
                            Vopt[:, :, :, 0] = Vd[:, 0, id_, :, :]
                        else:
                            ba = np.clip(b_grid + sp, b_grid[0], b_grid[-1])
                            Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            Vow = Vd[:, tn, id_, :, :]
                            if to == tn:
                                Vopt[:, :, :, tn] = Vow
                            elif to == 0:
                                bab = b_grid - hc
                                Vb = interp_on_grid(b_grid, Vow, bab)
                                for nn in range(npar):
                                    for cs in range(ncs):
                                        dpn = dp_choice[id_, tn, nn, cs]
                                        bmn = bmo[id_, tn, nn, cs]
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
                                        bmn = bmo[id_, tn, nn, cs]
                                        dpc = dpn - sp
                                        inf_m = (b_grid < dpc) | (bar < bmn)
                                        Vrs[inf_m, nn, cs] = -1e10
                                Vopt[:, :, :, tn] = Vrs
                        tc = np.argmax(Vopt, axis=3)
                        if use_tenure_logit:
                            ls, pr = logsumexp(Vopt / tenure_choice_kappa, axis=3)
                            pr[np.max(Vopt, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0
                            VH[:, to, id_, :, :] = tenure_choice_kappa * ls
                            tenure_probs[:, to, id_, j, zz, :, :, :] = pr.astype(np.float32)
                        else:
                            VH[:, to, id_, :, :] = np.max(Vopt, axis=3)
                        tcj[:, to, id_, :, :] = tc
            tenure_choice[:, :, :, j, zz, :, :] = tcj

            kl = P.kappa_loc
            if NUMBA_AVAILABLE and bool(getattr(P, "use_loc_kernel", True)):
                VI, lpj = location_logit_kernel(VH, iidx, iwt, loc_shift, kl)
            else:
                VI = np.zeros((Nb, nt, I, npar, ncs))
                lpj = np.zeros((Nb, nt, I, I, npar, ncs))
                for io in range(I):
                    for to in range(nt):
                        Va = np.zeros((Nb, I, npar, ncs))
                        Va[:, io, :, :] = VH[:, to, io, :, :]
                        idx = iidx[:, io, to]
                        wt = iwt[:, io, to]
                        for id_ in range(I):
                            if id_ == io:
                                continue
                            Vdst = VH[:, 0, id_, :, :]
                            Va[:, id_, :, :] = (1 - wt)[:, None, None] * Vdst[idx, :, :] + wt[:, None, None] * Vdst[idx + 1, :, :]
                        la = Va.copy()
                        for id_ in range(I):
                            la[:, id_, :, :] += loc_shift[io, id_]
                        la = la / kl
                        ls, pr = logsumexp(la, axis=1)
                        dead_loc = np.max(Va, axis=1) <= DEAD_VALUE_CUTOFF
                        pr = np.where(dead_loc[:, None, :, :], 0.0, pr)
                        VI[:, to, io, :, :] = kl * ls
                        lpj[:, to, io, :, :, :] = pr
            loc_probs[:, :, :, :, j, zz, :, :] = lpj

            if in_fert:
                pi_j = float(fec[j])
                if bool(getattr(P, "sequential_births", False)):
                    Vfa = np.empty((Nb, nt, I, 2))
                    settled_cs = readiness_settled_state(P)
                    Vfa[:, :, :, 0] = VI[:, :, :, 0, settled_cs]
                    Vfa[:, :, :, 1] = (
                        pi_j
                        * (
                            VI[:, :, :, 1, 1]
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
                            V2[:, :, :, 1] = (
                                pi_j * VI[:, :, :, nn + 1, destination_cs]
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
