@njit(cache=False, parallel=True)
def full_renter_block_kernel(
    Rv1d,           # (Nb,)
    Rvt1d,          # (Nb,) debt-blind means-test resources
    Vc_flat,        # (Nb, nc)
    bp_prev,        # (Nb, nc) or zeros (use has_prev to gate)
    has_prev,       # bool/int
    b_grid,         # (Nb,)
    cb_v,           # (nc,)
    hb_v,           # (nc,)
    psi_v,          # (nc,)
    gb_v,           # (nc,)
    alpha_v,        # (nc,)
    esc_v,          # (nc,)
    ri,
    hR_max,
    c_min,
    c_bar_0,
    h_bar_0,
    alpha,
    oms,
    beta,
    s_next,
    D_next,
    gs_alpha1,
    gs_alpha2,
    gs_tol,
    exhaustive_saving=0,
):
    # Full-Bellman renter block: golden-section search for bp + post-search
    # consumption / housing arithmetic, fused into one kernel per (i, j).
    # When bp_prev is given (j < J-1), the search interval is clamped to
    # [bp_prev - 2, bp_prev + 2] as a soft monotonicity prior — a
    # heuristic that mirrors the MATLAB implementation, not a strict
    # invariant of the model.
    if exhaustive_saving and has_prev:
        raise ValueError("Exhaustive saving requires the full feasible interval")
    Nb, nc = Vc_flat.shape
    Vo = np.empty((Nb, nc))
    bp_out = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    ho = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    inv_ri = 1.0 / ri
    bg0 = b_grid[0]
    for c in prange(nc):
        cbc = cb_v[c]
        hbc = hb_v[c]
        psic = psi_v[c]
        gc = gb_v[c]
        al = alpha_v[c]
        es = esc_v[c]
        # Keep the nested Stone--Geary arithmetic byte-for-byte in its own
        # branch; the eqscale branch is the only one that substitutes al/es.
        if es != 1.0 or al != alpha:
            Kr = (al ** al * ((1.0 - al) * inv_ri) ** (1.0 - al)) ** oms
            cap_c = ri * (hR_max - hbc) / (1.0 - al)
            ht_cap_pow = max(hR_max - hbc, 1e-10) ** (1.0 - al)
        else:
            Kr = (alpha ** alpha * ((1.0 - alpha) * inv_ri) ** (1.0 - alpha)) ** oms
            cap_c = ri * (hR_max - hbc) / (1.0 - alpha)
            ht_cap_pow = max(hR_max - hbc, 1e-10) ** (1.0 - alpha)
        dc = cbc + ri * hbc
        ht_cap_c = hR_max - hbc
        if ht_cap_c < 1e-10:
            ht_cap_c = 1e-10
        for b in range(Nb):
            Rvb = Rv1d[b]
            if gc > 0.0:
                Tb = gc - Rvt1d[b]
                if Tb > 0.0:
                    if Tb > gc:
                        Tb = gc
                    Rvb = Rvb + Tb
            current_b = b_grid[b]
            rollover_floor = s_next * (current_b if current_b < 0.0 else 0.0)
            line_floor = -D_next
            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor
            lo = unsecured_floor
            if bg0 > lo:
                lo = bg0
            hi = Rvb - dc - 1e-6
            if hi < lo:
                hi = lo
            if has_prev:
                lo_prev = bp_prev[b, c] - 2.0
                if lo_prev > lo:
                    lo = lo_prev
                hi_prev = bp_prev[b, c] + 2.0
                if hi_prev < hi:
                    hi = hi_prev
                if lo < unsecured_floor:
                    lo = unsecured_floor
                if lo < bg0:
                    lo = bg0
                if hi < lo:
                    hi = lo

            if exhaustive_saving:
                bp_best, v_best = exhaustive_saving_scalar(
                    lo, hi, Rvb, Vc_flat[:, c], b_grid, ri, hbc, cbc, psic,
                    hR_max, al, oms, beta, es, 0.0, 0.0, False)
            else:
                d = hi - lo
                x1 = lo + gs_alpha1 * d
                x2 = lo + gs_alpha2 * d
                f1 = eval_renter_scalar(x1, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es)
                f2 = eval_renter_scalar(x2, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es)
                d = gs_alpha1 * gs_alpha2 * d
                while d > gs_tol:
                    if f2 >= f1:
                        xe = x2 + d
                        if xe > hi:
                            xe = hi
                        fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es)
                        x1 = x2
                        f1 = f2
                        x2 = xe
                        f2 = fe
                    else:
                        xe = x1 - d
                        if xe < lo:
                            xe = lo
                        fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es)
                        x2 = x1
                        f2 = f1
                        x1 = xe
                        f1 = fe
                    d = d * gs_alpha2
                if f2 >= f1:
                    bp_best = x2
                    v_best = f2
                else:
                    bp_best = x1
                    v_best = f1
            bp_out[b, c] = bp_best
            Vo[b, c] = v_best

            surplus = Rvb - dc - bp_best
            if surplus <= 1e-10:
                co[b, c] = c_bar_0 + c_min
                ho[b, c] = h_bar_0 + 0.01
            else:
                ht_unc = (1.0 - al) * surplus * inv_ri
                if hbc + ht_unc > hR_max:
                    ct = Rvb - cbc - ri * hR_max - bp_best
                    if ct < 1e-10:
                        ct = 1e-10
                    ct_eff = ct if ct > c_min else c_min
                    ht_eff = ht_cap_c if ht_cap_c > 0.01 else 0.01
                    if exhaustive_saving:
                        # Report the intratemporal allocation used by the
                        # objective, without the legacy output-only floors.
                        ct_eff = ct
                        ht_eff = ht_cap_c
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
                else:
                    ct = al * surplus
                    ct_eff = ct if ct > c_min else c_min
                    ht_eff = ht_unc if ht_unc > 0.01 else 0.01
                    if exhaustive_saving:
                        ct_eff = ct
                        ht_eff = ht_unc
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
    return Vo, bp_out, co, ho

@njit(cache=False, parallel=True)
def full_owner_block_kernel(
    Rv1d,           # (Nb,)
    Rvt1d,          # (Nb,) debt-blind means-test resources
    Vco_flat,       # (Nb, nc)
    bp_prev,        # (Nb, nc) or zeros
    has_prev,
    b_grid,
    cb_v,           # (nc,)
    hb_v,           # (nc,)
    psi_v,          # (nc,)
    gb_v,           # (nc,)
    alpha_v,        # (nc,)
    esc_v,          # (nc,)
    bf_v,           # (nc,) — bmo[i, ten, nn, cs] flattened in F order
    oc,
    hsv,
    owner_h_bar_scale,
    owner_service_premium,
    c_min,
    alpha,
    oms,
    beta,
    s_next,
    D_next,
    gs_alpha1,
    gs_alpha2,
    gs_tol,
    strict_hbar_feasibility=0,
    exhaustive_saving=0,
):
    if exhaustive_saving and has_prev:
        raise ValueError("Exhaustive saving requires the full feasible interval")
    Nb, nc = Vco_flat.shape
    Vo = np.empty((Nb, nc))
    bp_out = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    bg0 = b_grid[0]
    for c in prange(nc):
        cbc = cb_v[c]
        hbc = hb_v[c]
        psic = psi_v[c]
        gc = gb_v[c]
        bf = bf_v[c]
        al = alpha_v[c]
        es = esc_v[c]
        ht_c = hsv - owner_h_bar_scale * hbc
        if strict_hbar_feasibility and ht_c <= 0.0:
            for b in range(Nb):
                Vo[b, c] = -1e10
                bp_out[b, c] = b_grid[0]
                co[b, c] = cbc + c_min
            continue
        if ht_c < 1e-10:
            ht_c = 1e-10
        ht_c = owner_service_premium * ht_c
        if es != 1.0 or al != alpha:
            Ko_c = ht_c ** ((1.0 - al) * oms)
        else:
            Ko_c = ht_c ** ((1.0 - alpha) * oms)
        for b in range(Nb):
            Rvb = Rv1d[b]
            if gc > 0.0:
                Tb = gc - Rvt1d[b]
                if Tb > 0.0:
                    if Tb > gc:
                        Tb = gc
                    Rvb = Rvb + Tb
            # b already contains secured mortgage debt.  Subtract the current
            # collateral floor before rolling only the unsecured component.
            current_unsecured = b_grid[b] - bf
            rollover_floor = s_next * (current_unsecured if current_unsecured < 0.0 else 0.0)
            line_floor = -D_next
            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor
            total_floor = bf + unsecured_floor
            lo = total_floor
            if bg0 > lo:
                lo = bg0
            hi = Rvb - oc - cbc - 1e-6
            if hi < lo:
                hi = lo
            if has_prev:
                lo_prev = bp_prev[b, c] - 2.0
                if lo_prev > lo:
                    lo = lo_prev
                hi_prev = bp_prev[b, c] + 2.0
                if hi_prev < hi:
                    hi = hi_prev
                if lo < total_floor:
                    lo = total_floor
                if lo < bg0:
                    lo = bg0
                if hi < lo:
                    hi = lo

            if exhaustive_saving:
                bp_best, v_best = exhaustive_saving_scalar(
                    lo, hi, Rvb, Vco_flat[:, c], b_grid, 1.0, 0.0, cbc, psic,
                    1.0, al, oms, beta, es, oc, Ko_c, True)
            else:
                d = hi - lo
                x1 = lo + gs_alpha1 * d
                x2 = lo + gs_alpha2 * d
                f1 = eval_owner_scalar(x1, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, al, oms, beta, es)
                f2 = eval_owner_scalar(x2, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, al, oms, beta, es)
                d = gs_alpha1 * gs_alpha2 * d
                while d > gs_tol:
                    if f2 >= f1:
                        xe = x2 + d
                        if xe > hi:
                            xe = hi
                        fe = eval_owner_scalar(xe, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, al, oms, beta, es)
                        x1 = x2
                        f1 = f2
                        x2 = xe
                        f2 = fe
                    else:
                        xe = x1 - d
                        if xe < lo:
                            xe = lo
                        fe = eval_owner_scalar(xe, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, al, oms, beta, es)
                        x2 = x1
                        f2 = f1
                        x1 = xe
                        f1 = fe
                    d = d * gs_alpha2
                if f2 >= f1:
                    bp_best = x2
                    v_best = f2
                else:
                    bp_best = x1
                    v_best = f1
            bp_out[b, c] = bp_best
            Vo[b, c] = v_best

            ct = Rvb - oc - cbc - bp_best
            ct_eff = ct if ct > c_min else c_min
            if exhaustive_saving and ct > 1e-10:
                ct_eff = ct
            co[b, c] = cbc + ct_eff
    return Vo, bp_out, co
