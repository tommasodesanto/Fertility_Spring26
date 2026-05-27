"""Compiled Numba kernels for the DT Python port.

Block kernels at the bottom (eval/full renter & owner, tenure choice,
location logit, forward distribution) replace the per-cell Python+NumPy
inner loops in `solver.py`. Most use `parallel=True` with `prange(nc)`;
the parallelism is race-free because each column `c` writes to disjoint
output cells.
"""

from __future__ import annotations

import numpy as np

try:  # pragma: no cover - availability depends on local environment
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def deco(fn):
            return fn

        return deco

    prange = range  # type: ignore


@njit(cache=True)
def interp_scalar(bg, V, x):
    if x <= bg[0]:
        idx = 0
    elif x >= bg[bg.size - 1]:
        idx = bg.size - 2
    else:
        lo = 0
        hi = bg.size - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if bg[mid] <= x:
                lo = mid
            else:
                hi = mid
        idx = lo
    wt = (x - bg[idx]) / (bg[idx + 1] - bg[idx])
    if wt < 0.0:
        wt = 0.0
    elif wt > 1.0:
        wt = 1.0
    return (1.0 - wt) * V[idx] + wt * V[idx + 1]


@njit(cache=True)
def eval_renter_scalar(bp, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta):
    surplus = Rv - dc - bp
    if surplus <= 1e-10:
        return -1e10
    if surplus > cc:
        ct = Rv - cb_c - ri * hRmax - bp
        if ct < 1e-10:
            ct = 1e-10
        u = (ct**alpha * ht_cap_c ** (1.0 - alpha)) ** oms / oms + pc
    else:
        ss = surplus
        if ss < 1e-10:
            ss = 1e-10
        u = Kr * ss**oms / oms + pc
    return u + beta * interp_scalar(bg, Vbar, bp)


@njit(cache=True)
def golden_renter_kernel(lo, hi, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, a1, a2, tol):
    n = Rv.size
    bp_out = np.empty(n)
    val_out = np.empty(n)
    for k in range(n):
        lok = lo[k]
        hik = hi[k]
        d = hik - lok
        x1 = lok + a1 * d
        x2 = lok + a2 * d
        f1 = eval_renter_scalar(x1, Rv[k], Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
        f2 = eval_renter_scalar(x2, Rv[k], Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
        d = a1 * a2 * d
        while d > tol:
            if f2 >= f1:
                xe = x2 + d
                fe = eval_renter_scalar(xe, Rv[k], Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
                x1 = x2
                f1 = f2
                x2 = xe
                f2 = fe
            else:
                xe = x1 - d
                fe = eval_renter_scalar(xe, Rv[k], Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
                x2 = x1
                f2 = f1
                x1 = xe
                f1 = fe
            d = d * a2
        if f2 >= f1:
            bp_out[k] = x2
            val_out[k] = f2
        else:
            bp_out[k] = x1
            val_out[k] = f1
    return bp_out, val_out


@njit(cache=True)
def eval_owner_scalar(bp, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta):
    ct = Rv - oc - cb_c - bp
    if ct <= 1e-10:
        return -1e10
    return Ko_c * ct ** (alpha * oms) / oms + pc + beta * interp_scalar(bg, Vbar, bp)


@njit(cache=True)
def golden_owner_kernel(lo, hi, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, a1, a2, tol):
    n = Rv.size
    bp_out = np.empty(n)
    val_out = np.empty(n)
    for k in range(n):
        lok = lo[k]
        hik = hi[k]
        d = hik - lok
        x1 = lok + a1 * d
        x2 = lok + a2 * d
        f1 = eval_owner_scalar(x1, Rv[k], Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta)
        f2 = eval_owner_scalar(x2, Rv[k], Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta)
        d = a1 * a2 * d
        while d > tol:
            if f2 >= f1:
                xe = x2 + d
                fe = eval_owner_scalar(xe, Rv[k], Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta)
                x1 = x2
                f1 = f2
                x2 = xe
                f2 = fe
            else:
                xe = x1 - d
                fe = eval_owner_scalar(xe, Rv[k], Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta)
                x2 = x1
                f2 = f1
                x1 = xe
                f1 = fe
            d = d * a2
        if f2 >= f1:
            bp_out[k] = x2
            val_out[k] = f2
        else:
            bp_out[k] = x1
            val_out[k] = f1
    return bp_out, val_out


@njit(cache=True)
def scatter_vec_kernel(idx, wt, mass, Nb):
    out = np.zeros(Nb)
    for r in range(mass.size):
        m = mass[r]
        if m == 0.0:
            continue
        k = idx[r]
        w = wt[r]
        out[k] += (1.0 - w) * m
        out[k + 1] += w * m
    return out


@njit(cache=True)
def scatter_cols_kernel(idx, wt, mass, Nb):
    nrow, ncol = mass.shape
    out = np.zeros((Nb, ncol))
    for c in range(ncol):
        for r in range(nrow):
            m = mass[r, c]
            if m == 0.0:
                continue
            k = idx[r, c]
            w = wt[r, c]
            out[k, c] += (1.0 - w) * m
            out[k + 1, c] += w * m
    return out


@njit(cache=True)
def scatter_cols_sameidx_kernel(idx, wt, mass, Nb):
    nrow, ncol = mass.shape
    out = np.zeros((Nb, ncol))
    for c in range(ncol):
        for r in range(nrow):
            m = mass[r, c]
            if m == 0.0:
                continue
            k = idx[r]
            w = wt[r]
            out[k, c] += (1.0 - w) * m
            out[k + 1, c] += w * m
    return out


@njit(cache=True)
def interp_cols_kernel(bg, V, bq):
    nrow, ncol = bq.shape
    out = np.empty((nrow, ncol))
    nb = bg.size
    for c in range(ncol):
        for r in range(nrow):
            x = bq[r, c]
            if x <= bg[0]:
                idx = 0
            elif x >= bg[nb - 1]:
                idx = nb - 2
            else:
                lo = 0
                hi = nb - 1
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if bg[mid] <= x:
                        lo = mid
                    else:
                        hi = mid
                idx = lo
            wt = (x - bg[idx]) / (bg[idx + 1] - bg[idx])
            if wt < 0.0:
                wt = 0.0
            elif wt > 1.0:
                wt = 1.0
            out[r, c] = (1.0 - wt) * V[idx, c] + wt * V[idx + 1, c]
    return out


@njit(cache=True)
def interp_cols_preidx_kernel(V, idx, wt):
    nrow, ncol = idx.shape
    out = np.empty((nrow, ncol))
    for c in range(ncol):
        for r in range(nrow):
            k = idx[r, c]
            w = wt[r, c]
            out[r, c] = (1.0 - w) * V[k, c] + w * V[k + 1, c]
    return out


@njit(cache=True)
def _advance_mass_component(
    g,
    mass,
    b_loc,
    ten_loc,
    loc,
    j,
    nn,
    cs,
    tenure_choice,
    bp_idx,
    bp_wt,
    tmx_idx,
    tmx_wt,
    Pi_child,
    K,
    ncs,
    use_stochastic_aging,
):
    # Advance one population cell from j -> j+1: tenure choice -> savings
    # interpolation -> child-state aging. `mass` is already conditioned on
    # the move/stay outcome by the caller. Returns mature entrants implied
    # by aging out of the youngest-child stage.
    if mass == 0.0:
        return 0.0

    tn = tenure_choice[b_loc, ten_loc, loc, j, nn, cs]
    kt = tmx_idx[loc, ten_loc, tn, nn, cs, b_loc]
    wt = tmx_wt[loc, ten_loc, tn, nn, cs, b_loc]
    entrant_units = 0.0

    for ot in range(2):
        if ot == 0:
            b_t = kt
            mt = (1.0 - wt) * mass
        else:
            b_t = kt + 1
            mt = wt * mass
        if mt == 0.0:
            continue

        ks = bp_idx[b_t, tn, loc, j, nn, cs]
        ws = bp_wt[b_t, tn, loc, j, nn, cs]
        for os in range(2):
            if os == 0:
                b_s = ks
                ms = (1.0 - ws) * mt
            else:
                b_s = ks + 1
                ms = ws * mt
            if ms == 0.0:
                continue

            if use_stochastic_aging:
                if cs == K and nn >= 1:
                    if nn == 1:
                        entrant_units += nn * Pi_child[cs, K + 1, nn] * ms
                    else:
                        entrant_units += nn * Pi_child[cs, K + 2, nn] * ms
                for csn in range(ncs):
                    pa = Pi_child[cs, csn, nn]
                    if pa != 0.0:
                        g[b_s, tn, loc, j + 1, nn, csn] += pa * ms
            else:
                if cs == 0:
                    csn = 0
                elif cs >= K + 1:
                    csn = cs
                elif cs < K:
                    csn = cs + 1
                else:
                    if nn == 0:
                        csn = 0
                    elif nn == 1:
                        csn = K + 1
                    else:
                        csn = K + 2
                if cs == K and csn >= K + 1 and nn >= 1:
                    entrant_units += nn * ms
                g[b_s, tn, loc, j + 1, nn, csn] += ms

    return entrant_units


@njit(cache=True)
def forward_distribution_fast_kernel(
    entry_by_loc,
    fert_probs,
    loc_probs,
    tenure_choice,
    bp_idx,
    bp_wt,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    Pi_child,
    Nb,
    nt,
    I,
    J,
    npar,
    ncs,
    entry_idx,
    A_f_start,
    A_f_end,
    K,
    use_stochastic_aging,
):
    # Fast-stat forward pass used during equilibrium iteration. Tracks
    # only price-relevant aggregates (population, births, mature entrants);
    # the full event-study moments are recomputed in Python on the final
    # accepted equilibrium, not here.
    g = np.zeros((Nb, nt, I, J, npar, ncs))
    births_by_loc = np.zeros(I)
    entrants_mature_by_loc = np.zeros(I)
    total_births = 0.0
    entrants_mature_total = 0.0

    for i in range(I):
        g[entry_idx, 0, i, 0, 0, 0] = entry_by_loc[i]

    for j in range(J - 1):
        age_idx = j + 1
        if age_idx >= A_f_start and age_idx <= A_f_end:
            for b in range(Nb):
                for ten in range(nt):
                    for loc in range(I):
                        m = g[b, ten, loc, j, 0, 0]
                        if m == 0.0:
                            continue
                        p0 = fert_probs[b, ten, loc, j, 0]
                        g[b, ten, loc, j, 0, 0] = m * p0
                        eb = 0.0
                        for nn in range(1, npar):
                            pn = fert_probs[b, ten, loc, j, nn]
                            mn = m * pn
                            if mn != 0.0:
                                g[b, ten, loc, j, nn, 1] += mn
                            eb += nn * pn
                        total_births += m * eb
                        births_by_loc[loc] += m * eb

        for b in range(Nb):
            for ten in range(nt):
                for io in range(I):
                    for nn in range(npar):
                        for cs in range(ncs):
                            m = g[b, ten, io, j, nn, cs]
                            if m == 0.0:
                                continue
                            for loc in range(I):
                                pl = loc_probs[b, ten, io, loc, j, nn, cs]
                                if pl == 0.0:
                                    continue
                                ml = m * pl
                                if loc == io:
                                    ent = _advance_mass_component(
                                        g,
                                        ml,
                                        b,
                                        ten,
                                        loc,
                                        j,
                                        nn,
                                        cs,
                                        tenure_choice,
                                        bp_idx,
                                        bp_wt,
                                        tmx_idx,
                                        tmx_wt,
                                        Pi_child,
                                        K,
                                        ncs,
                                        use_stochastic_aging,
                                    )
                                    if ent != 0.0:
                                        entrants_mature_by_loc[loc] += ent
                                        entrants_mature_total += ent
                                else:
                                    kl = lmm_idx[io, ten, b]
                                    wl = lmm_wt[io, ten, b]
                                    ent = _advance_mass_component(
                                        g,
                                        (1.0 - wl) * ml,
                                        kl,
                                        0,
                                        loc,
                                        j,
                                        nn,
                                        cs,
                                        tenure_choice,
                                        bp_idx,
                                        bp_wt,
                                        tmx_idx,
                                        tmx_wt,
                                        Pi_child,
                                        K,
                                        ncs,
                                        use_stochastic_aging,
                                    )
                                    if ent != 0.0:
                                        entrants_mature_by_loc[loc] += ent
                                        entrants_mature_total += ent
                                    ent = _advance_mass_component(
                                        g,
                                        wl * ml,
                                        kl + 1,
                                        0,
                                        loc,
                                        j,
                                        nn,
                                        cs,
                                        tenure_choice,
                                        bp_idx,
                                        bp_wt,
                                        tmx_idx,
                                        tmx_wt,
                                        Pi_child,
                                        K,
                                        ncs,
                                        use_stochastic_aging,
                                    )
                                    if ent != 0.0:
                                        entrants_mature_by_loc[loc] += ent
                                        entrants_mature_total += ent

    return g, total_births, births_by_loc, entrants_mature_by_loc, entrants_mature_total


@njit(cache=True, parallel=True)
def eval_renter_block_kernel(
    Rv1d,
    bpv,
    V,
    idx,
    wt,
    cb,
    hb,
    psi_v,
    ri,
    hR_max,
    c_min,
    c_bar_0,
    h_bar_0,
    alpha,
    oms,
    beta,
):
    # Howard-eval renter block for one (i, j): compute V, c, h policies at
    # the stored savings policy `bpv`. No optimization, just plug-in
    # evaluation. The `surplus > cap_c` branch is the renter housing-cap
    # corner: when desired h would exceed hR_max, switch to constrained
    # consumption with h pinned at hR_max - hb_c.
    Nb, nc = bpv.shape
    Vo = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    ho = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    inv_one_minus_alpha = 1.0 / (1.0 - alpha)
    inv_ri = 1.0 / ri
    Kr = (alpha ** alpha * ((1.0 - alpha) * inv_ri) ** (1.0 - alpha)) ** oms
    for c in prange(nc):
        cbc = cb[c]
        hbc = hb[c]
        psic = psi_v[c]
        dc = cbc + ri * hbc
        cap_c = ri * (hR_max - hbc) * inv_one_minus_alpha
        ht_cap_c = hR_max - hbc
        if ht_cap_c < 1e-10:
            ht_cap_c = 1e-10
        ht_cap_c_pow = ht_cap_c ** (1.0 - alpha)
        for b in range(Nb):
            bp = bpv[b, c]
            Rvb = Rv1d[b]
            surplus = Rvb - dc - bp
            k = idx[b, c]
            w = wt[b, c]
            Vinterp = (1.0 - w) * V[k, c] + w * V[k + 1, c]
            if surplus <= 1e-10:
                Vo[b, c] = -1e10
                co[b, c] = c_bar_0 + c_min
                ho[b, c] = h_bar_0 + 0.01
            elif surplus > cap_c:
                ct = Rvb - cbc - ri * hR_max - bp
                if ct < 1e-10:
                    ct = 1e-10
                comp = (ct ** alpha) * ht_cap_c_pow
                u = (comp ** oms) * inv_oms + psic
                Vo[b, c] = u + beta * Vinterp
                ct_eff = ct if ct > c_min else c_min
                co[b, c] = cbc + ct_eff
                ht_eff = ht_cap_c if ht_cap_c > 0.01 else 0.01
                ho[b, c] = hbc + ht_eff
            else:
                ss = surplus
                u = Kr * ss ** oms * inv_oms + psic
                Vo[b, c] = u + beta * Vinterp
                ct = alpha * surplus
                ht = (1.0 - alpha) * surplus * inv_ri
                ct_eff = ct if ct > c_min else c_min
                ht_eff = ht if ht > 0.01 else 0.01
                co[b, c] = cbc + ct_eff
                ho[b, c] = hbc + ht_eff
    return Vo, co, ho


@njit(cache=True)
def _interp_with_clip(bg, V, x):
    nb = bg.size
    if x <= bg[0]:
        return V[0]
    if x >= bg[nb - 1]:
        return V[nb - 1]
    lo = 0
    hi = nb - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bg[mid] <= x:
            lo = mid
        else:
            hi = mid
    w = (x - bg[lo]) / (bg[lo + 1] - bg[lo])
    if w < 0.0:
        w = 0.0
    elif w > 1.0:
        w = 1.0
    return (1.0 - w) * V[lo] + w * V[lo + 1]


@njit(cache=True)
def tenure_choice_kernel(
    Vd,                  # (Nb, nt, I, npar, ncs)
    b_grid,              # (Nb,)
    heq,                 # (I, nt)
    hcost,               # (I, nt)
    dp_arr,              # (I, nt, npar, ncs)
    bmo,                 # (I, nt, npar, ncs)
    birth_dp,            # (npar, ncs, nt, nt) bool
    birth_entry_grant,   # (I, nt, npar, ncs)
):
    # Discrete tenure-choice argmax over `tn` given conditional values Vd
    # for each (origin tenure `to`, location, parity, child-state, b).
    # Three branches handle: stay/move-as-renter (tn=0), buy-on-entry
    # (to=0 -> tn>=1, with optional birth grant or entry grant), and
    # sell-then-rebuy (to>=1 -> tn != to). Infeasible states (below
    # down-payment threshold dp or below borrowing limit bm) get -1e10.
    Nb, nt, I, npar, ncs = Vd.shape
    VH = np.empty((Nb, nt, I, npar, ncs))
    tcj = np.empty((Nb, nt, I, npar, ncs), dtype=np.int16)
    NEG_INF = -1e10
    for id_ in range(I):
        for to in range(nt):
            sp = heq[id_, to] if to > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    dpn0_unused = 0.0  # placeholder; per-tn dp/bm pulled inside tn loop
                    for b in range(Nb):
                        bg_b = b_grid[b]
                        best_v = NEG_INF
                        best_tn = 0
                        # tn = 0 (renter)
                        if to == 0:
                            v0 = Vd[b, 0, id_, nn, cs]
                        else:
                            ba = bg_b + sp
                            if ba < 0.0:
                                ba = 0.0
                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba)
                        if v0 > best_v:
                            best_v = v0
                            best_tn = 0
                        # tn >= 1 (owner tenures)
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            dpn = dp_arr[id_, tn, nn, cs]
                            bmn = bmo[id_, tn, nn, cs]
                            if to == tn:
                                v_tn = Vd[b, tn, id_, nn, cs]
                            elif to == 0:
                                bab = bg_b - hc
                                if birth_dp[nn, cs, to, tn]:
                                    bag = bab if bab > bmn else bmn
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bag)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], babg)
                                    if (bg_b + gfix) < dpn or babg < bmn:
                                        v_tn = NEG_INF
                                else:
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bab)
                                    if bg_b < dpn or bab < bmn:
                                        v_tn = NEG_INF
                            else:
                                bar = bg_b + sp - hc
                                v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bar)
                                dpc = dpn - sp
                                if bg_b < dpc or bar < bmn:
                                    v_tn = NEG_INF
                            if v_tn > best_v:
                                best_v = v_tn
                                best_tn = tn
                        VH[b, to, id_, nn, cs] = best_v
                        tcj[b, to, id_, nn, cs] = best_tn
    return VH, tcj


@njit(cache=True, parallel=True)
def full_renter_block_kernel(
    Rv1d,           # (Nb,)
    Vc_flat,        # (Nb, nc)
    bp_prev,        # (Nb, nc) or zeros (use has_prev to gate)
    has_prev,       # bool/int
    b_grid,         # (Nb,)
    cb_v,           # (nc,)
    hb_v,           # (nc,)
    psi_v,          # (nc,)
    ri,
    hR_max,
    c_min,
    c_bar_0,
    h_bar_0,
    alpha,
    oms,
    beta,
    gs_alpha1,
    gs_alpha2,
    gs_tol,
):
    # Full-Bellman renter block: golden-section search for bp + post-search
    # consumption / housing arithmetic, fused into one kernel per (i, j).
    # When bp_prev is given (j < J-1), the search interval is clamped to
    # [bp_prev - 2, bp_prev + 2] as a soft monotonicity prior — a
    # heuristic that mirrors the MATLAB implementation, not a strict
    # invariant of the model.
    Nb, nc = Vc_flat.shape
    Vo = np.empty((Nb, nc))
    bp_out = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    ho = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    inv_ri = 1.0 / ri
    Kr = (alpha ** alpha * ((1.0 - alpha) * inv_ri) ** (1.0 - alpha)) ** oms
    bg0 = b_grid[0]
    for c in prange(nc):
        cbc = cb_v[c]
        hbc = hb_v[c]
        psic = psi_v[c]
        dc = cbc + ri * hbc
        cap_c = ri * (hR_max - hbc) / (1.0 - alpha)
        ht_cap_c = hR_max - hbc
        if ht_cap_c < 1e-10:
            ht_cap_c = 1e-10
        ht_cap_pow = ht_cap_c ** (1.0 - alpha)
        for b in range(Nb):
            Rvb = Rv1d[b]
            lo = 0.0
            if bg0 > lo:
                lo = bg0
            hi = Rvb - dc - 1e-6
            if hi < 0.0:
                hi = 0.0
            if has_prev:
                lo_prev = bp_prev[b, c] - 2.0
                if lo_prev > lo:
                    lo = lo_prev
                hi_prev = bp_prev[b, c] + 2.0
                if hi_prev < hi:
                    hi = hi_prev
                if lo < 0.0:
                    lo = 0.0
                if hi < lo:
                    hi = lo

            d = hi - lo
            x1 = lo + gs_alpha1 * d
            x2 = lo + gs_alpha2 * d
            f1 = eval_renter_scalar(x1, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, alpha, oms, beta)
            f2 = eval_renter_scalar(x2, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, alpha, oms, beta)
            d = gs_alpha1 * gs_alpha2 * d
            while d > gs_tol:
                if f2 >= f1:
                    xe = x2 + d
                    fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, alpha, oms, beta)
                    x1 = x2
                    f1 = f2
                    x2 = xe
                    f2 = fe
                else:
                    xe = x1 - d
                    fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, alpha, oms, beta)
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
                ht_unc = (1.0 - alpha) * surplus * inv_ri
                if hbc + ht_unc > hR_max:
                    ct = Rvb - cbc - ri * hR_max - bp_best
                    if ct < 1e-10:
                        ct = 1e-10
                    ct_eff = ct if ct > c_min else c_min
                    ht_eff = ht_cap_c if ht_cap_c > 0.01 else 0.01
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
                else:
                    ct = alpha * surplus
                    ct_eff = ct if ct > c_min else c_min
                    ht_eff = ht_unc if ht_unc > 0.01 else 0.01
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
    return Vo, bp_out, co, ho


@njit(cache=True, parallel=True)
def full_owner_block_kernel(
    Rv1d,           # (Nb,)
    Vco_flat,       # (Nb, nc)
    bp_prev,        # (Nb, nc) or zeros
    has_prev,
    b_grid,
    cb_v,           # (nc,)
    hb_v,           # (nc,)
    psi_v,          # (nc,)
    bf_v,           # (nc,) — bmo[i, ten, nn, cs] flattened in F order
    oc,
    hsv,
    owner_h_bar_scale,
    c_min,
    alpha,
    oms,
    beta,
    gs_alpha1,
    gs_alpha2,
    gs_tol,
):
    Nb, nc = Vco_flat.shape
    Vo = np.empty((Nb, nc))
    bp_out = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    aoms = alpha * oms
    one_minus_alpha_oms = (1.0 - alpha) * oms
    bg0 = b_grid[0]
    for c in prange(nc):
        cbc = cb_v[c]
        hbc = hb_v[c]
        psic = psi_v[c]
        bf = bf_v[c]
        ht_c = hsv - owner_h_bar_scale * hbc
        if ht_c < 1e-10:
            ht_c = 1e-10
        Ko_c = ht_c ** one_minus_alpha_oms
        for b in range(Nb):
            Rvb = Rv1d[b]
            lo = bf
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
                if lo < bf:
                    lo = bf
                if hi < lo:
                    hi = lo

            d = hi - lo
            x1 = lo + gs_alpha1 * d
            x2 = lo + gs_alpha2 * d
            f1 = eval_owner_scalar(x1, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, alpha, oms, beta)
            f2 = eval_owner_scalar(x2, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, alpha, oms, beta)
            d = gs_alpha1 * gs_alpha2 * d
            while d > gs_tol:
                if f2 >= f1:
                    xe = x2 + d
                    fe = eval_owner_scalar(xe, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, alpha, oms, beta)
                    x1 = x2
                    f1 = f2
                    x2 = xe
                    f2 = fe
                else:
                    xe = x1 - d
                    fe = eval_owner_scalar(xe, Rvb, Vco_flat[:, c], b_grid, oc, cbc, psic, Ko_c, alpha, oms, beta)
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
            co[b, c] = cbc + ct_eff
    return Vo, bp_out, co


@njit(cache=True)
def location_logit_kernel(
    VH,         # (Nb, nt, I, npar, ncs)
    iidx,       # (Nb, I, nt) int
    iwt,        # (Nb, I, nt)
    loc_shift,  # (I, I)
    kappa_loc,
):
    Nb, nt, I, npar, ncs = VH.shape
    VI = np.empty((Nb, nt, I, npar, ncs))
    lpj = np.empty((Nb, nt, I, I, npar, ncs))
    inv_kl = 1.0 / kappa_loc
    for io in range(I):
        for to in range(nt):
            for nn in range(npar):
                for cs in range(ncs):
                    for b in range(Nb):
                        # Build value-to-go for each destination, then logit
                        # First find max for stable logsumexp
                        m = -1e300
                        for idd in range(I):
                            if idd == io:
                                v = VH[b, to, io, nn, cs]
                            else:
                                k = iidx[b, io, to]
                                w = iwt[b, io, to]
                                v = (1.0 - w) * VH[k, 0, idd, nn, cs] + w * VH[k + 1, 0, idd, nn, cs]
                            v_shift = (v + loc_shift[io, idd]) * inv_kl
                            if v_shift > m:
                                m = v_shift
                        # accumulate exp sum
                        se = 0.0
                        # store the shifted values temporarily in lpj
                        for idd in range(I):
                            if idd == io:
                                v = VH[b, to, io, nn, cs]
                            else:
                                k = iidx[b, io, to]
                                w = iwt[b, io, to]
                                v = (1.0 - w) * VH[k, 0, idd, nn, cs] + w * VH[k + 1, 0, idd, nn, cs]
                            v_shift = (v + loc_shift[io, idd]) * inv_kl
                            ex = np.exp(v_shift - m)
                            lpj[b, to, io, idd, nn, cs] = ex
                            se += ex
                        ls = m + np.log(se)
                        VI[b, to, io, nn, cs] = kappa_loc * ls
                        # normalize probs
                        for idd in range(I):
                            lpj[b, to, io, idd, nn, cs] = lpj[b, to, io, idd, nn, cs] / se
    return VI, lpj


@njit(cache=True, parallel=True)
def eval_owner_block_kernel(
    Rv1d,
    bpv_o,
    V,
    idx,
    wt,
    cb,
    hb,
    psi_v,
    oc,
    hsv,
    owner_h_bar_scale,
    c_min,
    alpha,
    oms,
    beta,
):
    Nb, nc = bpv_o.shape
    Vo = np.empty((Nb, nc))
    co = np.empty((Nb, nc))
    inv_oms = 1.0 / oms
    aoms = alpha * oms
    one_minus_alpha_oms = (1.0 - alpha) * oms
    for c in prange(nc):
        cbc = cb[c]
        hbc = hb[c]
        psic = psi_v[c]
        ht_c = hsv - owner_h_bar_scale * hbc
        if ht_c < 1e-10:
            ht_c = 1e-10
        Ko = ht_c ** one_minus_alpha_oms
        for b in range(Nb):
            bp = bpv_o[b, c]
            Rvb = Rv1d[b]
            ct_raw = Rvb - oc - cbc - bp
            ct = ct_raw if ct_raw > 1e-10 else 1e-10
            k = idx[b, c]
            w = wt[b, c]
            Vinterp = (1.0 - w) * V[k, c] + w * V[k + 1, c]
            if ct_raw <= 1e-10:
                Vo[b, c] = -1e10
            else:
                Vo[b, c] = Ko * (ct ** aoms) * inv_oms + psic + beta * Vinterp
            ct_eff = ct if ct > c_min else c_min
            co[b, c] = cbc + ct_eff
    return Vo, co
