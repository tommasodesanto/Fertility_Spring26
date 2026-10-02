"""Extracted from kernels.py (reviewed overlay) (sha256 379d179a8f80477e8d813f60bc526e700a270e37e0f411988aabd9a06fd74352).

Mechanical copy by refactor_lab/materialize.py: only reachable top-level
definitions, bodies byte-identical; import edits listed in the receipt.
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
def renter_wedge_flow(S_raw, cbc, hbc, ri, w0, w1, hk, hRmax, al, oms, es):
    """Intratemporal renter allocation under the size-dependent rental wedge.

    Total housing cost is C(h) = h*(ri + w0) + w1*h*max(0, h-hk); S_raw is
    resources net of saving but before committed (cbc, hbc) spending. Below
    the knee the allocation is Cobb-Douglas with effective rent ri + w0;
    above it the quadratic cost yields a closed-form root; at the knee
    itself there is a bunching interval (the marginal cost jumps), handled
    by comparing candidate utilities; beyond hRmax the household is capped
    at hRmax. Returns (u_flow, ct, ht) with the eqscale exponent applied to
    utility but the family shifter left to the caller. Mirrored in Python
    by solver.renter_wedge_flow_py; the two are cross-tested.
    """
    ri1 = ri + w0
    Cc = hbc * ri1 + w1 * hbc * (hbc - hk if hbc > hk else 0.0)
    S = S_raw - cbc - Cc
    if S <= 1e-10:
        return -1e10, 0.0, 0.0
    Kr1 = (al ** al * ((1.0 - al) / ri1) ** (1.0 - al)) ** oms
    # Candidate 1: interior below the knee.
    ht_a = (1.0 - al) * S / ri1
    h_a = hbc + ht_a
    u_a = Kr1 * S ** oms / oms
    if es != 1.0:
        u_a = es * u_a
    valid_a = (h_a <= hk)
    # Candidate 2: bunching at the knee (affordable and below the cap).
    ht_k = hk - hbc
    valid_k = (hbc < hk) and (hk <= hRmax) and (S >= ri1 * ht_k)
    if valid_k:
        ct_k = S - ri1 * ht_k
        if ct_k < 1e-10:
            ct_k = 1e-10
        u_k = (ct_k ** al * ht_k ** (1.0 - al)) ** oms / oms
        if es != 1.0:
            u_k = es * u_k
    else:
        u_k = -1e10
        ct_k = 0.0
    # Candidate 3: interior above the knee (quadratic root).
    if hbc < hk:
        S_eff = S + w1 * hbc * (hk - hbc)
    else:
        S_eff = S
    A2 = ri1 + w1 * (2.0 * hbc - hk)
    quad_b = w1 * (1.0 + al)
    if quad_b > 0.0:
        disc = A2 * A2 + 4.0 * quad_b * (1.0 - al) * S_eff
        ht_b = (-A2 + np.sqrt(disc)) / (2.0 * quad_b)
    else:
        ht_b = (1.0 - al) * S_eff / A2
    h_b = hbc + ht_b
    valid_b = (h_b > hk)
    if valid_b:
        ct_b = S_eff - A2 * ht_b - w1 * ht_b * ht_b
        if ct_b < 1e-10:
            ct_b = 1e-10
        u_b = (ct_b ** al * ht_b ** (1.0 - al)) ** oms / oms
        if es != 1.0:
            u_b = es * u_b
    else:
        u_b = -1e10
        ct_b = 0.0
    # Unconstrained pick among valid candidates.
    u_best = u_a if valid_a else -1e10
    ct_best = al * S if valid_a else 0.0
    ht_best = ht_a if valid_a else 0.0
    if valid_k and u_k > u_best:
        u_best = u_k
        ct_best = ct_k
        ht_best = ht_k
    if valid_b and u_b > u_best:
        u_best = u_b
        ct_best = ct_b
        ht_best = ht_b
    h_best = hbc + ht_best
    if h_best > hRmax:
        Ccap = hRmax * ri1 + w1 * hRmax * (hRmax - hk if hRmax > hk else 0.0)
        ct = S_raw - cbc - Ccap
        if ct < 1e-10:
            ct = 1e-10
        ht_use = hRmax - hbc
        if ht_use < 1e-10:
            ht_use = 1e-10
        u = (ct ** al * ht_use ** (1.0 - al)) ** oms / oms
        if es != 1.0:
            u = es * u
        return u, ct, ht_use
    if u_best <= -1e9:
        return -1e10, 0.0, 0.0
    return u_best, ct_best, ht_best


@njit(cache=True)
def eval_renter_scalar(bp, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, es,
                       hbc=0.0, w0=0.0, w1=0.0, hk=6.0, wedge_on=0):
    if wedge_on != 0:
        u_flow, _, _ = renter_wedge_flow(Rv - bp, cb_c, hbc, ri, w0, w1, hk, hRmax, alpha, oms, es)
        if u_flow <= -1e9:
            return -1e10
        return u_flow + pc + beta * interp_scalar(bg, Vbar, bp)
    surplus = Rv - dc - bp
    if surplus <= 1e-10:
        return -1e10
    if surplus > cc:
        ct = Rv - cb_c - ri * hRmax - bp
        if ct < 1e-10:
            ct = 1e-10
        if es != 1.0:
            u = es * (ct**alpha * ht_cap_c ** (1.0 - alpha)) ** oms / oms + pc
        else:
            u = (ct**alpha * ht_cap_c ** (1.0 - alpha)) ** oms / oms + pc
    else:
        ss = surplus
        if ss < 1e-10:
            ss = 1e-10
        if es != 1.0:
            u = es * Kr * ss**oms / oms + pc
        else:
            u = Kr * ss**oms / oms + pc
    return u + beta * interp_scalar(bg, Vbar, bp)


@njit(cache=True)
def eval_owner_scalar(bp, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es):
    ct = Rv - oc - cb_c - bp
    if ct <= 1e-10:
        return -1e10
    if es != 1.0:
        return es * Ko_c * ct ** (alpha * oms) / oms + pc + beta * interp_scalar(bg, Vbar, bp)
    return Ko_c * ct ** (alpha * oms) / oms + pc + beta * interp_scalar(bg, Vbar, bp)


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
def _interp_with_clip(bg, V, x, strict_interpolated_support=False, transaction_support=False):
    if transaction_support and (x < bg[0] or x > bg[-1]):
        return -1e10
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
    if strict_interpolated_support:
        # A finite infeasibility sentinel must never become a feasible
        # interpolant when its destination receives positive scatter mass.
        if ((1.0-w > 0.0 and V[lo] <= -1e9)
                or (w > 0.0 and V[lo+1] <= -1e9)):
            return -1e10
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
    Vd_stay,             # (Nb, nt, I, npar, ncs) values read at to == tn
    strict_interpolated_support=False,
    transaction_support=False,
    require_owner_sale_solvency=False,
):
    # Discrete tenure-choice argmax over `tn` given conditional values Vd
    # for each (origin tenure `to`, location, parity, child-state, b).
    # Three branches handle: stay/move-as-renter (tn=0), buy-on-entry
    # (to=0 -> tn>=1, with optional birth grant or entry grant), and
    # sell-then-rebuy (to>=1 -> tn != to). Infeasible states (below
    # down-payment threshold dp or below borrowing limit bm) get -1e10.
    # A stayer (to == tn) reads Vd_stay, which carries the stayer mortgage
    # floors when that switch is on and matches Vd bit for bit otherwise.
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
                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, strict_interpolated_support, transaction_support)
                            if require_owner_sale_solvency and bg_b + sp < 0.0:
                                v0 = NEG_INF
                        if v0 > best_v:
                            best_v = v0
                            best_tn = 0
                        # tn >= 1 (owner tenures)
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            dpn = dp_arr[id_, tn, nn, cs]
                            bmn = bmo[id_, tn, nn, cs]
                            if to == tn:
                                v_tn = Vd_stay[b, tn, id_, nn, cs]
                            elif to == 0:
                                bab = bg_b - hc
                                if birth_dp[nn, cs, to, tn]:
                                    bag = bab if bab > bmn else bmn
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bag, strict_interpolated_support, transaction_support)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], babg, strict_interpolated_support, transaction_support)
                                    if (bg_b + gfix) < dpn or babg < bmn:
                                        v_tn = NEG_INF
                                else:
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bab, strict_interpolated_support, transaction_support)
                                    if bg_b < dpn or bab < bmn:
                                        v_tn = NEG_INF
                            else:
                                bar = bg_b + sp - hc
                                v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bar, strict_interpolated_support, transaction_support)
                                dpc = dpn - sp
                                if bg_b < dpc or bar < bmn:
                                    v_tn = NEG_INF
                            if v_tn > best_v:
                                best_v = v_tn
                                best_tn = tn
                        VH[b, to, id_, nn, cs] = best_v
                        tcj[b, to, id_, nn, cs] = best_tn
    return VH, tcj


@njit(cache=True)
def tenure_logit_kernel(
    Vd,                  # (Nb, nt, I, npar, ncs)
    b_grid,              # (Nb,)
    heq,                 # (I, nt)
    hcost,               # (I, nt)
    dp_arr,              # (I, nt, npar, ncs)
    bmo,                 # (I, nt, npar, ncs)
    birth_dp,            # (npar, ncs, nt, nt) bool
    birth_entry_grant,   # (I, nt, npar, ncs)
    kappa,               # taste-shock scale
    Vd_stay,             # (Nb, nt, I, npar, ncs) values read at to == tn
    transaction_support=False,
    require_owner_sale_solvency=False,
):
    Nb, nt, I, npar, ncs = Vd.shape
    VH = np.empty((Nb, nt, I, npar, ncs))
    tcj = np.empty((Nb, nt, I, npar, ncs), dtype=np.int16)
    probs = np.zeros((Nb, nt, I, npar, ncs, nt), dtype=np.float32)
    NEG_INF = -1e10
    tiny_kappa = kappa if kappa > 1e-12 else 1e-12
    vals = np.empty(nt)
    for id_ in range(I):
        for to in range(nt):
            sp = heq[id_, to] if to > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    for b in range(Nb):
                        bg_b = b_grid[b]
                        best_v = NEG_INF
                        best_tn = 0
                        if to == 0:
                            v0 = Vd[b, 0, id_, nn, cs]
                        else:
                            ba = bg_b + sp
                            v0 = _interp_with_clip(b_grid, Vd[:, 0, id_, nn, cs], ba, False, transaction_support)
                            if require_owner_sale_solvency and bg_b + sp < 0.0:
                                v0 = NEG_INF
                        vals[0] = v0
                        if v0 > best_v:
                            best_v = v0
                            best_tn = 0
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            dpn = dp_arr[id_, tn, nn, cs]
                            bmn = bmo[id_, tn, nn, cs]
                            if to == tn:
                                v_tn = Vd_stay[b, tn, id_, nn, cs]
                            elif to == 0:
                                bab = bg_b - hc
                                if birth_dp[nn, cs, to, tn]:
                                    bag = bab if bab > bmn else bmn
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bag, False, transaction_support)
                                elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                    gfix = birth_entry_grant[id_, tn, nn, cs]
                                    babg = bab + gfix
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], babg, False, transaction_support)
                                    if (bg_b + gfix) < dpn or babg < bmn:
                                        v_tn = NEG_INF
                                else:
                                    v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bab, False, transaction_support)
                                    if bg_b < dpn or bab < bmn:
                                        v_tn = NEG_INF
                            else:
                                bar = bg_b + sp - hc
                                v_tn = _interp_with_clip(b_grid, Vd[:, tn, id_, nn, cs], bar, False, transaction_support)
                                dpc = dpn - sp
                                if bg_b < dpc or bar < bmn:
                                    v_tn = NEG_INF
                            vals[tn] = v_tn
                            if v_tn > best_v:
                                best_v = v_tn
                                best_tn = tn
                        tcj[b, to, id_, nn, cs] = best_tn
                        if best_v <= -1e9:
                            VH[b, to, id_, nn, cs] = best_v
                            for tn in range(nt):
                                probs[b, to, id_, nn, cs, tn] = 0.0
                        else:
                            denom = 0.0
                            for tn in range(nt):
                                ev = np.exp((vals[tn] - best_v) / tiny_kappa)
                                probs[b, to, id_, nn, cs, tn] = ev
                                denom += ev
                            if denom <= 0.0:
                                VH[b, to, id_, nn, cs] = best_v
                                probs[b, to, id_, nn, cs, best_tn] = 1.0
                            else:
                                for tn in range(nt):
                                    probs[b, to, id_, nn, cs, tn] = probs[b, to, id_, nn, cs, tn] / denom
                                VH[b, to, id_, nn, cs] = best_v + tiny_kappa * np.log(denom)
    return VH, tcj, probs


@njit(cache=True)
def _rank(bg, x):
    """Largest l with bg[l] <= x, or -1."""
    lo = -1
    hi = bg.size
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bg[mid] <= x:
            lo = mid
        else:
            hi = mid
    return lo


@njit(cache=True)
def _interp_ranked(bg, V, x, r):
    # interp_scalar with its binary-search index supplied: clamp(r, 0, nb-2).
    idx = r
    if idx < 0:
        idx = 0
    elif idx > bg.size - 2:
        idx = bg.size - 2
    wt = (x - bg[idx]) / (bg[idx + 1] - bg[idx])
    if wt < 0.0:
        wt = 0.0
    elif wt > 1.0:
        wt = 1.0
    return (1.0 - wt) * V[idx] + wt * V[idx + 1]


@njit(cache=True)
def _renter_value(bp, r, Rv, Vbar, bg, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, es):
    # eval_renter_scalar, wedge_on == 0 branch, same operation order.
    surplus = Rv - dc - bp
    if surplus <= 1e-10:
        return -1e10
    if surplus > cc:
        ct = Rv - cb_c - ri * hRmax - bp
        if ct < 1e-10:
            ct = 1e-10
        if es != 1.0:
            u = es * (ct**alpha * ht_cap_c ** (1.0 - alpha)) ** oms / oms + pc
        else:
            u = (ct**alpha * ht_cap_c ** (1.0 - alpha)) ** oms / oms + pc
    else:
        ss = surplus
        if ss < 1e-10:
            ss = 1e-10
        if es != 1.0:
            u = es * Kr * ss**oms / oms + pc
        else:
            u = Kr * ss**oms / oms + pc
    return u + beta * _interp_ranked(bg, Vbar, bp, r)


@njit(cache=True)
def _owner_value(bp, r, Rv, Vbar, bg, oc, cb_c, pc, Ko_c, alpha, oms, beta, es):
    ct = Rv - oc - cb_c - bp
    if ct <= 1e-10:
        return -1e10
    if es != 1.0:
        return es * Ko_c * ct ** (alpha * oms) / oms + pc + beta * _interp_ranked(bg, Vbar, bp, r)
    return Ko_c * ct ** (alpha * oms) / oms + pc + beta * _interp_ranked(bg, Vbar, bp, r)


@njit(cache=True)
def exhaustive_saving_scalar(lo, hi, resources, continuation, bg, rent, hb, cb, pc,
                              hmax, alpha, oms, beta, es, owner_cost, owner_K, owner):
    if not (0.0 < alpha < 1.0 and oms < 1.0 and oms != 0.0
            and beta > 0.0 and es > 0.0 and rent > 0.0 and hi >= lo):
        raise ValueError("Unsupported exhaustive-saving objective or interval")
    dc=cb+rent*hb
    cap=rent*(hmax-hb)/(1-alpha)
    hcap=max(hmax-hb,1e-10)
    Kr=(alpha**alpha*((1-alpha)/rent)**(1-alpha))**oms
    nb = bg.size
    # Grid nodes strictly inside (lo, hi) are the contiguous block [g0, g1).
    g0 = _rank(bg, lo) + 1
    r_hi = _rank(bg, hi)
    g1 = r_hi + 1
    if g1 > g0 and bg[g1 - 1] >= hi:
        g1 -= 1
    has_kink = False
    kink = 0.0
    r_kink = -1
    if not owner:
        kink = resources-dc-cap
        if lo < kink < hi:
            has_kink = True
            r_kink = _rank(bg, kink)
    points = np.empty(max(g1 - g0, 0) + 3)
    ranks = np.empty(points.size, dtype=np.int64)
    n = 0
    points[n] = lo; ranks[n] = g0 - 1; n += 1
    for i in range(g0, g1):
        if has_kink and kink < bg[i]:
            points[n] = kink; ranks[n] = r_kink; n += 1
            has_kink = False
        points[n] = bg[i]; ranks[n] = i; n += 1
    if has_kink:
        points[n] = kink; ranks[n] = r_kink; n += 1
    points[n] = hi; ranks[n] = r_hi; n += 1
    best=-1e300; bpbest=lo
    for k in range(n):
        x=points[k]
        r=ranks[k]
        if owner:
            value=_owner_value(x,r,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
        else:
            value=_renter_value(x,r,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
        if value>best:
            best=value;bpbest=x
        if k==n-1 or points[k+1]-x<1e-14:
            continue
        mid=(x+points[k+1])/2
        if mid<=bg[0] or mid>=bg[nb-1]:
            slope=0.0
        else:
            ix=r
            slope=(continuation[ix+1]-continuation[ix])/(bg[ix+1]-bg[ix])
        if slope<=0:
            continue
        if owner:
            optimal_c=(beta*slope/(es*owner_K*alpha))**(1/(alpha*oms-1))
            candidate=resources-owner_cost-cb-optimal_c
        elif resources-dc-mid>cap:
            K=hcap**((1-alpha)*oms)
            optimal_c=(beta*slope/(es*K*alpha))**(1/(alpha*oms-1))
            candidate=resources-cb-rent*hmax-optimal_c
        else:
            optimal_surplus=(beta*slope/(es*Kr))**(1/(oms-1))
            candidate=resources-dc-optimal_surplus
        if x<candidate<points[k+1]:
            if owner:
                value=_owner_value(candidate,r,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
            else:
                value=_renter_value(candidate,r,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
            if value>best:
                best=value;bpbest=candidate
    return bpbest,best


@njit(cache=True, parallel=True)
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
    yadj_v=None,
    pen_on=0,
    wedge_on=0,
    w0=0.0,
    w1=0.0,
    hk=6.0,
    exact_allocation_output=False,
    natural_floor_v=None,
    fixed_renter_floor=-np.inf,
):
    # Full-Bellman renter block: golden-section search for bp + post-search
    # consumption / housing arithmetic, fused into one kernel per (i, j).
    # yadj_v/pen_on carry the optional children-at-home earnings adjustment:
    # resources for family cell c rise by yadj_v[c]. With pen_on == 0 the
    # adjustment is skipped and the block is bit for bit the legacy one.
    # wedge_on/w0/w1/hk carry the optional size-dependent rental wedge (see
    # renter_wedge_flow). With wedge_on == 0 the block is bit for bit legacy.
    # The wedge is not supported under exhaustive saving (its kink/candidate
    # structure assumes linear rent); that combination raises below.
    # When bp_prev is given (j < J-1), the search interval is clamped to
    # [bp_prev - 2, bp_prev + 2] as a soft monotonicity prior — a
    # heuristic that mirrors the MATLAB implementation, not a strict
    # invariant of the model.
    if exhaustive_saving and has_prev:
        raise ValueError("Exhaustive saving requires the full feasible interval")
    if exhaustive_saving and wedge_on != 0:
        raise ValueError("Rental wedge requires the golden-section renter block")
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
            if pen_on != 0:
                if yadj_v is not None:
                    Rvb = Rv1d[b] + yadj_v[c]
                    Rvtb = Rvt1d[b] + yadj_v[c]
                else:
                    Rvb = Rv1d[b]
                    Rvtb = Rvt1d[b]
            else:
                Rvb = Rv1d[b]
                Rvtb = Rvt1d[b]
            if gc > 0.0:
                Tb = gc - Rvtb
                if Tb > 0.0:
                    if Tb > gc:
                        Tb = gc
                    Rvb = Rvb + Tb
            current_b = b_grid[b]
            rollover_floor = s_next * (current_b if current_b < 0.0 else 0.0)
            line_floor = -D_next
            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor
            if np.isfinite(fixed_renter_floor):
                unsecured_floor = fixed_renter_floor
            if natural_floor_v is not None:
                unsecured_floor = natural_floor_v[c]
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
                f1 = eval_renter_scalar(x1, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es, hbc, w0, w1, hk, wedge_on)
                f2 = eval_renter_scalar(x2, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es, hbc, w0, w1, hk, wedge_on)
                d = gs_alpha1 * gs_alpha2 * d
                while d > gs_tol:
                    if f2 >= f1:
                        xe = x2 + d
                        if xe > hi:
                            xe = hi
                        fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es, hbc, w0, w1, hk, wedge_on)
                        x1 = x2
                        f1 = f2
                        x2 = xe
                        f2 = fe
                    else:
                        xe = x1 - d
                        if xe < lo:
                            xe = lo
                        fe = eval_renter_scalar(xe, Rvb, Vc_flat[:, c], b_grid, dc, psic, cap_c, cbc, ri, hR_max, ht_cap_c, Kr, al, oms, beta, es, hbc, w0, w1, hk, wedge_on)
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

            if wedge_on != 0:
                uw, ctw, htw = renter_wedge_flow(Rvb - bp_best, cbc, hbc, ri, w0, w1, hk, hR_max, al, oms, es)
                if uw <= -1e9:
                    co[b, c] = c_bar_0 + c_min
                    ho[b, c] = h_bar_0 + 0.01
                else:
                    ct_eff = ctw if ctw > c_min else c_min
                    ht_eff = htw if htw > 0.01 else 0.01
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
                continue
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
                    if exhaustive_saving and (exact_allocation_output or v_best > -1e9):
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
                    if exhaustive_saving and (exact_allocation_output or v_best > -1e9):
                        ct_eff = ct
                        ht_eff = ht_unc
                    co[b, c] = cbc + ct_eff
                    ho[b, c] = hbc + ht_eff
    return Vo, bp_out, co, ho


@njit(cache=True, parallel=True)
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
    yadj_v=None,
    pen_on=0,
    stay_on=0,
    stay_orig=0,
    amort=0.0,
    exact_allocation_output=False,
    due_stayer=False,
    due_death_floor=-np.inf,
    purchase_saving_fraction=1.0,
):
    # yadj_v/pen_on carry the optional children-at-home earnings adjustment;
    # see full_renter_block_kernel. With pen_on == 0 the block matches the
    # legacy one bit for bit.
    # stay_on/stay_orig/amort carry the optional stayer mortgage floors: a
    # stayer (origin tenure == destination tenure) faces a no-cash-out rule
    # instead of the origination collateral floor. With stay_on == 0 the
    # floor block matches the legacy one bit for bit.
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
            if pen_on != 0:
                if yadj_v is not None:
                    Rvb = Rv1d[b] + yadj_v[c]
                    Rvtb = Rvt1d[b] + yadj_v[c]
                else:
                    Rvb = Rv1d[b]
                    Rvtb = Rvt1d[b]
            else:
                Rvb = Rv1d[b]
                Rvtb = Rvt1d[b]
            if gc > 0.0:
                Tb = gc - Rvtb
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
            if stay_on != 0:
                # Stayer mortgage floor: debt may not rise (no cash-out), and
                # with amortization it must fall by at least that share. The
                # taper/line rollover logic is bypassed: a stayer is never
                # forced below its own balance. With no debt, a first mortgage
                # on the owned house is an origination, so the collateral
                # floor applies when stay_orig is on.
                bb = b_grid[b]
                if bb < 0.0:
                    amort_floor = bb * (1.0 - amort)
                    if stay_orig != 0:
                        total_floor = amort_floor if amort_floor > bb else bb
                    elif total_floor < amort_floor:
                        total_floor = amort_floor
                elif stay_orig != 0:
                    total_floor = bf
            if due_stayer:
                # b is before interest; allowing R*b would capitalize interest.
                total_floor = min(b_grid[b], bf)
                total_floor = max(total_floor, due_death_floor)
            if purchase_saving_fraction < 1.0 and stay_on == 0 and not due_stayer:
                # Buyer rule A + lambda*S >= (1-phi)Q, where x=A-Q.
                # Its saving floor is bf + (1-lambda)/lambda * max(0, bf-x).
                x = b_grid[b]
                quarter_floor = bf + (1.0 - purchase_saving_fraction) / purchase_saving_fraction * max(0.0, bf - x)
                if quarter_floor > total_floor:
                    total_floor = quarter_floor
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
            if exact_allocation_output and exhaustive_saving and ct > 1e-10:
                ct_eff = ct
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
                        all_dead = True
                        for idd in range(I):
                            if idd == io:
                                v = VH[b, to, io, nn, cs]
                            else:
                                k = iidx[b, io, to]
                                w = iwt[b, io, to]
                                v = (1.0 - w) * VH[k, 0, idd, nn, cs] + w * VH[k + 1, 0, idd, nn, cs]
                            if v > -1e9:
                                all_dead = False
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
                            if all_dead:
                                lpj[b, to, io, idd, nn, cs] = 0.0
                            else:
                                lpj[b, to, io, idd, nn, cs] = lpj[b, to, io, idd, nn, cs] / se
    return VI, lpj
