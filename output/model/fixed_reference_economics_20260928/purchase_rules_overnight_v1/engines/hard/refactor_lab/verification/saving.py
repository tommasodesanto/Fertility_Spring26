"""Global saving choice on one continuation segment set.

The first five functions are verbatim extractions from
`intergen_eqscale_seq_optimized/kernels.py`; the acceptance suite asserts
source equality with the oracle, so any drift fails. The maintained Bellman
calls `exhaustive_saving_scalar` with `wedge_on == 0` (the wedge path raises
under exhaustive saving in the oracle).

`exhaustive_saving_indexed` is the proposed optimization, off unless selected.
It evaluates the same candidates in the same order with the same strict `>`
tie rule, but replaces per-candidate binary search, the `np.sort` of the merged
point list, and the per-segment `searchsorted` by rank bookkeeping (three
binary searches per call instead of about two per grid node). Justification:
between consecutive sorted points no grid node lies strictly inside, so the
interpolation index of a first-order candidate and the slope segment of a
midpoint both equal the rank of the left point. Bit identity must be shown by
the acceptance fixtures before it is used anywhere.
"""
from __future__ import annotations

import numpy as np

try:  # pragma: no cover
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*args, **kwargs):  # type: ignore
        def deco(fn):
            return fn
        return deco


# ---- verbatim extraction (do not edit; source-equality tested) -------------

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
def exhaustive_saving_scalar(lo, hi, resources, continuation, bg, rent, hb, cb, pc,
                   hmax, alpha, oms, beta, es, owner_cost, owner_K, owner):
    """Global maximum for the existing clipped-linear continuation objective.

    Utility is concave on each continuation-grid / renter-cap segment.
    Endpoints plus the interior first-order candidate therefore exhaust the
    segment maximum, even when continuation slopes are globally nonconcave.
    Formula and candidate ordering reproduce the September 5 audited oracle.
    """
    if not (0.0 < alpha < 1.0 and oms < 1.0 and oms != 0.0
            and beta > 0.0 and es > 0.0 and rent > 0.0 and hi >= lo):
        raise ValueError("Unsupported exhaustive-saving objective or interval")
    dc=cb+rent*hb
    cap=rent*(hmax-hb)/(1-alpha)
    hcap=max(hmax-hb,1e-10)
    Kr=(alpha**alpha*((1-alpha)/rent)**(1-alpha))**oms
    points=np.empty(bg.size+3)
    n=0
    points[n]=lo; n+=1
    for x in bg:
        if lo < x < hi:
            points[n]=x; n+=1
    if not owner:
        kink=resources-dc-cap
        if lo < kink < hi:
            points[n]=kink; n+=1
    points[n]=hi; n+=1
    points=np.sort(points[:n])
    best=-1e300; bpbest=lo
    for k in range(n):
        x=points[k]
        if owner:
            value=eval_owner_scalar(x,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
        else:
            value=eval_renter_scalar(x,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
        if value>best:
            best=value;bpbest=x
        if k==n-1 or points[k+1]-x<1e-14:
            continue
        mid=(x+points[k+1])/2
        if mid<=bg[0] or mid>=bg[-1]:
            slope=0.0
        else:
            ix=np.searchsorted(bg,mid)-1
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
                value=eval_owner_scalar(candidate,resources,continuation,bg,owner_cost,cb,pc,owner_K,alpha,oms,beta,es)
            else:
                value=eval_renter_scalar(candidate,resources,continuation,bg,dc,pc,cap,cb,rent,hmax,hcap,Kr,alpha,oms,beta,es)
            if value>best:
                best=value;bpbest=candidate
    return bpbest,best


# ---- proposed optimization (not yet verified) ------------------------------

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
def exhaustive_saving_indexed(lo, hi, resources, continuation, bg, rent, hb, cb, pc,
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
