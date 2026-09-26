"""Isolated, fixed-width housing-floor sensitivity for the compiled Markov path.

Install after the purchase/accounting adapter.  Frozen sources are never edited.
Only parent intratemporal housing services change; childless columns call the
installed native kernels unchanged.  Numerical work belongs on Torch.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path

import numpy as np
from numba import njit, prange


@njit(cache=False)
def soft_services(h, hb, softness):
    """Return g and g', using expm1 near zero to avoid cancellation."""
    if hb <= 0.0:
        return h, 1.0
    delta = softness * hb
    t = h / delta
    z0 = -1.0 / softness
    z = t + z0
    if t < 30.0:
        g = delta * math.log1p(math.exp(z0) * math.expm1(t) /
                             (1.0 + math.exp(z0)))
    else:
        sp = max(z, 0.0) + math.log1p(math.exp(-abs(z)))
        g = delta * (sp - math.log1p(math.exp(z0)))
    if z >= 0.0:
        gp = 1.0 / (1.0 + math.exp(-z))
    else:
        ez = math.exp(z)
        gp = ez / (1.0 + ez)
    return g, gp


@njit(cache=False)
def allocation(x, rent, hb, hmax, alpha, softness):
    """Exact budget, numerical monotone intratemporal FOC; x excludes c_bar."""
    if x <= 1e-10:
        return 0.0, 0.0
    if hb <= 0.0:
        h = min(hmax, (1.0 - alpha) * x / rent)
        return x - rent * h, h
    gcap, dcap = soft_services(hmax, hb, softness)
    ratio = alpha / (1.0 - alpha)
    cap_x = rent * (hmax + ratio * gcap / dcap)
    if x >= cap_x:
        return x - rent * hmax, hmax
    lo = 0.0
    hi = min(hmax, x / rent)
    h = 0.5 * (lo + hi)
    delta = softness * hb
    for _ in range(60):
        g, gp = soft_services(h, hb, softness)
        residual = h + ratio * g / gp - x / rent
        if abs(residual) <= 1e-12 * max(x / rent, 1e-12):
            break
        if residual > 0.0:
            hi = h
        else:
            lo = h
        derivative = 1.0 + ratio * (1.0 - g * (1.0 - gp) / (delta * gp))
        candidate = h - residual / derivative
        if not lo < candidate < hi:
            candidate = 0.5 * (lo + hi)
        h = candidate
    return x - rent * h, h


@njit(cache=False)
def linear_continuation(bg, values, x):
    if x <= bg[0]:
        return values[0]
    if x >= bg[-1]:
        return values[-1]
    ix = np.searchsorted(bg, x) - 1
    w = (x - bg[ix]) / (bg[ix + 1] - bg[ix])
    return (1.0 - w) * values[ix] + w * values[ix + 1]


@njit(cache=False)
def renter_value(bp, resources, values, bg, rent, hb, cb, pc, hmax,
                 alpha, beta, es, softness):
    x = resources - cb - bp
    if x <= 1e-10:
        return -1e10
    c, h = allocation(x, rent, hb, hmax, alpha, softness)
    g, _ = soft_services(h, hb, softness)
    return -es / (c ** alpha * g ** (1.0 - alpha)) + pc + beta * linear_continuation(bg, values, bp)


@njit(cache=False)
def demand_from_marginal(lam, rent, hb, hmax, alpha, es, softness):
    """Invert u_c=lambda without nesting an allocation solve inside a root."""
    gcap, dcap = soft_services(hmax, hb, softness)
    ccap = alpha * rent * gcap / ((1.0 - alpha) * dcap)
    mucap = es * alpha * ccap ** (-alpha - 1.0) * gcap ** (alpha - 1.0)
    if lam <= mucap:
        c = (es * alpha * gcap ** (alpha - 1.0) / lam) ** (1.0 / (alpha + 1.0))
        return c, hmax
    lo = 0.0
    hi = hmax
    for _ in range(65):
        h = 0.5 * (lo + hi)
        g, gp = soft_services(h, hb, softness)
        c = alpha * rent * g / ((1.0 - alpha) * gp)
        mu = es * alpha * c ** (-alpha - 1.0) * g ** (alpha - 1.0)
        if mu > lam:
            lo = h
        else:
            hi = h
    h = 0.5 * (lo + hi)
    g, gp = soft_services(h, hb, softness)
    return alpha * rent * g / ((1.0 - alpha) * gp), h


@njit(cache=False)
def exhaustive_soft(lo, hi, resources, values, bg, rent, hb, cb, pc,
                    hmax, alpha, beta, es, softness, demands=None):
    """Exhaust all continuation segments using concavity, including cap kink."""
    points = np.empty(bg.size + 3)
    n = 1
    points[0] = lo
    for x in bg:
        if lo < x < hi:
            points[n] = x
            n += 1
    gc, dc = soft_services(hmax, hb, softness)
    cap_x = rent * (hmax + alpha * gc / ((1.0 - alpha) * dc))
    kink = resources - cb - cap_x
    if lo < kink < hi:
        points[n] = kink
        n += 1
    points[n] = hi
    n += 1
    points = np.sort(points[:n])
    best = -1e300
    bpbest = lo
    for k in range(n):
        x = points[k]
        value = renter_value(x, resources, values, bg, rent, hb, cb, pc,
                             hmax, alpha, beta, es, softness)
        if value > best:
            best, bpbest = value, x
        if k == n - 1 or points[k + 1] - x < 1e-14:
            continue
        mid = 0.5 * (x + points[k + 1])
        if mid <= bg[0] or mid >= bg[-1]:
            continue
        ix = np.searchsorted(bg, mid) - 1
        slope = (values[ix + 1] - values[ix]) / (bg[ix + 1] - bg[ix])
        if slope <= 0.0:
            continue
        if demands is None:
            c, h = demand_from_marginal(beta * slope, rent, hb, hmax, alpha, es, softness)
            spending = c + rent * h
        else:
            spending = demands[ix]
        candidate = resources - cb - spending
        if x < candidate < points[k + 1]:
            value = renter_value(candidate, resources, values, bg, rent, hb, cb,
                                 pc, hmax, alpha, beta, es, softness)
            if value > best:
                best, bpbest = value, candidate
    return bpbest, best


@njit(cache=False, parallel=True)
def parent_renter_kernel(Rv, Rvt, Vc, previous, has_prev, bg, cbv, hbv,
                         psiv, gbv, alphav, esv, rent, hmax, cmin, cb0, hb0,
                         alpha, oms, beta, s_next, D_next, a1, a2, tol,
                         exhaustive, softness):
    nb, nc = Vc.shape
    vo = np.empty((nb, nc))
    bpout = np.empty((nb, nc))
    co = np.empty((nb, nc))
    ho = np.empty((nb, nc))
    for k in prange(nc):
        cb, hb, pc, al, es = cbv[k], hbv[k], psiv[k], alphav[k], esv[k]
        # Every wealth row shares this continuation column. Invert each positive
        # segment slope once, then shift optimal spending by row resources.
        demands = np.empty(bg.size - 1)
        if exhaustive:
            for ix in range(bg.size - 1):
                slope = (Vc[ix + 1, k] - Vc[ix, k]) / (bg[ix + 1] - bg[ix])
                demands[ix] = np.nan
                if slope > 0.0:
                    cc, hh = demand_from_marginal(beta * slope, rent, hb, hmax, al, es, softness)
                    demands[ix] = cc + rent * hh
        for b in range(nb):
            resources = Rv[b]
            if gbv[k] > 0.0:
                resources += min(max(gbv[k] - Rvt[b], 0.0), gbv[k])
            floor = min(s_next * min(bg[b], 0.0), -D_next)
            lo = max(floor, bg[0])
            hi = max(resources - cb - 1e-6, lo)
            if has_prev:
                lo = max(lo, previous[b, k] - 2.0)
                hi = min(hi, previous[b, k] + 2.0)
                lo = max(lo, floor, bg[0])
                hi = max(hi, lo)
            if exhaustive:
                bp, v = exhaustive_soft(lo, hi, resources, Vc[:, k], bg, rent,
                                        hb, cb, pc, hmax, al, beta, es, softness, demands)
            else:
                # Identical native golden-section update and tie ordering.
                d = hi - lo
                x1, x2 = lo + a1 * d, lo + a2 * d
                f1 = renter_value(x1, resources, Vc[:, k], bg, rent, hb, cb, pc, hmax, al, beta, es, softness)
                f2 = renter_value(x2, resources, Vc[:, k], bg, rent, hb, cb, pc, hmax, al, beta, es, softness)
                d *= a1 * a2
                while d > tol:
                    if f2 >= f1:
                        xe = min(x2 + d, hi)
                        fe = renter_value(xe, resources, Vc[:, k], bg, rent, hb, cb, pc, hmax, al, beta, es, softness)
                        x1, f1, x2, f2 = x2, f2, xe, fe
                    else:
                        xe = max(x1 - d, lo)
                        fe = renter_value(xe, resources, Vc[:, k], bg, rent, hb, cb, pc, hmax, al, beta, es, softness)
                        x2, f2, x1, f1 = x1, f1, xe, fe
                    d *= a2
                if f2 >= f1:
                    bp, v = x2, f2
                else:
                    bp, v = x1, f1
            vo[b, k], bpout[b, k] = v, bp
            x = resources - cb - bp
            if x <= 1e-10:
                co[b, k], ho[b, k] = cb0 + cmin, hb0 + 0.01
            else:
                c, h = allocation(x, rent, hb, hmax, al, softness)
                co[b, k], ho[b, k] = cb + c, h
    return vo, bpout, co, ho


def install(model, output_dir, *, softness=0.1):
    """Install into only the model's full Markov Bellman; return provenance."""
    if softness != 0.1:
        raise ValueError("Only the contracted fixed width delta=0.1*hP is supported")
    if getattr(model, "_soft_housing_installed", False):
        raise ValueError("Soft housing adapter must be installed exactly once")
    native_renter = model.full_renter_block_kernel
    native_owner = model.full_owner_block_kernel

    def renter(*args):
        if len(args) != 26:
            raise ValueError("Unexpected native renter signature")
        if args[18] != -1.0 or (args[25] and args[4]):
            raise ValueError("Soft renter requires sigma=2 and exhaustive full interval")
        hb = np.asarray(args[7])
        if np.any(hb < 0.0) or np.any(~np.isfinite(hb)):
            raise ValueError("Invalid parent housing requirement")
        if not np.any(hb > 0.0):
            return native_renter(*args)
        output = tuple(np.empty_like(args[2]) for _ in range(4))
        for mask, soft in ((hb == 0.0, False), (hb > 0.0, True)):
            if not np.any(mask):
                continue
            subset = list(args)
            for ix in (2, 3):
                subset[ix] = np.ascontiguousarray(args[ix][:, mask])
            for ix in range(6, 12):
                subset[ix] = np.ascontiguousarray(args[ix][mask])
            result = (parent_renter_kernel(*subset, softness) if soft else native_renter(*subset))
            for dest, source in zip(output, result):
                dest[:, mask] = source
        return output

    def owner(*args):
        if len(args) != 28:
            raise ValueError("Unexpected native owner signature")
        hb = np.asarray(args[7])
        if not np.any(hb > 0.0):
            return native_owner(*args)
        house, scale = float(args[14]), float(args[15])
        if scale != 1.0:
            raise ValueError("Contract requires owner_h_bar_scale=1")
        revised = list(args)
        effective = hb.copy()
        for k in np.flatnonzero(hb > 0.0):
            g, _ = soft_services(house, float(hb[k]), softness)
            if not g > 1e-10:
                raise ValueError("Native owner housing clamp would bind under soft services")
            effective[k] = (house - g) / scale
            reconstructed = house - scale * effective[k]
            if not np.isclose(reconstructed, g, rtol=1e-10, atol=0.0):
                raise ValueError("Owner effective floor loses numerical precision")
        revised[7] = effective
        return native_owner(*revised)

    source = inspect.getsource(model.solve_bellman_full_markov_income)
    if source.count("full_renter_block_kernel(") != 1 or source.count("full_owner_block_kernel(") != 1:
        raise ValueError("Unexpected full Markov kernel call sites")
    guard = '''    if (not NUMBA_AVAILABLE or not bool(getattr(P, "use_full_kernel", True))
            or str(getattr(P, "interp_method", "linear")) != "linear"
            or float(P.sigma) != 2.0
            or str(getattr(P, "preference_spec", "")) != "eqscale"
            or not bool(getattr(P, "child_room_floor", False))):
        raise ValueError("Unsupported path for isolated soft housing experiment")
'''
    anchor = "    t0 = time.perf_counter()\n"
    if source.count(anchor) != 1:
        raise ValueError("Unexpected Markov entry anchor")
    generated = source.replace(anchor, anchor + guard).replace(
        "full_renter_block_kernel(", "_soft_renter_block(").replace(
        "full_owner_block_kernel(", "_soft_owner_block(")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    generated_path = out / "soft_housing_bellman.generated.py"
    if generated_path.exists():
        raise FileExistsError(generated_path)
    generated_path.write_text(generated)
    model._soft_renter_block = renter
    model._soft_owner_block = owner
    exec(compile(generated, str(generated_path), "exec"), model.__dict__)
    model._soft_housing_installed = True
    receipt = {"softness": softness, "sigma": 2.0,
               "service_definition": "delta*(softplus((h-hP)/delta)-softplus(-hP/delta)); delta=.1*hP",
               "original_bellman_sha256": hashlib.sha256(source.encode()).hexdigest(),
               "generated_bellman_sha256": hashlib.sha256(generated.encode()).hexdigest(),
               "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "generated_source": str(generated_path),
               "childless": "native call-through", "policy_evaluation_path": "unsupported"}
    (out / "soft_housing_adapter.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def run_self_tests():
    """Torch-only independent scipy checks; no full model solve or import."""
    from scipy.optimize import minimize_scalar
    count = 0
    maximum_relative_error = 0.0
    for hb in (0.3, 1.2, 2.7):
        for rent in (0.08, 0.4):
            for x in (0.001, 0.05, 0.5, 2.0, 8.0):
                alpha, hmax, es = 0.72, 5.0, 1.234
                c, h = allocation(x, rent, hb, hmax, alpha, 0.1)
                g, gp = soft_services(h, hb, 0.1)
                assert abs(c + rent * h - x) < 1e-12
                assert c > 0.0 and 0.0 < h <= hmax
                # Independent scalar objective uses numpy logaddexp, not g().
                def objective(hh):
                    gg = 0.1 * hb * (np.logaddexp(0.0, (hh - hb) / (0.1 * hb)) - np.logaddexp(0.0, -10.0))
                    return es / ((x - rent * hh) ** alpha * gg ** (1.0 - alpha))
                upper = min(hmax, x / rent * (1.0 - 1e-12))
                ref = minimize_scalar(objective, bounds=(upper * 1e-12, upper), method="bounded",
                                      options={"xatol": 1e-13})
                exact = es / (c ** alpha * g ** (1.0 - alpha))
                reference = min(ref.fun, objective(upper))
                error = abs(exact - reference) / max(1.0, abs(exact))
                assert error < 2e-7, (hb, rent, x, error)
                if h < hmax * (1.0 - 1e-10):
                    assert abs(alpha / c - (1.0 - alpha) * gp / (rent * g)) / (alpha / c) < 1e-10
                maximum_relative_error = max(maximum_relative_error, error)
                count += 1
    bg = np.array([-1.0, 0.0, 0.7, 2.0, 5.0])
    for values in (np.array([-9., -5., -4., -1., 0.]), np.array([-9., -8., -1., -4., -2.])):
        for resources in (0.04, 0.8, 3.0, 8.0):
            lo, hi = 0.0, resources - 1e-6
            bp, value = exhaustive_soft(lo, hi, resources, values, bg, 0.3, 1.2,
                                        0.0, 0.8, 5.0, 0.72, 0.95, 1.234, 0.1)
            edges = sorted(set([lo, hi] + [float(x) for x in bg if lo < x < hi]))
            best = -np.inf
            for left, right in zip(edges[:-1], edges[1:]):
                objective = lambda xx: -renter_value(xx, resources, values, bg, 0.3, 1.2,
                                                     0.0, 0.8, 5.0, 0.72, 0.95, 1.234, 0.1)
                fit = minimize_scalar(objective, bounds=(left, right), method="bounded", options={"xatol": 1e-12})
                best = max(best, -fit.fun, -objective(left), -objective(right))
            error = abs(value - best) / max(abs(best), 1.0)
            assert error < 2e-8, (resources, bp, value, best)
            maximum_relative_error = max(maximum_relative_error, error)
            count += 1
    # Exercise the compiled parent block, cached segment inversions, borrowing
    # rollover, transfers, and exact reported allocations against scalar solves.
    resources = np.array([0.3, 0.7, 1.3, 2.0, 4.0])
    values = np.column_stack((np.array([-9., -5., -4., -1., 0.]),
                              np.array([-9., -8., -1., -4., -2.])))
    hb = np.array([0.3, 1.2])
    cb = np.zeros(2)
    pc = np.array([0.8, 0.9])
    grants = np.array([0.0, 0.5])
    al = np.full(2, 0.72)
    es = np.full(2, 1.234)
    result = parent_renter_kernel(resources, resources, values, np.zeros_like(values),
                                  0, bg, cb, hb, pc, grants, al, es, 0.3, 5.0,
                                  1e-4, 0.0, 0.0, 0.72, -1.0, 0.95, 0.6, 0.2,
                                  (3.0 - math.sqrt(5.0)) / 2.0,
                                  (math.sqrt(5.0) - 1.0) / 2.0, 1e-3, 1, 0.1)
    for k in range(2):
        for b in range(5):
            rr = resources[b] + min(max(grants[k] - resources[b], 0.0), grants[k])
            lo = max(min(0.6 * min(bg[b], 0.0), -0.2), bg[0])
            hi = max(rr - 1e-6, lo)
            bp, v = exhaustive_soft(lo, hi, rr, values[:, k], bg, 0.3, hb[k],
                                    0.0, pc[k], 5.0, 0.72, 0.95, 1.234, 0.1)
            assert abs(result[0][b, k] - v) < 1e-10
            assert abs(result[1][b, k] - bp) < 1e-10
            assert abs(result[2][b, k] + 0.3 * result[3][b, k] + bp - rr) < 1e-12
            count += 1
    return {"status": "passed", "checks": count, "maximum_relative_error": maximum_relative_error}
