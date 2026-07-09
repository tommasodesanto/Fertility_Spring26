"""Numerical utilities for the DT Python port."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np


def make_grid(P: SimpleNamespace) -> np.ndarray:
    """Four-segment liquid-wealth grid: optional sparse lower tail, a dense
    linear core, an optional linear mid band, and a power-spaced upper buffer.

    Configurable via P (defaults reproduce the frozen July dense grid):
      b_min, b_max, b_grid_power      : overall bounds and upper-buffer curvature
      b_core_lo, b_core_hi, b_mid_hi  : segment boundaries (defaults -5, 7, 15)
      b_frac_low, b_frac_core, b_frac_mid : node fractions for the first three
        segments (defaults 0.08, 0.72, 0.12); the upper buffer gets the rest.
    Set b_frac_low/b_frac_mid to 0 (and b_core_lo=b_min / b_mid_hi=b_core_hi) to
        collapse to a dense core + sparse buffer. The grid always pins nodes exactly
    at 0 and at the legacy scalar b_entry_fixed; external entry distributions
    are scattered onto the grid by the solver rather than pinned here.
    """
    Nb = int(P.Nb)
    core_lo = float(getattr(P, "b_core_lo", -5.0))
    core_hi = float(getattr(P, "b_core_hi", 7.0))
    mid_hi = float(getattr(P, "b_mid_hi", 15.0))
    f_low = float(getattr(P, "b_frac_low", 0.08))
    f_core = float(getattr(P, "b_frac_core", 0.72))
    f_mid = float(getattr(P, "b_frac_mid", 0.12))
    N1 = round(Nb * f_low)
    N2 = round(Nb * f_core)
    N3 = round(Nb * f_mid)
    N4 = Nb - N1 - N2 - N3
    segs = []
    if N1 > 0 and core_lo > P.b_min:
        segs.append(np.linspace(P.b_min, core_lo, N1 + 1)[:-1])
    else:
        N2 += N1
    if N2 > 0:
        segs.append(np.linspace(core_lo, core_hi, N2 + 1)[:-1])
    if N3 > 0 and mid_hi > core_hi:
        segs.append(np.linspace(core_hi, mid_hi, N3 + 1)[:-1])
    else:
        N4 += N3
        mid_hi = core_hi
    if N4 > 0 and P.b_max > mid_hi:
        u4 = np.linspace(0.0, 1.0, N4 + 1)[1:]
        segs.append(mid_hi + (P.b_max - mid_hi) * (u4 ** P.b_grid_power))
    b_grid = np.concatenate(segs).astype(float)
    b_grid = np.unique(b_grid)  # safety: sorted + dedup at segment seams
    b_grid[np.argmin(np.abs(b_grid))] = 0.0
    b_grid[np.argmin(np.abs(b_grid - P.b_entry_fixed))] = P.b_entry_fixed
    return b_grid


def _pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson monotone cubic Hermite node derivatives (== SciPy pchip)."""
    n = x.size
    h = np.diff(x)
    delta = np.diff(y) / h
    d = np.zeros(n)
    if n > 2:
        prod = delta[:-1] * delta[1:]
        w1 = 2.0 * h[1:] + h[:-1]
        w2 = h[1:] + 2.0 * h[:-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            harm = (w1 + w2) / (w1 / delta[:-1] + w2 / delta[1:])
        d[1:-1] = np.where(prod > 0.0, harm, 0.0)

    def _edge(h0, h1, m0, m1):
        dd = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
        if np.sign(dd) != np.sign(m0):
            dd = 0.0
        elif (np.sign(m0) != np.sign(m1)) and (abs(dd) > 3.0 * abs(m0)):
            dd = 3.0 * m0
        return dd

    if n == 2:
        d[0] = d[1] = delta[0]
    else:
        d[0] = _edge(h[0], h[1], delta[0], delta[1])
        d[-1] = _edge(h[-1], h[-2], delta[-1], delta[-2])
    return d


def make_value_interp(bg: np.ndarray, V: np.ndarray, method: str = "linear"):
    """Return a callable f(bq)->interpolated V along axis 0, by `method`.

    `linear` is the lottery interpolation used everywhere else (and what the
    compiled kernels use). `monotone_cubic` is a shape-preserving PCHIP
    (Fritsch-Carlson) Hermite spline -- no overshoot, preserves monotonicity, so
    it cannot create the spurious convex kinks that a natural cubic spline does.
    Intended for the Bellman continuation value only; the forward distribution
    stays linear for mass conservation. Pure NumPy, no SciPy dependency.
    """
    m = str(method or "linear").lower()
    if m in ("linear", "lin"):
        return lambda bq: interp_vector(bg, V, bq)
    if m in ("monotone_cubic", "pchip", "cubic"):
        bg = np.asarray(bg, dtype=float)
        y = np.asarray(V, dtype=float)
        d = _pchip_slopes(bg, y)
        h = np.diff(bg)
        lo, hi = float(bg[0]), float(bg[-1])

        def _f(bq):
            x = np.clip(np.asarray(bq, dtype=float), lo, hi)
            idx = np.clip(np.searchsorted(bg, x, side="right") - 1, 0, bg.size - 2)
            hk = h[idx]
            t = (x - bg[idx]) / hk
            t2 = t * t
            t3 = t2 * t
            h00 = 2.0 * t3 - 3.0 * t2 + 1.0
            h10 = t3 - 2.0 * t2 + t
            h01 = -2.0 * t3 + 3.0 * t2
            h11 = t3 - t2
            return h00 * y[idx] + h10 * hk * d[idx] + h01 * y[idx + 1] + h11 * hk * d[idx + 1]

        return _f
    raise ValueError(f"unknown interp method: {method!r}")


def interp_indices(bg: np.ndarray, bq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    bq_clip = np.clip(np.asarray(bq, dtype=float), bg[0], bg[-1])
    idx = np.searchsorted(bg, bq_clip, side="right") - 1
    idx = np.clip(idx, 0, len(bg) - 2)
    wt = (bq_clip - bg[idx]) / (bg[idx + 1] - bg[idx])
    wt = np.clip(wt, 0.0, 1.0)
    return idx.astype(np.int64), wt


def interp_vector(bg: np.ndarray, V: np.ndarray, bq: np.ndarray) -> np.ndarray:
    idx, wt = interp_indices(bg, bq)
    return (1.0 - wt) * V[idx] + wt * V[idx + 1]


def interp_on_grid(bg: np.ndarray, V: np.ndarray, bq: np.ndarray) -> np.ndarray:
    """Interpolate V along axis 0 using a common bq vector."""
    V = np.asarray(V)
    bq = np.asarray(bq, dtype=float)
    idx, wt = interp_indices(bg, bq)
    Vf = np.reshape(V, (len(bg), -1), order="F")
    out = (1.0 - wt[:, None]) * Vf[idx, :] + wt[:, None] * Vf[idx + 1, :]
    return np.reshape(out, bq.shape + V.shape[1:], order="F")


def interp_cols(bg: np.ndarray, V: np.ndarray, bq: np.ndarray) -> np.ndarray:
    """Column-wise interpolation for V and bq with shape (Nb, n_col)."""
    V = np.asarray(V)
    bq = np.asarray(bq, dtype=float)
    if V.shape != bq.shape:
        raise ValueError(f"interp_cols expects matching shapes, got {V.shape} and {bq.shape}.")
    idx, wt = interp_indices(bg, bq)
    cols = np.arange(V.shape[1])[None, :]
    return (1.0 - wt) * V[idx, cols] + wt * V[idx + 1, cols]


def scatter_redistribute(idx: np.ndarray, wt: np.ndarray, mass: np.ndarray, Nb: int) -> np.ndarray:
    return np.bincount(
        np.concatenate([idx, idx + 1]),
        weights=np.concatenate([(1.0 - wt) * mass, wt * mass]),
        minlength=Nb,
    )[:Nb]


def scatter_redistribute_cols(idx: np.ndarray, wt: np.ndarray, mass: np.ndarray, Nb: int) -> np.ndarray:
    out = np.zeros((Nb, mass.shape[1]))
    for col in range(mass.shape[1]):
        if np.sum(mass[:, col]) < 1e-15:
            continue
        out[:, col] = np.bincount(
            np.concatenate([idx[:, col], idx[:, col] + 1]),
            weights=np.concatenate([(1.0 - wt[:, col]) * mass[:, col], wt[:, col] * mass[:, col]]),
            minlength=Nb,
        )[:Nb]
    return out


def scatter_redistribute_cols_sameidx(idx: np.ndarray, wt: np.ndarray, mass: np.ndarray, Nb: int) -> np.ndarray:
    out = np.zeros((Nb, mass.shape[1]))
    rows = np.concatenate([idx, idx + 1])
    for col in range(mass.shape[1]):
        if np.sum(mass[:, col]) < 1e-15:
            continue
        out[:, col] = np.bincount(
            rows,
            weights=np.concatenate([(1.0 - wt) * mass[:, col], wt * mass[:, col]]),
            minlength=Nb,
        )[:Nb]
    return out


def logsumexp(a: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    m = np.max(a, axis=axis, keepdims=True)
    se = np.sum(np.exp(a - m), axis=axis, keepdims=True)
    ls = m + np.log(se)
    probs = np.exp(a - ls)
    return np.squeeze(ls, axis=axis), probs


# Fortran-order reshape between (Nb, n_parity, n_child_states) and (Nb, nc)
# matches the MATLAB reference, so the flattened column index is
# c = nn + n_parity * cs (parity varies fastest).
def flat_nc(x: np.ndarray, Nb: int, nc: int) -> np.ndarray:
    return np.reshape(x, (Nb, nc), order="F")


def unflat_nc(x: np.ndarray, Nb: int, n_parity: int, n_child_states: int) -> np.ndarray:
    return np.reshape(x, (Nb, n_parity, n_child_states), order="F")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probs: float | np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    probs_arr = np.asarray(probs, dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[keep]
    weights = weights[keep]
    if values.size == 0:
        return np.full(probs_arr.shape, np.nan)
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cw = np.cumsum(weights) / np.sum(weights)
    q = np.zeros(probs_arr.shape)
    for idx, prob in np.ndenumerate(probs_arr):
        pos = np.searchsorted(cw, prob, side="left")
        if pos >= values.size:
            pos = values.size - 1
        q[idx] = values[pos]
    if np.isscalar(probs):
        return float(q)
    return q


def weighted_median_from_cells(value_cells: list[np.ndarray], weight_cells: list[np.ndarray]) -> float:
    if not value_cells:
        return float("nan")
    return float(weighted_quantile(np.concatenate(value_cells), np.concatenate(weight_cells), 0.5))


def dict_to_namespace(d: dict) -> SimpleNamespace:
    return SimpleNamespace(**d)
