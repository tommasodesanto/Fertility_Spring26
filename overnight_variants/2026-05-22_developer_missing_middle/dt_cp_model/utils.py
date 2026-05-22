"""Numerical utilities for the DT Python port."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np


def make_grid(P: SimpleNamespace) -> np.ndarray:
    Nb = int(P.Nb)
    N1 = round(Nb * 0.15)
    N2 = round(Nb * 0.45)
    N3 = round(Nb * 0.15)
    N4 = Nb - N1 - N2 - N3
    s1 = np.linspace(P.b_min, -3.0, N1 + 1)[:-1]
    s2 = np.linspace(-3.0, 6.0, N2 + 1)[:-1]
    s3 = np.linspace(6.0, 20.0, N3 + 1)[:-1]
    u4 = np.linspace(0.0, 1.0, N4 + 1)[1:]
    s4 = 20.0 + (P.b_max - 20.0) * (u4 ** P.b_grid_power)
    b_grid = np.concatenate([s1, s2, s3, s4]).astype(float)
    b_grid[np.argmin(np.abs(b_grid))] = 0.0
    b_grid[np.argmin(np.abs(b_grid - P.b_entry_fixed))] = P.b_entry_fixed
    return b_grid


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
