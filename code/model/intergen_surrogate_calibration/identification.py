"""Smooth surrogate identification diagnostics.

The June-18 identification ledger's "Immediate Next Check" is a local Jacobian
J_{mp} = d m_m / d theta_p with rank, condition number, and near-collinear
parameter pairs. The existing audit computes it by finite differences over a
near-discrete objective, which is why it is noisy (discrete owner-median-room
moment goes locally flat -> rank deficiency).

This module reads the same Jacobian off the *smooth* emulator analytically and
applies the exact target-normalization used by the finite-difference audit, so
the two can be compared head-to-head:

    J~_{mp} = (d m_m / d theta_p) * s_p / d_m
    s_p   = max(|theta_p|, 1)         (beta uses max(|beta|, 0.05))
    d_m   = max(|model moment_m|, |target_m|, 1)

It also reports an ARD-lengthscale relevance table: a long lengthscale for
parameter p in moment m means the emulated moment barely responds to p -- a
model-derived weak-identification signal independent of the Jacobian SVD.
"""

from __future__ import annotations

import json
import os

import numpy as np

from . import data as D
from .emulator import MomentEmulator

FD_AUDIT_DIR = os.path.join(
    D.REPO_ROOT, "output/model/intergen_sensitivity_jacobian_20260618")


def _param_scale(theta_raw: np.ndarray, theta_params: list[str]) -> np.ndarray:
    s = np.maximum(np.abs(theta_raw), 1.0)
    if "beta" in theta_params:
        bi = theta_params.index("beta")
        s[bi] = max(abs(theta_raw[bi]), 0.05)
    return s


def _moment_scale(model_m: np.ndarray, target_m: np.ndarray) -> np.ndarray:
    return np.maximum.reduce([np.abs(model_m), np.abs(target_m),
                              np.ones_like(model_m)])


def surrogate_jacobian(emu: MomentEmulator, theta_raw: np.ndarray,
                       moment_model: np.ndarray | None = None):
    """Raw and target-normalized surrogate Jacobian at a raw-theta point.

    Returns dict with raw_matrix, scaled_matrix (moments x params), the model
    moments used for scaling, singular values, rank, condition number, and
    column (parameter) cosine-correlation matrix.
    """
    ds = emu.ds
    theta_raw = np.asarray(theta_raw, float)
    x = ds.scale(theta_raw)
    span = ds.hi - ds.lo

    names = ds.moment_names
    nm, npar = len(names), len(ds.theta_params)
    J_raw = np.zeros((nm, npar))
    pred_m = np.zeros(nm)
    for i, name in enumerate(names):
        gp = emu.gps[name]
        g_scaled = gp.predict_grad(x)          # d mean / d x_scaled
        J_raw[i] = g_scaled / span             # chain rule to raw theta
        pred_m[i] = float(gp.predict(x[None, :])[0])

    if moment_model is None:
        moment_model = pred_m
    target_vec = np.array([ds.targets[n] for n in names])
    s_p = _param_scale(theta_raw, ds.theta_params)
    d_m = _moment_scale(moment_model, target_vec)
    J_scaled = J_raw * s_p[None, :] / d_m[:, None]

    sv = np.linalg.svd(J_scaled, compute_uv=False)
    # Match the FD audit's rank tolerance EXACTLY (1e-8 * max(s_max, 1.0)); using
    # the matrix-dimension scaling instead would bias the surrogate rank vs FD.
    tol = 1e-8 * max(float(sv[0]) if sv.size else 0.0, 1.0)
    rank = int(np.sum(sv > tol))
    cond = float(sv[0] / sv[-1]) if sv[-1] > tol else float("inf")

    # Parameter-column collinearity (cosine similarity).
    cols = J_scaled
    norms = np.linalg.norm(cols, axis=0)
    safe = np.where(norms > 0, norms, 1.0)
    unit = cols / safe[None, :]
    colcorr = unit.T @ unit

    return {
        "theta_raw": theta_raw,
        "moment_names": names,
        "param_names": ds.theta_params,
        "raw_matrix": J_raw,
        "scaled_matrix": J_scaled,
        "moment_model": moment_model,
        "moment_pred": pred_m,
        "singular_values": sv,
        "rank": rank,
        "cond": cond,
        "smallest_sv": float(sv[-1]),
        "col_corr": colcorr,
    }


def ard_relevance(emu: MomentEmulator) -> dict:
    """Per-(moment, parameter) ARD relevance = span / lengthscale (scaled units).

    Inputs are in the unit cube, so a lengthscale >> 1 means the moment is
    essentially flat in that parameter over the whole searched range.
    """
    ds = emu.ds
    names = ds.moment_names
    params = ds.theta_params
    ls = np.zeros((len(names), len(params)))
    for i, name in enumerate(names):
        ls[i] = emu.gps[name].ls_
    relevance = 1.0 / ls  # in unit-cube units; larger = more relevant
    return {
        "moment_names": names,
        "param_names": params,
        "lengthscales": ls,
        "relevance": relevance,
        "param_relevance": relevance.mean(axis=0),  # avg over moments
    }


def load_fd_audit():
    """Load finite-difference audit anchors and per-point scaled Jacobians."""
    sp = json.load(open(os.path.join(FD_AUDIT_DIR, "source_points.json")))
    summary = json.load(open(os.path.join(FD_AUDIT_DIR, "audit_summary.json")))
    points = {}
    for p in sp["points"]:
        label = p["point_label"]
        jm_path = os.path.join(FD_AUDIT_DIR, label, "jacobian_matrix.json")
        jm = json.load(open(jm_path)) if os.path.exists(jm_path) else None
        points[label] = {
            "theta": p["theta"],
            "source_moments": p.get("source_moments", {}),
            "source_rank_loss": p.get("source_rank_loss"),
            "jacobian_matrix": jm,
        }
    # rank/cond per point from summary
    pt_summary = {pt.get("point_label", pt.get("label")): pt
                  for pt in summary.get("points", [])}
    return {"points": points, "summary": summary, "point_summary": pt_summary,
            "config": summary.get("config", {})}


def compare_to_fd(emu: MomentEmulator, fd, label: str):
    """Compare the surrogate Jacobian to the FD Jacobian at one anchor point.

    Uses the FD baseline (solved) moments as the moment-scale denominator for
    both, so only the derivative differs between the two scaled matrices.
    """
    ds = emu.ds
    pt = fd["points"][label]
    theta_raw = np.array([float(pt["theta"][p]) for p in ds.theta_params])

    fd_summary = fd["point_summary"].get(label, {})
    base_moments_raw = fd_summary.get("baseline_moments", pt["source_moments"])
    moment_model = np.array([float(base_moments_raw.get(n, np.nan))
                             for n in ds.moment_names])

    sur = surrogate_jacobian(emu, theta_raw, moment_model=moment_model)

    out = {
        "label": label,
        "surrogate_rank": sur["rank"],
        "surrogate_cond": sur["cond"],
        "surrogate_smallest_sv": sur["smallest_sv"],
        "surrogate_scaled_matrix": sur["scaled_matrix"],
        "moment_names": ds.moment_names,
        "param_names": ds.theta_params,
    }
    # FD reported rank/cond (audit_summary uses 'scaled_rank' / 'condition_number').
    out["fd_rank"] = fd_summary.get("scaled_rank", fd_summary.get("rank"))
    out["fd_cond"] = fd_summary.get("condition_number", fd_summary.get("cond"))

    # Entry-wise comparison if FD scaled matrix is available.
    jm = pt["jacobian_matrix"]
    if jm is not None and "scaled_matrix" in jm:
        fd_params = jm.get("parameters", ds.theta_params)
        fd_moments = jm.get("target_moments", ds.moment_names)
        fd_scaled = np.array(jm["scaled_matrix"], dtype=float)
        # Reorder FD matrix to our moment/param order.
        mi = [fd_moments.index(n) for n in ds.moment_names]
        pj = [fd_params.index(p) for p in ds.theta_params]
        fd_aligned = fd_scaled[np.ix_(mi, pj)]
        sur_s = sur["scaled_matrix"]
        a, b = fd_aligned.ravel(), sur_s.ravel()
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() > 3:
            pear = float(np.corrcoef(a[m], b[m])[0, 1])
            sign_agree = float(np.mean(np.sign(a[m]) == np.sign(b[m])))
        else:
            pear, sign_agree = np.nan, np.nan
        out["fd_vs_surrogate_pearson"] = pear
        out["fd_vs_surrogate_sign_agreement"] = sign_agree
        out["fd_scaled_matrix_aligned"] = fd_aligned
    return out
