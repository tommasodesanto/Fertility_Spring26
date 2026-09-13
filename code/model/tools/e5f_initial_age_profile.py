#!/usr/bin/env python3
"""Experimental initial fertility-stock age profiles and six added SMM rows.

``score_extra`` consumes the existing ``uniform_birth_time`` initial-fertility
observer packet.  It does not alter the retained twelve-row target system.  Its
diagonal weights use a synthetic standard deviation equal to five percent of
each target and are explicitly provisional rather than empirical uncertainty.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_PACKET = (
    ROOT
    / "output/model/e5f_matched_pf_20260909a/design_research/"
    "fertility_contract/age_profile/age_profile_candidates.json"
)
TOP_BIN_REPRESENTATIVE = 3.602359422009
AGE_WINDOWS = ((20, 24), (25, 29), (30, 34), (35, 39), (40, 44))
MEAN_SCORE_WINDOWS = ((25, 29), (30, 34), (35, 39), (40, 44))
CHILDLESS_SCORE_WINDOWS = ((25, 29), (35, 39))
SYNTHETIC_RELATIVE_SCALE = 0.05
TOLERANCE = 2.0e-12


def _window_label(lower: int, upper: int) -> str:
    return f"{int(lower)}_{int(upper)}"


def _array(value: Any, name: str, *, ndim: int) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim != ndim or result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a nonempty finite {ndim}-dimensional array")
    if np.any(result < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return result


def _model_profile(
    pre_parity: np.ndarray,
    post_parity: np.ndarray,
    age_starts: np.ndarray,
    *,
    lower: int,
    upper: int,
    cell_width: float,
) -> dict[str, float]:
    interval_right = float(upper) + 1.0
    overlap_left = np.maximum(age_starts, float(lower))
    overlap_right = np.minimum(age_starts + cell_width, interval_right)
    lengths = np.maximum(overlap_right - overlap_left, 0.0)
    if not math.isclose(float(np.sum(lengths)), interval_right - lower, abs_tol=TOLERANCE):
        raise ValueError(f"model ages do not cover the full [{lower},{interval_right:g}) window")
    overlap_weights = lengths / cell_width
    selected = lengths > 0.0
    post_shares = np.zeros_like(age_starts)
    post_shares[selected] = (
        (overlap_left[selected] + overlap_right[selected]) / 2.0
        - age_starts[selected]
    ) / cell_width
    projected = (1.0 - post_shares[:, None]) * pre_parity + post_shares[:, None] * post_parity
    window_mass = np.sum(overlap_weights[:, None] * projected, axis=0)
    total = float(np.sum(window_mass))
    if total <= 0.0:
        raise ValueError(f"model fertility-stock mass is zero for ages {lower}-{upper}")
    shares = window_mass / total
    if not math.isclose(float(np.sum(shares)), 1.0, abs_tol=TOLERANCE):
        raise RuntimeError("model fertility-stock shares do not sum to one")
    coded_mean = float(
        np.dot(shares, np.array([0.0, 1.0, 2.0, TOP_BIN_REPRESENTATIVE]))
    )
    return {
        "age_lower": int(lower),
        "age_upper": int(upper),
        "share_0": float(shares[0]),
        "share_1": float(shares[1]),
        "share_2": float(shares[2]),
        "share_3plus": float(shares[3]),
        "mean_model_coded_CEB": coded_mean,
        "window_model_mass": total,
        "overlap_weights": overlap_weights.tolist(),
        "post_fertility_interpolation_shares": [
            float(value) if selected[index] else None
            for index, value in enumerate(post_shares)
        ],
    }


def _model_profiles(packet: Mapping[str, Any]) -> dict[str, dict[str, float]]:
    metadata = packet.get("metadata", {})
    if metadata.get("age_projection") != "uniform_birth_time":
        raise ValueError("score_extra requires the existing uniform_birth_time observer packet")
    accounting = packet.get("accounting", {})
    pre = _array(accounting.get("pre_parity_mass_by_age"), "pre parity mass", ndim=2)
    post = _array(accounting.get("post_parity_mass_by_age"), "post parity mass", ndim=2)
    ages = _array(accounting.get("age_cell_start"), "age starts", ndim=1)
    if pre.shape != post.shape or pre.shape != (ages.size, 4):
        raise ValueError("observer packet must contain matching age-by-0/1/2/3+ arrays")
    if ages.size < 2:
        raise ValueError("at least two model age cells are required")
    widths = np.diff(ages)
    cell_width = float(widths[0])
    if cell_width <= 0.0 or not np.allclose(widths, cell_width, rtol=0.0, atol=TOLERANCE):
        raise ValueError("observer age cells must have a positive constant width")
    age_mass_error = float(np.max(np.abs(np.sum(pre, axis=1) - np.sum(post, axis=1))))
    if age_mass_error > 2.0e-10:
        raise RuntimeError(f"pre/post fertility age mass differs by {age_mass_error:.3e}")
    profiles = {
        _window_label(lower, upper): _model_profile(
            pre,
            post,
            ages,
            lower=lower,
            upper=upper,
            cell_width=cell_width,
        )
        for lower, upper in AGE_WINDOWS
    }
    original = packet.get("parity_shares_40_44")
    if not isinstance(original, Mapping):
        raise ValueError("observer packet is missing its original parity_shares_40_44 gate")
    generalized = profiles["40_44"]
    for old_key, new_key in (("0", "share_0"), ("1", "share_1"), ("2", "share_2"), ("3plus", "share_3plus")):
        if abs(float(original[old_key]) - float(generalized[new_key])) > TOLERANCE:
            raise RuntimeError("generalized age observer does not reproduce the original 40-44 observer")
    return profiles


def _empirical_profiles(path: Path = EMPIRICAL_PACKET) -> dict[str, dict[str, float]]:
    with path.open(encoding="utf-8") as handle:
        source = json.load(handle)
    if not math.isclose(
        float(source.get("top_bin_representative", math.nan)),
        TOP_BIN_REPRESENTATIVE,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("empirical packet uses a different 3+ representative")
    rows = {
        _window_label(int(row["age_lower"]), int(row["age_upper"])): row
        for row in source.get("rows", [])
        if row.get("window") == "pooled"
    }
    expected = {_window_label(*window) for window in AGE_WINDOWS}
    if set(rows) != expected:
        raise ValueError("empirical packet does not contain exactly the five pooled age profiles")
    for label, row in rows.items():
        shares = np.array([row["share_0"], row["share_1"], row["share_2"], row["share_3plus"]])
        if np.any(~np.isfinite(shares)) or np.any(shares < 0.0) or not math.isclose(
            float(np.sum(shares)), 1.0, abs_tol=TOLERANCE
        ):
            raise ValueError(f"empirical shares are invalid for {label}")
        reconstructed = float(
            np.dot(shares, np.array([0.0, 1.0, 2.0, TOP_BIN_REPRESENTATIVE]))
        )
        if abs(reconstructed - float(row["mean_model_coded_CEB"])) > TOLERANCE:
            raise ValueError(f"empirical model-coded mean is inconsistent for {label}")
    return rows


def score_extra(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Return six experimental score rows and all five age-profile shares."""

    model = _model_profiles(packet)
    target = _empirical_profiles()
    specifications = [
        (f"mean_model_coded_CEB_{_window_label(*window)}", window, "mean_model_coded_CEB")
        for window in MEAN_SCORE_WINDOWS
    ] + [
        (f"childless_rate_{_window_label(*window)}", window, "share_0")
        for window in CHILDLESS_SCORE_WINDOWS
    ]
    rows: list[dict[str, Any]] = []
    for name, window, field in specifications:
        label = _window_label(*window)
        target_value = float(target[label][field])
        model_value = float(model[label][field])
        synthetic_scale = SYNTHETIC_RELATIVE_SCALE * abs(target_value)
        if synthetic_scale <= 0.0:
            raise ValueError(f"synthetic five-percent scale is zero for {name}")
        weight = 1.0 / synthetic_scale**2
        gap = model_value - target_value
        rows.append(
            {
                "moment": name,
                "target": target_value,
                "model": model_value,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap**2,
                "synthetic_scale": synthetic_scale,
                "weight_status": "provisional_diagonal_synthetic_5pct_of_target_not_standard_error",
            }
        )
    profiles: list[dict[str, Any]] = []
    for lower, upper in AGE_WINDOWS:
        label = _window_label(lower, upper)
        row: dict[str, Any] = {"age_lower": lower, "age_upper": upper}
        for field in ("share_0", "share_1", "share_2", "share_3plus", "mean_model_coded_CEB"):
            row[f"target_{field}"] = float(target[label][field])
            row[f"model_{field}"] = float(model[label][field])
        row["model_window_mass"] = float(model[label]["window_model_mass"])
        row["overlap_weights"] = model[label]["overlap_weights"]
        row["post_fertility_interpolation_shares"] = model[label][
            "post_fertility_interpolation_shares"
        ]
        profiles.append(row)
    return {
        "rows": rows,
        "age_profiles": profiles,
        "extra_loss": float(sum(float(row["loss_contribution"]) for row in rows)),
        "metadata": {
            "status": "experimental augmented-initial pilot only",
            "existing_twelve_rows_changed": False,
            "empirical_source": str(EMPIRICAL_PACKET),
            "target_mean_definition": "pooled CPS 0/1/2/3+ shares with fixed 3+ representative",
            "top_bin_representative": TOP_BIN_REPRESENTATIVE,
            "age_projection": "uniform_birth_time pre/post stock interpolation",
            "weight_definition": "diagonal inverse variance using synthetic scale 0.05*abs(target)",
            "weight_is_empirical_standard_error": False,
            "overlap_covariance_treatment": "ignored in provisional experimental objective",
        },
    }
