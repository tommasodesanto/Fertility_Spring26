#!/usr/bin/env python3
"""Plot tenure-conditional consumption policies from a cached solution pickle."""

from __future__ import annotations

import os
import pickle
import sys
import tempfile
from pathlib import Path

import numpy as np

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "fertility_mplconfig"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


AUDIT_ROOT = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/"
    "Fertility_Spring26_fable_size_mapping_audit_20260701"
)
SOURCE_PICKLE = (
    AUDIT_ROOT
    / "output/model/diag_packet_nb120_best_20260707/solution_cache.pkl"
)
OUTPUT_FIGURE = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/"
    "latex/July_26_slides_eq_consumption_by_tenure.png"
)

REQUESTED_AGES = (30.0, 42.0)
ROOMS_TARGET = 4.0
LOCATION_INDEX = 0
CHILD_COUNT_INDEX = 0
CHILD_STATE_INDEX = 0


def add_unpickle_paths() -> None:
    """Add audit-copy model paths needed only for unpickling object graphs."""
    for rel in ("model_diag", "model_active_ref"):
        path = str(AUDIT_ROOT / rel)
        if path not in sys.path:
            sys.path.insert(0, path)


def add_numpy_pickle_compat() -> None:
    """Allow NumPy-2 pickles to load under the local NumPy-1.x environment."""
    sys.modules.setdefault("numpy._core", np.core)
    sys.modules.setdefault("numpy._core.multiarray", np.core.multiarray)
    sys.modules.setdefault("numpy._core.numeric", np.core.numeric)
    sys.modules.setdefault("numpy._core._multiarray_umath", np.core._multiarray_umath)


def load_payload() -> dict:
    add_numpy_pickle_compat()
    add_unpickle_paths()
    with SOURCE_PICKLE.open("rb") as handle:
        return pickle.load(handle)


def weighted_quantile_on_grid(
    grid: np.ndarray, weights: np.ndarray, probs: tuple[float, float]
) -> np.ndarray:
    grid = np.asarray(grid, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    if grid.shape != weights.shape:
        raise ValueError("grid and weights must have matching one-dimensional shapes")

    valid = np.isfinite(grid) & np.isfinite(weights) & (weights >= 0.0)
    grid = grid[valid]
    weights = weights[valid]
    if weights.sum() <= 0.0:
        raise ValueError("cannot compute weighted quantiles with zero mass")

    order = np.argsort(grid)
    grid = grid[order]
    weights = weights[order]
    cum_weights = np.cumsum(weights)
    targets = np.asarray(probs, dtype=float) * cum_weights[-1]
    return np.interp(targets, cum_weights, grid)


def visible_mask(grid: np.ndarray, xlo: float, xhi: float) -> np.ndarray:
    mask = (grid >= xlo) & (grid <= xhi)
    if int(mask.sum()) >= 3:
        return mask

    left = max(0, int(np.searchsorted(grid, xlo, side="left")) - 1)
    right = min(len(grid) - 1, int(np.searchsorted(grid, xhi, side="right")))
    while right - left + 1 < 3 and (left > 0 or right < len(grid) - 1):
        if left > 0:
            left -= 1
        if right < len(grid) - 1:
            right += 1
    expanded = np.zeros_like(grid, dtype=bool)
    expanded[left : right + 1] = True
    return expanded


def monotonicity_report(x: np.ndarray, y: np.ndarray) -> dict[str, float | int | bool]:
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    diffs = np.diff(y)
    tol = 1e-8
    negative = diffs < -tol
    second = np.diff(y, n=2)
    return {
        "n_points": int(y.size),
        "monotone_non_decreasing": bool(np.all(~negative)),
        "negative_segments": int(np.sum(negative)),
        "min_diff": float(np.min(diffs)) if diffs.size else float("nan"),
        "largest_drop": float(max(0.0, -np.min(diffs))) if diffs.size else 0.0,
        "max_abs_second_diff": float(np.max(np.abs(second))) if second.size else 0.0,
        "c_min": float(np.min(y)) if y.size else float("nan"),
        "c_max": float(np.max(y)) if y.size else float("nan"),
    }


def main() -> None:
    payload = load_payload()
    baseline = payload["baseline"]
    sol = baseline["sol"]
    P = baseline["P"]

    if not hasattr(sol, "c_pol"):
        raise AttributeError("sol has no direct consumption policy attribute c_pol")
    if not hasattr(sol, "g"):
        raise AttributeError("sol has no stationary distribution attribute g")

    b_grid = np.asarray(sol.b_grid, dtype=float)
    c_pol = np.asarray(sol.c_pol, dtype=float)
    g = np.asarray(sol.g, dtype=float)
    H_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    z_grid = np.asarray(P.z_grid, dtype=float).reshape(-1)
    age_grid = float(P.age_start) + np.arange(int(P.J)) * float(P.da)

    if c_pol.shape != g.shape:
        raise ValueError(f"c_pol shape {c_pol.shape} does not match g shape {g.shape}")
    if c_pol.ndim != 7:
        raise ValueError(f"expected 7D policies, found c_pol.ndim={c_pol.ndim}")

    owner_matches = np.flatnonzero(np.isclose(H_own, ROOMS_TARGET, rtol=0.0, atol=1e-10))
    if owner_matches.size != 1:
        raise ValueError(f"expected one owner rung with {ROOMS_TARGET:g} rooms, found {owner_matches}")
    owner_h_index = int(owner_matches[0])
    owner_tenure_index = owner_h_index + 1
    renter_tenure_index = 0

    age_indices = [int(np.argmin(np.abs(age_grid - age))) for age in REQUESTED_AGES]
    used_ages = [float(age_grid[idx]) for idx in age_indices]
    income_index = int(len(z_grid) // 2)

    tenure_specs = [
        ("Renter", renter_tenure_index),
        (f"Owner, {ROOMS_TARGET:g} rooms", owner_tenure_index),
    ]

    x_ranges: dict[int, tuple[float, float]] = {}
    for _, tenure_idx in tenure_specs:
        tenure_mass = g[:, tenure_idx, :, :, :, :, :]
        tenure_weights = tenure_mass.sum(axis=tuple(range(1, tenure_mass.ndim)))
        q01, q99 = weighted_quantile_on_grid(b_grid, tenure_weights, (0.01, 0.99))
        if not np.isfinite(q01) or not np.isfinite(q99) or q99 <= q01:
            raise ValueError(f"invalid tenure {tenure_idx} x-range: {(q01, q99)}")
        x_ranges[tenure_idx] = (float(q01), float(q99))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150, sharey=True)
    colors = ["#1f77b4", "#d62728"]
    line_reports: list[tuple[str, float, dict[str, float | int | bool]]] = []

    for ax, (panel_label, tenure_idx) in zip(axes, tenure_specs):
        xlo, xhi = x_ranges[tenure_idx]
        mask = visible_mask(b_grid, xlo, xhi)
        for color, age_idx, age in zip(colors, age_indices, used_ages):
            y = c_pol[
                :,
                tenure_idx,
                LOCATION_INDEX,
                age_idx,
                income_index,
                CHILD_COUNT_INDEX,
                CHILD_STATE_INDEX,
            ]
            ax.plot(b_grid[mask], y[mask], color=color, linewidth=2.0, label=f"Age {age:g}")
            report = monotonicity_report(b_grid[mask], y[mask])
            line_reports.append((panel_label, age, report))

        ax.set_xlim(xlo, xhi)
        ax.set_title(panel_label)
        ax.set_xlabel("Liquid wealth \\$")
        ax.set_ylabel("Consumption \\$")
        ax.grid(True, color="#d9d9d9", linewidth=0.8, alpha=0.8)
        ax.legend(frameon=False)

    fig.suptitle("Consumption policies by tenure", fontsize=14)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    OUTPUT_FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FIGURE)
    plt.close(fig)

    print(f"source_pickle={SOURCE_PICKLE}")
    print(f"output_figure={OUTPUT_FIGURE}")
    print(f"consumption_source=sol.c_pol")
    print(f"c_pol_shape={c_pol.shape}")
    print(f"stationary_distribution_source=sol.g")
    print(f"H_own={H_own.tolist()}")
    print(
        "owner_4_rooms_tenure_index="
        f"{owner_tenure_index} H_own_index={owner_h_index} rooms={H_own[owner_h_index]:g}"
    )
    print(f"age_grid={age_grid.tolist()}")
    print(f"requested_ages={list(REQUESTED_AGES)} used_ages={used_ages} age_indices={age_indices}")
    print(
        f"income_index={income_index} z_value={z_grid[income_index]:g} "
        f"z_grid={z_grid.tolist()}"
    )
    print(f"state_indices=n:{CHILD_COUNT_INDEX} s:{CHILD_STATE_INDEX} location:{LOCATION_INDEX}")
    for label, tenure_idx in tenure_specs:
        xlo, xhi = x_ranges[tenure_idx]
        print(f"x_range_{label.replace(' ', '_').replace(',', '').lower()}=({xlo:.10g}, {xhi:.10g})")
    print(f"panel_count={len(axes)}")
    print(f"plotted_line_count={len(line_reports)}")
    for label, age, report in line_reports:
        print(
            "line_check "
            f"tenure={label} age={age:g} "
            f"monotone_non_decreasing={report['monotone_non_decreasing']} "
            f"negative_segments={report['negative_segments']} "
            f"min_diff={report['min_diff']:.10g} "
            f"largest_drop={report['largest_drop']:.10g} "
            f"max_abs_second_diff={report['max_abs_second_diff']:.10g} "
            f"c_range=({report['c_min']:.10g}, {report['c_max']:.10g}) "
            f"n_points={report['n_points']}"
        )


if __name__ == "__main__":
    main()
