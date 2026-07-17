#!/usr/bin/env python3
"""Plot one intergen policy-function state slice from a saved solution cache."""

from __future__ import annotations

import argparse
import csv
import math
import pickle
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RUN_DIR = ROOT / "output" / "model" / "intergen_model_run_current"
DEFAULT_CACHE = DEFAULT_RUN_DIR / "solution_cache.pkl"
DEFAULT_OUTDIR = DEFAULT_RUN_DIR / "policy_slices"


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    payload = load_cache(args.cache.resolve())
    sol = payload["baseline"]["sol"]
    P = payload["baseline"]["P"]

    b_grid = np.asarray(sol.b_grid, dtype=float)
    ages = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    z_grid = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float)

    j = nearest_index(ages, args.age)
    zz = nearest_index(z_grid, args.z)
    ten = tenure_index(args.tenure)
    i = int(args.location)
    nn = int(args.parity)
    cs = int(args.child_state)

    rows = extract_rows(sol, P, b_grid, j, zz, ten, i, nn, cs, args.b_min, args.b_max)
    if not rows:
        raise ValueError("Selected state slice has no grid points after b-window filtering.")

    stem = args.name or (
        f"policy_slice_age{ages[j]:.0f}_z{z_grid[zz]:.2f}_"
        f"ten{ten}_loc{i}_n{nn}_cs{cs}"
    ).replace(".", "p")
    csv_path = outdir / f"{stem}.csv"
    png_path = outdir / f"{stem}.png"
    md_path = outdir / f"{stem}_readout.md"

    write_csv(csv_path, rows)
    plot_rows(rows, png_path, ages[j], z_grid[zz], ten, i, nn, cs)
    write_readout(md_path, csv_path, png_path, rows, ages[j], z_grid[zz], ten, i, nn, cs)

    print(f"Wrote {png_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--age", type=float, default=30.0)
    parser.add_argument("--z", type=float, default=1.0, help="Income state value; nearest grid value is used.")
    parser.add_argument("--tenure", default="renter", help="'renter' or integer tenure index; renter is 0.")
    parser.add_argument("--location", type=int, default=0)
    parser.add_argument("--parity", type=int, default=0)
    parser.add_argument("--child-state", type=int, default=0)
    parser.add_argument("--b-min", type=float, default=-2.0)
    parser.add_argument("--b-max", type=float, default=8.0)
    parser.add_argument("--name", default="", help="Optional output filename stem.")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} does not look like a solution cache.")
    return payload


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def tenure_index(value: str) -> int:
    text = str(value).strip().lower()
    if text in {"r", "rent", "renter"}:
        return 0
    if text in {"own", "owner"}:
        return 1
    return int(text)


def extract_rows(
    sol: Any,
    P: Any,
    b_grid: np.ndarray,
    j: int,
    zz: int,
    ten: int,
    i: int,
    nn: int,
    cs: int,
    b_min: float,
    b_max: float,
) -> list[dict[str, float]]:
    c_pol = np.asarray(sol.c_pol, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    hR_pol = np.asarray(sol.hR_pol, dtype=float)
    g = np.asarray(sol.g, dtype=float)
    tenure_probs = np.asarray(sol.tenure_probs, dtype=float)
    fert_probs = np.asarray(sol.fert_probs, dtype=float)
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float)

    rows: list[dict[str, float]] = []
    keep = (b_grid >= float(b_min)) & (b_grid <= float(b_max))
    for ib in np.where(keep)[0]:
        own_prob = float(np.sum(tenure_probs[ib, ten, i, j, zz, nn, cs, 1:]))
        if ten == 0:
            housing = float(hR_pol[ib, ten, i, j, zz, nn, cs])
        elif 0 < ten <= len(h_own):
            housing = float(h_own[ten - 1])
        else:
            housing = math.nan
        fert = float(np.sum(np.arange(fert_probs.shape[-1]) * fert_probs[ib, ten, i, j, zz, :]))
        rows.append(
            {
                "b": float(b_grid[ib]),
                "mass": float(g[ib, ten, i, j, zz, nn, cs]),
                "consumption": float(c_pol[ib, ten, i, j, zz, nn, cs]),
                "next_liquid_wealth": float(bp_pol[ib, ten, i, j, zz, nn, cs]),
                "housing_services": housing,
                "ownership_probability": own_prob,
                "expected_completed_children": fert,
            }
        )
    for prev, cur in zip(rows, rows[1:]):
        cur["delta_consumption"] = cur["consumption"] - prev["consumption"]
        cur["delta_next_liquid_wealth"] = cur["next_liquid_wealth"] - prev["next_liquid_wealth"]
        cur["delta_housing_services"] = cur["housing_services"] - prev["housing_services"]
    if rows:
        rows[0]["delta_consumption"] = math.nan
        rows[0]["delta_next_liquid_wealth"] = math.nan
        rows[0]["delta_housing_services"] = math.nan
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: list[dict[str, float]], path: Path, age: float, z: float, ten: int, i: int, nn: int, cs: int) -> None:
    b = arr(rows, "b")
    mass = arr(rows, "mass")
    occupied = mass > 1e-10
    sizes = 18.0 + 260.0 * mass / max(float(np.max(mass)), 1e-300)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), constrained_layout=True)
    fig.suptitle(f"Policy slice: age {age:.0f}, z={z:.2f}, tenure={ten}, loc={i}, n={nn}, cs={cs}")
    panels = [
        ("consumption", "consumption", "#2563eb"),
        ("next_liquid_wealth", "next liquid wealth b'", "#7c3aed"),
        ("housing_services", "housing services", "#059669"),
        ("expected_completed_children", "expected completed children", "#dc2626"),
    ]
    for ax, (key, title, color) in zip(axes.ravel(), panels):
        y = arr(rows, key)
        ax.plot(b, y, color=color, lw=2.0)
        if np.any(occupied):
            ax.scatter(b[occupied], y[occupied], s=sizes[occupied], color=color, edgecolor="white", zorder=3)
        if key == "expected_completed_children":
            own = arr(rows, "ownership_probability")
            ax.plot(b, own, color="#111827", lw=1.6, label="Pr(own)")
            ax.legend(frameon=False, loc="best")
        if key == "next_liquid_wealth":
            ax.plot(b, b, color="0.55", lw=1.0, ls="--")
        ax.axvline(0.0, color="0.45", lw=1.0)
        ax.set_title(title)
        ax.set_xlabel("liquid wealth before choice b")
        ax.grid(alpha=0.25)
    fig.text(0.01, 0.01, "Dots mark occupied grid points; dot area is proportional to mass in this exact state.", color="0.35")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_readout(
    path: Path,
    csv_path: Path,
    png_path: Path,
    rows: list[dict[str, float]],
    age: float,
    z: float,
    ten: int,
    i: int,
    nn: int,
    cs: int,
) -> None:
    occupied = [r for r in rows if r["mass"] > 1e-10]
    top = sorted(occupied, key=lambda r: r["mass"], reverse=True)[:8]
    dips = sorted(
        [r for r in rows if math.isfinite(r["delta_consumption"]) and r["delta_consumption"] < -0.05],
        key=lambda r: r["delta_consumption"],
    )[:8]
    lines = [
        f"# Policy slice: age {age:.0f}, z={z:.2f}, tenure={ten}, loc={i}, n={nn}, cs={cs}",
        "",
        f"- Figure: `{png_path}`",
        f"- CSV: `{csv_path}`",
        f"- Occupied mass: {sum(r['mass'] for r in occupied):.6g}",
    ]
    if occupied:
        lines.append(f"- Occupied b range: {min(r['b'] for r in occupied):.3f} to {max(r['b'] for r in occupied):.3f}")
        lines.append(f"- Max ownership probability on occupied points: {max(r['ownership_probability'] for r in occupied):.6g}")
    lines.extend(["", "## Top occupied points", "", "| b | mass | c | b_next | h | Pr(own) | E[n] |", "|---:|---:|---:|---:|---:|---:|---:|"])
    for r in top:
        lines.append(row_line(r))
    lines.extend(["", "## Largest consumption dips", "", "| b | delta_c | delta_b_next | delta_h | mass |", "|---:|---:|---:|---:|---:|"])
    for r in dips:
        lines.append(
            f"| {r['b']:.3f} | {r['delta_consumption']:.3f} | "
            f"{r['delta_next_liquid_wealth']:.3f} | {r['delta_housing_services']:.3f} | {r['mass']:.6g} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def row_line(r: dict[str, float]) -> str:
    return (
        f"| {r['b']:.3f} | {r['mass']:.6g} | {r['consumption']:.3f} | "
        f"{r['next_liquid_wealth']:.3f} | {r['housing_services']:.3f} | "
        f"{r['ownership_probability']:.3g} | {r['expected_completed_children']:.3f} |"
    )


def arr(rows: list[dict[str, float]], key: str) -> np.ndarray:
    return np.asarray([r[key] for r in rows], dtype=float)


if __name__ == "__main__":
    main()
