#!/usr/bin/env python3
"""Plot deterministic/modal policy-function overlays from a saved solution cache."""

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

from tools.build_intergen_mechanics_packet import (
    deterministic_branch_wealth,
    deterministic_fertility_choice,
    interp_policy_scalar,
    owner_asset_price_vector,
    precompute_shared,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PACKET = ROOT / "output/model/intergen_postaudit_sample3h_review_20260627/best_globalde_task1_packet"
DEFAULT_CACHE = DEFAULT_PACKET / "solution_cache.pkl"
DEFAULT_OUTDIR = DEFAULT_PACKET / "policy_overlays"


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    payload = load_cache(args.cache.resolve())
    baseline = payload["baseline"]
    sol = baseline["sol"]
    P = baseline["P"]

    rows = build_overlay_rows(sol, P, args)
    if not rows:
        raise ValueError("No policy rows generated for the requested slice.")

    stem = args.name or f"deterministic_policy_overlay_age{requested_age_label(rows)}"
    csv_path = outdir / f"{stem}.csv"
    png_path = outdir / f"{stem}.png"
    md_path = outdir / f"{stem}_readout.md"
    write_csv(csv_path, rows)
    plot_overlay(rows, png_path)
    write_readout(md_path, rows, csv_path, png_path, payload)

    print(f"Wrote {png_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--age", type=float, default=30.0)
    parser.add_argument("--z-low-index", type=int, default=0)
    parser.add_argument("--z-high-index", type=int, default=-1)
    parser.add_argument(
        "--owner-h",
        type=float,
        default=4.0,
        help="Current owner rung in room units; nearest H_own rung is used.",
    )
    parser.add_argument("--market", type=int, default=0)
    parser.add_argument("--parity", type=int, default=0)
    parser.add_argument("--child-state", type=int, default=0)
    parser.add_argument("--b-min", type=float, default=-2.5)
    parser.add_argument("--b-max", type=float, default=12.0)
    parser.add_argument("--name", default="")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} does not look like a mechanics-packet solution cache.")
    return payload


def build_overlay_rows(sol: Any, P: Any, args: argparse.Namespace) -> list[dict[str, Any]]:
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)

    j = int(np.argmin(np.abs(age_grid - float(args.age))))
    i = int(args.market)
    nn = int(args.parity)
    cs = int(args.child_state)
    owner_tenure = 1 + int(np.argmin(np.abs(h_own - float(args.owner_h))))
    z_low = normalize_index(int(args.z_low_index), len(z_grid))
    z_high = normalize_index(int(args.z_high_index), len(z_grid))

    state_specs = [
        ("low-income renter", 0, z_low),
        ("high-income renter", 0, z_high),
        (f"low-income owner H={h_own[owner_tenure - 1]:g}", owner_tenure, z_low),
        (f"high-income owner H={h_own[owner_tenure - 1]:g}", owner_tenure, z_high),
    ]
    c_pol = np.asarray(sol.c_pol, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    hR_pol = np.asarray(sol.hR_pol, dtype=float)
    tenure_choice = np.asarray(sol.tenure_choice)
    tenure_probs = getattr(sol, "tenure_probs", None)
    tenure_probs = None if tenure_probs is None else np.asarray(tenure_probs, dtype=float)
    fert_probs = np.asarray(sol.fert_probs, dtype=float)
    V = np.asarray(getattr(sol, "V", np.empty_like(c_pol)), dtype=float)
    price = owner_asset_price_vector(sol, P)
    shared = precompute_shared(P, b_grid)

    rows: list[dict[str, Any]] = []
    keep = (b_grid >= float(args.b_min)) & (b_grid <= float(args.b_max))
    for label, current_tenure, zz in state_specs:
        for bb in np.where(keep)[0]:
            if V.ndim >= 7 and V[bb, current_tenure, i, j, zz, nn, cs] <= -1e9:
                continue
            target_tenure = int(tenure_choice[bb, current_tenure, i, j, zz, nn, cs])
            branch_wealth = deterministic_branch_wealth(
                b=float(b_grid[bb]),
                origin_tenure=current_tenure,
                target_tenure=target_tenure,
                market=i,
                parity=nn,
                child_state=cs,
                P=P,
                price=price,
                shared=shared,
            )
            c_val = interp_policy_scalar(b_grid, c_pol[:, target_tenure, i, j, zz, nn, cs], branch_wealth)
            bp_val = interp_policy_scalar(b_grid, bp_pol[:, target_tenure, i, j, zz, nn, cs], branch_wealth)
            if target_tenure <= 0:
                h_val = interp_policy_scalar(b_grid, hR_pol[:, 0, i, j, zz, nn, cs], branch_wealth)
                target_h = 0.0
                target_label = "renter"
            else:
                target_h = float(h_own[target_tenure - 1])
                h_val = target_h
                target_label = f"owner H={target_h:g}"
            if tenure_probs is None:
                owner_prob = float(target_tenure > 0)
            else:
                owner_prob = float(np.sum(tenure_probs[bb, current_tenure, i, j, zz, nn, cs, 1:]))
            fert_choice, expected_children = deterministic_fertility_choice(
                fert_probs, bb, current_tenure, i, j, zz, nn
            )
            rows.append(
                {
                    "state": label,
                    "age": float(age_grid[j]),
                    "age_index": int(j),
                    "z": float(z_grid[zz]),
                    "z_index": int(zz),
                    "current_tenure_index": int(current_tenure),
                    "current_tenure": "renter" if current_tenure <= 0 else f"owner H={h_own[current_tenure - 1]:g}",
                    "liquid_wealth": float(b_grid[bb]),
                    "target_tenure_index": target_tenure,
                    "target_tenure": target_label,
                    "target_owner_rooms": target_h,
                    "branch_liquid_wealth_after_transaction": float(branch_wealth),
                    "consumption": float(c_val),
                    "next_liquid_wealth": float(bp_val),
                    "housing_services": float(h_val),
                    "owner_probability": owner_prob,
                    "deterministic_fertility_choice": float(fert_choice),
                    "expected_children_probability_weighted": float(expected_children),
                }
            )
    return rows


def normalize_index(index: int, n: int) -> int:
    if n <= 0:
        raise ValueError("No income states found.")
    if index < 0:
        index = n + index
    return int(np.clip(index, 0, n - 1))


def requested_age_label(rows: list[dict[str, Any]]) -> str:
    age = float(rows[0]["age"])
    return str(int(round(age)))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_overlay(rows: list[dict[str, Any]], path: Path) -> None:
    states = list(dict.fromkeys(str(r["state"]) for r in rows))
    colors = {
        "low-income renter": "#2563eb",
        "high-income renter": "#1d4ed8",
    }
    palette = ["#2563eb", "#1d4ed8", "#dc2626", "#991b1b", "#059669", "#7c3aed"]
    for idx, state in enumerate(states):
        colors.setdefault(state, palette[idx % len(palette)])

    age = float(rows[0]["age"])
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    fig.suptitle(f"Deterministic policy slices at age {age:g}: fixed current states, modal tenure choice")
    panels = [
        ("target_tenure_index", "chosen tenure/rung index"),
        ("consumption", "consumption"),
        ("next_liquid_wealth", "next liquid wealth b'"),
        ("housing_services", "housing services"),
        ("owner_probability", "logit Pr(owner)"),
        ("expected_children_probability_weighted", "expected children"),
    ]
    for ax, (key, title) in zip(axes.ravel(), panels):
        for state in states:
            state_rows = [r for r in rows if str(r["state"]) == state]
            x = np.asarray([float(r["liquid_wealth"]) for r in state_rows], dtype=float)
            y = np.asarray([float(r[key]) for r in state_rows], dtype=float)
            ax.plot(x, y, lw=2.0, label=state, color=colors[state])
        if key == "next_liquid_wealth":
            xmin, xmax = ax.get_xlim()
            lo, hi = min(xmin, xmax), max(xmin, xmax)
            ax.plot([lo, hi], [lo, hi], color="0.55", lw=1.0, ls="--")
        ax.axvline(0.0, color="0.45", lw=0.9)
        ax.set_title(title)
        ax.set_xlabel("current liquid wealth b")
        ax.grid(alpha=0.25)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_readout(
    path: Path,
    rows: list[dict[str, Any]],
    csv_path: Path,
    png_path: Path,
    payload: dict[str, Any],
) -> None:
    states = list(dict.fromkeys(str(r["state"]) for r in rows))
    lines = [
        "# Deterministic Policy Overlay",
        "",
        f"- Figure: `{png_path}`",
        f"- CSV: `{csv_path}`",
        f"- Source target set: `{payload.get('target_set', '')}`",
        f"- Source theta: `{payload.get('theta', {})}`",
        "",
        "This plot fixes the current state and plots the deterministic/modal policy over current liquid wealth.",
        "Consumption and saving are evaluated at branch wealth after the chosen tenure transaction.",
        "The owner-probability panel reports the logit tenure probability for reference only.",
        "",
        "## States",
        "",
    ]
    for state in states:
        sample = next(r for r in rows if str(r["state"]) == state)
        lines.append(
            f"- `{state}`: age `{sample['age']:.0f}`, z `{sample['z']:.3g}`, "
            f"current tenure `{sample['current_tenure']}`"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
