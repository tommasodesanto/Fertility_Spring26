"""Plot diagnostics for the population-closure mini runs."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    out_dir = ROOT / "benchmarks" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    mini_path = ROOT / "benchmarks" / "mini_parameter_solve_10_fast_2026_05_05.json"
    psi_path = ROOT / "benchmarks" / "scaled_equilibrium_psi_sensitivity_fixed_vo_2026_05_05.json"

    plot_mini_sweep(json.loads(mini_path.read_text()), out_dir / "mini_parameter_solve_10_fast_population.png")
    plot_scaled_psi(json.loads(psi_path.read_text()), out_dir / "scaled_psi_population_housing.png")


def plot_mini_sweep(payload: dict, out_path: Path) -> None:
    cases = [c for c in payload["cases"] if c.get("error") in (None, "None")]
    labels = [short_label(c["label"]) for c in cases]
    tfr = np.array([c["tfr"] for c in cases], dtype=float)
    scale = np.array([c["scale_factor"] for c in cases], dtype=float)
    own = np.array([c["own_rate"] for c in cases], dtype=float)
    err = np.array([c["best_eq_error"] for c in cases], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)

    sc = axes[0].scatter(tfr, scale, c=own, s=70, cmap="viridis", edgecolor="white", linewidth=0.8)
    axes[0].axhline(1.0, color="0.55", linewidth=0.9, linestyle="--")
    axes[0].set_yscale("log")
    axes[0].set_ylim(max(np.min(scale) * 0.45, 1e-5), np.max(scale) * 1.55)
    axes[0].set_xlabel("TFR")
    axes[0].set_ylabel("Implied city scale")
    axes[0].set_title("Population Scale")
    offsets = {
        "baseline": (6, -14),
        "low fertility": (-24, 8),
        "high fertility": (6, 7),
        "patient low child": (6, -3),
        "impatient high child": (6, -10),
        "big 1st child h": (6, 7),
        "flat child h": (6, 0),
        "low owner premium": (6, 8),
        "mobile": (6, 6),
        "sticky": (6, 7),
    }
    for x, y, lab in zip(tfr, scale, labels):
        axes[0].annotate(lab, (x, y), xytext=offsets.get(lab, (4, 4)), textcoords="offset points", fontsize=7)
    cbar = fig.colorbar(sc, ax=axes[0])
    cbar.set_label("Ownership rate")

    axes[1].scatter(tfr, own, color="#3366aa", edgecolor="white", linewidth=0.8, s=70)
    axes[1].set_ylim(-0.035, max(own) * 1.18)
    axes[1].set_xlabel("TFR")
    axes[1].set_ylabel("Ownership rate")
    axes[1].set_title("Tenure Response")
    for x, y, lab in zip(tfr, own, labels):
        axes[1].annotate(lab, (x, y), xytext=offsets.get(lab, (4, 4)), textcoords="offset points", fontsize=7)

    order = np.arange(len(cases))
    colors = ["#4c78a8" if e <= 0.01 else "#e45756" for e in err]
    axes[2].bar(order, err, color=colors)
    axes[2].axhline(0.01, color="0.35", linewidth=0.9, linestyle="--")
    axes[2].set_xticks(order)
    axes[2].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    axes[2].set_ylabel("Best GE error")
    axes[2].set_title("Solver Acceptance")

    fig.suptitle("Direct 10-Spec Mini Sweep: Fast GE Solves", fontsize=13)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_scaled_psi(payload: dict, out_path: Path) -> None:
    cases = payload["cases"]
    labels = [short_label(c["case"]) for c in cases]
    tfr = np.array([c["tfr"] for c in cases], dtype=float)
    scale = np.array([c["implied_total_population"] for c in cases], dtype=float)
    hd = np.array([c["implied_housing_demand"] for c in cases], dtype=float)
    prices = np.array([c["p_eq"] for c in cases], dtype=float)
    err = np.array([c["best_eq_error"] for c in cases], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(tfr, scale, marker="o", color="#2f6f4e")
    ax.axhline(1.0, color="0.55", linewidth=0.9, linestyle="--")
    ax.set_xlabel("TFR")
    ax.set_ylabel("Implied city scale")
    ax.set_title("Fertility Raises Scale")
    for x, y, lab in zip(tfr, scale, labels):
        ax.annotate(lab, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)

    ax = axes[0, 1]
    ax.plot(tfr, hd.sum(axis=1), marker="o", label="total", color="#2f6f4e")
    ax.plot(tfr, hd[:, 0], marker="s", label="periphery", color="#4c78a8")
    ax.plot(tfr, hd[:, 1], marker="^", label="center", color="#f58518")
    ax.set_xlabel("TFR")
    ax.set_ylabel("Aggregate housing demand")
    ax.set_title("Level Housing Demand")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(tfr, prices[:, 0], marker="s", label="periphery", color="#4c78a8")
    ax.plot(tfr, prices[:, 1], marker="^", label="center", color="#f58518")
    ax.set_xlabel("TFR")
    ax.set_ylabel("House price")
    ax.set_title("Equilibrium Prices")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.bar(np.arange(len(cases)), err, color="#4c78a8")
    ax.axhline(0.01, color="0.35", linewidth=0.9, linestyle="--")
    ax.set_xticks(np.arange(len(cases)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Best GE error")
    ax.set_title("Convergence Quality")

    fig.suptitle("Scaled-Housing Fertility Sensitivity", fontsize=13)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def short_label(label: str) -> str:
    replacements = {
        "baseline_x0": "baseline",
        "lower_fertility_value": "low fertility",
        "higher_fertility_value": "high fertility",
        "patient_low_child_value": "patient low child",
        "impatient_high_child_value": "impatient high child",
        "bigger_first_child_housing": "big 1st child h",
        "flatter_child_housing": "flat child h",
        "lower_owner_premium": "low owner premium",
        "mobile_high_loc_noise": "mobile",
        "sticky_low_loc_noise": "sticky",
        "low_psi": "low psi",
        "baseline_psi": "baseline",
        "high_psi": "high psi",
        "very_high_psi": "very high psi",
    }
    return replacements.get(label, label.replace("_", " "))


if __name__ == "__main__":
    main()
