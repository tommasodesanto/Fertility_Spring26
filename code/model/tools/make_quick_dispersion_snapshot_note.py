#!/usr/bin/env python3
"""Build a fast PDF snapshot from a direct-calibration record.

This avoids re-solving the model.  It is meant for live calibration triage:
saved record moments, saved parameters, target comparisons, and compact plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.direct_calibration import build_direct_calibration_setup  # noqa: E402


HARD_TARGET_ROWS = [
    ("tfr", "TFR"),
    ("childless_rate", "Childless"),
    ("mean_age_first_birth", "Age first birth"),
    ("tfr_gradient", "TFR gradient"),
    ("own_rate", "Ownership"),
    ("own_gradient", "Ownership gradient"),
    ("own_family_gap", "Family ownership gap"),
    ("housing_increment_0to1", "Rooms 0->1"),
    ("housing_increment_1to2", "Rooms 1->2"),
    ("young_liquid_wealth_to_income", "Young wealth/income"),
    ("center_share_nonparents", "Center nonparents"),
    ("center_share_newparents", "Center new parents"),
    ("migration_rate", "Migration rate"),
    ("old_age_own_rate", "Old ownership"),
    ("old_age_parent_childless_gap", "Old parent gap"),
    ("inv_pop_share_C", "Center share"),
    ("inv_rent_ratio_C_over_P", "Rent C/P"),
]

ROOM_TARGET_ROWS = [
    ("owner25_45_rooms_le6_share", "Owner <=6"),
    ("owner25_45_rooms_7to8_share", "Owner 7-8"),
    ("owner25_45_rooms_ge9_share", "Owner >=9"),
]

PARAMETER_ORDER = [
    "beta",
    "b_entry_fixed",
    "psi_child",
    "h_bar_jump",
    "h_bar_n",
    "c_bar_n",
    "kappa_fert",
    "chi",
    "kappa_loc",
    "mu_move",
    "theta0",
    "theta_n",
    "h_bar_0",
    "E_C",
    "r_bar_C",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--label", default="dispersion snapshot")
    parser.add_argument("--hR-max", type=float, default=8.0)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--note-tex", type=Path, required=True)
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    record = json.loads(args.record.read_text())
    setup = build_direct_calibration_setup(
        "benchmark",
        geo_weight=250.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=args.hR_max,
        alpha_cons=record.get("alpha_cons"),
        owner_h_bar_scale=record.get("owner_h_bar_scale"),
        owner_size_cost=record.get("owner_size_cost"),
        owner_size_cost_ref=record.get("owner_size_cost_ref"),
        owner_size_cost_power=record.get("owner_size_cost_power"),
        tenure_choice_kappa=record.get("tenure_choice_kappa"),
        weight_overrides=record.get("weight_overrides") or None,
        extra_targets=record.get("extra_targets") or None,
        H_own=record.get("H_own") or None,
    )
    moments = record.get("moments") or {}
    rows = comparison_rows(setup.targets, setup.weights, moments)
    write_csv(args.outdir / "target_comparison.csv", rows)
    write_csv(args.outdir / "parameters.csv", parameter_rows(record, args.hR_max))

    plot_key_fit(rows, args.outdir / "key_fit.png")
    plot_room_bins(rows, args.outdir / "owner_room_bins.png")
    plot_abs_errors(rows, args.outdir / "absolute_errors.png")
    write_note(args.note_tex, args.outdir, args.record, args.label, record, rows, args.hR_max)
    if args.compile:
        compile_latex(args.note_tex)
    print(f"Wrote {args.note_tex}")


def comparison_rows(targets: dict, weights: dict, moments: dict) -> list[dict]:
    rows = []
    for key, label in HARD_TARGET_ROWS + ROOM_TARGET_ROWS:
        if key not in targets:
            continue
        target = float(targets[key])
        model = float(moments.get(key, math.nan))
        rows.append(
            {
                "moment": key,
                "label": label,
                "target": target,
                "weight": float(weights.get(key, math.nan)),
                "model": model,
                "diff": model - target,
            }
        )
    return rows


def parameter_rows(record: dict, hR_max: float) -> list[dict]:
    params = record.get("parameters") or {}
    rows = [{"parameter": key, "value": float(params[key])} for key in PARAMETER_ORDER if key in params]
    controls = {
        "alpha_cons": record.get("alpha_cons"),
        "tenure_choice_kappa": record.get("tenure_choice_kappa"),
        "owner_h_bar_scale": record.get("owner_h_bar_scale"),
        "owner_size_cost": record.get("owner_size_cost"),
        "owner_size_cost_ref": record.get("owner_size_cost_ref"),
        "owner_size_cost_power": record.get("owner_size_cost_power"),
        "hR_max": hR_max,
    }
    prices = record.get("p_eq") or []
    if len(prices) >= 1:
        controls["p_P"] = prices[0]
    if len(prices) >= 2:
        controls["p_C"] = prices[1]
    for key, value in controls.items():
        if value is not None:
            rows.append({"parameter": key, "value": float(value)})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_key_fit(rows: list[dict], path: Path) -> None:
    keys = [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_gradient",
        "own_family_gap",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "inv_pop_share_C",
        "inv_rent_ratio_C_over_P",
    ]
    panel = [row for row in rows if row["moment"] in keys]
    labels = [row["label"] for row in panel]
    x = np.arange(len(panel))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.bar(x - width / 2, [row["target"] for row in panel], width, label="Target", color="#c9cdd3")
    ax.bar(x + width / 2, [row["model"] for row in panel], width, label="Model", color="#2d5f9a")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_title("Key Moment Fit")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_room_bins(rows: list[dict], path: Path) -> None:
    panel = [row for row in rows if row["moment"] in {key for key, _ in ROOM_TARGET_ROWS}]
    labels = [row["label"] for row in panel]
    x = np.arange(len(panel))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.bar(x - width / 2, [row["target"] for row in panel], width, label="ACS target", color="#c9cdd3")
    ax.bar(x + width / 2, [row["model"] for row in panel], width, label="Model", color="#7b3f8c")
    ax.set_ylim(0, 0.65)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Share")
    ax.set_title("Prime-Age Owner Room Distribution")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_abs_errors(rows: list[dict], path: Path) -> None:
    panel = sorted(rows, key=lambda row: abs(row["diff"]), reverse=True)[:12]
    labels = [row["label"] for row in panel]
    vals = [row["diff"] for row in panel]
    colors = ["#9d2f2f" if v > 0 else "#2e6d45" for v in vals]
    y = np.arange(len(panel))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(y, vals, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_title("Largest Raw Moment Gaps")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def write_note(
    note_tex: Path,
    outdir: Path,
    record_path: Path,
    label: str,
    record: dict,
    rows: list[dict],
    hR_max: float,
) -> None:
    note_tex.parent.mkdir(parents=True, exist_ok=True)
    rel = Path("../") / outdir.resolve().relative_to(ROOT)
    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=0.65in]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{graphicx}",
        r"\usepackage{float}",
        r"\title{Quick Dispersion-Wave Fit Snapshot}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Record}",
        r"\begin{itemize}",
        rf"\item Label: \textbf{{{tex(label)}}}.",
        rf"\item Source record: \texttt{{{tex(str(record_path))}}}.",
        rf"\item Run tag: \texttt{{{tex(record.get('run_tag', ''))}}}; job {record.get('job_id')}; eval {record.get('eval_id')}; loss {float(record.get('loss')):.3f}.",
        rf"\item $h_R^{{\max}}={hR_max:.1f}$; tenure/product logit $\kappa^T={float(record.get('tenure_choice_kappa') or math.nan):.3f}$; owner size cost {float(record.get('owner_size_cost') or 0.0):.4f}.",
        rf"\item Owner ladder: [{tex(', '.join(f'{float(x):.2g}' for x in (record.get('H_own') or [])))}].",
        r"\end{itemize}",
        r"\section*{Parameters}",
        table_from_rows(parameter_rows(record, hR_max), ["parameter", "value"]),
        r"\section*{Targets And Model Moments}",
        table_from_rows(rows, ["label", "target", "weight", "model", "diff"]),
        r"\section*{Plots}",
        figure_line(rel / "key_fit.png", "Key moment fit"),
        figure_line(rel / "owner_room_bins.png", "Room-bin discipline"),
        figure_line(rel / "absolute_errors.png", "Largest raw moment gaps"),
        r"\end{document}",
    ]
    note_tex.write_text("\n".join(lines) + "\n")


def table_from_rows(rows: list[dict], cols: list[str]) -> str:
    out = [
        r"\begin{table}[H]\centering\scriptsize",
        r"\begin{tabular}{l" + "r" * (len(cols) - 1) + r"}",
        r"\toprule",
        " & ".join(tex(col) for col in cols) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        vals = []
        for col in cols:
            value = row[col]
            if isinstance(value, (int, float)):
                vals.append("NA" if not math.isfinite(float(value)) else f"{float(value):.4g}")
            else:
                vals.append(tex(value))
        out.append(" & ".join(vals) + r" \\")
    out.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(out)


def figure_line(path: Path, caption: str) -> str:
    return "\n".join(
        [
            r"\begin{figure}[H]\centering",
            rf"\includegraphics[width=0.88\textwidth]{{{path.as_posix()}}}",
            rf"\caption{{{tex(caption)}}}",
            r"\end{figure}",
        ]
    )


def tex(value) -> str:
    return (
        str(value)
        .replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("$", r"\$")
        .replace("#", r"\#")
        .replace("_", r"\_")
        .replace("{", r"\{")
        .replace("}", r"\}")
    )


def compile_latex(note_tex: Path) -> None:
    for _ in range(2):
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", note_tex.name],
            cwd=note_tex.parent,
            check=True,
        )


if __name__ == "__main__":
    main()
