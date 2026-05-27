#!/usr/bin/env python3
"""Generate the slide-model figure bundle for direct-calibration best records."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]

import sys

if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.objective import extract_moments  # noqa: E402
from dt_cp_model.parameters import asdict  # noqa: E402
from dt_cp_model.solver import housing_demand_normalizer, run_model_cp_dt  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402
from dt_cp_model.utils import make_grid  # noqa: E402
from export_distributional_discipline import age_index, cell_arrays, current_children, room_bin_masks  # noqa: E402
from plot_lifecycle import lifecycle_arrays, plot_lifecycle  # noqa: E402


TARGET_ROWS = [
    ("tfr", "TFR"),
    ("childless_rate", "Childless"),
    ("mean_age_first_birth", "Age first birth"),
    ("own_rate", "Ownership"),
    ("own_gradient", "Ownership gradient"),
    ("own_family_gap", "Family ownership gap"),
    ("own_lifecycle_slope", "Lifecycle ownership slope"),
    ("prime_childless_renter_median_rooms", "Childless renter rooms"),
    ("prime_childless_owner_median_rooms", "Childless owner rooms"),
    ("housing_increment_0to1", "Rooms, 0 to 1 child"),
    ("housing_increment_1to2", "Rooms, 1 to 2 children"),
    ("old_age_own_rate", "Old-age ownership"),
    ("old_age_parent_childless_gap", "Old parent-childless gap"),
    ("inv_pop_share_C", "Center share"),
    ("inv_rent_ratio_C_over_P", "Rent ratio C/P"),
]

ROOM_BINS = ["<=4", "5", "6", "7-8", "9-10", "11+"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Build slide-style fit figures and a comparison note.")
    parser.add_argument("--case", action="append", required=True, help="Case spec LABEL=best.json,hR_max")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--note-tex", type=Path, required=True)
    parser.add_argument("--compile", action="store_true")
    parser.add_argument("--reuse-figures", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases = []
    for spec in args.case:
        label, record_path, hr_max = parse_case(spec)
        case = build_case(label, record_path, hr_max, args.outdir / label, reuse_figures=args.reuse_figures)
        cases.append(case)

    write_note(args.note_tex, args.outdir, cases)
    if args.compile:
        compile_latex(args.note_tex)
    print(f"Wrote note {args.note_tex}")


def parse_case(spec: str) -> tuple[str, Path, float]:
    if "=" not in spec:
        raise ValueError(f"Case must be LABEL=path,hR_max, got {spec!r}")
    label, rhs = spec.split("=", 1)
    path_str, hr_str = rhs.rsplit(",", 1)
    return label, Path(path_str), float(hr_str)


def build_case(label: str, record_path: Path, hr_max: float, outdir: Path, *, reuse_figures: bool = False) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    record = json.loads(record_path.read_text())
    setup = build_direct_calibration_setup(
        "benchmark",
        geo_weight=100.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=hr_max,
    )
    moments = dict(record.get("moments", {}))
    figure_names = [
        "benchmark_lifecycle_compare.png",
        "market_equilibrium.png",
        "housing_choice_heatmap.png",
        "housing_size_fit.png",
    ]
    have_figures = all((outdir / name).exists() for name in figure_names)
    if not (reuse_figures and have_figures):
        sol, P, p_eq = solve_record(record, setup)
        b_grid = make_grid(P)
        fresh = extract_moments(sol, P, p_eq, 0.0, 0.0, 0.0, True)
        moments["fresh_own_rate"] = float(getattr(fresh, "own_rate", np.nan))
        moments["fresh_tfr"] = float(getattr(fresh, "tfr", np.nan))

        plot_lifecycle(
            lifecycle_arrays(sol, P, b_grid),
            P,
            sol,
            outdir / "benchmark_lifecycle_compare.png",
            f"Lifecycle Fit: {label}",
        )
        plot_market_equilibrium(sol, P, p_eq, setup, outdir / "market_equilibrium.png", label)
        plot_housing_choice_heatmap(sol, P, b_grid, outdir / "housing_choice_heatmap.png", label)
        plot_housing_size_fit(sol, P, b_grid, setup, moments, outdir / "housing_size_fit.png", label)
    write_case_summary(outdir / "target_comparison.csv", setup.targets, moments, label)

    return {
        "label": label,
        "outdir": outdir,
        "record_path": record_path,
        "hR_max": hr_max,
        "record": record,
        "targets": setup.targets,
        "moments": moments,
        "prices": [float(x) for x in np.asarray(record.get("p_eq", [])).reshape(-1)],
    }


def solve_record(record: dict, setup) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    theta = np.asarray(record["theta"], dtype=float)
    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    P = apply_theta(P, [theta_dict[name] for name in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if OUTSIDE_VALUE_NAME in theta_dict:
        P.outside_value = theta_dict[OUTSIDE_VALUE_NAME]
        P.outside_value_is_calibrated = True
        P.allow_uncalibrated_outside_value = False
    if RENEWAL_FLOW_NAME in theta_dict:
        P.outside_entry_flow = theta_dict[RENEWAL_FLOW_NAME]
    return run_model_cp_dt(P, verbose=False)


def plot_market_equilibrium(sol, P, p_eq, setup, out_path: Path, label: str) -> None:
    locs = ["Periphery", "Center"]
    prices = np.asarray(p_eq, dtype=float)
    rents = P.user_cost_rate * prices
    demand = np.asarray(sol.housing_demand, dtype=float) * housing_demand_normalizer(P)
    supply = np.asarray(P.H0, dtype=float) * (rents / np.asarray(P.r_bar, dtype=float)) ** np.asarray(P.eta_supply, dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    fig.suptitle(f"Market Equilibrium: {label}", fontsize=13)

    x = np.arange(len(locs))
    width = 0.35
    ax = axes[0, 0]
    ax.bar(x - width / 2, demand, width, label="Demand", color="tab:blue")
    ax.bar(x + width / 2, supply, width, label="Supply", color="tab:orange")
    ax.set_xticks(x, locs)
    ax.set_title("Housing quantities")
    ax.set_ylabel("Room-equivalent units")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[0, 1]
    ax.bar(x, prices, color=["tab:blue", "tab:red"])
    ax.set_xticks(x, locs)
    ax.set_title("House prices")
    ax.set_ylabel("Price per room-equivalent unit")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1, 0]
    pop_share = np.asarray(sol.pop_share, dtype=float)
    ax.bar(x, pop_share, color=["tab:blue", "tab:red"])
    ax.axhline(setup.targets["inv_pop_share_C"], color="black", linestyle="--", linewidth=1, label="Center target")
    ax.set_xticks(x, locs)
    ax.set_title(f"Population share, center={pop_share[1]:.3f}")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1, 1]
    ratio = float(rents[1] / rents[0])
    ax.bar([0, 1], [ratio, setup.targets["inv_rent_ratio_C_over_P"]], color=["tab:red", "gray"])
    ax.set_xticks([0, 1], ["Model", "Target"])
    ax.set_title("Center/periphery rent-per-room ratio")
    ax.set_ylim(0, max(1.8, ratio * 1.15))
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_housing_choice_heatmap(sol, P, b_grid: np.ndarray, out_path: Path, label: str) -> None:
    age_axis = np.arange(P.J) + P.age_start
    states = [(0, 0, "Childless"), (1, 1, "One young child")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
    fig.suptitle(f"Conditional Renter Housing Choice: {label}", fontsize=13)
    vmax = max(float(P.hR_max), 1.0)
    im = None
    for row, loc in enumerate([0, 1]):
        for col, (nn, cs, state_label) in enumerate(states):
            ax = axes[row, col]
            if nn >= P.n_parity or cs >= P.n_child_states:
                data = np.full((len(b_grid), P.J), np.nan)
            else:
                data = np.asarray(sol.hR_pol[:, 0, loc, :, nn, cs], dtype=float)
            im = ax.imshow(
                data,
                origin="lower",
                aspect="auto",
                extent=[age_axis[0], age_axis[-1], b_grid[0], b_grid[-1]],
                vmin=0,
                vmax=vmax,
                cmap="viridis",
            )
            ax.set_title(f"{['Periphery', 'Center'][loc]}, {state_label}")
            ax.set_xlabel("Age")
            ax.set_ylabel("Liquid wealth")
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85)
    cbar.set_label("Rooms")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_housing_size_fit(sol, P, b_grid: np.ndarray, setup, moments: dict, out_path: Path, label: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(f"Housing Size Fit: {label}", fontsize=13)

    moment_keys = [
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_1to2",
    ]
    labels = ["Renter med.", "Owner med.", "0->1 inc.", "1->2 inc."]
    target_vals = [float(setup.targets[k]) for k in moment_keys]
    model_vals = [float(moments.get(k, np.nan)) for k in moment_keys]
    x = np.arange(len(labels))
    width = 0.36
    ax = axes[0]
    ax.bar(x - width / 2, target_vals, width, label="Target", color="gray")
    ax.bar(x + width / 2, model_vals, width, label="Model", color="tab:blue")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("Rooms")
    ax.set_title("Targeted room moments")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    bins = ROOM_BINS
    acs = load_acs_room_bins()
    model_r = model_room_bin_shares(sol, P, b_grid, tenure="Renter", child_bin="0")
    acs_r = [acs.get(("Renter", "0", b), np.nan) for b in bins]
    x = np.arange(len(bins))
    ax.plot(x, acs_r, marker="o", color="black", label="ACS renter childless")
    ax.plot(x, model_r, marker="s", color="tab:blue", label="Model renter childless")
    ax.set_xticks(x, bins)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Share")
    ax.set_title("Prime-age childless renter room bins")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def load_acs_room_bins() -> dict[tuple[str, str, str], float]:
    path = ROOT / "output/model/acs_2023_rooms_bins.csv"
    out: dict[tuple[str, str, str], float] = {}
    if not path.exists():
        return out
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("age_window") != "25_45":
                continue
            out[(row["tenure"], row["child_bin"], row["bin"])] = float(row["share"])
    return out


def model_room_bin_shares(sol, P, b_grid: np.ndarray, *, tenure: str, child_bin: str) -> list[float]:
    ten_filter = (lambda ten: ten == 0) if tenure == "Renter" else (lambda ten: ten > 0)
    if child_bin == "0":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 0
    elif child_bin == "1":
        child_filter = lambda nn, cs: current_children(nn, cs, P) == 1
    else:
        child_filter = lambda nn, cs: current_children(nn, cs, P) >= 2
    vals, wts = cell_arrays(
        sol,
        P,
        b_grid,
        ages=range(age_index(P, 25), age_index(P, 45) + 1),
        mask_fn=lambda ten, loc, nn, cs: ten_filter(ten) and child_filter(nn, cs),
    )
    if wts.size == 0 or float(np.sum(wts)) <= 0:
        return [float("nan")] * len(ROOM_BINS)
    masks = room_bin_masks(vals["rooms"])
    total = float(np.sum(wts))
    return [float(np.sum(wts[masks[name]]) / total) for name in ROOM_BINS]


def write_case_summary(path: Path, targets: dict, moments: dict, label: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "moment", "target", "model", "diff"])
        for key, _ in TARGET_ROWS:
            if key not in targets:
                continue
            target = float(targets[key])
            model = float(moments.get(key, np.nan))
            writer.writerow([label, key, target, model, model - target])


def write_note(note_tex: Path, outdir: Path, cases: list[dict]) -> None:
    note_tex.parent.mkdir(parents=True, exist_ok=True)
    rel_base = Path("../") / outdir.resolve().relative_to(ROOT)
    rows = []
    for key, name in TARGET_ROWS:
        vals = []
        target = None
        for case in cases:
            if key in case["targets"]:
                target = float(case["targets"][key])
                vals.append(float(case["moments"].get(key, np.nan)))
            else:
                vals.append(np.nan)
        if target is None:
            continue
        rows.append((name, target, vals))

    lines = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[margin=0.7in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{float}",
        r"\usepackage{hyperref}",
        r"\title{Current Best-Fit Slide Diagnostics}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Source Runs}",
        r"\begin{itemize}",
    ]
    for case in cases:
        rec = case["record"]
        lines.append(
            rf"\item \textbf{{{latex_escape(case['label'])}}}: loss {float(rec['loss']):.3f}, "
            rf"job {rec.get('job_id')}, eval {rec.get('eval_id')}, $h_R^{{\max}}={case['hR_max']:.1f}$."
        )
    lines.extend(
        [
            r"\end{itemize}",
            r"\section*{Moment Table}",
            r"\begin{table}[H]\centering\small",
            r"\begin{tabular}{l" + "r" * (2 + len(cases)) + r"}",
            r"\toprule",
            "Moment & Target & " + " & ".join(latex_escape(c["label"]) for c in cases) + r" \\",
            r"\midrule",
        ]
    )
    for name, target, vals in rows:
        lines.append(
            latex_escape(name)
            + f" & {target:.3f} & "
            + " & ".join("NA" if not np.isfinite(v) else f"{v:.3f}" for v in vals)
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])

    lines.append(r"\section*{Slide Figure Bundle}")
    for case in cases:
        case_rel = rel_base / case["label"]
        lines.extend(
            [
                rf"\subsection*{{{latex_escape(case['label'])}}}",
                figure_line(case_rel / "benchmark_lifecycle_compare.png", "Lifecycle fit"),
                figure_line(case_rel / "market_equilibrium.png", "Market equilibrium"),
                figure_line(case_rel / "housing_choice_heatmap.png", "Housing choice heatmap"),
                figure_line(case_rel / "housing_size_fit.png", "Appendix housing-size fit"),
            ]
        )
    lines.append(r"\end{document}")
    note_tex.write_text("\n".join(lines) + "\n")


def figure_line(path: Path, caption: str) -> str:
    return (
        r"\begin{figure}[H]\centering"
        + "\n"
        + rf"\includegraphics[width=0.92\textwidth]{{{path.as_posix()}}}"
        + "\n"
        + rf"\caption{{{latex_escape(caption)}}}"
        + "\n"
        + r"\end{figure}"
    )


def latex_escape(text: str) -> str:
    return (
        str(text)
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
