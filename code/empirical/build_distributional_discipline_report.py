"""Build the data-vs-model distributional discipline PDF source.

This script assumes the empirical and model exporters have already run:

1. PSID income/wealth master sample and CSV.
2. Current Python model distributional exports.

It writes a LaTeX file and supporting PDF figures. The LaTeX file is compiled
separately from `latex/`.
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_ROOT = REPO_ROOT / "code" / "model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.direct_calibration import build_direct_calibration_setup  # noqa: E402


MODEL_OUT = REPO_ROOT / "output" / "model" / "distributional_discipline_current"
FIGDIR = MODEL_OUT / "report_figures"
LATEX = REPO_ROOT / "latex" / "distributional_empirics_report.tex"
BEST_JSON = (
    REPO_ROOT
    / "code"
    / "cluster"
    / "results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506"
    / "direct_geometry_best.json"
)
PSID_MASTER = (
    REPO_ROOT
    / "code"
    / "data"
    / "psid_followup_mar2026"
    / "output"
    / "income_wealth_fertility_master_v1"
    / "income_wealth_fertility_master_v1.csv"
)
PSID_WITHIN = (
    REPO_ROOT
    / "code"
    / "data"
    / "psid_followup_mar2026"
    / "output"
    / "fertility_wealth_v1"
    / "within_wealth_gap_v9.csv"
)
ACS_ROOMS = REPO_ROOT / "output" / "model" / "acs_2023_rooms_summary.csv"
MMS_LOCATION = (
    REPO_ROOT
    / "code"
    / "data"
    / "mms_center_periphery"
    / "output_middle_center"
    / "mms_location_by_parent_summary.csv"
)
MMS_TRANSITION = (
    REPO_ROOT
    / "code"
    / "data"
    / "mms_center_periphery"
    / "output_middle_center"
    / "mms_origin_transition_summary_allparent.csv"
)
PSID_MOVE_WINDOWS = (
    REPO_ROOT
    / "code"
    / "data"
    / "psid_followup_mar2026"
    / "output"
    / "moving_first_birth_window_summary_v1"
    / "moving_first_birth_window_summary_v1.csv"
)
PSID_MOVE_SA = (
    REPO_ROOT
    / "code"
    / "data"
    / "psid_followup_mar2026"
    / "output"
    / "sa_moving_first_birth_window_v1"
    / "moving_first_birth_sa_window_summary_v1.csv"
)


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    figures = build_figures(data)
    tex = build_tex(data, figures)
    LATEX.write_text(tex)
    print(f"Wrote {LATEX}")


def load_inputs() -> dict[str, pd.DataFrame | dict]:
    out: dict[str, pd.DataFrame | dict] = {}
    out["psid_master"] = read_csv(PSID_MASTER)
    out["psid_within"] = read_csv(PSID_WITHIN)
    out["acs_rooms"] = read_csv(ACS_ROOMS)
    out["mms_location"] = read_csv(MMS_LOCATION)
    out["mms_transition"] = read_csv(MMS_TRANSITION)
    out["psid_move_windows"] = read_csv(PSID_MOVE_WINDOWS)
    out["psid_move_sa"] = read_csv(PSID_MOVE_SA)
    out["model_metadata"] = read_csv(MODEL_OUT / "model_metadata.csv")
    out["model_income"] = read_csv(MODEL_OUT / "model_income_bins.csv")
    out["model_wealth"] = read_csv(MODEL_OUT / "model_wealth_bins.csv")
    out["model_rooms"] = read_csv(MODEL_OUT / "model_rooms_summary.csv")
    out["model_location"] = read_csv(MODEL_OUT / "model_location_tenure_by_parent.csv")
    out["model_mobility"] = read_csv(MODEL_OUT / "model_mobility_by_parent.csv")
    out["model_transition"] = read_csv(MODEL_OUT / "model_origin_transition_by_parent.csv")
    out["best_record"] = json.loads(BEST_JSON.read_text())[0]
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        geo_weight=100.0,
        population_closure="renewal_valve_calibrated",
        scale_target=1.0,
        renewal_retention=1.0,
    )
    out["targets"] = setup.targets
    return out


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in df.columns:
        if col not in {"wealth_concept", "group", "outcome", "bin_type", "domain", "age_bin", "concept", "parent_status", "tenure", "child_bin", "age_window", "origin_label", "dest_label", "parent_compare_all", "window", "label", "interpretation", "key", "value", "transition"}:
            df[col] = pd.to_numeric(df[col], errors="ignore")
    return df


def build_figures(data: dict) -> dict[str, Path]:
    figures = {
        "target_errors": FIGDIR / "target_errors.pdf",
        "model_income_profile": FIGDIR / "model_lifecycle_income_profile.pdf",
        "income_parent": FIGDIR / "parenthood_by_income_bins.pdf",
        "wealth_parent": FIGDIR / "parenthood_by_wealth_bins.pdf",
        "tenure_wealth_childless": FIGDIR / "childlessness_by_tenure_wealth.pdf",
        "rooms_mean": FIGDIR / "rooms_mean_by_children.pdf",
        "rooms_share7": FIGDIR / "rooms_share_ge7_by_children.pdf",
        "location_center": FIGDIR / "location_center_share.pdf",
        "location_owner": FIGDIR / "location_owner_rate.pdf",
        "location_rooms": FIGDIR / "location_rooms_by_parent.pdf",
        "transition": FIGDIR / "center_periphery_transitions.pdf",
        "mobility": FIGDIR / "first_birth_mobility.pdf",
    }

    make_target_error_plot(data, figures["target_errors"])
    make_model_income_profile(figures["model_income_profile"])

    psid = data["psid_master"]
    model_income = data["model_income"]
    model_wealth = data["model_wealth"]

    data_income = psid_slice(psid, "income_k", "all", "parent_by_45", "quintile")
    data_wealth = psid_slice(psid, "total_nw_k", "all", "parent_by_45", "quintile")
    model_inc = model_income.query("age_bin == 'ages25_30' and concept == 'current_income' and group == 'all'")
    model_w = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'all'")

    x = np.arange(1, 6)
    line_plot(
        figures["income_parent"],
        x,
        [
            ("PSID: income at ages 25--30", 100 * data_income["mean_outcome"].to_numpy()),
            ("Model: deterministic income rank, ages 25--30", 100 * model_wide(model_inc, "parent_share")),
        ],
        "Parent by 45 / parent share (percent)",
        "Young-adult income quintile",
        ylim=(0, 105),
    )
    line_plot(
        figures["wealth_parent"],
        x,
        [
            ("PSID: total wealth at ages 25--30", 100 * data_wealth["mean_outcome"].to_numpy()),
            ("Model: sale-net wealth at age 45", 100 * model_wide(model_w, "parent_share")),
        ],
        "Parent by 45 / parent share (percent)",
        "Wealth quintile",
        ylim=(0, 105),
    )

    psid_within = data["psid_within"]
    d_stay = psid_within[psid_within["transition"].eq("Stayed renter")].sort_values("wealth_q")
    d_owner = psid_within[psid_within["transition"].eq("Became owner")].sort_values("wealth_q")
    m_renter = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'renters'")
    m_owner = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'owners'")
    line_plot(
        figures["tenure_wealth_childless"],
        x,
        [
            ("Data: stayed renter", 100 * d_stay["childless_by_45"].to_numpy()),
            ("Data: became owner", 100 * d_owner["childless_by_45"].to_numpy()),
            ("Model: renter at age 45", 100 * model_wide(m_renter, "childless_share")),
            ("Model: owner at age 45", 100 * model_wide(m_owner, "childless_share")),
        ],
        "Childless share (percent)",
        "Wealth quintile",
        ylim=(0, 85),
    )

    acs = data["acs_rooms"]
    model_rooms = data["model_rooms"]
    child_order = ["0", "1", "2+"]
    rooms_facet_plot(acs, model_rooms, "mean", figures["rooms_mean"], "Mean rooms", child_order, ylim=(3, 8))
    rooms_facet_plot(acs, model_rooms, "share_ge_7", figures["rooms_share7"], "Share with at least 7 rooms (percent)", child_order, scale=100, ylim=(0, 100))

    mms = data["mms_location"].set_index("parent_status")
    model_loc = data["model_location"].set_index("parent_status")
    statuses = ["Non-Parents", "New Parents", "Older Parents"]
    grouped_bar_plot(
        figures["location_center"],
        statuses,
        [("MMS data", [100 * mms.loc[s, "center_share"] for s in statuses]), ("Model", [100 * model_loc.loc[s, "center_share"] for s in statuses])],
        "Center share (percent)",
        ylim=(0, 65),
    )
    grouped_bar_plot(
        figures["location_owner"],
        statuses,
        [("MMS data", [100 * mms.loc[s, "owner_rate"] for s in statuses]), ("Model", [100 * model_loc.loc[s, "owner_rate"] for s in statuses])],
        "Owner rate (percent)",
        ylim=(0, 75),
    )
    grouped_bar_plot(
        figures["location_rooms"],
        statuses,
        [("MMS data", [mms.loc[s, "mean_rooms"] for s in statuses]), ("Model", [model_loc.loc[s, "mean_rooms"] for s in statuses])],
        "Mean rooms",
        ylim=(4.5, 7.8),
    )

    transition_plot(data["mms_transition"], data["model_transition"], figures["transition"])
    mobility_plot(data["psid_move_windows"], data["model_mobility"], figures["mobility"])

    return figures


def build_tex(data: dict, figures: dict[str, Path]) -> str:
    calibration = calibration_table(data)
    meta = dict(zip(data["model_metadata"]["key"], data["model_metadata"]["value"]))

    return rf"""\documentclass[11pt]{{article}}
\usepackage[margin=0.8in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{float}}
\usepackage{{graphicx}}
\usepackage{{adjustbox}}
\usepackage{{caption}}
\usepackage{{array}}
\usepackage{{hyperref}}

\title{{Distributional Empirical Discipline: Data and Current Model}}
\author{{}}
\date{{Built May 8, 2026}}

\begin{{document}}
\maketitle

\noindent This note is a discipline file, not a validation claim. For each empirical object that is currently available, it shows the data object and the closest object in the active Python model. The model is the current direct-geometry benchmark, job {latex_escape(str(meta.get("job_id", "")))}, evaluation {latex_escape(str(meta.get("eval_id", "")))}, with the renewal-valve stationary closure. The benchmark loss in the cluster record is {fmt_num(meta.get("loss"), 3)} and the resolved equilibrium error is {fmt_num(meta.get("final_eq_error"), 4)}.

\medskip
\noindent The organizing economic references are standard: Becker and Becker--Lewis for the fertility quantity-quality tradeoff; Sommer and Doepke--Kindermann style lifecycle fertility models for income and lifecycle discipline; Rosen--Roback, Henderson--Ioannides, Sommer--Sullivan--Verbrugge, and Gyourko--Mayer--Sinai for spatial sorting, tenure, and housing-price discipline. The point of this file is to make clear which of those cross-sectional facts the current model can speak to now.

\medskip
\noindent Income needs special care. The model has a deterministic lifecycle income profile, so income moves with age. It does not have persistent idiosyncratic income risk or a within-age income distribution. Therefore the model income-bin plot below uses the stationary mass generated by deterministic age/location income; it is a lifecycle-rank diagnostic, not a clean analog of PSID young-adult income quintiles.

{calibration}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.86\textwidth]{{\detokenize{{{rel_to_latex(figures["target_errors"])}}}}}
\caption{{Calibration target deviations. Positive values mean the model is above the target.}}
\end{{figure}}

\section*{{Income}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.74\textwidth]{{\detokenize{{{rel_to_latex(figures["model_income_profile"])}}}}}
\caption{{Model lifecycle income profile. This is the deterministic age component that creates income dynamics in the model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.82\textwidth]{{\detokenize{{{rel_to_latex(figures["income_parent"])}}}}}
\caption{{Parenthood by young-adult income bins. Data use PSID family income averaged over ages 25--30. The model uses ages 25--30 current income, which varies with age and location but not with idiosyncratic permanent income. Missing model points are empty rank bins.}}
\end{{figure}}

\section*{{Wealth And Tenure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.82\textwidth]{{\detokenize{{{rel_to_latex(figures["wealth_parent"])}}}}}
\caption{{Parenthood by wealth bins. Data use young-adult PSID total wealth; model uses cross-sectional sale-net wealth at age 45.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.82\textwidth]{{\detokenize{{{rel_to_latex(figures["tenure_wealth_childless"])}}}}}
\caption{{Childlessness by wealth and tenure. The data object follows young renters by whether they become owners by age 35; the model object is current tenure at age 45.}}
\end{{figure}}

\section*{{Housing Space}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.86\textwidth]{{\detokenize{{{rel_to_latex(figures["rooms_mean"])}}}}}
\caption{{Mean rooms by tenure and children, ACS versus model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.86\textwidth]{{\detokenize{{{rel_to_latex(figures["rooms_share7"])}}}}}
\caption{{Share of households with at least 7 rooms by tenure and children, ACS versus model.}}
\end{{figure}}

\section*{{Space And Mobility}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\textwidth]{{\detokenize{{{rel_to_latex(figures["location_center"])}}}}}
\caption{{Center residence shares by parent status, MMS data versus model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\textwidth]{{\detokenize{{{rel_to_latex(figures["location_owner"])}}}}}
\caption{{Owner rates by parent status, MMS data versus model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\textwidth]{{\detokenize{{{rel_to_latex(figures["location_rooms"])}}}}}
\caption{{Mean rooms by parent status, MMS data versus model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.90\textwidth]{{\detokenize{{{rel_to_latex(figures["transition"])}}}}}
\caption{{Center/periphery origin-destination shares. The model analog is the one-period location-choice probability averaged over the stationary distribution.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.78\textwidth]{{\detokenize{{{rel_to_latex(figures["mobility"])}}}}}
\caption{{First-birth moving facts. The current model has an analog for any location change, but no move-for-size or moved-to-ownership reason code.}}
\end{{figure}}

\section*{{Read}}

\noindent The model has genuine lifecycle income dynamics, but not the cross-sectional income heterogeneity needed to match PSID income quintile fertility facts. The wealth comparison is closer, but still not exact because the data object is young-adult wealth predicting later fertility while the model object is a stationary cross-section. The room and spatial plots are the cleanest diagnostics: renter housing services are too high, the model does not generate enough room dispersion by tenure and children, and center/periphery transition rates are much too low. These are useful failures because they identify exactly which empirical margins the current model cannot yet claim to validate.

\end{{document}}
"""


def make_target_error_plot(data: dict, path: Path) -> None:
    targets: dict = data["targets"]
    moments = data["best_record"]["moments"]
    order = [
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "tfr_gradient",
        "own_rate",
        "own_gradient",
        "own_family_gap",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "young_liquid_wealth_to_income",
        "center_share_nonparents",
        "center_share_newparents",
        "migration_rate",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "inv_pop_share_C",
        "inv_rent_ratio_C_over_P",
    ]
    labels = [short_moment_name(k) for k in order]
    vals = []
    for k in order:
        target = float(targets[k])
        model = float(moments[k])
        vals.append(100 * (model - target) / max(abs(target), 1e-12))
    fig, ax = plt.subplots(figsize=(8.3, 6.1))
    y = np.arange(len(order))
    colors = ["#d95f02" if v > 0 else "#1b9e77" for v in vals]
    ax.barh(y, vals, color=colors, alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Percent deviation from target")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_model_income_profile(path: Path) -> None:
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        geo_weight=100.0,
        population_closure="renewal_valve_calibrated",
        scale_target=1.0,
        renewal_retention=1.0,
    )
    P = setup.P_base
    ages = np.arange(18, 76)
    income = np.interp(ages, P.income_age_breaks, P.income_age_values)
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(ages, income, marker="o", markevery=5, linewidth=2.0)
    ax.set_xlabel("Age")
    ax.set_ylabel("Model income, relative units")
    ax.set_ylim(0, max(income) * 1.15)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def line_plot(path: Path, x: np.ndarray, series: list[tuple[str, np.ndarray]], ylabel: str, xlabel: str, ylim: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.1))
    for label, y in series:
        y = np.asarray(y, dtype=float)
        ax.plot(x, y, marker="o", linewidth=2.0, label=label)
        for xi, yi in zip(x, y):
            if np.isfinite(yi):
                ax.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def rooms_facet_plot(
    acs: pd.DataFrame,
    model_rooms: pd.DataFrame,
    value_col: str,
    path: Path,
    ylabel: str,
    child_order: list[str],
    *,
    scale: float = 1.0,
    ylim: tuple[float, float] | None = None,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), sharey=True)
    x = np.arange(len(child_order))
    for ax, tenure in zip(axes, ["Renter", "Owner"]):
        d = acs[(acs["age_window"].eq("25_45")) & (acs["tenure"].eq(tenure)) & (acs["child_bin"].isin(child_order))]
        m = model_rooms[(model_rooms["age_window"].eq("25_45")) & (model_rooms["tenure"].eq(tenure)) & (model_rooms["child_bin"].isin(child_order))]
        d = d.set_index("child_bin").loc[child_order]
        m = m.set_index("child_bin").loc[child_order]
        ax.plot(x, scale * d[value_col].to_numpy(), marker="o", linewidth=2.0, label="Data")
        ax.plot(x, scale * m[value_col].to_numpy(), marker="s", linewidth=2.0, label="Model")
        ax.set_title(tenure)
        ax.set_xticks(x, child_order)
        ax.set_xlabel("Children in household")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(ylabel)
    if ylim is not None:
        axes[0].set_ylim(*ylim)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def grouped_bar_plot(path: Path, xlabels: list[str], series: list[tuple[str, list[float]]], ylabel: str, ylim: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    x = np.arange(len(xlabels))
    width = 0.75 / max(len(series), 1)
    start = -0.5 * width * (len(series) - 1)
    for idx, (label, vals) in enumerate(series):
        vals = np.asarray(vals, dtype=float)
        ax.bar(x + start + idx * width, vals, width, label=label)
        for xi, yi in zip(x + start + idx * width, vals):
            if np.isfinite(yi):
                ax.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7)
    ax.set_xticks(x, xlabels)
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def transition_plot(data_transition: pd.DataFrame, model_transition: pd.DataFrame, path: Path) -> None:
    groups = ["All Parents", "Non-Parents"]
    pairs = [("Center", "Center"), ("Center", "Periphery"), ("Periphery", "Center"), ("Periphery", "Periphery")]
    labels = ["C->C", "C->P", "P->C", "P->P"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), sharey=True)
    for ax, group in zip(axes, groups):
        dvals = []
        mvals = []
        for origin, dest in pairs:
            d = data_transition[
                data_transition["parent_compare_all"].eq(group)
                & data_transition["origin_label"].eq(origin)
                & data_transition["dest_label"].eq(dest)
            ]
            m = model_transition[
                model_transition["parent_compare_all"].eq(group)
                & model_transition["origin_label"].eq(origin)
                & model_transition["dest_label"].eq(dest)
            ]
            dvals.append(100 * float(d.iloc[0]["share"]))
            mvals.append(100 * float(m.iloc[0]["share"]))
        x = np.arange(len(labels))
        width = 0.36
        ax.bar(x - width / 2, dvals, width, label="MMS data")
        ax.bar(x + width / 2, mvals, width, label="Model")
        ax.set_title(group)
        ax.set_xticks(x, labels)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Transition share (percent)")
    axes[0].set_ylim(0, 105)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def mobility_plot(windows: pd.DataFrame, model_mobility: pd.DataFrame, path: Path) -> None:
    outcomes = [("move_any", "Any move"), ("moved_for_size", "Moved for size"), ("moved_to_own", "Moved to own")]
    data_vals = []
    for outcome, _ in outcomes:
        row = windows[(windows["outcome"].eq(outcome)) & (windows["window"].eq("post0to3"))]
        data_vals.append(100 * float(row.iloc[0]["mean"]))
    model = model_mobility.set_index("group")
    model_vals = [100 * float(model.loc["Current Parents", "model_move_probability"]), np.nan, np.nan]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    x = np.arange(len(outcomes))
    width = 0.36
    ax.bar(x - width / 2, data_vals, width, label="PSID data")
    ax.bar(x + width / 2, np.nan_to_num(model_vals, nan=0.0), width, label="Model")
    for xi, val in zip(x - width / 2, data_vals):
        ax.annotate(f"{val:.1f}", (xi, val), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7)
    for xi, val in zip(x + width / 2, model_vals):
        if np.isfinite(val):
            ax.annotate(f"{val:.1f}", (xi, val), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7)
        else:
            ax.annotate("no analog", (xi, 2), rotation=90, ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x, [label for _, label in outcomes])
    ax.set_ylabel("Post-birth move rate / model location-change probability (percent)")
    ax.set_ylim(0, 80)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def short_moment_name(name: str) -> str:
    replacements = {
        "childless_rate": "childless",
        "mean_age_first_birth": "age first birth",
        "tfr_gradient": "TFR gradient",
        "own_rate": "ownership",
        "own_gradient": "own gradient",
        "own_family_gap": "own family gap",
        "prime_childless_renter_median_rooms": "renter rooms",
        "prime_childless_owner_median_rooms": "owner rooms",
        "housing_increment_0to1": "rooms 0->1",
        "housing_increment_1to2": "rooms 1->2",
        "young_liquid_wealth_to_income": "young liq wealth/inc",
        "center_share_nonparents": "center nonparents",
        "center_share_newparents": "center new parents",
        "migration_rate": "migration",
        "old_age_own_rate": "old ownership",
        "old_age_parent_childless_gap": "old own parent gap",
        "inv_pop_share_C": "center pop",
        "inv_rent_ratio_C_over_P": "rent ratio",
    }
    return replacements.get(name, name.upper() if name == "tfr" else name.replace("_", " "))


def calibration_table(data: dict) -> str:
    targets: dict = data["targets"]
    moments = data["best_record"]["moments"]
    order = [
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "tfr_gradient",
        "own_rate",
        "own_gradient",
        "own_family_gap",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "young_liquid_wealth_to_income",
        "center_share_nonparents",
        "center_share_newparents",
        "migration_rate",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "inv_pop_share_C",
        "inv_rent_ratio_C_over_P",
    ]
    rows = [[k, fmt_num(targets.get(k), 3), fmt_num(moments.get(k), 3)] for k in order]
    return latex_table(
        "Benchmark calibration moments, with targets next to model moments.",
        "tab:calibration",
        ["Moment", "Target", "Model"],
        rows,
        "lrr",
        note="These are the live benchmark moments from the current direct-geometry calibration record. Distributional tables below are not additional calibration targets.",
    )


def income_wealth_table(data: dict) -> str:
    psid = data["psid_master"]
    model_income = data["model_income"]
    model_wealth = data["model_wealth"]
    d_income = psid_slice(psid, "income_k", "all", "parent_by_45", "quintile").set_index("bin")
    d_wealth = psid_slice(psid, "total_nw_k", "all", "parent_by_45", "quintile").set_index("bin")
    d_childless = psid_slice(psid, "total_nw_k", "all", "childless_by_45", "quintile").set_index("bin")
    m_income = model_income.query("age_bin == 'ages25_45' and concept == 'current_income' and group == 'all'").set_index("bin")
    m_wealth = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'all'").set_index("bin")
    rows = []
    for b in range(1, 6):
        rows.append(
            [
                str(b),
                fmt_pct(d_income.loc[b, "mean_outcome"]),
                fmt_pct(d_wealth.loc[b, "mean_outcome"]),
                fmt_pct(d_childless.loc[b, "mean_outcome"]),
                fmt_pct(m_income.loc[b, "parent_share"]) if b in m_income.index else "--",
                fmt_pct(m_wealth.loc[b, "parent_share"]) if b in m_wealth.index else "--",
                fmt_num(m_wealth.loc[b, "mean_completed_fertility"], 2) if b in m_wealth.index else "--",
            ]
        )
    return latex_table(
        "Fertility by income and wealth bins: data and model.",
        "tab:income-wealth",
        [
            "Bin",
            "Data income parent by 45",
            "Data wealth parent by 45",
            "Data wealth childless by 45",
            "Model income parent",
            "Model wealth parent",
            "Model wealth children",
        ],
        rows,
        "lrrrrrr",
        note="Data are PSID bins measured at ages 25--30. The model income column is a stationary cross-section over ages 25--45; most bins are empty because the current model does not contain idiosyncratic income heterogeneity. The model wealth column uses sale-net wealth at age 45.",
    )


def tenure_wealth_table(data: dict) -> str:
    psid = data["psid_within"].copy()
    model = data["model_wealth"]
    m_renter = model.query("age_bin == 'age45' and concept == 'sale_net' and group == 'renters'").set_index("bin")
    m_owner = model.query("age_bin == 'age45' and concept == 'sale_net' and group == 'owners'").set_index("bin")
    rows = []
    for b in range(1, 6):
        d_r = psid[(psid["wealth_q"].eq(b)) & (psid["transition"].eq("Stayed renter"))]
        d_o = psid[(psid["wealth_q"].eq(b)) & (psid["transition"].eq("Became owner"))]
        rows.append(
            [
                str(b),
                fmt_pct(d_r.iloc[0]["childless_by_45"]) if len(d_r) else "--",
                fmt_pct(d_o.iloc[0]["childless_by_45"]) if len(d_o) else "--",
                fmt_pct(m_renter.loc[b, "childless_share"]) if b in m_renter.index else "--",
                fmt_pct(m_owner.loc[b, "childless_share"]) if b in m_owner.index else "--",
            ]
        )
    return latex_table(
        "Childlessness by wealth and tenure status.",
        "tab:tenure-wealth",
        ["Wealth bin", "Data stayed renter", "Data became owner", "Model renter", "Model owner"],
        rows,
        "lrrrr",
        note="The PSID data classify young renters by whether they become owners by age 35. The model columns are current age-45 tenure cross-sections, so this is not yet a true transition analog.",
    )


def rooms_table(data: dict) -> str:
    acs = data["acs_rooms"]
    model = data["model_rooms"]
    rows = []
    for tenure in ["Renter", "Owner"]:
        for child in ["0", "1", "2+"]:
            d = acs[(acs["age_window"].eq("25_45")) & (acs["tenure"].eq(tenure)) & (acs["child_bin"].eq(child))]
            m = model[(model["age_window"].eq("25_45")) & (model["tenure"].eq(tenure)) & (model["child_bin"].eq(child))]
            rows.append(
                [
                    tenure,
                    child,
                    fmt_num(d.iloc[0]["mean"], 2),
                    fmt_num(d.iloc[0]["p50"], 1),
                    fmt_pct(d.iloc[0]["share_ge_7"]),
                    fmt_num(m.iloc[0]["mean"], 2),
                    fmt_num(m.iloc[0]["p50"], 1),
                    fmt_pct(m.iloc[0]["share_ge_7"]),
                ]
            )
    return latex_table(
        "Rooms by tenure and children, ages 25--45.",
        "tab:rooms",
        ["Tenure", "Children", "Data mean", "Data p50", "Data share >=7", "Model mean", "Model p50", "Model share >=7"],
        rows,
        "llrrrrrr",
        note="Data are ACS 2023 room counts. Model rooms are room-equivalent services in the current benchmark.",
    )


def location_table(data: dict) -> str:
    d = data["mms_location"].set_index("parent_status")
    m = data["model_location"].set_index("parent_status")
    rows = []
    for status in ["Non-Parents", "New Parents", "Older Parents"]:
        rows.append(
            [
                status,
                fmt_pct(d.loc[status, "center_share"]),
                fmt_pct(m.loc[status, "center_share"]),
                fmt_pct(d.loc[status, "owner_rate"]),
                fmt_pct(m.loc[status, "owner_rate"]),
                fmt_num(d.loc[status, "mean_rooms"], 2),
                fmt_num(m.loc[status, "mean_rooms"], 2),
            ]
        )
    return latex_table(
        "Location, tenure, and rooms by parent status.",
        "tab:location",
        ["Status", "Data center", "Model center", "Data owner", "Model owner", "Data rooms", "Model rooms"],
        rows,
        "lrrrrrr",
        note="Data are the MMS center/periphery summary with middle locations assigned to center in the binary output. The model table uses ages 22--45.",
    )


def transition_table(data: dict) -> str:
    d = data["mms_transition"]
    m = data["model_transition"]
    rows = []
    for group in ["All Parents", "Non-Parents"]:
        for origin, dest in [("Center", "Center"), ("Center", "Periphery"), ("Periphery", "Center"), ("Periphery", "Periphery")]:
            dd = d[(d["parent_compare_all"].eq(group)) & (d["origin_label"].eq(origin)) & (d["dest_label"].eq(dest))]
            mm = m[(m["parent_compare_all"].eq(group)) & (m["origin_label"].eq(origin)) & (m["dest_label"].eq(dest))]
            rows.append([group, origin, dest, fmt_pct(dd.iloc[0]["share"]), fmt_pct(mm.iloc[0]["share"])])
    return latex_table(
        "Center/periphery transition shares by parent status.",
        "tab:transition",
        ["Group", "Origin", "Destination", "Data share", "Model share"],
        rows,
        "lllrr",
        note="Data are MMS origin-destination shares. The model analog is the one-period location-choice probability averaged over the stationary distribution.",
    )


def mobility_table(data: dict) -> str:
    windows = data["psid_move_windows"]
    sa = data["psid_move_sa"].set_index("outcome")
    model = data["model_mobility"].set_index("group")
    rows = []
    for outcome, label in [
        ("move_any", "Any move"),
        ("moved_for_size", "Moved for size"),
        ("moved_to_own", "Moved to ownership"),
    ]:
        w = windows[(windows["outcome"].eq(outcome)) & (windows["window"].eq("post0to3"))]
        rows.append(
            [
                label,
                fmt_pct(w.iloc[0]["mean"]) if len(w) else "--",
                fmt_pp(sa.loc[outcome, "coef_0"]) if outcome in sa.index else "--",
                fmt_pp(sa.loc[outcome, "se_0"]) if outcome in sa.index else "--",
                fmt_pct(model.loc["Current Parents", "model_move_probability"]) if outcome == "move_any" else "No direct analog",
            ]
        )
    return latex_table(
        "First-birth moving facts and model analog.",
        "tab:mobility",
        ["Object", "Data post 0--3 mean", "SA coef at birth", "SA s.e.", "Model analog"],
        rows,
        "lrrrl",
        note="The current model records location changes, not PSID move reasons. Therefore moved-for-size and moved-to-ownership are empirical discipline objects but not exact model moments yet.",
    )


def psid_slice(df: pd.DataFrame, concept: str, group: str, outcome: str, bin_type: str) -> pd.DataFrame:
    out = df[
        df["wealth_concept"].eq(concept)
        & df["group"].eq(group)
        & df["outcome"].eq(outcome)
        & df["bin_type"].eq(bin_type)
    ].copy()
    out.sort_values("bin", inplace=True)
    return out


def model_wide(df: pd.DataFrame, col: str) -> np.ndarray:
    vals = np.full(5, np.nan)
    for _, row in df.iterrows():
        b = int(row["bin"])
        if 1 <= b <= 5:
            vals[b - 1] = row[col]
    return vals


def latex_table(caption: str, label: str, headers: list[str], rows: list[list[str]], align: str, note: str = "") -> str:
    head = " & ".join(latex_escape(h) for h in headers) + r" \\"
    body = "\n".join(" & ".join(latex_escape(str(x)) for x in row) + r" \\" for row in rows)
    note_block = ""
    if note:
        note_block = rf"\begin{{minipage}}{{0.97\textwidth}}\footnotesize \emph{{Note:}} {latex_escape(note)}\end{{minipage}}"
    return rf"""
\begin{{table}}[H]
\centering
\caption{{{latex_escape(caption)}}}
\label{{{label}}}
\scriptsize
\begin{{adjustbox}}{{max width=\textwidth}}
\begin{{tabular}}{{{align}}}
\toprule
{head}
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
{note_block}
\end{{table}}
"""


def fmt_num(x, digits: int = 2) -> str:
    try:
        value = float(x)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(value):
        return "--"
    return f"{value:.{digits}f}"


def fmt_pct(x) -> str:
    try:
        value = float(x)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(value):
        return "--"
    return f"{100 * value:.1f}%"


def fmt_pp(x) -> str:
    try:
        value = float(x)
    except (TypeError, ValueError):
        return "--"
    if not math.isfinite(value):
        return "--"
    return f"{100 * value:.1f}"


def latex_escape(s: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in s)


def rel_to_latex(path: Path) -> str:
    return os.path.relpath(path, LATEX.parent).replace("\\", "/")


if __name__ == "__main__":
    main()
