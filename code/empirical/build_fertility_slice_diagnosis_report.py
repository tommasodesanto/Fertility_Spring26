"""Build a data-vs-model fertility-slice diagnostic PDF.

This is a diagnostic report, not a calibration target file. It compares the
active Python benchmark to available empirical fertility cuts by age, location,
tenure, income, and wealth.
"""

from __future__ import annotations

import json
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
OUTDIR = REPO_ROOT / "output" / "model" / "fertility_slice_diagnosis_current"
FIGDIR = OUTDIR / "report_figures"
LATEX = REPO_ROOT / "latex" / "fertility_slice_diagnosis_report.tex"
BEST_JSON = (
    REPO_ROOT
    / "code"
    / "cluster"
    / "results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506"
    / "direct_geometry_best.json"
)
MMS = REPO_ROOT / "code" / "data" / "mms_center_periphery" / "output_middle_center"
PSID = REPO_ROOT / "code" / "data" / "psid_followup_mar2026" / "output"
MMS_AGE = MMS / "mms_age_fertility_profiles.csv"
MMS_LOCATION = MMS / "mms_location_summary.csv"
MMS_TENURE = MMS / "mms_tenure_parent_summary.csv"
MMS_TFR = MMS / "reconstructed_fertility_targets.csv"
PSID_MASTER = PSID / "income_wealth_fertility_master_v1" / "income_wealth_fertility_master_v1.csv"
PSID_TENURE_WEALTH = PSID / "fertility_wealth_v1" / "within_wealth_gap_v9.csv"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    figures = build_figures(data)
    summary = build_summary(data)
    summary.to_csv(OUTDIR / "fertility_slice_gap_summary.csv", index=False)
    LATEX.write_text(build_tex(data, figures, summary))
    print(f"Wrote {LATEX}")
    print(f"Wrote {OUTDIR / 'fertility_slice_gap_summary.csv'}")


def load_inputs() -> dict[str, pd.DataFrame | dict]:
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        geo_weight=100.0,
        population_closure="renewal_valve_calibrated",
        scale_target=1.0,
        renewal_retention=1.0,
    )
    return {
        "mms_age": read_csv(MMS_AGE),
        "mms_location": read_csv(MMS_LOCATION),
        "mms_tenure": read_csv(MMS_TENURE),
        "mms_tfr": read_csv(MMS_TFR),
        "psid_master": read_csv(PSID_MASTER),
        "psid_tenure_wealth": read_csv(PSID_TENURE_WEALTH),
        "model_fertility": read_csv(MODEL_OUT / "model_fertility_slices.csv"),
        "model_income": read_csv(MODEL_OUT / "model_income_bins.csv"),
        "model_wealth": read_csv(MODEL_OUT / "model_wealth_bins.csv"),
        "model_metadata": read_csv(MODEL_OUT / "model_metadata.csv"),
        "best_record": json.loads(BEST_JSON.read_text())[0],
        "targets": setup.targets,
    }


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    text_cols = {
        "dest_label",
        "mms_location",
        "tenure",
        "parent_status",
        "wealth_concept",
        "group",
        "outcome",
        "bin_type",
        "domain",
        "age_bin",
        "age_window",
        "concept",
        "slice",
        "location",
        "key",
        "value",
        "transition",
    }
    for col in df.columns:
        if col not in text_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def build_figures(data: dict) -> dict[str, Path]:
    figures = {
        "age": FIGDIR / "fertility_by_age_location.pdf",
        "location": FIGDIR / "fertility_by_location.pdf",
        "tenure": FIGDIR / "fertility_by_tenure.pdf",
        "income": FIGDIR / "fertility_by_income.pdf",
        "wealth": FIGDIR / "fertility_by_wealth.pdf",
        "tenure_wealth": FIGDIR / "fertility_by_wealth_within_tenure.pdf",
    }
    make_age_plot(data["mms_age"], data["model_fertility"], figures["age"])
    make_location_plot(data["mms_tfr"], data["mms_location"], data["model_fertility"], figures["location"])
    make_tenure_plot(data["mms_tenure"], data["model_fertility"], figures["tenure"])
    make_income_plot(data["psid_master"], data["model_income"], figures["income"])
    make_wealth_plot(data["psid_master"], data["model_wealth"], figures["wealth"])
    make_tenure_wealth_plot(data["psid_tenure_wealth"], data["model_wealth"], figures["tenure_wealth"])
    return figures


def make_age_plot(mms_age: pd.DataFrame, model: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True)
    for ax, loc in zip(axes, ["Center", "Periphery"]):
        data_loc = mms_age[mms_age["dest_label"].eq(loc)].sort_values("age")
        model_loc = model[
            model["slice"].eq("age") & model["location"].eq(loc) & model["tenure"].eq("All")
        ].sort_values("age")
        ax.plot(data_loc["age"], 100 * data_loc["has_children_rate"], linewidth=2.2, label="MMS data")
        ax.plot(model_loc["age"], 100 * model_loc["current_parent_share"], linewidth=2.2, label="Model")
        ax.set_title(loc)
        ax.set_xlabel("Age")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Current parent share (percent)")
    axes[0].set_ylim(0, 85)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_location_plot(mms_tfr: pd.DataFrame, mms_loc: pd.DataFrame, model: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
    locs = ["Center", "Periphery"]
    data_tfr = [float(mms_tfr.iloc[0]["tfr_center"]), float(mms_tfr.iloc[0]["tfr_periphery"])]
    model_tfr = [
        model_value(model, "location", "age45", loc, "All", "tfr") for loc in locs
    ]
    grouped_bar(
        axes[0],
        locs,
        [("MMS target", data_tfr), ("Model age 45", model_tfr)],
        "TFR / TFR-equivalent",
        ylim=(0, 2.4),
    )
    mms_by_loc = mms_loc.assign(location=mms_loc["mms_location"].map({"center": "Center", "periphery": "Periphery"}))
    data_parent = [float(mms_by_loc[mms_by_loc["location"].eq(loc)].iloc[0]["has_children_rate"]) * 100 for loc in locs]
    model_parent = [
        model_value(model, "location", "ages22_45", loc, "All", "current_parent_share") * 100 for loc in locs
    ]
    grouped_bar(
        axes[1],
        locs,
        [("MMS data", data_parent), ("Model", model_parent)],
        "Current parent share (percent)",
        ylim=(0, 65),
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_tenure_plot(mms_tenure: pd.DataFrame, model: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    tenure_order = ["Renter", "Owner"]
    data_overall = tenure_parent_share(mms_tenure, location=None)
    model_overall = {
        ten: model_value(model, "tenure", "ages22_45", "All", ten, "current_parent_share")
        for ten in tenure_order
    }
    grouped_bar(
        axes[0],
        tenure_order,
        [
            ("MMS data", [100 * data_overall[ten] for ten in tenure_order]),
            ("Model", [100 * model_overall[ten] for ten in tenure_order]),
        ],
        "Current parent share (percent)",
        ylim=(0, 65),
    )
    axes[0].set_title("Overall")

    categories = [("Center", "Renter"), ("Center", "Owner"), ("Periphery", "Renter"), ("Periphery", "Owner")]
    labels = ["C Rent", "C Own", "P Rent", "P Own"]
    data_cells = []
    model_cells = []
    for loc, ten in categories:
        data_cells.append(100 * tenure_parent_share(mms_tenure, loc.lower())[ten])
        model_cells.append(100 * model_value(model, "location_tenure", "ages22_45", loc, ten, "current_parent_share"))
    grouped_bar(axes[1], labels, [("MMS data", data_cells), ("Model", model_cells)], "", ylim=(0, 70))
    axes[1].set_title("By Location")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_income_plot(psid: pd.DataFrame, model_income: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.7))
    x = np.arange(1, 6)
    data_parent = psid_slice(psid, "income_k", "all", "parent_by_45").set_index("bin")
    data_childless = psid_slice(psid, "income_k", "all", "childless_by_45").set_index("bin")
    data_children = psid_slice(psid, "income_k", "all", "children_43_50").set_index("bin")
    model = model_income.query("age_bin == 'age45' and concept == 'current_income' and group == 'all'").set_index("bin")
    panels = [
        ("Parent by 45 / parent share", data_parent["mean_outcome"] * 100, model["parent_share"] * 100, (0, 105), "Percent"),
        ("Childless by 45 / childless share", data_childless["mean_outcome"] * 100, model["childless_share"] * 100, (0, 45), "Percent"),
        ("Completed children", data_children["mean_outcome"], 2.0 * model["mean_completed_fertility"], (0, 2.8), "Count"),
    ]
    for ax, (title, data_y, model_y, ylim, ylabel) in zip(axes, panels):
        ax.plot(x, data_y.reindex(x), marker="o", linewidth=2.0, label="PSID data")
        ax.plot(x, model_y.reindex(x), marker="s", linewidth=2.0, label="Model age 45")
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xlabel("Income quintile")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_wealth_plot(psid: pd.DataFrame, model_wealth: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.7))
    x = np.arange(1, 6)
    data_parent = psid_slice(psid, "total_nw_k", "all", "parent_by_45").set_index("bin")
    data_childless = psid_slice(psid, "total_nw_k", "all", "childless_by_45").set_index("bin")
    data_children = psid_slice(psid, "total_nw_k", "all", "children_43_50").set_index("bin")
    model = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'all'").set_index("bin")
    panels = [
        ("Parent by 45 / parent share", data_parent["mean_outcome"] * 100, model["parent_share"] * 100, (0, 105), "Percent"),
        ("Childless by 45 / childless share", data_childless["mean_outcome"] * 100, model["childless_share"] * 100, (0, 45), "Percent"),
        ("Completed children", data_children["mean_outcome"], 2.0 * model["mean_completed_fertility"], (0, 2.8), "Count"),
    ]
    for ax, (title, data_y, model_y, ylim, ylabel) in zip(axes, panels):
        ax.plot(x, data_y.reindex(x), marker="o", linewidth=2.0, label="PSID data")
        ax.plot(x, model_y.reindex(x), marker="s", linewidth=2.0, label="Model age 45")
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xlabel("Wealth quintile")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make_tenure_wealth_plot(data: pd.DataFrame, model_wealth: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    x = np.arange(1, 6)
    stayed = data[data["transition"].eq("Stayed renter")].set_index("wealth_q")["childless_by_45"] * 100
    became = data[data["transition"].eq("Became owner")].set_index("wealth_q")["childless_by_45"] * 100
    renters = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'renters'").set_index("bin")["childless_share"] * 100
    owners = model_wealth.query("age_bin == 'age45' and concept == 'sale_net' and group == 'owners'").set_index("bin")["childless_share"] * 100
    ax.plot(x, stayed.reindex(x), marker="o", linewidth=2.0, label="PSID: stayed renter")
    ax.plot(x, became.reindex(x), marker="o", linewidth=2.0, label="PSID: became owner")
    ax.plot(x, renters.reindex(x), marker="s", linewidth=2.0, label="Model: renter at 45")
    ax.plot(x, owners.reindex(x), marker="s", linewidth=2.0, label="Model: owner at 45")
    ax.set_xticks(x)
    ax.set_xlabel("Wealth quintile")
    ax.set_ylabel("Childless share (percent)")
    ax.set_ylim(0, 75)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def build_summary(data: dict) -> pd.DataFrame:
    rows = []
    targets = data["targets"]
    moments = data["best_record"]["moments"]
    for key in ["tfr", "childless_rate", "mean_age_first_birth", "tfr_gradient"]:
        rows.append(summary_row("Target", key, targets[key], moments[key]))

    mms_loc = data["mms_location"].assign(location=data["mms_location"]["mms_location"].map({"center": "Center", "periphery": "Periphery"}))
    model = data["model_fertility"]
    for loc in ["Center", "Periphery"]:
        data_val = float(mms_loc[mms_loc["location"].eq(loc)].iloc[0]["has_children_rate"])
        model_val = model_value(model, "location", "ages22_45", loc, "All", "current_parent_share")
        rows.append(summary_row("Location current parent", loc, data_val, model_val))

    data_tenure = tenure_parent_share(data["mms_tenure"], location=None)
    for ten in ["Renter", "Owner"]:
        model_val = model_value(model, "tenure", "ages22_45", "All", ten, "current_parent_share")
        rows.append(summary_row("Tenure current parent", ten, data_tenure[ten], model_val))

    psid = data["psid_master"]
    model_income = data["model_income"].query("age_bin == 'age45' and concept == 'current_income' and group == 'all'").set_index("bin")
    d_inc = psid_slice(psid, "income_k", "all", "parent_by_45").set_index("bin")
    rows.append(summary_row("Income parent Q1", "parent_by_45", d_inc.loc[1, "mean_outcome"], np.nan))
    rows.append(summary_row("Income parent Q5", "parent_by_45", d_inc.loc[5, "mean_outcome"], model_income.loc[5, "parent_share"]))

    model_wealth = data["model_wealth"].query("age_bin == 'age45' and concept == 'sale_net' and group == 'all'").set_index("bin")
    d_w = psid_slice(psid, "total_nw_k", "all", "parent_by_45").set_index("bin")
    for q in [1, 5]:
        rows.append(summary_row("Wealth parent", f"Q{q}", d_w.loc[q, "mean_outcome"], model_wealth.loc[q, "parent_share"]))
    return pd.DataFrame(rows)


def summary_row(domain: str, statistic: str, data_value: float, model_value: float) -> dict[str, float | str]:
    return {
        "domain": domain,
        "statistic": statistic,
        "data_or_target": data_value,
        "model": model_value,
        "model_minus_data": model_value - data_value if np.isfinite(model_value) else np.nan,
    }


def build_tex(data: dict, figures: dict[str, Path], summary: pd.DataFrame) -> str:
    meta = dict(zip(data["model_metadata"]["key"], data["model_metadata"]["value"]))
    calib = fertility_target_table(data)
    slice_table = summary_table(summary)
    return rf"""\documentclass[11pt]{{article}}
\usepackage[margin=0.82in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{float}}
\usepackage{{graphicx}}
\usepackage{{caption}}
\usepackage{{array}}
\usepackage{{hyperref}}

\title{{Fertility Slice Diagnostics: Data and Current Model}}
\author{{}}
\date{{Built May 8, 2026}}

\begin{{document}}
\maketitle

\noindent This note isolates fertility moments by age, location, tenure, income, and wealth. It uses the current Python direct-geometry benchmark, job {latex_escape(str(meta.get("job_id", "")))}, evaluation {latex_escape(str(meta.get("eval_id", "")))}, with the renewal-valve stationary closure. The cluster-record loss is {fmt_num(meta.get("loss"), 3)} and the final equilibrium error is {fmt_num(meta.get("final_eq_error"), 4)}.

\medskip
\noindent The objects are deliberately separated. MMS moments are current-household fertility objects: whether an adult currently has children in the household, by age, location, and tenure. PSID moments are completed or late-life fertility objects: parenthood by 45, childlessness by 45, and children observed at ages 43--50, by young-adult income or wealth. The model analog is exact for location, tenure, and age in the sense that it averages the stationary distribution over the same state variables. For income, the comparison is intentionally diagnostic because the current model has deterministic lifecycle income and no persistent within-age income distribution.

{calib}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.94\textwidth]{{\detokenize{{{rel_to_latex(figures["age"])}}}}}
\caption{{Fertility by age and location. Data are MMS shares with children in the household; model is the stationary current-parent share.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.94\textwidth]{{\detokenize{{{rel_to_latex(figures["location"])}}}}}
\caption{{Fertility by location. Left panel uses reconstructed MMS TFR targets and the model TFR-equivalent completed fertility at age 45. Right panel uses current parenthood for adults up to age 45.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.96\textwidth]{{\detokenize{{{rel_to_latex(figures["tenure"])}}}}}
\caption{{Fertility by tenure. Data use MMS current parenthood by owner/renter status; model uses the stationary current-parent share over ages 22--45.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\textwidth]{{\detokenize{{{rel_to_latex(figures["income"])}}}}}
\caption{{Fertility by income. PSID bins use young-adult family income. The model line is the age-45 current-income bin. Because model income has no within-age idiosyncratic dispersion, all age-45 mass falls in one rank bin; this is a failure of the current model as an income-distribution diagnostic, not a substantive income-gradient result.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\textwidth]{{\detokenize{{{rel_to_latex(figures["wealth"])}}}}}
\caption{{Fertility by wealth. PSID bins use young-adult total net worth. The model uses sale-net wealth at age 45. Completed children in the model are reported in TFR-equivalent units, $2n$.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.78\textwidth]{{\detokenize{{{rel_to_latex(figures["tenure_wealth"])}}}}}
\caption{{Childlessness by wealth within tenure. The PSID classification follows young renters by later ownership status; the model classification is current tenure at age 45, so this is an imperfect but useful tenure--wealth diagnostic.}}
\end{{figure}}

{slice_table}

\section*{{Diagnosis}}

\noindent The model gets the broad sign of the location fertility gradient: periphery fertility is above center fertility. It still misses timing: the target mean age at first birth is 26.0 and the model is 33.5, so the age profile is too delayed.

\medskip
\noindent The tenure facts are sharper. In MMS data, current parenthood is higher among owners than renters, and this holds in both center and periphery cells. The model reproduces that ordering, but it overshoots owner parenthood and undershoots renter parenthood.

\medskip
\noindent The income comparison should not be oversold. The PSID income gradient is a within-cohort distributional fact. The current model has only deterministic lifecycle income. Pooling ages would mechanically create an income-fertility pattern through age composition; the clean age-45 plot shows the right limitation: there is no model income distribution to compare with PSID income quintiles.

\medskip
\noindent The wealth comparison is the most informative cross-sectional distributional cut. In the PSID, parenthood rises mildly through the middle of young-adult wealth and falls at the top. In the model, age-45 sale-net wealth produces a much steeper negative fertility gradient at the top, so wealth and tenure selection are doing too much work relative to the data.

\end{{document}}
"""


def fertility_target_table(data: dict) -> str:
    targets = data["targets"]
    moments = data["best_record"]["moments"]
    rows = []
    for key, label in [
        ("tfr", "TFR"),
        ("childless_rate", "Childless rate"),
        ("mean_age_first_birth", "Mean age first birth"),
        ("tfr_gradient", "TFR periphery minus center"),
    ]:
        rows.append([label, fmt_num(targets[key], 3), fmt_num(moments[key], 3)])
    return latex_table(
        "Fertility moments in the live benchmark calibration.",
        ["Moment", "Target", "Model"],
        rows,
        "lrr",
    )


def summary_table(summary: pd.DataFrame) -> str:
    rows = []
    keep = summary[
        summary["domain"].isin(
            ["Location current parent", "Tenure current parent", "Wealth parent", "Income parent Q5"]
        )
    ]
    for _, row in keep.iterrows():
        rows.append(
            [
                latex_escape(str(row["domain"])),
                latex_escape(str(row["statistic"])),
                fmt_num(row["data_or_target"], 3),
                fmt_num(row["model"], 3),
                fmt_num(row["model_minus_data"], 3),
            ]
        )
    return latex_table(
        "Selected slice gaps. Shares are in levels, not percent.",
        ["Slice", "Statistic", "Data", "Model", "Gap"],
        rows,
        "llrrr",
    )


def psid_slice(df: pd.DataFrame, concept: str, group: str, outcome: str) -> pd.DataFrame:
    out = df[
        df["wealth_concept"].eq(concept)
        & df["group"].eq(group)
        & df["outcome"].eq(outcome)
        & df["bin_type"].eq("quintile")
    ].sort_values("bin")
    if len(out) != 5:
        raise ValueError(f"Expected 5 rows for {concept}, {group}, {outcome}; found {len(out)}")
    return out


def tenure_parent_share(df: pd.DataFrame, location: str | None) -> dict[str, float]:
    sub = df if location is None else df[df["mms_location"].eq(location)]
    out: dict[str, float] = {}
    for tenure in ["Renter", "Owner"]:
        d = sub[sub["tenure"].eq(tenure)]
        total = float(d["pop_weight"].sum())
        parent = float(d[d["parent_status"].eq("With children")]["pop_weight"].sum())
        out[tenure] = parent / total if total > 0 else np.nan
    return out


def model_value(model: pd.DataFrame, slice_name: str, age_window: str, location: str, tenure: str, column: str) -> float:
    row = model[
        model["slice"].eq(slice_name)
        & model["age_window"].eq(age_window)
        & model["location"].eq(location)
        & model["tenure"].eq(tenure)
    ]
    if len(row) != 1:
        raise ValueError(f"Expected one model row for {slice_name}, {age_window}, {location}, {tenure}; found {len(row)}")
    return float(row.iloc[0][column])


def grouped_bar(ax, xlabels: list[str], series: list[tuple[str, list[float]]], ylabel: str, ylim: tuple[float, float] | None = None) -> None:
    x = np.arange(len(xlabels))
    width = 0.75 / max(len(series), 1)
    start = -0.5 * width * (len(series) - 1)
    for idx, (label, vals) in enumerate(series):
        vals = np.asarray(vals, dtype=float)
        pos = x + start + idx * width
        ax.bar(pos, vals, width, label=label)
        for xi, yi in zip(pos, vals):
            if np.isfinite(yi):
                ax.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7)
    ax.set_xticks(x, xlabels)
    if ylabel:
        ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def latex_table(caption: str, columns: list[str], rows: list[list[str]], align: str) -> str:
    body = "\n".join(" & ".join(row) + r" \\" for row in rows)
    header = " & ".join(columns) + r" \\"
    return rf"""\begin{{table}}[H]
\centering
\caption{{{caption}}}
\begin{{tabular}}{{{align}}}
\toprule
{header}
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def rel_to_latex(path: Path) -> str:
    return os.path.relpath(path, LATEX.parent)


def latex_escape(value: str) -> str:
    return (
        value.replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("$", r"\$")
        .replace("#", r"\#")
        .replace("_", r"\_")
        .replace("{", r"\{")
        .replace("}", r"\}")
    )


def fmt_num(value, digits: int = 3) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "--"
    if not np.isfinite(x):
        return "--"
    return f"{x:.{digits}f}"


if __name__ == "__main__":
    main()
