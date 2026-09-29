"""Reviewed v2: corrections to two prototype figures per lead review.

Zero model/lifecycle solves; reads the same two saved small CSVs staged in
./inputs/. Writes plotted_data.csv (audit trail), two PNG+PDF pairs, and
manifest.json (input/output SHA-256) to ./actual_output/.
"""
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, "inputs")
OUT = os.path.join(HERE, "actual_output")
os.makedirs(OUT, exist_ok=True)

REF_LABEL = "2007 stationary reference — block0506, September 28 verified export"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Figure 1: who contributes to the credit-driven birth response (corrected)
# ---------------------------------------------------------------------------
decomp_path = os.path.join(IN, "impact_birth_decomposition.csv")
decomp = pd.read_csv(decomp_path)

decomp["age_label"] = (
    decomp["age_left"].astype(int).astype(str) + "-" + decomp["age_right"].astype(int).astype(str)
)
decomp["cell"] = decomp["inherited_tenure"] + ", " + decomp["inherited_net_financial_wealth"] + " wealth"

first = decomp[decomp["birth_order"] == "first"].copy()
total_all_orders = decomp["birth_difference"].sum()
first_share_of_total = first["birth_difference"].sum() / total_all_orders

pivot = first.pivot_table(index="age_label", columns="cell", values="birth_difference", aggfunc="sum")
age_order = sorted(pivot.index, key=lambda s: int(s.split("-")[0]))
pivot = pivot.loc[age_order]
cell_order = [c for c in [
    "renter, nonpositive wealth", "renter, positive wealth",
    "owner, nonpositive wealth", "owner, positive wealth",
] if c in pivot.columns]
pivot = pivot[cell_order]

colors = {
    "renter, nonpositive wealth": "#1f77b4",
    "renter, positive wealth": "#89b8dd",
    "owner, nonpositive wealth": "#d95f02",
    "owner, positive wealth": "#f2b57f",
}

fig, ax = plt.subplots(figsize=(8.2, 5.2))
bottom = pd.Series(0.0, index=pivot.index)
for cell in cell_order:
    vals = pivot[cell].fillna(0.0) * 1000.0
    ax.bar(pivot.index, vals, bottom=bottom, label=cell, color=colors.get(cell))
    bottom = bottom + vals

ax.set_ylabel("Additional first births per 1,000 initial households\n(credit minus baseline, four-year impact; fixed prices)")
ax.set_xlabel("Inherited age at shock (years)")
ax.set_title("Who drives the credit-induced first-birth increase")
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.text(
    0.01, 0.10,
    f"First births are {first_share_of_total*100:.1f}% of the total impact birth increase (all birth orders).\n"
    "Denominator: all initial households (occupied-state contributions), not at-risk or within-cell shares.",
    fontsize=7.5, va="bottom",
)
fig.text(0.01, 0.015, REF_LABEL, fontsize=7.5, va="bottom")
fig.tight_layout(rect=(0, 0.16, 0.98, 1))
fig1_png = os.path.join(OUT, "birth_response_contributors.png")
fig1_pdf = os.path.join(OUT, "birth_response_contributors.pdf")
fig.savefig(fig1_png, dpi=200, bbox_inches="tight")
fig.savefig(fig1_pdf, bbox_inches="tight")
plt.close(fig)

pivot_out = pivot.copy()
pivot_out.insert(0, "age_label", pivot_out.index)
pivot_out["first_share_of_total_all_orders"] = first_share_of_total

# ---------------------------------------------------------------------------
# Figure 2: three-case dot plot, corrected labels (recomputed-cohort values)
# ---------------------------------------------------------------------------
cmp_path = os.path.join(IN, "borrowing_comparison.csv")
cmp_df = pd.read_csv(cmp_path).set_index("outcome")

rows = [
    ("Completed fertility", "Completed fertility"),
    ("Homeownership (%)", "Homeownership (%)"),
    ("Mean first-birth age (years)", "Mean first-birth age (years)"),
]
cases = [
    ("matched_grid_baseline_fixed_prices", "Baseline (matched grid)"),
    ("solvency_only_fixed_prices", "Credit: fixed prices"),
    ("solvency_only_ge", "Credit: stationary GE"),
]

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.6))
plotted_rows = []
for ax, (row_key, title) in zip(axes, rows):
    vals = [cmp_df.loc[row_key, col] for col, _ in cases]
    labels = [lab for _, lab in cases]
    y = list(range(len(vals)))[::-1]
    ax.scatter(vals, y, s=90, color=["#4c72b0", "#dd8452", "#55a868"], zorder=3)
    ax.axvline(vals[0], color="black", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    span = max(vals) - min(vals)
    pad = max(span * 0.35, 1e-6)
    ax.set_xlim(min(vals) - pad, max(vals) + pad)
    for v, yy, lab in zip(vals, y, labels):
        plotted_rows.append({"panel": title, "case": lab, "value": v})
fig.suptitle("Recomputed-cohort outcomes: fixed prices vs. closed stationary GE", fontsize=12)
fig.text(
    0.5, 0.03,
    "Dashed line = Baseline (matched 262-node grid, original credit limits; not the original 160-node frozen\n"
    "reference, which differs numerically). All values are recomputed-cohort outcomes; not a transition path.",
    fontsize=7.5, ha="center", va="bottom",
)
fig.text(0.5, -0.005, REF_LABEL, fontsize=7.5, ha="center", va="bottom")
fig.tight_layout(rect=(0, 0.14, 1, 0.90))
fig2_png = os.path.join(OUT, "fixed_price_vs_ge.png")
fig2_pdf = os.path.join(OUT, "fixed_price_vs_ge.pdf")
fig.savefig(fig2_png, dpi=200, bbox_inches="tight")
fig.savefig(fig2_pdf, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
import csv

plotted_data_path = os.path.join(OUT, "plotted_data.csv")
fig1_unit = "additional_first_births_per_1000_initial_households"
fig2_units = {
    "Completed fertility": "children",
    "Homeownership (%)": "percent",
    "Mean first-birth age (years)": "years",
}
with open(plotted_data_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["figure", "panel_or_cell", "x_or_age", "value", "unit"])
    w.writeheader()
    for age_label, row in pivot.iterrows():
        for cell in cell_order:
            v = row[cell]
            if pd.notna(v):
                w.writerow({
                    "figure": "birth_response_contributors",
                    "panel_or_cell": cell,
                    "x_or_age": age_label,
                    "value": float(v) * 1000.0,
                    "unit": fig1_unit,
                })
    for r in plotted_rows:
        w.writerow({
            "figure": "fixed_price_vs_ge",
            "panel_or_cell": r["panel"],
            "x_or_age": r["case"],
            "value": r["value"],
            "unit": fig2_units.get(r["panel"], ""),
        })

# Validate: re-read with csv, check field counts and compare values.
with open(plotted_data_path, newline="") as f:
    reread = list(csv.DictReader(f))
assert all(len(row) == 5 for row in reread), "malformed row(s) in plotted_data.csv"
n_expected = sum(pd.notna(pivot[c]).sum() for c in cell_order) + len(plotted_rows)
assert len(reread) == n_expected, f"row count mismatch: {len(reread)} vs {n_expected}"

manifest = {
    "inputs": {
        "impact_birth_decomposition.csv": sha256(decomp_path),
        "borrowing_comparison.csv": sha256(cmp_path),
    },
    "outputs": {
        os.path.basename(p): sha256(p)
        for p in [fig1_png, fig1_pdf, fig2_png, fig2_pdf, plotted_data_path]
    },
    "first_birth_share_of_total_impact_birth_increase": float(first_share_of_total),
    "zero_solves": True,
    "corrections_v2_to_v3": [
        "Fig2 baseline case corrected from frozen_reference back to matched_grid_baseline_fixed_prices (262-node grid), matching v1 and the lead's 'exact same three cases/values' instruction; label now 'Baseline (matched grid)' with a footer noting it differs numerically from the original 160-node frozen reference.",
        "plotted_data.csv rewritten with csv.DictWriter (proper escaping of comma-containing cell names), explicit unit column, fig1 values in additional-births-per-1000-initial-households, fig2 values in children/percent/years; re-read and validated (field counts, row counts) after writing.",
    ],
}
with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)

print("done")
print(json.dumps(manifest, indent=2))
