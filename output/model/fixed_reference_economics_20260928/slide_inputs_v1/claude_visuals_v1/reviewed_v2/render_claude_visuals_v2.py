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
    ("frozen_reference", "Baseline"),
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
    "Dashed line = Baseline (frozen reference, original credit limits). All values are recomputed-cohort\n"
    "outcomes at fixed household preferences; \"fixed prices\" and \"stationary GE\" are not a transition path.",
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
plotted_data_path = os.path.join(OUT, "plotted_data.csv")
with open(plotted_data_path, "w") as f:
    f.write("figure,panel_or_cell,x_or_age,value\n")
    for age_label, row in pivot.iterrows():
        for cell in cell_order:
            v = row[cell]
            if pd.notna(v):
                f.write(f"birth_response_contributors,{cell},{age_label},{v}\n")
    for r in plotted_rows:
        f.write(f"fixed_price_vs_ge,{r['panel']},{r['case']},{r['value']}\n")

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
    "corrections_v1_to_v2": [
        "Removed invented 83.7% figure; single derived share only, computed from manifest total (all orders denominator).",
        "Labeled denominator as all initial households, not at-risk/within-cell.",
        "Added explicit four-year impact / fixed-price labels.",
        "Fig2: corrected y-axis label from pct.pts. to %; values are recomputed-cohort, not immediate impact.",
        "Fig2: replaced truncated bar chart with dot plot at true scale, three cases (baseline, credit fixed-price, credit stationary GE).",
        "Fig2: dropped 'transition' framing; first-birth age described as flow-weighted mean per recomputed population.",
    ],
}
with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)

print("done")
print(json.dumps(manifest, indent=2))
