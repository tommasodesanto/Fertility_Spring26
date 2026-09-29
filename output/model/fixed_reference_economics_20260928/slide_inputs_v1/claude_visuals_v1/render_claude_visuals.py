"""Render two supplemental prototype figures from already-saved small CSVs.

Zero model/lifecycle solves. Inputs are the existing occupied-state birth
decomposition (credit_v1/summary_v1/impact_birth_decomposition.csv) and the
existing fixed-price-vs-GE outcome table (borrowing_comparison.csv). Both are
staged read-only in ./inputs/. Outputs: two PNG+PDF pairs in ./actual_output/.
"""
import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, "inputs")
OUT = os.path.join(HERE, "actual_output")
os.makedirs(OUT, exist_ok=True)

CAPTION = (
    "2007 stationary reference — block0506, September 28 verified export. "
    "Fixed preferences, earnings, entry endowments, fiscal and housing-supply "
    "primitives; not a general-equilibrium or transition result unless labeled GE."
)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Figure 1: who contributes to the credit-driven birth response
# ---------------------------------------------------------------------------
decomp_path = os.path.join(IN, "impact_birth_decomposition.csv")
decomp = pd.read_csv(decomp_path)

decomp["age_label"] = (
    decomp["age_left"].astype(int).astype(str) + "-" + decomp["age_right"].astype(int).astype(str)
)
decomp["cell"] = decomp["inherited_tenure"] + ", " + decomp["inherited_net_financial_wealth"] + " wealth"

first = decomp[decomp["birth_order"] == "first"].copy()
total_diff = decomp["birth_difference"].sum()
first_share = first["birth_difference"].sum() / total_diff

pivot = first.pivot_table(
    index="age_label", columns="cell", values="birth_difference", aggfunc="sum"
)
age_order = sorted(pivot.index, key=lambda s: int(s.split("-")[0]))
pivot = pivot.loc[age_order]

cell_order = [
    "renter, nonpositive wealth",
    "renter, positive wealth",
    "owner, nonpositive wealth",
    "owner, positive wealth",
]
cell_order = [c for c in cell_order if c in pivot.columns]
pivot = pivot[cell_order]

colors = {
    "renter, nonpositive wealth": "#1f77b4",
    "renter, positive wealth": "#89b8dd",
    "owner, nonpositive wealth": "#d95f02",
    "owner, positive wealth": "#f2b57f",
}

fig, ax = plt.subplots(figsize=(8.0, 5.0))
bottom = pd.Series(0.0, index=pivot.index)
for cell in cell_order:
    vals = pivot[cell].fillna(0.0) * 1000.0
    ax.bar(pivot.index, vals, bottom=bottom, label=cell, color=colors.get(cell))
    bottom = bottom + vals

ax.set_ylabel("Additional first births per 1,000 households\n(credit minus baseline, impact cohort)")
ax.set_xlabel("Inherited age at shock (years)")
ax.set_title("Who drives the credit-induced first-birth increase")
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.text(
    0.01, 0.01,
    "First births are 83.4% of the total impact birth increase, shown here (83.7% within this decomposition).\n"
    "Renters before the shock contribute 95.6% and nonpositive inherited net financial wealth 67.6% of the total increase (all orders, README).\n"
    + CAPTION,
    fontsize=7, va="bottom", ha="left", wrap=True,
)
fig.tight_layout(rect=(0, 0.14, 0.98, 1))
fig1_png = os.path.join(OUT, "birth_response_contributors.png")
fig1_pdf = os.path.join(OUT, "birth_response_contributors.pdf")
fig.savefig(fig1_png, dpi=200, bbox_inches="tight")
fig.savefig(fig1_pdf, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 2: fixed-price impact behavior vs closed stationary GE endpoint
# ---------------------------------------------------------------------------
cmp_path = os.path.join(IN, "borrowing_comparison.csv")
cmp_df = pd.read_csv(cmp_path).set_index("outcome")

rows = [
    ("Completed fertility", "Completed fertility", None),
    ("Homeownership (%)", "Homeownership (pct. pts.)", None),
    ("Mean first-birth age (years)", "Mean first-birth age (years)", None),
]
regimes = [
    ("matched_grid_baseline_fixed_prices", "Baseline credit,\nfixed price"),
    ("solvency_only_fixed_prices", "Solvency-only credit,\nfixed price"),
    ("solvency_only_ge", "Solvency-only credit,\nclosed GE"),
]

fig, axes = plt.subplots(1, 3, figsize=(12.5, 5.2))
for ax, (row_key, title, _) in zip(axes, rows):
    ref = cmp_df.loc[row_key, "frozen_reference"]
    vals = [cmp_df.loc[row_key, col] for col, _ in regimes]
    labels = [lab for _, lab in regimes]
    x = range(len(vals))
    colors_bar = ["#4c72b0", "#dd8452", "#55a868"]
    ax.bar(x, vals, color=colors_bar)
    ax.axhline(ref, color="black", linewidth=1.0, linestyle="--")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8, rotation=20, ha="right")
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ymin = min(vals + [ref]) * 0.985
    ymax = max(vals + [ref]) * 1.015
    ax.set_ylim(ymin, ymax)
axes[0].set_ylabel("Level")
fig.suptitle("Fixed-price credit response vs. closed stationary GE endpoint", fontsize=13)
fig.text(
    0.5, 0.005,
    "Dashed line = frozen reference (original borrowing limits). Fixed-price panels hold house price and rent fixed;\n"
    "closed GE lets price and population clear so completed fertility returns to replacement. " + CAPTION,
    fontsize=7, ha="center", va="bottom",
)
fig.tight_layout(rect=(0, 0.12, 1, 0.92))
fig2_png = os.path.join(OUT, "fixed_price_vs_ge.png")
fig2_pdf = os.path.join(OUT, "fixed_price_vs_ge.pdf")
fig.savefig(fig2_png, dpi=200, bbox_inches="tight")
fig.savefig(fig2_pdf, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
manifest = {
    "inputs": {
        "impact_birth_decomposition.csv": sha256(decomp_path),
        "borrowing_comparison.csv": sha256(cmp_path),
    },
    "outputs": {
        os.path.basename(p): sha256(p)
        for p in [fig1_png, fig1_pdf, fig2_png, fig2_pdf]
    },
    "first_order_share_of_total_impact_birth_increase": float(first_share),
    "zero_solves": True,
}
with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=2)

print("done")
print(json.dumps(manifest, indent=2))
