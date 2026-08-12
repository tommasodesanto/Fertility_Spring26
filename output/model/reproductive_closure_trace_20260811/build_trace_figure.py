#!/usr/bin/env python3
"""Two-panel loci figure from the fixed-price trace (smoke grade)."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRATCH = Path(__file__).resolve().parent
ROWS = [
    row
    for row in csv.DictReader((SCRATCH / "trace_points.csv").open())
    if row["status"] == "ok"
]
ROWS.sort(key=lambda row: float(row["price"]))

rent = [float(row["rent"]) for row in ROWS]
b_over_e = [float(row["B_over_E"]) for row in ROWS]
scale = [float(row["clearing_scale_S"]) for row in ROWS]
nu = [float(row["lifetime_births_nu"]) for row in ROWS]

base = min(ROWS, key=lambda row: abs(float(row["price"]) - 0.776083))
base_rent = float(base["rent"])
base_be = float(base["B_over_E"])
base_scale = float(base["clearing_scale_S"])
leak = 1.0 - float(base["mature_entrants_B"]) / (
    0.5 * float(base["births_per_household"])
)
data_line = 0.5 * (1.0 - leak) * 1.918
ceiling = max(b_over_e)
print(
    f"baseline rent={base_rent:.4f} B/E={base_be:.4f} S={base_scale:.4f} "
    f"leak={leak:.4f} data_line={data_line:.4f} ceiling={ceiling:.4f} "
    f"nu(0)={max(nu):.4f} nu(base)={float(base['lifetime_births_nu']):.4f}"
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))

# Panel (a): reproduction against the rent.
ax1.plot(rent, b_over_e, color="black", lw=1.8)
ax1.axhline(1.0, color="black", ls="--", lw=1.0)
ax1.axhline(data_line, color="0.45", ls=":", lw=1.0)
ax1.plot([base_rent], [base_be], "o", color="black", ms=5)
ax1.annotate(
    "benchmark", (base_rent, base_be), textcoords="offset points",
    xytext=(6, -14), fontsize=9,
)
ax1.text(0.395, 1.005, "replacement", fontsize=9, ha="right", va="bottom")
ax1.text(0.395, data_line + 0.005, "US completed fertility", fontsize=9,
         ha="right", va="bottom", color="0.35")
ax1.set_xlabel("rental user cost per room $r$")
ax1.set_ylabel("reproduction $B(r)/E$")
ax1.set_xlim(0, 0.41)
ax1.set_ylim(0.79, 1.05)
ax1.set_title("(a) The reproduction schedule", fontsize=10)

# Panel (b): clearing locus, population relative to the benchmark.
rel = [s / base_scale for s in scale]
keep = [(x, y) for x, y in zip(rel, rent) if x <= 4.05]
ax2.plot([x for x, _ in keep], [y for _, y in keep], color="black", lw=1.8)
ax2.plot([1.0], [base_rent], "o", color="black", ms=5)
ax2.annotate(
    "benchmark", (1.0, base_rent), textcoords="offset points",
    xytext=(8, -4), fontsize=9,
)
ax2.set_xlabel("population supported, relative to benchmark $S/S_0$")
ax2.set_ylabel("rental user cost per room $r$")
ax2.set_xlim(0, 4.05)
ax2.set_title("(b) The clearing locus", fontsize=10)

for ax in (ax1, ax2):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(SCRATCH / "reproduction_and_clearing_loci.pdf")
fig.savefig(SCRATCH / "reproduction_and_clearing_loci.png", dpi=200)
print("figure written")
