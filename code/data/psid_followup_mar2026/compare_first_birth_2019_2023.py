"""Compare the A2h first-birth event study on PSID 1968-2019 with 1968-2023 (Oct 20 2026).

Reads the household-clustered coefficient files of both runs, writes a side-by-side CSV and one PDF
with four panels (rooms full sample, rooms common support, ownership, moved since last interview).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parent / "output"
OLD = OUT / "sa_first_birth_hhcluster_v1" / "coefficients_household.csv"
NEW = OUT / "sa_first_birth_hhcluster_2023" / "coefficients_household.csv"
DEST = OUT / "sa_first_birth_hhcluster_2023"

WINDOWS = ["Wleft", "Wm2", "Wm1", "ref", "Wp1", "Wp2", "Wp3", "Wp4", "Wp5", "Wp6", "Wright"]
LABELS = ["≤-8", "-7/-6", "-5/-4", "-3/-2", "-1/0", "1/2", "3/4", "5/6", "7/8", "9/10", "≥11"]
RUNS = [("rooms_full", "Rooms, full sample"), ("rooms_cs", "Rooms, common support"),
        ("own", "Owner (share)"), ("moved", "Moved since last interview (share)")]

old = pd.read_csv(OLD)
new = pd.read_csv(NEW)
m = old.merge(new, on=["run", "window"], suffixes=("_1968_2019", "_1968_2023"))
m["change"] = m["estimate_1968_2023"] - m["estimate_1968_2019"]
m.to_csv(DEST / "comparison_1968_2019_vs_2023.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
for ax, (run, title) in zip(axes.flat, RUNS):
    for df, lab, col, off in [(old, "1968–2019", "#7a7a7a", -0.12), (new, "1968–2023", "#1f5fa8", 0.12)]:
        d = df[df.run == run].set_index("window")
        x, y, e = [], [], []
        for i, w in enumerate(WINDOWS):
            if w == "ref":
                x.append(i + off); y.append(0.0); e.append(0.0)
            elif w in d.index:
                x.append(i + off); y.append(d.loc[w, "estimate"]); e.append(1.96 * d.loc[w, "se_household"])
        ax.errorbar(x, y, yerr=e, fmt="o", ms=4, capsize=2, color=col, label=lab)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(3, color="black", lw=0.4, ls=":")
    ax.set_xticks(range(len(WINDOWS)))
    ax.set_xticklabels(LABELS, fontsize=7)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Years relative to first birth", fontsize=8)
axes[0, 0].legend(fontsize=8, frameon=False)
fig.suptitle("First-birth event study (A2h), PSID 1968–2019 vs 1968–2023; 95% CI, household-clustered", fontsize=10)
fig.tight_layout()
fig.savefig(DEST / "first_birth_event_study_1968_2019_vs_2023.pdf")
print(m[m.window.isin(["Wp1", "Wp2", "Wp3", "Wp4"])].to_string(index=False))
