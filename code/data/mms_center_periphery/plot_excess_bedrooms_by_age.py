"""Age-conditional spare-space measure over time, to net out mechanical
cohort aging from the raw 'intergenerational allocation' series.

excess_bedrooms = bedrooms - household size, owner-occupied heads.
If a given age band holds more spare rooms over time, that is real
intensification net of cohort flow.
"""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt

OUTDIR = Path(__file__).resolve().parent / "output_big_homes_timeseries"

# ---- load age x year ----
by_age = defaultdict(dict)
with open(OUTDIR / "excess_bedrooms_by_age_year.csv") as f:
    for r in csv.DictReader(f):
        by_age[int(r["age5"])][int(r["year"])] = float(r["excess_bedrooms"])

# ---- load overall (old vs young) ----
by_old = defaultdict(dict)
with open(OUTDIR / "excess_bedrooms_overall_by_year.csv") as f:
    for r in csv.DictReader(f):
        by_old[int(r["year"])][int(r["old"])] = float(r["excess_bedrooms"])

ages = sorted(by_age)
years = sorted(by_old)
y0, y1 = years[0], years[-1]

lines = []
lines.append("Excess bedrooms (rooms - household size), owner-occupied heads")
lines.append(f"Age profile, {y0} vs {y1}:")
lines.append(f"{'age':>5} {str(y0):>8} {str(y1):>8} {'delta':>8}")
for a in ages:
    d = by_age[a]
    if y0 in d and y1 in d:
        lines.append(f"{a:>5} {d[y0]:>8.3f} {d[y1]:>8.3f} {d[y1]-d[y0]:>+8.3f}")
lines.append("")
lines.append(f"Overall owner-head excess bedrooms by year ({y0}-{y1}):")
lines.append(f"{'year':>5} {'<60':>8} {'60+':>8}")
for yr in years:
    lines.append(f"{yr:>5} {by_old[yr].get(0,float('nan')):>8.3f} {by_old[yr].get(1,float('nan')):>8.3f}")
summary = "\n".join(lines)
(OUTDIR / "excess_bedrooms_summary.txt").write_text(summary + "\n")
print(summary)

# ---- Plot 1: age profile, first vs last year ----
fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
ax.plot(ages, [by_age[a].get(y0) for a in ages], marker="o", lw=2.5, color="#5276A6", label=str(y0))
ax.plot(ages, [by_age[a].get(y1) for a in ages], marker="o", lw=2.5, color="#C73E3A", label=str(y1))
ax.set_xlabel("Head age (5-year band, lower edge)")
ax.set_ylabel("Excess bedrooms (rooms - household size)")
ax.grid(True, alpha=0.3); ax.legend(frameon=False)
fig.savefig(OUTDIR / "ts_excess_bedrooms_ageprofile.png", dpi=120, bbox_inches="tight")

# ---- Plot 2: time series for selected bands ----
fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
sel = {30: "#5276A6", 45: "#E1812C", 60: "#7B9C76", 70: "#C73E3A"}
for a, c in sel.items():
    if a in by_age:
        ys = sorted(by_age[a])
        ax.plot(ys, [by_age[a][y] for y in ys], marker="o", lw=2.5, color=c, label=f"age {a}-{a+4}")
ax.set_xlabel("Year"); ax.set_ylabel("Excess bedrooms (rooms - household size)")
ax.grid(True, alpha=0.3); ax.legend(frameon=False)
fig.savefig(OUTDIR / "ts_excess_bedrooms_byband.png", dpi=120, bbox_inches="tight")
print("\nwrote ts_excess_bedrooms_ageprofile.png and ts_excess_bedrooms_byband.png")
