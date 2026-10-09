"""Slide figures from saved results only (no solves): fertility fit, full transition four-panel (H100 two-shock path at
the 14.402 base, psi 0.1311/0.110) and a simplified property-tax-by-supply-regime figure (policy battery h16).
Usage: plot_jmp_slides_transition_figures_20261009.py [OUTDIR]"""
import json, sys, glob, pickle
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "output/model/jmp_slides_transition_figures_20261009"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.fontsize": 10})
C, R, B = "#d97a1c", "#c0392b", "#2c5f99"
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(OUT / f"{name}.{ext}", dpi=170)
    plt.close(fig)

# ---------------- shared transition path
T = ROOT / "output/model/transition_h100_14402_20261007"
r = json.load(open(T / "pairs_h100r_20261007.json"))["0.1311_0.11"]
path = r["stage1"]["path"][:2] + r["stage2"]["path"]
yrs = [p["year"] for p in path]; tfr = [p["tfr"] for p in path]
Q0, xi = 0.77941391535061, 0.63
ep = r["stage2"]["endpoint"]
H_init = r["stage2"]["path"][0]["stock"] * (Q0 / r["stage2"]["path"][0]["price"]) ** xi
H_end = H_init * (ep["price"] / Q0) ** xi

# ---------------- housing-stock-fixed overlay (saved GE result, stock frozen at its 2023 level from 2023)
FX = ROOT / "output/model/fixed_supply_2023/h100_pair12_v2_20261008/supply_ge_result.json"
fx = json.load(open(FX))["last"]
fy = [x["calendar_year"] for x in fx["rows"]]
fx_tfr = [x["period_tfr_topcode_adjusted"] for x in fx["fertility"]]
fx_pop = [x["adult_population"] for x in fx["rows"]]
fx_stock = [x["housing_supply"] for x in fx["rows"]]
fx_price = [x["asset_price"] for x in fx["rows"]]
FXLAB1, FXLAB2 = "Model, housing stock fixed from 2023", "Housing stock fixed from 2023"

# ---------------- Figure 1: fertility fit, points at period midpoints
wdi = {int(x["date"]): float(x["value"]) for x in json.load(open(
    ROOT / "output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence/source/stock_forecast/us_period_fertility_wdi.json"))[1]}
k = [i for i, y in enumerate(yrs) if y <= 2063]
mx = [2005] + [yrs[i] + 2 for i in k]; my = [2.1] + [tfr[i] for i in k]
fig, ax = plt.subplots(figsize=(9, 5))
ax.axhline(2.1, color="grey", ls=":", lw=1)
ay = sorted(y for y in wdi if y >= 2000)
ax.plot(ay, [wdi[y] for y in ay], "-o", color=B, lw=1, ms=3, alpha=.8, label="U.S. data (annual)")
# four-year averages as used in the slide (period labelled t averages data years t+1..t+4), plotted at model midpoints
dp = [2007, 2011, 2015, 2019]
ax.plot([t + 2 for t in dp], [np.mean([wdi[y] for y in range(t + 1, t + 5)]) for t in dp], "s", color=B, ms=8, label="U.S. data (four-year averages)")
ax.plot(mx, my, "-o", color=C, lw=2.2, ms=5, label="Model")
# author (Oct 9): no fixed-stock overlay on the shock-fitting slide
ax.set(xlim=(2000, 2066), ylim=(1.5, 2.25), xlabel="Year", ylabel="Births per woman"); ax.grid(axis="y", color="#eee")
ax.legend(frameon=False, loc="lower right")
save(fig, "fertility_fit_2007_2063")

# ---------------- Figure 2: four-panel
idx = lambda v, v0: [100 * x / v0 for x in v]
panels = [("Period fertility", tfr, 2.1, 2.1, "Births per woman"),
          ("Household heads", idx([p["population"] for p in path], 1.0), 100, 100 * ep["population_scale"], "Initial steady state = 100"),
          ("Total housing", idx([p["stock"] for p in path], H_init), 100, 100 * H_end / H_init, "Initial steady state = 100"),
          ("House price", idx([p["price"] for p in path], Q0), 100, 100 * ep["price"] / Q0, "Initial steady state = 100")]
fv = {"Period fertility": fx_tfr, "Household heads": idx(fx_pop, 1.0), "Total housing": idx(fx_stock, H_init),
      "House price": idx(fx_price, Q0)}
# two builds for a two-step slide: baseline only, then with the fixed-stock overlay (same axes)
for with_fx, fname in ((False, "full_transition_four_panel"), (True, "full_transition_four_panel_fixed")):
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5)); Tl = yrs[-1]
    for ax, (title, v, init, new, yl) in zip(axes.ravel(), panels):
        ax.plot([2003, 2007], [init, init], color=C, lw=2.5); ax.plot([2007, 2007], [init, v[0]], color=C, lw=2.5); ax.plot(yrs, v, color=C, lw=2.5)
        ax.axhline(init, color="grey", ls=":", lw=1); ax.axvline(2007, color="grey", ls=":", lw=1)
        # invisible overlay in the first build keeps identical axis limits across the two steps
        ax.plot(fy, fv[title], color="#2a7f62", ls="--", lw=2.2, label=FXLAB2 if with_fx else None, alpha=1.0 if with_fx else 0.0)
        ax.axhline(new, color=R, ls="--", lw=1.2); ax.plot([Tl], [new], "D", color=R, ms=7, label="New steady state")
        ax.set(title=title, ylabel=yl, xlabel="Start of four-year period", xlim=(2003, Tl + 4), xticks=[2007, 2103, 2203, 2303, 2403]); ax.grid(axis="y", color="#eee")
    axes[0, 0].legend(frameon=False, loc="lower right")
    fig.tight_layout(); save(fig, fname)

# ---------------- Figure 3: property tax by housing-supply regime (h16 battery)
RUNS = Path.home() / "fertility_runs/policy_battery_h16_20261008/runs"
def L(n):
    s = pickle.load(open(sorted(glob.glob(f"{RUNS}/{n}/map_*/summary.pkl"))[-1], "rb")); rr = s["rows"]
    g = lambda k: np.array([x[k] for x in rr])
    return dict(y=g("calendar_year"), b=g("birth_children_topcode_adjusted"), q=g("asset_price"), r=g("renter_price"), n=g("adult_population"))
reg = [("elastic", "base_elastic", "Supply responds", "-", "#1f5fa8"), ("slow", "slow_base", "Slowly adjusting stock", "--", "#2a9d8f"),
       ("frozen", "frozen_base", "Fixed stock", "-", "#c2452d")]
# author (Oct 9): four panels (births, house price, rent, population)
fig, axs = plt.subplots(2, 2, figsize=(10, 6.4)); axs = axs.ravel()
keys = [("b", "Births (% change)"), ("q", "House price (% change)"), ("r", "Rent per room (% change)"), ("n", "Population (% change)")]
for a, base, lab, ls, c in reg:
    x = L(f"ptax_{a}"); bb = L(base); m = (x["y"] >= 2023) & (x["y"] <= 2063)
    for ax, (k, _) in zip(axs, keys):
        v = 100 * (x[k] / bb[k] - 1)
        ax.plot(x["y"][m], v[m], ls=ls, color=c, lw=2.2, label=lab)
        print(a, k, [round(float(t), 2) for t in v[m][[0, 2, 5]]])
for ax, (_, t) in zip(axs, keys):
    ax.axhline(0, color="0.5", lw=.7); ax.set_title(t, loc="left"); ax.set_xlabel("Year"); ax.set_xlim(2023, 2063); ax.grid(axis="y", color="#eee")
axs[0].legend(frameon=False, loc="best"); fig.tight_layout(); save(fig, "ptax_by_supply_regime")
print("years model", yrs[:6], "last", yrs[-1], "tfr", [round(t, 3) for t in tfr[:5]])
