#!/usr/bin/env python3
"""Render established transition panels from saved dated arrays; no model solves."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--native", type=Path, required=True, help="Saved candidate/latest_completed_full.json")
parser.add_argument("--receipt", type=Path, required=True, help="Saved run/latest_completed.json")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--targets", type=Path, help="Optional JSON mapping from birth-window end year to observed fertility")
parser.add_argument("--allow-incomplete", action="store_true", help="Label and plot saved paths that have not passed root/replay gates")
args = parser.parse_args()
NATIVE, REVIEW, HERE = args.native, args.receipt, args.output
HERE.mkdir(parents=True, exist_ok=True)
TARGETS = {int(k): float(v) for k, v in json.loads(args.targets.read_text()).items()} if args.targets else {}

native = json.loads(NATIVE.read_text())
review = json.loads(REVIEW.read_text())
rows, fertility = native["rows"], native["fertility"]
assert len(rows) == len(fertility) == review["horizon"]
if not args.allow_incomplete and not (review["root_pass"] and review["root_gate_passed"] and review["replay_pass"]):
    raise ValueError("Saved path has not passed root and replay gates; use --allow-incomplete for a diagnostic plot")
assert all(r["calendar_year"] == f["calendar_year"] for r, f in zip(rows, fertility))
assert rows == review["full_dated_rows"] and fertility == review["fertility_rows"]
if not review.get("two_shock_result", False):
    assert all(abs(r["psi_child"] - review["psi"]) < 1e-15 for r in rows)
x = np.array([r["calendar_year"] for r in rows], dtype=int)
assert np.array_equal(x, np.arange(x[0], x[0] + 4*len(x), 4))
get = lambda k: np.array([float(r[k]) for r in rows])
fy = np.array([float(f["period_tfr_topcode_adjusted"]) for f in fertility])
assert np.isfinite(fy).all()
for key in ("psi_child", "birth_children_topcode_adjusted", "birth_children", "adult_population", "asset_price", "renter_price", "owner_rate", "housing_demand"):
    if not np.isfinite(get(key)).all():
        raise ValueError(f"Nonfinite saved series: {key}")
if np.any(get("adult_population") <= 0):
    raise ValueError("Adult population must be positive for rates")
stage = "two-shock" if review.get("two_shock_result", False) else "first-shock"
gate_note = "root/replay pass" if review["root_pass"] and review["root_gate_passed"] and review["replay_pass"] else "root/replay incomplete"
horizon_note = "horizon comparison pending" if review.get("pending_horizon_comparison", False) else "horizon comparison status in receipt"
terminal_note = "terminal diagnostic pass" if review.get("terminal_diagnostic_pass", False) else "terminal certification pending"
scope_note = "second shock included" if review.get("two_shock_result", False) else "no second shock"
status_note = "Production-ready receipt" if review.get("production_ready", False) else "Diagnostic only"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

def layout():
    fig, axs = plt.subplots(2, 2, figsize=(11.8, 7.3))
    for ax in axs.flat:
        ax.grid(alpha=.18)
        ax.set_xlabel("Start of four-year period")
    return fig, axs

def finish(fig, stem):
    fig.suptitle(f"Saved {stage} transition: H{review['horizon']}", fontsize=15, y=.99)
    fig.text(.5, .925, f"First-shock preference = {review['psi']:.6g} | four-year dates from {x[0]} through {x[-1]}",
             ha="center", fontsize=10)
    fig.text(.5, .025, f"{status_note}: {scope_note}; {gate_note}; {horizon_note}; {terminal_note}.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .08, 1, .90))
    fig.savefig(HERE / f"{stem}.png", dpi=160, facecolor="white")
    fig.savefig(HERE / f"{stem}.pdf", facecolor="white")
    plt.close(fig)

fig, ax = layout()
ax[0, 0].plot(x, get("psi_child"), "o-", ms=3)
ax[0, 0].set(title="Preference input", ylabel="Child-preference coefficient")
ax[0, 1].plot(x+4, fy, "o-", ms=3, label="Model: topcode-adjusted period fertility")
if TARGETS:
    data_x = np.array(list(TARGETS)); data_y = np.array(list(TARGETS.values()))
    ax[0, 1].plot(data_x, data_y, "ks--", ms=4, label="US data: four-year average")
ax[0, 1].set(title="Fertility: model path and historical data", xlabel="End of four-year birth window", ylabel="Period fertility")
ax[0, 1].legend(fontsize=8)
ax[1, 0].plot(x, get("birth_children_topcode_adjusted"), "o-", ms=3, label="Topcode-adjusted")
ax[1, 0].plot(x, get("birth_children"), "s--", ms=3, label="Explicit-state")
ax[1, 0].set(title="Births over each four-year interval", ylabel="Births per initial model household")
ax[1, 0].legend(fontsize=8)
ax[1, 1].plot(x, 100*get("adult_population")/get("adult_population")[0], "o-", ms=3)
ax[1, 1].set(title=f"Adult household population ({x[0]} = 100)", ylabel="Index")
finish(fig, "fertility_demography")

fig, ax = layout()
for a, series, title, ylabel in zip(ax.flat,
        [get("asset_price"), get("renter_price"), 100*get("owner_rate"), get("housing_demand")/get("adult_population")],
        ["House asset price", "Rent per physical room", "Ownership", "Occupied physical rooms"],
        ["Model asset-price units", "Model rent units", "Percent of household heads", "Rooms per household head"]):
    a.plot(x, series, "o-", ms=3)
    a.set(title=title, ylabel=ylabel)
finish(fig, "housing")

# Match the established 2007–2063 single-panel fertility slide layout.
if x[0] == 2007 and x[-1] >= 2063:
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 12})
    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    use = x <= 2063
    ax.plot(x[use], fy[use], color="#d97815", lw=2.6, label=f"Model: {stage}")
    if TARGETS:
        dx = np.array(sorted(TARGETS), dtype=int) - 4
        dy = np.array([TARGETS[int(year)] for year in dx + 4])
        ax.plot(dx, dy, "s--", color="#245f99", lw=1.6, ms=5, label="US data")
    ax.set(xlim=(2003, 2065), xticks=[2007, 2023, 2043, 2063],
           xlabel="Start of four-year period", ylabel="Period fertility", title="Period fertility: 2007–2063")
    ax.grid(axis="y", alpha=.16)
    ax.legend(frameon=False)
    fig.text(.5, .02,
             f"{status_note}: first-shock preference {review['psi']:.6g}; H{review['horizon']}; {scope_note}; {horizon_note}; {terminal_note}.",
             ha="center", fontsize=10, color=".3")
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(HERE / "fertility_2007_2063.png", dpi=180, facecolor="white")
    fig.savefig(HERE / "fertility_2007_2063.pdf", facecolor="white")
    plt.close(fig)

with (HERE / "series.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["calendar_year", "preference", "period_fertility", "births_topcode_adjusted", "births_explicit", "adult_population", "asset_price", "rent", "ownership_percent", "rooms_per_head"])
    for i, year in enumerate(x):
        writer.writerow([year, rows[i]["psi_child"], fy[i], rows[i]["birth_children_topcode_adjusted"], rows[i]["birth_children"], rows[i]["adult_population"], rows[i]["asset_price"], rows[i]["renter_price"], 100*rows[i]["owner_rate"], rows[i]["housing_demand"]/rows[i]["adult_population"]])
receipt = {
    "classification": f"{status_note}: saved {stage} transition; {scope_note}; {gate_note}; {horizon_note}; {terminal_note}",
    "candidate": review.get("task_id"), "preference": review["psi"], "horizon": review["horizon"],
    "source_files": [str(NATIVE), str(REVIEW)],
    "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (NATIVE, REVIEW, *([args.targets] if args.targets else []))},
    "root_pass": review["root_pass"], "root_gate_passed": review["root_gate_passed"],
    "replay_pass": review["replay_pass"], "pending_horizon_comparison": review["pending_horizon_comparison"],
    "production_ready": review["production_ready"], "empirical_fitted": review["empirical_fitted"],
    "two_shock_result": review.get("two_shock_result", False), "terminal_diagnostic_pass": review.get("terminal_diagnostic_pass", False),
    "series_rows": len(x), "target_data": TARGETS,
    "layout_source": "code/model/tools/build_e5f_saved_transition_panels.py; values and captions replaced for c03_h24",
    "model_solves": 0,
}
(HERE / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
