"""Table and PDF figures for the ACS first-birth headship designs (2026-11-12).

Reads the collected Torch outputs under
output/national_acs_comparison/headship_20261112/collected/ and writes the
table (CSV + markdown) and one PDF per comparison next to it. Read-only on
inputs; writes only new files in the headship folder.
"""
from pathlib import Path
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output/national_acs_comparison/headship_20261112"
COL = OUT / "collected"
DESIGNS = ["A", "B", "C", "D", "E", "A_age3045", "B_age3045", "A_fullpre", "B_fullpre"]
LABEL = {
    "A": "A baseline (all positions)",
    "B": "B householder/spouse/partner",
    "C": "C householder/spouse",
    "D": "D households only (no group quarters)",
    "E": "E donor-only household-position restriction",
    "A_age3045": "A, first birth at 30-45",
    "B_age3045": "B, first birth at 30-45",
    "A_fullpre": "A, anchors with donors at every pre-birth lag",
    "B_fullpre": "B, anchors with donors at every pre-birth lag",
}
OUTCOMES = [("rooms9", "Rooms (cap 9)", 1.0), ("bedrooms5", "Bedrooms (cap 5)", 1.0),
            ("ownership_lw", "Owner-occupied dwelling (pp)", 100.0)]


def load(design, name):
    p = COL / "estimates" / design / name
    return pd.read_csv(p) if p.exists() else None


def table():
    rows = []
    for d in DESIGNS:
        con, raw, cur, pc = (load(d, n) for n in ["national_contrasts.csv", "national_raw_baselines.csv",
                                                   "national_event_curves.csv", "pool_counts.csv"])
        if con is None:
            continue
        for o, lab, sc in OUTCOMES:
            c = con[(con.outcome == o) & (con.specification == "full")].iloc[0]
            r = raw[(raw.outcome == o) & (raw.event_time == -1)].iloc[0]
            n = cur[(cur.outcome == o) & (cur.specification == "full")].nobs.iloc[0]
            rows.append(dict(design=d, outcome=o, raw_level_minus1=r.weighted_mean * sc,
                             contrast_p3_m1=c.estimate * sc, se=c["std.error"] * sc, regression_n=int(n),
                             persons=int(pc.persons.iloc[0]), households=int(pc.households.iloc[0]),
                             anchors=int(pc.anchors.iloc[0])))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "headship_contrast_table.csv", index=False)
    return t


def curves(designs, fname, title):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    cols = ["#1f4e79", "#c0504d", "#4f8a3c", "#7f6084"]
    for j, (o, lab, sc) in enumerate(OUTCOMES):
        ax = axes[j]
        for i, d in enumerate(designs):
            cur = load(d, "national_event_curves.csv")
            if cur is None:
                continue
            z = cur[(cur.outcome == o) & (cur.specification == "full")].sort_values("event_time")
            x = z.event_time + (i - (len(designs) - 1) / 2) * 0.12
            ax.errorbar(x, z.estimate * sc, yerr=1.96 * z["std.error"] * sc, fmt="o-", ms=3, lw=1.2,
                        capsize=0, color=cols[i], label=LABEL[d])
        ax.axhline(0, color="grey", lw=0.6); ax.axvline(-0.5, color="grey", lw=0.6, ls=":")
        ax.set_title(lab, fontsize=10); ax.set_xlabel("Years relative to first birth")
        ax.set_xticks(range(-5, 11))
    axes[0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle(title + " (ACS 2005-19, women; reference -2; state, age, year FE; 95% CI, household-clustered)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / fname); plt.close(fig)


def positions():
    files = glob.glob(str(COL / "states" / "statefip_*" / "diag_positions.csv"))
    if not files:
        return None
    d = pd.concat([pd.read_csv(f) for f in files])
    g = d.groupby(["t", "role", "pos"], as_index=False)[["n", "w", "w_parent", "w_rooms", "wy_rooms", "w_own", "wy_own"]].sum()
    tot = g.groupby(["t", "role"]).w.transform("sum")
    g["share"] = g.w / tot
    g["rooms_mean"] = g.wy_rooms / g.w_rooms
    g["own_mean"] = g.wy_own / g.w_own
    g.to_csv(OUT / "headship_position_shares_by_event.csv", index=False)
    par = d.groupby(["t", "role"], as_index=False)[["w", "w_parent"]].sum()
    par["parent_coresident_share"] = par.w_parent / par.w
    par.to_csv(OUT / "headship_parent_coresidence_by_event.csv", index=False)
    mar = d.copy(); mar["head"] = mar.pos.isin(["householder", "spouse", "unmarried_partner"])
    mm = mar.groupby(["t", "married_flag", "head"], as_index=False).w.sum()
    mm["share"] = mm.w / mm.groupby(["t", "married_flag"]).w.transform("sum")
    mm[mm["head"]].to_csv(OUT / "headship_head_share_by_marital_status.csv", index=False)

    order = ["householder", "spouse", "unmarried_partner", "child_of_householder", "other_relative", "other_nonrelative", "group_quarters"]
    colors = ["#1f4e79", "#5b9bd5", "#9dc3e6", "#c0504d", "#f4b183", "#a5a5a5", "#404040"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    ax = axes[0]
    for p, c in zip(order, colors):
        z = g[g.pos == p].sort_values("t")
        ax.plot(z.t, z.share, "o-", ms=3, color=c, label=p.replace("_", " "))
    ax.axvline(-0.5, color="grey", lw=0.6, ls=":")
    ax.set_title("Household position, weighted share\n(t<0 matched childless donors, t>=0 mothers)", fontsize=9)
    ax.set_xlabel("Years relative to first birth"); ax.legend(fontsize=7, frameon=False)
    ax = axes[1]
    z = par.sort_values("t"); ax.plot(z.t, z.parent_coresident_share, "o-", color="#c0504d", ms=3)
    ax.axvline(-0.5, color="grey", lw=0.6, ls=":")
    ax.set_title("Living with own mother or father (MOMLOC/POPLOC)", fontsize=9); ax.set_xlabel("Years relative to first birth")
    ax = axes[2]
    for p, c in zip(["householder", "spouse", "unmarried_partner", "child_of_householder"], colors[:3] + [colors[3]]):
        z = g[g.pos == p].sort_values("t")
        ax.plot(z.t, z.rooms_mean, "o-", ms=3, color=c, label=p.replace("_", " "))
    ax.axvline(-0.5, color="grey", lw=0.6, ls=":")
    ax.set_title("Rooms (cap 9), weighted mean within position", fontsize=9); ax.set_xlabel("Years relative to first birth")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Baseline A matched panel: household position by event time (ACS 2005-19, women)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "headship_positions_by_event.pdf"); plt.close(fig)
    return g


def support():
    files = glob.glob(str(COL / "states" / "statefip_*" / "support_by_lag.csv"))
    if not files:
        return None
    s = pd.concat([pd.read_csv(f) for f in files])
    s = s.groupby(["design", "lag", "cat"], as_index=False)[["anchors", "anchor_wgt"]].sum()
    s.to_csv(OUT / "headship_support_by_lag.csv", index=False)
    return s


def markdown(t):
    fmt = lambda v, k=3: f"{v:.{k}f}"
    lines = ["| Design | Outcome | Raw level at -1 | +3 minus -1 (SE) | Regression N | Persons | Households | Anchors |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    for _, r in t.iterrows():
        lines.append(f"| {LABEL[r.design]} | {r.outcome} | {fmt(r.raw_level_minus1, 2)} | {fmt(r.contrast_p3_m1)} ({fmt(r.se)}) | "
                     f"{r.regression_n:,} | {r.persons:,} | {r.households:,} | {r.anchors:,} |")
    lines += ["", "PSID references (not like for like; sample and windows stated): women-only household sample "
              "(women reference person or spouse/partner each interview) rooms 1.025 (0.073), +3/+4 vs -3/-2; "
              "main sample (men and women, reference person or spouse 2-3 years before the birth) rooms 1.465 (0.050), "
              "+3/+4 vs -3/-2. Source: code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/summary.csv (H, A2h)."]
    (OUT / "headship_contrast_table.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    t = table()
    markdown(t)
    curves(["A", "B"], "headship_curves_A_vs_B.pdf", "A baseline vs B householder/spouse/partner")
    curves(["A", "E", "B"], "headship_curves_A_E_B.pdf", "A baseline, E donor-only restriction, B")
    curves(["A", "C", "D"], "headship_curves_A_C_D.pdf", "A baseline, C householder/spouse, D households only")
    curves(["B", "B_age3045", "B_fullpre"], "headship_curves_B_sensitivity.pdf", "B and its sensitivities")
    positions(); support()
    print(t.to_string())
