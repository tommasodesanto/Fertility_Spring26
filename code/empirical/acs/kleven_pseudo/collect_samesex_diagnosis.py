"""Collect the national same-sex diagnosis receipts into a compact report.

Reads only the small files copied from the diagnosis output directory (case
receipts JSON, CSV summaries, identity JSON, logs) under
output/acs_fertility_iv/samesex_diagnosis/remote_collection/. No regression is
fitted. The only derived statistic is the BB-GG contrast, computed from each
fit's saved coefficients and clustered covariance (Var = V_bb + V_gg - 2 C).

Writes, in output/acs_fertility_iv/samesex_diagnosis/:
  samesex_diagnosis_all_cases.csv, summary.json, samesex_diagnosis.{pdf,png},
  report.md, report_tables.tex (+ report.pdf via pdflatex).
Run: code/model/.venv/bin/python code/empirical/acs/kleven_pseudo/collect_samesex_diagnosis.py
"""
import csv
import glob
import json
import math
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
OUT = os.path.join(ROOT, "output", "acs_fertility_iv", "samesex_diagnosis")
COL = os.path.join(OUT, "remote_collection")
REC = os.path.join(COL, "outdir", "case_receipts")
BASE = os.path.join(ROOT, "output", "acs_fertility_iv", "national_128g_results", "national_18row_table.csv")
OUTCOMES = ["ROOMS_out", "BEDROOMS_out", "OWNERSHP_out"]
LABEL = {"ROOMS_out": "Rooms", "BEDROOMS_out": "Bedrooms", "OWNERSHP_out": "Owner"}
Z = 1.959963984540054


def load(name):
    with open(os.path.join(REC, name + ".json")) as fh:
        return json.load(fh)


def row(case, spec, outcome, event_age, term, rf, rf_se, fs, fs_se, n, hh, nz, status, nwarn, extra=""):
    ci = lambda b, s: (b - Z * s, b + Z * s) if b is not None and s is not None else (None, None)
    rl, ru = ci(rf, rf_se)
    fl, fu = ci(fs, fs_se)
    return dict(case=case, spec=spec, outcome=outcome, event_age=event_age, term=term,
                rf_coef=rf, rf_se=rf_se, rf_ci_lower=rl, rf_ci_upper=ru,
                fs_coef=fs, fs_se=fs_se, fs_ci_lower=fl, fs_ci_upper=fu,
                first_stage_F=(fs / fs_se) ** 2 if fs is not None and fs_se else None,
                n_usable=n, n_households=hh, n_instrument_positive=nz,
                status=status, n_warnings=nwarn, note=extra)


def main():
    checks, failures, rows = {}, [], []
    files = sorted(glob.glob(os.path.join(REC, "*.json")))
    checks["n_case_receipts"] = len(files)
    if len(files) != 31:
        failures.append(f"expected 31 receipts, found {len(files)}")

    # A: exact baseline reproduction
    a = load("A_reproduction_gate")
    checks["A_pass"] = bool(a["pass"]) and all(a["gate"].values()) and all(a["checks"].values())
    if not checks["A_pass"]:
        failures.append("A reproduction gate")

    # C: coverage and audit
    cg, cs, ca = load("C_full_coverage_gate"), load("C_smoke"), load("C_full_audit")
    checks["C_full_coverage_ok"] = cg["keys_ok"] and cg["one_to_one_and_recomputed_equal"] and \
        cg["n_cached"] == cg["n_metadata_eligible"] == cg["n_matched_join_rows"] == 656986 and \
        cg["n_recompute_mismatches"] == 0
    checks["C_smoke_pass"] = bool(cs["smoke_pass"])
    checks["C_n_states"] = ca["n_states"]
    checks["C_pooled_audit"] = ca["pooled_audit"]
    for k, ok in [("C full coverage", checks["C_full_coverage_ok"]), ("C smoke", checks["C_smoke_pass"]),
                  ("C 51 states", ca["n_states"] == 51)]:
        if not ok:
            failures.append(k)

    # Baseline pooled rows (the national table; A certifies exact reproduction)
    with open(BASE) as fh:
        base = {r["outcome"]: r for r in csv.DictReader(fh) if r["design"] == "SameSex2_pooled0_5"}
    for oc in OUTCOMES:
        r = base[oc]
        rows.append(row("A", "baseline pooled 0-5", oc, "0-5", "samesex", float(r["rf_coef"]), float(r["rf_se"]),
                        float(r["fs_coef"]), float(r["fs_se"]), int(r["n_usable_prefit"]),
                        int(r["n_households_prefit"]), int(r["n_instrument_positive"]), r["status"], 0))

    statuses = {}
    # B: each event age
    for e in range(6):
        for oc in OUTCOMES:
            nm = f"B_event_age_{e}__{oc}"
            d = load(nm)
            statuses[nm] = d["status"]
            rows.append(row("B", f"event age {e}", oc, str(e), "samesex", d["rf_coef"], d["rf_se"], d["fs_coef"],
                            d["fs_se"], d["n_usable"], d["n_households"], d["n_instrument_positive"],
                            d["status"], len(d["warnings"])))
            assert d["rf_nobs"] == d["fs_nobs"] == d["n_usable"], nm
    # D additive Boy1/Boy2 controls
    for oc in OUTCOMES:
        nm = f"D_additive_sex_control__{oc}"
        d = load(nm)
        statuses[nm] = d["status"]
        assert "firstchildboy" in d["rf_b"] and "secondchildboy" in d["rf_b"] and d["constant_added_controls"] == []
        rows.append(row("D", "+ first-child boy, second-child boy", oc, "0-5", "samesex", d["rf_coef"], d["rf_se"],
                        d["fs_coef"], d["fs_se"], d["n_usable"], d["n_households"], d["n_instrument_positive"],
                        d["status"], len(d["warnings"]),
                        f"rf firstchildboy={d['rf_b']['firstchildboy']:.5f} secondchildboy={d['rf_b']['secondchildboy']:.5f}"))
    # D joint BB/GG vs mixed
    for oc in OUTCOMES:
        nm = f"D_joint_bb_gg__{oc}"
        d = load(nm)
        statuses[nm] = d["status"]
        rf, fs = d["rf"], d["fs"]
        assert rf["nobs"] == fs["nobs"] == d["n_usable"], nm
        common = dict(n=d["n_usable"], hh=d["n_households"])
        rows.append(row("D", "joint BB, GG vs mixed", oc, "0-5", "bothboys", rf["bb_coef"], rf["bb_se"], fs["bb_coef"],
                        fs["bb_se"], common["n"], common["hh"], d["n_bothboys"], d["status"], len(d["warnings"])))
        rows.append(row("D", "joint BB, GG vs mixed", oc, "0-5", "bothgirls", rf["gg_coef"], rf["gg_se"], fs["gg_coef"],
                        fs["gg_se"], common["n"], common["hh"], d["n_bothgirls"], d["status"], len(d["warnings"])))
        cse = lambda x: math.sqrt(x["bb_se"] ** 2 + x["gg_se"] ** 2 - 2 * x["bb_gg_cov"])
        rows.append(row("D", "joint BB, GG vs mixed", oc, "0-5", "BB minus GG", rf["bb_coef"] - rf["gg_coef"], cse(rf),
                        fs["bb_coef"] - fs["gg_coef"], cse(fs), common["n"], common["hh"], None, d["status"],
                        len(d["warnings"]), "SE from saved clustered V (bb, gg, cov)"))
    # E second/third tie exclusion
    for oc in OUTCOMES:
        nm = f"E_ambiguity_sensitivity__{oc}"
        d = load(nm)
        statuses[nm] = d["status"]
        rows.append(row("E", "exclude second = third child age", oc, "0-5", "samesex", d["rf_coef"], d["rf_se"],
                        d["fs_coef"], d["fs_se"], d["n_usable"], d["n_households"], d["n_instrument_positive"],
                        d["status"], len(d["warnings"]),
                        f"excluded {d['n_excluded_a2_eq_a3']} of {d['n_before']}"))
    checks["regression_case_statuses"] = statuses
    bad = [k for k, v in statuses.items() if v != "full_fit"]
    checks["n_regression_cases"] = len(statuses)
    checks["all_regression_cases_full_fit"] = not bad
    checks["total_warnings"] = sum(r["n_warnings"] for r in rows)
    failures += [f"status {k}" for k in bad]

    # Group means (weighted) from the sex-cell file
    with open(os.path.join(COL, "outdir", "D_sex_cell_counts.csv")) as fh:
        cells = list(csv.DictReader(fh))

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "samesex_diagnosis_all_cases.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    summary = dict(checks=checks, failures=failures, sex_cells=cells,
                   source=dict(remote_outdir="/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_samesex_diagnosis_20260922a",
                               job="18284841", identity=json.load(open(os.path.join(COL, "outdir", "samesex_roster_metadata_identity.json")))),
                   note="No IV or AR fits. Baseline pooled rows are the national table, certified identical by gate A.")
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    figure(rows)
    tables_tex(rows, cells)
    print("failures:", failures)
    print("checks:", {k: v for k, v in checks.items() if k not in ("regression_case_statuses", "C_pooled_audit")})


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 3, figsize=(9.5, 7.6))
    for j, oc in enumerate(OUTCOMES):
        sc = 100.0 if oc == "OWNERSHP_out" else 1.0
        b = [r for r in rows if r["case"] == "B" and r["outcome"] == oc]
        xs = [int(r["event_age"]) for r in b]
        for i, (k, lab) in enumerate([("rf", "Reduced form"), ("fs", "First stage (pp)")]):
            s = sc if k == "rf" else 100.0
            y = [r[f"{k}_coef"] * s for r in b]
            e = [Z * r[f"{k}_se"] * s for r in b]
            ax[i, j].errorbar(xs, y, yerr=e, fmt="o", color="black", ms=4, capsize=3, lw=1)
            ax[i, j].axhline(0, color="0.6", lw=0.8, ls="--")
            ax[i, j].set_xticks(range(6))
            ax[i, j].tick_params(labelsize=8)
            if j == 0:
                ax[i, j].set_ylabel(lab, fontsize=9)
            if i == 0:
                ax[i, j].set_title(LABEL[oc] + (" (pp)" if sc == 100 else ""), fontsize=9)
            if i == 1:
                ax[i, j].set_xlabel("Second child's age (years)", fontsize=8)
        specs = [("A", "samesex", "Baseline"), ("D", "samesex", "+Boy1, Boy2"), ("D", "bothboys", "BB vs mixed"),
                 ("D", "bothgirls", "GG vs mixed"), ("D", "BB minus GG", "BB - GG"), ("E", "samesex", "Excl. 2nd=3rd age")]
        y, e, labs = [], [], []
        for c, t, lab in specs:
            r = next(r for r in rows if r["case"] == c and r["term"] == t and r["outcome"] == oc)
            y.append(r["rf_coef"] * sc); e.append(Z * r["rf_se"] * sc); labs.append(lab)
        yy = list(range(len(specs)))[::-1]
        ax[2, j].errorbar(y, yy, xerr=e, fmt="o", color="black", ms=4, capsize=3, lw=1)
        ax[2, j].axvline(0, color="0.6", lw=0.8, ls="--")
        ax[2, j].set_yticks(yy)
        ax[2, j].set_yticklabels(labs if j == 0 else [""] * len(labs), fontsize=8)
        ax[2, j].tick_params(labelsize=8)
        ax[2, j].set_xlabel("Reduced form, ages 0-5 pooled", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "samesex_diagnosis.pdf"))
    fig.savefig(os.path.join(OUT, "samesex_diagnosis.png"), dpi=180)


def fmt(x, d):
    return "" if x is None else f"{x:,.{d}f}".replace("-", "$-$")


def tables_tex(rows, cells):
    L = []
    L.append("\\begin{tabular}{@{}llrrrrrrr@{}}\\toprule")
    L.append("Case & Outcome & RF & (SE) & FS pp & (SE) & $F$ & N & HH \\\\\\midrule")
    for r in rows:
        if r["case"] == "B":
            continue
        sc = 100 if r["outcome"] == "OWNERSHP_out" else 1
        d = 2 if sc == 100 else 4
        L.append(f"{r['case']}: {r['spec']} [{r['term']}] & {LABEL[r['outcome']]} & {fmt(r['rf_coef']*sc,d)} & ({fmt(r['rf_se']*sc,d)}) & "
                 f"{fmt(r['fs_coef']*100,2)} & ({fmt(r['fs_se']*100,2)}) & {fmt(r['first_stage_F'],0)} & {r['n_usable']:,} & {r['n_households']:,} \\\\")
    L.append("\\bottomrule\\end{tabular}")
    spec_tab = "\n".join(L).replace("_", "\\_")
    L = ["\\begin{tabular}{@{}lrrrrrrr@{}}\\toprule",
         "Outcome, age & RF & (SE) & FS pp & (SE) & $F$ & N & HH \\\\\\midrule"]
    for r in rows:
        if r["case"] != "B":
            continue
        sc = 100 if r["outcome"] == "OWNERSHP_out" else 1
        d = 2 if sc == 100 else 4
        L.append(f"{LABEL[r['outcome']]}, {r['event_age']} & {fmt(r['rf_coef']*sc,d)} & ({fmt(r['rf_se']*sc,d)}) & {fmt(r['fs_coef']*100,2)} & "
                 f"({fmt(r['fs_se']*100,2)}) & {fmt(r['first_stage_F'],1)} & {r['n_usable']:,} & {r['n_households']:,} \\\\")
    L.append("\\bottomrule\\end{tabular}")
    age_tab = "\n".join(L)
    L = ["\\begin{tabular}{@{}lrrrrrr@{}}\\toprule",
         "Cell & N & Weighted N & $D$ (w) & Rooms (w) & Bedrooms (w) & Owner (w) \\\\\\midrule"]
    for c in cells:
        L.append(f"{c['cell']} & {int(c['n_unweighted']):,} & {float(c['n_weighted']):,.0f} & {float(c['D_mean_weighted']):.4f} & "
                 f"{float(c['ROOMS_out_Y_mean_weighted']):.4f} & {float(c['BEDROOMS_out_Y_mean_weighted']):.4f} & {float(c['OWNERSHP_out_Y_mean_weighted']):.4f} \\\\")
    L.append("\\bottomrule\\end{tabular}")
    with open(os.path.join(OUT, "report_tables.tex"), "w") as fh:
        fh.write("\\newcommand{\\spectab}{" + spec_tab + "}\n\\newcommand{\\agetab}{" + age_tab +
                 "}\n\\newcommand{\\celltab}{" + "\n".join(L) + "}\n")


if __name__ == "__main__":
    main()
