"""Build the data-appendix tables and figure for the national ACS twin-like and
sibling-sex housing designs from the saved 18-row estimate table only.

Reads output/acs_fertility_iv/national_128g_results/national_18row_table.csv
(native units: rooms, bedrooms, ownership probability) and writes LaTeX table
fragments plus the reduced-form figure under output/acs_fertility_iv/data_appendix/.
No estimation is performed. Ownership coefficients and all first stages are
converted to percentage points for display only.

Run: python3 code/empirical/acs/kleven_pseudo/build_acs_iv_appendix_tables.py
"""
import csv
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
SRC = os.path.join(ROOT, "output", "acs_fertility_iv", "national_128g_results", "national_18row_table.csv")
OUT = os.path.join(ROOT, "output", "acs_fertility_iv", "data_appendix")
TAB = os.path.join(OUT, "tables")
FIG = os.path.join(OUT, "figures")

DESIGNS = [("Twin1", "A. Same-age oldest pair"), ("SameSex2", "B. Same-sex oldest pair")]
OUTCOMES = [("ROOMS_out", "Rooms", 1.0), ("BEDROOMS_out", "Bedrooms", 1.0), ("OWNERSHP_out", "Owner (pp)", 100.0)]
WINDOWS = [("pooled0_5", "0--5"), ("event3", "3"), ("event5", "5")]


def load():
    with open(SRC) as fh:
        return {(r["design"], r["outcome"]): r for r in csv.DictReader(fh)}


def f(x):
    return float(x)


def num(x, d):
    s = f"{x:,.{d}f}"
    return s.replace("-", "$-$")


def ar_cell(r, scale, d):
    lo, hi = r["ar_summary_lower"], r["ar_summary_upper"]
    if lo in ("", "NA") or hi in ("", "NA"):
        return "none"
    lo, hi = f(lo) * scale, f(hi) * scale
    if lo == hi:
        return "\\{" + num(lo, d) + "\\}"
    return "[" + num(lo, d) + ", " + num(hi, d) + "]"


def rows_for(by, design, window):
    out = []
    for oc, label, scale in OUTCOMES:
        r = by[(f"{design}_{window}", oc)]
        d_out = 2 if scale == 100.0 else 3
        rf, rf_se = f(r["rf_coef"]) * scale, f(r["rf_se"]) * scale
        fs, fs_se = f(r["fs_coef"]) * 100, f(r["fs_se"]) * 100
        iv, iv_se = f(r["iv_coef"]) * scale, f(r["iv_se"]) * scale
        n, hh, zpos = int(r["n_usable_prefit"]), int(r["n_households_prefit"]), int(r["n_instrument_positive"])
        assert int(r["rf_nobs_fit"]) == n == int(r["fs_nobs_fit"]) == int(r["iv_nobs_fit"]), (design, window, oc)
        assert r["status"] == "full_fit" and int(r["ar_n_errors"]) == 0, (design, window, oc)
        ar_d = 0 if scale == 100.0 else 1
        out.append(
            f"{label} & {n:,} & {hh:,} & {zpos:,} & {num(rf, d_out)} & {num(fs, 2)} & "
            f"{f(r['first_stage_F']):,.0f} & {num(iv, d_out)} & {ar_cell(r, scale, ar_d)} \\\\\n"
            f" & & & & ({num(rf_se, d_out)}) & ({num(fs_se, 2)}) & & ({num(iv_se, d_out)}) & \\\\\n"
        )
    return "".join(out)


def table(by, windows, caption, label):
    head = (
        "\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{@{}lrrrrrrrl@{}}\n\\toprule\n"
        " & & & Instrument & Reduced & First stage & First-stage & & AR accepted \\\\\n"
        "Outcome & Mothers & Households & equal to one & form & (pp) & $F$ & 2SLS & grid values \\\\\n"
        "\\midrule\n"
    )
    body = ""
    for design, dlabel in DESIGNS:
        for w, wlabel in windows:
            title = f"{dlabel}, event age {wlabel}"
            body += f"\\multicolumn{{9}}{{@{{}}l}}{{\\textit{{{title}}}}} \\\\\n"
            body += rows_for(by, design, w)
            body += "\\addlinespace\n"
    return head + body.removesuffix("\\addlinespace\n") + "\\bottomrule\n\\end{tabular}\n"


def figure(by):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(8.0, 4.8), sharex=True)
    xs = [0, 1, 2]
    xt = ["0–5", "3", "5"]
    for i, (design, dlabel) in enumerate(DESIGNS):
        for j, (oc, olabel, scale) in enumerate(OUTCOMES):
            ax = axes[i, j]
            y, lo, hi = [], [], []
            for w, _ in WINDOWS:
                r = by[(f"{design}_{w}", oc)]
                y.append(f(r["rf_coef"]) * scale)
                lo.append(f(r["rf_ci_lower"]) * scale)
                hi.append(f(r["rf_ci_upper"]) * scale)
            ax.errorbar(xs, y, yerr=[[a - b for a, b in zip(y, lo)], [b - a for a, b in zip(y, hi)]],
                        fmt="o", color="black", ms=4, capsize=3, lw=1)
            ax.axhline(0, color="0.6", lw=0.8, ls="--")
            ax.set_xticks(xs)
            ax.set_xticklabels(xt)
            ax.set_xlim(-0.5, 2.5)
            ax.tick_params(labelsize=8)
            if i == 0:
                ax.set_title({"Rooms": "Rooms", "Bedrooms": "Bedrooms", "Owner (pp)": "Owner (percentage points)"}[olabel], fontsize=9)
            if j == 0:
                ax.set_ylabel(dlabel.split(". ")[1], fontsize=9)
            if i == 1:
                ax.set_xlabel("Event age (years)", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "acs_iv_reduced_forms.pdf"))
    fig.savefig(os.path.join(FIG, "acs_iv_reduced_forms.png"), dpi=200)


DIAG = os.path.join(ROOT, "output", "acs_fertility_iv", "samesex_diagnosis")
DIAG_CASES = os.path.join(DIAG, "samesex_diagnosis_all_cases.csv")
DIAG_CELLS = os.path.join(DIAG, "remote_collection", "outdir", "D_sex_cell_counts.csv")
SPECS = [("A", "samesex", "Baseline"),
         ("D", "samesex", "With first- and second-child sex"),
         ("E", "samesex", "Excluding tied second and third ages"),
         ("D", "bothboys", "Two boys vs.\\ mixed"),
         ("D", "bothgirls", "Two girls vs.\\ mixed"),
         ("D", "BB minus GG", "Two boys minus two girls")]


def load_diag():
    with open(DIAG_CASES) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        assert r["status"] == "full_fit" and r["n_warnings"] == "0", r
    return rows


def pick(rows, case, term, oc, age="0-5"):
    m = [r for r in rows if r["case"] == case and r["term"] == term and r["outcome"] == oc and r["event_age"] == age]
    assert len(m) == 1, (case, term, oc, age)
    return m[0]


def diag_spec_table(rows):
    s = ("\\setlength{\\tabcolsep}{5pt}\n\\begin{tabular}{@{}lrrrrr@{}}\n\\toprule\n"
         " & \\multicolumn{3}{c}{Reduced form} & First stage & \\\\\n\\cmidrule(lr){2-4}\n"
         "Specification & Rooms & Bedrooms & Owner (pp) & (pp) & Mothers \\\\\n\\midrule\n")
    for case, term, lab in SPECS:
        rr = [pick(rows, case, term, oc) for oc, _, _ in OUTCOMES]
        sc = [s_ for _, _, s_ in OUTCOMES]
        b = [num(f(r["rf_coef"]) * k, 2 if k == 100 else 4) for r, k in zip(rr, sc)]
        se = [num(f(r["rf_se"]) * k, 2 if k == 100 else 4) for r, k in zip(rr, sc)]
        fs, fse = num(f(rr[1]["fs_coef"]) * 100, 2), num(f(rr[1]["fs_se"]) * 100, 2)
        s += f"{lab} & {b[0]} & {b[1]} & {b[2]} & {fs} & {int(rr[1]['n_usable']):,} \\\\\n"
        s += f" & ({se[0]}) & ({se[1]}) & ({se[2]}) & ({fse}) & \\\\\n"
    return s + "\\bottomrule\n\\end{tabular}\n"


def diag_age_table(rows):
    s = ("\\setlength{\\tabcolsep}{5pt}\n\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule\n"
         "Age of second & First stage & First-stage & \\multicolumn{3}{c}{Reduced form} & \\\\\n\\cmidrule(lr){4-6}\n"
         "child (years) & (pp) & $F$ & Rooms & Bedrooms & Owner (pp) & Mothers \\\\\n\\midrule\n")
    for e in range(6):
        rr = [pick(rows, "B", "samesex", oc, str(e)) for oc, _, _ in OUTCOMES]
        b = [num(f(r["rf_coef"]) * k, 2 if k == 100 else 4) for r, (_, _, k) in zip(rr, OUTCOMES)]
        se = [num(f(r["rf_se"]) * k, 2 if k == 100 else 4) for r, (_, _, k) in zip(rr, OUTCOMES)]
        r1 = rr[1]
        s += (f"{e} & {num(f(r1['fs_coef'])*100, 2)} & {f(r1['first_stage_F']):,.1f} & {b[0]} & {b[1]} & {b[2]} & "
              f"{int(r1['n_usable']):,} \\\\\n")
        s += f" & ({num(f(r1['fs_se'])*100, 2)}) & & ({se[0]}) & ({se[1]}) & ({se[2]}) & \\\\\n"
    return s + "\\bottomrule\n\\end{tabular}\n"


def diag_cell_table():
    with open(DIAG_CELLS) as fh:
        cells = {r["cell"]: r for r in csv.DictReader(fh)}
    names = [("BB", "Two boys"), ("GG", "Two girls"), ("BG", "Boy, then girl"), ("GB", "Girl, then boy")]
    s = ("\\setlength{\\tabcolsep}{5pt}\n\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule\n"
         "Two oldest children & Mothers & Weighted & Three or more & Rooms & Bedrooms & Owner \\\\\n"
         " & & (millions) & children & & & \\\\\n\\midrule\n")
    for pan, suf in [("A. Weighted means", "weighted"), ("B. Unweighted means", "unweighted")]:
        s += f"\\multicolumn{{7}}{{@{{}}l}}{{\\textit{{{pan}}}}} \\\\\n"
        for k, lab in names:
            c = cells[k]
            s += (f"{lab} & {int(c['n_unweighted']):,} & {f(c['n_weighted'])/1e6:.2f} & {f(c['D_mean_' + suf]):.3f} & "
                  f"{f(c['ROOMS_out_Y_mean_' + suf]):.3f} & {f(c['BEDROOMS_out_Y_mean_' + suf]):.3f} & "
                  f"{f(c['OWNERSHP_out_Y_mean_' + suf]):.3f} \\\\\n")
        s += "\\addlinespace\n"
    return s.removesuffix("\\addlinespace\n") + "\\bottomrule\n\\end{tabular}\n"


def diag_figures(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Figure: first stage and reduced forms by age of the second child
    fig, ax = plt.subplots(1, 4, figsize=(9.0, 2.6))
    panels = [("fs", "BEDROOMS_out", 100.0, "First stage (pp)")] + \
             [("rf", oc, k, lab if k == 1 else "Owner (pp)") for oc, lab, k in OUTCOMES]
    for a, (kind, oc, k, title) in zip(ax, panels):
        rr = [pick(rows, "B", "samesex", oc, str(e)) for e in range(6)]
        y = [f(r[f"{kind}_coef"]) * k for r in rr]
        e = [1.959964 * f(r[f"{kind}_se"]) * k for r in rr]
        a.errorbar(range(6), y, yerr=e, fmt="o", color="black", ms=3.5, capsize=2.5, lw=1)
        a.axhline(0, color="0.6", lw=0.8, ls="--")
        a.set_xticks(range(6))
        a.set_title(title, fontsize=9)
        a.set_xlabel("Age of second child", fontsize=8)
        a.tick_params(labelsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "acs_iv_samesex_by_age.pdf"))
    fig.savefig(os.path.join(FIG, "acs_iv_samesex_by_age.png"), dpi=200)
    # Figure: pooled reduced forms across specifications
    fig, ax = plt.subplots(1, 3, figsize=(9.0, 2.9), sharey=True)
    labs = [lab.replace("\\ ", " ") for _, _, lab in SPECS]
    yy = list(range(len(SPECS)))[::-1]
    for a, (oc, lab, k) in zip(ax, OUTCOMES):
        rr = [pick(rows, c, t, oc) for c, t, _ in SPECS]
        x = [f(r["rf_coef"]) * k for r in rr]
        e = [1.959964 * f(r["rf_se"]) * k for r in rr]
        a.errorbar(x, yy, xerr=e, fmt="o", color="black", ms=3.5, capsize=2.5, lw=1)
        a.axvline(0, color="0.6", lw=0.8, ls="--")
        a.set_title(lab if k == 1 else "Owner (pp)", fontsize=9)
        a.tick_params(labelsize=7.5)
    ax[0].set_yticks(yy)
    ax[0].set_yticklabels(labs, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "acs_iv_samesex_specs.pdf"))
    fig.savefig(os.path.join(FIG, "acs_iv_samesex_specs.png"), dpi=200)


def companion_csv(by):
    """All 18 baseline rows: conventional 95% CIs and tested-grid AR summaries, display units."""
    path = os.path.join(OUT, "acs_iv_18row_estimates_companion.csv")
    cols = ["design", "window", "outcome", "units", "mothers", "households", "instrument_positive",
            "rf_coef", "rf_se", "rf_ci95_lower", "rf_ci95_upper", "fs_coef_pp", "fs_se_pp", "fs_ci95_lower_pp",
            "fs_ci95_upper_pp", "first_stage_F_cluster_wald", "iv_coef", "iv_se", "iv_ci95_lower", "iv_ci95_upper",
            "ar_grid", "ar_accepted_grid_min", "ar_accepted_grid_max", "ar_n_accepted_runs"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for design, _ in DESIGNS:
            for win, _ in WINDOWS:
                for oc, lab, k in OUTCOMES:
                    r = by[(f"{design}_{win}", oc)]
                    g = "-50:50 by 2 (pp)" if k == 100 else "-3:3 by 0.1"
                    lo = "" if r["ar_summary_lower"] in ("", "NA") else f(r["ar_summary_lower"]) * k
                    hi = "" if r["ar_summary_upper"] in ("", "NA") else f(r["ar_summary_upper"]) * k
                    fse = f(r["fs_se"]) * 100
                    w.writerow([design, win, oc, "pp" if k == 100 else "count", r["n_usable_prefit"],
                                r["n_households_prefit"], r["n_instrument_positive"],
                                f(r["rf_coef"]) * k, f(r["rf_se"]) * k, f(r["rf_ci_lower"]) * k, f(r["rf_ci_upper"]) * k,
                                f(r["fs_coef"]) * 100, fse, f(r["fs_coef"]) * 100 - 1.959964 * fse,
                                f(r["fs_coef"]) * 100 + 1.959964 * fse, r["first_stage_F"],
                                f(r["iv_coef"]) * k, f(r["iv_se"]) * k, f(r["iv_ci_lower"]) * k, f(r["iv_ci_upper"]) * k,
                                g, lo, hi, r["ar_n_components"]])


def main():
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    by = load()
    assert len(by) == 18
    with open(os.path.join(TAB, "tab_acs_iv_pooled.tex"), "w") as fh:
        fh.write(table(by, [WINDOWS[0]], "Pooled estimates, event ages 0--5", "tab:acsiv_pooled"))
    with open(os.path.join(TAB, "tab_acs_iv_ages.tex"), "w") as fh:
        fh.write(table(by, WINDOWS[1:], "Estimates at event ages 3 and 5", "tab:acsiv_ages"))
    figure(by)
    companion_csv(by)
    rows = load_diag()
    for name, body in [("tab_acs_iv_diag_specs.tex", diag_spec_table(rows)),
                       ("tab_acs_iv_diag_ages.tex", diag_age_table(rows)),
                       ("tab_acs_iv_diag_cells.tex", diag_cell_table())]:
        with open(os.path.join(TAB, name), "w") as fh:
            fh.write(body)
    diag_figures(rows)
    print("wrote", TAB, FIG)


if __name__ == "__main__":
    main()
