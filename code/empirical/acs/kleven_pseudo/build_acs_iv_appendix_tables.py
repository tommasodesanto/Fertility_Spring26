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
    print("wrote", TAB, FIG)


if __name__ == "__main__":
    main()
