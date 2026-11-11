"""LaTeX table from summary_stats.csv written by sa_first_birth_summary_stats_20261111.do (Nov 11, 2026)."""
import csv, pathlib, subprocess, sys
out = pathlib.Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/empirical_psid_first_birth_summary_stats_20261111")
d = {(r["sample"], r["grp"], r["variable"]): r for r in csv.DictReader(open(out / "summary_stats.csv"))}
cols = [("A2h", "parents"), ("A2h", "controls"), ("A2h", "ctrl_age"), ("H", "parents"), ("H", "controls"), ("H", "ctrl_age")]
# Nov 11: six-column version with age-matched controls; the four-column files are kept.
f = lambda x, k=2: f"{float(x):,.{k}f}"
n = lambda x: f"{int(float(x)):,}"
def row(label, v, k=2, sd=True):
    m = " & ".join(f(d[(s, g, v)]["mean"], k) for s, g in cols)
    lines = [f"{label} & {m} \\\\"]
    if sd:
        lines.append(" & " + " & ".join(f"({f(d[(s, g, v)]['sd'], k)})" for s, g in cols) + " \\\\")
    return lines
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Summary Statistics, First-Birth Event-Study Samples (PSID)}",
     r"\label{tab:app-first-birth-sumstats}", r"\small",
     r"\begin{tabular}{@{}l ccc ccc@{}}", r"\toprule",
     r" & \multicolumn{3}{c}{Adults who were reference} & \multicolumn{3}{c}{Women who are reference} \\",
     r" & \multicolumn{3}{c}{person or spouse two to three} & \multicolumn{3}{c}{person or spouse} \\",
     r" & \multicolumn{3}{c}{years before the birth} & \multicolumn{3}{c}{in each year} \\",
     r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
     r" & First-time & Childless & Childless, & First-time & Childless & Childless, \\",
     r" & parents & controls & parents' ages & parents & controls & parents' ages \\", r"\midrule"]
L += row("Age", "age", 1)
L += row("Woman", "woman", 2, sd=False)
L += row("Married or cohabiting", "married", 2, sd=False)
L += row("Years of education", "edu", 1)
L += row("Family income (thousands of 2022 dollars)", "inc", 1)
L.append("Family income, median & " + " & ".join(f(d[(s, g, 'inc')]["p50"], 1) for s, g in cols) + r" \\")
L += row("Rooms", "rooms", 2)
L += row("Owner", "own", 2, sd=False)
L += [r"\addlinespace", r"\multicolumn{7}{@{}l}{\textit{Observations}} \\"]
L.append("Persons & " + " & ".join(n(d[(s, g, 'age' if g == 'ctrl_age' else 'py_all')]["persons"]) for s, g in cols) + r" \\")
L.append("Person-years & " + " & ".join(n(d[(s, g, 'py_all')]["n"]) for s, g in cols) + r" \\")
pr = lambda s: n(d[(s, "parents", "age")]["persons"])
py = lambda s: n(d[(s, "parents", "py_ref")]["n"])
L.append(f"Parents observed in years $-3$ and $-2$ & {pr('A2h')} & & & {pr('H')} & & \\\\")
L.append(f"Person-years in years $-3$ and $-2$ & {py('A2h')} & & & {py('H')} & & \\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\par\smallskip", r"\begin{minipage}{0.97\linewidth}",
 r"\footnotesize Note: This table reports means, weighted by the PSID longitudinal individual weight, with",
 r"standard deviations in parentheses. The samples are those of the rooms regressions in",
 r"Table~\ref{tab:app-rooms-path}, after dropping cohorts without an observation in the reference",
 r"window and people observed only once. First-time parents are measured in years $-3$ and $-2$",
 r"before the first birth; childless controls are measured in all their person-years in the sample, and",
 r"in the third column of each sample only at ages between the 5th and 95th percentiles of the parents'",
 r"age in years $-3$ and $-2$ (20 to 35 in the first sample, 20 to 34 in the second). Person fixed effects",
 r"absorb level differences between the groups; the controls identify the age and survey-year effects.",
 r"Married or cohabiting equals one when the person is the reference person or spouse and the",
 r"reference person is married or permanently cohabiting; other household members are coded zero.",
 r"Family income is total family-unit income in the previous tax year, deflated to 2022 dollars;",
 r"it is not reported in the 1994 and 1995 waves. Owner is computed on observations with a reported",
 r"tenure of owning or renting. Persons and person-years cover the full sample period of each group.",
 r"\end{minipage}", r"\end{table}"]
(out / "first_birth_summary_stats_agematched.tex").write_text("\n".join(L) + "\n")
doc = ("\\documentclass[11pt]{article}\n\\usepackage[margin=0.8in]{geometry}\n\\usepackage{booktabs}\n"
       "\\begin{document}\n\\pagestyle{empty}\n\\input{first_birth_summary_stats_agematched.tex}\n\\end{document}\n")
(out / "first_birth_summary_stats_agematched_standalone.tex").write_text(doc)
subprocess.run(["latexmk", "-pdf", "-quiet", "first_birth_summary_stats_agematched_standalone.tex"], cwd=out, check=True)
subprocess.run(["latexmk", "-c", "first_birth_summary_stats_agematched_standalone.tex"], cwd=out, check=True)
