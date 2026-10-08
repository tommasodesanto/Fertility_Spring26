#!/usr/bin/env python3
"""Build v8 deck LaTeX table fragments from the official saved-result reader (slot base_2007)."""
import datetime, hashlib, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "latex/JMP_slides/prod_tables"
OUT.mkdir(parents=True, exist_ok=True)

def read(cmd):
    r = subprocess.run(["python3", "code/model/tools/read_results.py", cmd, "base_2007"],
                       cwd=ROOT, capture_output=True, text=True, check=True)
    return json.loads(r.stdout)

fit, par, show = read("fit"), read("parameters"), read("show")
F = {r["moment"]: r for r in fit["rows"]}
P = {r["parameter"]: r for r in par["rows"]}
summ = fit["summary"]

def f(x, d): return f"{float(x):.{d}f}"
def pct(x, d): return f"{100 * float(x):.{d}f}"
def num(r, key): return float(r[key])

# (moment, label, kind, decimals); None = group break
SPEC = [
    ("initial_normalization", "Completed fertility (normalization)", "n", 2),
    ("cps_childlessness", r"Childless, women 40--44 (\%)", "p", 1),
    ("cps_exactly_one", r"Exactly one child among mothers, 40--44 (\%)", "p", 1),
    ("nchs_mean_age", "Mean age at first birth (years)", "n", 2),
    ("nchs_share30", r"First births at age 30+ (\%)", "p", 1),
    ("early_fertility", "Children ever born by age 25", "n", 3),
    None,
    ("mean_rooms", "Mean rooms", "n", 2),
    ("ownership_30_55", r"Ownership, heads 30--55 (\%)", "p", 1),
    ("first_birth_rooms", "Rooms added at first birth", "n", 3),
    ("family_rooms", "Rooms: 3+ versus 1--2 children", "n", 3),
    ("recent_parent_ownership", "Recent-parent ownership gap", "n", 3),
    None,
    ("wealth_earnings", "Wealth / annual earnings", "n", 2),
    ("bequest_wealth", r"Annual bequests / wealth (\%)", "p", 2),
    ("old_dispersion", "Wealth p90/median, old age", "n", 2),
]
names = [s[0] for s in SPEC if s]
if set(names) != set(F) or len(names) != len(F):
    raise SystemExit(f"moment mismatch: only in spec {set(names)-set(F)}, only in JSON {set(F)-set(names)}")

def label(s): return s[1] + (r"$^{v}$" if F[s[0]]["role"] == "validation" else "")
def val(r, key, s): return (pct if s[2] == "p" else f)(r[key], s[3])
def fmt_weight(r):
    w = r["weight"]
    if w == "": return "norm."
    w = float(w)
    if w == 0: return "val."
    return f"{w:,.0f}" if w >= 100 else f"{w:.2f}"

def table(cols, header, extra_cells, total=None):
    L = [rf"\begin{{tabularx}}{{\textwidth}}{{{cols}}}", r"\toprule", header + r" \\", r"\midrule"]
    for s in SPEC:
        if s is None:
            L.append(r"\addlinespace[4pt]"); continue
        r = F[s[0]]
        L.append(" & ".join([label(s), val(r, "target", s), val(r, "model", s)] + extra_cells(r)) + r" \\")
    if total is not None:
        L += [r"\midrule", total]
    L += [r"\bottomrule", r"\end{tabularx}"]
    return "\n".join(L) + "\n"

def w(name, text): (OUT / name).write_text(text); print(f"wrote {OUT / name}")

# 1. base point
w("base_point.tex", "\n".join([
    rf"\def\baseloss{{{summ['loss']:.2f}}}", rf"\def\basechain{{{int(summ['chain'])}}}",
    rf"\def\baseprice{{{summ['price']:.3f}}}", rf"\def\baseHzero{{{summ['H0']:.2f}}}"]) + "\n")

# 2-3. fit tables
w("fit.tex", table(r"@{}Xrr@{}", "Moment & Target & Model", lambda r: []))
contrib = [float(r["loss_contribution"]) for r in F.values() if r["loss_contribution"] != ""]
tot = sum(contrib)
assert abs(tot - float(summ["loss"])) < 1e-6, (tot, summ["loss"])
w("fit_full.tex", table(
    r"@{}Xrrrr@{}", "Moment & Target & Model & Weight & Loss",
    lambda r: [fmt_weight(r), f(r["loss_contribution"], 3) if r["loss_contribution"] != "" else ""],
    total=rf"Total & & & & {tot:.3f} \\"))

# 4. parameters
PSPEC = [
    ("beta_annual", r"$\beta_{\rm annual}$", "Annual discount factor"),
    ("kappa_fert", r"$\kappa_1$", "First-birth taste dispersion"),
    ("kappa_fert_continuation", r"$\kappa_C$", "Later-birth taste dispersion"),
    ("first_birth_fixed_cost", r"$\xi$", "First-birth utility cost"),
    ("psi_child", r"$\psi_0$", "Preference for children (2007)"),
    ("child_benefit_curvature", r"$\gamma$", "Curvature of the benefit of children"),
    ("h_P", r"$h_P$", "Parenthood space requirement (rooms)"),
    ("chi", r"$\chi$", "Owner housing-service multiplier"),
    ("tenure_choice_kappa", r"$\kappa_T$", "Tenure taste dispersion"),
    ("theta0", r"$\theta_0$", "Bequest-motive scale"),
    ("H0", r"$H_0$", "Housing supply scale (derived)"),
]
def need(k):
    if k not in P: raise SystemExit(f"missing parameter {k}")
    return P[k]
L = [r"\begin{tabularx}{\textwidth}{@{}lXrl@{}}", r"\toprule",
     r"Parameter & Economic interpretation & Value & Bounds \\", r"\midrule"]
for k, sym, desc in PSPEC:
    r = need(k)
    v = f(r["estimate"], 3) + (r"$^{*}$" if r["near_bound"] == "True" else "")
    b = "--" if k == "H0" else f"[{float(r['lower']):g}, {float(r['upper']):g}]"
    L.append(f"{sym} & {desc} & {v} & {b} \\\\")
L += [r"\bottomrule", r"\end{tabularx}"]
w("params.tex", "\n".join(L) + "\n")

# 5. external parameters
E = lambda k: float(need(k)["estimate"])
EXT = [
    ("Risk aversion $\\sigma$", f"{E('sigma'):g}", ""),
    ("Consumption share $\\alpha_0$", f(E("alpha_cons"), 3), "CEX 2019--23, childless renters"),
    ("Equivalence scale $(\\omega,\\,\\nu)$".replace(",\\,", ","), "$(0.7,\\,0.7)$", "Scholz et al.\\ (2006)"),
    ("Housing-supply elasticity $\\eta$", f(E("housing_supply_elasticity"), 2), "Baum-Snow \\& Han (2024)"),
    ("Annual real interest rate", pct(E("q_annual"), 1) + "\\%", "Greaney et al.\\ (2025)"),
    ("Payroll tax rate $\\tau^{w}$", pct(E("payroll_tax"), 1) + "\\%", "balances the pay-as-you-go pension"),
    ("Financed share $\\phi$", f(E("financed_share"), 2), "20\\% down payment"),
    ("Depreciation $\\delta_H$", pct(E("annual_depreciation"), 2) + "\\% per year", ""),
    ("Property tax $\\tau^p$", pct(E("annual_property_tax"), 2) + "\\% per year", "rebated equally"),
    ("Selling cost $\\psi^s$", pct(E("selling_cost"), 0) + "\\%", "Greaney et al.\\ (2025)"),
    ("Rental size cap $\\bar h^R$", f"{E('rental_cap'):g} rooms", "Greaney et al.\\ (2025)"),
    ("Bequest wealth shift $\\theta_1$", f(E("theta1"), 4), "fixed externally"),
    ("Earnings persistence, innovation s.d.", "0.735, 0.484", "four-year AR(1)"),
    ("Ages", "18, 66, 45", "entry, retirement, last fertile age"),
]
L = [r"\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.36\textwidth}"
     r">{\raggedright\arraybackslash}p{0.20\textwidth}>{\raggedright\arraybackslash}X@{}}",
     r"\toprule", r"Parameter & Value & Source \\", r"\midrule"]
L += [f"{a} & {b} & {c} \\\\" for a, b, c in EXT]
L += [r"\bottomrule", r"\end{tabularx}"]
w("external.tex", "\n".join(L) + "\n")

# 6. receipt
paths = [ROOT / p if not Path(p).is_absolute() else Path(p)
         for k in ("fit", "parameters") for p in show["readiness"][k]["paths"]]
(OUT / "receipt.json").write_text(json.dumps({
    "timestamp": datetime.datetime.now().astimezone().isoformat(),
    "reference": show["reference"], "status": show["status"], "epoch": show["epoch"],
    "summary": show["summary"], "csv_paths": [str(p) for p in paths],
    "csv_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}, indent=2) + "\n")
print(f"wrote {OUT / 'receipt.json'}")

# 7. policy numbers -> prod_tables/policy.tex, used in the deck as \pol{case.field} (default precision) or \pol{case.field.1} (1 dp).
# Source: the reconciliation at the current base (not yet a production readout slot); point POLICY_SRC at the new run after a re-run.
import csv
POLICY_SRC = ROOT / "output/model/reconciliation_14p40_20261007"
def sgn(x, d): return f"{x:+.{d}f}"
pol = {}
for f in ("credit/results.csv", "property_tax/results.csv"):
    for r in csv.DictReader(open(POLICY_SRC / f)):
        for fld, d in (("births_pct", 2), ("own_18_29_pp", 1), ("price_pct", 2)):
            v = float(r[fld]); pol[f"{r['case']}.{fld}"] = sgn(v, d); pol[f"{r['case']}.{fld}.1"] = sgn(v, 1)
ep = json.load(open(POLICY_SRC / "psi2023/results.json"))
b = ep["psi110_new_base"]
for arm in ("phi095", "ptax2"):
    a = ep[f"psi110_new_{arm}"]
    for fld, v in (("price_pct", 100 * (a["price"] / b["price"] - 1)), ("pop_pct", 100 * (a["pop"] / b["pop"] - 1)),
                   ("own_18_29_pp", 100 * (a["own_18_29"] - b["own_18_29"]))):
        pol[f"end_{arm}.{fld}"] = sgn(v, 2); pol[f"end_{arm}.{fld}.1"] = sgn(v, 1)
lines = [r"% generated by build_prod_tables.py from " + str(POLICY_SRC.relative_to(ROOT))]
lines += [r"\expandafter\def\csname pol:%s\endcsname{%s}" % kv for kv in sorted(pol.items())]
w("policy.tex", "\n".join(lines) + "\n")
