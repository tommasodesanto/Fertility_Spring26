"""Generate the CES-appendix tables of the JMP deck from saved CES calibration results.

Writes latex/JMP_slides/ces_tables/{fit.tex, params.tex, mechanisms.tex, numbers.tex, receipt.json}.
numbers.tex defines \cesLTVhundred, \cesLowerRisk and \cesCapTen (completed-fertility responses quoted in the frame text).
The frame's sentence on which rows the CES fits better or worse than the benchmark is written by hand; recheck it at a new leader.
Run from anywhere:  python3 latex/JMP_slides/build_ces_tables.py [--leader PATH] [--cd-table PATH]
                     [--ces-cells DIR] [--ces-rental PATH] [--cd-cells DIR]

Defaults point at the Oct 8 2026 final CES leader (Torch chain 2 case 0038_nm, loss 64.313, Mac repeat to 1e-13):
  leader      output/model/ces_calibration_20261008/readout_20261008/final/leader_best_so_far.json
              (target_fit: every scored moment with target and model value; parameters; H0, price, rebate)
  cd-table    output/model/ces_calibration_20261008/readout_20261008/table.json
              (column 7 of each row = the production Cobb-Douglas 14.402 point measured with the same observers)
  ces-cells   output/model/ces_calibration_20261008/start_diagnostics/final_leader/{base,LTV100,riskB}.json
  ces-rental  output/model/ces_calibration_20261008/postcal_rental_risk/final_leader/results.json  (cap6, cap10)
  cd-cells    output/model/ces_calibration_20261008/comparison_20261008/cd_cells/{base,LTV100,riskB,cap10}.json
To refresh at a new leader, pass --leader <its best_so_far.json> and the folders of the mechanism cells run at that
point; the mechanism rows must come from cells solved at the same parameter point as the fit table.
"""
import argparse, json, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CES = ROOT / "output/model/ces_calibration_20261008"
OUT = ROOT / "latex/JMP_slides/ces_tables"

# (moment key, row label, formatter)  in the order of the deck's main calibration table
def pct1(x): return f"{100*x:.1f}"
def pct2(x): return f"{100*x:.2f}"
def f2(x): return f"{x:.2f}"
def f3(x): return f"{x:.3f}"
FIT_ROWS = [
    ("cps_childlessness",            r"Childless, women 40--44 (\%)",                  pct1),
    ("cps_exactly_one",              r"One child among mothers, 40--44 (\%)",          pct1),
    ("three_plus_among_mothers",     r"Mothers with 3+ children, 40--44 (\%)",          pct1),
    ("nchs_mean_age",                r"Mean age at first birth (years)",               f2),
    ("early_fertility",              r"Children ever born by age 25",                  f3),
    None,
    ("mean_rooms",                   r"Mean occupied rooms",                           f2),
    ("ownership_30_55",              r"Ownership, heads 30--55 (\%)",                  pct1),
    ("renter_to_owner_4yr",          r"Renter to owner within 4 years (\%)",            pct1),
    ("rent_share_childless_renters", r"Rent share, childless renters (\%)",             pct1),
    ("first_birth_rooms",            r"First-birth room response",                     f3),
    ("family_rooms",                 r"Room gap, 3+ vs 1--2 children at home",         f3),
    None,
    ("wealth_earnings",              r"Wealth / annual gross earnings",                f2),
    ("bequest_wealth",               r"Annual bequests / wealth (\%)",                 pct2),
]
PARAM_ROWS = [  # (key, symbol, interpretation)  main-deck order first, then the CES block
    ("beta_annual",             r"$\beta_{\rm annual}$", "Annual discount factor"),
    ("kappa_fert",              r"$\kappa_1$",           "First-birth taste scale"),
    ("kappa_fert_continuation", r"$\kappa_C$",           "Later-birth taste scale"),
    ("chi",                     r"$\chi$",               "Owner housing-service multiplier"),
    ("H0",                      r"$H_0$",                "Housing supply scale"),
    ("theta0",                  r"$\theta_0$",           "Bequest-motive scale"),
    ("first_birth_fixed_cost",  r"$\xi$",                "First-birth utility cost"),
    ("psi_child",               r"$\psi_0$",             "Initial preference for children"),
    ("child_benefit_curvature", r"$\gamma$",             "Curvature of the benefit of children"),
    ("tenure_choice_kappa",     r"$\kappa_T$",           "Tenure taste scale"),
    ("omega",                   r"$\alpha_0$",           "Consumption weight in the composite"),
    ("eta1",                    r"$\eta_1$",             "Housing need of the first child"),
    ("eta2",                    r"$\eta_2$",             "Housing need of each later child"),
]

def pctchange(new, base):
    v = 100 * (new / base - 1)
    return f"${'+' if v >= 0 else '-'}{abs(v):.1f}\\%$"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leader", default=CES / "readout_20261008/final/leader_best_so_far.json")
    ap.add_argument("--cd-table", default=CES / "readout_20261008/table.json")
    ap.add_argument("--ces-cells", default=CES / "start_diagnostics/final_leader")
    ap.add_argument("--ces-rental", default=CES / "postcal_rental_risk/final_leader/results.json")
    ap.add_argument("--cd-cells", default=CES / "comparison_20261008/cd_cells")
    a = ap.parse_args()
    leader = json.load(open(a.leader))["best"]
    fit = {r["moment"]: r for r in leader["target_fit"]}
    cd = {row[0]: row[6] for row in json.load(open(a.cd_table))["rows"]}
    OUT.mkdir(exist_ok=True)

    lines = [r"\begin{tabularx}{\textwidth}{@{}Xrrr@{}}", r"\toprule",
             r"Moment & Target & CES & Cobb--Douglas \\", r"\midrule"]
    for row in FIT_ROWS:
        if row is None:
            lines.append(r"\addlinespace[3pt]"); continue
        key, label, fmt = row
        r = fit[key]
        assert r["role"] == "scored", key
        lines.append(f"{label} & {fmt(r['target'])} & {fmt(r['model'])} & {fmt(cd[key])} \\\\")
    lines += [r"\bottomrule", r"\end{tabularx}"]
    (OUT / "fit.tex").write_text("\n".join(lines) + "\n")

    vals = dict(leader["parameters"]); vals["H0"] = leader["H0"]
    lines = [r"\begin{tabular}{@{}lr@{}}", r"\toprule", r"Parameter & Value \\", r"\midrule"]
    for key, sym, text in PARAM_ROWS:   # interpretation kept as a LaTeX comment
        lines.append(f"{sym} & {vals[key]:.3f} \\\\ % {text}")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (OUT / "params.tex").write_text("\n".join(lines) + "\n")

    cell = lambda d, name: json.load(open(Path(d) / f"{name}.json"))["completed_fertility_46_49"]
    rental = json.load(open(a.ces_rental))
    ces = {"base": cell(a.ces_cells, "base"), "LTV100": cell(a.ces_cells, "LTV100"),
           "riskB": cell(a.ces_cells, "riskB"), "cap6": rental["cap6"]["F"], "cap10": rental["cap10"]["F"]}
    cdm = {k: cell(a.cd_cells, k) for k in ("base", "LTV100", "riskB", "cap10")}
    assert abs(ces["base"] - ces["cap6"]) < 1e-9, (ces["base"], ces["cap6"])
    rows = [(r"Financed share $80\%\to100\%$", pctchange(ces["LTV100"], ces["base"]), pctchange(cdm["LTV100"], cdm["base"])),
            (r"Lower earnings risk, mean-preserving", pctchange(ces["riskB"], ces["base"]), pctchange(cdm["riskB"], cdm["base"])),
            (r"Rental cap $6\to10$ rooms", pctchange(ces["cap10"], ces["cap6"]), pctchange(cdm["cap10"], cdm["base"]))]
    lines = [r"\begin{tabular}{@{}lrr@{}}", r"\toprule", r"Experiment & CES & Cobb--Douglas \\", r"\midrule"]
    lines += [f"{x} & {y} & {z} \\\\" for x, y, z in rows]
    lines += [r"\bottomrule", r"\end{tabular}"]
    (OUT / "mechanisms.tex").write_text("\n".join(lines) + "\n")
    (OUT / "numbers.tex").write_text(
        "% completed fertility at 46--49, CES at the calibrated point, fixed price and rebate (generated)\n"
        f"\\def\\cesLTVhundred{{{pctchange(ces['LTV100'], ces['base'])}}}\n"
        f"\\def\\cesLowerRisk{{{pctchange(ces['riskB'], ces['base'])}}}\n"
        f"\\def\\cesCapTen{{{pctchange(ces['cap10'], ces['cap6'])}}}\n")

    receipt = {"generated": datetime.datetime.now().isoformat(timespec="seconds"),
               "leader_file": str(a.leader), "leader_label": leader.get("label"), "leader_loss": leader["loss"],
               "price": leader["price"], "H0": leader["H0"], "rebate": leader.get("rebate"),
               "cd_table": str(a.cd_table), "ces_cells": str(a.ces_cells), "ces_rental": str(a.ces_rental),
               "cd_cells": str(a.cd_cells), "ces_completed_fertility": ces, "cd_completed_fertility": cdm}
    (OUT / "receipt.json").write_text(json.dumps(receipt, indent=1) + "\n")
    print("wrote", OUT, "leader", leader.get("label"), "loss", leader["loss"])

if __name__ == "__main__":
    main()
