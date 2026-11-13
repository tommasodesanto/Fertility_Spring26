"""Easy-to-read tables for the twins / same-sex IV with the household-position rule (2026-11-12).

Inputs (read-only):
  output/acs_fertility_iv/national_128g_results/rf_fs_iv_ar_table.csv   (all mothers, job 18274624)
  output/national_acs_comparison/headship_iv_20261112/collected/iv_restrict/rf_fs_iv_ar_table.csv
  output/national_acs_comparison/headship_iv_20261112/collected/iv_outcome/rf_fs_iv_ar_table.csv
Outputs: headship_iv_table.md / .csv in output/national_acs_comparison/headship_iv_20261112/.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output/national_acs_comparison/headship_iv_20261112"
SRC = {
    "All mothers (saved run 18274624)": ROOT / "output/acs_fertility_iv/national_128g_results/rf_fs_iv_ar_table.csv",
    "All mothers (rerun today)": OUT / "collected/iv_outcome/rf_fs_iv_ar_table.csv",
    "Householder/spouse/partner": OUT / "collected/iv_restrict/rf_fs_iv_ar_table.csv",
}
DESIGNS = {"Twin1_pooled0_5": "Twins at first birth, child ages 0-5 pooled",
           "SameSex2_pooled0_5": "First two children same sex, child ages 0-5 pooled",
           "Twin1_event3": "Twins, child age 3", "Twin1_event5": "Twins, child age 5",
           "SameSex2_event3": "Same-sex, child age 3", "SameSex2_event5": "Same-sex, child age 5"}


def fmt(c, s, sc=1.0, k=2):
    return f"{c * sc:.{k}f} ({s * sc:.{k}f})"


def main():
    rows, head = [], []
    for ver, p in SRC.items():
        if not p.exists():
            continue
        t = pd.read_csv(p)
        for d, dl in DESIGNS.items():
            z = t[t.design == d].set_index("outcome")
            if "ROOMS_out" not in z.index:
                continue
            r = z.loc["ROOMS_out"]
            row = {"design": dl, "version": ver,
                   "first stage, pp more likely extra child (SE)": fmt(r.fs_coef, r.fs_se, 100, 1),
                   "F": f"{r.first_stage_F:,.0f}",
                   "IV rooms (SE)": fmt(r.iv_coef, r.iv_se),
                   "IV bedrooms (SE)": fmt(z.loc["BEDROOMS_out"].iv_coef, z.loc["BEDROOMS_out"].iv_se),
                   "IV ownership, pp (SE)": fmt(z.loc["OWNERSHP_out"].iv_coef, z.loc["OWNERSHP_out"].iv_se, 100, 1),
                   "N mothers": f"{int(r.iv_nobs_fit):,}"}
            rows.append(row)
            if "HEAD_out" in z.index:
                h = z.loc["HEAD_out"]
                head.append({"design": dl,
                             "instrument effect on headship, pp (SE)": fmt(h.rf_coef, h.rf_se, 100, 2),
                             "IV: extra child effect on headship, pp (SE)": fmt(h.iv_coef, h.iv_se, 100, 1),
                             "N mothers": f"{int(h.rf_nobs_fit):,}"})
    t1 = pd.DataFrame(rows); t2 = pd.DataFrame(head)
    t1.to_csv(OUT / "headship_iv_table.csv", index=False)
    t2.to_csv(OUT / "headship_iv_headship_effect.csv", index=False)
    md = ["## IV effect of an extra child on housing (ACS 2005-19, mothers 21-35)", "",
          t1[t1.design.str.contains("pooled")].to_markdown(index=False), "",
          "## Instrument effect on the mother being householder/spouse/partner (all mothers)", "",
          t2[t2.design.str.contains("pooled")].to_markdown(index=False) if len(t2) else "(pending)", "",
          "## By child age (3 and 5)", "", t1[~t1.design.str.contains("pooled")].to_markdown(index=False)]
    (OUT / "headship_iv_table.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
