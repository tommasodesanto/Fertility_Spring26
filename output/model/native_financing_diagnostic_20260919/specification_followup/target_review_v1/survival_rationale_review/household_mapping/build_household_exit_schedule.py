#!/usr/bin/env python3
"""Small standard-library calculation: 2007 qx + SCF 2007 household mix."""
from __future__ import annotations
import csv, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
BASE = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/bequest_flow_2007"
QX_PATH = BASE / "mortality/nchs_2007_qx_by_age_sex.csv"
SCF_PATH = BASE / "calculation/scf2007_keys_ages.csv"
OUT = Path(__file__).resolve().parent
AGES = list(range(18, 82, 4))  # transition start ages, 18 through 78; age-82 is terminal cell
J_R = 12
J = 17
PERIODS = 4
CURRENT_POSTRET_SURVIVAL = {66: 0.9391263063710125, 70: 0.9184976343249724,
                            74: 0.8849521927812863, 78: 0.8300468061015381}


def q4(qx: dict[int, tuple[float, float]], age: int, sex: int, topcode: int) -> float:
    # SCF public spouse ages are top-coded at 95. Main maps 95+ to qx95;
    # sensitivity maps to qx99 or terminal 100+.
    start = min(max(age, 0), topcode)
    col = 0 if sex == 1 else 1
    survival = 1.0
    for a in range(start, start + PERIODS):
        qa = qx.get(a, qx[100])[col]
        survival *= 1.0 - qa
    return 1.0 - survival


def load_qx() -> dict[int, tuple[float, float]]:
    with QX_PATH.open(newline="") as f:
        return {int(r["age_start"]): (float(r["male_qx"]), float(r["female_qx"])) for r in csv.DictReader(f)}


def load_peus() -> list[dict[str, float | int]]:
    # Official SCF simple-estimate convention: X42001/5 across five implicates.
    # Head age/sex and living arrangement are invariant; spouse age may be imputed.
    grouped: dict[str, list[dict[str, str]]] = {}
    with SCF_PATH.open(newline="") as f:
        for r in csv.DictReader(f):
            grouped.setdefault(r["YY1"], []).append(r)
    demographic = ("X8021", "X8023", "X103", "X14")
    out = []
    mismatch = []
    for key, rows in grouped.items():
        first = rows[0]
        if any(any(r[k] != first[k] for k in demographic) for r in rows[1:]):
            mismatch.append(key)
        # Retain all five SCF implicates with official simple-statistic weight
        # X42001/5; X19 (partner age) is sometimes imputed and can vary.
        for imp, r in enumerate(rows):
            out.append({k: int(float(r[k])) for k in (*demographic, "X19")} |
                       {"peu": key, "weight": float(r["X42001"])/5.0,
                        "n_implicates": len(rows), "implicate": imp+1})
    if mismatch:
        raise ValueError(f"Head demographics/living arrangement vary across implicates for {len(mismatch)} PEUs")
    return out


def weighted_mean(rows, key):
    den = sum(float(r["weight"]) for r in rows)
    return sum(float(r["weight"]) * float(r[key]) for r in rows) / den if den else math.nan


def run(topcode: int, suffix: str):
    qx = load_qx()
    peus = load_peus()
    schedule = []
    for age0 in AGES:
        rows = [r for r in peus if age0 <= int(r["X14"]) < age0 + 4 and 1 <= int(r["X8021"]) <= 2]
        wgt = sum(float(r["weight"]) for r in rows)
        q_head = []
        q_first = []
        q_exit = []
        q_single = []
        q_couple = []
        cp = []
        for r in rows:
            headage, headsex = int(r["X14"]), int(r["X8021"])
            qh = q4(qx, headage, headsex, topcode)
            partnered = int(r["X8023"]) in (1, 2) and int(r["X103"]) in (1, 2) and int(r["X19"]) >= 18
            wt = float(r["weight"])
            q_head.append((wt, qh))
            if partnered:
                qp = q4(qx, int(r["X19"]), int(r["X103"]), topcode)
                first = 1.0 - (1.0-qh)*(1.0-qp)
                last = qh*qp  # both adults die within the same four-year window, independence
                q_first.append((wt, first)); q_exit.append((wt, last)); q_couple.append((wt, last)); cp.append((wt, 1.0))
            else:
                q_first.append((wt, qh)); q_exit.append((wt, qh)); q_single.append((wt, qh)); cp.append((wt, 0.0))
        def wm(items):
            den=sum(w for w,_ in items)
            return sum(w*v for w,v in items)/den if den else math.nan
        schedule.append({
            "age_start": age0, "age_interval": f"{age0}-{age0+4}",
            "scf_unweighted_peus": len({r["peu"] for r in rows}), "scf_weighted_peus": wgt,
            "couple_share": wm(cp), "single_head_q4": wm(q_single),
            "couple_last_survivor_exit_q4": wm(q_couple),
            "reference_adult_q4_all": wm(q_head),
            "first_adult_death_q4_mix": wm(q_first),
            "household_exit_q4_mix": wm(q_exit),
            "head_age_bins": [age0, age0+1, age0+2, age0+3],
        })
    if suffix == "main":
        outpath = OUT / "candidate_schedule.csv"
        with outpath.open("w", newline="") as f:
            w=csv.DictWriter(f,fieldnames=[k for k in schedule[0] if k != "head_age_bins"])
            w.writeheader(); w.writerows({k:v for k,v in row.items() if k != "head_age_bins"} for row in schedule)
        global LAST_SCHEDULE
        LAST_SCHEDULE = schedule
    return schedule


def profile(hazards):
    mass=[1.0]
    for h in hazards:
        mass.append(mass[-1]*(1.0-h))
    worker=sum(mass[:J_R]); retiree=sum(mass[J_R:J])
    return {"age_cell_mass_from_one_age18_entrant":mass,
            "worker_mass_18_66":worker,"retiree_mass_66_82_terminal_cell":retiree,
            "worker_to_retiree_mass_ratio":worker/retiree,
            "survival_to_age66":mass[J_R],"survival_to_age82_cell":mass[J-1],
            "total_age_cell_person_year_equivalent_periods":sum(mass)}


def main():
    global LAST_SCHEDULE
    qx=load_qx(); peus=load_peus()
    n=len({r["peu"] for r in peus})
    main_schedule=run(95,"main"); sens99=run(99,"sens99"); sens100=run(100,"sens100")
    # All hazards are conditional four-year probabilities for each age cell.
    candidate=[float(r["household_exit_q4_mix"]) for r in main_schedule]
    first=[float(r["first_adult_death_q4_mix"]) for r in main_schedule]
    head=[float(r["reference_adult_q4_all"]) for r in main_schedule]
    current=[0.0]*len(AGES)
    for k,a in enumerate(AGES):
        current[k]=1.0-CURRENT_POSTRET_SURVIVAL[a] if a in CURRENT_POSTRET_SURVIVAL else 0.0
    sens_candidate=[float(r["household_exit_q4_mix"]) for r in sens99]
    # An independent q4 check across ages 18-66 using 2007 official source.
    def survival_sex(sex):
        s=1.0
        for a in range(18,66): s*=1.0-qx[a][0 if sex==1 else 1]
        return s
    m_surv,f_surv=survival_sex(1),survival_sex(2)
    result={
      "scope":"cross-sectional 2007 age-specific household-exit proxy, not dynamic couple/head succession or a GE/model solve",
      "sources":{"mortality":"NCHS United States Life Tables, 2007; exact-age one-year qx by sex","household_mix":"Federal Reserve 2007 SCF public PEU keys/ages; X8023=married or living with partner; weight X42001/5 across five implicates"},
      "sample":{"unique_peus":n,"official_weight_rule":"X42001/5 on each of five implicate rows; head age/sex and living arrangement checked invariant; spouse age X19 is retained by implicate because it can vary"},
      "mapping":{"single_household_exit":"head dies within next four years","partnered_household_exit":"both adults die within next four years; assumes independent spouse death times and household persists after first death","first_adult_death":"1-(1-q_head_4)(1-q_partner_4) for partnered units; q_head_4 for singles","candidate_age_hazard":"SCF-weighted single/couple mix at reference-person age; partner age/sex used in paired hazard","age_groups":"each start age x pools SCF reference-person ages x through x+3","topcode":"main maps SCF public spouse age 95+ to qx95; sensitivities map to qx99 and qx100plus"},
      "life_table_scale":{"male_survival_18_to_66":m_surv,"male_death_risk_18_to_66":1-m_surv,"female_survival_18_to_66":f_surv,"female_death_risk_18_to_66":1-f_surv,"equal_sex_count_survival_18_to_66":(m_surv+f_surv)/2,"equal_sex_count_death_risk_18_to_66":1-(m_surv+f_surv)/2},
      "demographic_profiles":{"selected_checkpoint_current_effective_2023_postretirement_only":profile(current),"2007_candidate_household_exit_mix":profile(candidate),"2007_head_death_upper_bound":profile(head),"2007_first_adult_death_as_exit_upper_bound":profile(first),"2007_candidate_topcode_qx99_sensitivity":profile(sens_candidate)},
      "schedule_sensitivities":{"topcode_qx95":main_schedule,"topcode_qx99":sens99,"topcode_qx100plus":sens100},
      "limitations":["Cross-sectional SCF marital/living-arrangement shares are held fixed at each age; this is an age-specific one-state exit schedule, not a dynamically updated widow/single process.","First-adult death exposure is distinct from household exit; the model has no spouse status, and a surviving spouse may change head age and earnings.","No household children-at-home exposure is calculated here; checkpoint distribution exposure is a separate Torch-only step.","No births/2.1 demographic closure is inferred. Candidate hazards change age masses, and fertility/entry must be recomputed before claiming a closed population."]
    }
    (OUT/"demographic_profiles.json").write_text(json.dumps(result,indent=2)+"\n")
    with (OUT/"household_composition_by_age.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=[k for k in main_schedule[0] if k!="head_age_bins"]); w.writeheader();w.writerows({k:v for k,v in r.items() if k!="head_age_bins"} for r in main_schedule)

if __name__=="__main__": main()
