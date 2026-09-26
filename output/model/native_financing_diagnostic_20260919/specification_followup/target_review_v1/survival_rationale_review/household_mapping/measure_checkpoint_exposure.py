#!/usr/bin/env python3
"""Torch-only fixed-distribution exposure calculation; performs no model solve."""
from __future__ import annotations
import argparse, csv, gzip, json, pickle
from pathlib import Path
import numpy as np

EXPECTED_CHECKPOINT = "d5ef71bdaf9960273035c722a2428a55f14bab160e0596881c8a981e67b8ead1"
EXPECTED_SOURCE = "237904131d159f775c7ae89d1bbf1e8d1dd70c79ad012c80a36d948658d6f9c6"


def read_json(path):
    return json.loads(Path(path).read_text())


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--packet", required=True)
    ap.add_argument("--receipt", required=True)
    ap.add_argument("--lifecycle", required=True)
    ap.add_argument("--schedule", required=True)
    ap.add_argument("--output", required=True)
    a=ap.parse_args()
    receipt=read_json(a.receipt)
    if receipt.get("selected_reference") != EXPECTED_CHECKPOINT:
        raise RuntimeError("control packet does not point to author-specified selected checkpoint")
    if receipt.get("source_manifest_sha256") != EXPECTED_SOURCE:
        raise RuntimeError("control packet source manifest differs")
    with gzip.open(a.packet,"rb") as stream: packet=pickle.load(stream)
    P=packet["parameters"]
    g=np.asarray(packet["stationary_g_pre"],dtype=np.float64)
    J=int(P.J); ncs=int(P.n_child_states)
    if g.shape[-1] != ncs: raise RuntimeError(f"child state axis mismatch: shape={g.shape}, ncs={ncs}")
    axes=[i for i,size in enumerate(g.shape) if size==J]
    if len(axes)!=1: raise RuntimeError(f"age axis is ambiguous: shape={g.shape}, J={J}")
    age_axis=axes[0]
    by_age=np.moveaxis(g,age_axis,0)
    if by_age.shape[0] != J: raise RuntimeError("age move failed")
    age_mass=by_age.sum(axis=tuple(range(1,by_age.ndim)))
    # With K=1, cs=1 is the model's active dependent-child stage. The separate
    # parity axis nn carries the model child-count category used by child costs;
    # this exposure calculation sums over nn and therefore reports households,
    # not the corresponding m child-units. Biological child ages are not stored.
    dependent_mass=by_age[...,1].sum(axis=tuple(range(1,by_age[...,1].ndim)))
    schedule={int(r["age_start"]):r for r in csv.DictReader(open(a.schedule,newline=""))}
    if len(schedule)!=J-1: raise RuntimeError(f"schedule length {len(schedule)} does not match J-1={J-1}")
    if bool(P.use_age_survival) is not True: raise RuntimeError("serialized checkpoint survival switch is not enabled")
    current_s=np.asarray(P.survival_probs,dtype=float).reshape(-1)
    if current_s.size!=J-1: raise RuntimeError("serialized survival schedule length differs")
    entry_age=int(P.age_start); period=int(round(float(P.period_years)))
    candidate_exit=[]; candidate_first=[]; ages=[]
    for j in range(J-1):
        age=entry_age+j*period
        if age not in schedule: raise RuntimeError(f"candidate schedule misses age {age}")
        row=schedule[age]; ages.append(age)
        candidate_exit.append(float(row["household_exit_q4_mix"]))
        candidate_first.append(float(row["first_adult_death_q4_mix"]))
    candidate_exit=np.asarray(candidate_exit); candidate_first=np.asarray(candidate_first)
    current_exit=1.0-current_s
    # Reweight the saved age-cell distribution for the candidate survival path
    # while holding choices, child states and age-18 entry mass fixed.
    rel=np.ones(J,dtype=float)
    for j in range(J-1): rel[j+1]=rel[j]*(1-candidate_exit[j])/current_s[j]
    cand_age_mass=age_mass*rel
    cand_dep_mass=dependent_mass*rel
    worker_j=int(P.J_R)
    def stats(m):
        work=float(m[:worker_j].sum()); ret=float(m[worker_j:].sum())
        return {"working_household_mass_age18_62":work,
                "retiree_household_mass_age66_82":ret,
                "working_to_retiree_household_mass_ratio":work/ret,
                "total_adult_household_mass_age18_82":float(m.sum()),
                "mass_at_age66":float(m[worker_j]),
                "mass_at_terminal_cell_age82":float(m[-1])}
    rows=[]
    for j in range(J):
        age=entry_age+j*period
        if j<J-1:
            eh=float(candidate_exit[j]); fh=float(candidate_first[j]); ch=float(current_exit[j])
            first_event=float(cand_dep_mass[j]*fh)
            last_event=float(cand_dep_mass[j]*eh)
            old_event=float(dependent_mass[j]*ch)
        else:
            # Terminal age 82 has an age-support exit in the lifecycle model,
            # not another four-year mortality transition. Leave hazards undefined.
            eh=fh=ch=None
            first_event=last_event=old_event=0.0
        rows.append({"age_start":age,"saved_checkpoint_household_mass":float(age_mass[j]),
          "saved_checkpoint_dependent_state_mass_cs1":float(dependent_mass[j]),
          "candidate_fixed_distribution_household_mass":float(cand_age_mass[j]),
          "candidate_fixed_distribution_dependent_state_mass_cs1":float(cand_dep_mass[j]),
          "current_exit_probability_next_four_years":ch,
          "candidate_first_adult_death_probability_next_four_years":fh,
          "candidate_household_exit_probability_next_four_years":eh,
          "terminal_transition_defined":j<J-1,
          "terminal_scope_note":"" if j<J-1 else "No four-year mortality transition defined; terminal age-support exit is separate.",
          "candidate_dependent_households_exposed_to_nonterminal_first_adult_death_per_period":first_event,
          "candidate_dependent_households_exiting_nonterminal_per_period":last_event,
          "current_dependent_households_exiting_nonterminal_per_period":old_event})
    lifecycle={int(float(r["age_node"])):float(r["mass"]) for r in csv.DictReader(open(a.lifecycle,newline=""))}
    margdiff={str(entry_age+j*period):float(age_mass[j]-lifecycle.get(entry_age+j*period,float("nan"))) for j in range(J)}
    maxdiff=max(abs(v) for v in margdiff.values() if np.isfinite(v))
    if maxdiff>5e-8: raise RuntimeError(f"saved distribution age mass does not match compact control lifecycle (max gap {maxdiff})")
    out={"status":"completed_fixed_distribution_exposure_no_solve",
      "checkpoint_reference_sha256":EXPECTED_CHECKPOINT,"source_manifest_sha256":EXPECTED_SOURCE,
      "control_receipt_case":receipt.get("case"),"target_system_sha256":receipt.get("target_system_sha256"),
      "parameters":{"age_start":entry_age,"period_years":float(P.period_years),"J":J,"J_R":worker_j,
        "n_child_stages":int(P.n_child_stages),"n_child_states":ncs,
        "use_age_survival":bool(P.use_age_survival),"current_survival_probs":[float(x) for x in current_s]},
      "distribution":{"field":"stationary_g_pre","shape":[int(x) for x in g.shape],"age_axis":age_axis,"parity_axis_size":int(by_age.shape[-3]),"child_state_axis_size":ncs,
        "child_state_axis":"last axis; cs=1 is the active dependent-child stage; the separate parity axis nn carries model child-count category used in child costs; current age-only exposure marginalizes over nn and counts households, not model m child-units; biological child ages/minor status are absent",
        "age_marginal_max_abs_gap_vs_control_lifecycle_csv":maxdiff},
      "fixed_entry_demographic_aggregates":{"saved_checkpoint_current":stats(age_mass),"candidate_reweighted_fixed_distribution":stats(cand_age_mass),
        "candidate_minus_current_household_stock_percent":100*(float(cand_age_mass.sum())/float(age_mass.sum())-1),
        "interpretation":"age-cell reweighting under fixed age-18 entrant mass and fixed household choices; not a closed population, equilibrium, or behavioral response"},
      "dependent_state_exposure":{"saved_checkpoint_cs1_households_by_age":[float(x) for x in dependent_mass],
        "candidate_fixed_distribution_cs1_households_by_age":[float(x) for x in cand_dep_mass],
        "candidate_expected_cs1_households_exposed_to_nonterminal_first_adult_death_total_per_period":sum(r["candidate_dependent_households_exposed_to_nonterminal_first_adult_death_per_period"] for r in rows),
        "candidate_expected_cs1_nonterminal_household_exits_total_per_period":sum(r["candidate_dependent_households_exiting_nonterminal_per_period"] for r in rows),
        "current_expected_cs1_nonterminal_household_exits_total_per_period":sum(r["current_dependent_households_exiting_nonterminal_per_period"] for r in rows),
        "interpretation":"cs=1 marks model dependent stage; separate parity n can recover model-coded m from the joint n,cs state (subject to top-bin mapping); this saved exposure marginalizes over n and therefore reports household events, not m child-unit events. Biological child ages/minor status are not encoded."},
      "normalization":{"current_birth_based_adult_entry":float(receipt.get("birth_based_adult_entry")),
        "current_normalized_entry_rate":float(receipt.get("normalized_entrant_flow")),
        "birth_to_household_conversion":1/2.1,
        "candidate":"not recomputed or reclosed; adding mortality changes age-weighted fertility and deaths, so births/2.1 cannot silently be reused as a stationary replacement law"},
      "rows":rows,
      "terminal_scope":"reported mortality events are nonterminal transitions only (age starts 18 through 78); age-82 terminal population exit is not a mortality event in these totals. The age-82 terminal cell contains positive cs=1 mass and is not zero terminal exposure.",
      "not_done":["no household or equilibrium solve","no policy or fertility reoptimization","no spouse-to-child genealogical link or widow/head-age succession","no child-unit m exposure aggregation (parity axis was collapsed); no biological child-age/minor or orphan count"]}
    Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with Path(a.output).with_name("checkpoint_age_exposure.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({k:out[k] for k in ("status","distribution","fixed_entry_demographic_aggregates","dependent_state_exposure","normalization")},indent=2))

if __name__=="__main__": main()
