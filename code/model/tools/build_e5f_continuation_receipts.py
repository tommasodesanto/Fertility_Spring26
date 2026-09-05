#!/usr/bin/env python3
"""Read-only measurement and conditional sensitivity audit receipts."""
from __future__ import annotations
import argparse
import csv
import gc
import json
from pathlib import Path
import numpy as np
import run_e5f_independent_numerical_audit as audit


def read_csv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def age_receipts(panel, lifecycle, age_profile, out):
    annual = read_csv(panel / "ownership_annual_age_receipts.csv")
    original = json.loads((panel / "summary.json").read_text())["ownership_measurement"]
    expected = original["lifecycle_sha256"]
    if audit.digest(lifecycle) != expected:
        raise RuntimeError("Production lifecycle receipt differs")
    selected = [r for r in read_csv(lifecycle) if 26 <= float(r["age_node"]) <= 34]
    if [float(r["age_node"]) for r in selected] != [26.,30.,34.]:
        raise RuntimeError("Production ownership age-node selection differs")
    production = sum(float(r["mass"])*float(r["owner_rate"]) for r in selected)/sum(float(r["mass"]) for r in selected)
    weighted = {}
    for window in ("25_34", "26_37"):
        rows = [r for r in annual if r["age_window"] == window]
        weights = np.array([float(r["hhwt"]) for r in rows])
        weighted[window] = {k:float(weights @ np.array([float(r[k]) for r in rows]) / weights.sum())
            for k in ("acs_owner_rate", "model_cell_rate", "model_interpolated_rate")}
    specs = [
        ("stored comparison: model nodes 26,30,34; ACS ages 25-34", production, weighted["25_34"]["acs_owner_rate"]),
        ("whole-cell window aligned to ACS ages 26-37; original model weights", production, weighted["26_37"]["acs_owner_rate"]),
        ("ages 26-37; common ACS annual-age weights; constant model cell rates", weighted["26_37"]["model_cell_rate"], weighted["26_37"]["acs_owner_rate"]),
        ("ages 25-34; common ACS annual-age weights; constant model cell rates", weighted["25_34"]["model_cell_rate"], weighted["25_34"]["acs_owner_rate"]),
        ("ages 25-34; common ACS annual-age weights; linear between labelled nodes", weighted["25_34"]["model_interpolated_rate"], weighted["25_34"]["acs_owner_rate"])]
    rows = [dict(measurement=s, model=m, acs=t, gap_pp=100*(m-t), status="unchanged stored comparison" if i==0 else "diagnostic measurement only") for i,(s,m,t) in enumerate(specs)]
    audit.baseline.write_csv(out / "ownership_age_alignment.csv", rows)
    decomposition = dict(original_gap_pp=rows[0]["gap_pp"], full_cell_acs_weighted_gap_pp=rows[2]["gap_pp"],
        acs_window_change_contribution_pp=-100*(weighted["26_37"]["acs_owner_rate"]-weighted["25_34"]["acs_owner_rate"]),
        model_age_weight_change_contribution_pp=100*(weighted["26_37"]["model_cell_rate"]-production),
        interpretation="Ordered arithmetic decomposition: first align ACS window to the three full model cells, then apply common metropolitan ACS age weights. Not a newly approved calibration moment.")
    if abs(decomposition["original_gap_pp"]+decomposition["acs_window_change_contribution_pp"]+decomposition["model_age_weight_change_contribution_pp"]-decomposition["full_cell_acs_weighted_gap_pp"]) > 1e-12:
        raise RuntimeError("Age-gap decomposition does not add up")
    audit.save_json(out / "ownership_age_decomposition.json", decomposition)
    if audit.digest(age_profile) != original["source_sha256"]:
        raise RuntimeError("Authoritative annual-age ACS source differs")
    profile = [r for r in read_csv(age_profile) if r["sample"] == "household_heads_hhwt_due_housing"]
    by_node = {int(float(r["age_node"])):r for r in read_csv(lifecycle)}
    prime = [by_node[node] for node in range(30,55,4)]
    prime_production = sum(float(r["mass"])*float(r["owner_rate"]) for r in prime)/sum(float(r["mass"]) for r in prime)
    prime_rows = []
    for lo,hi in ((30,55),(30,57)):
        annual_rows = [r for r in profile if lo <= int(r["age"]) <= hi]
        den = sum(float(r["weight_sum"]) for r in annual_rows)
        target = sum(float(r["weight_sum"])*float(r["owner_rate"]) for r in annual_rows)/den
        model_rate = sum(float(r["weight_sum"])*float(by_node[18+4*((int(r["age"])-18)//4)]["owner_rate"]) for r in annual_rows)/den
        prime_rows.append(dict(acs_age_min=lo,acs_age_max=hi,acs_owner_rate=target,
            original_model_rate=prime_production,common_acs_weights_constant_cell_model_rate=model_rate,
            original_model_gap_pp=100*(prime_production-target),common_weights_gap_pp=100*(model_rate-target),
            status="diagnostic measurement only; active target and score unchanged"))
    audit.baseline.write_csv(out / "prime_ownership_age_alignment.csv",prime_rows)
    return dict(lifecycle_sha256=expected, acs_age_profile_sha256=original["source_sha256"],
        annual_age_receipt_sha256=audit.digest(panel / "ownership_annual_age_receipts.csv"), rows=rows,
        decomposition=decomposition, prime_age_alignment=prime_rows)


def exposure_receipts(saving_dir, out):
    rows = read_csv(saving_dir / "case_summary.csv")
    selected = {(r["method"],r["policy"]):r for r in rows if r["pricing"] == "separately_cleared"}
    summaries, detailed, hashes, replays = [], [], {}, {}
    for name in ("baseline", "supply-plus-20", "dependent-child-ltv95"):
        packets = {}
        for method in ("local", "global"):
            case = saving_dir / "cases" / f"{method}_{name}_separately_cleared" / "case_state.pkl.gz"
            expected = selected[method,name]["case_checkpoint_sha256"]
            if audit.digest(case) != expected:
                raise RuntimeError(f"Saving checkpoint differs: {case}")
            hashes[f"{method}_{name}"] = expected
            packets[method] = audit.load_checkpoint(case)
        local, global_ = packets["local"], packets["global"]
        P, e, e2, bg = local["parameters"], local["evaluation"], global_["evaluation"], local["b_grid"]
        if abs(e.policy.price[0]-e2.policy.price[0]) > 1e-14 or np.max(np.abs(e.g_pre-e2.g_pre)) > 2e-12:
            raise RuntimeError("Irregularity comparison requires a common price and pre-choice measure")
        settled = audit.model.readiness_settled_state(P)
        mass = e.g_pre[..., 0, settled]
        p_local, p_global = e.policy.fert_probs[..., 1], e2.policy.fert_probs[..., 1]
        delta = p_global-p_local
        pi = audit.model.get_fecundity_by_age(P).copy()
        fertile = (np.arange(P.J)+1 >= P.A_f_start) & (np.arange(P.J)+1 <= P.A_f_end)
        pi[~fertile] = 0
        birth_delta = mass * delta * pi[None,None,None,:,None]
        first_birth_change_from_stocks = float(e.g_post_fertility[...,0,settled].sum()-e2.g_post_fertility[...,0,settled].sum())
        replay_gap = float(birth_delta.sum())-first_birth_change_from_stocks
        if abs(replay_gap) > 2e-12:
            raise RuntimeError("First-birth probability attribution disagrees with post-birth childless stocks")
        occupied_fertile = (mass > 1e-12) & fertile[None,None,None,:,None]
        replays[name] = dict(probability_first_birth_change=float(birth_delta.sum()),
            post_fertility_childless_stock_first_birth_change=first_birth_change_from_stocks,
            replay_gap=replay_gap, maximum_occupied_fertile_probability_change=float(np.max(np.abs(delta[occupied_fertile]))))
        valid = (e.policy.V[...,0,settled] > -1e9) & (e2.policy.V[...,0,settled] > -1e9)
        large = valid & (np.abs(delta) > .1) & fertile[None,None,None,:,None]
        spike = valid & (p_local > .9) & (p_global < .5) & (bg[:,None,None,None,None] < 0)
        age42 = audit.model.age_to_index(P,42)
        selected_spike = np.zeros_like(spike)
        selected_spike[:,0,0,age42,:] = spike[:,0,0,age42,:]
        for scope, mask in (("all fertile childless states: absolute probability change above 0.1",large),
                            ("age-42 childless renter, negative wealth: local above 0.9 and global below 0.5",selected_spike)):
            summaries.append(dict(policy=name, scope=scope, grid_states=int(mask.sum()), occupied_grid_states=int(np.count_nonzero(mask & (mass > 1e-12))),
                exposed_household_mass=float(mass[mask].sum()), share_all_households=float(mass[mask].sum()/e.g_pre.sum()),
                signed_first_birth_change=float(birth_delta[mask].sum()),
                maximum_probability_change=float(np.max(np.abs(delta[mask]))) if np.any(mask) else 0.,
                occupancy_threshold=1e-12, exposure_mass_uses_all_positive_mass=True))
        for idx in np.argwhere(selected_spike):
            b,t,i,j,z = map(int,idx)
            detailed.append(dict(policy=name, wealth=float(bg[b]), age=float(P.age_start+j*P.da), income_index=z,
                household_mass=float(mass[tuple(idx)]), local_attempt_probability=float(p_local[tuple(idx)]),
                global_attempt_probability=float(p_global[tuple(idx)]), fecundity=float(pi[j]),
                first_birth_change=float(birth_delta[tuple(idx)])))
        summaries.append(dict(policy=name,scope="all first births: exact frozen-state probability difference",grid_states=None,occupied_grid_states=None,
            exposed_household_mass=float(mass.sum()),share_all_households=float(mass.sum()/e.g_pre.sum()),
            signed_first_birth_change=float(birth_delta.sum()),maximum_probability_change=float(np.max(np.abs(delta))),
            occupancy_threshold=None,exposure_mass_uses_all_positive_mass=True))
        del packets, local, global_, e, e2, mass, p_local, p_global, delta, birth_delta, valid, large, spike, selected_spike
        gc.collect()
    audit.baseline.write_csv(out / "first_birth_irregularity_exposure.csv", summaries)
    audit.baseline.write_csv(out / "age42_first_birth_spike_states.csv", detailed)
    return dict(saving_summary_sha256=audit.digest(saving_dir / "case_summary.csv"), case_checkpoint_hashes=hashes,
        first_birth_stock_replay=replays, thresholds="Screening definitions stated in each row; not inference or welfare bounds", rows=summaries)


def grid_receipts(stages, saving_dir, out):
    rows = [dict(r,nodes=120,subdivisions=1) for r in read_csv(saving_dir / "case_summary.csv") if r["method"]=="global" and r["pricing"]=="separately_cleared"]
    receipts = {}
    for stage in stages:
        summary_path = stage / "summary.json"
        summary = json.loads(summary_path.read_text())
        if summary.get("status") != "complete" or not summary.get("all_case_gates_passed"):
            raise RuntimeError(f"Grid stage incomplete: {stage}")
        receipts[str(stage)] = audit.digest(summary_path)
        for row in summary["rows"]:
            if not row["market_cleared"]:
                continue
            folder = stage / "cases" / row["case"]
            if audit.digest(folder / "case_state.pkl.gz") != row["case_checkpoint_sha256"]:
                raise RuntimeError("Grid checkpoint differs")
            if len(list((folder / "standard_diagnostics").glob("*.png"))) != 17:
                raise RuntimeError("Standard graph packet incomplete")
            rows.append(row)
    effects = []
    levels = sorted(set(int(r["nodes"]) for r in rows))
    for nodes in levels:
        subset = {r["policy"]:r for r in rows if int(r["nodes"]) == nodes}
        if set(subset) != {"baseline", "supply-plus-20", "dependent-child-ltv95"}:
            raise RuntimeError(f"Incomplete three-policy comparison at {nodes} nodes")
        b = subset["baseline"]
        for name in ("supply-plus-20", "dependent-child-ltv95"):
            r = subset[name]
            effects.append(dict(nodes=nodes,policy=name,
                births_pct=100*(float(r["births_per_household"])/float(b["births_per_household"])-1),
                rooms_pct=100*(float(r["rooms_per_household"])/float(b["rooms_per_household"])-1),
                ownership_pp=100*(float(r["ownership"])-float(b["ownership"])),
                price_pct=100*(float(r["asset_price"])/float(b["asset_price"])-1)))
    fields=("nodes","policy","births_per_household","ownership","rooms_per_household","asset_price","market_residual","budget_excess_share","case_checkpoint_sha256")
    audit.baseline.write_csv(out / "grid_levels.csv", [{k:r[k] for k in fields} for r in rows])
    audit.baseline.write_csv(out / "grid_policy_effects.csv",effects)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(9,3.5))
    for ax,name,title in zip(axes,("supply-plus-20","dependent-child-ltv95"),("Housing supply +20%","Dependent-child 95% LTV")):
        selected=[r for r in effects if r["policy"]==name]
        values=[r["births_pct"] for r in selected]
        bars=ax.bar([str(r["nodes"]) for r in selected],values,color="#183f59",width=.6)
        ax.bar_label(bars,labels=[f"{v:.6f}%" for v in values],padding=4,fontsize=9)
        ax.set_ylim(0,max(values)*1.22)
        ax.set_xlabel("Nested liquid-wealth grid nodes")
        ax.set_ylabel("Impact births per household (%)");ax.set_title(title);ax.grid(alpha=.2)
    fig.suptitle("Supplemental: fixed history, global saving, separately cleared prices",fontsize=10)
    fig.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(out / f"supplemental_nested_grid.{ext}",dpi=180)
    plt.close(fig)
    return dict(stage_summary_hashes=receipts, effects=effects,
        interpretation="Conditional impact sensitivity only. Common original pre-choice measure and unchanged parameters. Neither history nor calibration is recomputed.")


def rental_receipts(stages, out, grid_stages):
    """Compare the economic rental-cap variant with the same-grid reference."""
    rows = [dict(r,rental_cap=6.0) for r in read_csv(out / "grid_levels.csv") if int(r["nodes"]) == 239]
    if len(rows) != 3:
        raise RuntimeError("Complete 239-node reference required for rental-cap comparison")
    references = []
    for stage in grid_stages:
        for row in json.loads((stage / "summary.json").read_text())["rows"]:
            if int(row["nodes"]) == 239 and row["policy"] == "baseline":
                references.append((stage / "cases" / row["case"] / "case_state.pkl.gz",row["case_checkpoint_sha256"]))
    if len(references) != 1 or audit.digest(references[0][0]) != references[0][1]:
        raise RuntimeError("Unique hash-verified fine-grid baseline required")
    reference = audit.load_checkpoint(references[0][0])
    reference_g = reference["state"].g_pre.copy()
    reference_grid = reference["b_grid"].copy()
    del reference
    gc.collect()
    hashes = {}
    inherited_checks = []
    for stage in stages:
        path = stage / "summary.json";summary = json.loads(path.read_text())
        if summary.get("status") != "complete" or not summary.get("all_case_gates_passed"):
            raise RuntimeError(f"Rental-cap stage incomplete: {stage}")
        if summary.get("economic_parameter_changes") != {"hR_max":{"reference":6.0,"diagnostic":8.0}}:
            raise RuntimeError("Unexpected cap-variant economic changes")
        hashes[str(stage)] = audit.digest(path)
        for row in summary["rows"]:
            folder = stage / "cases" / Path(row["diagnostic_directory"]).name
            if audit.digest(folder / "case_state.pkl.gz") != row["case_checkpoint_sha256"]:
                raise RuntimeError("Rental-cap checkpoint differs")
            if len(list((folder / "standard_diagnostics").glob("*.png"))) != 17:
                raise RuntimeError("Rental-cap standard graph packet incomplete")
            if not row["all_gates_passed"] or row["nodes"] != 239 or row["rental_cap"] != 8.0:
                raise RuntimeError("Rental-cap case contract failed")
            packet = audit.load_checkpoint(folder / "case_state.pkl.gz")
            if not np.array_equal(packet["b_grid"],reference_grid) or not np.array_equal(packet["state"].g_pre,reference_g):
                raise RuntimeError("Rental-cap inherited input distribution or wealth grid changed")
            evaluation=packet["evaluation"]
            processed,projected=audit.calendar.gate_pre_fertility_distribution(reference_g,
                evaluation.policy,packet["parameters"],packet["b_grid"],packet["shared"])
            if not np.array_equal(processed,evaluation.g_pre) or projected != evaluation.feasibility_projection_mass:
                raise RuntimeError("Saved pre-choice preparation does not exactly replay")
            inherited_checks.append(dict(policy=row["policy"],case_checkpoint_sha256=row["case_checkpoint_sha256"],
                exact_grid_identity=True,exact_inherited_input_measure_identity=True,
                exact_processed_measure_identity=np.array_equal(evaluation.g_pre,reference_g),
                exact_pre_choice_preparation_replay=True,feasibility_projection_mass=float(projected),
                processed_measure_l1_difference=float(np.sum(np.abs(evaluation.g_pre-reference_g))),
                maximum_initial_input_mass_difference=0.0,reference_checkpoint_sha256=references[0][1]))
            del packet,evaluation,processed
            gc.collect()
            rows.append(row)
    effects = [];baselines = {}
    for cap in (6.,8.):
        subset = {r["policy"]:r for r in rows if float(r["rental_cap"]) == cap}
        if set(subset) != {"baseline","supply-plus-20","dependent-child-ltv95"}:
            raise RuntimeError("Incomplete three-policy rental-cap comparison")
        b = subset["baseline"];baselines[cap] = b
        for name in ("supply-plus-20","dependent-child-ltv95"):
            r = subset[name]
            effects.append(dict(rental_cap=cap,nodes=239,policy=name,
                births_pct=100*(float(r["births_per_household"])/float(b["births_per_household"])-1),
                rooms_pct=100*(float(r["rooms_per_household"])/float(b["rooms_per_household"])-1),
                ownership_pp=100*(float(r["ownership"])-float(b["ownership"])),
                price_pct=100*(float(r["asset_price"])/float(b["asset_price"])-1)))
    fields=("rental_cap","nodes","policy","births_per_household","ownership","rooms_per_household","asset_price","market_residual","budget_excess_share","case_checkpoint_sha256")
    audit.baseline.write_csv(out / "rental_cap_levels.csv",[{k:r[k] for k in fields} for r in rows])
    audit.baseline.write_csv(out / "rental_cap_policy_effects.csv",effects)
    b6,b8=baselines[6.],baselines[8.]
    baseline_change=dict(births_pct=100*(float(b8["births_per_household"])/float(b6["births_per_household"])-1),
        rooms_pct=100*(float(b8["rooms_per_household"])/float(b6["rooms_per_household"])-1),
        ownership_pp=100*(float(b8["ownership"])-float(b6["ownership"])),
        price_pct=100*(float(b8["asset_price"])/float(b6["asset_price"])-1))
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,4.8))
    for column,(name,title) in enumerate(zip(("supply-plus-20","dependent-child-ltv95"),("Housing supply +20%","Dependent-child 95% LTV"))):
        selected=[r for r in effects if r["policy"]==name]
        for index,(field,label,suffix) in enumerate((("births_pct","Births per household (%)","%"),("ownership_pp","Ownership (percentage points)"," pp"))):
            ax=axes[index,column];values=[r[field] for r in selected]
            bars=ax.bar(["Six rooms","Eight rooms"],values,color=["#183f59","#b76532"])
            ax.bar_label(bars,labels=[f"{v:.4f}{suffix}" for v in values],padding=4,fontsize=9)
            ax.margins(y=.25)
            ax.axhline(0,color="black",linewidth=.6);ax.set_ylabel(label)
            if index == 0: ax.set_title(title)
            else: ax.set_xlabel("Rental upper bound")
            ax.grid(axis="y",alpha=.2)
    fig.suptitle("Supplemental economic sensitivity: fixed parameters and history, 239 nodes",fontsize=10)
    fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(out / f"supplemental_rental_cap.{ext}",dpi=180)
    plt.close(fig)
    audit.save_json(out / "rental_inherited_measure_checks.json",inherited_checks)
    return dict(stage_summary_hashes=hashes,effects=effects,baseline_change_from_six_to_eight=baseline_change,
        inherited_measure_checks=inherited_checks,
        interpretation="An uncalibrated economic opportunity change at fixed parameters and inherited history. Each policy is compared with its own cap-specific baseline. Not a numerical repair, first-birth room-target remeasurement, long-run or welfare estimate.")


def housing_stock_receipts(path, out):
    """Independently aggregate the existing descriptive ACS/MMS stock cells."""
    expected = "2c9ae796b6bc8aeb6f37df6b640c5a52e8a811b40fa1d6a9eac4144fbffb052b"
    if audit.digest(path) != expected:
        raise RuntimeError("Reviewed housing-stock cells differ")
    cells = read_csv(path)
    keys = [tuple(r[k] for k in ("met2013","mms_location","size_bin","tenure")) for r in cells]
    if len(set(keys)) != len(keys) or sum(int(r["n_hh"]) for r in cells) != 4103889:
        raise RuntimeError("Housing-stock cell keys or reported sample count differ")
    rows = []
    for tenure in ("renter","owner"):
        selected = [r for r in cells if r["tenure"] == tenure]
        denominator = sum(float(r["hh_weight"]) for r in selected)
        for size in ("S_1_4","M_5_6","L_7plus"):
            group = [r for r in selected if r["size_bin"] == size]
            weight = sum(float(r["hh_weight"]) for r in group)
            if not group or min(float(r["hh_weight"]) for r in group) <= 0:
                raise RuntimeError("Invalid housing-stock weights")
            rows.append(dict(tenure=tenure,size_bin=size,weighted_share=weight/denominator,
                n_hh=sum(int(r["n_hh"]) for r in group),household_weight=weight,
                mean_rooms=sum(float(r["mean_rooms"])*float(r["hh_weight"]) for r in group)/weight))
        if abs(sum(r["weighted_share"] for r in rows if r["tenure"] == tenure)-1) > 1e-12:
            raise RuntimeError("Housing-stock shares do not add up")
    audit.baseline.write_csv(out / "descriptive_acs_housing_stock.csv",rows)
    return dict(source=str(path),source_sha256=expected,rows=rows,
        authoritative_builder="code/data/mms_center_periphery/analyze_family_size_supply_menu.R",
        estimator="Household-weighted stock shares, aggregating mutually exclusive metro-location-size-tenure cells",
        sample="Existing 2012-2023 packet: matched MMS metropolitan household heads age 18+, GQ codes 1 or 2, positive HHWT and rooms; middle PUMAs assigned to center",
        uncertainty="No empirical standard error computed from saved aggregate cells",
        fixed_effects="Not applicable to these descriptive shares",clustering="Not applicable; no inference performed",
        scope="Saved-cell arithmetic and builder inspection, not a raw-IPUMS rebuild or newly approved calibration target. This sample does not impose the active ownership target's standard-structure filter.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--saving-dir",type=Path,required=True)
    parser.add_argument("--panel-dir",type=Path,required=True)
    parser.add_argument("--lifecycle",type=Path,required=True)
    parser.add_argument("--acs-age-profile",type=Path,required=True)
    parser.add_argument("--grid-stage",type=Path,action="append",default=[])
    parser.add_argument("--rental-stage",type=Path,action="append",default=[])
    parser.add_argument("--housing-stock-cells",type=Path)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--reuse-exposure",action="store_true")
    args=parser.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    _,configured=audit.transition.configure_sequential_model()
    live=audit.baseline.calibration.code_fingerprint_contract(configured)["bundle_sha256"]
    if live != audit.BUNDLE or audit.digest(audit.__file__) != "ba0f94f52705f43fd0bf7a18ea8c2053f9897674abeb20b030eaf37c57de8623":
        raise RuntimeError("Read-only receipts require the frozen production source and audit helper")
    ages=age_receipts(args.panel_dir,args.lifecycle,args.acs_age_profile,out)
    if args.reuse_exposure:
        exposure=json.loads((out / "exposure_receipt.json").read_text())
        if exposure["saving_summary_sha256"] != audit.digest(args.saving_dir / "case_summary.csv"):
            raise RuntimeError("Reused exposure receipt references different cases")
    else:
        exposure=exposure_receipts(args.saving_dir,out)
        audit.save_json(out / "exposure_receipt.json",exposure)
    grid=grid_receipts(args.grid_stage,args.saving_dir,out) if args.grid_stage else None
    if args.rental_stage and grid is None:
        raise RuntimeError("Rental-cap comparison requires verified grid reference stages")
    rental=rental_receipts(args.rental_stage,out,args.grid_stage) if args.rental_stage else None
    stock=housing_stock_receipts(args.housing_stock_cells,out) if args.housing_stock_cells else None
    audit.save_json(out / "summary.json",dict(status="complete",driver_sha256=audit.digest(__file__),source_bundle_sha256=live,
        age_alignment=ages,irregularity_exposure=exposure,grid_comparison=grid,rental_cap_comparison=rental,
        descriptive_housing_stock=stock))


if __name__=="__main__":
    main()
