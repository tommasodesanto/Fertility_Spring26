#!/usr/bin/env python3
"""Rebuild complete panel tables and compare 12 versus 13 moment sensitivity.

The strict ridge planner remains untouched. This records its rejection and
separately reconstructs exact flow identities and the original solver gates.
No parameter proposals, model solves, target changes, or promotion occur.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
import analyze_e5f_transition_calibration_panel as panel
import build_e5f_transition_ridge_refinement as ridge

ROOT=Path(__file__).resolve().parents[3]
PANEL=ROOT/"output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_youngown_overid_coord_20260904a"
OUT=ROOT/"output/model/e5f_overnight_independent_verification_20260905a/panel_receipts"


def ownership_measurement(out,aggregate_room_target):
    lifecycle=out.parent/"numerical_full/lifecycle_2023.csv"
    if not lifecycle.exists():
        return dict(status="dated lifecycle not yet collected")
    ages=panel.read_csv(ROOT/"code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv")
    ages=[r for r in ages if r["source"]=="ACS" and r["sample"]=="household_heads_hhwt_due_housing"]
    model_rows=panel.read_csv(lifecycle)
    model={int(float(r["age_node"])):float(r["owner_rate"]) for r in model_rows}
    rows=[];detail=[]
    for lo,hi in ((25,34),(26,37),(65,75),(80,84)):
        sample=[r for r in ages if lo<=int(r["age"])<=hi]
        denominator=sum(float(r["weight_sum"]) for r in sample)
        target=sum(float(r["weight_sum"])*float(r["owner_rate"]) for r in sample)/denominator
        start=sum(float(r["weight_sum"])*model[min(82,18+4*((int(r["age"])-18)//4))] for r in sample)/denominator
        linear=sum(float(r["weight_sum"])*float(np.interp(int(r["age"]),list(model),list(model.values()))) for r in sample)/denominator
        rows.append(dict(age_min=lo,age_max=hi,acs_hhwt=denominator,target=target,
            model_cell_start=start,cell_start_gap_pp=100*(start-target),
            model_linear_between_nodes=linear,linear_gap_pp=100*(linear-target)))
        for r in sample:
            age=int(r["age"]);node=min(82,18+4*((age-18)//4))
            detail.append(dict(age_window=f"{lo}_{hi}",age=age,hhwt=float(r["weight_sum"]),
                acs_owner_rate=float(r["owner_rate"]),model_node=node,model_cell_rate=model[node],
                model_interpolated_rate=float(np.interp(age,list(model),list(model.values())))))
    panel.write_csv(out/"ownership_aligned_age_windows.csv",rows)
    panel.write_csv(out/"ownership_annual_age_receipts.csv",detail)
    import matplotlib.pyplot as plt
    figure_dir=out.parent/"supplemental_figures"
    figure_dir.mkdir(exist_ok=True)
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False}):
        fig,axes=plt.subplots(1,2,figsize=(7.2,2.7),layout="constrained")
        axes[0].plot(list(model),100*np.asarray(list(model.values())),"-o",color="#173E55",lw=1.6,markersize=2.5,label="Model, 2023 nodes")
        axes[0].plot([int(r["age"]) for r in ages],[100*float(r["owner_rate"]) for r in ages],color="#B35A32",lw=1.5,label="ACS, pooled 2012-23")
        axes[0].set(xlabel="Age",ylabel="Ownership (%)",ylim=(0,105),title="Tenure over the lifecycle")
        axes[1].plot([float(r["age_node"]) for r in model_rows],[float(r["mean_rooms"]) for r in model_rows],"-o",color="#173E55",lw=1.6,markersize=2.5,label="Model, 2023 nodes")
        axes[1].axhline(aggregate_room_target,color="#B35A32",ls="--",label="Aggregate rooms target")
        axes[1].set(xlabel="Age",ylabel="Mean occupied rooms",title="Housing use over the lifecycle")
        for ax in axes:
            ax.legend(frameon=False,fontsize=7.5,loc="lower right")
            ax.grid(axis="y",alpha=.17)
        fig.savefig(figure_dir/"lifecycle_validation.pdf",bbox_inches="tight")
        fig.savefig(figure_dir/"lifecycle_validation.png",dpi=200,bbox_inches="tight")
        plt.close(fig)
    return dict(status="complete diagnostic remeasurement",rows=rows,
        assumptions="Same ACS annual-age HHWT on both sides. Cell-start mapping holds each model rate over [node,node+4); linear variant interpolates between age nodes and holds the last node at ages 82-84. Neither changes the production target contract.",
        source_sha256=panel.sha256(ROOT/"code/data/mms_center_periphery/output_ownership_audit/acs_ownership_age_profiles.csv"),
        lifecycle_sha256=panel.sha256(lifecycle))


def smallstep_receipts(root,out,original,moments,targets,weights,Jlarge):
    """Compare three columns at two step sizes; never impute missing columns."""
    ids=(1,2,3,8,9,14,15)
    actual=sorted(p.name for p in root.glob("task_*") if p.is_dir())
    if actual!=[f"task_{i:03d}" for i in ids]:
        raise RuntimeError("Small-step readout requires precisely the seven declared tasks")
    out.mkdir(parents=True,exist_ok=True)
    expected_contract=ridge.scientific_contract(original[0].summary)
    tasks={};fits=[];parameters=[];case_rows=[];hashes={}
    mask=np.asarray([m!="own_rate_2534" for m in moments])
    for i in ids:
        directory=root/f"task_{i:03d}"
        summary=panel.read_json(directory/"summary.json")
        design=summary["panel_design"];candidate=summary["best_candidate"]
        expected_design="anchor" if i==1 else f"coordinate_{(i-2)//2}_{'minus' if i%2==0 else 'plus'}"
        if summary.get("status")!="complete_transition_calibration_panel_task" or design.get("design")!=expected_design:
            raise RuntimeError(f"Small-step task {i} is incomplete or has the wrong design")
        if design["panel_seed"]!=2026090505 or not math.isclose(design["local_radius"],.005,abs_tol=1e-14):
            raise RuntimeError("Small-step radius/seed contract failed")
        for name in ("center_sha256","domain","panel_size","panel_design"):
            if design[name]!=original[0].panel[name]: raise RuntimeError(f"Small-step {name} differs")
        if json.dumps(ridge.scientific_contract(summary),sort_keys=True)!=json.dumps(expected_contract,sort_keys=True):
            raise RuntimeError("Small-step scientific contract differs from original panel")
        ridge.validate_candidate_eligibility(summary)
        o=summary["renewal_accounting_old_state"];g=summary["stationary_operator_gates"]
        E=o["old_entry_flow_E"];B=o["old_queue_mature_flow_B"]
        if max(abs(B-g["stationary_topcode_adjusted_birth_flow"]/2.1),
               abs(o["old_queue_B_over_E"]-B/E),
               abs(o["outside_flow_M"]-o["outside_origin_entry_share"]*E),
               abs(o["retention_rho"]-(1-o["outside_origin_entry_share"])*E/B),
               abs(E-o["outside_flow_M"]-o["retention_rho"]*B))>2e-12:
            raise RuntimeError("Small-step exact queue-flow identity failed")
        measurement_gap=o["old_queue_B_over_E"]-summary["old_model_completed_fertility"]/2.1
        if abs(measurement_gap)>2e-8: raise RuntimeError("Small-step original measurement gate failed")
        names,ts,ms,ws,residual,loss=panel.parse_target_rows(panel.read_csv(directory/"target_fit_long.csv"),candidate["candidate"])
        if names!=moments or not np.array_equal(ts,targets) or not np.array_equal(ws,weights):
            raise RuntimeError("Small-step target table differs")
        if not math.isclose(loss,candidate["transition_loss"],abs_tol=1e-9,rel_tol=1e-10):
            raise RuntimeError("Small-step objective does not reconstruct")
        unit=np.asarray(design["unit_vector"],float)
        tasks[i]=dict(unit=unit,residual=residual,models=ms,loss=loss)
        case_rows.append(dict(task=i,design=expected_design,loss_original12=float(residual[mask]@residual[mask]),
            loss_augmented13=loss,young_ownership=float(ms[moments.index("own_rate_2534")]),
            first_birth_rooms=float(ms[moments.index("housing_increment_0to1")]),
            flow_stock_measurement_gap=measurement_gap))
        for k,name in enumerate(moments):
            fits.append(dict(task=i,moment=name,target=ts[k],model=ms[k],gap=ms[k]-ts[k],
                weight=ws[k],loss_contribution=residual[k]**2,production_target=name!="own_rate_2534"))
        parameters.extend(dict(task=i,**r) for r in panel.read_csv(directory/"parameter_table.csv"))
        hashes[str(i)]={name:panel.sha256(directory/name) for name in ("summary.json","target_fit_long.csv","parameter_table.csv")}
    anchor=tasks[1]
    if not np.allclose(anchor["unit"],original[0].unit,rtol=0,atol=2e-12) or not np.allclose(anchor["models"],original[0].model,rtol=0,atol=2e-10):
        raise RuntimeError("Small-step anchor does not reproduce original panel")
    sides=[];derivatives=[]
    for dimension,minus_id,plus_id in ((0,2,3),(3,8,9),(6,14,15)):
        minus,plus=tasks[minus_id],tasks[plus_id]
        for t,sign in ((minus,-1),(plus,1)):
            expected=anchor["unit"].copy();expected[dimension]+=sign*.005
            if not np.allclose(t["unit"],expected,rtol=0,atol=2e-12):
                raise RuntimeError("Small-step pair changes an undeclared coordinate")
        backward=(anchor["residual"]-minus["residual"])/.005
        forward=(plus["residual"]-anchor["residual"])/.005
        central=(forward+backward)/2
        name=original[0].panel["domain"][dimension]["name"]
        norms=np.linalg.norm(backward),np.linalg.norm(forward)
        large=Jlarge[:,dimension]
        sides.append(dict(parameter=name,step=.005,
            forward_backward_cosine=float(forward@backward/max(norms[0]*norms[1],1e-30)),
            side_disagreement_ratio=float(np.linalg.norm(forward-backward)/max(sum(norms),1e-30)),
            central_norm=float(np.linalg.norm(central)),large_step_central_norm=float(np.linalg.norm(large)),
            central_small_large_cosine=float(central@large/max(np.linalg.norm(central)*np.linalg.norm(large),1e-30)),
            central_relative_change=float(np.linalg.norm(central-large)/max(np.linalg.norm(large),1e-30))))
        for k,moment in enumerate(moments):
            derivatives.append(dict(parameter=name,moment=moment,small_step=.005,large_step=.02,
                small_central=central[k],small_backward=backward[k],small_forward=forward[k],large_central=large[k],
                small_raw_moment_derivative=central[k]/np.sqrt(weights[k]),large_raw_moment_derivative=large[k]/np.sqrt(weights[k])))
    for name,rows in (("all_case_objectives",case_rows),("all_target_fits",fits),("all_parameters",parameters),
                      ("step_comparison",sides),("moment_derivatives",derivatives)):
        panel.write_csv(out/f"{name}.csv",rows)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(8,4.4),layout="constrained")
    for column,(dimension,minus_id,plus_id) in enumerate(((0,2,3),(3,8,9),(6,14,15))):
        observations=[original[minus_id-1].model,tasks[minus_id]["models"],anchor["models"],
                      tasks[plus_id]["models"],original[plus_id-1].model]
        for row,(moment,factor,label) in enumerate((("own_rate_2534",100,"Young ownership (%)"),
                ("housing_increment_0to1",1,"First-birth rooms response"))):
            k=moments.index(moment);ax=axes[row,column]
            ax.plot([-.02,-.005,0,.005,.02],[float(v[k])*factor for v in observations],"-o",color="#173E55",markersize=3)
            ax.axhline(float(targets[k])*factor,ls="--",color="#B35A32",lw=1)
            ax.set_ylabel(label,fontsize=8);ax.tick_params(labelsize=7)
            ax.grid(alpha=.15)
            if row==0:ax.set_title(original[0].panel["domain"][dimension]["name"],fontsize=10)
            if row==1:ax.set_xlabel("Transformed-coordinate change",fontsize=8)
    fig.savefig(out/"supplemental_step_stability.pdf",bbox_inches="tight")
    fig.savefig(out/"supplemental_step_stability.png",dpi=180,bbox_inches="tight");plt.close(fig)
    summary=dict(status="complete_seven_case_diagnostic_subset",task_ids=list(ids),anchor_reproduced=True,
        step_comparison=sides,source_hashes=hashes,production_unchanged=True,
        scope="Three paired sensitivity columns at radius .005 compared with .02. This is not a complete new Jacobian, a rank test, or a calibration search. No missing columns are imputed and no candidate is promoted.")
    panel.write_json(out/"summary.json",summary)
    return summary


def global_saving_receipts(root,out):
    """Read the completed conditional policy experiment without model imports."""
    summary=panel.read_json(root/"summary.json")
    if summary.get("status")!="complete" or summary.get("mode")!="equilibrium" or summary.get("production_changed") is not False:
        raise RuntimeError("Global-saving policy panel is not a complete diagnostic comparison")
    if summary["checkpoint_sha256"]!="bbe10a21a843facaf2bceed56e89281e00992d632cddab640ff0c572d3eb494f":
        raise RuntimeError("Global-saving inherited state differs")
    if summary["input_hashes"]["code_bundle_sha256"]!="630ba20bca6a1b54eb4c46aca904c4a087afb8c808b9c7f4660d5fcd316a970e":
        raise RuntimeError("Global-saving production source differs")
    for key,path in (("driver_sha256",ROOT/"code/model/tools/run_e5f_global_saving_quantification.py"),
                     ("oracle_driver_sha256",ROOT/"code/model/tools/run_e5f_independent_numerical_audit.py")):
        if summary[key]!=panel.sha256(path):raise RuntimeError(f"Global-saving {key} differs from reviewed source")
    rows=summary["rows"];policies=("baseline","supply-plus-20","dependent-child-ltv95")
    expected={(method,p,"separately_cleared") for method in ("local","global") for p in policies}
    expected|={("global",p,"fixed_local_policy_price") for p in policies}
    if len(rows)!=9 or {(r["method"],r["policy"],r["pricing"]) for r in rows}!=expected:
        raise RuntimeError("Global-saving case set differs")
    if max(r["adult_households"] for r in rows)-min(r["adult_households"] for r in rows)>2e-10:
        raise RuntimeError("Global-saving common population differs")
    comparisons=[];checkpoints={}
    for r in rows:
        case=f'{r["method"]}_{r["policy"]}_{r["pricing"]}'
        path=root/"cases"/case/"case_state.pkl.gz"
        if panel.sha256(path)!=r["case_checkpoint_sha256"]:raise RuntimeError(f"Case checkpoint {case} differs")
        checkpoints[case]=r["case_checkpoint_sha256"]
        if r["market_cleared"] and r["market_residual"]>2e-4:raise RuntimeError("Market gate failed")
        if r["policy"]=="baseline":continue
        base=next(b for b in rows if b["method"]==r["method"] and b["pricing"]==r["pricing"] and b["policy"]=="baseline")
        comparisons.append(dict(method=r["method"],policy=r["policy"],pricing=r["pricing"],
            births_percent=100*(r["births_per_household"]/base["births_per_household"]-1),
            rooms_percent=100*(r["rooms_per_household"]/base["rooms_per_household"]-1),
            ownership_pp=100*(r["ownership"]-base["ownership"]),
            price_percent=100*(r["asset_price"]/base["asset_price"]-1)))
    out.mkdir(parents=True,exist_ok=True)
    panel.write_csv(out/"all_case_quantities.csv",rows)
    panel.write_csv(out/"policy_effects.csv",comparisons)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.7),layout="constrained")
    for ax,p,title in zip(axes,policies[1:],("Supply +20%","Dependent-child LTV 95%")):
        values=[next(r["births_percent"] for r in comparisons if r["method"]==m and r["policy"]==p and r["pricing"]=="separately_cleared") for m in ("local","global")]
        ax.bar([0,1],values,color=["#173E55","#B35A32"],width=.6)
        ax.set_xticks([0,1],["Production optimizer","Global saving"],fontsize=8)
        ax.set_ylabel("Impact births per household (%)",fontsize=8);ax.set_title(title,fontsize=10)
        ax.axhline(0,color=".3",lw=.6);ax.spines[['top','right']].set_visible(False)
        for x,value in enumerate(values):ax.annotate(f"{value:.5f}%",(x,value),xytext=(0,5 if value>=0 else -12),textcoords="offset points",ha="center",fontsize=8)
        ax.margins(y=.25)
    fig.savefig(out/"supplemental_global_saving.pdf",bbox_inches="tight")
    fig.savefig(out/"supplemental_global_saving.png",dpi=180,bbox_inches="tight");plt.close(fig)
    result=dict(status="complete_verified_conditional_policy_readout",comparisons=comparisons,
        summary_sha256=panel.sha256(root/"summary.json"),checkpoint_hashes=checkpoints,
        max_cleared_market_residual=max(r["market_residual"] for r in rows if r["market_cleared"]),
        max_budget_excess_share=max(r["budget_excess_share"] for r in rows),production_unchanged=True,
        scope="Common inherited 2023 pre-choice distribution, unchanged calibrated parameters. Each method clears current markets separately. Fixed_local_policy_price rows use the corresponding local-method policy and baseline prices; those markets need not clear under global saving.")
    panel.write_json(out/"summary.json",result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel",type=Path,default=PANEL)
    parser.add_argument("--outdir",type=Path,default=OUT)
    parser.add_argument("--smallstep-panel",type=Path,help="Explicit completed seven-case subset; omitted by default.")
    parser.add_argument("--global-saving-dir",type=Path,help="Explicit completed nine-case numerical comparison.")
    args=parser.parse_args();out=args.outdir
    out.mkdir(parents=True,exist_ok=True)
    tasks,moments,targets,weights,provenance=panel.load_panel(args.panel)
    assert provenance["contract"]["source_sha256"]=="0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6"
    if provenance["contract"]["target_fingerprint"]!="186d692167d2d0905b2621f5dc31a9f9edca2fe3e7f9c9a3c20d35201442f0ac":
        raise RuntimeError("Unexpected target fingerprint")
    if provenance["contract"]["code_bundle_sha256"]!="b033604b0b647f200bb03c5260fa476bcf61cccbb1689b04e76c745aa8f216a6":
        raise RuntimeError("Unexpected source bundle")
    J,side,back,fwd=panel.central_jacobian(tasks[0],[(tasks[1+2*i],tasks[2+2*i]) for i in range(11)])
    iy=moments.index("own_rate_2534");mask=np.arange(len(moments))!=iy
    comparison=[];fit=[];parameters=[];renewal=[];hash_contract=None
    for t in tasks:
        s=t.summary;o=s["renewal_accounting_old_state"];g=s["stationary_operator_gates"]
        contract=ridge.scientific_contract(s)
        serialized=json.dumps(contract,sort_keys=True)
        if hash_contract is None: hash_contract=serialized
        elif serialized!=hash_contract: raise RuntimeError("Complete scientific contracts differ")
        ridge.validate_candidate_eligibility(s)
        E=o["old_entry_flow_E"];B=o["old_queue_mature_flow_B"]
        gaps=dict(queue_from_recorded_birth_flow=B-g["stationary_topcode_adjusted_birth_flow"]/2.1,
            reported_ratio=o["old_queue_B_over_E"]-B/E,
            outside_flow=o["outside_flow_M"]-o["outside_origin_entry_share"]*E,
            retention=o["retention_rho"]-(1-o["outside_origin_entry_share"])*E/B,
            renewal=E-o["outside_flow_M"]-o["retention_rho"]*B)
        if max(map(abs,gaps.values()))>2e-12: raise RuntimeError("Exact flow accounting failed")
        measurement_gap=o["old_queue_B_over_E"]-s["old_model_completed_fertility"]/2.1
        if abs(measurement_gap)>2e-8: raise RuntimeError("Original driver measurement gate failed")
        try:
            ridge.validate_renewal_accounting(s); rejection=""
        except ridge.ContractError as error: rejection=str(error)
        renewal.append(dict(task=t.task_id,**gaps,flow_stock_measurement_gap=measurement_gap,
             stationary_measurement_max_abs_gap=s["stationary_measurement_max_abs_gap"],
             one_step_stationary_nesting_l1=g["one_step_constant_path_nesting_l1"],
             strict_planner_passed=not bool(rejection),strict_planner_rejection=rejection))
        loss12=float(t.residual[mask]@t.residual[mask])
        comparison.append(dict(task=t.task_id,design=t.panel["design"],loss_original12=loss12,
            loss_augmented13=t.loss,young_ownership=t.model[iy],young_gap_pp=100*(t.model[iy]-targets[iy]),
            first_birth_rooms=t.model[moments.index("housing_increment_0to1")],
            mean_rooms=t.model[moments.index("aggregate_mean_occupied_rooms_18_85")],
            own_rate=t.model[moments.index("own_rate")],young_loss=t.residual[iy]**2))
        for k,name in enumerate(moments):
            fit.append(dict(task=t.task_id,moment=name,target=targets[k],model=t.model[k],
                 gap=t.model[k]-targets[k],weight=weights[k],loss_contribution=t.residual[k]**2,
                 production_target=name!="own_rate_2534"))
        for row in panel.read_csv(t.summary_path.parent/"parameter_table.csv"):
            parameters.append(dict(task=t.task_id,**row))
    ranks=[];directions=[]
    for label,A in (("original12",J[mask]),("augmented13",J)):
        u,s,v=np.linalg.svd(A,full_matrices=False)
        ranks.append(dict(target_system=label,**panel.rank_report(A)))
        for i,row in enumerate(tasks[0].panel["domain"]):
            directions.append(dict(target_system=label,parameter=row["name"],
                weakest_loading=v[-1,i],squared_loading=v[-1,i]**2))
    theta1_index=next(i for i,d in enumerate(tasks[0].panel["domain"]) if d["name"]=="theta1")
    conditional_rank=panel.rank_report(np.delete(J[mask],theta1_index,axis=1))
    # Both systems use precisely the same parameter coordinate metric. Adding
    # a row cannot reduce J'J in Loewner order; a larger condition number alone
    # must not be interpreted as information destruction.
    info_increment=J.T@J-J[mask].T@J[mask]
    if not np.allclose(info_increment,np.outer(J[iy],J[iy]),rtol=1e-10,atol=1e-9):
        raise RuntimeError("Information add-up failed")
    for name,rows in (("all_case_objectives",comparison),("all_target_fits",fit),
                      ("all_parameters",parameters),("renewal_identity_receipts",renewal),
                      ("weakest_directions",directions),("coordinate_side_consistency",side)):
        panel.write_csv(out/(name+".csv"),rows)
    for label,A,names in (("original12",J[mask],[m for m in moments if m!="own_rate_2534"]),
                           ("augmented13",J,moments)):
        panel.write_csv(out/("jacobian_"+label+".csv"),[dict(moment=name,**{d["name"]:A[k,i] for i,d in enumerate(tasks[0].panel["domain"])}) for k,name in enumerate(names)])
    summary=dict(status="complete_read_only_receipts",rank_comparison=ranks,
        conditional_rank_if_theta1_externally_fixed=conditional_rank,
        conditional_rank_scope="Column deletion at the current point only, not an external restriction chosen or a new calibration.",
        ownership_measurement=ownership_measurement(out,targets[moments.index("aggregate_mean_occupied_rooms_18_85")]),
        strict_planner_pass_count=sum(r["strict_planner_passed"] for r in renewal),
        strict_planner_fail_count=sum(not r["strict_planner_passed"] for r in renewal),
        max_exact_flow_gap=max(abs(r[k]) for r in renewal for k in gaps),
        max_flow_stock_measurement_gap=max(abs(r["flow_stock_measurement_gap"]) for r in renewal),
        best_original12=min(comparison,key=lambda r:r["loss_original12"]),
        best_augmented13=min(comparison,key=lambda r:r["loss_augmented13"]),
        production_unchanged=True,proposals_written=False,
        interpretation="The augmented row adds information mostly to an already strong ownership-premium direction; the weakest bequest-shifter direction hardly changes. Both matrices have algebraic rank 11 and relative rank 10 at 1e-3. This is weak local identification, not a proof of structural underidentification.",
        provenance=provenance)
    if args.smallstep_panel:
        summary["smallstep"]=smallstep_receipts(args.smallstep_panel,out.parent/"smallstep_receipts",tasks,moments,targets,weights,J)
    if args.global_saving_dir:
        summary["global_saving"]=global_saving_receipts(args.global_saving_dir,out.parent/"global_saving_receipts")
    panel.write_json(out/"summary.json",summary)
    print(json.dumps({k:v for k,v in summary.items() if k!="provenance"},indent=2))


if __name__=="__main__":main()
