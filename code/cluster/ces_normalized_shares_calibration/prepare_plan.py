"""Prepare four deterministic CES-share starts; this never initializes the model."""
import hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/experiments/ces_normalized_shares/overnight_v1"
ANCHOR = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json"

def canonical(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    anchor=json.loads(ANCHOR.read_text())
    if anchor.get("status") != "selected_numerically_verified" or anchor.get("chain") != 13 or anchor.get("arm") != "alternative":
        raise SystemExit("Canonical post-interest chain-13 source mismatch")
    selected=dict(anchor["selected"]["parameters"])
    bounds={r["parameter"]:[float(r["lower"]),float(r["upper"])] for r in anchor["parameters"] if r.get("parameter") in selected and r.get("lower") not in ("",None)}
    if selected.pop("h_P", None) is None or "h_P" not in bounds: raise SystemExit("Anchor h_P coordinate missing")
    bounds.pop("h_P")
    bounds["delta_alpha_jump"]=[0.0,0.25]
    bounds["delta_alpha"]=[0.0,0.25]
    coordinates=list(selected)+["delta_alpha_jump","delta_alpha"]
    if len(coordinates)!=11 or len(bounds)!=11: raise SystemExit("Expected nine inherited coordinates plus jump and slope")
    target=[{k:r[k] for k in ("moment","target","weight","role")} for r in anchor["target_fit"]]
    wealth=[r for r in target if r["moment"]=="wealth_earnings"]
    family=[r for r in target if r["moment"]=="family_rooms"]
    if family != [{"moment":"family_rooms","target":"0.38509964969278165","weight":"0.0","role":"validation"}]: raise SystemExit("Canonical family_rooms row drift")
    family[0].update(weight="280.52808370152104",role="scored")
    if len(target)!=14 or len([r for r in target if r["role"]=="scored"])!=11 or wealth != [{"moment":"wealth_earnings","target":"6.92658379107299","weight":"7.595098472533724","role":"scored"}]:
        raise SystemExit("Retained same-14 old-wealth target contract drift")
    starts=[]
    for jump,slope in ((.07780442689806688,.03897536437154123),(.04,.02),(.12,.04),(.18,.08)):
        p=dict(selected, delta_alpha_jump=jump, delta_alpha=slope)
        if not all(bounds[k][0]<=p[k]<=bounds[k][1] for k in coordinates): raise SystemExit("start outside bound")
        starts.append(p)
    plan=dict(status="prepared_not_launched", name="ces_normalized_shares_overnight_20261003_v1",
      reference="canonical post-interest soft chain 13; retained 6.926584 wealth target; family_rooms promoted to scored (11 of 14)",
      source_checkpoint=str(ANCHOR.relative_to(ROOT)),source_checkpoint_sha256=sha(ANCHOR),selected_source_sha256=sha(ANCHOR),
      coordinates=coordinates,bounds=bounds,starts=starts,start_provenance=[{"chain":i,"delta_alpha_jump":p["delta_alpha_jump"],"delta_alpha":p["delta_alpha"],"base":"chain13; first start is historical lastweek estimate, not adoption"} for i,p in enumerate(starts)],
      target_contract=target,target_fingerprint=canonical(target),weight_fingerprint=canonical(dict(base_contract=target,multipliers={})),
      contract_id="ces_normalized_jump_slope_family_rooms_v1", economics=dict(utility="Q=c^a*s^(1-a)/(a^a(1-a)^(1-a)); no r* or alpha0 numerator", alpha="author-selected isolated test: alpha(m)=.733 if m=0, otherwise clip(.733-delta_alpha_jump-delta_alpha*m,.05,.95)", unchanged="all other chain-13 economics, grid, fixed inputs and raw-unit utility costs; no estateA and no new birth menu"),
      family_rooms_provenance=dict(definition="E[min(ROOMS,9)|NCHILD>=3,YNGCH<18]-E[min(ROOMS,9)|NCHILD=1-2,YNGCH<18]",sample="ACS 2005/06 national heads 30-55, HHWT>0, no DUE restriction",estimator="weighted mean difference; no FE; national SE unavailable",weight_note="280.52808370152104 is inherited 42-metro bootstrap weight, not national inverse variance",mapping_warning="Native observer remains a model-dependent proxy: model m can differ from resident own-child count and does not guarantee ACS minor-child restriction.",builder="output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py; output/model/e5f_matched_pf_20260909a/design_research/housing/summarize_early_housing.py",pointer="specification_followup/housing_profiles_v1/full/target_recomputed.json#/recomputed/national/family_rooms"),
      resource_budget=dict(chains=4,cpus_per_chain=1,memory_GiB_per_chain=24,wall_seconds_per_chain=21600,max_objective_calls_per_chain=500,max_lifecycle_per_GE=32,final_native_reserve_seconds=1800,smoke_objective_calls=2,smoke_wall_seconds=5400,total_planned_max_GE=2000,total_planned_max_lifecycle=64000,observed_reference="two evaluations plus postcheck: about 350--440 seconds; CES shape requires smoke measurement"),
      method="bounded Nelder-Mead; selected postcheck is a fresh child and native exact repeat", no_auto_retry=True,no_adoption=True,production_release_required=True)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"start_plan.json").write_text(json.dumps(plan,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":plan["status"],"starts":len(starts),"fingerprint":plan["target_fingerprint"]}))
if __name__=="__main__": main()
