"""Build the 13-row validation table from a strict native 2023 readout.

This is deliberately table-only: it does not import or invoke the historical
patch driver, stationary solvers, or sequence-fertility routines.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

VINTAGES = {
    "completed_fertility": "CPS 2024 (nearest fertility supplement)",
    "childlessness": "CPS 2024 (nearest fertility supplement)",
    "exactly_one": "CPS 2024 (nearest fertility supplement)",
    "first_birth_age": "NCHS 2023",
    "first_birth_share30": "NCHS 2023",
    "mean_rooms": "ACS 2023",
    "ownership_30_55": "ACS 2023",
    "first_birth_rooms": "PSID pooled; retained Sun–Abraham receipt",
    "family_rooms": "ACS 2023",
    "recent_parent": "ACS 2023",
    "wealth_earnings": "PSID 2005–2019",
    "bequest_wealth": "External benchmark",
    "old_dispersion": "PSID 1984–2019",
}

def _rows(model, empirical):
    f, h = model["fertility_stock_timing"], model["housing_wealth"]["moments"]
    flow = model["fertility"]
    third, explicit = flow["birth_flow_third_bin_entry"], flow["birth_flow_explicit"]
    adjusted = flow["birth_flow_topcode_adjusted"]
    reps = [3 + (a-e)/t for a,e,t in zip(adjusted, explicit, third) if t > 1e-12]
    if not reps or max(reps)-min(reps) > 1e-10: raise ValueError("inconsistent top-bin representative")
    shares = f["parity_shares_40_44"]
    completed = shares["1"] + 2*shares["2"] + reps[0]*shares["3plus"]
    def add(key, label, data, value, scale=1, decimals=2, note=""):
        if not all(math.isfinite(float(x)) for x in (data,value)):
            raise ValueError('Nonfinite validation moment: '+key)
        rows.append(dict(moment_key=key, moment=label, data=scale*float(data),
            model=scale*float(value), gap=scale*(float(value)-float(data)),
            data_vintage=VINTAGES[key], model_year=2023, decimals=decimals,
            measurement_note=note))
    rows=[]
    cps, nchs, acs, psid, wealth, beq = (empirical[k] for k in ("cps","nchs","acs","psid","wealth","bequest"))
    add("completed_fertility", "Completed fertility, ages 40–44", cps["tfr"], completed, note="Children ever born stock; distinct from period fertility.")
    add("childlessness", "Childless, ages 40–44 (%)", cps["childless_rate"], f["moments"]["childless_rate_40_44"], 100, 1)
    add("exactly_one", "Exactly one child among mothers, 40–44 (%)", cps["parity_share_1"]/(1-cps["childless_rate"]), f["moments"]["exactly_one_among_mothers_40_44"], 100, 1)
    add("first_birth_age", "Mean age at first birth (years)", nchs["mean_age"], f["moments"]["period_mean_age_first_birth"], note="Model 2024–2027 flow versus annual 2023 NCHS births.")
    add("first_birth_share30", "First births at age 30+ (%)", nchs["share30"], f["moments"]["period_share_first_births_age30plus"], 100, 1, note="Model 2024–2027 flow versus annual 2023 NCHS births.")
    add("mean_rooms", "Mean occupied rooms (capped at 9)", acs["mean_rooms"], h["aggregate_mean_occupied_rooms_capped9_18_85"])
    add("ownership_30_55", "Ownership, heads 30–55 (%)", acs["ownership"], h["own_rate_30_55"], 100, 1)
    add("first_birth_rooms", "First-birth room response, −1 to +3", psid["first_birth_rooms"], model["dated_first_birth_rooms"]["housing_response"], 1, 3, note="Sun–Abraham empirical contrast; model matched branch from 2019 into 2023.")
    add("family_rooms", "Rooms: 3+ versus 1–2 resident children", acs["family_rooms"], h["prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9"], 1, 3, note="Model dependent counts proxy resident own children under 18.")
    add("recent_parent", "Recent-parent ownership gap (pp)", acs["recent_parent"], model["recent_parent"]["model_value"], 100, 1, note="Model current birth into an empty dependent home versus currently empty home; retained flow proxy, not exact ACS oldest-child-age reconstruction.")
    add("wealth_earnings", "Wealth / annual gross earnings", wealth, h["aggregate_wealth_to_annual_gross_labor_earnings"])
    add("bequest_wealth", "Annual bequests / wealth (%)", beq, h["annual_bequest_flow_to_aggregate_wealth"], 100, 2, note="External historical restriction; not a 2023 observation.")
    add("old_dispersion", "Wealth/income p90 / median, ages 76–84", psid["old_dispersion"], h["old_total_wealth_to_annual_income_p90_p50_7684"], note="Beginning-period wealth of living households; model pension-income proxy versus empirical family income.")
    if len(rows) != 13: raise AssertionError("expected exactly 13 rows")
    return rows

def load_empirical():
    def records(path):
        with Path(path).open() as f: return list(csv.DictReader(f))
    d=ROOT/'output/model/e5f_matched_pf_20260909a/design_research'
    cps_path=ROOT/'code/data/cps_fertility/output/cps_fertility_targets.csv'; cps={r['moment_key']:float(r['estimate']) for r in records(cps_path)}
    nchs_path=ROOT/'code/data/nchs_natality_timing/first_birth_counts_year_age.csv'; births=[(int(r['age']),float(r['n_first_births'])) for r in records(nchs_path) if int(r['year'])==2023]
    midpoint=lambda a:20 if a<22 else 44 if a>=42 else 20+4*((a-18)//4); total=sum(n for _,n in births)
    acs_path=d/'housing/early_housing_target_candidates.csv'; acs={r['moment']:float(r['point']) for r in records(acs_path) if r['window']=='2023'}
    room_path=ROOT/'code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv'; room_rows=records(room_path); room=room_rows[0]
    wealth_path=d/'wealth/aggregate_wealth_results.csv'; wealth=next(float(r['estimate']) for r in records(wealth_path) if r['window']=='pooled_2005_2019')
    old_path=d/'wealth/old_wealth_results.csv'; old=next(float(r['estimate']) for r in records(old_path) if r['window']=='pooled_1984_2019' and r['moment']=='old_p90_p50')
    beq_path=ROOT/'output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_weights.csv'; beq=next(float(r['target']) for r in records(beq_path) if r['restriction_id']=='bequest_wealth')
    paths={"cps":cps_path,"nchs":nchs_path,"acs":acs_path,"psid_room":room_path,"wealth":wealth_path,"old":old_path,"bequest":beq_path}
    empirical={"cps":cps,"nchs":{"mean_age":sum(n*midpoint(a) for a,n in births)/total,"share30":sum(n for a,n in births if a>=30)/total},"acs":{"mean_rooms":acs['aggregate_mean_occupied_rooms_capped9_18_85'],"ownership":acs['own_rate_30_55'],"family_rooms":acs['prime30_55_resident_3plus_minus_1to2_rooms_capped9'],"recent_parent":acs['recent_parent_minus_no_resident_child_ownership_30_55']},"psid":{"first_birth_rooms":float(room['estimate']),"old_dispersion":old},"wealth":wealth,"bequest":beq}
    meta={"paths":{k:{"path":str(v),"sha256":hashlib.sha256(Path(v).read_bytes()).hexdigest()} for k,v in paths.items()},"predicates":{"cps":"all rows; moment_key lookup","nchs":"year == 2023; age-cell midpoint weighting","acs":"window == 2023; moment lookup","psid_room":"first row of retained target_receipt.csv; estimator/sample defined by receipt","wealth":"window == pooled_2005_2019","old":"window == pooled_1984_2019 and moment == old_p90_p50","bequest":"restriction_id == bequest_wealth"},"psid_room_receipt":room}
    return empirical,meta

def build_table(model_path: Path, out_dir: Path, empirical=None, empirical_meta=None, fixture=False):
    model_path = Path(model_path).resolve(); model = json.loads(model_path.read_text())
    if model.get("calendar_year") != 2023: raise ValueError("strict reader model must be calendar year 2023")
    verification=model_path.with_name('verification.json')
    check={}
    if not fixture:
        if empirical is not None:raise ValueError('Empirical overrides require explicit fixture mode')
        if not verification.is_file():raise ValueError('Strict-reader verification.json is required')
        check=json.loads(verification.read_text())
        if (check.get('status')!='PASS' or check.get('finite_converged') is not True
                or check.get('verification_method')!='native_saved_snapshot_aggregate_match'
                or check.get('year')!=2023 or not check.get('historical_fit_status',{}).get('complete')):
            raise ValueError('Complete finite historical readout verification is required')
        if not model.get('forecast_receipt_sha256') or model['forecast_receipt_sha256']!=check.get('root_receipt_sha256'):
            raise ValueError('Model/readout receipt linkage is missing or mismatched')
    if empirical is None: empirical, empirical_meta = load_empirical()
    rows = _rows(model, empirical); out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir/"validation_2023.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    manifest={"rows":13,"model_source":str(model_path),"model_source_sha256":hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "data_vintages":VINTAGES,"sequence_windows":{"2007":"2008–2011","2011":"2012–2015","2015":"2016–2019","2019":"2020–2023"},
        "forecast_window":"2024–2027","property_tax":"1% annual; equally rebated; PAYGO balanced","empirical_values_changed":False,
        "horizon_verified":False,"deterministic_fixture":bool(fixture),"empirical_provenance":empirical_meta or {},
        "readout_verification":check,
        "readout_verification_sha256":hashlib.sha256(verification.read_bytes()).hexdigest() if check else None}
    (out_dir/"validation_2023_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return rows

def main():
    p=argparse.ArgumentParser(); p.add_argument("--model", type=Path, required=True); p.add_argument("--out", type=Path, required=True)
    a=p.parse_args(); build_table(a.model,a.out)
if __name__ == "__main__": main()
