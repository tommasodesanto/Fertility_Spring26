#!/usr/bin/env python3
"""Private staging and aggregate-only verification for the PSID A2h income sensitivity."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, shutil, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_REFERENCE = HERE / "sa_rooms_first_birth_v2.do"
DEFAULT_OUTPUT = HERE / "output" / "first_birth_rooms_income_sensitivity_v1"
REFERENCE_SHA = "49ebdb6c4c780eaa6eb00be0c66c0f8523af85ed0f2c049c48d286f2281d8b8e"
PREPARED_SHA = "9fd181c48e1951a7a0052836926260a89d056d666aad3e49e2ae9cca05f1ba00"
FITS = ("baseline_full", "common_no_income", "common_income")
EVENTS = (("Wleft", "≤−8"), ("Wm2", "−7/−6"), ("Wm1", "−5/−4"),
          ("baseline", "−3/−2"), ("Wp1", "−1/0"), ("Wp2", "+1/+2"),
          ("Wp3", "+3/+4"), ("Wp4", "+5/+6"), ("Wp5", "+7/+8"),
          ("Wp6", "+9/+10"), ("Wright", "≥+11"))

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def identical(left: Path, right: Path) -> bool:
    if left.stat().st_size != right.stat().st_size: return False
    with left.open("rb") as a, right.open("rb") as b:
        while True:
            x, y = a.read(1024 * 1024), b.read(1024 * 1024)
            if x != y: return False
            if not x: return True

def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding="utf-8")

def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f: return list(csv.DictReader(f))

def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1: raise RuntimeError(f"required literal hook occurs {text.count(old)} times")
    return text.replace(old, new, 1)

def extract_income(a: argparse.Namespace) -> None:
    source, private = Path(a.source_dta).resolve(), Path(a.private_dir).resolve()
    if not source.is_file(): raise SystemExit(f"source microdata not found: {source}")
    private.mkdir(parents=True, exist_ok=True); os.chmod(private, 0o700)
    income, do, log = private / "income.dta", private / "extract_income.do", private / "income_extraction.log"
    before = source.stat()
    write(do, f'''version 17.0
clear all
set more off
set processors 1
log using "{log}", replace text
use ID year INCFAMR using "{source}", clear
isid ID year
save "{income}", replace
di "INCOME_EXTRACTION_PASS"
log close
''')
    meta = {"source_dta": str(source), "source_sha256": sha(source), "source_stat_before": [before.st_size, before.st_mtime_ns],
            "income_dta": str(income), "income_clock": "FU total family income, PCE-adjusted to 2022 USD, tax year interview_year-1"}
    write(private / "income_source_metadata.json", json.dumps(meta, indent=2) + "\n")
    if not a.run: return
    subprocess.run([a.stata, "-bq", "do", str(do)], cwd=private, check=True, timeout=300)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns): raise RuntimeError("source stat changed during extraction")
    if not income.is_file() or "INCOME_EXTRACTION_PASS" not in log.read_text(errors="replace"): raise RuntimeError("income extraction marker missing")
    os.chmod(income, 0o600); meta.update({"income_sha256": sha(income), "source_stat_after": [after.st_size, after.st_mtime_ns]})
    write(private / "income_source_metadata.json", json.dumps(meta, indent=2) + "\n")

def income_hook(with_income: bool) -> str:
    cov = "covariates(i.AGEREP i.EDUYEAR log_family_income)" if with_income else "covariates(i.AGEREP i.EDUYEAR)"
    return f'''* INCOME_SENSITIVITY_HOOK: canonical status/support already constructed.
unab master_variables : _all
tempfile originalmaster
preserve
    sort ID year
    keep `master_variables'
    save `originalmaster'
restore
local master_rows = _N
merge 1:1 ID year using "__INCOME_DTA__", keep(master match)
assert inlist(_merge,1,3)
assert _N == `master_rows'
drop _merge
isid ID year
sort ID year
cf `master_variables' using `originalmaster'
gen byte income_missing = missing(INCFAMR)
gen byte income_nonpositive = !missing(INCFAMR) & INCFAMR <= 0
gen byte income_topcoded = INCFAMR == 9999999
gen byte income_valid = INCFAMR > 0 & INCFAMR < . & INCFAMR != 9999999
gen double log_family_income = ln(INCFAMR) if income_valid
quietly count if income_missing
local income_missing_rows = r(N)
quietly count if income_nonpositive
local income_nonpositive_rows = r(N)
quietly count if income_topcoded
local income_topcoded_rows = r(N)
quietly count if income_valid
local income_valid_rows = r(N)
quietly count if !income_valid
local income_excluded_rows = r(N)
* Restriction is intentionally after canonical support and immediately before estimation.
keep if income_valid
bysort f_c_y: egen long common_reference_rows = total(inrange(K,`baseline_lo',`baseline_hi'))
quietly count if !missing(f_c_y) & control == 0 & common_reference_rows == 0
assert r(N) == 0
drop common_reference_rows
eventstudyinteract rooms `dummies' `weightspec', ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(control) ///
    {cov}
'''

def staged(source: str, fit: str, income: Path) -> str:
    if fit == "baseline_full": return source
    anchor = "eventstudyinteract rooms `dummies' `weightspec', ///\n    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(control) ///\n    covariates(i.AGEREP i.EDUYEAR)"
    out = replace_once(source, anchor, income_hook(fit == "common_income").replace("__INCOME_DTA__", str(income)))
    receipt = "gen double runtime_seconds = `seconds'"
    extra = receipt + '''
gen long income_missing_rows = `income_missing_rows'
gen long income_nonpositive_rows = `income_nonpositive_rows'
gen long income_topcoded_rows = `income_topcoded_rows'
gen long income_valid_rows = `income_valid_rows'
gen long income_excluded_rows = `income_excluded_rows'
gen str100 income_definition = "ln(INCFAMR); valid iff 0<INCFAMR<. and INCFAMR!=9999999"'''
    return replace_once(out, receipt, extra)

def stage(a: argparse.Namespace) -> None:
    task, ref, income = Path(a.task_root).resolve(), Path(a.reference_estimator).resolve(), Path(a.income_dta).resolve()
    sample, ado = Path(a.analysis_sample).resolve(), Path(a.ado_dir).resolve()
    if not all(x.is_file() for x in (ref, income, sample)): raise SystemExit("reference, prepared sample, or private income.dta is missing")
    if sha(ref) != a.reference_sha256 or sha(sample) != a.analysis_sample_sha256: raise SystemExit("reference or prepared-data SHA256 pin mismatch")
    stage_dir, results = task / "staged", task / "results"; stage_dir.mkdir(parents=True, exist_ok=True); results.mkdir(parents=True, exist_ok=True)
    source = ref.read_text(encoding="utf-8")
    pins = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha(Path(__file__)), "reference_estimator": str(ref), "reference_sha256": sha(ref), "analysis_sample": str(sample), "analysis_sample_sha256": sha(sample), "income_dta": str(income), "income_sha256": sha(income), "ado_dir": str(ado), "entry": str(task / "entry.do"), "processors": 8, "fits": FITS}
    for fit in FITS:
        (results / fit).mkdir(parents=True, exist_ok=True); write(stage_dir / f"{fit}.do", staged(source, fit, income))
    if not identical(stage_dir / "baseline_full.do", ref): raise RuntimeError("baseline stage is not byte-identical")
    entry = f'''clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "{ado}"
adopath ++ "{ado}"
mata: mata mlib index
log using "{task / 'stata_run.log'}", replace text
'''
    for fit in FITS:
        entry += f'''use "{sample}", clear
capture noisily do "{stage_dir / (fit + '.do')}" A2h "{results / fit}"
if _rc exit _rc
file open done using "{task / 'latest_completed_case.txt'}", write replace
file write done "{fit}" _n
file close done
di "INCOME_FIT_PASS {fit}"
'''
    entry += 'di "ROOMS_INCOME_THREE_FIT_PASS"\nlog close\nexit 0\n'
    write(task / "entry.do", entry); pins["entry_sha256"] = sha(task / "entry.do")
    write(task / "run_config.json", json.dumps(pins, indent=2) + "\n")
    sums = [(sha(stage_dir / f"{f}.do"), f"staged/{f}.do") for f in FITS] + [(sha(task / "entry.do"), "entry.do"), (sha(Path(__file__)), str(Path(__file__).resolve()))]
    write(task / "SHA256SUMS", "".join(f"{h}  {n}\n" for h,n in sums))

def validate_stage(a: argparse.Namespace) -> None:
    task = Path(a.task_root).resolve(); cfg = json.loads((task / "run_config.json").read_text())
    for key in ("driver", "reference_estimator", "analysis_sample", "income_dta", "entry"):
        p = Path(cfg[key]); expect = cfg.get("reference_sha256" if key=="reference_estimator" else key + "_sha256")
        if not p.is_file() or (expect and sha(p) != expect): raise RuntimeError(f"pin failed: {key}")
    if not identical(task / "staged/baseline_full.do", Path(cfg["reference_estimator"])): raise RuntimeError("baseline source changed")

def verify_fit(directory: Path) -> None:
    import numpy as np
    c, v = rows(directory / "coefficients.csv"), rows(directory / "covariance.csv")
    if not c or any(not math.isfinite(float(x[k])) for x in c for k in ("estimate", "standard_error")): raise RuntimeError(f"{directory}: nonfinite coefficients")
    names = [x["coefficient"] for x in c]; idx = {n:i for i,n in enumerate(names)}; m = np.full((len(names),len(names)), np.nan)
    for x in v: m[idx[x["coefficient_i"]], idx[x["coefficient_j"]]] = float(x["covariance"])
    if not np.isfinite(m).all() or not np.allclose(m,m.T,atol=1e-10) or np.linalg.eigvalsh((m+m.T)/2).min() < -1e-8: raise RuntimeError(f"{directory}: covariance gate failed")
    for i,x in enumerate(c):
        if abs(float(x["standard_error"]) - math.sqrt(max(0,m[i,i]))) > 1e-7: raise RuntimeError(f"{directory}: SE/diagonal mismatch")
    r = rows(directory / "fit_receipt.csv")[0]
    if int(float(r["fitted_unsupported_rows"])) != 0: raise RuntimeError(f"{directory}: unsupported fitted cohort")
    if not math.isfinite(float(r["headline_effect"])) or not math.isfinite(float(r["headline_se"])): raise RuntimeError(f"{directory}: headline not finite")
    point=next(x for x in c if x["coefficient"]==r["headline"])
    for receipt_key,point_key in (("headline_effect","estimate"),("headline_se","standard_error")):
        if abs(float(r[receipt_key])-float(point[point_key]))>1e-7: raise RuntimeError(f"{directory}: headline/coefficients mismatch")

def finalize_private(a: argparse.Namespace) -> None:
    task, result = Path(a.task_root).resolve(), Path(a.task_root).resolve() / "results"
    for fit in FITS: verify_fit(result / fit)
    hashes = {}
    for fit in FITS:
        key = result / fit / "private_sample_keys.csv"
        if not key.is_file(): raise RuntimeError(f"{fit}: private keys absent before finalization")
        hashes[fit] = sha(key); write(result / fit / "sample_key_hashes.json", json.dumps({"sorted_ID_year_sha256":hashes[fit]},indent=2)+"\n")
    if hashes["common_no_income"] != hashes["common_income"]: raise RuntimeError("common final e(sample) key digests differ")
    for fit in FITS:
        (result / fit / "private_sample_keys.csv").unlink(); write(result / fit / "completion.txt", "pass\n")

def same_receipt(actual: dict[str,str], canonical: dict[str,str]) -> dict[str, object]:
    skip={"runtime_seconds"}; diff={}
    if set(actual) != set(canonical): diff["keys"]=[sorted(actual),sorted(canonical)]
    for k in set(actual)&set(canonical)-skip:
        try: ok=abs(float(actual[k])-float(canonical[k]))<=1e-7
        except ValueError: ok=actual[k]==canonical[k]
        if not ok: diff[k]=[actual[k],canonical[k]]
    return diff

def collect_verify(a: argparse.Namespace) -> None:
    import matplotlib.pyplot as plt
    task, result, dest, canonical = Path(a.task_root).resolve(), Path(a.task_root).resolve()/"results", Path(a.destination).resolve(), Path(a.canonical_dir).resolve()
    run=json.loads((task/"run_receipt.json").read_text())
    if run.get("status")!="pass" or run.get("exit_code")!=0: raise RuntimeError("run receipt failed")
    verification={"fits":{},"baseline_differences":{},"run_receipt":run}
    for fit in FITS:
        d=result/fit
        for n in ("completion.txt","coefficients.csv","covariance.csv","fit_receipt.csv","sample_key_hashes.json","fitted_support.csv","input_support.csv"):
            if not (d/n).is_file(): raise RuntimeError(f"missing {fit}/{n}")
        if (d/"completion.txt").read_text().strip()!="pass": raise RuntimeError(f"{fit}: completion gate failed")
        verify_fit(d); verification["fits"][fit]="pass"
    log=(task/"stata_run.log").read_text(errors="replace")
    if log.count("ROOMS_V2_ARM_PASS A2h") != 3:
        raise RuntimeError("expected exactly three canonical Stata completion markers")
    for fit in FITS:
        if f"INCOME_FIT_PASS {fit}" not in log: raise RuntimeError("Stata completion marker missing")
    if "ROOMS_INCOME_THREE_FIT_PASS" not in log: raise RuntimeError("three-fit marker missing")
    h={f:json.loads((result/f/"sample_key_hashes.json").read_text())["sorted_ID_year_sha256"] for f in FITS}
    if h["baseline_full"] != a.expected_sample_key_sha or h["common_no_income"] != h["common_income"]: raise RuntimeError("sample-key digest gate failed")
    aa={x["coefficient"]:x for x in rows(result/"baseline_full"/"coefficients.csv")}; bb={x["coefficient"]:x for x in rows(canonical/"coefficients.csv")}
    if aa.keys()!=bb.keys(): raise RuntimeError("baseline coefficient keys differ")
    for k in aa:
        for col in ("estimate","standard_error"):
            if abs(float(aa[k][col])-float(bb[k][col]))>1e-7: verification["baseline_differences"][f"{k}:{col}"]=[aa[k][col],bb[k][col]]
    rd=same_receipt(rows(result/"baseline_full"/"fit_receipt.csv")[0],rows(canonical/"fit_receipt.csv")[0]); verification["baseline_differences"].update(rd)
    if verification["baseline_differences"]: raise RuntimeError("baseline differs from canonical receipt/coefficients")
    dest.mkdir(parents=True,exist_ok=True)
    comp=[]; curves=[]; readable={"baseline_full":"Original sample","common_no_income":"Income-observed sample","common_income":"Income-observed sample + income"}
    for fit in FITS:
        d=result/fit; shutil.copy2(d/"coefficients.csv",dest/f"{fit}_coefficients.csv"); shutil.copy2(d/"covariance.csv",dest/f"{fit}_covariance.csv"); shutil.copy2(d/"fit_receipt.csv",dest/f"{fit}_fit_receipt.csv"); shutil.copy2(d/"sample_key_hashes.json",dest/f"{fit}_sample_key_hashes.json"); shutil.copy2(d/"fitted_support.csv",dest/f"{fit}_fitted_support.csv"); shutil.copy2(d/"input_support.csv",dest/f"{fit}_input_support.csv")
        comp.append({"fit":fit,**rows(d/"fit_receipt.csv")[0]}); p={x["coefficient"]:x for x in rows(d/"coefficients.csv")}
        for event,label in EVENTS:
            est,se=(0.,0.) if event=="baseline" else (float(p[event]["estimate"]),float(p[event]["standard_error"]))
            curves.append({"fit":fit,"event":event,"event_label":label,"event_time_years_relative_first_birth":label,"estimate":est,"standard_error":se,"ci_lo":est-1.96*se,"ci_hi":est+1.96*se})
    with (dest/"comparison.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=sorted({k for x in comp for k in x}),lineterminator="\n");w.writeheader();w.writerows(comp)
    with (dest/"full_window.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=list(curves[0]),lineterminator="\n");w.writeheader();w.writerows(curves)
    fig,ax=plt.subplots(figsize=(10,5))
    for fit in FITS:
        d=[x for x in curves if x["fit"]==fit]; x=list(range(len(d))); y=[z["estimate"] for z in d]; lo=[z["ci_lo"] for z in d]; hi=[z["ci_hi"] for z in d]
        z=next(z for z in d if z["event"]=="Wp3")
        ax.plot(x,y,marker="o",label=f"{readable[fit]}: {z['estimate']:.3f}"); ax.fill_between(x,lo,hi,alpha=.15)
    ax.axhline(0,color="black",lw=.7); ax.set_xticks(range(len(EVENTS)),[x[1] for x in EVENTS]); ax.set_xlabel("Years relative to first birth"); ax.set_ylabel("Change in rooms"); ax.set_title("First-birth housing response (95% confidence bands)"); ax.legend(title="+3/+4 estimates"); fig.tight_layout(); fig.savefig(dest/"full_window.png",dpi=180); plt.close(fig)
    verification.update({"status":"pass","sample_key_hashes":h,"completion":"pass"}); write(dest/"verification.json",json.dumps(verification,indent=2)+"\n")

def toy(a: argparse.Namespace) -> None:
    """Run the same staged three-fit machinery on a master intentionally lacking INCFAMR."""
    task=Path(a.task_root).resolve(); task.mkdir(parents=True,exist_ok=True)
    ref=Path(a.reference_estimator).resolve()
    master,income=task/"toy_analysis_sample.dta",task/"income.dta"
    make=task/"make_toy.do"
    write(make, f'''version 17.0
clear all
set more off
set seed 24092026
set obs 1680
gen long ID=ceil(_n/21)
bysort ID: gen int year=1979+2*(_n-1)
gen byte current=1
gen byte AGEREP=25+mod(ID,12)
gen byte EDUYEAR=12+mod(ID,5)
gen double iw=1+mod(ID,4)
gen double bio_first_year=1987+2*mod(ID,7)
replace bio_first_year=. if mod(ID,5)==0
gen double rooms=4+.15*(year>=bio_first_year)+rnormal()
gen double relchirep_max=cond(missing(bio_first_year),0,1)
gen double year_entry_adult=1979
gen byte rel=cond(mod(ID,3)==0,3,cond(mod(ID,2),1,2))
gen double INCFAMR=25000+1000*mod(ID+year,31)
replace INCFAMR=. if mod(ID,17)==0
replace INCFAMR=0 if mod(ID,19)==0
replace INCFAMR=9999999 if mod(ID,23)==0
isid ID year
preserve
 keep ID year INCFAMR
 save "{income}", replace
restore
drop INCFAMR
isid ID year
save "{master}", replace
''')
    subprocess.run([a.stata,"-bq","do",str(make)],cwd=task,check=True)
    actual_ref=sha(ref); actual_master=sha(master)
    stage(argparse.Namespace(task_root=str(task), income_dta=str(income), reference_estimator=str(ref), analysis_sample=str(master), ado_dir=str(task/"ado"), reference_sha256=actual_ref, analysis_sample_sha256=actual_master))
    canonical=task/"canonical_toy.do"
    write(canonical, f'''clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "{task / 'ado'}"
adopath ++ "{task / 'ado'}"
mata: mata mlib index
use "{master}", clear
do "{ref}" A2h "{task / 'results/canonical_toy'}"
exit 0
''')
    subprocess.run([a.stata,"-bq","do",str(canonical)],cwd=task,check=True)
    subprocess.run([a.stata,"-bq","do",str(task/"entry.do")],cwd=task,check=True)
    base={x["coefficient"]:x for x in rows(task/"results/baseline_full/coefficients.csv")}
    canon={x["coefficient"]:x for x in rows(task/"results/canonical_toy/coefficients.csv")}
    if base != canon: raise RuntimeError("toy baseline is not exactly invariant")
    finalize_private(argparse.Namespace(task_root=str(task)))
    if sha(task/"results/common_no_income/sample_key_hashes.json") != sha(task/"results/common_income/sample_key_hashes.json"):
        # File bytes include the same digest but retain this explicit toy gate without using private keys.
        raise RuntimeError("toy common hash receipts differ")
    write(task/"TOY_PASS", "TOY_PASS\n")

def main() -> None:
    p=argparse.ArgumentParser(); q=p.add_subparsers(dest="cmd",required=True)
    x=q.add_parser("extract-income"); x.add_argument("--source-dta",required=True);x.add_argument("--private-dir",required=True);x.add_argument("--stata",default="stata-mp");x.add_argument("--run",action="store_true");x.set_defaults(fn=extract_income)
    s=q.add_parser("stage");s.add_argument("task_root");s.add_argument("--income-dta",required=True);s.add_argument("--reference-estimator",default=str(DEFAULT_REFERENCE));s.add_argument("--analysis-sample");s.add_argument("--ado-dir");s.add_argument("--reference-sha256",default=REFERENCE_SHA);s.add_argument("--analysis-sample-sha256",default=PREPARED_SHA);s.set_defaults(fn=stage)
    v=q.add_parser("validate-stage");v.add_argument("task_root");v.set_defaults(fn=validate_stage)
    f=q.add_parser("finalize-private");f.add_argument("task_root");f.set_defaults(fn=finalize_private)
    t=q.add_parser("toy");t.add_argument("task_root");t.add_argument("--stata",default="stata-mp");t.add_argument("--reference-estimator",default=str(DEFAULT_REFERENCE));t.set_defaults(fn=toy)
    c=q.add_parser("collect-verify");c.add_argument("task_root");c.add_argument("--canonical-dir",required=True);c.add_argument("--expected-sample-key-sha",required=True);c.add_argument("--destination",default=str(DEFAULT_OUTPUT));c.set_defaults(fn=collect_verify)
    a=p.parse_args()
    if a.cmd=="stage":
        t=Path(a.task_root).resolve(); a.analysis_sample=a.analysis_sample or str(t/"analysis_sample.dta"); a.ado_dir=a.ado_dir or str(t/"ado")
    a.fn(a)
if __name__=="__main__": main()
