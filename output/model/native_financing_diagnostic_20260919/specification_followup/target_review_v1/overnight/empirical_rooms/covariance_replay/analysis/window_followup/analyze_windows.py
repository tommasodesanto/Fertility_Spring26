"""Mechanical event-window contrasts from the frozen reviewed replay; no refit."""
from __future__ import annotations
import csv, json, math, hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPLAY = HERE.parents[1]
ANALYSIS = REPLAY / "analysis"
OUT = REPLAY / "output"

def read_csv(p):
    with p.open(newline="") as f: return list(csv.DictReader(f))
def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()
def event(k): return f"F{-k}event" if k < 0 else f"L{k}event"

# Reconfirm reviewed terminal gate and replay analysis; never invoke Stata or raw data.
execution=json.loads((REPLAY/"execution.json").read_text())
old=json.loads((ANALYSIS/"replay_analysis.json").read_text())
if execution.get("status") != "process_completed_pending_scientific_review" or execution.get("exit_code") != 0:
    raise SystemExit("STOP: source replay is not terminal-successful")
if old.get("status", "").find("analyzed completed exact replay") < 0:
    raise SystemExit("STOP: reviewed replay analysis status not established")
coefs=read_csv(OUT/"cohort_interaction_coefficients.csv")
q=len(coefs)
if q != 900 or old.get("cohort_interaction_matrix_dimensions") != [900,900]:
    raise SystemExit(f"STOP: full interaction covariance dimension mismatch ({q})")
byidx={int(x["matrix_index"]):x for x in coefs}
idx={(int(float(x["cohort"])),x["event"]):int(x["matrix_index"]) for x in coefs}
if sorted(byidx) != list(range(1,q+1)): raise SystemExit("STOP: coefficient index alignment invalid")
# Read all 810,000 cells of full e(V), preserving exported matrix indices.
V=[[0.0]*q for _ in range(q)]
rows=read_csv(OUT/"cohort_interaction_covariance.csv")
if len(rows)!=q*q: raise SystemExit(f"STOP: covariance has {len(rows)} rows, expected {q*q}")
for r in rows:
    i,j=int(r["row_index"])-1,int(r["column_index"])-1
    V[i][j]=float(r["covariance"])
for i in range(q):
    if abs(V[i][i]-float(byidx[i+1]["marginal_variance"]))>1e-14: raise SystemExit("STOP: full e(V) diagonal alignment failed")
    for j in range(i+1,q):
        if abs(V[i][j]-V[j][i])>1e-10: raise SystemExit("STOP: full e(V) symmetry failed")
# only estimation-sample cohort/event support; ignore never-treated rows.
support={}
for r in read_csv(OUT/"estimation_cohort_event_support.csv"):
    g=r["cohort_group"].strip(); k=r["event_time"].strip()
    if not g or g=="never_treated" or not k: continue
    support[(int(float(g)),int(float(k)))]=(int(r["estimation_observations"]),float(r["longitudinal_weight"]))

def cohorts_at(k): return {g for (g,t),(n,w) in support.items() if t==k and n>0}
def calc(name,start_k,end_k,weight_k,start_coeff=True, require_ref_support=False, omitted_end=False, restrict_cohorts=None):
    cs=sorted(cohorts_at(start_k)&cohorts_at(end_k))
    if restrict_cohorts is not None: cs=[g for g in cs if g in set(restrict_cohorts)]
    if require_ref_support: cs=[g for g in cs if g in cohorts_at(-2)]
    if not cs: raise ValueError(f"no common support: {name}")
    raw={g:support[g,weight_k][1] for g in cs}
    if any(not math.isfinite(w) or w<=0 for w in raw.values()): raise ValueError("invalid weights")
    total=math.fsum(raw.values()); weights={g:raw[g]/total for g in cs}
    c=[0.0]*q; bdiff={}
    for g in cs:
        end_i=None if omitted_end else idx.get((g,event(end_k)))
        if not omitted_end and end_i is None: raise ValueError(f"missing endpoint coefficient cohort={g}, event={event(end_k)}")
        if omitted_end and g not in cohorts_at(end_k): raise ValueError(f"unsupported omitted endpoint reference for {g}")
        endb=0.0 if omitted_end else float(byidx[end_i]["estimate"])
        if start_coeff:
            start_i=idx.get((g,event(start_k)))
            if start_i is None: raise ValueError(f"missing start coefficient cohort={g}, event={event(start_k)}")
            startb=float(byidx[start_i]["estimate"])
            c[start_i-1]-=weights[g]
        else:
            # The omitted -2 reference is zero only on supported cohorts.
            if g not in cohorts_at(-2): raise ValueError(f"unsupported omitted -2 reference for {g}")
            startb=0.0
        if end_i is not None: c[end_i-1]+=weights[g]
        bdiff[g]=endb-startb
    estimate=math.fsum(weights[g]*bdiff[g] for g in cs)
    variance=math.fsum(c[i]*V[i][j]*c[j] for i in range(q) for j in range(q))
    if variance < -1e-12: raise ValueError(f"negative variance for {name}: {variance}")
    # Algebraic normalization check: arbitrary cohort-specific shifts added to both endpoints cancel.
    shifts={g:math.sin(g*0.37)+g/10000 for g in cs}
    shifted_est=math.fsum(weights[g]*(bdiff[g]+shifts[g]-shifts[g]) for g in cs)
    cohort_loading_sums={str(g):((-weights[g] if start_coeff else 0.0)
       +(weights[g] if not omitted_end else 0.0)
       +(weights[g] if omitted_end else 0.0)) for g in cs}
    # The final +w for an omitted endpoint is the shifted omitted-reference loading.
    if not omitted_end and not start_coeff:
        cohort_loading_sums={str(g):weights[g]-weights[g] for g in cs}
    shifted_gap=abs(shifted_est-estimate)
    loading=math.fsum(c)
    return {"name":name,"start_event_time":start_k,"end_event_time":end_k,
      "start_coefficient":event(start_k) if start_coeff else "omitted F2event = 0, supported cohorts only",
      "end_coefficient":"omitted F2event = 0, supported cohorts only" if omitted_end else event(end_k),"weight_endpoint":weight_k,
      "cohorts":cs,"cohort_count":len(cs),
      "start_endpoint_supported_cohort_count":len(cohorts_at(start_k)),
      "end_endpoint_supported_cohort_count":len(cohorts_at(end_k)),
      "retained_start_weight_mass":math.fsum(support[g,start_k][1] for g in cs)/math.fsum(support[g,start_k][1] for g in cohorts_at(start_k)),
      "retained_end_weight_mass":math.fsum(support[g,end_k][1] for g in cs)/math.fsum(support[g,end_k][1] for g in cohorts_at(end_k)),
      "estimate":estimate,"variance":max(0.0,variance),"standard_error":math.sqrt(max(0.0,variance)),
      "weight_shares":{str(g):weights[g] for g in cs},
      "normalization_invariance":{"cohort_specific_endpoint_shift_test":"add a distinct arbitrary constant to both endpoints within each cohort; omitted reference is shifted algebraically when applicable, while its empirical support requirement remains in force","max_abs_estimate_gap_after_cohort_specific_endpoint_shifts":shifted_gap,"per_cohort_loading_sums":cohort_loading_sums,"max_abs_per_cohort_loading_sum":max(abs(v) for v in cohort_loading_sums.values()),"aggregate_coefficient_loading_sum":loading,"weight_scale_check_factor":13.7},
      "se_interpretation":"conditional on these fixed cohort weights; excludes sampling variation in estimated cohort shares"}

# Reproduce all four previously reviewed contrasts before adding requested windows.
checks=[]
for name, spec in [
 ("common_cohorts_fixed_prebirth_IW",(-1,3,-1,True)),
 ("common_cohorts_fixed_postbirth_IW",(-1,3,3,True)),
 ("common_cohorts_minus2_to_plus2_fixed_minus2_IW",(-2,2,-2,False)),
 ("common_cohorts_minus2_to_plus2_fixed_plus2_IW",(-2,2,2,False))]:
    a=next(x for x in old["contrasts"] if x["name"]==name)
    specarg=spec
    # Historical -2 +2 uses both endpoints' support, and -2 itself is observed for included cohorts.
    got=calc(name,*specarg,require_ref_support=(not specarg[3]))
    gap=abs(got["estimate"]-a["estimate"])
    if gap>1e-12: raise SystemExit(f"STOP: old contrast {name} mismatch {gap:.3g}")
    variance_gap=abs(got["variance"]-a["variance"])
    se_gap=abs(got["standard_error"]-a["standard_error"])
    if variance_gap>1e-12 or se_gap>1e-12: raise SystemExit(f"STOP: old variance/SE mismatch {name}: {variance_gap:.3g}/{se_gap:.3g}")
    got.update({"validation_reference_estimate":a["estimate"],"validation_abs_gap":gap,"validation_reference_variance":a["variance"],"validation_variance_abs_gap":variance_gap,"validation_reference_se":a["standard_error"],"validation_se_abs_gap":se_gap})
    checks.append(got)

# Requested -4 to +4 using fixed pre (-4) shares.
minus4_plus4=calc("minus4_to_plus4_fixed_minus4_weights",-4,4,-4,True)
minus4_plus4_post=calc("minus4_to_plus4_fixed_plus4_weights",-4,4,4,True)
# Requested prebirth contrast on exactly that support cannot assign omitted -2 zero:
# cohorts 1984 and 1985 have no e(sample) observations at -2. Do not silently include them.
shared=minus4_plus4["cohorts"]
unsupported=[g for g in shared if g not in cohorts_at(-2)]
same_support_prebirth={"name":"minus4_to_minus2_same_minus4_plus4_support_fixed_minus4_weights",
 "status":"not_calculated_source_support_failure","required_cohorts":shared,"unsupported_omitted_minus2_reference_cohorts":unsupported,
 "reason":"The omitted -2 coefficient may be set to zero only for cohorts with positive e(sample) support at -2; exact -4/+4 support includes cohorts lacking that support."}
# Separate own common support, still fixed at -4 endpoint weights.
minus4_minus2=calc("minus4_to_minus2_own_common_support_fixed_minus4_weights",-4,-2,-4,True,omitted_end=True)
# Three-way shared support makes both comparisons use identical cohorts and fixed -4 weights.
threeway=sorted(cohorts_at(-4)&cohorts_at(4)&cohorts_at(-2))
threeway_post=calc("threeway_common_support_minus4_to_plus4_fixed_minus4_weights",-4,4,-4,True,restrict_cohorts=threeway)
threeway_pre=calc("threeway_common_support_minus4_to_minus2_fixed_minus4_weights",-4,-2,-4,True,omitted_end=True,restrict_cohorts=threeway)
# Include retained-weight summaries and cohort lists for all contrasts in compact table.
results={"status":"completed_with_one_requested_same_support_contrast_unavailable",
 "scope":"Mechanical postprocessing of the single reviewed replay; no refit, rawdata, or model solve.",
 "source_covariance":"full 900x900 cohort-interaction block from e(V); not e(V_interact).",
 "covariance_cells":len(rows),"max_symmetry_gap":max(abs(V[i][j]-V[j][i]) for i in range(q) for j in range(q)),
 "old_contrast_validation_tolerance":1e-12,"old_contrast_validation":checks,
 "requested_contrasts":{"minus4_to_plus4_fixed_minus4_weights":minus4_plus4,"minus4_to_plus4_fixed_plus4_weights":minus4_plus4_post,"minus4_to_minus2_same_minus4_plus4_support":same_support_prebirth,"minus4_to_minus2_own_common_support":minus4_minus2,"threeway_common_support_23_cohorts":{"cohorts":threeway,"minus4_to_plus4":threeway_post,"minus4_to_minus2":threeway_pre}},
 "limitations":["All standard errors are conditional on fixed cohort weights.","No causal, preferred-window, or model-timing interpretation is made."]}
(HERE/"window_results.json").write_text(json.dumps(results,indent=2)+"\n")
rowsout=[]
for x in checks+[minus4_plus4,minus4_plus4_post,minus4_minus2,threeway_post,threeway_pre]:
    rowsout.append({"name":x["name"],"start_event_time":x["start_event_time"],"end_event_time":x["end_event_time"],"weight_endpoint":x["weight_endpoint"],"cohort_count":x["cohort_count"],"cohorts":";".join(map(str,x["cohorts"])),"retained_start_weight_mass":x["retained_start_weight_mass"],"retained_end_weight_mass":x["retained_end_weight_mass"],"estimate":x["estimate"],"variance":x["variance"],"standard_error":x["standard_error"],"max_abs_cohort_loading_sum":x["normalization_invariance"]["max_abs_per_cohort_loading_sum"],"validation_abs_gap":x.get("validation_abs_gap","")})
with (HERE/"window_results.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rowsout[0]),lineterminator="\n");w.writeheader();w.writerows(rowsout)
inputs=[OUT/"cohort_interaction_coefficients.csv",OUT/"cohort_interaction_covariance.csv",OUT/"estimation_cohort_event_support.csv",ANALYSIS/"replay_analysis.json",REPLAY/"execution.json"]
receipt={"status":results["status"],"script_sha256":sha256(Path(__file__)),"input_sha256":{p.name:sha256(p) for p in inputs},"summary":{"old_validation_max_abs_estimate_gap":max(x["validation_abs_gap"] for x in checks),"old_validation_max_abs_variance_gap":max(x["validation_variance_abs_gap"] for x in checks),"old_validation_max_abs_se_gap":max(x["validation_se_abs_gap"] for x in checks),"minus4_plus4_preweighted":{"estimate":minus4_plus4["estimate"],"se":minus4_plus4["standard_error"]},"minus4_plus4_postweighted":{"estimate":minus4_plus4_post["estimate"],"se":minus4_plus4_post["standard_error"]},"threeway_minus4_plus4":{"estimate":threeway_post["estimate"],"se":threeway_post["standard_error"]},"threeway_minus4_minus2":{"estimate":threeway_pre["estimate"],"se":threeway_pre["standard_error"]},"same_support_minus2_unavailable_cohorts":unsupported}}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
