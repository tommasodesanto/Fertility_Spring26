#!/usr/bin/env bash
# Credit-policy retention replay: three household solves (phi=.8, cap=6, lambda in {0,1,5})
# against the immutable finance_dose_refit_v2 runtime.  Dry-run by default (nothing durable
# is created locally or remotely).  SUBMIT=1 stages remote-to-remote copies of the runtime
# helpers plus the frozen core into a fresh experiment root and submits smoke -> production.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TAG="${TAG:-credit_policy_retention_v1}"; SUBMIT="${SUBMIT:-0}"; SSH_HOST="${SSH_HOST:-torch}"
[[ "$TAG" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "invalid TAG: $TAG" >&2; exit 2; }
REMOTE_ROOT="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/$TAG"
LOCAL_ROOT="$ROOT/output/model/native_financing_diagnostic_20260919/specification_followup/$TAG"
RUNTIME="/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/finance_dose_refit_v2"
FROZEN="/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2"
CHECKPOINT="/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/production/selected_verification/evaluation/raw/repetition_02/initial_state.pkl.gz"
CHECKPOINT_SHA=b3491eedcee6250cf94833067646d3e6463496cbf5a64bdcabe7b13bc7e89eb2
SOLVER_SHA=2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da
PARAMS_SHA=c0c1c18500fba069152659eaf588c3c895993cdcecb104cee7d6edca6bfae6a5
REMOTE_SUMMARY="$RUNTIME/refit_new_income_production/refit_new_income/summary.json"
LOCAL_SUMMARY="$ROOT/output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms/refit_new_income/summary.json"
LOCAL_PLAN="$ROOT/output/model/native_financing_diagnostic_20260919/overnight/finance_dose_refit_v2/plan.remote.json"
WRAPPER_REL=code/model/tools/run_e5f_credit_policy_retention.py
FACTORIAL_REL=code/model/tools/run_e5f_financing_factorial.py
for f in "$ROOT/$WRAPPER_REL" "$ROOT/$FACTORIAL_REL" "$LOCAL_SUMMARY" "$LOCAL_PLAN"; do [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 2; }; done
if [[ "$SUBMIT" == 1 ]]; then
  [[ ! -e "$LOCAL_ROOT" ]] || { echo "local receipt root exists; choose a fresh TAG: $LOCAL_ROOT" >&2; exit 2; }
  STAGE="$LOCAL_ROOT"; mkdir -p "$STAGE"
else
  STAGE="$(mktemp -d "${TMPDIR:-/tmp}/credit_policy_retention_dryrun.XXXXXX")"
fi
MANIFEST="$STAGE/launch_manifest.json"; W="$REMOTE_ROOT/$WRAPPER_REL"; RESULTS="$REMOTE_ROOT/results"

# ---- local syntax review (pure parse; writes no bytecode) --------------------------------------
python3 - "$ROOT/$WRAPPER_REL" "$ROOT/$FACTORIAL_REL" <<'PY'
import ast, sys
for p in sys.argv[1:]: ast.parse(open(p).read(), p)
PY
bash -n "${BASH_SOURCE[0]}"

# ---- pins: every used Python source, summary, plan and the new wrapper -------------------------
python3 - "$ROOT" "$MANIFEST" "$REMOTE_ROOT" "$RUNTIME" "$FROZEN" "$CHECKPOINT" "$CHECKPOINT_SHA" "$SOLVER_SHA" "$PARAMS_SHA" "$REMOTE_SUMMARY" "$LOCAL_SUMMARY" "$LOCAL_PLAN" "$WRAPPER_REL" "$FACTORIAL_REL" "$SSH_HOST" "$STAGE/pins.env" <<'PY'
import hashlib, json, subprocess, sys, time
from pathlib import Path
root, out, exp, runtime, frozen, ckpt, ckpt_sha, solver_sha, params_sha, rsummary, lsummary, lplan, wrapper_rel, factorial_rel, host, pins_env = sys.argv[1:17]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan = json.loads(Path(lplan).read_text())
helpers = {}
for k in ("constructor", "adapter", "controller", "overnight_controller", "factorial_controller"):
    path = plan[k + "_path"]
    if not path.startswith(runtime + "/code/model/tools/"): raise SystemExit(f"plan {k}_path is not under the runtime: {path}")
    helpers[Path(path).name] = plan[k + "_sha256"]
helpers.update(plan["factorial_helper_sha256"])
if len(helpers) != 9: raise SystemExit(f"expected 9 runtime helpers, found {sorted(helpers)}")
if plan["source_root"] != frozen or plan["plan_path"] != runtime + "/plan.json": raise SystemExit("plan source_root/plan_path do not match the runtime contract")
summary = json.loads(Path(lsummary).read_text()); cases = summary["cases"]
if summary.get("status") != "complete" or summary.get("design") != "dose" or len(cases) != 50: raise SystemExit("local summary is not the complete 50-row dose summary")
c0 = cases[0]["contract"]
for key, want in (("checkpoint", ckpt), ("checkpoint_sha256", ckpt_sha), ("source_root", frozen), ("population_source", "saved_evaluation"), ("family", "refit_new_income")):
    if c0[key] != want: raise SystemExit(f"summary contract {key}={c0[key]!r} != {want!r}")
targets = [(0.8, 0.0, 6.0, "case_0_baseline_lambda0"), (0.8, 1.0, 6.0, "case_1_credit_lambda1"), (0.8, 5.0, 6.0, "case_2_credit_lambda5")]
old = {label: [r["label"] for r in cases if (r["phi"], r["lambda"], r["rental_cap"]) == (p, l, c)] for p, l, c, label in targets}
if not all(old.values()): raise SystemExit(f"old receipts missing for some arms: {old}")
try: commit = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
except Exception: commit = "unknown"
manifest = {"schema": "credit_policy_retention_launch_v1", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "local_commit": commit, "ssh_host": host,
            "experiment_root": exp, "results_root": exp + "/results", "runtime_root": runtime, "frozen_root": frozen, "frozen_copy": exp + "/source",
            "checkpoint": {"path": ckpt, "sha256": ckpt_sha, "verification": "sha256 before every case; never copied"},
            "summary": {"remote_path": rsummary, "sha256": sha(lsummary), "local_path": lsummary, "old_labels": old},
            "plan": {"remote_path": runtime + "/plan.json", "copy": exp + "/plan.json", "sha256": sha(lplan), "local_path": lplan},
            "core": {"solver.py": solver_sha, "parameters.py": params_sha}, "runtime_helpers": helpers,
            "wrapper": {"path": exp + "/" + wrapper_rel, "sha256": sha(Path(root) / wrapper_rel), "local_path": str(Path(root) / wrapper_rel)},
            "local_factorial_informational": {"path": str(Path(root) / factorial_rel), "sha256": sha(Path(root) / factorial_rel), "note": "not staged; the runtime copy pinned above is what runs"},
            "design": {"cases": [dict(phi=p, lam=l, cap=c, label=label) for p, l, c, label in targets], "household_solves": 3, "population_source": "saved_evaluation",
                       "policy_identity_gate": False, "reproduction_tolerance": 1e-10, "smoke_minutes": 10, "production_minutes": 20, "case_seconds": 600, "total_seconds": 1800,
                       "threads": 1, "numba_disable_jit": 0, "no_active_core": True, "scope": "fixed-price partial-equilibrium replay; no new target, parameter, entry, source, population, GE or calibration change"}}
Path(out).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
Path(pins_env).write_text(f"SUMMARY_SHA={manifest['summary']['sha256']}\nPLAN_SHA={manifest['plan']['sha256']}\nWRAPPER_SHA={manifest['wrapper']['sha256']}\n"
                          + "HELPER_PINS='" + " ".join(f"{k}:{v}" for k, v in sorted(helpers.items())) + "'\n")
PY
# shellcheck disable=SC1090
source "$STAGE/pins.env"

# ---- remote staging script (read-only pre-copy hash checks, then fresh root + copies) -----------
cat >"$STAGE/remote_stage.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
E="$REMOTE_ROOT"; R="$RUNTIME"; F="$FROZEN"
check() { local got; got="\$(sha256sum "\$1" | cut -d' ' -f1)"; [[ "\$got" == "\$2" ]] || { echo "hash mismatch: \$1 got \$got expected \$2" >&2; exit 4; }; }
check "$CHECKPOINT" "$CHECKPOINT_SHA"
check "\$F/code/model/intergen_eqscale_seq_optimized/solver.py" "$SOLVER_SHA"
check "\$F/code/model/intergen_eqscale_seq_optimized/parameters.py" "$PARAMS_SHA"
check "$REMOTE_SUMMARY" "$SUMMARY_SHA"
check "\$R/plan.json" "$PLAN_SHA"
for pin in $HELPER_PINS; do check "\$R/code/model/tools/\${pin%%:*}" "\${pin##*:}"; done
[[ "\$(ls "\$R/code/model/tools/"*.py | wc -l)" == 9 ]] || { echo "runtime tools dir does not hold exactly the 9 pinned helpers; refusing to stage" >&2; exit 4; }
[[ -f "\$F/code/model/tools/run_e5f_matched_pf_smoke.py" && -f "\$F/code/model/tools/run_e5f_independent_numerical_audit.py" ]] || { echo "frozen tools missing matched_pf_smoke/independent_numerical_audit" >&2; exit 4; }
[[ ! -e "\$E" ]] || { echo "experiment root exists: \$E" >&2; exit 4; }
mkdir -p "\$E/code/model/tools" "\$E/source/code/model" "\$E/results" "\$E/logs" "\$E/numba_cache"
rsync -a --include='*.py' --exclude='*' "\$R/code/model/tools/" "\$E/code/model/tools/"
rsync -a --exclude='__pycache__' --exclude='*.pyc' --exclude='.pytest_cache' --exclude='numba_cache' --exclude='.numba_cache' --exclude='output' --exclude='outputs' --exclude='tmp' --exclude='*.log' --exclude='*.out' --exclude='*.png' --exclude='*.pdf' "\$F/code/model/" "\$E/source/code/model/"
cp -p "\$R/plan.json" "\$E/plan.json"
python3 - "\$F/code/model" "\$E/source/code/model" "\$E/source_manifest.json" <<'PY'
import hashlib,json,pathlib,sys
original,copy,out=map(pathlib.Path,sys.argv[1:])
pins={str(p.relative_to(original)):hashlib.sha256(p.read_bytes()).hexdigest() for p in original.rglob('*.py') if '__pycache__' not in p.parts}
for rel,pin in pins.items():
 assert hashlib.sha256((copy/rel).read_bytes()).hexdigest()==pin,rel
out.write_text(json.dumps(pins,indent=2,sort_keys=True)+'\n')
PY
[[ "\$(find "\$E/source" -name __pycache__ | wc -l)" == 0 ]] || { echo "cache directories copied" >&2; exit 4; }
echo "staged \$E"
EOF

# ---- remote finalize: verify copied files against pins, then make the staged inputs read-only --
cat >"$STAGE/remote_finalize.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
E="$REMOTE_ROOT"
check() { local got; got="\$(sha256sum "\$1" | cut -d' ' -f1)"; [[ "\$got" == "\$2" ]] || { echo "hash mismatch: \$1 got \$got expected \$2" >&2; exit 4; }; }
for pin in $HELPER_PINS; do check "\$E/code/model/tools/\${pin%%:*}" "\${pin##*:}"; done
check "\$E/code/model/tools/run_e5f_credit_policy_retention.py" "$WRAPPER_SHA"
check "\$E/source/code/model/intergen_eqscale_seq_optimized/solver.py" "$SOLVER_SHA"
check "\$E/source/code/model/intergen_eqscale_seq_optimized/parameters.py" "$PARAMS_SHA"
check "\$E/plan.json" "$PLAN_SHA"
[[ "\$(ls "\$E/code/model/tools/"*.py | wc -l)" == 10 ]] || { echo "unexpected helper count in \$E/code/model/tools" >&2; exit 4; }
chmod -R a-w "\$E/code" "\$E/source" "\$E/plan.json" "\$E/launch_manifest.json" "\$E/smoke.sbatch" "\$E/production.sbatch"
echo "finalized \$E"
EOF

# ---- sbatch files: fresh root, explicit logs/chdir, fail-closed hashes before each job ----------
sbatch_header() {
  cat <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cpr_$1
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=$2
#SBATCH --chdir=$REMOTE_ROOT
#SBATCH --output=$REMOTE_ROOT/logs/%x_%j.out
#SBATCH --error=$REMOTE_ROOT/logs/%x_%j.err
#SBATCH --kill-on-invalid-dep=yes
set -euo pipefail
module load anaconda3/2025.06
unset PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR="$REMOTE_ROOT/numba_cache" MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd "$REMOTE_ROOT"
W="$W"; M="$REMOTE_ROOT/launch_manifest.json"; RES="$RESULTS"
python -B "\$W" --mode verify --manifest "\$M" --results "\$RES" --receipt "\$RES/verify_$1_\${SLURM_JOB_ID}.json"
EOF
}
{ sbatch_header smoke 00:10:00; cat <<EOF
DEADLINE=\$(( \$(date +%s) + 1800 ))
timeout --signal=TERM --kill-after=30s 600s python -B "\$W" --mode case --case 0 --manifest "\$M" --results "\$RES" --deadline-epoch "\$DEADLINE"
EOF
} >"$STAGE/smoke.sbatch"
{ sbatch_header production 00:20:00; cat <<EOF
python -B "\$W" --mode check-smoke --manifest "\$M" --results "\$RES"
DEADLINE=\$(( \$(date +%s) + 1800 ))
for i in 1 2; do
  timeout --signal=TERM --kill-after=30s 600s python -B "\$W" --mode case --case "\$i" --manifest "\$M" --results "\$RES" --deadline-epoch "\$DEADLINE"
done
python -B "\$W" --mode compare --manifest "\$M" --results "\$RES"
EOF
} >"$STAGE/production.sbatch"
for f in "$STAGE/remote_stage.sh" "$STAGE/remote_finalize.sh" "$STAGE/smoke.sbatch" "$STAGE/production.sbatch"; do bash -n "$f"; done
python3 - "$MANIFEST" "$STAGE" <<'PY'
import hashlib, json, sys
from pathlib import Path
m = json.loads(Path(sys.argv[1]).read_text()); stage = Path(sys.argv[2])
m["generated_files"] = {n: hashlib.sha256((stage / n).read_bytes()).hexdigest() for n in ("remote_stage.sh", "remote_finalize.sh", "smoke.sbatch", "production.sbatch")}
Path(sys.argv[1]).write_text(json.dumps(m, indent=2, sort_keys=True) + "\n")
PY
if [[ "$SUBMIT" != 1 ]]; then
  echo "dry-run: syntax review passed; nothing staged or submitted."
  echo "review: $MANIFEST $STAGE/remote_stage.sh $STAGE/remote_finalize.sh $STAGE/smoke.sbatch $STAGE/production.sbatch"
  echo "remote root (untouched): $REMOTE_ROOT"; exit 0
fi

# ---- SUBMIT=1: stage remote-to-remote, upload wrapper/manifest/sbatch, verify, submit ------------
SSH=(ssh -o BatchMode=yes "$SSH_HOST")
"${SSH[@]}" bash -s <"$STAGE/remote_stage.sh"
scp -q "$SSH_HOST:$REMOTE_ROOT/source_manifest.json" "$STAGE/source_manifest.json"
python3 - "$MANIFEST" "$STAGE/source_manifest.json" "$REMOTE_ROOT/source_manifest.json" <<'PY'
import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]);m=json.loads(p.read_text())
m['source_manifest']={'path':sys.argv[3],'sha256':hashlib.sha256(pathlib.Path(sys.argv[2]).read_bytes()).hexdigest()}
p.write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
PY
scp -q "$ROOT/$WRAPPER_REL" "$SSH_HOST:$REMOTE_ROOT/code/model/tools/run_e5f_credit_policy_retention.py"
scp -q "$MANIFEST" "$STAGE/smoke.sbatch" "$STAGE/production.sbatch" "$STAGE/remote_stage.sh" "$STAGE/remote_finalize.sh" "$SSH_HOST:$REMOTE_ROOT/"
"${SSH[@]}" bash -s <"$STAGE/remote_finalize.sh"
SMOKE="$("${SSH[@]}" "sbatch --parsable '$REMOTE_ROOT/smoke.sbatch'")"
[[ "$SMOKE" =~ ^[0-9]+$ ]] || { echo "smoke submission failed: $SMOKE" >&2; exit 5; }
python3 - "$STAGE/submission.json" "$SMOKE" "$REMOTE_ROOT" "$MANIFEST" <<'PY'
import json, sys, time
json.dump({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "remote_root": sys.argv[3], "manifest": sys.argv[4], "jobs": {"smoke": sys.argv[2], "production": None},
           "status": "smoke submitted; production pending afterok"}, open(sys.argv[1], "w"), indent=2); print(file=open(sys.argv[1], "a"))
PY
PROD="$("${SSH[@]}" "sbatch --parsable --dependency=afterok:$SMOKE '$REMOTE_ROOT/production.sbatch'")"
[[ "$PROD" =~ ^[0-9]+$ ]] || { echo "production submission failed after smoke=$SMOKE: $PROD" >&2; exit 5; }
python3 - "$STAGE/submission.json" "$PROD" <<'PY'
import json, sys
p = sys.argv[1]; d = json.load(open(p)); d["jobs"]["production"] = sys.argv[2]; d["status"] = "submitted; production afterok smoke"
json.dump(d, open(p, "w"), indent=2); print(file=open(p, "a"))
PY
echo "submitted smoke=$SMOKE production=$PROD root=$REMOTE_ROOT receipts=$STAGE"
