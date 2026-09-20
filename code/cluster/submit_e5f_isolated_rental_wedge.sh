#!/usr/bin/env bash
# Isolated rental-size wedge: six fixed-price household solves on the original checkpoint
# (cap6 slope0, cap10 slope0, cap10 slopes .05/.2/1, plus cap10 slope .2 at phi=1) from an immutable experiment root that
# combines the reviewed rental-wedge source snapshot (frozen Sep-14 core + explicit isolated
# patch, uploaded from tmp/e5f_rental_wedge_exhaustive_v1) with the nine finance_dose_v1 runtime
# helpers (copied remote-to-remote).  Dry-run by default: SUBMIT=0 parses, pins, writes the
# launch manifest, the driver plan, and the sbatch scripts into a temporary directory and
# submits nothing.  SUBMIT=1 requires source_manifest.rental_wedge_port_reviewed=true, stages a
# fresh remote root, verifies every hash, makes the staged code read-only, and submits
# smoke (cap6zero + cap10s02) -> production (afterok: cap10zero + cap10s005 + cap10s1 + cap10s02_phi1 + combine).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TAG="${TAG:-rental_wedge_v1}"; SUBMIT="${SUBMIT:-0}"; SSH_HOST="${SSH_HOST:-torch}"
[[ "$TAG" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "invalid TAG: $TAG" >&2; exit 2; }
REMOTE_ROOT="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/$TAG"
LOCAL_ROOT="$ROOT/output/model/native_financing_diagnostic_20260919/specification_followup/$TAG"
RUNTIME="/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/finance_dose_v1"
REMOTE_SUMMARY="$RUNTIME/original_production/original/summary.json"
LOCAL_SUMMARY="$ROOT/output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms/original/summary.json"
LOCAL_PLAN="$ROOT/output/model/native_financing_diagnostic_20260919/overnight/finance_dose_v1/plan.remote.json"
REFIT_PLAN="$ROOT/output/model/native_financing_diagnostic_20260919/overnight/finance_dose_refit_v2/plan.remote.json"
SNAPSHOT="$ROOT/tmp/e5f_rental_wedge_exhaustive_v1"
PORT_MANIFEST="$SNAPSHOT/source_manifest.json"
LOCAL_CHECKPOINT="$ROOT/output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_02/initial_state.pkl.gz"
CHECKPOINT_SHA=3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993
DRIVER_REL=code/model/tools/run_e5f_isolated_rental_wedge.py
FACTORIAL_REL=code/model/tools/run_e5f_financing_factorial.py
for f in "$ROOT/$DRIVER_REL" "$LOCAL_SUMMARY" "$LOCAL_PLAN" "$PORT_MANIFEST"; do [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 2; }; done
[[ -d "$SNAPSHOT/code/model" ]] || { echo "missing ported snapshot: $SNAPSHOT/code/model" >&2; exit 2; }
if [[ "$SUBMIT" == 1 ]]; then
  [[ ! -e "$LOCAL_ROOT" ]] || { echo "local receipt root exists; choose a fresh TAG: $LOCAL_ROOT" >&2; exit 2; }
  STAGE="$LOCAL_ROOT"; mkdir -p "$STAGE"
else
  STAGE="$(mktemp -d "${TMPDIR:-/tmp}/rental_wedge_dryrun.XXXXXX")"
fi
MANIFEST="$STAGE/launch_manifest.json"; W="$REMOTE_ROOT/$DRIVER_REL"; RESULTS="$REMOTE_ROOT/results"

# ---- local syntax review (pure parse; writes no bytecode) --------------------------------------
python3 - "$ROOT/$DRIVER_REL" <<'PY'
import ast, sys
for p in sys.argv[1:]: ast.parse(open(p).read(), p)
PY
bash -n "${BASH_SOURCE[0]}"

# ---- driver plan: verifies every snapshot hash locally; never requires the review marker -------
PLAN_ARGS=(--mode plan --source-root "$SNAPSHOT" --manifest "$PORT_MANIFEST" --output "$STAGE/driver_plan.json")
if [[ -f "$LOCAL_CHECKPOINT" ]]; then PLAN_ARGS+=(--checkpoint "$LOCAL_CHECKPOINT"); else echo "note: local checkpoint copy absent; the remote stage verifies the retained checkpoint hash"; fi
python3 -B "$ROOT/$DRIVER_REL" "${PLAN_ARGS[@]}"

# ---- pins: checkpoint/frozen root from the old summary contract, nine runtime helpers, port ----
python3 -B - "$ROOT" "$MANIFEST" "$REMOTE_ROOT" "$RUNTIME" "$CHECKPOINT_SHA" "$REMOTE_SUMMARY" "$LOCAL_SUMMARY" "$LOCAL_PLAN" "$REFIT_PLAN" "$SNAPSHOT" "$PORT_MANIFEST" "$DRIVER_REL" "$FACTORIAL_REL" "$SSH_HOST" "$STAGE" "$SUBMIT" "$LOCAL_CHECKPOINT" <<'PY'
import hashlib, json, subprocess, sys, time
from pathlib import Path
(root, out, exp, runtime, ckpt_sha, rsummary, lsummary, lplan, refit_plan, snapshot, port_manifest, driver_rel, factorial_rel, host, stage, submit, local_ckpt) = sys.argv[1:18]
sys.path.insert(0, str(Path(root) / "code/model/tools"))
import run_e5f_isolated_rental_wedge as drv
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
port = json.loads(Path(port_manifest).read_text())
source = drv.load_source_manifest(Path(snapshot), Path(port_manifest))
if port.get("schema") != "e5f_isolated_rental_wedge_source_v2": raise SystemExit(f"unexpected port manifest schema: {port.get('schema')}")
if port.get("checkpoint_sha256") != ckpt_sha: raise SystemExit("port manifest checkpoint pin differs")
reviewed = bool(port.get("rental_wedge_port_reviewed", False))
if submit == "1" and not reviewed: raise SystemExit("SUBMIT=1 refused: source_manifest.rental_wedge_port_reviewed is false (lead review marker required; dry-run remains available)")
relevant = list(port["relevant_files"]); base = port["base_hashes_frozen"]; ported = port["ported_hashes"]
if sorted(ported) != sorted(relevant) or not (set(relevant) - set(port.get("owned_tests", []))) <= set(base): raise SystemExit("port manifest base/ported maps do not cover changed production files")
changed = drv.changed_source_hashes(port)
plan = json.loads(Path(lplan).read_text()); helpers = {}
for k in ("constructor", "adapter", "controller", "overnight_controller", "factorial_controller"):
    path = plan[k + "_path"]
    if not path.startswith(runtime + "/code/model/tools/"): raise SystemExit(f"plan {k}_path is not under the runtime: {path}")
    helpers[Path(path).name] = plan[k + "_sha256"]
helpers.update(plan["factorial_helper_sha256"])
if sorted(helpers) != sorted(n + ".py" for n in drv.RUNTIME_HELPERS): raise SystemExit(f"runtime helper set {sorted(helpers)} differs from the driver's pinned nine")
summary = json.loads(Path(lsummary).read_text()); cases = summary["cases"]
if summary.get("status") != "complete" or summary.get("design") != "dose" or len(cases) != 50: raise SystemExit("local summary is not the complete 50-row original dose summary")
c0 = cases[0]["contract"]; frozen = c0["source_root"]; ckpt = c0["checkpoint"]
if c0["checkpoint_sha256"] != ckpt_sha or c0["family"] != "original": raise SystemExit("summary contract is not the original-family retained checkpoint")
if not ckpt.startswith("/scratch/") or not frozen.startswith("/scratch/"): raise SystemExit("summary contract paths are not shared scratch paths")
if plan["source_root"] != frozen or plan["plan_path"] != runtime + "/plan.json" or plan["source_manifest_path"] != c0["source_manifest"]: raise SystemExit("plan source_root/plan_path/source_manifest_path do not match the summary contract")
for r in cases:
    c = r["contract"]
    if (c["checkpoint"], c["checkpoint_sha256"], c["source_root"], c["family"], c["price"], c["initial_population_shape"], c["source_manifest"], c.get("candidate_fingerprint")) != (ckpt, ckpt_sha, frozen, "original", c0["price"], c0["initial_population_shape"], c0["source_manifest"], None):
        raise SystemExit(f"old summary rows carry mixed contracts: {r['label']}")
    if r["population"] != drv.OLD_POPULATION_LABEL or r["preferences"] != drv.OLD_PREFERENCES_LABEL: raise SystemExit(f"old row labels differ: {r['label']}")
old = {f"cap{cap:g}": [r["label"] for r in drv.old_rows(summary, cap)] for cap in (6.0, 10.0)}
refit = None
if Path(refit_plan).is_file():
    rp = json.loads(Path(refit_plan).read_text())
    rh = {Path(rp[k + "_path"]).name: rp[k + "_sha256"] for k in ("constructor", "adapter", "controller", "overnight_controller", "factorial_controller")}; rh.update(rp["factorial_helper_sha256"])
    refit = {"plan": refit_plan, "identical_to_finance_dose_v1": rh == helpers, "differing": sorted(k for k in helpers if rh.get(k) != helpers[k]), "note": "informational; finance_dose_v1 produced the original summary and is the pinned runtime"}
try: commit = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
except Exception: commit = "unknown"
local_ckpt_row = {"path": local_ckpt, "sha256": sha(local_ckpt), "matches": sha(local_ckpt) == ckpt_sha} if Path(local_ckpt).is_file() else {"path": local_ckpt, "present": False}
manifest = {"schema": "e5f_isolated_rental_wedge_launch_v1", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "local_commit": commit, "ssh_host": host,
            "experiment_root": exp, "results_root": exp + "/results", "runtime_root": runtime, "frozen_root": frozen, "source_copy": exp + "/source",
            "checkpoint": {"path": ckpt, "sha256": ckpt_sha, "verification": "sha256 before staging and before every case; never copied", "local_copy_informational": local_ckpt_row},
            "summary": {"remote_path": rsummary, "sha256": sha(lsummary), "local_path": lsummary, "old_labels": old, "population_label": drv.OLD_POPULATION_LABEL},
            "plan": {"remote_path": runtime + "/plan.json", "copy": exp + "/plan.json", "sha256": sha(lplan), "local_path": lplan},
            "runtime_helpers": helpers, "runtime_helpers_used_by_driver": list(drv.RUNTIME_HELPERS_USED), "refit_v2_helpers_informational": refit,
            "driver": {"path": exp + "/" + driver_rel, "sha256": sha(Path(root) / driver_rel), "local_path": str(Path(root) / driver_rel)},
            "port_manifest": {"local_path": port_manifest, "copy": exp + "/source_manifest.json", "sha256": sha(port_manifest), "snapshot": snapshot, "source_files": source["source_files"],
                              "rental_wedge_port_reviewed": reviewed, "frozen_source_root_local": port["frozen_source_root"], "diff_path_local": port.get("diff_path"),
                              "remaining_review_blockers": port.get("remaining_review_blockers"), "port_scope": port.get("port_scope")},
            "changed_source_vs_frozen": changed,
            "local_factorial_informational": {"path": str(Path(root) / factorial_rel), "sha256": sha(Path(root) / factorial_rel), "note": "not staged and not imported; the runtime copy pinned above is staged for provenance only"},
            "design": {"cases": [c.__dict__ for c in drv.CASES], "stages": {k: list(v) for k, v in drv.STAGE_CASES.items()}, "household_solves": drv.HOUSEHOLD_SOLVES_PLANNED, "max_lifecycle_solves": drv.MAX_LIFECYCLE_SOLVES,
                       "case_timeout_seconds": dict(drv.CASE_TIMEOUT_SECONDS), "stage_seconds": dict(drv.STAGE_SECONDS), "slurm_minutes": {"smoke": 45, "production": 60}, "cpus": 1, "mem_gb": 24,
                       "reproduction_tolerance": drv.TOL, "saving_audit": {"draws": drv.SAVING_DRAWS, "gain_tolerance": drv.SAVING_GAIN_TOLERANCE, "retained_definition": drv.RETAINED_SAVING_AUDIT_DEFINITION},
                       "cost_function": "C(h) = rent*h + slope*h*max(h-6,0); intercept 0; knee 6", "threads": 1, "numba_disable_jit": 0, "no_active_core": True, "auto_retries": False,
                       "identity_gates": {"population_bitwise": True, "entry_cohort_bitwise": True, "policy_identity": False},
            "scope": "fixed-price partial equilibrium on the original checkpoint; five fixed-phi=.8 arms plus one matched cap10/slope=.2 financed-share contrast at phi=1; chi, preferences, prices, lambda=0, raw stationary_g_pre, entry and fiscal contracts otherwise unchanged; no GE or stationary calibration; positive slopes are findings"}}
Path(out).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
frozen_pins = " ".join(f"{rel}:{base[rel]}" for rel in sorted(base))
Path(stage, "pins.env").write_text(f"CHECKPOINT={ckpt}\nFROZEN={frozen}\nSUMMARY_SHA={manifest['summary']['sha256']}\nPLAN_SHA={manifest['plan']['sha256']}\nDRIVER_SHA={manifest['driver']['sha256']}\nPORT_SHA={manifest['port_manifest']['sha256']}\nREVIEWED={'1' if reviewed else '0'}\n"
                                   + "HELPER_PINS='" + " ".join(f"{k}:{v}" for k, v in sorted(helpers.items())) + "'\n" + f"FROZEN_PINS='{frozen_pins}'\n")
cli = drv.planned_cli(exp)
for stage_name in ("smoke", "production"):
    Path(stage, f"{stage_name}_body.sh").write_text("\n".join(cli[stage_name]) + "\n")
PY
# shellcheck disable=SC1090
source "$STAGE/pins.env"

# ---- remote staging: read-only pre-copy hash checks, fresh root, helpers + frozen non-Python ---
cat >"$STAGE/remote_stage.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
E="$REMOTE_ROOT"; R="$RUNTIME"; F="$FROZEN"
check() { local got; got="\$(sha256sum "\$1" | cut -d' ' -f1)"; [[ "\$got" == "\$2" ]] || { echo "hash mismatch: \$1 got \$got expected \$2" >&2; exit 4; }; }
check "$CHECKPOINT" "$CHECKPOINT_SHA"
check "$REMOTE_SUMMARY" "$SUMMARY_SHA"
check "\$R/plan.json" "$PLAN_SHA"
for pin in $HELPER_PINS; do check "\$R/code/model/tools/\${pin%%:*}" "\${pin##*:}"; done
[[ "\$(ls "\$R/code/model/tools/"*.py | wc -l)" == 9 ]] || { echo "runtime tools dir does not hold exactly the 9 pinned helpers; refusing to stage" >&2; exit 4; }
for pin in $FROZEN_PINS; do check "\$F/\${pin%%:*}" "\${pin##*:}"; done
[[ ! -e "\$E" ]] || { echo "experiment root exists: \$E" >&2; exit 4; }
mkdir -p "\$E/code/model/tools" "\$E/source/code/model" "\$E/results" "\$E/logs" "\$E/numba_cache"
rsync -a --include='*.py' --exclude='*' "\$R/code/model/tools/" "\$E/code/model/tools/"
rsync -a --exclude='__pycache__' --exclude='*.pyc' --exclude='.pytest_cache' --exclude='numba_cache' --exclude='.numba_cache' --exclude='output' --exclude='outputs' --exclude='tmp' --exclude='*.log' --exclude='*.out' --exclude='*.png' --exclude='*.pdf' --exclude='*.py' "\$F/code/model/" "\$E/source/code/model/"
cp -p "\$R/plan.json" "\$E/plan.json"
echo "staged \$E: runtime helpers and frozen non-Python files; awaiting the ported snapshot upload"
EOF

# ---- remote finalize: verify the uploaded snapshot against its manifest and the frozen root ----
cat >"$STAGE/remote_finalize.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
E="$REMOTE_ROOT"; F="$FROZEN"
check() { local got; got="\$(sha256sum "\$1" | cut -d' ' -f1)"; [[ "\$got" == "\$2" ]] || { echo "hash mismatch: \$1 got \$got expected \$2" >&2; exit 4; }; }
check "\$E/source_manifest.json" "$PORT_SHA"
check "\$E/plan.json" "$PLAN_SHA"
check "\$E/code/model/tools/run_e5f_isolated_rental_wedge.py" "$DRIVER_SHA"
for pin in $HELPER_PINS; do check "\$E/code/model/tools/\${pin%%:*}" "\${pin##*:}"; done
[[ "\$(ls "\$E/code/model/tools/"*.py | wc -l)" == 10 ]] || { echo "unexpected helper count in \$E/code/model/tools" >&2; exit 4; }
python3 - "\$E" "\$F" <<'PY'
import hashlib, json, pathlib, sys
E, F = map(pathlib.Path, sys.argv[1:3])
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
port = json.loads((E / "source_manifest.json").read_text()); launch = json.loads((E / "launch_manifest.json").read_text())
pins = port["source_files"]; relevant = set(port["relevant_files"]); base = port["base_hashes_frozen"]
present = {str(p.relative_to(E / "source")) for p in (E / "source/code/model").rglob("*.py") if "__pycache__" not in p.parts}
if present != set(pins): raise SystemExit(f"staged Python set differs from the port manifest: extra={sorted(present - set(pins))} missing={sorted(set(pins) - present)}")
for rel, pin in pins.items():
    if sha(E / "source" / rel) != pin: raise SystemExit(f"staged source hash mismatch: {rel}")
rows = {}
for rel, pin in pins.items():
    fpath = F / rel
    if rel in relevant and rel not in port.get("owned_tests", []):
        if not fpath.is_file(): raise SystemExit(f"frozen root lacks relevant file {rel}")
        fh = sha(fpath)
        if fh != base[rel]: raise SystemExit(f"frozen root {rel} hash {fh} != port base {base[rel]}")
        rows[rel] = {"frozen": fh, "ported": pin, "changed": fh != pin}
    elif fpath.is_file():
        fh = sha(fpath)
        if fh != pin: raise SystemExit(f"non-relevant snapshot file differs from the frozen root: {rel}")
        rows[rel] = {"frozen": fh, "ported": pin, "changed": False}
    else:
        if not rel.startswith("code/model/intergen_eqscale_seq_optimized/tests/"): raise SystemExit(f"snapshot-only file outside tests/: {rel}")
        rows[rel] = {"frozen": None, "ported": pin, "changed": True, "snapshot_only": True}
changed = sorted(r for r, v in rows.items() if v["changed"] and not v.get("snapshot_only"))
expected = sorted(r for r, v in launch["changed_source_vs_frozen"].items() if v["changed"] and v["frozen"] is not None)
if changed != expected: raise SystemExit(f"changed-file set {changed} differs from the launch manifest {expected}")
for rel in changed:
    if launch["changed_source_vs_frozen"][rel] != {"frozen": rows[rel]["frozen"], "ported": rows[rel]["ported"], "changed": True}: raise SystemExit(f"changed-file hashes differ from the launch manifest: {rel}")
nonpy = {str(p.relative_to(E / "source")): sha(p) for p in sorted((E / "source").rglob("*")) if p.is_file() and p.suffix != ".py"}
(E / "frozen_comparison.json").write_text(json.dumps({"frozen_root": str(F), "files": rows, "changed_vs_frozen": changed, "snapshot_only": sorted(r for r, v in rows.items() if v.get("snapshot_only")), "python_files": len(rows), "nonpython_files": len(nonpy)}, indent=2, sort_keys=True) + "\n")
(E / "source_nonpython_manifest.json").write_text(json.dumps(nonpy, indent=2, sort_keys=True) + "\n")
print(f"verified {len(rows)} Python files ({len(changed)} patched vs frozen) and {len(nonpy)} non-Python frozen files")
PY
[[ "\$(find "\$E/source" -name __pycache__ | wc -l)" == 0 ]] || { echo "cache directories present in staged source" >&2; exit 4; }
chmod -R a-w "\$E/code" "\$E/source" "\$E/plan.json" "\$E/launch_manifest.json" "\$E/source_manifest.json" "\$E/smoke.sbatch" "\$E/production.sbatch" "\$E/remote_stage.sh" "\$E/remote_finalize.sh" "\$E/driver_plan.json" "\$E/frozen_comparison.json" "\$E/source_nonpython_manifest.json"
echo "finalized \$E"
EOF

# ---- sbatch files: fresh root, explicit logs/chdir, threads=1, private Numba cache -------------
sbatch_header() {
  cat <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=rw_$1
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
echo "job \${SLURM_JOB_ID:-none} stage $1 start \$(date -u +%FT%TZ) host \$(hostname)"
EOF
}
{ sbatch_header smoke 00:45:00; cat "$STAGE/smoke_body.sh"; echo 'echo "smoke stage complete $(date -u +%FT%TZ)"'; } >"$STAGE/smoke.sbatch"
{ sbatch_header production 01:00:00; cat "$STAGE/production_body.sh"; echo 'echo "production stage complete $(date -u +%FT%TZ)"'; } >"$STAGE/production.sbatch"
for f in "$STAGE/remote_stage.sh" "$STAGE/remote_finalize.sh" "$STAGE/smoke.sbatch" "$STAGE/production.sbatch"; do bash -n "$f"; done
python3 - "$MANIFEST" "$STAGE" <<'PY'
import hashlib, json, sys
from pathlib import Path
m = json.loads(Path(sys.argv[1]).read_text()); stage = Path(sys.argv[2])
m["generated_files"] = {n: hashlib.sha256((stage / n).read_bytes()).hexdigest() for n in ("remote_stage.sh", "remote_finalize.sh", "smoke.sbatch", "production.sbatch", "driver_plan.json")}
Path(sys.argv[1]).write_text(json.dumps(m, indent=2, sort_keys=True) + "\n")
PY
if [[ "$SUBMIT" != 1 ]]; then
  echo "dry-run: syntax review and pins passed; nothing staged or submitted."
  echo "port reviewed marker: $REVIEWED (SUBMIT=1 requires 1)"
  echo "household solves planned: 6 (smoke: cap6zero, cap10s02; production: cap10zero, cap10s005, cap10s1, cap10s02_phi1; combine reuses the verified smoke cap10s02)"
  echo "review: $MANIFEST $STAGE/driver_plan.json $STAGE/remote_stage.sh $STAGE/remote_finalize.sh $STAGE/smoke.sbatch $STAGE/production.sbatch"
  echo "remote root (untouched): $REMOTE_ROOT"; exit 0
fi

# ---- SUBMIT=1: stage remote-to-remote, upload snapshot/driver/manifests, verify, submit ---------
[[ "$REVIEWED" == 1 ]] || { echo "refusing to stage: port manifest is not marked reviewed" >&2; exit 3; }
SSH=(ssh -o BatchMode=yes "$SSH_HOST")
"${SSH[@]}" bash -s <"$STAGE/remote_stage.sh"
rsync -a --exclude='__pycache__' --include='*/' --include='*.py' --exclude='*' --prune-empty-dirs "$SNAPSHOT/code/model/" "$SSH_HOST:$REMOTE_ROOT/source/code/model/"
scp -q "$ROOT/$DRIVER_REL" "$SSH_HOST:$REMOTE_ROOT/code/model/tools/run_e5f_isolated_rental_wedge.py"
scp -q "$PORT_MANIFEST" "$SSH_HOST:$REMOTE_ROOT/source_manifest.json"
scp -q "$MANIFEST" "$STAGE/smoke.sbatch" "$STAGE/production.sbatch" "$STAGE/remote_stage.sh" "$STAGE/remote_finalize.sh" "$STAGE/driver_plan.json" "$SSH_HOST:$REMOTE_ROOT/"
"${SSH[@]}" bash -s <"$STAGE/remote_finalize.sh"
scp -q "$SSH_HOST:$REMOTE_ROOT/frozen_comparison.json" "$STAGE/frozen_comparison.json"
scp -q "$SSH_HOST:$REMOTE_ROOT/source_nonpython_manifest.json" "$STAGE/source_nonpython_manifest.json"
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
