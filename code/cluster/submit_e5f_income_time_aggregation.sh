#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TAG="${TAG:-income_aggregation_v1}"; SUBMIT="${SUBMIT:-0}"
REMOTE_ROOT="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/$TAG"
LOCAL_ROOT="$ROOT/output/model/native_financing_diagnostic_20260919/specification_followup/$TAG"
[[ "$TAG" =~ ^[A-Za-z0-9._-]+$ ]] || { echo invalid_TAG >&2; exit 2; }
[[ "$TAG" == income_aggregation_v1 || "${ALLOW_TAG_OVERRIDE:-0}" == 1 ]] || { echo TAG_override_requires_ALLOW_TAG_OVERRIDE >&2; exit 2; }
[[ ! -e "$LOCAL_ROOT" ]] || { echo "refusing existing local output: $LOCAL_ROOT" >&2; exit 2; }
DRIVER=code/model/tools/diagnose_e5f_income_time_aggregation.py
CONSTRUCTOR=code/model/tools/build_persistent_transitory_income_candidate.py
CANDIDATE=output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json
AUTOCOV=code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/md_autocovariance_fit.csv
PACKAGE=code/model/intergen_eqscale_seq_optimized
for f in "$DRIVER" "$CONSTRUCTOR" "$CANDIDATE" "$AUTOCOV"; do [[ -f "$ROOT/$f" ]] || { echo "missing input: $f" >&2; exit 2; }; done
if [[ "$SUBMIT" == 1 ]]; then mkdir -p "$LOCAL_ROOT"; else LOCAL_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/e5f_income_time_aggregation.XXXXXX")"; fi
MANIFEST="$LOCAL_ROOT/launch_manifest.json"
python3 - "$ROOT" "$MANIFEST" "$DRIVER" "$CONSTRUCTOR" "$CANDIDATE" "$AUTOCOV" "$PACKAGE" <<'PY'
import hashlib,json,pathlib,sys
root,out,*items=map(pathlib.Path,sys.argv[1:]); files=[]
for x in items: files.extend(sorted((root/x).rglob('*.py')) if (root/x).is_dir() else [root/x])
rows=[{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
payload={'design':{'smoke':[2,2000,40],'full':[20,20000,120],'seed':20260920,'internal_seconds':{'smoke':60,'full':720},'process_seconds':{'smoke':540,'full':840},'zero_household_model_or_equilibrium_solves':True,'hypothesis':'Quantify differences between four-year averages of annual lognormal income and the current endpoint approximation; no process adoption or recalibration.','stop_rule':'Stop at failed validity gate or budget; never retry automatically.'},'inputs':rows}
pathlib.Path(out).write_text(json.dumps(payload,indent=2)+'\n')
PY
if [[ "$SUBMIT" != 1 ]]; then printf 'dry-run manifest: %s\nremote root: %s\n' "$MANIFEST" "$REMOTE_ROOT"; exit 0; fi
REMOTE='ssh -o BatchMode=yes torch'; $REMOTE "test ! -e '$REMOTE_ROOT' && mkdir -p '$REMOTE_ROOT'" || { echo remote_root_exists_or_unavailable >&2; exit 3; }
ARCHIVE="$LOCAL_ROOT/source_snapshot.tar"
python3 - "$ROOT" "$MANIFEST" "$ARCHIVE" <<'PY'
import json,pathlib,sys,tarfile
root=pathlib.Path(sys.argv[1])
with tarfile.open(sys.argv[3],'w') as archive:
 for row in json.load(open(sys.argv[2]))['inputs']:
  archive.add(root/row['path'],arcname=row['path'],recursive=False)
PY
scp -q "$ARCHIVE" "$MANIFEST" "torch:$REMOTE_ROOT/"
$REMOTE "cd '$REMOTE_ROOT' && tar -xf source_snapshot.tar && rm source_snapshot.tar && python3 - <<'PY'
import hashlib,json,pathlib
for r in json.load(open('launch_manifest.json'))['inputs']:
 p=pathlib.Path(r['path']); assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']
PY"
cat > "$LOCAL_ROOT/smoke.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name=income_agg_smoke
#SBATCH --chdir=$REMOTE_ROOT
#SBATCH --output=$REMOTE_ROOT/%x-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$REMOTE_ROOT"
python3 - <<'PY'
import hashlib,json,pathlib
for r in json.load(open('launch_manifest.json'))['inputs']:
 p=pathlib.Path(r['path']); assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']
PY
timeout 540 python3 -B $DRIVER --mode smoke --candidate-json $CANDIDATE --output output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/smoke
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/smoke/receipt.json
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/smoke/level_covariance_comparison.png
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/smoke/block_average_quantiles.png
EOF
scp -q "$LOCAL_ROOT/smoke.sbatch" "torch:$REMOTE_ROOT/smoke.sbatch"
SMOKE="$($REMOTE "sbatch --parsable '$REMOTE_ROOT/smoke.sbatch'")"
python3 - "$LOCAL_ROOT/submission.json" "$SMOKE" "$REMOTE_ROOT" "$MANIFEST" <<'PY'
import json,sys,pathlib
pathlib.Path(sys.argv[1]).write_text(json.dumps({'smoke_job':sys.argv[2],'remote_root':sys.argv[3],'manifest':sys.argv[4]},indent=2)+'\n')
PY
cat > "$LOCAL_ROOT/full.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=income_agg_full
#SBATCH --chdir=$REMOTE_ROOT
#SBATCH --output=$REMOTE_ROOT/%x-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$REMOTE_ROOT"
python3 - <<'PY'
import hashlib,json,pathlib
for r in json.load(open('launch_manifest.json'))['inputs']:
 p=pathlib.Path(r['path']); assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']
s=pathlib.Path('output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/smoke'); receipt=json.loads((s/'receipt.json').read_text()); assert receipt.get('status')=='completed'
for key in ('all_mean_one_checks_pass','all_exact_level_moment_checks_pass','all_dependencies_unchanged','all_batches_completed','all_required_plots_present'): assert receipt['validity'][key] is True,key
for n in ('level_covariance_comparison.png','block_average_quantiles.png'): assert (s/n).is_file()
PY
timeout 840 python3 -B $DRIVER --mode full --candidate-json $CANDIDATE --output output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/full
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/full/receipt.json
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/full/progress.json
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/full/level_covariance_comparison.png
test -s output/model/native_financing_diagnostic_20260919/specification_followup/$TAG/results/full/block_average_quantiles.png
EOF
scp -q "$LOCAL_ROOT/full.sbatch" "torch:$REMOTE_ROOT/full.sbatch"
FULL="$($REMOTE "sbatch --parsable --dependency=afterok:$SMOKE '$REMOTE_ROOT/full.sbatch'")"
python3 - "$LOCAL_ROOT/submission.json" "$SMOKE" "$FULL" <<'PY'
import json,sys,pathlib
p=pathlib.Path(sys.argv[1]); d=json.loads(p.read_text()); d['full_job']=sys.argv[3]; p.write_text(json.dumps(d,indent=2)+'\n')
PY
printf 'submitted smoke=%s full=%s\n' "$SMOKE" "$FULL"
