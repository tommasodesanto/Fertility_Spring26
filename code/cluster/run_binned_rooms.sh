#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
umask 077
module load stata/19.0
phase=${1:?Specify toy or full}
task_root=${2:?Specify frozen task directory}
if [[ "$phase" == full ]]; then
    arms=(original_binned aligned_binned)
    arm=${arms[${SLURM_ARRAY_TASK_ID:?Array index required}]}
else
    [[ "$phase" == toy ]]
    arm=toy
fi
outdir="$task_root/results/$arm"
if [[ "$phase" == toy ]]; then
    outdir="$task_root/results/smoke_${SLURM_JOB_ID}/toy"
fi
mkdir -p "$outdir"
cd "$outdir"
started_epoch=$(date +%s)
write_receipt() {
    result_code=$?
    trap - EXIT
    python3 - "$task_root" "$outdir" "$arm" "$started_epoch" "$result_code" <<'PY'
import hashlib,json,pathlib,sys,time
root,out,arm,started,code=sys.argv[1:]
root,out=pathlib.Path(root),pathlib.Path(out)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
passed=int(code)==0 and (out/'completion.txt').exists()
r=dict(arm=arm,route='Torch Stata 19',status='pass' if passed else 'failed',
       exit_code=int(code),elapsed_seconds=time.time()-int(started),
       estimator_do_sha256=sha(root/'audit_binned_rooms.do'),
       data_sha256=sha(root/'analysis_sample.dta'),microdata_uploaded=True)
hp=out/'sample_key_hashes.json'
if hp.exists():
    hashes=json.loads(hp.read_text())
    if 'private_sample_keys.csv' in hashes:r['sample_keys_sha256']=hashes['private_sample_keys.csv']
(out/'run_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
PY
    exit "$result_code"
}
trap write_receipt EXIT
(cd "$task_root" && sha256sum -c SHA256SUMS > "$outdir/source_verification.txt")
cp "$task_root/audit_binned_rooms.do" .
cat > entry.do <<EODO
clear all
set more off
set processors 8
version 17.0
sysdir set PLUS "$task_root/ado/"
adopath ++ "$task_root/ado/"
log using "$outdir/estimation.log", replace text
mata: mata mlib index
EODO
if [[ "$phase" == full ]]; then
    printf 'use "%s/analysis_sample.dta", clear\n' "$task_root" >> entry.do
fi
printf 'do "audit_binned_rooms.do" %s "%s"\n' "$arm" "$outdir" >> entry.do
stata-mp -bq do entry.do &
stata_pid=$!
trap 'kill "$stata_pid" 2>/dev/null || true' TERM INT
while kill -0 "$stata_pid" 2>/dev/null; do
    date -u '+%Y-%m-%dT%H:%M:%SZ' > heartbeat.txt
    sleep 10
done
wait "$stata_pid"
if [[ "$phase" == toy ]]; then
    grep -q '^BINNED_ROOMS_TOY_PASS$' estimation.log
else
    grep -q "^BINNED_ROOMS_ARM_PASS $arm$" estimation.log
fi
python3 - "$outdir" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1])
keys=list(root.rglob('private_sample_keys.csv'))
assert keys
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in keys}
if root.name=='toy':
    assert hashes['original_binned/private_sample_keys.csv']==hashes['aligned_binned/private_sample_keys.csv']
(root/'sample_key_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
for p in keys:
    p.unlink()
PY
printf 'pass\n' > completion.txt
