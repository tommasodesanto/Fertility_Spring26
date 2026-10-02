#!/usr/bin/env bash
# Submit the native replay gate, then the four fixed-coordinate evaluations.
set -euo pipefail
remote=/scratch/td2248/projects/normalized_floor_extension_v2
ssh -o BatchMode=yes torch "cd '$remote' && test -f launch_torch.sh && test ! -e submission_receipt.txt && FLOOR_MODE=gate sbatch --parsable launch_torch.sh" | tee /tmp/normalized_floor_extension_v2_gate_jobid.txt
gate=$(cat /tmp/normalized_floor_extension_v2_gate_jobid.txt)
ssh -o BatchMode=yes torch "cd '$remote' && FLOOR_MODE=point sbatch --parsable --array=0-3%4 --dependency=afterok:$gate launch_torch.sh" | tee /tmp/normalized_floor_extension_v2_array_jobid.txt
printf 'gate=%s\narray=%s\n' "$gate" "$(cat /tmp/normalized_floor_extension_v2_array_jobid.txt)" | tee /tmp/normalized_floor_extension_v2_submission_receipt.txt
scp /tmp/normalized_floor_extension_v2_submission_receipt.txt "torch:$remote/submission_receipt.txt"
