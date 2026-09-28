#!/usr/bin/env bash
# Run locally only after final freeze; transfers small reviewed source files.
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
scp code/model/tools/run_e5f_night_calibration.py code/model/tools/prepare_e5f_night_contract.py code/model/tools/test_e5f_night_calibration.py "torch:$stage/project/code/model/tools/"
scp code/cluster/submit_e5f_night_calibration.sh code/cluster/submit_e5f_night_gated.sh "torch:$stage/project/code/cluster/"
scp output/model/overnight_calibration_20260928/cluster/test_and_prepare.sh "torch:$stage/test_and_prepare.sh"
ssh torch 'stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1; sha256sum "$stage/project/code/model/tools/run_e5f_night_calibration.py" "$stage/project/code/model/tools/prepare_e5f_night_contract.py" "$stage/project/code/model/tools/test_e5f_night_calibration.py" "$stage/project/code/cluster/submit_e5f_night_calibration.sh" "$stage/project/code/cluster/submit_e5f_night_gated.sh" "$stage/test_and_prepare.sh" > "$stage/refresh_manifest.sha256"; cat "$stage/refresh_manifest.sha256"; sbatch --parsable --output="$stage/test_prepare_%j.log" "$stage/test_and_prepare.sh"'
