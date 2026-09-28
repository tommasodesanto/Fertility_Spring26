#!/usr/bin/env bash
set -euo pipefail
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
cd "$original"
python=/share/apps/anaconda3/2025.06/bin/python
"$python" output/model/fertility_identification_20260928/prepare.py
contract="$original/output/model/fertility_identification_20260928/contract_v1/contract.json"
export EXPECTED_E5F_IDENTIFICATION_SHA256="$(sha256sum "$contract" | cut -d' ' -f1)"
"$python" code/model/tools/run_e5f_fertility_identification.py --stage prepare --contract "$contract" --output output/model/fertility_identification_20260928/verified_v1
