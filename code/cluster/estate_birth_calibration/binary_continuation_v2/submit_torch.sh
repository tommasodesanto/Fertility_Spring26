#!/usr/bin/env bash
# Submit fixed 20-chain one-birth array. No retries or extensions are configured.
set -euo pipefail
remote=/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2
ssh -o BatchMode=yes torch "/share/apps/anaconda3/2025.06/bin/python '$remote/submit_production.py'"
