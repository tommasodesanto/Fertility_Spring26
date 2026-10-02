#!/usr/bin/env bash
set -euo pipefail
ssh -4 -oConnectTimeout=20 -oBatchMode=yes torch 'bash /scratch/td2248/projects/transition_readiness_v1/normalized_resumed_fit_v2/floor_launch.sh --mode restore --seconds 300 --label restore_smoke_v2 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/v2/resume_plan.json'
