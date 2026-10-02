#!/usr/bin/env bash
set -euo pipefail
folder=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp "$folder/progress.XXXXXX")
trap 'rm -f "$tmp"' EXIT
ssh -4 -oConnectTimeout=20 -oBatchMode=yes torch '/share/apps/anaconda3/2025.06/bin/python -' < "$folder/compact_progress.py" > "$tmp"
python3 -m json.tool "$tmp" >/dev/null
mv "$tmp" "$folder/latest_progress.json"
cat "$folder/latest_progress.json"
