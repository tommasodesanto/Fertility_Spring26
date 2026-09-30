#!/usr/bin/env bash
# Stage only the isolated driver and copied small engine; does not submit jobs.
set -euo pipefail
packet="$(cd "$(dirname "$0")" && pwd)"
remote=/scratch/td2248/projects/small_credit_v1
staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT
mkdir -p "$staging/source"
rsync -a --exclude='__pycache__' --exclude='*.pyc' "$packet/source/" "$staging/source/"
cp "$packet/driver.py" "$packet/phase_a.py" "$packet/single_price.py" \
  "$packet/phase_b_ge.py" "$packet/plan.json" "$packet/launch_torch.sh" "$staging/"
(cd "$staging" && find . -type f | sort | xargs sha256sum > source.sha256)
ssh torch "test ! -e '$remote/source' && mkdir -p '$remote/source' '$remote/results' '$remote/cache' '$remote/logs'" || {
  echo 'remote source already exists; refusing overwrite' >&2; exit 2;
}
rsync -a "$staging/" "torch:$remote/source/"
ssh torch "cd '$remote/source' && sha256sum -c source.sha256 && test \"\$(sha256sum /scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs/bundle.json | cut -d' ' -f1)\" = 427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7 && chmod -R a-w '$remote/source'"
echo "Staged and checked $remote/source; no job submitted."
