#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:50:00
#SBATCH --job-name=refactor_lab
# Slurm may spool this script, so the source checkout must be explicit.
: "${LAB_SRC:?set LAB_SRC to the directory containing refactor_lab}"
DRIVER="$LAB_SRC/refactor_lab/verification/verify.sh"
if [ ! -f "$DRIVER" ]; then
  echo "verification driver not found: $DRIVER" >&2
  exit 2
fi
MODE=torch exec bash "$DRIVER" "$@"
