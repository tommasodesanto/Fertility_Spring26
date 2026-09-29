#!/usr/bin/env bash
#SBATCH --job-name=e5f_credit_no_taper_check
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fixed_reference_economics_20260928/credit_no_taper_v1/slurm_%j.log
set -euo pipefail

# Prepared only; submit after lead review.  This job applies an isolated overlay
# and runs pure parameter/floor tests.  It does not import a checkpoint or solve.
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/credit_no_taper_v1
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'submit with sbatch after lead review' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 NUMBA_DISABLE_JIT=1
cd "$stage/project/$packet"
sha256sum -c source.sha256
export FROZEN_PROJECT_ROOT="$original"
export RENTER_NO_TAPER_OVERLAY="$original/$packet/overlay_v1"
export REFERENCE_LABEL='2007 stationary reference — block0506, September 28 verified export'
python=(apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
"${python[@]}" -c 'import hashlib,json,os,pathlib; root=pathlib.Path(os.environ["FROZEN_PROJECT_ROOT"]); packet=root/"output/model/fixed_reference_economics_20260928/credit_no_taper_v1"; rel=pathlib.Path("code/model/intergen_eqscale_seq_optimized"); job=os.environ["SLURM_JOB_ID"]; files={"parameters.py":root/rel/"parameters.py", "solver.py":root/rel/"solver.py", "kernels.py":root/rel/"kernels.py", "patch_driver.py":packet/"patch_driver.py", "test_renter_no_taper.py":packet/"test_renter_no_taper.py", "run.sh":packet/"run.sh"}; data={"status":"PINNED_PRE_EXECUTION","reference":os.environ["REFERENCE_LABEL"],"hashes":{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in files.items()}}; (packet/("preflight_pins_"+job+".json")).write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")'
"${python[@]}" "$packet/patch_driver.py" --apply --frozen-root "$original" \
  --overlay "$RENTER_NO_TAPER_OVERLAY" --test-file "$original/$packet/test_renter_no_taper.py"
"${python[@]}" "$packet/test_renter_no_taper.py"
