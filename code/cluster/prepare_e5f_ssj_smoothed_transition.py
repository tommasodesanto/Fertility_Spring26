"""Run ON TORCH: write manifest + sbatch for one smoothed-transition experiment batch.

Usage: python prepare_e5f_ssj_smoothed_transition.py <batch_dir> <kappa> [seconds] [slurm_minutes]
The batch_dir must already contain source/ with the three driver files.
"""
import hashlib, json, sys
from pathlib import Path

R = "/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches"
ANN = f"{R}/announced_original_queue_20260913c"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    batch = Path(sys.argv[1]); kappa = float(sys.argv[2])
    seconds = int(sys.argv[3]) if len(sys.argv) > 3 else 31000
    minutes = int(sys.argv[4]) if len(sys.argv) > 4 else 540
    import os
    start_env = os.environ.get("E5F_SMOOTH_STATIONARY_START", "")
    stationary_start = [float(x) for x in start_env.split(",")] if start_env else None
    stationary_evaluations = int(os.environ.get("E5F_SMOOTH_STATIONARY_EVALS", "16"))
    receipt = f"{ANN}/output/run/root_receipt.json"
    ckpt = f"{R}/announced_original_queue_20260913c_ssj_rescue_toeplitz/run/best_so_far.json"
    toep = f"{R}/afternoon_original_queue_20260913a_ssj_toeplitz_10/derivative/derivative_receipt.json"
    sources = ["run_e5f_ssj_smoothed_transition.py", "run_e5f_ssj_announced_rescue.py", "e5f_ssj_toeplitz_jacobian.py",
               "e5f_ssj_scaled_step_root.py"]
    for s in sources:
        if not (batch / "source" / s).exists():
            raise SystemExit("missing source file: " + s)
    pins = [str(batch / "source" / s) for s in sources] + [f"{ANN}/manifest.json", receipt, ckpt, toep,
                                                           f"{ANN}/source/run_e5f_announced_original_queue.py"]
    pins = {k: sha(k) for k in pins}
    m = dict(announced_manifest=f"{ANN}/manifest.json", announced_source_dir=f"{ANN}/source",
             tenure_choice_kappa=kappa, warm_start_checkpoint=ckpt, toeplitz_receipt=toep,
             jacobian_mode="toeplitz", step_rule="scaled", mapping_budget=12, trim_count=0,
             skip_mapping_plots=True, output=str(batch), seconds=seconds, stationary_seconds=3600,
             stationary_start=stationary_start, stationary_evaluations=stationary_evaluations,
             author_request=f"Experiment only, author-approved 2026-09-16: re-solve stationary and terminal equilibria at tenure_choice_kappa={kappa}, then the 104-date announced root with the measured-Jacobian start at the retained gates. No production change.",
             file_sha256=pins)
    json.dump(m, open(batch / "manifest.json", "w"), indent=2, sort_keys=True)
    (batch / "run.sbatch").write_text(f"""#!/bin/bash
#SBATCH --job-name=e5f_ssj_smooth_k{kappa}
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time={minutes}
#SBATCH --output={batch}/run_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/share/apps/anaconda3/2025.06/bin/python {batch}/source/run_e5f_ssj_smoothed_transition.py --manifest {batch}/manifest.json
""")
    print(json.dumps(dict(batch=str(batch), kappa=kappa, manifest_sha256=sha(batch / "manifest.json"),
                          script_sha256=sha(batch / "run.sbatch"))))


if __name__ == "__main__":
    main()
