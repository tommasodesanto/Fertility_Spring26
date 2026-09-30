"""Identity of a same-machine ORIGINAL-engine baseline (stdlib + numpy only).

A baseline directory (written by `acceptance_oracle old-fixed-price`) is usable
for exact comparison only if ALL hold:
  * baseline_receipt.json matches the caller's pin (sha256);
  * engine == "old_frozen_original_sources" and the original solver.py under
    the reference ROOT still has the recorded hash;
  * the recorded checkpoint, manifest and frozen-driver hashes match the files
    under ROOT now;
  * every recorded artifact hash matches (packet, flat solution+shared arrays,
    target_fit.csv, parameters.csv, 17 plots);
  * the recorded runtime identity equals the CURRENT process (OS, architecture,
    processor, Python, NumPy, Numba, SciPy, NumPy BLAS/LAPACK build info).
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

ENGINE = "old_frozen_original_sources"
CHECKPOINT_REL = "output/model/fertility_identification_20260928/resume_v1/selected_export/primary/initial_state.pkl.gz"
MANIFEST_REL = "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
DRIVER_REL = "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py"
SOLVER_REL = "code/model/intergen_eqscale_seq_optimized/solver.py"
ARTIFACTS = ("baseline_packet.pkl.gz", "solution_arrays.npz", "certificate/target_fit.csv", "certificate/parameters.csv")


def sha(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as s:
        for b in iter(lambda: s.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def runtime_identity() -> dict:
    import numba
    import numpy
    import scipy
    try:
        blas = {k: str(v) for k, v in numpy.show_config(mode="dicts").get("Build Dependencies", {}).items()}
    except TypeError:   # NumPy < 1.26 has no dict mode
        blas = {}
    return dict(system=platform.system(), release=platform.release(), machine=platform.machine(),
                processor=platform.processor(), python=platform.python_version(),
                python_implementation=platform.python_implementation(),
                numpy=numpy.__version__, numba=numba.__version__, scipy=scipy.__version__, numpy_build=blas)


def record(root: Path, directory: Path) -> dict:
    """Fields the baseline writer must put into baseline_receipt.json."""
    return dict(engine=ENGINE, model_file=str(root / SOLVER_REL), model_sha256=sha(root / SOLVER_REL),
                checkpoint_sha256=sha(root / CHECKPOINT_REL), manifest_sha256=sha(root / MANIFEST_REL),
                frozen_driver_sha256=sha(root / DRIVER_REL), runtime=runtime_identity(),
                artifacts={a: sha(directory / a) for a in ARTIFACTS},
                plots={p.name: sha(p) for p in sorted((directory / "certificate/standard_diagnostics").glob("*.png"))})


def verify(root: Path, directory: Path, pin: str) -> dict:
    receipt_path = Path(directory) / "baseline_receipt.json"
    if sha(receipt_path) != pin:
        raise RuntimeError("baseline_receipt.json differs from its pin")
    r = json.loads(receipt_path.read_text())
    now = record(Path(root), Path(directory))
    problems = [k for k in ("engine", "model_file", "model_sha256", "checkpoint_sha256", "manifest_sha256",
                            "frozen_driver_sha256", "artifacts", "plots") if r.get(k) != now[k]]
    if r.get("runtime") != now["runtime"]:
        problems.append("runtime (not the same machine/runtime: " + json.dumps(
            {k: (r.get("runtime", {}).get(k), v) for k, v in now["runtime"].items() if r.get("runtime", {}).get(k) != v}) + ")")
    if len(now["plots"]) != 17:
        problems.append(f"plots: {len(now['plots'])} != 17")
    if problems:
        raise RuntimeError("Baseline identity failed: " + "; ".join(problems))
    return r
