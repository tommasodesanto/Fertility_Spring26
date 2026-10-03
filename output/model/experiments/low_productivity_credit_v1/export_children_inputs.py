"""Export unchanged saved distributions for the shared child-count plotter.

No model solve or alteration to the native result is performed here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "code/model"))
from production.storage import load_case


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(name: str, expected_phi: float):
    case = HERE / name
    saved, _ = load_case(case)
    P, sol = saved.P, saved.solution
    if not np.allclose(P.phi, expected_phi, rtol=0, atol=1e-12):
        raise RuntimeError(f"Wrong saved financing share in {case}")
    if not np.isfinite(sol.g).all() or not np.isclose(sol.g.sum(), 1, rtol=0, atol=1e-8):
        raise RuntimeError(f"Invalid saved mass in {case}")
    dest = HERE / "children_by_age" / name
    dest.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest / "solution_arrays.npz",
                        g=sol.g, g_beginning_distribution=sol.g_beginning_distribution)
    # The shared plotter's streaming reader expects the phi array to be followed
    # by another field (a closing `],` line in the pretty JSON representation).
    fields = {"phi": np.asarray(P.phi).tolist(), "J": int(P.J), "age_start": float(P.age_start),
              "da": float(P.da), "n_parity": int(P.n_parity),
              "use_postdecision_current_distribution": bool(P.use_postdecision_current_distribution)}
    (dest / "executed_P.json").write_text(json.dumps(fields, indent=2) + "\n")
    (dest / "provenance.json").write_text(json.dumps({
        "source_case": str(case.resolve()),
        "source_native_result_sha256": sha(case / "native_result.npz"),
        "exported_solution_arrays_sha256": sha(dest / "solution_arrays.npz"),
        "population_mass": float(sol.g.sum()),
        "outside_low_income_mass": float(sol.g.sum(axis=(0, 1, 2, 3, 5, 6))[1:].sum()),
        "scope": "unchanged native distribution; no buyer-policy or aggregate interpolation",
    }, indent=2) + "\n")
    return dest


if __name__ == "__main__":
    for name, phi in (("phi_08", 0.8), ("phi_10", 1.0)):
        print(export(name, phi))
