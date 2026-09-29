#!/usr/bin/env python3
"""Torch-only fail-closed overlay for the frozen renter-debt rule.

This program never writes the frozen source.  It verifies the three supplied
native hashes and writes only ``overlay_v1/.../parameters.py``.  The overlay
adds a default-off flag.  With the flag on it accepts only lambda_d == 0 and
maps the native floor to min(b, 0) before a possible-death decision and to 0
at a possible-death (or terminal) decision.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path


REL = Path("code/model/intergen_eqscale_seq_optimized")
EXPECTED = {
    "solver.py": "b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1",
    "parameters.py": "66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464",
    "kernels.py": "639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27",
}
FLAG = "renter_no_taper_estate_bound"
REFERENCE_LABEL = "2007 stationary reference — block0506, September 28 verified export"
MANIFEST_REL = Path("output/model/fertility_identification_20260928/fixed_reference_manifest.json")
PHYSICAL_TORCH_ROOT = "/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def torch_guard() -> None:
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID", "").isdigit():
        raise RuntimeError("Torch Slurm execution required; local execution is prohibited")


def verify_source(root: Path) -> dict[str, str]:
    manifest_path = root / MANIFEST_REL
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    physical_root = manifest.get("torch_project_root")
    container_root = manifest.get("container_project_root")
    if manifest.get("label") != REFERENCE_LABEL:
        raise RuntimeError("frozen reference manifest label mismatch")
    if physical_root != PHYSICAL_TORCH_ROOT:
        raise RuntimeError("frozen reference manifest physical Torch root mismatch")
    if str(root) not in (physical_root, container_root):
        raise RuntimeError("frozen root is neither the authenticated physical root nor its container mount")
    got = {}
    for name, expected in EXPECTED.items():
        path = root / REL / name
        if not path.is_file():
            raise FileNotFoundError(path)
        got[name] = sha(path)
        if got[name] != expected:
            raise RuntimeError(f"frozen source hash mismatch for {name}: {got[name]}")
    return got


def patch_text(source: str) -> str:
    """Apply one anchored edit to the authenticated parameters source."""
    old_keys = '    "use_tenure_kernel", "w_fixed",\n'
    new_keys = '    "renter_no_taper_estate_bound", "use_tenure_kernel", "w_fixed",\n'
    old_default = "    P.lambda_d = 0.0\n    P.debt_taper_start_age = 42.0\n"
    new_default = (
        "    P.lambda_d = 0.0\n"
        "    # Default-off isolated renter rule.  When enabled, only zero new\n"
        "    # unsecured borrowing is admitted; existing renter debt rolls over\n"
        "    # until the current decision has positive death probability.\n"
        "    P.renter_no_taper_estate_bound = False\n"
        "    P.debt_taper_start_age = 42.0\n"
    )
    old_taper = (
        "    taper = np.ones(J, dtype=float)\n"
        "    middle = (ages > start) & (ages < end)\n"
        "    taper[middle] = (end - ages[middle]) / (end - start)\n"
        "    taper[ages >= end] = 0.0\n\n"
        "    income = np.asarray(P.income, dtype=float)\n"
    )
    new_taper = (
        "    taper = np.ones(J, dtype=float)\n"
        "    if bool(getattr(P, \"renter_no_taper_estate_bound\", False)):\n"
        "        if lam != 0.0:\n"
        "            raise ValueError(\"renter_no_taper_estate_bound requires lambda_d == 0\")\n"
        "        if bool(getattr(P, \"use_age_survival\", False)):\n"
        "            survival = np.asarray(getattr(P, \"survival_probs\"), dtype=float).reshape(-1)\n"
        "            if survival.shape != (J - 1,) or not np.all(np.isfinite(survival)):\n"
        "                raise ValueError(\"survival_probs must be finite with shape (J - 1)\")\n"
        "            if np.any((survival < 0.0) | (survival > 1.0)):\n"
        "                raise ValueError(\"survival_probs must lie in [0, 1]\")\n"
        "            # Array entry j+1 controls the decision at j.  The appended\n"
        "            # terminal entry below is always zero, so this is exactly the\n"
        "            # native condition j == J-1 or survival_probs[j] < 1.\n"
        "            for j in range(J - 1):\n"
        "                if survival[j] < 1.0:\n"
        "                    taper[j + 1] = 0.0\n"
        "    else:\n"
        "        middle = (ages > start) & (ages < end)\n"
        "        taper[middle] = (end - ages[middle]) / (end - start)\n"
        "        taper[ages >= end] = 0.0\n\n"
        "    income = np.asarray(P.income, dtype=float)\n"
    )
    for old, new, label in ((old_keys, new_keys, "dynamic-key"),
                            (old_default, new_default, "default"),
                            (old_taper, new_taper, "taper")):
        if source.count(old) != 1:
            raise RuntimeError(f"expected one {label} anchor, found {source.count(old)}")
        source = source.replace(old, new, 1)
    return source


def apply(root: Path, overlay: Path, test_file: Path | None) -> Path:
    torch_guard()
    source_hashes = verify_source(root)
    target = overlay / REL / "parameters.py"
    if overlay.exists():
        raise FileExistsError(f"overlay must be new: {overlay}")
    target.parent.mkdir(parents=True, exist_ok=False)
    changed = patch_text((root / REL / "parameters.py").read_text())
    target.write_text(changed)
    manifest = {
        "status": "APPLIED_NOT_APPROVED",
        "reference": REFERENCE_LABEL,
        "frozen_root": str(root),
        "overlay": str(overlay),
        "source_hashes": source_hashes,
        "effective_parameters_sha256": sha(target),
        "patch_driver_sha256": sha(Path(__file__)),
        "test_sha256": sha(test_file) if test_file else None,
        "economic_change": "Default-off renter estate-bound rule; lambda_d must equal zero when active.",
        "uncomputed": ["No lifecycle solve", "No KFE/GE solve", "No recalibration", "No diagnostic plots regenerated"],
    }
    (overlay / "overlay_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return target


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--frozen-root", type=Path, required=True)
    ap.add_argument("--overlay", type=Path, required=True)
    ap.add_argument("--test-file", type=Path)
    ns = ap.parse_args()
    if not ns.apply:
        raise SystemExit("only --apply is supported")
    target = apply(ns.frozen_root.resolve(), ns.overlay.resolve(),
                   ns.test_file.resolve() if ns.test_file else None)
    print(json.dumps({"status": "APPLIED_NOT_APPROVED", "parameters": str(target)}, sort_keys=True))


if __name__ == "__main__":
    main()
