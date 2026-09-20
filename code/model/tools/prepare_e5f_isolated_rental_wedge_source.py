#!/usr/bin/env python3
"""Build the hash-pinned frozen-source rental-wedge snapshot.

The builder copies only the retained September source, applies the committed
review patch, and writes a manifest/diff. It never imports active model code
or runs a lifecycle solve. Existing output is never overwritten unless
``--force`` is explicit; a separate scratch directory is recommended for
verification.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SNAPSHOT = ROOT / "tmp/e5f_rental_wedge_runtime_v4"
DEFAULT_FROZEN = ROOT / "tmp/e5f_complete_frozen_source"
DEFAULT_ACTIVE = ROOT
DEFAULT_PATCH = ROOT / "code/model/tools/patches/e5f_isolated_rental_wedge.patch"
PACKAGE = "code/model/intergen_eqscale_seq_optimized"
RELEVANT = (f"{PACKAGE}/parameters.py", f"{PACKAGE}/kernels.py", f"{PACKAGE}/solver.py")
HELPERS = (
    "code/model/tools/run_e5f_matched_pf_smoke.py",
    "code/model/tools/run_dynamic_population_transition.py",
    "code/model/tools/run_e5f_open_population_transition.py",
    "code/model/tools/run_e5f_perfect_foresight_transition.py",
    "code/model/tools/run_e5f_independent_numerical_audit.py",
)
OWNED_TESTS = (f"{PACKAGE}/tests/test_switch_wedge.py", f"{PACKAGE}/tests/test_rental_wedge_port_tiny.py")
DIFF_FILES = RELEVANT + HELPERS + OWNED_TESTS
BASE_MANIFEST = ROOT / "output/model/paper_baseline_sep14/main_expected_source_manifest.json"
FULL_FROZEN_MANIFEST = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/credit_policy_retention_v2/source_manifest.json"
EXCLUDED_INACTIVE_BUILDERS = {
    "code/model/tools/build_e5f_bounded_refinement_plan.py",
    "code/model/tools/build_simplified_olg_theory_slides.py",
}
CHECKPOINT_SHA256 = "3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} missing: {path}")


def hash_files(root: Path, relative_paths: tuple[str, ...]) -> dict[str, str]:
    result = {}
    for relative in relative_paths:
        path = root / relative
        require_file(path, relative)
        result[relative] = sha256(path)
    return result


def all_snapshot_python_hashes(snapshot: Path) -> dict[str, str]:
    paths = sorted((snapshot / "code/model").rglob("*.py"))
    if not paths:
        raise FileNotFoundError(f"snapshot has no Python source: {snapshot}")
    return {str(path.relative_to(snapshot)): sha256(path) for path in paths}


def copied_python_inventory(frozen: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Use the complete previously verified frozen snapshot, including pickle helpers."""
    pins = json.loads(FULL_FROZEN_MANIFEST.read_text())
    relative = tuple(sorted("code/model/" + name for name in pins))
    present = {str(p.relative_to(frozen)) for p in (frozen / "code/model").rglob("*.py")}
    if present != set(relative):
        raise RuntimeError("complete frozen Python inventory differs from its retained manifest")
    # These unused archive builders differ between the historical main manifest
    # and the verified full snapshot. Do not overwrite either pin or stage them.
    relative = tuple(name for name in relative if name not in EXCLUDED_INACTIVE_BUILDERS)
    tests = tuple(name for name in relative if "/tests/" in name or Path(name).name.startswith("test_"))
    production = tuple(name for name in relative if name not in tests)
    return production, tests


def frozen_base_hashes(frozen: Path) -> dict[str, str]:
    """Verify the retained manifest against the actual frozen files first."""
    production, _ = copied_python_inventory(frozen)
    actual = hash_files(frozen, production)
    if BASE_MANIFEST.is_file():
        expected_source = json.loads(BASE_MANIFEST.read_text()).get("source_files", {})
        complete_pins = json.loads(FULL_FROZEN_MANIFEST.read_text())
        for relative in production:
            short = relative.removeprefix("code/model/")
            if short in complete_pins:
                expected_source.setdefault(relative, complete_pins[short])
        missing = [relative for relative in production if relative not in expected_source]
        if missing:
            raise RuntimeError(f"retained frozen manifest omits source files: {missing}")
        mismatch = {relative: {"manifest": str(expected_source[relative]), "actual": actual[relative]}
                    for relative in actual if str(expected_source[relative]) != actual[relative]}
        if mismatch:
            raise RuntimeError(f"retained frozen manifest disagrees with source: {mismatch}")
    return actual


def apply_patch(snapshot: Path, patch_path: Path) -> None:
    require_file(patch_path, "rental-wedge patch")
    result = subprocess.run(
        ["patch", "--batch", "--forward", "--fuzz=0", "-p1", "-i", str(patch_path)],
        cwd=snapshot, text=True, capture_output=True,
    )
    if result.returncode:
        raise RuntimeError(f"patch failed ({result.returncode}):\n{result.stdout}\n{result.stderr}")


def build_snapshot(snapshot: Path, frozen: Path, patch_path: Path, *, force: bool) -> None:
    frozen_base_hashes(frozen)
    active = ROOT.resolve()
    frozen = frozen.resolve()
    snapshot = snapshot.resolve()
    if snapshot == active or snapshot == frozen or active.is_relative_to(snapshot) or frozen.is_relative_to(snapshot) or snapshot.is_relative_to(active / "code") or snapshot.is_relative_to(frozen):
        raise RuntimeError(f"refusing snapshot target at/above active or frozen root: {snapshot}")
    if snapshot.exists():
        if not force:
            raise FileExistsError(f"refusing to overwrite existing snapshot: {snapshot}; use --force or a fresh path")
        shutil.rmtree(snapshot)
    (snapshot / "code/model").mkdir(parents=True)
    production, tests = copied_python_inventory(frozen)
    for relative in (*production, *tests):
        destination = snapshot / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(frozen / relative, destination)
    apply_patch(snapshot, patch_path)


def unified_diff(snapshot: Path, frozen: Path) -> str:
    lines: list[str] = []
    for relative in DIFF_FILES:
        frozen_path, snapshot_path = frozen / relative, snapshot / relative
        require_file(snapshot_path, f"snapshot {relative}")
        old = frozen_path.read_text().splitlines(keepends=True) if frozen_path.exists() else []
        lines.extend(difflib.unified_diff(
            old, snapshot_path.read_text().splitlines(keepends=True),
            fromfile=f"frozen/{relative}", tofile=f"snapshot/{relative}",
        ))
    return "".join(lines)


def assert_zero_slope_legacy_fragments(snapshot: Path) -> list[str]:
    text = "".join((snapshot / relative).read_text() for relative in RELEVANT)
    required = {
        "parameters rental defaults": "P.rental_wedge_slope = 0.0",
        "parameters active helper": "def rental_wedge_active",
        "kernel legacy renter candidate": "surplus = Rv - dc - bp",
        "kernel legacy cap arithmetic": "ct = Rvb - cbc - ri * hR_max - bp_best",
        "solver legacy golden branch": "bp, val = golden_renter(",
        "solver wedge rejection": "Rental wedge requires the golden-section renter block",
    }
    missing = [label for label, fragment in required.items() if fragment not in text]
    if missing:
        raise RuntimeError(f"zero-slope legacy fragments missing: {missing}")
    return sorted(required)


def verify_against_manifest(snapshot: Path, reference: Path) -> None:
    expected = json.loads(reference.read_text()).get("source_files", {})
    actual = all_snapshot_python_hashes(snapshot)
    mismatch = {
        relative: {"expected": expected[relative], "actual": actual.get(relative)}
        for relative in expected if expected[relative] != actual.get(relative)
    }
    if mismatch or set(actual) != set(expected):
        raise RuntimeError(f"built snapshot hashes differ from reference: {mismatch}")


def write_manifest(snapshot: Path, frozen: Path, active: Path, patch_path: Path,
                   base_hash: dict[str, str]) -> dict[str, Any]:
    source_hashes = all_snapshot_python_hashes(snapshot)
    port_hashes = hash_files(snapshot, RELEVANT + HELPERS + OWNED_TESTS)
    active_hashes = hash_files(active, RELEVANT)
    copied_production, frozen_tests = copied_python_inventory(frozen)
    frozen_test_hashes = {relative: sha256(frozen / relative) for relative in frozen_tests}
    retained_tests = set(json.loads(BASE_MANIFEST.read_text()).get("source_files", {})) if BASE_MANIFEST.is_file() else set()
    diff_path = snapshot / "rental_wedge_port.diff"
    diff_path.write_text(unified_diff(snapshot, frozen))
    manifest = {
        "schema": "e5f_isolated_rental_wedge_source_v2",
        "source_snapshot": str(snapshot), "frozen_source_root": str(frozen),
        "active_checkout": str(active), "patch_path": str(patch_path),
        "patch_sha256": sha256(patch_path), "rental_wedge_port_reviewed": False,
        "complete_frozen_pin_manifest": str(FULL_FROZEN_MANIFEST),
        "complete_frozen_pin_manifest_sha256": sha256(FULL_FROZEN_MANIFEST),
        "excluded_inactive_builders": sorted(EXCLUDED_INACTIVE_BUILDERS),
        "control_required_before_wedge": True, "checkpoint_sha256": CHECKPOINT_SHA256,
        "base_hashes_frozen": base_hash, "ported_hashes": port_hashes,
        "copied_production_python_files": list(copied_production),
        "frozen_test_hashes_informational": frozen_test_hashes,
        "frozen_tests_missing_from_retained_manifest": sorted(set(frozen_tests) - retained_tests),
        "active_hashes_relevant": active_hashes, "source_files": source_hashes,
        "required_helpers": list(HELPERS), "owned_tests": list(OWNED_TESTS),
        "relevant_files": list(DIFF_FILES),
        "zero_slope_legacy_fragments": assert_zero_slope_legacy_fragments(snapshot),
        "port_scope": {
            "parameters": "rental_wedge_intercept/slope/knee defaults, validation, active predicate, total cost",
            "kernels": "renter candidate branches, saving objective, committed cost, cap output",
            "solver": "Markov-income renter branch only; core path rejects active wedge",
            "audits": "snapshot-local dated and independent audits use full rental_wedge_total_cost; active-wedge saving audit enumerates every b_grid segment with endpoint checks and an independent bounded solve per segment",
        },
        "remaining_review_blockers": [
            "Line-by-line verify quadratic above-knee candidate and cap consumption against C(h).",
            "Review snapshot-local dated and independent budget gates and array broadcasting before any lifecycle run.",
            "Verify Python and Numba branches on cap, knee, and infeasible committed-housing states.",
            "Run the original frozen control and record mandatory policy/mass comparisons before wedge arms.",
        ],
        "diff_path": str(diff_path),
    }
    manifest_path = snapshot / "source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {"manifest": str(manifest_path), "diff": str(diff_path), "files": len(source_hashes),
            "diff_lines": len(diff_path.read_text().splitlines())}


def prepare(snapshot: Path = DEFAULT_SNAPSHOT, frozen: Path = DEFAULT_FROZEN,
            active: Path = DEFAULT_ACTIVE, patch_path: Path = DEFAULT_PATCH,
            *, build: bool = False, force: bool = False,
            verify_against: Path | None = None) -> dict[str, Any]:
    snapshot, frozen, active, patch_path = map(Path.resolve, (snapshot, frozen, active, patch_path))
    base_hash = frozen_base_hashes(frozen)
    if build:
        build_snapshot(snapshot, frozen, patch_path, force=force)
    else:
        require_file(snapshot / RELEVANT[0], "existing snapshot")
        require_file(patch_path, "rental-wedge patch")
    if verify_against is not None:
        verify_against_manifest(snapshot, Path(verify_against).resolve())
    return write_manifest(snapshot, frozen, active, patch_path, base_hash)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--frozen", type=Path, default=DEFAULT_FROZEN)
    parser.add_argument("--active", type=Path, default=DEFAULT_ACTIVE)
    parser.add_argument("--patch", type=Path, default=DEFAULT_PATCH)
    parser.add_argument("--build", action="store_true", help="copy frozen source and apply the patch")
    parser.add_argument("--force", action="store_true", help="allow replacing an existing snapshot during --build")
    parser.add_argument("--verify-against", type=Path, help="compare the fresh build with this prior source manifest")
    parser.add_argument("--dry-run", action="store_true", help="verify frozen hashes and patch presence without writing")
    args = parser.parse_args(argv)
    frozen = args.frozen.resolve()
    base_hash = frozen_base_hashes(frozen)
    require_file(args.patch.resolve(), "rental-wedge patch")
    if args.dry_run:
        print(json.dumps({"dry_run": True, "frozen_hashes": base_hash,
                          "patch_sha256": sha256(args.patch.resolve())}, indent=2, sort_keys=True))
        return 0
    print(json.dumps(prepare(args.snapshot, args.frozen, args.active, args.patch,
                             build=args.build, force=args.force,
                             verify_against=args.verify_against), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
