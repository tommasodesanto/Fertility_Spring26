#!/usr/bin/env python3
"""Small, fail-closed before/after comparator for E5F refactor evidence.

The command compares two explicitly supplied artifacts, or two directories whose
``manifest.json`` enumerates the artifacts.  It is deliberately format-limited:
JSON, CSV, NPY, and NPZ are supported; pickle based formats are rejected.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


class ComparisonError(AssertionError):
    """Raised when evidence does not satisfy the exact comparison contract."""


# These are the only administrative fields exempted inside the manifest's
# top-level run_metadata object. Artifact JSON is always compared exactly.
MANIFEST_RUN_METADATA_IGNORED_FIELDS = frozenset(
    {"created_at", "timestamp", "started_at", "finished_at", "elapsed_seconds",
     "wall_seconds", "runtime_seconds", "duration_seconds", "path", "source_path",
     "output_path", "working_directory", "hostname", "pid"}
)
SUPPORTED_SCHEMAS = frozenset({"e5f_refactor_evidence_v1"})


def _sha256(path: Path) -> str:
    try:
        h = hashlib.sha256()
        with path.open("rb") as f:
            for block in iter(lambda: f.read(1024 * 1024), b""):
                h.update(block)
        return h.hexdigest()
    except OSError as exc:
        raise ComparisonError(f"cannot hash {path}: {exc}") from exc


def _json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ComparisonError(f"cannot read JSON {path}: {exc}") from exc


def _normalise(value: Any, *, manifest: bool = False, path: tuple[str, ...] = ()) -> Any:
    """Canonicalize JSON while preserving scalar types.

    In manifest mode, only provenance.source_identity, direct artifact path and
    sha256 claims, and allowlisted keys in top-level run_metadata are omitted.
    """
    if isinstance(value, dict):
        out = {}
        for key, item in sorted(value.items()):
            if manifest and path == ("provenance",) and key == "source_identity":
                continue
            if manifest and len(path) == 2 and path[0] == "artifacts" and key in {"path", "sha256"}:
                continue
            if manifest and path == ("run_metadata",) and key in MANIFEST_RUN_METADATA_IGNORED_FIELDS:
                continue
            out[key] = _normalise(item, manifest=manifest, path=path + (key,))
        return out
    if isinstance(value, list):
        return [_normalise(v, manifest=manifest, path=path + (str(i),)) for i, v in enumerate(value)]
    # JSON's native equality conflates True with 1 and 1 with 1.0. Preserve
    # scalar types because gates and discrete decisions are scientific data.
    if isinstance(value, bool):
        return ("bool", value)
    if isinstance(value, int):
        return ("int", value)
    if isinstance(value, float):
        return ("float", repr(value))
    if isinstance(value, str):
        return ("str", value)
    if value is None:
        return ("null", None)
    return (type(value).__name__, repr(value))


def _compare_json(a: Path, b: Path, label: str) -> list[str]:
    left, right = _json(a), _json(b)
    failures: list[str] = []
    if isinstance(left, dict) and isinstance(right, dict):
        if "provenance" in left or "provenance" in right:
            if not isinstance(left.get("provenance"), dict) or not isinstance(right.get("provenance"), dict):
                failures.append(f"{label}: provenance must be an object on both sides")
    if _normalise(left) != _normalise(right):
        failures.append(f"{label}: JSON values differ")
    return failures


def _compare_array(a: Path, b: Path, label: str) -> list[str]:
    try:
        x = np.load(a, allow_pickle=False)
        y = np.load(b, allow_pickle=False)
    except Exception as exc:  # includes object arrays and malformed files
        raise ComparisonError(f"{label}: NPY/NPZ load failed or pickle/object data forbidden: {exc}") from exc
    if isinstance(x, np.lib.npyio.NpzFile) or isinstance(y, np.lib.npyio.NpzFile):
        if not isinstance(x, np.lib.npyio.NpzFile) or not isinstance(y, np.lib.npyio.NpzFile):
            return [f"{label}: one side is NPZ and the other is not"]
        try:
            if set(x.files) != set(y.files):
                return [f"{label}: NPZ keys differ: {x.files} vs {y.files}"]
            out: list[str] = []
            for key in sorted(x.files):
                try:
                    out.extend(_array_values(x[key], y[key], f"{label}[{key}]"))
                except Exception as exc:
                    raise ComparisonError(f"{label}[{key}]: object/pickle NPZ data forbidden: {exc}") from exc
            return out
        finally:
            x.close(); y.close()
    return _array_values(x, y, label)


def _array_values(x: np.ndarray, y: np.ndarray, label: str) -> list[str]:
    if x.shape != y.shape or x.dtype != y.dtype:
        return [f"{label}: shape/dtype differ: {x.shape}/{x.dtype} vs {y.shape}/{y.dtype}"]
    if not np.array_equal(x, y, equal_nan=True):
        return [f"{label}: array values differ (NaN positions are compared equal)"]
    return []


def _compare_csv(a: Path, b: Path, label: str) -> list[str]:
    try:
        with a.open(newline="", encoding="utf-8") as fa, b.open(newline="", encoding="utf-8") as fb:
            equal = list(csv.reader(fa)) == list(csv.reader(fb))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ComparisonError(f"{label}: CSV read failed: {exc}") from exc
    return [] if equal else [f"{label}: CSV cells differ"]


def _safe_artifact(root: Path, raw: Any) -> Path:
    if not isinstance(raw, str) or not raw:
        raise ComparisonError("manifest artifact path must be a nonempty string")
    root_resolved = root.resolve()
    p = (root / raw).resolve()
    if root_resolved not in p.parents and p != root_resolved:
        raise ComparisonError(f"manifest artifact escapes bundle: {raw}")
    if not p.is_file():
        raise ComparisonError(f"manifest artifact missing: {raw}")
    return p


def _manifest_pairs(reference: Path, candidate: Path, expected_source_ids: tuple[str, str] | None) -> list[tuple[str, Path, Path]]:
    rm, cm = reference / "manifest.json", candidate / "manifest.json"
    if not rm.is_file() or not cm.is_file():
        raise ComparisonError("directory inputs require manifest.json on both sides")
    r, c = _json(rm), _json(cm)
    if (not isinstance(expected_source_ids, (tuple, list)) or len(expected_source_ids) != 2
            or not all(isinstance(item, str) and item for item in expected_source_ids)):
        raise ComparisonError("bundle comparison requires explicit reference and candidate source identities")
    for side, obj in (("reference", r), ("candidate", c)):
        if not isinstance(obj, dict):
            raise ComparisonError(f"{side} manifest must be a JSON object")
        if not isinstance(obj.get("schema"), str) or obj["schema"] not in SUPPORTED_SCHEMAS:
            raise ComparisonError(f"{side} manifest schema is missing or unsupported")
        if not isinstance(obj.get("provenance"), dict):
            raise ComparisonError(f"{side} manifest requires a provenance object")
        if not isinstance(obj.get("artifacts"), dict) or not obj["artifacts"]:
            raise ComparisonError(f"{side} manifest requires nonempty artifacts mapping")
        provenance = obj["provenance"]
        if (not isinstance(provenance.get("source_identity"), str) or not provenance["source_identity"]
                or not isinstance(provenance.get("input_fingerprint"), str) or not provenance["input_fingerprint"]):
            raise ComparisonError(f"{side} manifest provenance IDs must be nonempty strings")
        if "run_metadata" in obj and not isinstance(obj["run_metadata"], dict):
            raise ComparisonError(f"{side} manifest run_metadata must be an object")
        expected = expected_source_ids[0 if side == "reference" else 1]
        if provenance["source_identity"] != expected:
            raise ComparisonError(f"{side} source identity does not match declared comparison contract")
    if r["provenance"]["input_fingerprint"] != c["provenance"]["input_fingerprint"]:
        raise ComparisonError("reference and candidate input fingerprints differ")
    for side, artifacts in (("reference", r["artifacts"]), ("candidate", c["artifacts"])):
        if any(not isinstance(key, str) or not key for key in artifacts):
            raise ComparisonError(f"{side} manifest artifact logical names must be nonempty strings")
        if "manifest.json" in artifacts:
            raise ComparisonError("artifact logical name manifest.json is reserved")
    if set(r["artifacts"]) != set(c["artifacts"]):
        raise ComparisonError("manifest artifact keys differ")
    pairs = [("manifest.json", rm, cm)]
    for key in sorted(r["artifacts"]):
        rv, cv = r["artifacts"][key], c["artifacts"][key]
        if not isinstance(key, str) or not key:
            raise ComparisonError("manifest artifact logical names must be nonempty strings")
        if not isinstance(rv, dict) or not isinstance(cv, dict):
            raise ComparisonError(f"artifact {key} metadata must be an object on both manifests")
        for side, item in (("reference", rv), ("candidate", cv)):
            claimed = item.get("sha256")
            if (not isinstance(claimed, str) or len(claimed) != 64
                    or any(ch not in "0123456789abcdef" for ch in claimed.lower())):
                raise ComparisonError(f"artifact {key} {side} sha256 must be a 64-character hexadecimal string")
        ra, ca = _safe_artifact(reference, rv.get("path")), _safe_artifact(candidate, cv.get("path"))
        if _sha256(ra) != rv["sha256"] or _sha256(ca) != cv["sha256"]:
            raise ComparisonError(f"artifact {key} sha256 claim does not match bytes")
        pairs.append((key, ra, ca))
    return pairs


def compare(reference: Path, candidate: Path, *, expected_source_ids: tuple[str, str] | None = None) -> dict[str, Any]:
    reference, candidate = Path(reference), Path(candidate)
    self_check = reference.resolve() == candidate.resolve()
    if self_check:
        raise ComparisonError("reference and candidate are the same path: self-check cannot certify cross-version equivalence")
    if not reference.exists() or not candidate.exists():
        raise ComparisonError("reference and candidate paths must exist")
    if reference.is_dir() != candidate.is_dir():
        raise ComparisonError("reference and candidate must both be files or both be directories")
    bundle = reference.is_dir()
    pairs = _manifest_pairs(reference, candidate, expected_source_ids) if bundle else [(reference.name, reference, candidate)]
    failures: list[str] = []
    checks = []
    for label, a, b in pairs:
        suffix = a.suffix.lower()
        if suffix == ".json":
            if bundle and label == "manifest.json":
                if _normalise(_json(a), manifest=True) != _normalise(_json(b), manifest=True):
                    failures.append(f"{label}: manifest scientific metadata differ")
            else:
                failures.extend(_compare_json(a, b, label))
        elif suffix in {".npy", ".npz"}: failures.extend(_compare_array(a, b, label))
        elif suffix == ".csv": failures.extend(_compare_csv(a, b, label))
        else: raise ComparisonError(f"{label}: unsupported artifact format {suffix!r}")
        checks.append({"label": label, "reference_sha256": _sha256(a), "candidate_sha256": _sha256(b)})
    return {"status": "pass" if not failures else "fail", "self_check": self_check,
            "scope": "artifact_bundle" if bundle else "artifact_only",
            "reference": str(reference.resolve()), "candidate": str(candidate.resolve()),
            "ignored_manifest_run_metadata_fields": sorted(MANIFEST_RUN_METADATA_IGNORED_FIELDS),
            "checks": checks, "failures": failures}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference", type=Path); p.add_argument("candidate", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--reference-source-id")
    p.add_argument("--candidate-source-id")
    args = p.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        ids = (args.reference_source_id, args.candidate_source_id) if args.reference_source_id and args.candidate_source_id else None
        result = compare(args.reference, args.candidate, expected_source_ids=ids)
    except Exception as exc:
        result = {"status": "error", "failures": [f"{type(exc).__name__}: {exc}"], "self_check": False}
    (args.output_dir / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [f"E5F refactor verification: {result['status']}", f"reference: {result.get('reference', args.reference)}", f"candidate: {result.get('candidate', args.candidate)}"]
    if result.get("failures"): lines.append("failures:\n- " + "\n- ".join(result["failures"]))
    (args.output_dir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
