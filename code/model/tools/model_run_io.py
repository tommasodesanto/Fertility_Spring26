"""Validated, pickle-free storage for local fixed-price model runs."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import uuid

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RUNS_ROOT = ROOT / "tmp/model_runs"
FORMAT_VERSION = 1


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _encode(value, arrays: dict[str, np.ndarray], *, label: str):
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            raise TypeError(f"Cannot safely serialize object array {label}")
        key = f"array_{len(arrays):05d}"
        arrays[key] = np.ascontiguousarray(value)
        return {"__ndarray__": key}
    if isinstance(value, np.generic):
        return _encode(value.item(), arrays, label=label)
    if isinstance(value, Path):
        return {"__path__": str(value)}
    if hasattr(value, "__dict__"):
        return {"__namespace__": _encode(vars(value), arrays, label=label)}
    if isinstance(value, float) and not np.isfinite(value):
        return {"__float__": "nan" if np.isnan(value) else ("inf" if value > 0 else "-inf")}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise TypeError(f"Only string-keyed mappings are supported ({label})")
        return {key: _encode(item, arrays, label=f"{label}.{key}") for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        tag = "__tuple__" if isinstance(value, tuple) else "__list__"
        return {tag: [_encode(item, arrays, label=f"{label}[{index}]")
                      for index, item in enumerate(value)]}
    raise TypeError(f"Unsupported value in native result at {label}: {type(value).__name__}")


def _decode(value, arrays):
    if isinstance(value, list):
        return [_decode(item, arrays) for item in value]
    if isinstance(value, dict):
        if set(value) == {"__ndarray__"}:
            return arrays[value["__ndarray__"]].copy()
        if set(value) == {"__path__"}:
            return Path(value["__path__"])
        if set(value) == {"__namespace__"}:
            from types import SimpleNamespace
            return SimpleNamespace(**_decode(value["__namespace__"], arrays))
        if set(value) == {"__float__"}:
            return {"nan": float("nan"), "inf": float("inf"), "-inf": -float("inf")}[
                value["__float__"]
            ]
        if set(value) == {"__tuple__"}:
            return tuple(_decode(item, arrays) for item in value["__tuple__"])
        if set(value) == {"__list__"}:
            return [_decode(item, arrays) for item in value["__list__"]]
        return {key: _decode(item, arrays) for key, item in value.items()}
    return value


def _namespace_values(value, label: str) -> dict:
    if not hasattr(value, "__dict__"):
        raise TypeError(f"Expected a native object with fields for {label}")
    return vars(value)


def _array_inventory(values: dict, prefix: str) -> dict[str, np.ndarray]:
    found = {}

    def visit(item, name):
        if isinstance(item, np.ndarray):
            found[name] = item
        elif hasattr(item, "__dict__"):
            visit(vars(item), name)
        elif isinstance(item, dict):
            for key, child in item.items():
                visit(child, f"{name}.{key}")
        elif isinstance(item, (list, tuple)):
            for index, child in enumerate(item):
                visit(child, f"{name}[{index}]")

    visit(values, prefix)
    return found


def _assert_equal_tree(left, right, label: str) -> None:
    if isinstance(left, np.ndarray):
        if not isinstance(right, np.ndarray) or left.dtype != right.dtype or left.shape != right.shape:
            raise RuntimeError(f"Round-trip array metadata changed: {label}")
        equal = (np.array_equal(left, right, equal_nan=True)
                 if left.dtype.kind in "fc" else np.array_equal(left, right))
        if not equal:
            raise RuntimeError(f"Round-trip array values changed: {label}")
        return
    if hasattr(left, "__dict__"):
        if not hasattr(right, "__dict__"):
            raise RuntimeError(f"Round-trip object type changed: {label}")
        _assert_equal_tree(vars(left), vars(right), label)
        return
    if isinstance(left, dict):
        if not isinstance(right, dict) or set(left) != set(right):
            raise RuntimeError(f"Round-trip mapping changed: {label}")
        for key in left:
            _assert_equal_tree(left[key], right[key], f"{label}.{key}")
        return
    if isinstance(left, (list, tuple)):
        if type(left) is not type(right) or len(left) != len(right):
            raise RuntimeError(f"Round-trip sequence changed: {label}")
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            _assert_equal_tree(left_item, right_item, f"{label}[{index}]")
        return
    if isinstance(left, float) and np.isnan(left):
        if not (isinstance(right, float) and np.isnan(right)):
            raise RuntimeError(f"Round-trip scalar changed: {label}")
        return
    if left != right:
        raise RuntimeError(f"Round-trip scalar changed: {label}")


def _validate_roundtrip(original, loaded) -> dict:
    if original.parameters != loaded.parameters or original.price != loaded.price:
        raise RuntimeError("Saved parameters or fixed price changed during round-trip")
    original_groups = {
        "solution": _namespace_values(original.solution, "solution"),
        "P": _namespace_values(original.P, "P"),
    }
    loaded_groups = {
        "solution": _namespace_values(loaded.solution, "solution"),
        "P": _namespace_values(loaded.P, "P"),
    }
    counts = {}
    for group, before in original_groups.items():
        after = loaded_groups[group]
        if set(before) != set(after):
            raise RuntimeError(f"Round-trip field inventory changed for {group}")
        _assert_equal_tree(before, after, group)
        before_arrays = _array_inventory(before, group)
        after_arrays = _array_inventory(after, group)
        if set(before_arrays) != set(after_arrays):
            raise RuntimeError(f"Round-trip array inventory changed for {group}")
        for key, left in before_arrays.items():
            right = after_arrays[key]
            if left.dtype != right.dtype or left.shape != right.shape:
                raise RuntimeError(f"Round-trip array metadata changed: {key}")
            equal = (np.array_equal(left, right, equal_nan=True)
                     if left.dtype.kind in "fc" else np.array_equal(left, right))
            if not equal:
                raise RuntimeError(f"Round-trip array values changed: {key}")
        counts[group] = {"fields": len(before), "arrays": len(before_arrays)}
    return counts


def _load_directory(run_directory: Path, *, allow_unvalidated: bool = False):
    run_directory = Path(run_directory).resolve()
    metadata_path = run_directory / "metadata.json"
    archive_path = run_directory / "native_result.npz"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    valid_status = metadata.get("status") == "complete" and "roundtrip_validation" in metadata
    private_status = allow_unvalidated and metadata.get("status") == "validating"
    if metadata.get("format_version") != FORMAT_VERSION or not (valid_status or private_status):
        raise RuntimeError("Run is incomplete or uses an unsupported format")
    if _sha256(archive_path) != metadata.get("native_result_sha256"):
        raise RuntimeError("Native result archive SHA-256 mismatch")
    with np.load(archive_path, allow_pickle=False) as archive:
        encoded = metadata["serialized"]
        decoded = _decode(encoded, archive)
    from types import SimpleNamespace
    tools_directory = str(Path(__file__).resolve().parent)
    if tools_directory not in sys.path:
        sys.path.insert(0, tools_directory)
    from model_playground import ModelResult
    solution = SimpleNamespace(**decoded["solution"])
    P = SimpleNamespace(**decoded["P"])
    result = ModelResult(solution, P, metadata["parameters"], metadata["fixed_price"],
                         label=metadata.get("label", "saved fixed-price run"))
    return result, metadata


def reserve_run_directory(runs_root: Path | None = None) -> Path:
    """Create a unique directory before solving so diagnostics can be written there."""
    runs_root = RUNS_ROOT if runs_root is None else Path(runs_root)
    runs_root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_directory = runs_root / f"{stamp}_{uuid.uuid4().hex[:8]}"
    run_directory.mkdir(parents=False, exist_ok=False)
    return run_directory


def save_run(result, *, run_metadata: dict, update_latest: bool = True,
             runs_root: Path | None = None, run_directory: Path | None = None,
             metadata_finalizer=None) -> Path:
    """Write, reopen, and validate a complete result before advancing latest.json."""
    runs_root = RUNS_ROOT if runs_root is None else Path(runs_root)
    runs_root.mkdir(parents=True, exist_ok=True)
    if run_directory is None:
        run_directory = reserve_run_directory(runs_root)
    else:
        run_directory = Path(run_directory).resolve()
        if not run_directory.is_dir() or run_directory.parent != runs_root.resolve():
            raise ValueError("run_directory must be a reserved directory directly under runs_root")
        if (run_directory / "metadata.json").exists() or (run_directory / "native_result.npz").exists():
            raise FileExistsError("Refusing to overwrite an existing run result")
    run_id = run_directory.name

    arrays: dict[str, np.ndarray] = {}
    serialized = {
        "solution": _encode(_namespace_values(result.solution, "solution"), arrays, label="solution"),
        "P": _encode(_namespace_values(result.P, "P"), arrays, label="P"),
    }
    archive_path = run_directory / "native_result.npz"
    temporary_archive = run_directory / ".native_result.npz.tmp"
    with temporary_archive.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary_archive, archive_path)

    metadata = dict(run_metadata)
    metadata.update(
        format_version=FORMAT_VERSION,
        status="validating",
        run_id=run_id,
        created_utc=datetime.now(timezone.utc).isoformat(),
        fixed_price=float(result.price),
        label=str(result.label),
        parameters={key: float(value) for key, value in result.parameters.items()},
        native_result_sha256=_sha256(archive_path),
        native_array_count=len(arrays),
        serialized=serialized,
    )
    metadata_path = run_directory / "metadata.json"
    temporary_metadata = run_directory / ".metadata.json.tmp"
    temporary_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True,
                                             allow_nan=False) + "\n", encoding="utf-8")
    with temporary_metadata.open("rb") as stream:
        os.fsync(stream.fileno())
    os.replace(temporary_metadata, metadata_path)

    loaded, loaded_metadata = _load_directory(run_directory, allow_unvalidated=True)
    validation = _validate_roundtrip(result, loaded)
    loaded_metadata["roundtrip_validation"] = validation
    if metadata_finalizer is not None:
        finalized_fields = metadata_finalizer()
        if not isinstance(finalized_fields, dict):
            raise TypeError("metadata_finalizer must return a dictionary")
        loaded_metadata.update(finalized_fields)
    loaded_metadata["status"] = "complete"
    temporary_metadata.write_text(json.dumps(loaded_metadata, indent=2, sort_keys=True,
                                             allow_nan=False) + "\n", encoding="utf-8")
    with temporary_metadata.open("rb") as stream:
        os.fsync(stream.fileno())
    os.replace(temporary_metadata, metadata_path)

    if update_latest:
        pointer = {"run_id": run_id, "directory": run_id,
                   "metadata_sha256": _sha256(metadata_path),
                   "updated_utc": datetime.now(timezone.utc).isoformat()}
        pointer_tmp = runs_root / ".latest.json.tmp"
        pointer_tmp.write_text(json.dumps(pointer, indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")
        with pointer_tmp.open("rb") as stream:
            os.fsync(stream.fileno())
        os.replace(pointer_tmp, runs_root / "latest.json")
    return run_directory


def load_run(RUN_DIRECTORY=None):
    """Load a saved run without solving; prefer the canonical local GE cache."""
    production_root = ROOT / "output/model/local_solution"
    if RUN_DIRECTORY is None:
        if (production_root / "latest").exists():
            from production.storage import load_latest
            result, run_directory = load_latest(production_root)
            return result, run_directory
        pointer_path = RUNS_ROOT / "latest.json"
        pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
        run_directory = (RUNS_ROOT / pointer["directory"]).resolve()
        if RUNS_ROOT.resolve() not in run_directory.parents:
            raise RuntimeError("latest.json points outside the model-run directory")
        metadata_path = run_directory / "metadata.json"
        if _sha256(metadata_path) != pointer.get("metadata_sha256"):
            raise RuntimeError("latest.json metadata SHA-256 mismatch")
    else:
        run_directory = Path(RUN_DIRECTORY).expanduser().resolve()
        if production_root.resolve() in run_directory.parents or run_directory == (production_root / "latest").resolve():
            from production.storage import load_case
            return load_case(run_directory)
    result, metadata = _load_directory(run_directory)
    return result, run_directory
