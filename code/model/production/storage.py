"""Pickle-free storage for completed local production GE cases.

The public ``latest`` pointer is advanced only after the native arrays reopen
cleanly and every cached artifact has been written.  Failed attempts therefore
cannot hide the preceding completed case.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace
import uuid

import numpy as np


FORMAT_VERSION = 1


class StoredResult:
    """Cached GE result with the convenience API expected by old plotters."""

    def __init__(self, solution, P, b_grid, price, *, parameters=None, label="stationary GE"):
        self.solution, self.P = solution, P
        self.b_grid = np.asarray(b_grid).copy()
        self.price = float(price)
        self.parameters = dict(parameters or {})
        self.label = str(label)

    def __getattr__(self, name):
        return getattr(self.solution, name)

    @staticmethod
    def _policy_tools():
        import sys
        model_tools = Path(__file__).resolve().parents[1] / "tools"
        if str(model_tools) not in sys.path:
            sys.path.insert(0, str(model_tools))
        import model_policy_tools
        return model_policy_tools

    def aggregates(self):
        aggregate_solution = self._policy_tools().aggregate_solution
        return aggregate_solution(self.solution, houses=np.asarray(self.P.H_own, dtype=float),
                                  age_start=int(self.P.age_start), period_years=float(self.P.period_years))

    def plot_policy(self, **kwargs):
        plot_policy = self._policy_tools().plot_policy
        return plot_policy(self.solution, houses=np.asarray(self.P.H_own, dtype=float),
                           age_start=int(self.P.age_start), period_years=float(self.P.period_years), **kwargs)

    def plot_aggregates(self, *, wealth_range="central"):
        plot_aggregates = self._policy_tools().plot_aggregates
        return plot_aggregates(self.solution, houses=np.asarray(self.P.H_own, dtype=float),
                               age_start=int(self.P.age_start), period_years=float(self.P.period_years),
                               wealth_range=wealth_range)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _encode(value, arrays: dict[str, np.ndarray], label: str):
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            raise TypeError(f"Object array cannot be stored safely: {label}")
        key = f"array_{len(arrays):05d}"
        arrays[key] = np.ascontiguousarray(value)
        return {"__array__": key}
    if isinstance(value, np.generic):
        return _encode(value.item(), arrays, label)
    if hasattr(value, "__dict__"):
        return {"__namespace__": _encode(vars(value), arrays, label)}
    if isinstance(value, Path):
        return {"__path__": str(value)}
    if isinstance(value, float) and not np.isfinite(value):
        return {"__float__": "nan" if np.isnan(value) else ("inf" if value > 0 else "-inf")}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise TypeError(f"Only string mapping keys are supported: {label}")
        return {key: _encode(item, arrays, f"{label}.{key}") for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return {"__tuple__" if isinstance(value, tuple) else "__list__":
                [_encode(item, arrays, f"{label}[]") for item in value]}
    raise TypeError(f"Unsupported native value {type(value).__name__} at {label}")


def _decode(value, arrays):
    if isinstance(value, dict):
        if set(value) == {"__array__"}: return arrays[value["__array__"]].copy()
        if set(value) == {"__namespace__"}: return SimpleNamespace(**_decode(value["__namespace__"], arrays))
        if set(value) == {"__path__"}: return Path(value["__path__"])
        if set(value) == {"__float__"}: return {"nan": float("nan"), "inf": float("inf"), "-inf": -float("inf")}[value["__float__"]]
        if set(value) == {"__tuple__"}: return tuple(_decode(item, arrays) for item in value["__tuple__"])
        if set(value) == {"__list__"}: return [_decode(item, arrays) for item in value["__list__"]]
        return {key: _decode(item, arrays) for key, item in value.items()}
    if isinstance(value, list): return [_decode(item, arrays) for item in value]
    return value


def reserve_case(output_root: Path) -> Path:
    cases = Path(output_root) / "cases"
    cases.mkdir(parents=True, exist_ok=True)
    name = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "_" + uuid.uuid4().hex[:8]
    case = cases / name
    case.mkdir()
    return case


def save_case(result: StoredResult, case: Path, *, metadata: dict) -> Path:
    """Save and reopen a case; this intentionally does not publish ``latest``."""
    arrays: dict[str, np.ndarray] = {}
    serialized = {"solution": _encode(result.solution, arrays, "solution"),
                  "P": _encode(result.P, arrays, "P"),
                  "b_grid": _encode(result.b_grid, arrays, "b_grid")}
    archive = case / "native_result.npz"
    temporary = case / ".native_result.npz.tmp"
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays); stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, archive)
    payload = dict(metadata, format_version=FORMAT_VERSION, status="complete",
                   created_utc=datetime.now(timezone.utc).isoformat(), price=result.price,
                   parameters=result.parameters, label=result.label, serialized=serialized,
                   native_result_sha256=_sha256(archive))
    metadata_path = case / "metadata.json"
    temporary = case / ".metadata.json.tmp"
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    with temporary.open("rb") as stream: os.fsync(stream.fileno())
    os.replace(temporary, metadata_path)
    reopened, _ = load_case(case)
    _assert_roundtrip(result, reopened)
    return case


def _assert_roundtrip(original: StoredResult, reopened: StoredResult) -> None:
    """Require exact recursive equality of every stored native field/value."""
    for label, left, right in (("solution", original.solution, reopened.solution),
                               ("P", original.P, reopened.P),
                               ("b_grid", original.b_grid, reopened.b_grid),
                               ("parameters", original.parameters, reopened.parameters),
                               ("price", original.price, reopened.price),
                               ("label", original.label, reopened.label)):
        _assert_equal(left, right, label)
    for field in ("closure", "report_directory"):
        if hasattr(original, field):
            _assert_equal(getattr(original, field), getattr(reopened, field), field)


def _assert_equal(left, right, label: str) -> None:
    if isinstance(left, np.generic):
        left = left.item()
    if isinstance(right, np.generic):
        right = right.item()
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
            raise RuntimeError(f"Round-trip native object changed: {label}")
        _assert_equal(vars(left), vars(right), label)
        return
    if isinstance(left, dict):
        if not isinstance(right, dict) or set(left) != set(right):
            raise RuntimeError(f"Round-trip mapping changed: {label}")
        for key, value in left.items():
            _assert_equal(value, right[key], f"{label}.{key}")
        return
    if isinstance(left, (list, tuple)):
        if type(left) is not type(right) or len(left) != len(right):
            raise RuntimeError(f"Round-trip sequence changed: {label}")
        for index, (a, b) in enumerate(zip(left, right)):
            _assert_equal(a, b, f"{label}[{index}]")
        return
    if isinstance(left, float) and np.isnan(left):
        if not (isinstance(right, float) and np.isnan(right)):
            raise RuntimeError(f"Round-trip scalar changed: {label}")
        return
    if left != right:
        raise RuntimeError(f"Round-trip scalar changed: {label}")


def load_case(case: Path) -> tuple[StoredResult, Path]:
    case = Path(case).resolve()
    metadata = json.loads((case / "metadata.json").read_text())
    archive = case / "native_result.npz"
    if metadata.get("format_version") != FORMAT_VERSION or metadata.get("status") != "complete":
        raise RuntimeError("Case is incomplete or uses an unsupported storage format")
    if _sha256(archive) != metadata.get("native_result_sha256"):
        raise RuntimeError("Native result archive SHA-256 mismatch")
    with np.load(archive, allow_pickle=False) as arrays:
        decoded = _decode(metadata["serialized"], arrays)
    result = StoredResult(decoded["solution"], decoded["P"], decoded["b_grid"], metadata["price"],
                          parameters=metadata.get("parameters"), label=metadata.get("label", "stationary GE"))
    result.closure = metadata.get("closure")
    result.report_directory = metadata.get("report_directory")
    return result, case


def publish_latest(case: Path, output_root: Path) -> Path:
    """Atomically replace the public latest symlink after a completed case."""
    root, case = Path(output_root).resolve(), Path(case).resolve()
    if root not in case.parents or case.parent.name != "cases":
        raise ValueError("latest can only point to a case directly under output_root/cases")
    load_case(case)
    latest = root / "latest"
    temporary = root / ".latest.tmp"
    temporary.unlink(missing_ok=True)
    temporary.symlink_to(Path("cases") / case.name)
    os.replace(temporary, latest)
    return latest


def load_latest(output_root: Path) -> tuple[StoredResult, Path]:
    latest = Path(output_root) / "latest"
    if not latest.exists():
        raise FileNotFoundError(f"No completed local GE case at {latest}")
    return load_case(latest)
