"""Explicit, authenticated model inputs.

The package never builds parameters from constructor defaults. Inputs come
from a bundle exported once from the authenticated block0506 checkpoint
(`export_inputs.py`, run once on the pinned checkpoint). Loading re-applies the
manifest's serialization convention and requires exact equality with the
manifest's `actual_serialized_parameters` for every loaded field. The caller
supplies the bundle.json SHA from an independent receipt.

Bundle layout (directory):
  bundle.json     typed fields, reference price, primitive/artifact lists, arrays pin
  arrays.npz      every ndarray field plus `b_grid` and `reference_price`
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from . import LABEL

MANIFEST_RELATIVE = "output/model/fertility_identification_20260928/fixed_reference_manifest.json"
MANIFEST_SHA256 = "147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4"
CHECKPOINT_SHA256 = "b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d"
PSI_CHILD = 0.1355551166583114
BUNDLE_SCHEMA = "refactor_lab_inputs_v1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def serialized(value: Any) -> Any:
    """Manifest identity convention (verbatim from run_fixed_price.py)."""
    if hasattr(value, "shape") and getattr(value, "size", 0) > 1024:
        h = hashlib.sha256()
        for start in range(0, value.size, 65536):
            h.update(value.flat[start:start + 65536].tobytes())
        return dict(serialized_array=True, shape=list(value.shape), dtype=str(value.dtype),
                    size=int(value.size), sha256_c_order_bytes=h.hexdigest())
    if hasattr(value, "tolist"):
        return serialized(value.tolist())
    if isinstance(value, dict):
        return {str(k): serialized(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialized(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite_float": repr(value)}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError("Unsupported serialized value: " + str(type(value)))


# Typed encoding so arrays, tuples and numpy scalars survive the round trip.
def encode(value: Any, key: str, arrays: dict[str, np.ndarray]) -> Any:
    if isinstance(value, np.ndarray):
        name = f"a{len(arrays):04d}"
        arrays[name] = value
        return {"__ndarray__": name, "field": key}
    if isinstance(value, np.generic):
        return {"__npscalar__": value.dtype.str, "value": value.item()}
    if isinstance(value, tuple):
        return {"__tuple__": [encode(v, key, arrays) for v in value]}
    if isinstance(value, list):
        return [encode(v, key, arrays) for v in value]
    if isinstance(value, dict):
        if not all(isinstance(k, str) for k in value):
            raise TypeError(f"Non-string dict key in {key}")
        return {"__dict__": {k: encode(v, key, arrays) for k, v in value.items()}}
    if isinstance(value, float) and not math.isfinite(value):
        return {"__float__": repr(value)}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"Unsupported input type in {key}: {type(value)}")


def decode(value: Any, arrays: Any) -> Any:
    if isinstance(value, list):
        return [decode(v, arrays) for v in value]
    if isinstance(value, dict):
        if "__ndarray__" in value:
            return np.array(arrays[value["__ndarray__"]])
        if "__npscalar__" in value:
            return np.dtype(value["__npscalar__"]).type(value["value"])
        if "__tuple__" in value:
            return tuple(decode(v, arrays) for v in value["__tuple__"])
        if "__dict__" in value:
            return {k: decode(v, arrays) for k, v in value["__dict__"].items()}
        if "__float__" in value:
            return float(value["__float__"])
        raise ValueError("Untagged dict in bundle")
    return value


def is_artifact(name: str) -> bool:
    """Underscore fields are outputs of a previous solve stored on P.

    solve_bellman_full_markov_income overwrites _fert2_probs, _joint_choice,
    _bp_pol_stay and _c_pol_stay before the KFE reads them; the KFE overwrites
    _entry_*, _g_stay_distribution and the *_by_age flows. They are therefore
    verification artifacts, not economic inputs, and are omitted from the
    normal runtime P (equivalence is tested by the acceptance suite).
    """
    return name.startswith("_")


@dataclass(frozen=True)
class ModelInputs:
    """Everything the stationary calculation reads; nothing is defaulted."""

    parameters: SimpleNamespace   # primitive checkpoint fields (no '_' artifacts)
    b_grid: np.ndarray            # exact checkpoint wealth grid (Nb,)
    reference_price: np.ndarray   # checkpoint solution.p_eq (fixed-price replay / GE start)
    artifacts: dict               # '_' fields, only when load_artifacts=True
    identity: dict                # manifest/checkpoint/bundle pins


def load_manifest(root: Path) -> dict:
    path = Path(root) / MANIFEST_RELATIVE
    if sha256_file(path) != MANIFEST_SHA256:
        raise RuntimeError("Reference manifest SHA differs")
    manifest = json.loads(path.read_text())
    if manifest["label"] != LABEL or manifest["checkpoint"]["sha256"] != CHECKPOINT_SHA256:
        raise RuntimeError("Manifest does not name the block0506 repeat_0212 checkpoint")
    return manifest


def load_inputs(bundle_dir: Path, root: Path, bundle_sha256: str, *, load_artifacts: bool = False) -> ModelInputs:
    """Load a bundle pinned by an external bundle.json hash and prove identity.

    Chain: caller pin -> bundle.json -> arrays.npz. Primitive fields must equal
    the manifest's serialized identity; artifacts are checked only if loaded.
    """
    bundle_dir = Path(bundle_dir)
    manifest = load_manifest(root)
    if sha256_file(bundle_dir / "bundle.json") != bundle_sha256:
        raise RuntimeError("bundle.json differs from the external pin")
    meta = json.loads((bundle_dir / "bundle.json").read_text())
    if meta["schema"] != BUNDLE_SCHEMA or meta["checkpoint_sha256"] != CHECKPOINT_SHA256:
        raise RuntimeError("Bundle schema or checkpoint pin differs")
    if sha256_file(bundle_dir / "arrays.npz") != meta["arrays_sha256"]:
        raise RuntimeError("Bundle arrays differ from their recorded hash")
    expected = manifest["actual_serialized_parameters"]
    if set(meta["parameters"]) != set(expected):
        raise RuntimeError("Bundle field set differs from manifest: "
                           + str(sorted(set(meta["parameters"]) ^ set(expected))))
    wanted = [k for k in expected if load_artifacts or not is_artifact(k)]
    with np.load(bundle_dir / "arrays.npz", allow_pickle=False) as arrays:
        fields = {k: decode(meta["parameters"][k], arrays) for k in wanted}
        b_grid = np.array(arrays["b_grid"])
        price = np.array(arrays["reference_price"])
    mismatched = [k for k in wanted if serialized(fields[k]) != expected[k]]
    if mismatched:
        raise RuntimeError("Bundle fields differ from manifest identity: " + ", ".join(mismatched))
    if float(fields["psi_child"]) != PSI_CHILD:
        raise RuntimeError("psi_child differs from the frozen reference value")
    if not np.array_equal(b_grid, fields["fixed_reference_entry_grid"]):
        raise RuntimeError("Wealth grid differs from the fixed entry grid")
    if meta["reference_price"] != price.tolist():
        raise RuntimeError("Reference price differs from the pinned bundle record")
    primitives = {k: v for k, v in fields.items() if not is_artifact(k)}
    artifacts = {k: v for k, v in fields.items() if is_artifact(k)}
    identity = dict(label=LABEL, manifest_sha256=MANIFEST_SHA256, checkpoint_sha256=CHECKPOINT_SHA256,
                    bundle_json_sha256=bundle_sha256, arrays_sha256=meta["arrays_sha256"],
                    primitive_fields=len(primitives), artifact_fields_loaded=len(artifacts),
                    reference_price=price.tolist())
    return ModelInputs(SimpleNamespace(**primitives), b_grid, price, artifacts, identity)
