"""One-time input export from the authenticated checkpoint (verification side; local or Torch).

The only lab file that needs the old classes, and only to unpickle the
checkpoint. Paths are explicit; nothing is resolved from the manifest's
recorded Mac/Torch locations. Refuses an existing destination.

Writes
  <out>/inputs/bundle.json, arrays.npz    economic inputs (+ '_' artifacts, tagged)
  <out>/verification/reference_solution.npz   every ndarray attribute of the
      saved stationary solution plus stationary_g_pre (comparison only)
  <out>/export_receipt.json   pins; the lead copies bundle_json_sha256 into
      the launcher, so bundle edits cannot self-authenticate.

    python -m refactor_lab.export_inputs --reference-root R --checkpoint C \
        --oracle-code R/code/model --out DIR
"""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
import sys
from pathlib import Path

import numpy as np

from .inputs import (BUNDLE_SCHEMA, CHECKPOINT_SHA256, encode, is_artifact, load_inputs,
                     load_manifest, serialized, sha256_file)


def export(reference_root: Path, checkpoint: Path, oracle_code: Path, out: Path) -> dict:
    manifest = load_manifest(reference_root)
    if sha256_file(checkpoint) != CHECKPOINT_SHA256:
        raise RuntimeError("Checkpoint SHA differs from repeat_0212 pin")
    sys.path[:0] = [str(oracle_code / "tools"), str(oracle_code)]  # class resolution only
    with gzip.open(checkpoint, "rb") as stream:
        packet = pickle.load(stream)
    fields = vars(packet["parameters"])
    if serialized(fields) != manifest["actual_serialized_parameters"]:
        raise RuntimeError("Checkpoint parameters differ from manifest identity")
    arrays: dict[str, np.ndarray] = {}
    encoded = {k: encode(v, k, arrays) for k, v in sorted(fields.items())}
    price = np.asarray(packet["solution"].p_eq, dtype=float).reshape(-1)
    arrays["b_grid"] = np.asarray(packet["b_grid"], dtype=float)
    arrays["reference_price"] = price
    inputs_dir, verify_dir = out / "inputs", out / "verification"
    inputs_dir.mkdir(parents=True, exist_ok=False)
    verify_dir.mkdir()
    np.savez(inputs_dir / "arrays.npz", **arrays)
    meta = dict(schema=BUNDLE_SCHEMA, checkpoint_sha256=CHECKPOINT_SHA256,
                arrays_sha256=sha256_file(inputs_dir / "arrays.npz"), parameters=encoded,
                reference_price=price.tolist(),
                primitive_fields=sorted(k for k in fields if not is_artifact(k)),
                artifact_fields=sorted(k for k in fields if is_artifact(k)))
    (inputs_dir / "bundle.json").write_text(json.dumps(meta, indent=1, sort_keys=True, allow_nan=False) + "\n")
    solution = {k: v for k, v in vars(packet["solution"]).items()
                if isinstance(v, np.ndarray) and v.dtype != object}
    solution["stationary_g_pre"] = np.asarray(packet["stationary_g_pre"])
    np.savez(verify_dir / "reference_solution.npz", **solution)
    bundle_sha = sha256_file(inputs_dir / "bundle.json")
    loaded = load_inputs(inputs_dir, reference_root, bundle_sha, load_artifacts=True)
    receipt = dict(loaded.identity, checkpoint_path=str(checkpoint), bundle_json_sha256=bundle_sha,
                   reference_solution_sha256=sha256_file(verify_dir / "reference_solution.npz"),
                   reference_solution_arrays=sorted(solution),
                   skipped_solution_attributes=sorted(k for k, v in vars(packet["solution"]).items()
                                                      if k not in solution))
    (out / "export_receipt.json").write_text(json.dumps(receipt, indent=1, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-root", type=Path, required=True)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--oracle-code", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    print(json.dumps(export(a.reference_root, a.checkpoint, a.oracle_code, a.out), indent=1))


if __name__ == "__main__":
    main()
