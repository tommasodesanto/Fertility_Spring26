#!/usr/bin/env python3
"""Aggregate saved E5F housing arrays without solving or loading raw data.

The collector is deliberately a local/remote-array utility.  It expects a
saved packet containing ``g_current`` and ``hR_pol``; policy arrays are never
reconstructed by solving.  The tenure axis is conventionally 0=renter and
1..n_house=owner rungs, matching the E5F saved-state contract.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import sys
from pathlib import Path

import numpy as np


AGE_COUNT = 17


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def _aggregate(g: np.ndarray, h_r: np.ndarray, h_own: np.ndarray) -> dict:
    if g.ndim != 7 or h_r.shape != g.shape:
        raise ValueError(f"expected g/hR shape [wealth,tenure,location,age,income,n,cs], got {g.shape} and {h_r.shape}")
    if g.shape[3] != AGE_COUNT:
        raise ValueError(f"expected 17 model age cells, got {g.shape[3]}")
    if g.shape[1] - 1 != len(h_own):
        raise ValueError(f"H_own length {len(h_own)} does not match owner tenure states {g.shape[1]-1}")
    if not (np.isfinite(g).all() and np.isfinite(h_own).all()):
        raise ValueError("nonfinite mass or housing policy")
    if np.any(g < -1e-12):
        raise ValueError("negative g_current mass")
    occupied = g[:, 0] > 0
    if not np.isfinite(h_r[:, 0][occupied]).all() or np.any(h_r[:, 0][occupied] < 0):
        raise ValueError("negative occupied renter policy")

    # Cap every state before weighting and age aggregation.  Owners use the
    # pinned physical rung, while renters use the saved hR policy.
    room = np.zeros_like(g)
    room[:, 0] = np.where(occupied, np.minimum(h_r[:, 0], 9.0), 0.0)
    for ten, h in enumerate(h_own, start=1):
        room[:, ten, ...] = min(float(h), 9.0)
    age_mass = g.sum(axis=(0, 1, 2, 4, 5, 6))
    owner_mass = g[:, 1:, ...].sum(axis=(0, 1, 2, 4, 5, 6))
    room_mass = (g * room).sum(axis=(0, 1, 2, 4, 5, 6))
    renter_room_mass = (g[:, 0:1, ...] * room[:, 0:1, ...]).sum(axis=(0, 1, 2, 4, 5, 6))
    owner_room_mass = (g[:, 1:, ...] * room[:, 1:, ...]).sum(axis=(0, 1, 2, 4, 5, 6))
    def series(m, denominator=age_mass):
        return [None if x == 0 else float(y / x) for x, y in zip(denominator, m)]
    renter_mass = age_mass - owner_mass
    raw_renter_sum = float(np.sum(g[:, 0] * np.where(occupied, h_r[:, 0], 0.0)))
    raw_owner_sum = sum(float(g[:, ten].sum()) * float(h) for ten, h in enumerate(h_own, 1))
    ages = 18.0 + 4.0 * np.arange(AGE_COUNT)
    overlap = np.maximum(0.0, np.minimum(ages + 4.0, 56.0) - np.maximum(ages, 30.0)) / 4.0
    return {
        "age_cells": list(range(AGE_COUNT)),
        "total_mass_by_age": age_mass.astype(float).tolist(),
        "ownership_by_age": series(owner_mass),
        "mean_rooms_capped9_by_age": series(room_mass),
        "renter_mean_rooms_capped9_by_age": series(renter_room_mass, renter_mass),
        "owner_mean_rooms_capped9_by_age": series(owner_room_mass, owner_mass),
        "scalar_total_mass": float(g.sum()),
        "scalar_mean_rooms_capped9": float(room_mass.sum() / g.sum()),
        "scalar_ownership": float(g[:, 1:, ...].sum() / g.sum()),
        "scalar_mean_rooms_uncapped": (raw_renter_sum + raw_owner_sum) / float(g.sum()),
        "ownership_30_55_uniform_age_cells": float(overlap @ owner_mass / (overlap @ age_mass)) if overlap @ age_mass else None,
    }


def selftest() -> None:
    g = np.zeros((2, 2, 1, AGE_COUNT, 1, 1, 1))
    h = np.zeros_like(g)
    g[0, 0, 0, 0, 0, 0, 0] = 2.0; h[0, 0, 0, 0, 0, 0, 0] = 12.0
    g[1, 1, 0, 0, 0, 0, 0] = 1.0
    out = _aggregate(g, h, np.array([11.0]))
    assert out["scalar_total_mass"] == 3.0
    assert abs(out["scalar_mean_rooms_capped9"] - 9.0) < 1e-12
    assert abs(out["scalar_ownership"] - 1/3) < 1e-12
    assert out["renter_mean_rooms_capped9_by_age"][0] == 9.0
    assert out["owner_mean_rooms_capped9_by_age"][0] == 9.0
    h[0, 1, 0, 5, 0, 0, 0] = np.nan
    assert _aggregate(g, h, np.array([11.0]))["scalar_mean_rooms_capped9"] == 9.0


def collect(packet: Path, output: Path, family: str, h_own_arg: Path | None) -> None:
    with np.load(packet, allow_pickle=False) as z:
        keys = set(z.files)
        required = {"g_current", "hR_pol"}
        missing = required - keys
        if missing:
            raise FileNotFoundError(f"packet lacks {sorted(missing)}; keys={sorted(keys)}")
        g = np.asarray(z["g_current"], dtype=float)
        h_r = np.asarray(z["hR_pol"], dtype=float)
        h_own = np.asarray(z["H_own"], dtype=float) if "H_own" in keys else None
    if h_own is None and h_own_arg is not None:
        h_own = np.asarray(json.loads(h_own_arg.read_text())["H_own"], dtype=float)
    if h_own is None:
        raise FileNotFoundError("H_own absent from packet and --h-own was not supplied")
    result = {"family": family, "packet": str(packet), "packet_sha256": sha256(packet),
              "array_keys": sorted(keys), "g_current_shape": list(g.shape),
              "hR_pol_shape": list(h_r.shape), "H_own": h_own.tolist(),
              "occupied_policy_check": True, "child_at_home": "not reported: saved n/cs semantics not established",
              "profiles": _aggregate(g, h_r, h_own)}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")


def _get(obj, key):
    if isinstance(obj, dict):
        return obj[key]
    return getattr(obj, key)


def collect_checkpoint(checkpoint: Path, output: Path, family: str, source_root: Path, summary_path: Path) -> None:
    """Read the approved pinned checkpoint and emit only small sufficient stats."""
    if not source_root.is_dir():
        raise FileNotFoundError(f"frozen source root missing: {source_root}")
    summary = json.loads(summary_path.read_text())
    rows = [r for r in summary['cases'] if (r.get('phi'), r.get('lambda'), r.get('rental_cap')) == (0.8, 0.0, 6.0)]
    baseline = rows[0]
    if summary.get('status') != 'complete' or baseline.get('status') != 'completed':
        raise ValueError('baseline receipt is incomplete')
    contract = baseline['contract']
    if sha256(checkpoint) != contract['checkpoint_sha256']:
        raise ValueError('checkpoint hash mismatch')
    if (source_root / 'code/model').is_dir():
        source_root = source_root / 'code/model'
    sys.path[0:0] = [str(source_root / 'tools'), str(source_root)]
    with gzip.open(checkpoint, "rb") as stream:
        state = pickle.load(stream)
    evaluation = _get(state, "evaluation")
    policy = _get(evaluation, "policy")
    g = np.asarray(_get(evaluation, "g_current"), dtype=float)
    h_r = np.asarray(_get(policy, "hR_pol"), dtype=float)
    params = _get(state, "parameters")
    h_own = np.asarray(_get(params, "H_own"), dtype=float)
    if g.shape != h_r.shape or g.ndim != 7:
        raise ValueError(f"checkpoint g_current/hR_pol shape mismatch: {g.shape} {h_r.shape}")
    modules = {type(evaluation).__module__, type(policy).__module__, type(params).__module__}
    origins = {}
    for name in sorted(modules):
        mod = sys.modules.get(name)
        origin = getattr(mod, "__file__", None)
        origins[name] = origin
        if name in {'types', 'builtins'}:
            continue
        if origin and source_root.resolve() not in Path(origin).resolve().parents:
            raise RuntimeError(f"module origin outside frozen root: {name} -> {origin}")
    if (float(params.age_start), float(params.da), int(params.J)) != (18.0, 4.0, 17):
        raise ValueError('unexpected age grid')
    profiles = _aggregate(g, h_r, h_own)
    checks = {'mean_rooms': profiles['scalar_mean_rooms_uncapped'], 'ownership': profiles['scalar_ownership']}
    replay = {key: {'saved': baseline['metrics'][key], 'observed': value, 'gap': value-baseline['metrics'][key]} for key, value in checks.items()}
    if any(abs(row['gap']) > 1e-10 for row in replay.values()):
        raise ValueError(f'baseline scalar replay failed: {replay}')
    result = {"family": family, "checkpoint": str(checkpoint), "checkpoint_sha256": sha256(checkpoint),
              "source_root": str(source_root), "module_origins": origins,
              "g_current_shape": list(g.shape), "hR_pol_shape": list(h_r.shape),
              "H_own": h_own.tolist(), "child_at_home": "not reported: saved n/cs semantics not established",
              "profiles": profiles, "baseline_scalar_replay": replay,
              "summary_sha256": sha256(summary_path), "status": "completed"}
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--packet", type=Path)
    ap.add_argument("--checkpoint", type=Path)
    ap.add_argument("--source-root", type=Path)
    ap.add_argument("--summary", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--family", default="unknown")
    ap.add_argument("--h-own", type=Path)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest(); print("selftest: PASS"); return
    if not a.output:
        ap.error("--output is required unless --selftest")
    if a.checkpoint:
        if not a.source_root or not a.summary:
            ap.error("--source-root and --summary are required with --checkpoint")
        collect_checkpoint(a.checkpoint, a.output, a.family, a.source_root, a.summary)
    elif a.packet:
        collect(a.packet, a.output, a.family, a.h_own)
    else:
        ap.error("provide --checkpoint or --packet")


if __name__ == "__main__":
    main()
