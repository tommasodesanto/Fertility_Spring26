#!/usr/bin/env python3
"""Build a compact report from the native fixed-price financing packets.

This is a reporting-only script.  It never solves the model.  Each case is
read one at a time from ``arrays.npz`` and is measured against the same saved
beginning-of-period distribution ``g_pre``.  The first-birth calculation is
the native sequential-transition calculation, not a rescaled attempt rate.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

CHECKPOINT_SHA256 = "3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993"

POLICY_NAMES = (
    "V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs",
    "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price",
)
CASE_ORDER = ("baseline", "mortgage_only", "unsecured_only", "both")


def _jsonable(x: Any) -> Any:
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def _load_metadata(checkpoint: Path, source_root: Path) -> tuple[dict[str, Any], str | None]:
    if not checkpoint.exists() or not source_root.exists():
        raise ValueError("checkpoint and frozen source root are required")
    try:
        if source_root and source_root.exists():
            sys.path.insert(0, str(source_root))
            sys.path.insert(0, str(source_root / "tools"))
        import hashlib
        h = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        if h != CHECKPOINT_SHA256:
            raise ValueError(f"checkpoint hash mismatch: {h}")
        with gzip.open(checkpoint, "rb") as f:
            packet = pickle.load(f)
        p = packet.get("parameters", packet) if isinstance(packet, dict) else packet
        b = packet.get("b_grid") if isinstance(packet, dict) else None
        if b is None:
            b = getattr(p, "b_grid", None)
        meta = {
            "b_grid": np.asarray(b, dtype=float),
            "age_start": float(p.age_start), "da": float(p.da),
            "period_years": float(p.period_years), "H_own": np.asarray(p.H_own, dtype=float),
            "g_pre_reference": np.asarray(packet["stationary_g_pre"]),
        }
        if meta["b_grid"].ndim != 1 or meta["b_grid"].size < 2 or meta["H_own"].ndim != 1:
            raise ValueError("checkpoint has invalid wealth grid or owner-room sizes")
        meta["checkpoint_sha256"] = h
        return meta, None
    except Exception as exc:
        raise RuntimeError(f"checkpoint metadata unavailable: {type(exc).__name__}: {exc}") from exc


def _case_dirs(input_root: Path) -> dict[str, Path]:
    found: dict[str, Path] = {}
    for p in sorted((input_root / "cases").glob("*/arrays.npz")):
        name = p.parent.name
        for case in CASE_ORDER:
            if name == case or name.startswith(case + "_"):
                found.setdefault(case, p.parent)
    return found


def _load_case(path: Path) -> dict[str, np.ndarray]:
    receipt = json.loads((path / "receipt.json").read_text())
    if receipt.get("status") != "completed":
        raise ValueError(f"case receipt is not completed: {receipt.get('status')}")
    actual = receipt.get("contract", {}).get("checkpoint_sha256")
    if actual != CHECKPOINT_SHA256:
        raise ValueError(f"case checkpoint hash mismatch: {actual}")
    with np.load(path / "arrays.npz", allow_pickle=False) as z:
        missing = [x for x in ("g_pre", "g_post_fertility", "g_current", "births", *POLICY_NAMES) if x not in z]
        if missing:
            raise ValueError(f"missing mandatory arrays: {missing}")
        return {x: np.asarray(z[x]) for x in z.files}


def _age(meta: dict[str, Any], j: int) -> float:
    return float(meta.get("age_start", 18.0) + j * meta.get("da", meta.get("period_years", 4.0)))


def _first_birth_flow(a: dict[str, np.ndarray], meta: dict[str, Any], age_lo: float | None = None, age_hi: float | None = None) -> float:
    """Exact n=0 loss at the fertility stage, before housing/current choice."""
    g, post = a["g_pre"], a["g_post_fertility"]
    if g.shape != post.shape or g.ndim != 7:
        raise ValueError(f"unexpected native fertility distributions: {g.shape}, {post.shape}")
    total = 0.0
    for j in range(g.shape[3]):
        age = _age(meta, j)
        if age_lo is not None and age < age_lo:
            continue
        if age_hi is not None and age > age_hi:
            continue
        loss = float(np.sum(g[:, :, :, j, :, 0, :]) - np.sum(post[:, :, :, j, :, 0, :]))
        if loss < -2e-11:
            raise ValueError(f"n=0 mass increases at age index {j}: {loss}")
        total += max(loss, 0.0)
    return total


def _topcode_births(a: dict[str, np.ndarray], meta: dict[str, Any]) -> tuple[float, float, float]:
    explicit = float(np.asarray(a["births"]).sum())
    return explicit, float("nan"), explicit


def _physical_housing(a: dict[str, np.ndarray], meta: dict[str, Any]) -> tuple[float, float, float]:
    gc = a["g_current"]
    hown = np.asarray(meta.get("H_own", []), dtype=float)
    total = float(gc.sum())
    if total <= 0 or hown.size != gc.shape[1] - 1:
        raise ValueError("owner-size metadata unavailable")
    owner_mass = float(gc[:, 1:, ...].sum())
    owner_rooms = float(sum(gc[:, t, ...].sum() * hown[t - 1] for t in range(1, gc.shape[1])))
    renter_rooms = float(np.sum(gc[:, 0, ...] * a["hR_pol"][:, 0, ...]))
    rooms = (owner_rooms + renter_rooms) / total
    return owner_mass / total, rooms, owner_rooms / max(owner_mass, 1e-15)


def _censoring(a: dict[str, np.ndarray], meta: dict[str, Any]) -> tuple[float, float]:
    gc = a["g_current"]
    realized_boundary = float(gc[0, ...].sum()) / max(float(gc.sum()), 1e-15)
    b = np.asarray(meta.get("b_grid", []), dtype=float)
    if b.size == 0:
        return realized_boundary, float("nan")
    bp = a["bp_pol"]
    bp_boundary = float(np.sum(gc * ((bp <= b[0] + 2e-11) | (bp >= b[-1] - 2e-11)))) / max(float(gc.sum()), 1e-15)
    return realized_boundary, bp_boundary


def _probability_postcheck(a: dict[str, np.ndarray]) -> dict[str, Any]:
    """Check normalized choice probabilities on occupied post-fertility states."""
    g = np.asarray(a["g_post_fertility"], dtype=float)
    checks: dict[str, Any] = {}
    for name, axis, aligned in (("tenure_probs", -1, g), ("loc_probs", 3, g)):
        raw = np.asarray(a[name])
        if not np.issubdtype(raw.dtype, np.floating):
            raise ValueError(f"probability array is not floating: {name}, {raw.dtype}")
        p = raw.astype(np.float64)
        if name == "loc_probs": expected = (g.shape[0], g.shape[1], g.shape[2], p.shape[3], *g.shape[3:])
        else: expected = aligned.shape + (p.shape[-1],)
        if p.shape != expected:
            checks[name] = {"status": "unavailable", "reason": f"shape {p.shape} does not align with g_post {g.shape}"}
            continue
        occ = aligned > 2e-11
        sums = p.sum(axis=axis)
        gap = np.abs(sums - 1.0)
        tolerance = 2e-11 + float(np.finfo(raw.dtype).eps)
        bad = occ & (~np.isfinite(sums) | (gap > tolerance))
        weighted = float(np.sum(aligned * gap) / max(float(np.sum(aligned)), 1e-300))
        checks[name] = {"status": "invalid" if np.any(bad) else "ok", "dtype": str(raw.dtype), "tolerance": tolerance, "occupied_cells": int(occ.sum()), "bad_cells": int(bad.sum()), "max_abs_gap": float(np.max(gap[occ])) if np.any(occ) else 0.0, "mass_weighted_abs_gap": weighted}
        if np.any(bad):
            raise ValueError(f"probability normalization failed: {name}, dtype={raw.dtype}, tolerance={tolerance:.17g}, max_gap={checks[name]['max_abs_gap']:.17g}, weighted_gap={weighted:.17g}, bad_cells={int(bad.sum())}")
    return checks


def _measure(case: str, a: dict[str, np.ndarray], meta: dict[str, Any]) -> dict[str, Any]:
    if not np.array_equal(a["g_pre"], meta["g_pre_reference"]):
        raise ValueError("case g_pre differs from frozen checkpoint stationary_g_pre")
    mass = float(a["g_pre"].sum())
    explicit, top_entry, adjusted = _topcode_births(a, meta)
    owner, rooms, owner_rooms = _physical_housing(a, meta)
    boundary, bp_boundary = _censoring(a, meta)
    probability_check = _probability_postcheck(a)
    first = _first_birth_flow(a, meta)
    if case == "baseline" and not np.isclose(first, 0.05064983385601245, rtol=0, atol=2e-11):
        raise ValueError(f"baseline first-birth flow mismatch: {first}")
    return {
        "case": case, "initial_households": mass, "births_per_initial_household": adjusted / mass,
        "explicit_births_per_initial_household": explicit / mass, "top_bin_entry_flow": top_entry / mass if np.isfinite(top_entry) else float("nan"),
        "first_birth_flow_initial_n0": first / mass,
        "first_birth_flow_age26_38_initial_n0": _first_birth_flow(a, meta, 26.0, 38.0) / mass,
        "ownership_rate_realized": owner, "mean_physical_rooms": rooms,
        "mean_owner_physical_rooms": owner_rooms, "saving_boundary_mass_realized": boundary,
        "saving_boundary_mass_bp_policy": bp_boundary,
        "probability_postcheck": json.dumps(probability_check, sort_keys=True),
    }


def _write_figures(rows: list[dict[str, Any]], out: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        (out / "figures_unavailable.txt").write_text(f"matplotlib unavailable: {exc}\n")
        return
    labels = [r["case"] for r in rows]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(2, 2, figsize=(9, 6), constrained_layout=True)
    ax[0, 0].bar(x, [r["first_birth_flow_age26_38_initial_n0"] for r in rows]); ax[0, 0].set_title("First birth flow, ages 26–38"); ax[0, 0].set_ylabel("per initial household")
    ax[0, 1].bar(x, [r["ownership_rate_realized"] for r in rows]); ax[0, 1].set_title("Realized ownership rate")
    ax[1, 0].bar(x, [r["mean_physical_rooms"] for r in rows]); ax[1, 0].set_title("Mean physical rooms")
    ax[1, 1].bar(x, [r["saving_boundary_mass_realized"] for r in rows]); ax[1, 1].set_title("Realized saving-boundary mass")
    for axes in ax.flat: axes.set_xticks(x, labels, rotation=25, ha="right")
    fig.savefig(out / "native_financing_supplement.png", dpi=160); fig.savefig(out / "native_financing_supplement.pdf"); plt.close(fig)


def build(args: argparse.Namespace) -> None:
    args.output.mkdir(parents=True, exist_ok=True)
    meta, meta_reason = _load_metadata(args.checkpoint, args.source_root)
    dirs = _case_dirs(args.input)
    rows, missing = [], {}
    for case in CASE_ORDER:
        if case not in dirs:
            missing[case] = "case arrays.npz not available"
            continue
        try:
            rows.append(_measure(case, _load_case(dirs[case]), meta))
        except Exception as exc:
            missing[case] = f"measurement unavailable: {type(exc).__name__}: {exc}"
            print(f"{case}: {missing[case]}", file=sys.stderr)
    (args.output / "metadata.json").write_text(json.dumps({"missing_cases": missing, "metadata_note": meta_reason, "checkpoint_sha256": meta.get("checkpoint_sha256"), "normalization_status": "fertility normalization not separately audited; tenure and location normalization audited", "cases": [r["case"] for r in rows]}, indent=2, default=_jsonable) + "\n")
    with (args.output / "comparisons.csv").open("w", newline="") as f:
        fields = sorted({k for r in rows for k in r} | {"case", "status", "reason"})
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n"); w.writeheader()
        for r in rows: w.writerow({**r, "status": "available", "reason": ""})
        for case, reason in missing.items(): w.writerow({"case": case, "status": "missing", "reason": reason})
    lines = ["# Native financing diagnostic", "", "Fixed-price partial-equilibrium comparison using the saved `g_pre` for every arm.", "", "First births equal the exact loss of `n=0` mass between `g_pre` and `g_post_fertility`, before housing and current-tenure choices. Total births are the native `births` scalar. The report does not infer an adjusted top-code total unless native accounting supplies it. Tenure and location probability checks sum in float64 and allow only the stored array dtype's single machine epsilon in addition to the reporting arithmetic tolerance of 2e-11; they never renormalize probabilities. Fertility normalization is not separately audited.", "", "## Comparison", "", "| Case | Births / initial HH | First births / initial HH | First births age 26–38 / initial HH | Ownership | Mean physical rooms | Saving boundary (realized) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append("| {case} | {births_per_initial_household:.6g} | {first_birth_flow_initial_n0:.6g} | {first_birth_flow_age26_38_initial_n0:.6g} | {ownership_rate_realized:.6g} | {mean_physical_rooms:.6g} | {saving_boundary_mass_realized:.6g} |".format(**r))
    for case, reason in missing.items():
        lines.append(f"| {case} | missing | missing | missing | missing | missing | missing |")
        lines.append(f"Reason for {case}: {reason}")
    lines += ["", "All grouped birth outcomes use initial states in `g_pre`; current-tenure redistribution is not reclassified as initial tenure. Renter rooms use realized renter mass and `hR_pol`; owner rooms use realized owner mass and `H_own`. Saving-boundary censoring uses realized current mass and both saving-grid boundaries.", "", "Supplementary figures: `native_financing_supplement.png` and `.pdf`."]
    if meta_reason: lines += ["", f"Metadata note: {meta_reason}."]
    (args.output / "report.md").write_text("\n".join(lines) + "\n")
    if rows: _write_figures(rows, args.output)
    if missing:
        raise RuntimeError(f"report incomplete; reasons saved in metadata.json: {missing}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--source-root", type=Path, required=True)
    build(p.parse_args(argv))


if __name__ == "__main__":
    main()
