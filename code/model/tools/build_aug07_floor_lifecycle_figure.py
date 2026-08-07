#!/usr/bin/env python3
"""Build the Aug. 7 floor-arm lifecycle model-versus-data figure.

The baseline reconstruction follows ``run_e5_repaired_policy_with_entry``'s
E5F routing exactly.  The packet's nine free parameters are loaded from its
parameter table, converted to solver names/units, and checked against the raw
run metadata before the one permitted model solve.  No output artifact is
written unless every solved target moment reproduces the packet value.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import signal
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterator


# Set deterministic single-thread execution before importing NumPy/Numba.
for _name in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np  # noqa: E402


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

PACKET = ROOT / "output/model/eqscale_seq_e5_policy_quota_closure_20260807"
PARAMETER_TABLE = PACKET / "calibration/floor/parameter_table_full.csv"
TARGET_FIT = PACKET / "calibration/floor/target_fit_full.csv"
RAW_METADATA = PACKET / "raw_runs/floor/quota/metadata.json"
OUTDIR = PACKET / "diagnostics/floor/standard"
FIGURE = OUTDIR / "lifecycle_model_vs_data.png"
VERIFICATION = OUTDIR / "lifecycle_model_vs_data_verification.json"
SEED = 20260807
PERIOD_YEARS = 4.0
EXPECTED_PARAMETERS = {
    "beta_annual",
    "psi_child",
    "kappa_fert",
    "kappa_fert_continuation",
    "chi",
    "H0",
    "theta0",
    "theta1",
    "hbar_child_rooms",
}

# Provenance for the producing-chain reconstruction below. Line references are
# deliberately kept in the verification artifact as well as this source.
PROVENANCE = {
    "packet_theta": "output/model/eqscale_seq_e5_policy_quota_closure_20260807/calibration/floor/parameter_table_full.csv:1",
    "packet_raw_metadata": "output/model/eqscale_seq_e5_policy_quota_closure_20260807/raw_runs/floor/quota/metadata.json:3",
    "policy_profile_routing": "code/model/tools/run_e5_repaired_policy_with_entry.py:89",
    "chain_common_overrides": "code/model/intergen_eqscale_seq_optimized/run_e1_chain.py:289",
    "e5_externals": "code/model/intergen_eqscale_seq_optimized/e5_profile.py:96",
    "maturation_overrides": "code/model/intergen_eqscale_seq_optimized/e5_maturation_repair.py:17",
    "floor_overrides": "code/model/intergen_eqscale_seq_optimized/e5f_floor_profile.py:28",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-solve-seconds",
        type=float,
        default=900.0,
        help="Abort the one model solve after this many seconds (default: 900).",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=1.0e-6,
        help="Relative tolerance for every packet model-moment comparison.",
    )
    parser.add_argument(
        "--allow-overwrite",
        action="store_true",
        help="Replace only this driver's two output artifacts if they already exist.",
    )
    args = parser.parse_args()
    if not math.isfinite(args.max_solve_seconds) or args.max_solve_seconds <= 0.0:
        parser.error("--max-solve-seconds must be finite and positive")
    if not math.isfinite(args.rtol) or args.rtol <= 0.0:
        parser.error("--rtol must be finite and positive")
    return args


def load_floor_theta(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Load the nine packet estimates and convert annual beta to period beta."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    annual = {str(row["parameter"]): float(row["estimate"]) for row in rows}
    if len(rows) != len(annual) or set(annual) != EXPECTED_PARAMETERS:
        raise ValueError(
            f"Expected exactly the nine floor parameters in {path}; found {sorted(annual)}"
        )
    solver = dict(annual)
    solver["beta"] = solver.pop("beta_annual") ** PERIOD_YEARS
    return annual, solver


def configure_floor_profile() -> dict[str, str]:
    """Set the import-time gates used by the packet's E5F producing chain."""
    for name in (
        "E4_SPLIT",
        "E5",
        "E5_MATURATION_REPAIR",
        "E5F",
        "E5_PROBE_FIX_KE",
        "E6A",
        "E6B",
        "E6C",
    ):
        os.environ.pop(name, None)
    gates = {
        "E3_L4": "1",
        "E3_TFR_TOP_BIN_WEIGHT": "3.602359422009",
        "E5": "1",
        "E5_MATURATION_REPAIR": "1",
        "E5F": "1",
        "E5_PSI_MIN": "0.0",
        "E5_PSI_BOUND": "3.0",
    }
    os.environ.update(gates)
    return gates


def reconstruct_overrides(theta: dict[str, float]) -> tuple[dict[str, Any], dict[str, str]]:
    """Mirror the packet-producing policy driver's baseline construction."""
    gates = configure_floor_profile()
    from intergen_eqscale_seq_optimized import run_e1_chain

    run_e1_chain.load_runtime()
    runtime = SimpleNamespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5)
    overrides = run_e1_chain.common_overrides(runtime)
    overrides.update(theta)
    # These are the final baseline assignments in load_repaired_profile().
    overrides.update(
        child_state_mode="independent_count",
        property_tax_lump_sum_transfer=0.0,
    )
    return overrides, gates


def verify_packet_theta(overrides: dict[str, Any], metadata_path: Path) -> dict[str, float]:
    """Require the reconstructed solver values to equal the packet raw metadata."""
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    if payload.get("source_arm") != "E5F":
        raise ValueError(f"{metadata_path} is not the floor E5F arm")
    packet_theta = payload.get("theta")
    if not isinstance(packet_theta, dict):
        raise ValueError(f"{metadata_path} does not contain a theta dictionary")
    checked: dict[str, float] = {}
    mismatches: list[str] = []
    for name, expected_raw in packet_theta.items():
        expected = float(expected_raw)
        if name not in overrides:
            mismatches.append(f"{name}: absent from reconstructed overrides")
            continue
        actual = float(overrides[name])
        checked[str(name)] = actual
        if actual != expected:
            mismatches.append(f"{name}: reconstructed={actual:.17g}, packet={expected:.17g}")
    if mismatches:
        raise RuntimeError(
            "Floor theta/externals do not match packet raw metadata:\n" + "\n".join(mismatches)
        )
    return checked


@contextmanager
def solve_deadline(seconds: float) -> Iterator[None]:
    """Interrupt the single solve if it exceeds the requested wall-clock cap."""
    if not hasattr(signal, "SIGALRM"):
        yield
        return
    previous = signal.getsignal(signal.SIGALRM)

    def expired(_signum: int, _frame: Any) -> None:
        raise TimeoutError(f"The single model solve exceeded {seconds:.0f} seconds")

    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous)


def certified_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            {
                "moment": str(row["moment"]),
                "target_value": float(row["target"]),
                "packet_model_value": float(row["model"]),
            }
            for row in csv.DictReader(handle)
        ]
    names = [str(row["moment"]) for row in rows]
    if len(rows) != 12 or len(set(names)) != 12:
        raise ValueError(f"Expected 12 distinct packet target-fit rows in {path}")
    return rows


def compare_moments(
    moments: dict[str, Any], rows: list[dict[str, float | str]], rtol: float
) -> list[dict[str, Any]]:
    comparison: list[dict[str, Any]] = []
    for row in rows:
        name = str(row["moment"])
        if name not in moments:
            raise ValueError(f"Solved model does not report packet moment {name!r}")
        expected = float(row["packet_model_value"])
        solved = float(moments[name])
        absolute = abs(solved - expected)
        relative = absolute / max(abs(expected), 1.0e-300)
        comparison.append(
            {
                **row,
                "model_value": solved,
                "abs_deviation": absolute,
                "rel_deviation": relative,
                "passed": bool(absolute <= 1.0e-12 + rtol * abs(expected)),
            }
        )
    return comparison


def atomic_write_json(path: Path, payload: dict[str, Any]) -> Path:
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        return Path(handle.name)


def main() -> None:
    args = parse_args()
    existing = [path for path in (FIGURE, VERIFICATION) if path.exists()]
    if existing and not args.allow_overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing lifecycle artifact(s): "
            + ", ".join(str(path) for path in existing)
            + ". Pass --allow-overwrite only to replace this driver's own files."
        )

    random.seed(SEED)
    np.random.seed(SEED)
    annual_theta, solver_theta = load_floor_theta(PARAMETER_TABLE)
    overrides, gates = reconstruct_overrides(solver_theta)
    checked_packet_theta = verify_packet_theta(overrides, RAW_METADATA)

    from intergen_eqscale_seq_optimized.calibration import extract_moments
    from intergen_eqscale_seq_optimized.solver import run_model_cp_dt

    print("Starting the single deterministic floor-baseline solve.", flush=True)
    started = time.perf_counter()
    with solve_deadline(args.max_solve_seconds):
        solution, parameters, _price = run_model_cp_dt(overrides, verbose=False)
    solve_seconds = time.perf_counter() - started
    print(f"Single solve completed in {solve_seconds:.3f} seconds.", flush=True)

    moments = extract_moments(solution, parameters)
    comparison = compare_moments(moments, certified_rows(TARGET_FIT), args.rtol)
    for row in comparison:
        print(
            f"{row['moment']}: solved={row['model_value']:.15g} "
            f"packet={row['packet_model_value']:.15g} "
            f"abs={row['abs_deviation']:.3e} rel={row['rel_deviation']:.3e} "
            f"{'PASS' if row['passed'] else 'FAIL'}",
            flush=True,
        )
    failures = [row for row in comparison if not row["passed"]]
    if failures:
        mismatch_table = "\n".join(
            f"  {row['moment']}: abs={row['abs_deviation']:.6e}, "
            f"rel={row['rel_deviation']:.6e}"
            for row in failures
        )
        raise RuntimeError(
            "Packet moment verification failed; no figure or verification JSON was written.\n"
            + mismatch_table
            + "\nMost likely causes are a missing producing-chain override, an incorrect "
            "annual-to-period beta conversion, or a stale packet target-fit artifact."
        )

    max_abs = max(float(row["abs_deviation"]) for row in comparison)
    max_rel = max(float(row["rel_deviation"]) for row in comparison)
    # Import the established July 24 figure functions only after verification
    # passes, keeping the no-figure-on-mismatch contract explicit.
    from build_eqscale_note_draft_figures import (
        acs_age_profiles,
        model_age_profiles,
        plot_lifecycle,
    )

    model_profiles = model_age_profiles(solution, parameters)
    acs_profiles = acs_age_profiles()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    verification_payload = {
        "status": "pass",
        "relative_tolerance": float(args.rtol),
        "absolute_tolerance_floor": 1.0e-12,
        "max_abs_deviation": max_abs,
        "max_rel_deviation": max_rel,
        "solve_seconds": solve_seconds,
        "deterministic_seed": SEED,
        "parameter_table": str(PARAMETER_TABLE),
        "target_fit_table": str(TARGET_FIT),
        "figure": str(FIGURE),
        "override_reconstruction": {
            "named_baseline": (
                "intergen_eqscale_seq_optimized.run_e1_chain.common_overrides"
                "(SimpleNamespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5))"
            ),
            "environment_gates": gates,
            "applied_solver_theta": solver_theta,
            "annual_parameter_table_values": annual_theta,
            "final_explicit_assignments": {
                "child_state_mode": "independent_count",
                "property_tax_lump_sum_transfer": 0.0,
            },
            "packet_theta_metadata_values_checked": checked_packet_theta,
            "provenance": PROVENANCE,
        },
        "moment_comparison": comparison,
    }

    temporary_figure: Path | None = None
    temporary_json: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=OUTDIR, prefix=f".{FIGURE.name}.", suffix=".png", delete=False
        ) as handle:
            temporary_figure = Path(handle.name)
        plot_lifecycle(model_profiles, acs_profiles, temporary_figure)
        temporary_json = atomic_write_json(VERIFICATION, verification_payload)
        os.replace(temporary_figure, FIGURE)
        temporary_figure = None
        os.replace(temporary_json, VERIFICATION)
        temporary_json = None
    finally:
        for temporary in (temporary_figure, temporary_json):
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    print(f"Verification PASS: max_abs={max_abs:.3e}, max_rel={max_rel:.3e}", flush=True)
    print(f"Figure: {FIGURE}", flush=True)
    print(f"Verification JSON: {VERIFICATION}", flush=True)


if __name__ == "__main__":
    main()
