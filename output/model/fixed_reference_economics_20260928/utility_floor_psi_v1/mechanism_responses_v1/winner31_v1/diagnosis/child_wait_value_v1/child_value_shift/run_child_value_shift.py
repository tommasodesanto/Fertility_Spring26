#!/usr/bin/env python3
"""Saved-policy odds shifts from a fixed current first-child flow benefit."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np

HERE = Path(__file__).resolve().parent
WIN = HERE.parents[2]
BASE = WIN / "purchase_ltv_v1/local_run/retry5/results/baseline_80_80"
PERM = WIN / "purchase_ltv_v1/local_run/retry9/results/both_100_100"
PRE = BASE.parent / "q0_reference_inherited_states.npz"
BASE_PSI = 0.17156192800028292
FLOW_LEVELS = (0.01, 0.10, BASE_PSI, 0.25, 0.50)
EXPECTED_ZERO_SHIFT_PP = -0.051680757225407


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shifted_probability(p, shift):
    """Apply a finite log-odds shift, preserving exact saved 0/1 endpoints."""
    out = p.copy()
    interior = (p > 0.0) & (p < 1.0)
    z = np.log(p[interior]) - np.log1p(-p[interior]) + shift[interior]
    positive = z >= 0.0
    vals = np.empty_like(z)
    vals[positive] = 1.0 / (1.0 + np.exp(-z[positive]))
    ez = np.exp(z[~positive])
    vals[~positive] = ez / (1.0 + ez)
    out[interior] = vals
    return out


def main():
    meta = json.loads((HERE.parent / "summary.json").read_text())
    params = {r["parameter"]: float(r["estimate"])
              for r in csv.DictReader((BASE / "parameters.csv").open())}
    kappa = params["kappa_fert"]
    if not np.isclose(params["psi_child"], BASE_PSI, rtol=0.0, atol=1e-14):
        raise ValueError("Saved baseline psi_child differs from the pinned diagnostic value")
    pi = np.asarray(meta["fecundity_by_age"], dtype=float)[None, :, None]
    with np.load(PRE) as a:
        weights = a["g_pre"][:, 0, 0, :, :, 0, 0].copy()
    weights[:, 7:] = 0.0
    mass = float(weights.sum())
    if mass <= 0.0:
        raise ValueError("Young childless renter PRE mass is empty")

    arms = {}
    for name, folder in (("baseline_80_80", BASE), ("permanent_both_100_100", PERM)):
        with np.load(folder / "solution_arrays.npz") as a:
            probs = a["fert_probs"][:, 0, 0, :, :, :2]
        p = probs[..., 1].astype(float)
        if np.any((p < 0.0) | (p > 1.0)):
            raise ValueError(f"Invalid attempt probabilities in {name}")
        if p.shape != weights.shape:
            raise ValueError(f"Saved PRE/policy shape mismatch in {name}: {weights.shape} vs {p.shape}")
        arms[name] = p

    rows = []
    for flow in FLOW_LEVELS:
        delta = float(flow - BASE_PSI)
        shifted = {}
        aggregate = {}
        for name, p in arms.items():
            shift = np.broadcast_to(pi * delta / kappa, p.shape)
            pnew = shifted_probability(p, shift)
            shifted[name] = pnew
            aggregate[name] = weights * np.broadcast_to(pi, p.shape) * pnew
        diff = aggregate["permanent_both_100_100"] - aggregate["baseline_80_80"]
        pos_mass = float(np.maximum(diff, 0.0).sum())
        neg_mass = float(np.minimum(diff, 0.0).sum())
        net_mass = float(diff.sum())
        rows.append({
            "current_child_flow_value": float(flow),
            "delta_from_baseline_psi": delta,
            "baseline_psi_child": BASE_PSI,
            "baseline_firstbirth_probability_per_young_n0_renter_PRE": float(aggregate["baseline_80_80"].sum() / mass),
            "credit_firstbirth_probability_per_young_n0_renter_PRE": float(aggregate["permanent_both_100_100"].sum() / mass),
            "credit_minus_baseline_pp": float(net_mass / mass * 100.0),
            "positive_firstbirth_mass": pos_mass,
            "negative_firstbirth_mass": neg_mass,
            "net_firstbirth_mass": net_mass,
            "zero_probability_mass_baseline": float(weights[arms["baseline_80_80"] == 0.0].sum()),
            "zero_probability_mass_credit": float(weights[arms["permanent_both_100_100"] == 0.0].sum()),
            "one_probability_mass_baseline": float(weights[arms["baseline_80_80"] == 1.0].sum()),
            "one_probability_mass_credit": float(weights[arms["permanent_both_100_100"] == 1.0].sum()),
        })

    zero = next(r for r in rows if r["current_child_flow_value"] == BASE_PSI)
    zero_gap = zero["credit_minus_baseline_pp"] - EXPECTED_ZERO_SHIFT_PP
    if abs(zero_gap) > 1e-10:
        raise ValueError(f"delta=0 fails published response check: error={zero_gap:.12g} pp")

    with (HERE / "child_value_shift.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "saved_policy_current_child_value_shift_complete_no_model_solves",
        "interpretation": "Diagnostic shift in current successful-first-child value in both arms; continuation values/policies and all other model objects held fixed. Not psi recalibration, not GE, not adoption.",
        "flow_benefit_definition": "child_preferences.py: benefit[n,m]=psi_child*m^(1-child_benefit_curvature); m=1 gives psi_child exactly",
        "logit_shift_definition": "logit(p_try_new)=logit(p_try)+pi*delta/kappa_fert for interior probabilities",
        "fixed_baseline_psi_child": BASE_PSI,
        "kappa_fert": kappa,
        "young_n0_renter_PRE_mass": mass,
        "zero_shift_validation": {
            "expected_credit_minus_baseline_pp": EXPECTED_ZERO_SHIFT_PP,
            "computed_credit_minus_baseline_pp": zero["credit_minus_baseline_pp"],
            "absolute_error_pp": abs(zero_gap),
            "tolerance_pp": 1e-10,
        },
        "rows": len(rows),
        "new_model_solves": 0,
        "probability_endpoint_rule": "Exact saved p=0 and p=1 remain 0 and 1 under finite log-odds shifts; any underlying saved-probability underflow is retained, not repaired.",
        "sources": {
            "driver_sha256": sha256(Path(__file__)),
            "diagnostic_parent_summary_sha256": sha256(HERE.parent / "summary.json"),
            "baseline_solution_arrays_sha256": sha256(BASE / "solution_arrays.npz"),
            "credit_solution_arrays_sha256": sha256(PERM / "solution_arrays.npz"),
            "common_PRE_sha256": sha256(PRE),
            "baseline_parameters_sha256": sha256(BASE / "parameters.csv"),
        },
        "rows_data": rows,
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
