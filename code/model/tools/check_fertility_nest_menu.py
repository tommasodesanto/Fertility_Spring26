"""Pure arithmetic check of a contingent-plan GEV representation.

This is a specification diagnostic, not a model implementation or calibration.
The additional outcome-dependent subnests are NOT adopted economic assumptions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def logsum(q, scale, axis=-1):
    q = np.asarray(q, dtype=float)
    maximum = np.max(q, axis=axis, keepdims=True)
    finite = np.isfinite(maximum)
    delta = np.full_like(q, -np.inf)
    np.subtract(q, maximum, out=delta, where=np.broadcast_to(finite, q.shape))
    exponent = np.exp(delta / scale)
    total = exponent.sum(axis=axis, keepdims=True)
    probability = np.divide(exponent, total, out=np.zeros_like(exponent), where=total > 0)
    logged = np.zeros_like(total)
    np.log(total, out=logged, where=total > 0)
    return np.squeeze(np.where(finite, maximum + scale * logged, -np.inf), axis), probability


def check_case(q0, q1, pi, housing, fertility, cost):
    """Compare independent nested enumeration with the sequential recursion."""
    if not 0 <= pi <= 1 or not 0 < housing <= fertility:
        raise ValueError("Require valid conception and GEV scales")
    l0, p0 = logsum(q0, housing)
    l1, p1 = logsum(q1, housing)
    attempt = 0.0
    if pi < 1:
        attempt += (1 - pi) * l0
    if pi > 0:
        attempt += pi * (l1 - cost)
    v_seq, a_seq = logsum(np.stack([l0, attempt], axis=-1), fertility)

    # Independently enumerate every complete contingent housing plan.
    if pi == 0:
        attempt_nested, conditional0 = logsum(q0, housing)
        conditional1 = np.zeros_like(p1)
    elif pi == 1:
        attempt_nested, conditional1 = logsum(q1 - cost, housing)
        conditional0 = np.zeros_like(p0)
    else:
        plans = (1 - pi) * q0[..., :, None] + pi * (q1[..., None, :] - cost)
        if pi <= 0.5:
            inner, inner_probability = logsum(plans, pi * housing, axis=-1)
            attempt_nested, outer_probability = logsum(inner, (1 - pi) * housing)
            joint = outer_probability[..., :, None] * inner_probability
        else:
            inner, inner_probability = logsum(plans, (1 - pi) * housing, axis=-2)
            attempt_nested, outer_probability = logsum(inner, pi * housing)
            joint = inner_probability * outer_probability[..., None, :]
        conditional0 = joint.sum(axis=-1)
        conditional1 = joint.sum(axis=-2)
    v_nested, a_nested = logsum(np.stack([l0, attempt_nested], axis=-1), fertility)
    # Unconditional realized product/family probabilities drive the forward law.
    no_seq = (a_seq[..., 0] + (1 - pi) * a_seq[..., 1])[..., None] * p0
    yes_seq = pi * a_seq[..., 1, None] * p1
    no_nested = a_nested[..., 0, None] * p0 + (1 - pi) * a_nested[..., 1, None] * conditional0
    yes_nested = pi * a_nested[..., 1, None] * conditional1
    return {
        "value_error": float(np.max(abs(v_seq - v_nested))),
        "action_probability_error": float(np.max(abs(a_seq - a_nested))),
        "realized_family_product_error": float(max(np.max(abs(no_seq - no_nested)), np.max(abs(yes_seq - yes_nested)))),
        "total_probability_error": float(np.max(abs(no_nested.sum(axis=-1) + yes_nested.sum(axis=-1) - 1))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20260907)
    cases = []
    pis = [0.0, 1e-9, 0.01, 0.2, 0.499999, 0.5, 0.500001, 0.8, 0.99, 1 - 1e-9, 1.0]
    # Six alternatives: renter plus the five maintained owner products.
    for spread in [0.001, 0.1, 5.0, 1000.0]:
        for housing, fertility in [(0.005, 2.1681730392479377), (0.005, 1.7364706586958831), (0.5, 0.5)]:
            for pi in pis:
                for repeat in range(10):
                    q0 = rng.normal(size=6) * spread
                    q1 = rng.normal(size=6) * spread
                    q0[rng.random(6) < 0.3] = -np.inf
                    q1[rng.random(6) < 0.3] = -np.inf
                    q0[0] = 0.0
                    q1[-1] = spread * 0.2
                    cases.append(check_case(q0, q1, pi, housing, fertility, 0.3))
    errors = {key: max(case[key] for case in cases) for key in cases[0]}
    # Values at utility magnitudes up to thousands accumulate floating roundoff.
    assert errors["value_error"] < 2e-12, errors
    assert errors["action_probability_error"] < 2e-12, errors
    assert errors["realized_family_product_error"] < 2e-11, errors
    assert errors["total_probability_error"] < 2e-12, errors
    summary = {
        "status": "PASS",
        "scope": "Pure contingent-plan GEV arithmetic; no Bellman, equilibrium or objective run",
        "candidate_status": "Additional outcome-dependent subnests NOT adopted",
        "seed": 20260907,
        "cases": len(cases),
        "housing_products": 6,
        "maximum_errors": errors,
        "limitations": [
            "Structural plan-level GEV shocks differ from original reused housing shocks.",
            "Full source integration and event-cohort coupling are not certified.",
            "No empirical fit is measured by this diagnostic.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
