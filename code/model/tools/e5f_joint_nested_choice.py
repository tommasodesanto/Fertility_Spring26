"""Pure experimental four-plan choice rules; no imports or edits to the model.

The last two axes always mean (tenure: rent/own, action: wait/attempt).
Joint nesting groups alternatives by tenure. The sequential rule is an
explicitly different random-utility control on the same four payoffs.
"""
from __future__ import annotations

import numpy as np


def logsum_prob(values, scale, axis=-1):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Choice scale must be finite and positive")
    if np.isnan(values).any() or np.isposinf(values).any():
        raise ValueError("Utilities may contain negative infinity, not NaN/+inf")
    maximum = np.max(values, axis=axis, keepdims=True)
    feasible = np.isfinite(maximum)
    centered = np.full_like(values, -np.inf)
    np.subtract(values, maximum, out=centered, where=np.broadcast_to(feasible, values.shape))
    exponent = np.exp(centered / scale)
    total = np.sum(exponent, axis=axis, keepdims=True)
    probabilities = np.divide(exponent, total, out=np.zeros_like(exponent), where=total > 0)
    logged = np.zeros_like(total)
    np.log(total, out=logged, where=total > 0)
    result = np.where(feasible, maximum + scale * logged, -np.inf)
    return np.squeeze(result, axis=axis), probabilities


def plan_values(no_birth, birth, conception, first_cost=0.0, available=True):
    """Values of success-contingent housing plans within a committed tenure."""
    q0 = np.asarray(no_birth, dtype=float)
    if q0.shape[-1] != 2 or not 0 <= conception <= 1:
        raise ValueError("Require two tenures and conception probability in [0,1]")
    if not np.isfinite(first_cost):
        raise ValueError("Cost must be finite")
    q0 = np.where(q0 > -1e9, q0, -np.inf)
    result = np.full(q0.shape + (2,), -np.inf)
    result[..., 0] = q0
    if not available:
        return result
    q1 = np.asarray(birth, dtype=float)
    if q1.shape != q0.shape:
        raise ValueError("Birth/no-birth dimensions differ")
    q1 = np.where(q1 > -1e9, q1, -np.inf)
    # Never let an impossible zero-probability outcome eliminate a plan.
    attempt = np.zeros_like(q0)
    if conception < 1:
        attempt += (1 - conception) * q0
    if conception > 0:
        attempt += conception * (q1 - first_cost)
    result[..., 1] = attempt
    return result


def choose(plans, outer_scale, dissimilarity, rule="joint"):
    plans = np.asarray(plans, dtype=float)
    if plans.shape[-2:] != (2, 2):
        raise ValueError("Last axes must be tenure, action")
    if not np.isfinite(dissimilarity) or not 0 < dissimilarity <= 1:
        raise ValueError("Nested-GEV dissimilarity must lie in (0,1]")
    inner = float(outer_scale) * float(dissimilarity)
    if rule == "joint":
        inclusive, action_given_tenure = logsum_prob(plans, inner, axis=-1)
        value, tenure = logsum_prob(inclusive, outer_scale, axis=-1)
        probabilities = tenure[..., :, None] * action_given_tenure
    elif rule == "sequential":
        inclusive, tenure_given_action = logsum_prob(plans, outer_scale, axis=-2)
        value, action = logsum_prob(inclusive, inner, axis=-1)
        probabilities = action[..., None, :] * tenure_given_action
    else:
        raise ValueError("Unknown experimental rule")
    return value, probabilities


def scatter_joint_block(mass, probabilities, choices_no_birth, choices_birth,
                        conception, destination_family, origin_family,
                        transaction_index, transaction_weight):
    """Scatter an origin's four plans directly, preserving joint selection.

    mass: (wealth, inherited product); choices: (..., tenure); maps have axes
    (inherited product, destination product, parity, child count, wealth).
    Returns current (wealth, product, parity, child count), post-birth mass
    before transactions, and explicit births. No age/income transitions occur.
    """
    nb, nt = mass.shape
    npar, ncs = transaction_index.shape[2:4]
    current = np.zeros((nb, nt, npar, ncs))
    post = np.zeros_like(current)
    born = 0.0
    for action in range(2):
        outcomes = [(1.0, origin_family, choices_no_birth)] if action == 0 else [
            (1 - conception, origin_family, choices_no_birth),
            (conception, destination_family, choices_birth),
        ]
        for chance, family, choices in outcomes:
            if chance == 0 or family is None:
                continue
            nn, cs = family
            for tenure in range(2):
                selected_mass = mass * probabilities[..., tenure, action] * chance
                post[:, :, nn, cs] += selected_mass
                if action == 1 and family != origin_family:
                    born += float(selected_mass.sum())
                for old in range(nt):
                    for new in range(nt):
                        weights = selected_mass[:, old] * (choices[:, old, tenure] == new)
                        if not np.any(weights):
                            continue
                        idx = transaction_index[old, new, nn, cs]
                        wt = transaction_weight[old, new, nn, cs]
                        np.add.at(current[:, new, nn, cs], idx, weights * (1 - wt))
                        np.add.at(current[:, new, nn, cs], idx + 1, weights * wt)
    return current, post, born
