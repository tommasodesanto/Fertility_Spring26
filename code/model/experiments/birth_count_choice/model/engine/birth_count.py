"""Direct intended-birth count menu and one-draw distribution transition.

Actions k and realized increments x each use the final axis, padded to length 4.
A state tensor has axes (..., children ever born n, children at home m).
Policy tensors have axes (..., n, m, k) or (..., n, m, x). No new Gumbel
scale or fertility cost is introduced. Birth-order tags count each crossed
order once: a 0->2 household appears in both first- and second-birth tags.
"""
from __future__ import annotations
import math
import json
import numpy as np
from .shared import DEAD_VALUE_CUTOFF, DEAD_MASS_TOL
from .parameters import independent_child_maturation_active, parent_age_maturation_active, readiness_gate_active


def enabled(P):
    return bool(getattr(P, 'birth_count_choice_enabled', False))


def validate_contract(P):
    if not enabled(P):
        return
    if int(P.n_parity) != 4 or int(P.n_child_states) != 4:
        raise ValueError('Birth-count experiment requires n=0..3 and m=0..3.')
    if not independent_child_maturation_active(P) or not bool(getattr(P, 'sequential_births', False)):
        raise ValueError('Birth-count experiment requires the retained independent-count child architecture.')
    if (bool(getattr(P, 'joint_nested_choice', False)) or readiness_gate_active(P)
            or parent_age_maturation_active(P) or bool(getattr(P, 'two_shock_choice', False))
            or bool(getattr(P, 'fertility_nest_choice', False))):
        raise ValueError('Birth-count experiment does not support joint/readiness/parent-age choice variants.')
    cap = getattr(P, 'birth_count_choice_cap', 3)
    if isinstance(cap, bool) or int(cap) != cap or int(cap) not in (1, 2, 3):
        raise ValueError('birth_count_choice_cap must be 1, 2, or 3.')


def binomial_probabilities(k, pi):
    if not 0 <= int(k) <= 3 or int(k) != k or not np.isfinite(pi) or not 0 <= pi <= 1:
        raise ValueError('Expected k=0..3 and finite pi in [0,1].')
    out = np.zeros(4)
    for x in range(int(k) + 1):
        out[x] = math.comb(int(k), x) * pi**x * (1 - pi)**(int(k) - x)
    return out


def birth_count_menu(VI, n, m, pi, kappa, first_birth_fixed_cost=0., cap=3):
    """Return value, padded action probabilities, realized probabilities, U(k).

    VI is the existing housing/consumption optimized value (...,n,m), with no
    changes to its budget or maturation. The first-birth fixed cost is paid
    exactly once iff n==0 and x>0. An entirely dead menu has zero probabilities.
    """
    VI = np.asarray(VI, dtype=float)
    if not (0 <= m <= n <= 3) or VI.shape[-2:] != (4, 4):
        raise ValueError('Expected supported n,m and VI ending in (4,4).')
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError('The retained fertility noise scale must be positive.')
    kmax = min(int(cap), 3 - n)
    utilities = np.full(VI.shape[:-2] + (4,), -np.inf)
    kernels = np.zeros((4, 4))
    for k in range(kmax + 1):
        kernels[k] = binomial_probabilities(k, pi)
        u = np.zeros(VI.shape[:-2])
        for x in range(k + 1):
            if kernels[k, x] > 0:
                u += kernels[k, x] * (VI[..., n + x, m + x]
                    - (first_birth_fixed_cost if n == 0 and x > 0 else 0.))
        utilities[..., k] = u
    # Include precisely the existing alternatives at cap=1. Padding must not
    # enter the logsumexp or alter its scale.
    scaled = utilities[..., :kmax + 1] / kappa
    maximum = np.max(scaled, axis=-1, keepdims=True)
    exponentials = np.exp(scaled - maximum)
    denominator = np.sum(exponentials, axis=-1, keepdims=True)
    ls = np.squeeze(maximum + np.log(denominator), axis=-1)
    # Computing exp(scaled-ls) loses normalization when the utility level is
    # very negative. Normalize in centered units, without changing the noise.
    action_short = exponentials / denominator
    dead = np.max(utilities[..., :kmax + 1], axis=-1) <= DEAD_VALUE_CUTOFF
    action_short[dead, :] = 0.
    action = np.zeros_like(utilities)
    action[..., :kmax + 1] = action_short
    realized = action @ kernels
    return kappa * ls, action, realized, utilities


def identity_realized_probabilities(shape):
    """P(X=0)=1 at valid family states, including non-fertile periods."""
    out = np.zeros(tuple(shape) + (4,))
    for n in range(4):
        for m in range(n + 1):
            out[..., n, m, 0] = 1.
    return out


class BirthCountProbabilityError(ValueError):
    """Strict normalization failure with occupied-state diagnostic evidence."""
    def __init__(self, evidence):
        self.evidence = evidence
        super().__init__('Birth-count occupied-state normalization failure: ' + json.dumps(evidence, sort_keys=True))


def _probability_failure(pre, probabilities, action_probs, *, kind, age_index, state_values):
    sums = np.sum(probabilities, axis=-1)
    bad = (pre > 0) & (np.abs(sums - 1.) > 1e-12)
    values = None if state_values is None else np.asarray(state_values, dtype=float)
    if values is not None and values.shape != pre.shape:
        raise ValueError('Diagnostic state_values must have pre-distribution shape.')
    zero = bad & (sums == 0.)
    dead = np.zeros_like(bad) if values is None else bad & (values <= DEAD_VALUE_CUTOFF)
    evidence = dict(kind=kind, age_index=age_index, tolerance=1e-12,
        invalid_cell_count=int(np.sum(bad)), invalid_mass=float(np.sum(pre[bad])),
        probability_sum_min=float(np.min(sums[bad])), probability_sum_max=float(np.max(sums[bad])),
        zero_menu_cell_count=int(np.sum(zero)), zero_menu_mass=float(np.sum(pre[zero])),
        dead_value_cutoff=float(DEAD_VALUE_CUTOFF), dead_status_available=values is not None,
        dead_invalid_mass=float(np.sum(pre[dead])) if values is not None else None,
        living_invalid_mass=float(np.sum(pre[bad & ~dead])) if values is not None else None,
        cells=[])
    action_sums = None if action_probs is None else np.sum(action_probs, axis=-1)
    for idx in np.argwhere(bad)[:8]:
        key = tuple(int(x) for x in idx)
        evidence['cells'].append(dict(index=list(key), mass=float(pre[key]),
            probability_sum=float(sums[key]),
            action_probability_sum=None if action_sums is None else float(action_sums[key]),
            value=None if values is None else float(values[key]),
            dead=None if values is None else bool(values[key] <= DEAD_VALUE_CUTOFF)))
    raise BirthCountProbabilityError(evidence)


def birth_count_transition(pre, realized_probs, action_probs=None, cap=3, *, age_index=None, state_values=None):
    """Propagate one immutable pre-birth snapshot; return auditable flow objects.

    post and first_birth_tagged_post have pre.shape. order_tagged_post has
    (order=3,)+pre.shape; a crossing of several orders is tagged in each.
    births_by_order counts children at each order; expected_births sums these.
    any_birth_mass counts households with X>0. at_risk_by_order includes n<q
    for which the available intended-count menu can reach order q. attempts
    counts P(k>=q-n), not expected k, and shares that denominator. This reduces
    to the existing binary order-specific risk/attempt flow when cap=1.
    """
    pre = np.asarray(pre, dtype=float)
    rp = np.asarray(realized_probs, dtype=float)
    if pre.shape[-2:] != (4, 4) or rp.shape != pre.shape + (4,):
        raise ValueError('Expected pre (...,4,4) and realized_probs (...,4,4,4).')
    if not np.all(np.isfinite(pre)) or np.any(pre < 0) or not np.all(np.isfinite(rp)) or np.any(rp < 0):
        raise ValueError('Mass and realized probabilities must be finite and nonnegative.')
    if isinstance(cap, bool) or int(cap) != cap or int(cap) not in (1, 2, 3):
        raise ValueError('The action cap must be 1, 2, or 3.')
    if action_probs is not None:
        action_probs = np.asarray(action_probs, dtype=float)
        if action_probs.shape != rp.shape:
            raise ValueError('Action probabilities must have realized-probability shape.')
        if not np.all(np.isfinite(action_probs)) or np.any(action_probs < 0):
            raise ValueError('Action probabilities must be finite and nonnegative.')
    # The native occupied-dead gate already admits at most DEAD_MASS_TOL
    # numerical tail mass. The legacy binary zero-choice sentinel leaves that
    # mass in its origin with zero births; reproduce that identity exactly.
    # This is neither dropping mass nor renormalizing a nonzero policy.
    zero_menu = (pre > 0) & (np.sum(rp, axis=-1) == 0.)
    zero_menu_mass = float(np.sum(pre[zero_menu]))
    if zero_menu_mass > DEAD_MASS_TOL:
        _probability_failure(pre, rp, action_probs, kind='zero_menu_mass_exceeds_native_dead_tolerance',
            age_index=age_index, state_values=state_values)
    if state_values is not None:
        values = np.asarray(state_values, dtype=float)
        if values.shape != pre.shape:
            raise ValueError('Diagnostic state_values must have pre-distribution shape.')
        if np.any(zero_menu & (values > DEAD_VALUE_CUTOFF)):
            _probability_failure(pre, rp, action_probs, kind='zero_menu_living_source',
                age_index=age_index, state_values=state_values)
    if action_probs is not None and np.any(zero_menu & (np.sum(action_probs, axis=-1) != 0.)):
        _probability_failure(pre, rp, action_probs, kind='zero_realized_menu_with_nonzero_action_menu',
            age_index=age_index, state_values=state_values)
    post = np.zeros_like(pre)
    tagged = np.zeros((3,) + pre.shape)
    born = np.zeros(3)
    risk = np.zeros(3)
    attempts = np.zeros(3)
    any_birth_mass = 0.
    expected_by_cell = np.zeros(pre.shape[:-2])
    for n in range(4):
        for m in range(4):
            source = pre[..., n, m]
            if m > n:
                if np.any(source > 0):
                    raise ValueError('Positive household mass at m>n.')
                continue
            probs = rp[..., n, m, :]
            if np.any(probs[..., min(int(cap), 3 - n) + 1:] > 0):
                raise ValueError('Birth increment exceeds remaining child cap.')
            sentinel = zero_menu[..., n, m]
            occupied = (source > 0) & ~sentinel
            post[..., n, m] += np.where(sentinel, source, 0.)
            if np.any(np.abs(np.sum(probs, axis=-1)[occupied] - 1) > 1e-12):
                _probability_failure(pre, rp, action_probs, kind='realized',
                    age_index=age_index, state_values=state_values)
            if action_probs is not None:
                actions = action_probs[..., n, m, :]
                if np.any(actions[..., min(int(cap), 3 - n) + 1:] > 0):
                    raise ValueError('Intended birth count exceeds the available menu.')
                if np.any(np.abs(np.sum(actions, axis=-1)[occupied] - 1) > 1e-12):
                    _probability_failure(pre, action_probs, action_probs, kind='action',
                        age_index=age_index, state_values=state_values)
            for q in range(n + 1, min(3, n + int(cap)) + 1):
                risk[q - 1] += float(source.sum())
                if action_probs is not None:
                    attempts[q - 1] += float(np.sum(source * np.sum(action_probs[..., n, m, q - n:], axis=-1)))
            for x in range(4 - n):
                flow = source * probs[..., x]
                post[..., n + x, m + x] += flow
                if x:
                    any_birth_mass += float(flow.sum())
                    expected_by_cell += x * flow
                    for q in range(n + 1, n + x + 1):
                        tagged[q - 1, ..., n + x, m + x] += flow
                        born[q - 1] += float(flow.sum())
    return dict(post=post, births_by_order=born, expected_births=float(born.sum()),
        expected_births_by_cell=expected_by_cell, any_birth_mass=any_birth_mass,
        first_birth_tagged_post=tagged[0].copy(), order_tagged_post=tagged,
        at_risk_by_order=risk, attempts_by_order=attempts,
        zero_menu_identity_mass=zero_menu_mass, zero_menu_identity_cell_count=int(np.sum(zero_menu)))


def transition_at_age(pre, P, j, *, state_values=None):
    """Apply the saved menu outcome distribution at one pre-birth age."""
    return birth_count_transition(pre, P.birth_count_realized_probs[:, :, :, j],
        P.birth_count_action_probs[:, :, :, j], getattr(P, 'birth_count_choice_cap', 3),
        age_index=int(j), state_values=state_values)
