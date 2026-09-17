"""Simple simultaneous fertility GEV over complete housing plans.

Wait alternatives contain one product; attempt alternatives contain one product
per positive-probability conception outcome. There is ONE housing scale within
each fertility nest and no additional outcome-dependent subnests.
"""
import numpy as np


def logsum(values, scale):
    values = np.asarray(values, dtype=float)
    scale = np.broadcast_to(np.asarray(scale, dtype=float), values.shape[:-1])
    if np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("Choice scales must be finite and positive")
    if np.isnan(values).any() or np.isposinf(values).any():
        raise ValueError("Invalid deterministic choice value")
    maximum = np.max(values, axis=-1, keepdims=True)
    delta = np.full_like(values, -np.inf)
    finite = np.isfinite(maximum)
    np.subtract(values, maximum, out=delta, where=np.broadcast_to(finite, values.shape))
    exponent = np.exp(delta / scale[..., None])
    total = exponent.sum(axis=-1, keepdims=True)
    probability = np.divide(exponent, total, out=np.zeros_like(exponent), where=total > 0)
    logged = np.zeros_like(total)
    np.log(total, out=logged, where=total > 0)
    value = np.where(finite[..., 0], maximum[..., 0] + scale * logged[..., 0], -np.inf)
    return value, probability


def choose(q0, q1, pi, kappa, sigma, cost=0.0, available=True):
    q0 = np.asarray(q0, dtype=float)
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), q0.shape[:-1])
    if not np.isfinite(pi) or not 0 <= pi <= 1:
        raise ValueError("Invalid conception probability")
    if not np.isfinite(kappa) or kappa <= 0 or np.any(~np.isfinite(sigma)) or np.any(sigma < kappa):
        raise ValueError("Simple fertility GEV requires sigma >= kappa > 0")
    if not np.isfinite(cost):
        raise ValueError("Invalid birth cost")
    wait_value, wait_housing = logsum(q0, kappa)
    failure_housing = np.zeros_like(q0)
    success_housing = np.zeros_like(q0)
    attempt_value = np.full(q0.shape[:-1], -np.inf)
    if available:
        attempt_value = np.zeros(q0.shape[:-1])
        if pi < 1:
            failure_value, failure_housing = logsum((1 - pi) * q0, kappa)
            attempt_value += failure_value
        if pi > 0:
            q1 = np.asarray(q1, dtype=float)
            if q1.shape != q0.shape:
                raise ValueError("Outcome housing menus must have matching product axes")
            success_value, success_housing = logsum(pi * q1, kappa)
            attempt_value += success_value - pi * cost
    value, action = logsum(np.stack((wait_value, attempt_value), axis=-1), sigma)
    feasible_attempt = np.isfinite(attempt_value)[..., None]
    failure_housing = np.where(feasible_attempt, failure_housing, 0.0)
    success_housing = np.where(feasible_attempt, success_housing, 0.0)
    occupied = np.isfinite(value)
    if np.any(abs(action.sum(axis=-1)[occupied] - 1) > 2e-11):
        raise RuntimeError("Fertility-nest probability accounting failed")
    return dict(value=value, action_probability=action, wait_housing=wait_housing,
                failure_housing=failure_housing, success_housing=success_housing)


def bellman_block(Vd, kernel_args, P, j, fecundity, deterministic_kernel):
    """Retain every original product, budget and post-conception saving policy."""
    if np.any(kernel_args[-1]):
        raise NotImplementedError("Fertility-nest entry grants are outside this experiment")
    nt = Vd.shape[1]
    q = np.empty(Vd.shape + (nt,))
    for product in range(nt):
        restricted = np.full_like(Vd, -1e10)
        restricted[:, product] = Vd[:, product]
        value, _ = deterministic_kernel(restricted, *kernel_args, True)
        q[..., product] = np.where(value > -1e9, value + float(P.E_loc[0] - P.mu_stay), -np.inf)
    kappa = float(P.tenure_choice_kappa)
    value, wait = logsum(q, kappa)
    probabilities = np.zeros(q.shape + (2,))
    probabilities[..., 0] = wait
    failure = np.zeros_like(q)
    fertile = int(P.A_f_start) <= j + 1 <= int(P.A_f_end)
    pi = float(fecundity[j])
    if fertile:
        for nn in range(P.n_parity - 1):
            raw = P.kappa_fert if nn == 0 else getattr(P, "kappa_fert_continuation", None)
            sigma = float(P.kappa_fert if raw is None else raw)
            for cs in range(nn + 1):
                no = q[..., nn, cs, :]
                yes = q[..., nn + 1, cs + 1, :]
                result = choose(no, yes, pi, kappa, sigma, float(P.first_birth_fixed_cost) if nn == 0 else 0.0)
                value[..., nn, cs] = result["value"]
                action = result["action_probability"]
                probabilities[..., nn, cs, :, 0] = action[..., 0, None] * result["wait_housing"]
                # This branch stores the action marginal even at pi=0, where
                # the success housing coordinate is economically absent.
                housing_success = result["success_housing"] if pi > 0 else result["failure_housing"]
                probabilities[..., nn, cs, :, 1] = action[..., 1, None] * housing_success
                failure[..., nn, cs, :] = action[..., 1, None] * result["failure_housing"]
    products = np.broadcast_to(np.arange(nt, dtype=np.int16), q.shape)
    return np.where(np.isfinite(value), value, -1e10), probabilities, products, wait, failure


def factor_age(g_pre, joint, P, j, mode='natural'):
    """Realize the independent conception outcome after joint plan selection."""
    from .parameters import get_fecundity_by_age
    if mode not in ('natural', 'wait', 'first_birth_treated', 'first_birth_control'):
        raise ValueError('Unknown fertility-nest population mode')
    if g_pre.shape[2] != 1 or bool(getattr(P, 'birth_entry_grant', False)):
        raise NotImplementedError('Fertility-nest experiment requires one location and no entry grants')
    pi = float(get_fecundity_by_age(P)[j])
    fertile = int(P.A_f_start) <= j + 1 <= int(P.A_f_end)
    nt = g_pre.shape[1]
    post = np.zeros_like(g_pre)
    weighted = np.zeros(g_pre.shape + (nt,))
    births = np.zeros(P.n_parity)
    attempts = np.zeros_like(births)
    risk = np.zeros_like(births)
    selected_cohort = mode.startswith('first_birth_')
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            mass = g_pre[..., nn, cs]
            if not np.any(mass):
                continue
            if cs > nn:
                raise RuntimeError('Occupied impossible child-count state')
            eligible = fertile and nn < P.n_parity - 1
            if selected_cohort and (nn != 0 or cs != 0 or not eligible):
                continue
            pr = joint.probabilities[:, :, :, j, :, nn, cs]
            wait = joint.wait_probabilities[:, :, :, j, :, nn, cs]
            attempt_mass = pr[..., :, 1].sum(axis=-1)
            if mode == 'natural':
                outcomes = [(nn, cs, pr[..., :, 0])]
                if eligible:
                    failure = joint.failure_probabilities[:, :, :, j, :, nn, cs]
                    outcomes.extend([(nn, cs, (1 - pi) * failure), (nn + 1, cs + 1, pi * pr[..., :, 1])])
                    risk[nn] += mass.sum()
                    attempts[nn] += np.sum(mass * attempt_mass)
                    births[nn] += pi * np.sum(mass * attempt_mass)
            elif mode == 'wait':
                outcomes = [(nn, cs, wait)]
            else:
                births[0] += pi * np.sum(mass * attempt_mass)
                if mode == 'first_birth_treated':
                    outcomes = [(1, 1, pi * pr[..., :, 1])]
                else:
                    # Preserve production measurement: selected origin states
                    # held childless use the restricted-wait housing policy,
                    # not the different unsuccessful-attempt plan policy.
                    outcomes = [(0, 0, pi * attempt_mass[..., None] * wait)]
            for dn, dc, product_probability in outcomes:
                selected = mass[..., None] * product_probability
                post[..., dn, dc] += selected.sum(axis=-1)
                weighted[..., dn, dc, :] += selected
    effective = np.divide(weighted, post[..., None], out=np.zeros_like(weighted), where=post[..., None] > 0)
    expected = births.sum() if selected_cohort else float(g_pre.sum())
    if abs(float(post.sum()) - expected) > 2e-10 * max(1.0, expected):
        raise RuntimeError('Fertility-nest product factorization lost occupied mass')
    return post, effective, births, attempts, risk
