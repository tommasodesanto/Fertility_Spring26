"""Experimental simultaneous tenure/attempt GEV and exact mass factorization.

Two tenure nests, a common dissimilarity, deterministic product choice within
committed tenure after conception. No location choice or readiness states.
Effective product probabilities below depend on the supplied distribution;
they are a lossless computational compression, never a household policy.
"""
from types import SimpleNamespace
import numpy as np


def logsum_prob(values, scale, axis=-1):
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Invalid positive GEV scale')
    maximum = np.max(values, axis=axis, keepdims=True)
    feasible = np.isfinite(maximum)
    centered = np.full_like(values, -np.inf)
    np.subtract(values, maximum, out=centered, where=np.broadcast_to(feasible, values.shape))
    exponent = np.exp(centered / scale)
    total = exponent.sum(axis=axis, keepdims=True)
    probability = np.divide(exponent, total, out=np.zeros_like(exponent), where=total > 0)
    logged = np.zeros_like(total)
    np.log(total, out=logged, where=total > 0)
    return np.squeeze(np.where(feasible, maximum + scale * logged, -np.inf), axis), probability


def allocate(shape, P):
    if P.I != 1 or not bool(getattr(P, 'sequential_births', False)):
        raise NotImplementedError('Joint experiment requires one market and sequential parity states')
    if str(getattr(P, 'child_maturation_mode', '')) not in ('independent', 'independent_count'):
        from .parameters import independent_child_maturation_active
        if not independent_child_maturation_active(P):
            raise NotImplementedError('Joint experiment requires independent child-count maturation')
    from .parameters import readiness_gate_active
    if readiness_gate_active(P):
        raise NotImplementedError('Readiness extension is outside this joint experiment')
    kappa, lam = float(P.tenure_choice_kappa), float(P.joint_nest_lambda)
    if not (np.isfinite(kappa) and kappa > 0 and np.isfinite(lam) and 0 < lam <= 1):
        raise ValueError('Require kappa>0 and 0<lambda<=1')
    return SimpleNamespace(probabilities=np.zeros(shape + (2, 2)),
                           products=np.zeros(shape + (2,), dtype=np.int16),
                           wait_probabilities=np.zeros(shape + (2,)))


def bellman_block(Vd, kernel_args, P, j, fecundity, deterministic_kernel):
    """Return EV, p(tenure,attempt), product plans, restricted-wait tenure law."""
    if np.any(kernel_args[-1]):
        raise NotImplementedError("Joint experimental entry grants are not implemented")
    rental = Vd.copy(); rental[:, 1:] = -1e10
    qr, cr = deterministic_kernel(rental, *kernel_args)
    owner = Vd.copy(); owner[:, 0] = -1e10
    qo, co = deterministic_kernel(owner, *kernel_args)
    q = np.stack((qr, qo), axis=-1)
    q = np.where(q > -1e9, q + float(P.E_loc[0] - P.mu_stay), -np.inf)
    products = np.stack((cr, co), axis=-1)
    kappa, lam = float(P.tenure_choice_kappa), float(P.joint_nest_lambda)
    _, wait = logsum_prob(q, kappa)
    plans = np.full(q.shape + (2,), -np.inf)
    plans[..., 0] = q
    fertile = int(P.A_f_start) <= j + 1 <= int(P.A_f_end)
    pi = float(fecundity[j])
    if fertile:
        for nn in range(P.n_parity - 1):
            for cs in range(nn + 1):
                no = q[..., nn, cs, :]
                yes = q[..., nn + 1, cs + 1, :]
                attempt = np.zeros_like(no)
                if pi < 1:
                    attempt += (1 - pi) * no
                if pi > 0:
                    cost = float(P.first_birth_fixed_cost) if nn == 0 else 0.0
                    attempt += pi * (yes - cost)
                plans[..., nn, cs, :, 1] = attempt
    inclusive, action = logsum_prob(plans, kappa * lam)
    value, tenure = logsum_prob(inclusive, kappa)
    probability = tenure[..., :, None] * action
    return np.where(np.isfinite(value), value, -1e10), probability, products, wait


def factor_age(g_pre, joint, P, j, mode='natural'):
    """Compress selected joint-plan mass before transaction at a single age."""
    if bool(getattr(P, 'birth_entry_grant', False)):
        raise NotImplementedError('Joint experimental entry grants are not implemented')
    modes = ('natural', 'wait', 'first_birth_treated', 'first_birth_control')
    if mode not in modes:
        raise ValueError(f'Unknown joint transition mode {mode}')
    if g_pre.shape[2] != 1:
        raise NotImplementedError('One market only')
    from .parameters import get_fecundity_by_age
    pi = float(get_fecundity_by_age(P)[j])
    fertile = int(P.A_f_start) <= j + 1 <= int(P.A_f_end)
    nt = g_pre.shape[1]
    post = np.zeros_like(g_pre)
    weighted = np.zeros(g_pre.shape + (nt,))
    births = np.zeros(P.n_parity)
    attempts = np.zeros_like(births); risk = np.zeros_like(births)
    selection = mode.startswith('first_birth_')
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            mass = g_pre[..., nn, cs]
            if not np.any(mass):
                continue
            if cs > nn:
                raise RuntimeError('Occupied impossible child-count state')
            eligible = fertile and nn < P.n_parity - 1
            if selection and (nn != 0 or cs != 0 or not eligible):
                continue
            pr = joint.probabilities[:, :, :, j, :, nn, cs, :, :]
            if eligible and mode == 'natural':
                risk[nn] += mass.sum()
                attempts[nn] += np.sum(mass * pr[..., :, 1].sum(axis=-1))
                births[nn] += pi * np.sum(mass * pr[..., :, 1].sum(axis=-1))
            if mode == 'natural':
                outcomes = [(nn, cs, pr[..., :, 0])]
                if eligible:
                    outcomes += [(nn, cs, (1 - pi) * pr[..., :, 1]),
                                 (nn + 1, cs + 1, pi * pr[..., :, 1])]
            elif mode == 'wait':
                outcomes = [(nn, cs, joint.wait_probabilities[:, :, :, j, :, nn, cs, :])]
            else:
                dest = (1, 1) if mode == 'first_birth_treated' else (0, 0)
                outcomes = [(dest[0], dest[1], pi * pr[..., :, 1])]
                births[0] += pi * np.sum(mass * pr[..., :, 1].sum(axis=-1))
            for dn, dc, tenure_mass in outcomes:
                products = joint.products[:, :, :, j, :, dn, dc, :]
                for tau in range(2):
                    selected = mass * tenure_mass[..., tau]
                    invalid_product = (products[..., tau] != 0) if tau == 0 else (products[..., tau] <= 0)
                    if np.any((selected > 0) & invalid_product):
                        raise RuntimeError("Infeasible forced joint-plan outcome within committed tenure")
                    post[..., dn, dc] += selected
                    for product in range(nt):
                        weighted[..., dn, dc, product] += selected * (products[..., tau] == product)
    effective = np.divide(weighted, post[..., None], out=np.zeros_like(weighted), where=post[..., None] > 0)
    expected = births.sum() if selection else float(g_pre.sum())
    if abs(float(post.sum()) - expected) > 2e-10 * max(1.0, expected):
        raise RuntimeError('Joint factorization lost occupied mass or used infeasible plans')
    return post, effective, births, attempts, risk


def factor_distribution(g_pre, joint, P, mode='natural'):
    if joint is None:
        raise RuntimeError('Joint policy missing matching plan object')
    post = np.empty_like(g_pre)
    effective = np.empty(g_pre.shape + (g_pre.shape[1],))
    births = np.zeros((P.J, P.n_parity)); attempts = np.zeros_like(births); risk = np.zeros_like(births)
    for j in range(P.J):
        result = factor_age(g_pre[:, :, :, j], joint, P, j, mode)
        post[:, :, :, j], effective[:, :, :, j], births[j], attempts[j], risk[j] = result
    return post, effective, births, attempts, risk


def stationary_first_birth_response(g_pre, joint, P, bg, SD, lp, tc, bp, hr, maps):
    """Fixed-price dated matched branch, independently of calibration wrappers."""
    from . import solver as model
    lidx, lwt, tidx, twt = maps
    _, _, pi_z = model.income_transition_values(P)
    stochastic = bool(P.use_stochastic_aging and hasattr(P, 'Pi_child'))
    means = []
    masses = []
    for mode in ('first_birth_treated', 'first_birth_control'):
        post, effective, *_ = factor_distribution(g_pre, joint, P, mode)
        future = np.zeros_like(post)
        for j in range(P.J - 1):
            survival = float(P.survival_probs[j]) if P.use_age_survival else 1.0
            future[:, :, :, j + 1] = model.advance_cohort_one_period_markov_income(
                survival * post[:, :, :, j], j, lp, tc, effective, bp, P, bg, SD,
                lidx, lwt, tidx, twt, stochastic, P.Pi_child if stochastic else None, pi_z)
        post, effective, *_ = factor_distribution(
            future, joint, P, 'natural' if mode.endswith('treated') else 'wait')
        current = model.realize_current_cross_section(
            post, lp, tc, effective, lidx, lwt, tidx, twt,
            use_compiled_scatter=bool(getattr(P, 'use_numba_scatter', False)))
        mass = float(current.sum()); masses.append(mass)
        housing = float(np.sum(current[:, 0] * hr[:, 0]))
        for ten in range(1, 1 + P.n_house):
            housing += float(current[:, ten].sum()) * float(P.H_own[ten - 1])
        means.append(housing / max(mass, 1e-300))
    if abs(masses[0] - masses[1]) > 2e-10 or min(masses) <= 1e-14:
        raise RuntimeError('Invalid stationary matched joint branch mass')
    return means[0] - means[1]


def action_marginals(probability):
    """Binary marginals with complementary rounding and exact [0,1] support."""
    raw = probability.sum(axis=-2)
    if not np.isfinite(raw).all() or raw.min() < -1e-14 or raw.max() > 1 + 1e-14:
        raise RuntimeError("Invalid GEV action marginal before roundoff correction")
    attempt = np.clip(raw[..., 1], 0.0, 1.0)
    wait = np.where(raw.sum(axis=-1) > 0, 1.0 - attempt, 0.0)
    return np.stack((wait, attempt), axis=-1)
