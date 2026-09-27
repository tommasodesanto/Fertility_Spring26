"""Experimental natural-credit adapter; NOT a production model default.

Install only after the authenticated purchase-income adapter. Preserve all
parameters, including psi_child, and use a direct stationary equilibrium call.
The admissible savings set is the grid-supported upper interval on which every
positive-probability continuation is feasible and estates repay at death.
This conservative nodal boundary is a discretization, not a continuum proof.
"""
from __future__ import annotations
import inspect
import difflib
from pathlib import Path
import numpy as np

DEAD = -1e10
CUTOFF = -1e9


def natural_support_floor(values, grid):
    """First node in an upper feasible interval; reject holes, preserve dead cols."""
    v = np.asarray(values)
    bg = np.asarray(grid)
    if v.ndim != 2 or v.shape[0] != bg.size or np.any(np.diff(bg) <= 0):
        raise ValueError('Invalid continuation/grid shape')
    if not np.isfinite(v).all():
        raise ValueError('Nonfinite continuation value')
    feasible = v > CUTOFF
    if np.any(feasible[:-1] & ~feasible[1:]):
        raise ValueError('Continuation feasibility is not an upper interval')
    dead = ~feasible.any(axis=0)
    first = np.argmax(feasible, axis=0)
    floors = bg[first].astype(float)
    floors[dead] = bg[-1]
    return floors, dead


def net_estate_value(liquid, house_cost, sale_cost):
    return np.asarray(liquid) + (1.0 - float(sale_cost)) * np.asarray(house_cost)


def strict_expectation(values, probabilities):
    """Independent reference: any reachable bankruptcy rules out saving."""
    v = np.asarray(values, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    if v.shape[0] != p.size or np.any(p < 0) or not np.isclose(p.sum(), 1):
        raise ValueError('Invalid probability support')
    out = np.tensordot(p, v, axes=(0, 0))
    bad = np.any((v <= CUTOFF) & (p > 0).reshape((-1,) + (1,) * (v.ndim-1)), axis=0)
    return np.where(bad, DEAD, out)


def _once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('Unexpected source anchor: ' + old[:100])
    return source.replace(old, new)


def _kernel_wrapper(original, renter):
    """Replace artificial floor by the strict continuation-domain floor."""
    def wrapped(*args):
        args = list(args)
        floor, dead = natural_support_floor(args[2], args[5])
        if renter:
            # New experimental kernel argument is after exhaustive_saving.
            args.extend([floor])
        else:
            args[12] = floor
            args[21] = 0.0
            args[22] = 0.0
        result = list(original(*args))
        if np.any(dead):
            result[0][:, dead] = DEAD
        return tuple(result)
    return wrapped


def install(model, output, *, enabled=False):
    """Explicit opt-in, source-auditable patch. Baseline disabled is a true no-op.

    No model solve is performed. Lead must review generated source and test a
    full baseline replay before enabling this on the economic benchmark.
    """
    if not enabled:
        return {'enabled': False, 'baseline_noop': True}
    if getattr(model, '_solvency_credit_benchmark_installed', False):
        raise ValueError('Already installed')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    original = model.solve_bellman_full_markov_income
    source = inspect.getsource(original)
    source_before = source
    if 'bmo_purchase' not in source:
        raise ValueError('Authenticated purchase-income adapter must run first')
    source = _once(source, '    t0 = time.perf_counter()\n', '''    t0 = time.perf_counter()
    if continuation_V is not None or not bool(getattr(P, "exhaustive_saving_control", False)):
        raise ValueError("Natural-credit draft supports exhaustive stationary Bellman only")
''')
    # Death wealth convention follows the provisional net-liquidation ledger.
    source = _once(source,
        '                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)',
        '''                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)
                    Vbq[b_grid + (1.0 - P.psi) * hv < 0.0, ten, i, nn, cs] = -1e10''')
    source = _once(source, '                Vnr = np.zeros((Nb, nt, I, npar, ncs))',
        '                Vnr = np.zeros((Nb, nt, I, npar, ncs))\n                continuation_bad = np.zeros_like(Vnr, dtype=bool)')
    anchor = '''                        Vnr += transition_weight * next_values[
                            :, :, :, j + 1, znext, :, :
                        ]'''
    source = _once(source, anchor, anchor + '''
                        continuation_bad |= next_values[:, :, :, j + 1, znext, :, :] <= DEAD_VALUE_CUTOFF''')
    anchor = '''                    Vnr = survival * Vnr + (1.0 - survival) * Vbq'''
    source = _once(source, anchor, anchor + '''
                    continuation_bad = ((survival > 0.0) & continuation_bad) | ((survival < 1.0) & (Vbq <= DEAD_VALUE_CUTOFF))
                Vnr[continuation_bad] = -1e10''')
    anchor = '            Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs, age_index=j)'
    source = _once(source, anchor, anchor + '''
            continuation_bad = apply_child_aging((Vnr <= DEAD_VALUE_CUTOFF).astype(float), P, Nb, nt, I, npar, ncs, age_index=j) > 0.0
            Vc[continuation_bad] = -1e10''')
    source = _once(source, '            dp_choice = dp_arr - income_for_purchase',
        '            dp_choice = np.full_like(dp_arr, -np.inf)')
    source = _once(source, '            bmo_purchase = np.maximum(bmo - income_for_purchase, b_grid[0])',
        '            bmo_purchase = np.full_like(bmo, b_grid[0])')
    # Add one floor vector to already allocation-corrected renter kernel.
    kernel = model.full_renter_block_kernel
    py = getattr(kernel, 'py_func', kernel)
    ks = inspect.getsource(py)
    kernel_before = ks
    ks = _once(ks, '    exhaustive_saving=0,\n', '    exhaustive_saving=0,\n    natural_floor_v=None,\n')
    ks = _once(ks, '            unsecured_floor = rollover_floor if rollover_floor < line_floor else line_floor', '            unsecured_floor = natural_floor_v[c]')
    kp = out / 'natural_credit_renter.generated.py'
    kp.write_text(ks)
    ns = dict(py.__globals__)
    exec(compile(ks, str(kp), 'exec'), ns)
    model.full_renter_block_kernel = _kernel_wrapper(ns['full_renter_block_kernel'], True)
    model.full_owner_block_kernel = _kernel_wrapper(model.full_owner_block_kernel, False)
    bp = out / 'natural_credit_bellman.generated.py'
    bp.write_text(source)
    (out / "natural_credit.diff").write_text("".join(difflib.unified_diff(source_before.splitlines(True), source.splitlines(True), fromfile="authenticated_bellman", tofile="natural_credit_bellman")) + "".join(difflib.unified_diff(kernel_before.splitlines(True), ks.splitlines(True), fromfile="authenticated_renter", tofile="natural_credit_renter")))
    exec(compile(source, str(bp), 'exec'), model.__dict__)
    model._solvency_credit_benchmark_installed = True
    return dict(enabled=True, source=str(bp), renter_kernel=str(kp),
                status='DRAFT_REQUIRES_LEAD_REVIEW_AND_FULL_OPERATOR_AUDITS',
                grid_convention='strict nodal continuation support; retain original grid',
                classification_limitation='inherits native V > -1e9 feasibility classifier; very negative finite utility is not distinguished from bankruptcy')


def audit_solvency_arrays(mass, saving, house_cost, death_probability, sale_cost, grid):
    """Independent broadcast-array audit, after all policies are realized.

    Caller must provide current/post-tenure mass, with age-aligned death
    probabilities and terminal probability exactly one. No continuation values
    or borrowing parameters enter this test.
    """
    g, bp, hc, dp = np.broadcast_arrays(mass, saving, house_cost, death_probability)
    if np.any(g < 0) or np.any(dp < 0) or np.any(dp > 1):
        raise ValueError('Invalid masses/death probabilities')
    if not all(np.isfinite(x).all() for x in (g, bp, hc, dp)):
        raise ValueError('Nonfinite audit input')
    estate = net_estate_value(bp, hc, sale_cost)
    tolerance = 1e-10
    bad = (dp > 0) & (estate < -tolerance)
    return dict(negative_estate_exposure_mass=float(g[bad].sum()),
                negative_estate_death_mass=float((g * dp)[bad].sum()),
                negative_estate_liability=float(np.sum(g * dp * np.maximum(-estate, 0))),
                saving_at_lower_grid_mass=float(g[bp <= grid[0] + tolerance].sum()),
                saving_at_upper_grid_mass=float(g[bp >= grid[-1] - tolerance].sum()),
                saving_outside_grid_mass=float(g[(bp < grid[0]-tolerance) | (bp > grid[-1]+tolerance)].sum()))


def death_probabilities(P):
    """J-1 survival transitions followed by certain terminal death."""
    J = int(P.J)
    if J < 1:
        raise ValueError('Positive lifecycle length required')
    deaths = np.zeros(J)
    if bool(getattr(P, 'use_age_survival', False)):
        survival = np.asarray(P.survival_probs, dtype=float)
        if (survival.shape != (J-1,) or not np.isfinite(survival).all()
                or np.any(survival < 0) or np.any(survival > 1)):
            raise ValueError('Expected J-1 finite survival probabilities in [0,1]')
        deaths[:-1] = 1.0 - survival
    deaths[-1] = 1.0
    return deaths


def audit_purchase_accounting(evaluation, P, shared, grid, model):
    """Audit exact transactions and death solvency without any LTV assumption.

    Same calling convention as the inherited audit. Native household budget and
    forward-operator audits must ALSO run; this does not replace them.
    """
    if int(P.I) != 1 or evaluation.policy.tenure_probs is None:
        raise ValueError('Requires one market and probabilistic tenure')
    policy = evaluation.policy
    bg = np.asarray(grid)
    costs = np.r_[0., float(policy.price[0]) * np.asarray(P.H_own)]
    sale = (1.0 - float(P.psi)) * costs
    outside_mass = max_error = 0.
    for age in range(P.J):
        for zz in range(len(P.z_grid)):
            for old in range(len(costs)):
                mass = evaluation.g_post_fertility[:, old, 0, age, zz] * policy.loc_probs[:, old, 0, 0, age, zz]
                probs = np.asarray(policy.tenure_probs[:, old, 0, age, zz], dtype=float)
                sums = probs.sum(axis=-1, keepdims=True)
                probs = np.divide(probs, sums, out=np.zeros_like(probs), where=sums > 0)
                for new in range(len(costs)):
                    branch_mass = mass * probs[..., new]
                    x = bg if old == new else bg + sale[old] - costs[new]
                    outside = (x < bg[0]-1e-12) | (x > bg[-1]+1e-12)
                    outside_mass += float(branch_mass[outside].sum())
                    idx = policy.maps.tmx_idx[0, old, new]
                    wt = policy.maps.tmx_wt[0, old, new]
                    mapped = (1-wt)*bg[idx] + wt*bg[idx+1]
                    # Map arrays are (ever born, at home, wealth).
                    error = np.abs(mapped - x[None,None,:]).transpose(2,0,1)
                    occupied = branch_mass > 0
                    if np.any(occupied):
                        max_error = max(max_error, float(error[occupied].max()))
    g = np.asarray(evaluation.g_current)
    bp = np.asarray(policy.bp_pol)
    deaths = death_probabilities(P)
    death_shape = [1]*g.ndim; death_shape[3] = P.J
    cost_shape = [1]*g.ndim; cost_shape[1] = len(costs)
    result = audit_solvency_arrays(g, bp, costs.reshape(cost_shape), deaths.reshape(death_shape), P.psi, bg)
    result.update(transaction_outside_grid_mass=outside_mass,
                  maximum_occupied_transaction_wealth_error=max_error,
                  audit_id='natural_credit_net_estate_solvency_v1',
                  purchase_threshold='none; actual transaction support and conditional budget feasibility only',
                  borrowing_constraint='strict continuation solvency and net estates at possible death',
                  native_value_cutoff_limitation=True)
    if outside_mass > 1e-10 or max_error > 1e-8 or result['negative_estate_exposure_mass'] > 1e-10:
        raise RuntimeError('Natural-credit transaction or solvency audit failed: '+str(result))
    return result
