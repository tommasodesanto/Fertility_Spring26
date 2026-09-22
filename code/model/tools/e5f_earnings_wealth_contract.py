"""Explicit diagnostic entry and purchase-timing adapters for a frozen E5F source.

Nothing is installed by import.  The original source tree remains immutable;
the wrapper records this module's hash and the exact generated source diff.
"""
from __future__ import annotations

import difflib
import inspect
from pathlib import Path
from typing import Any

import numpy as np


def rank_coupled_entry(model: Any, reference: Any, grid: np.ndarray,
                       persistent_weights: np.ndarray, iid_weights: np.ndarray):
    """Preserve reference wealth and its income-rank ordering; iid is independent.

    This is a declared diagnostic coupling, not an estimated joint distribution.
    Probability-interval overlaps map old total-income ranks to new persistent
    ranks.  Conditional transitory draws never rescale wealth.
    """
    if int(reference.I) != 1:
        raise ValueError("entry mapping currently requires the one-market reference")
    z, zw, _ = model.income_transition_values(reference)
    order = np.argsort(z, kind="stable")
    old_joint = np.zeros((len(grid), len(z)))
    for k in range(len(z)):
        idx, wt = model.entry_wealth_grid_weights(grid, reference, i=0, j=0,
                                                z_value=float(z[k]))
        old_joint[idx, k] = np.asarray(wt) * zw[k]
    wp = np.asarray(persistent_weights, dtype=float)
    we = np.asarray(iid_weights, dtype=float)
    for weights in (zw, wp, we):
        if np.any(weights <= 0) or not np.isclose(np.sum(weights), 1., atol=1e-12):
            raise ValueError("entry weights must be positive probabilities")
    old_edges = np.r_[0., np.cumsum(zw[order])]
    new_edges = np.r_[0., np.cumsum(wp)]
    joint = np.zeros((len(grid), len(wp)))
    for k, old_k in enumerate(order):
        conditional = old_joint[:, old_k] / zw[old_k]
        for j in range(len(wp)):
            overlap = max(0., min(old_edges[k+1], new_edges[j+1])
                          - max(old_edges[k], new_edges[j]))
            joint[:, j] += overlap * conditional
    conditional = np.repeat(joint / wp[None, :], len(we), axis=1)
    new_weights = np.kron(wp, we)
    marginal = conditional @ new_weights
    np.testing.assert_allclose(marginal, old_joint.sum(axis=1), atol=2e-14, rtol=0)
    np.testing.assert_allclose(conditional.sum(axis=0), 1., atol=2e-14, rtol=0)
    return conditional, {
        "rule": "reference wealth marginal; old total-income rank coupled to new persistent rank; iid independent conditional on persistent rank",
        "empirical_joint_distribution_estimated": False,
        "wealth_marginal_max_gap": float(np.max(np.abs(marginal-old_joint.sum(axis=1)))),
        "reference_wealth_marginal": old_joint.sum(axis=1).tolist(),
        "candidate_wealth_marginal": marginal.tolist(),
        "reference_wealth_mean": float(grid @ old_joint.sum(axis=1)),
        "candidate_wealth_mean": float(grid @ marginal),
    }


def install_fixed_entry(model: Any):
    original = model.entry_wealth_grid_weights

    def fixed_entry(b_grid, P, *, i=0, j=0, z_value=1.0):
        if not hasattr(P, "fixed_reference_entry_conditional"):
            return original(b_grid, P, i=i, j=j, z_value=z_value)
        if i != 0 or j != 0:
            raise ValueError("fixed entry mapping is defined only at age-zero in one market")
        np.testing.assert_array_equal(b_grid, P.fixed_reference_entry_grid)
        zz = np.flatnonzero(np.asarray(P.z_grid) == z_value)
        if len(zz) != 1:
            raise ValueError("entry lookup requires a unique exact income node")
        weights = np.asarray(P.fixed_reference_entry_conditional)[:, int(zz[0])]
        idx = np.flatnonzero(weights > 0)
        return idx, weights[idx]

    model.entry_wealth_grid_weights = fixed_entry


def _replace_once(source: str, before: str, after: str) -> str:
    if source.count(before) != 1:
        raise ValueError(f"source patch anchor count != 1: {before[:80]}")
    return source.replace(before, after)


def extend_upper_transaction_grid(core, *, upper, extra_nodes):
    """Preserve every original knot and append a declared geometric tail.

    This expands the numerical wealth domain, not an economic saving limit.
    Full transaction-support and occupied-value gates remain necessary; this
    constructor alone does not establish convergence or a reachability bound.
    """
    core = np.asarray(core, dtype=float)
    if (core.ndim != 1 or len(core) != 120 or not np.isfinite(core).all()
            or not np.all(np.diff(core) > 0) or core[-1] != 30.
            or not np.isfinite(upper) or upper <= core[-1]
            or int(extra_nodes) != extra_nodes or extra_nodes < 1):
        raise ValueError("unsupported explicit upper transaction-grid contract")
    result = np.r_[core, np.geomspace(core[-1], upper, int(extra_nodes) + 1)[1:]]
    np.testing.assert_array_equal(result[:len(core)], core)
    return result


def install_explicit_transaction_grid(model):
    original = model.make_grid

    def make_grid(P):
        if not hasattr(P, "earnings_transaction_grid"):
            return original(P)
        grid = np.asarray(P.earnings_transaction_grid, dtype=float)
        if len(grid) != int(P.Nb) or grid[-1] != P.b_max:
            raise ValueError("explicit transaction-grid metadata mismatch")
        return grid.copy()

    model.make_grid = make_grid


def rewrite_probe_grid(source, expected_nodes):
    """Change only the declared numerical geometry; retain all native gates."""
    if expected_nodes <= 120:
        raise ValueError("transaction-grid probe requires an explicit extension")
    source = _replace_once(source, "    grid=np.asarray(old.b_grid).copy()\n",
                           "    grid=np.asarray(model.make_grid(base)).copy()\n")
    return _replace_once(source, "if (len(grid)!=120 or", f"if (len(grid)!={int(expected_nodes)} or")


def _ordinary_forward_maps(P, b_grid, hc, he, phi_choice, birth_dp, birth_entry_grant):
    """Exact transaction wealth, without inventing a collateral-limit transfer.

    Maps outside numerical support are placeholders only.  The backward choice
    gate forbids those purchases; an independent occupied-branch audit verifies
    that their realized probability mass is zero.
    """
    from intergen_eqscale_seq_optimized.utils import interp_indices
    if np.any(birth_dp) or np.any(birth_entry_grant):
        raise ValueError("purchase-income diagnostic does not combine grant/waiver policies")
    bg = np.asarray(b_grid)
    I, nt = hc.shape
    npar, ncs = phi_choice.shape[2:]
    idx = np.zeros((I, nt, nt, npar, ncs, len(bg)), dtype=np.int64)
    wt = np.zeros_like(idx, dtype=float)
    for i in range(I):
        for old in range(nt):
            for new in range(nt):
                x = bg if old == new else bg + he[i, old] - hc[i, new]
                ii, ww = interp_indices(bg, np.clip(x, bg[0], bg[-1]))
                idx[i, old, new, :, :, :] = ii
                wt[i, old, new, :, :, :] = ww
    return idx, wt


def _install_transaction_support(model: Any, diff_path: Path):
    """Exclude unsupported transaction wealth before tenure choice/logit.

    Conditional values are solved only on b_grid. Clipping a sale or purchase
    to an endpoint changes resources. Keep the original interpolation exactly
    within support, and give unsupported transactions the native infeasible
    value. The forward map may retain placeholders only on zero-probability
    branches; the independent occupied-branch audit still checks this.
    """
    from numba import njit
    support = '''def _interp_on_transaction_grid(bg, values, x, strict_interpolated_support=False):
    if x < bg[0] or x > bg[-1]:
        return -1e10
    return _native_transaction_interp(bg, values, x, strict_interpolated_support)
'''
    sources = []
    diffs = []
    for name in ("tenure_choice_kernel", "tenure_logit_kernel"):
        original = getattr(model, name)
        python_function = getattr(original, "py_func", original)
        namespace = dict(python_function.__globals__)
        namespace["_native_transaction_interp"] = namespace["_interp_with_clip"]
        exec(compile(support, str(diff_path.with_suffix(".tenure.py")), "exec"), namespace)
        namespace["_interp_on_transaction_grid"] = njit(cache=False)(namespace["_interp_on_transaction_grid"])
        source = inspect.getsource(python_function)
        revised = _replace_once(source, "@njit(cache=True)", "@njit(cache=False)")
        if "_interp_with_clip(" not in revised:
            raise ValueError("missing native tenure interpolation call sites")
        revised = revised.replace("_interp_with_clip(", "_interp_on_transaction_grid(")
        exec(compile(revised, str(diff_path.with_suffix(".tenure.py")), "exec"), namespace)
        setattr(model, name, namespace[name])
        sources.append(revised)
        diffs.extend(difflib.unified_diff(source.splitlines(True), revised.splitlines(True),
            fromfile=f"frozen/{name}", tofile=f"diagnostic/{name}"))
    diff_path.with_suffix(".tenure.py").write_text(support + "\n".join(sources))
    return "".join(diffs)


def _install_exact_allocation_output(model: Any, diff_path: Path):
    """Clone native full kernels with exact exhaustive-saving output arithmetic.

    The optimizer and its returned value/saving arrays are left untouched.  In
    exhaustive-saving mode only, realized consumption and housing report the
    allocation used by that optimizer even when continuation value is deeply
    negative; genuinely infeasible and nonexhaustive branches retain native
    placeholders/floors.
    """
    from numba import njit

    replacements = {
        "full_renter_block_kernel": (
            "if exhaustive_saving and v_best > -1e9:",
            "if exhaustive_saving:",
            2,
        ),
        "full_owner_block_kernel": (
            "            ct_eff = ct if ct > c_min else c_min\n"
            "            co[b, c] = cbc + ct_eff",
            "            ct_eff = ct if ct > c_min else c_min\n"
            "            if exhaustive_saving and ct > 1e-10:\n"
            "                ct_eff = ct\n"
            "            co[b, c] = cbc + ct_eff",
            1,
        ),
    }
    generated = []
    diffs = []
    for name in ("full_renter_block_kernel", "full_owner_block_kernel"):
        original = getattr(model, name)
        python_function = getattr(original, "py_func", original)
        source = inspect.getsource(python_function)
        before, after, expected = replacements[name]
        if source.count(before) != expected:
            raise ValueError(f"allocation patch anchor count != {expected}: {name}")
        revised = source.replace(before, after)
        revised = _replace_once(
            revised, "@njit(cache=True, parallel=True)",
            "@njit(cache=False, parallel=True)")
        namespace = dict(python_function.__globals__)
        exec(compile(revised, str(diff_path.with_suffix(".allocation.py")), "exec"), namespace)
        setattr(model, name, namespace[name])
        generated.append(revised)
        diffs.extend(difflib.unified_diff(
            source.splitlines(True), revised.splitlines(True),
            fromfile=f"frozen/{name}", tofile=f"diagnostic/{name}.allocation"))
    diff_path.with_suffix(".allocation.py").write_text("\n".join(generated))
    return "".join(diffs)


def install_purchase_income(model: Any, diff_path: Path):
    """Apply the reviewed timing change only to the supported Markov branch.

    x = b + sale proceeds - purchase price is unchanged. Conditional spending
    remains c + owner costs + b' = R*x + y. The purchase threshold admits y/R,
    while the final mortgage floor is still b' >= -phi*pH, without rollover of
    the temporary purchase-stage shortfall. No income is added to x itself.
    """
    original = inspect.getsource(model.solve_bellman_full_markov_income)
    revised = _replace_once(original, "    t0 = time.perf_counter()\n", '''    t0 = time.perf_counter()
    if (bool(getattr(P, "joint_nested_choice", False))
            or bool(getattr(P, "use_pti_constraint", False))
            or float(P.lambda_d) != 0.0
            or not bool(getattr(P, "use_tenure_kernel", True))
            or not bool(getattr(P, "use_full_kernel", True))
            or not NUMBA_AVAILABLE
            or str(getattr(P, "interp_method", "linear")) != "linear"
            or np.any(SD.birth_dp) or np.any(SD.birth_entry_grant)):
        raise ValueError("unsupported mechanism combined with purchase-income diagnostic")
''')
    revised = _replace_once(revised, '''                            alpha, oms, beta, s_next, D_next, gs_alpha1, gs_alpha2, gs_tol,
                            strict_owner_hbar_feasibility, int(exhaustive_saving),''', '''                            alpha, oms, beta, 0.0, 0.0, gs_alpha1, gs_alpha2, gs_tol,
                            strict_owner_hbar_feasibility, int(exhaustive_saving),''')
    revised = _replace_once(revised, "            dp_choice = dp_arr\n", '''            income_for_purchase = np.array([
                income_at_state(P, i, j, float(z_value)) for i in range(I)
            ], dtype=float).reshape(I, 1, 1, 1) / Rg
            dp_choice = dp_arr - income_for_purchase
            # Require actual post-transaction wealth to lie on the solved grid.
            bmo_purchase = np.maximum(bmo - income_for_purchase, b_grid[0])
''')
    # Only the two compiled tenure calls are enabled by the guard above.
    before = "Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant"
    if revised.count(before) != 2:
        raise ValueError("unexpected compiled tenure call sites")
    revised = revised.replace(before, before.replace("dp_choice, bmo,", "dp_choice, bmo_purchase,"))
    diff = "".join(difflib.unified_diff(original.splitlines(True), revised.splitlines(True),
                                      fromfile="frozen/solve_bellman_full_markov_income",
                                      tofile="diagnostic/solve_bellman_full_markov_income"))
    diff_path.parent.mkdir(parents=True, exist_ok=True)
    diff += _install_transaction_support(model, diff_path)
    diff += _install_exact_allocation_output(model, diff_path)
    diff_path.write_text(diff)
    diff_path.with_suffix(".generated.py").write_text(revised)
    # Use the original module globals so every calendar/stationary caller sees
    # the same function; source files themselves remain byte-for-byte frozen.
    exec(compile(revised, str(diff_path.with_suffix(".generated.py")), "exec"), model.__dict__)
    model.build_forward_tenure_transition_maps = _ordinary_forward_maps
    return diff


def audit_purchase_accounting(evaluation, P, shared, grid, model):
    """Check occupied transactions and final debt against independent formulas.

    This supplements, and does not replace, the native dated household budget
    gate. The support test catches interpolation clipping before accepting a
    scored solution. End-of-period debt cannot retain the purchase shortfall.
    """
    if int(P.I) != 1 or evaluation.policy.tenure_probs is None:
        raise ValueError("purchase audit requires one market and probabilistic tenure")
    policy = evaluation.policy
    grid = np.asarray(grid, dtype=float)
    price = float(policy.price[0])
    costs = np.r_[0., price * np.asarray(P.H_own)]
    sale = (1. - float(P.psi)) * costs
    outside_mass = threshold_mass = debt_mass = 0.
    largest_wealth_error = largest_debt_shortfall = 0.
    for age in range(P.J):
        for zz, z in enumerate(P.z_grid):
            y = model.income_at_state(P, 0, age, float(z))
            for old in range(len(costs)):
                mass = (evaluation.g_post_fertility[:, old, 0, age, zz]
                        * policy.loc_probs[:, old, 0, 0, age, zz])
                probs = np.asarray(policy.tenure_probs[:, old, 0, age, zz], dtype=float)
                sums = probs.sum(axis=-1, keepdims=True)
                probs = np.divide(probs, sums, out=np.zeros_like(probs), where=sums > 0)
                for new in range(len(costs)):
                    branch_mass = mass * probs[..., new]
                    x = grid if old == new else grid + sale[old] - costs[new]
                    outside = (x < grid[0] - 1e-12) | (x > grid[-1] + 1e-12)
                    outside_mass += float(branch_mass[outside].sum())
                    for nn in range(P.n_parity):
                        for cs in range(P.n_child_states):
                            occupied = branch_mass[:, nn, cs] > 1e-12
                            idx = policy.maps.tmx_idx[0, old, new, nn, cs]
                            wt = policy.maps.tmx_wt[0, old, new, nn, cs]
                            mapped = (1. - wt) * grid[idx] + wt * grid[idx + 1]
                            if np.any(occupied):
                                largest_wealth_error = max(largest_wealth_error,
                                    float(np.max(np.abs(mapped[occupied] - x[occupied]))))
                            if new > 0 and new != old:
                                floor = -float(shared.phi_choice[0, new, nn, cs]) * costs[new]
                                invalid = x + y / float(P.R_gross) < floor - 1e-10
                                threshold_mass += float(branch_mass[invalid, nn, cs].sum())
            for new in range(1, len(costs)):
                final = policy.bp_pol[:, new, 0, age, zz]
                mass = evaluation.g_current[:, new, 0, age, zz]
                floor = -np.asarray(shared.phi_choice[0, new]) * costs[new]
                shortfall = floor[None, :, :] - final
                debt_mass += float(mass[shortfall > 1e-9].sum())
                occupied = mass > 1e-12
                if np.any(occupied):
                    largest_debt_shortfall = max(largest_debt_shortfall,
                        float(shortfall[occupied].max()))
    receipt = dict(transaction_outside_grid_mass=outside_mass,
                   purchase_threshold_violation_mass=threshold_mass,
                   end_mortgage_floor_violation_mass=debt_mass,
                   maximum_occupied_transaction_wealth_error=largest_wealth_error,
                   maximum_occupied_end_debt_shortfall=largest_debt_shortfall,
                   mass_tolerance=2e-10, wealth_tolerance=1e-9)
    if (max(outside_mass, threshold_mass, debt_mass) > 2e-10
            or max(largest_wealth_error, largest_debt_shortfall) > 1e-9):
        raise RuntimeError(f"purchase accounting gate failed: {receipt}")
    return receipt
