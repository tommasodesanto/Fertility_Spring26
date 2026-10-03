"""Distribution stage of the extracted solver (bodies byte-identical; see split_receipt.json)."""
from __future__ import annotations
import copy
import math
import time
from types import SimpleNamespace
from typing import Any
import numpy as np
from . import joint_nested, birth_count
from .adult_entry import adjusted_births, potential_entry_households
from .warm_price import search_warm_price
from .child_preferences import apply_child_preferences
from .parameters import (bequest_utility_net_active, child_earnings_multiplier, child_earnings_penalty_active, children_at_home_count, estate_flow_net_active, estate_housing_value, estate_receiver_active, estate_transfer_at_age, get_fecundity_by_age, independent_child_maturation_active, mortgage_stay_floor_active, parent_age_maturation_active, readiness_childless_states, readiness_cumulative_probability, readiness_gate_active, readiness_settled_state, readiness_transition_hazard, rental_wedge_active, unsecured_debt_floor)
from .kernels import (NUMBA_AVAILABLE, full_owner_block_kernel, full_renter_block_kernel, location_logit_kernel, scatter_cols_kernel, scatter_cols_sameidx_kernel, scatter_vec_kernel, tenure_choice_kernel, tenure_logit_kernel)
from .utils import (decode_flat_family_state, flat_nc, interp_indices, interp_on_grid, interp_vector, logsumexp, make_grid, make_value_interp, scatter_redistribute, scatter_redistribute_cols, scatter_redistribute_cols_sameidx, unflat_nc, weighted_median_from_cells, weighted_quantile)
from .shared import (DEAD_MASS_TOL, DEAD_VALUE_CUTOFF, ENTRY_WEALTH_INCOME_RATIO_MODES, InfeasibleThetaError, _linear_grid_weights_for_points, add_aggregate_wealth_gross_labor_diagnostics, annual_gross_income_at_state, housing_demand_normalizer, income_at_state, income_transition_values, markov_grant_outlays, normalize_population_mass, penalized_income_at_state, property_tax_revenue_from_distribution)
from .household import (_blended_child_Pi_for_cell, _child_Pa_exempt_for_age, _child_Pa_for_age, birth_destination_child_state, build_forward_tenure_transition_maps, current_child_bin_dt, effective_owner_collateral_floor, owner_borrowing_floor, renter_borrowing_floor)


def uses_entry_wealth_income_ratio(P: SimpleNamespace) -> bool:
    return str(getattr(P, "entry_wealth_mode", "scalar")).lower() in ENTRY_WEALTH_INCOME_RATIO_MODES


def entry_wealth_ratio_distribution(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    ratios = np.asarray(getattr(P, "entry_wealth_ratio_nodes", []), dtype=float).reshape(-1)
    weights = np.asarray(getattr(P, "entry_wealth_ratio_weights", []), dtype=float).reshape(-1)
    if ratios.size == 0:
        raise ValueError("entry_wealth_mode uses income ratios but entry_wealth_ratio_nodes is empty.")
    if weights.size == 0:
        weights = np.ones(ratios.size, dtype=float) / float(ratios.size)
    if weights.size != ratios.size:
        raise ValueError("entry_wealth_ratio_weights must have the same length as entry_wealth_ratio_nodes.")
    weights = np.maximum(weights, 0.0)
    if float(np.sum(weights)) <= 0.0:
        weights = np.ones(ratios.size, dtype=float) / float(ratios.size)
    else:
        weights = weights / float(np.sum(weights))
    return ratios, weights


def _scalar_entry_wealth_grid_weights(b_grid: np.ndarray, P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    """Grid-node distribution for entrant liquid wealth.

    The default is the historical point-mass injection at the first grid node
    weakly above ``b_entry_fixed``. If ``entry_wealth_spread_nodes`` is larger
    than one, entrants are distributed over the nearest grid nodes with weights
    exponentially tilted to preserve mean wealth at ``b_entry_fixed``.
    """
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    if bg.size == 0:
        raise ValueError("b_grid must contain at least one node")
    b_entry = float(np.clip(getattr(P, "b_entry_fixed", 0.0), bg[0], bg[-1]))
    spread_nodes = max(int(getattr(P, "entry_wealth_spread_nodes", 1)), 1)
    if spread_nodes <= 1 or bg.size == 1:
        idx = int(np.argmax(bg >= b_entry)) if np.any(bg >= b_entry) else int(bg.size - 1)
        return np.array([idx], dtype=np.int64), np.array([1.0], dtype=float)

    n_nodes = min(spread_nodes, int(bg.size))
    nearest = np.argsort(np.abs(bg - b_entry), kind="mergesort")[:n_nodes]
    idx = np.sort(nearest).astype(np.int64)
    x = bg[idx]
    if b_entry <= x[0] + 1e-12:
        wt = np.zeros_like(x, dtype=float)
        wt[0] = 1.0
        return idx, wt
    if b_entry >= x[-1] - 1e-12:
        wt = np.zeros_like(x, dtype=float)
        wt[-1] = 1.0
        return idx, wt

    scale = max(float(x[-1] - x[0]), 1e-12)
    x_scaled = (x - b_entry) / scale

    def tilted(lam: float) -> tuple[np.ndarray, float]:
        z = np.clip(lam * x_scaled, -700.0, 700.0)
        z = z - float(np.max(z))
        wt0 = np.exp(z)
        wt0 = wt0 / float(np.sum(wt0))
        return wt0, float(np.sum(wt0 * x))

    lo = -1.0
    hi = 1.0
    _, mlo = tilted(lo)
    _, mhi = tilted(hi)
    for _ in range(80):
        if mlo <= b_entry <= mhi:
            break
        if mlo > b_entry:
            hi = lo
            lo *= 2.0
            _, mlo = tilted(lo)
        elif mhi < b_entry:
            lo = hi
            hi *= 2.0
            _, mhi = tilted(hi)

    wt = np.ones_like(x, dtype=float) / float(x.size)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        wt_mid, mean_mid = tilted(mid)
        wt = wt_mid
        if mean_mid < b_entry:
            lo = mid
        else:
            hi = mid
    wt = wt / float(np.sum(wt))
    return idx, wt


def entry_wealth_grid_weights(
    b_grid: np.ndarray,
    P: SimpleNamespace,
    *,
    i: int = 0,
    j: int = 0,
    z_value: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Grid-node distribution for entrant liquid wealth.

    In the legacy scalar mode this preserves the historical `b_entry_fixed`
    injection. In `income_ratio_distribution` mode, entrants draw empirical
    wealth/income ratios and those ratios are converted to model wealth using
    annual gross income at the entrant state.
    """
    if bool(getattr(P, "native_fixed_reference_entry", False)):
        if i != 0 or j != 0 or int(P.I) != 1:
            raise ValueError("fixed entry mapping requires age-zero in one market")
        np.testing.assert_array_equal(b_grid, P.fixed_reference_entry_grid)
        zz = np.flatnonzero(np.asarray(P.z_grid) == z_value)
        if len(zz) != 1:
            raise ValueError("entry lookup requires a unique exact income node")
        conditional = np.asarray(P.fixed_reference_entry_conditional, dtype=float)
        if (conditional.shape != (len(b_grid), len(P.z_grid))
                or not np.isfinite(conditional).all() or np.any(conditional < 0)
                or not np.allclose(conditional.sum(axis=0), 1.0, rtol=0, atol=2e-12)):
            raise ValueError("invalid fixed conditional entry distribution")
        weights = conditional[:, int(zz[0])]
        idx = np.flatnonzero(weights > 0)
        return idx, weights[idx]
    if not uses_entry_wealth_income_ratio(P):
        return _scalar_entry_wealth_grid_weights(b_grid, P)
    ratios, weights = entry_wealth_ratio_distribution(P)
    y_entry = annual_gross_income_at_state(P, int(i), int(j), float(z_value))
    points = ratios * y_entry
    return _linear_grid_weights_for_points(b_grid, points, weights)


def aggregate_entry_wealth_grid_weights(b_grid: np.ndarray, P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    """Unconditional entrant wealth distribution over grid nodes for reporting."""
    if not uses_entry_wealth_income_ratio(P) and not bool(getattr(P, "native_fixed_reference_entry", False)):
        return entry_wealth_grid_weights(b_grid, P)
    z_grid, z_weights, _ = income_transition_values(P)
    entry_by_loc = np.maximum(np.asarray(getattr(P, "entry_by_loc", np.ones(P.I)), dtype=float).reshape(-1), 0.0)
    if entry_by_loc.size != int(P.I) or float(np.sum(entry_by_loc)) <= 0.0:
        entry_by_loc = np.ones(int(P.I), dtype=float) / max(int(P.I), 1)
    else:
        entry_by_loc = entry_by_loc / float(np.sum(entry_by_loc))
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    mass = np.zeros(bg.size, dtype=float)
    for i in range(int(P.I)):
        for zz, wz in enumerate(z_weights):
            idx, wt = entry_wealth_grid_weights(bg, P, i=i, j=0, z_value=float(z_grid[zz]))
            mass[idx] += float(entry_by_loc[i]) * float(wz) * wt
    total = float(np.sum(mass))
    if total <= 0.0:
        return entry_wealth_grid_weights(bg, P)
    idx = np.flatnonzero(mass > 0.0).astype(np.int64)
    return idx, mass[idx] / total


def attach_entry_wealth_stats(
    stats: SimpleNamespace,
    b_grid: np.ndarray,
    entry_idx: np.ndarray,
    entry_wt: np.ndarray,
    P: SimpleNamespace,
) -> SimpleNamespace:
    x = np.asarray(b_grid, dtype=float).reshape(-1)[np.asarray(entry_idx, dtype=np.int64)]
    wt = np.asarray(entry_wt, dtype=float).reshape(-1)
    stats.entry_wealth_spread_nodes = int(getattr(P, "entry_wealth_spread_nodes", 1))
    stats.entry_wealth_grid_indices = np.asarray(entry_idx, dtype=np.int64).copy()
    stats.entry_wealth_grid_values = x.copy()
    stats.entry_wealth_grid_weights = wt.copy()
    stats.entry_wealth_grid_mean = float(np.sum(x * wt))
    stats.entry_wealth_target = float(np.clip(getattr(P, "b_entry_fixed", 0.0), x[0], x[-1])) if x.size else np.nan
    stats.entry_wealth_mode = str(getattr(P, "entry_wealth_mode", "scalar"))
    if uses_entry_wealth_income_ratio(P):
        ratios, ratio_wt = entry_wealth_ratio_distribution(P)
        stats.entry_wealth_ratio_nodes = ratios.copy()
        stats.entry_wealth_ratio_weights = ratio_wt.copy()
        stats.entry_wealth_ratio_mean = float(np.sum(ratios * ratio_wt))
        stats.entry_wealth_ratio_source = str(getattr(P, "entry_wealth_ratio_source", ""))
        stats.b_entry_fixed_legacy = float(getattr(P, "b_entry_fixed", np.nan))
        stats.entry_wealth_target = stats.entry_wealth_grid_mean
    return stats


def upgrade_fast_markov_solution(
    fast_solution: SimpleNamespace, P: SimpleNamespace, b_grid: np.ndarray, SD: SimpleNamespace,
) -> SimpleNamespace:
    """Build full reporting objects by reusing the accepted Bellman payload."""
    payload = getattr(fast_solution, "_model_payload", None)
    if payload is None:
        raise ValueError("fast solution does not retain a household-solution payload")
    V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, fert2_probs = payload
    # Must precede the KFE: forward_distribution reads P._fert2_probs only as
    # an input; its other P._* fertility fields are KFE outputs.
    P._fert2_probs = fert2_probs
    if birth_count.enabled(P):
        P.birth_count_action_probs = fast_solution.birth_count_action_probs.copy()
        P.birth_count_realized_probs = fast_solution.birth_count_realized_probs.copy()
        P.birth_count_policy_axes = fast_solution.birth_count_policy_axes
    P._joint_choice = getattr(fast_solution, "joint_choice", None)
    P._bp_pol_stay = getattr(fast_solution, "_bp_pol_stay", None)
    P._c_pol_stay = getattr(fast_solution, "_c_pol_stay", None)
    start = time.perf_counter()
    g, stats = forward_distribution_markov_income(
        bp_pol, hR_pol, tc, lp_j, fp, V, r, p, P, b_grid, SD,
        fast_stats=False, tenure_probs=tp, bp_pol_stay=P._bp_pol_stay,
    )
    solution = pack_solution_markov_income(
        V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, g, stats, P.w_hat, p, P,
    )
    solution.b_grid = np.asarray(b_grid, dtype=float).copy()
    solution.timings = {
        "bellman_full": 0.0, "distribution": time.perf_counter() - start,
        "n_full": 0, "n_eval": 0, "n_dist": 1, "income_process": "markov",
        "bellman_mode": "reused_accepted_price",
        "reused_bellman_seconds": float(getattr(fast_solution, "timings", {}).get("bellman_full", 0.0)),
    }
    return solution


def realize_current_choices(
    cohort: np.ndarray,
    j: int,
    loc_probs: np.ndarray,
    tenure_choice: np.ndarray,
    tenure_probs: np.ndarray | None,
    lmm_idx: np.ndarray,
    lmm_wt: np.ndarray,
    tmx_idx: np.ndarray,
    tmx_wt: np.ndarray,
    *,
    use_compiled_scatter: bool = False,
    mass_pruning_tolerance: float = 1e-15,
) -> np.ndarray:
    """Apply location, tenure, and housing transactions without aging.

    ``cohort`` is the age-j mass after any current family-size choice. The
    returned mass is indexed by realized current location, tenure/rung, and
    post-transaction liquid wealth. Saving, income, family-stage, and age
    transitions are deliberately excluded.
    """

    Nb, nt, I, npar, ncs = cohort.shape
    nc = npar * ncs
    after_location = np.zeros_like(cohort)
    for io in range(I):
        for to in range(nt):
            origin = flat_nc(cohort[:, to, io, :, :], Nb, nc)
            if np.sum(origin) == 0.0 or np.sum(origin) < mass_pruning_tolerance:
                continue
            probs = np.reshape(loc_probs[:, to, io, :, j, :, :], (Nb, I, nc), order="F")
            after_location[:, to, io, :, :] += unflat_nc(origin * probs[:, io, :], Nb, npar, ncs)
            idx = lmm_idx[io, to, :]
            wt = lmm_wt[io, to, :]
            for id_ in range(I):
                if id_ == io:
                    continue
                moved = origin * probs[:, id_, :]
                if use_compiled_scatter:
                    moved = scatter_cols_sameidx_kernel(idx, wt, moved, Nb)
                else:
                    moved = scatter_redistribute_cols_sameidx(idx, wt, moved, Nb, mass_pruning_tolerance=mass_pruning_tolerance)
                after_location[:, 0, id_, :, :] += unflat_nc(moved, Nb, npar, ncs)

    realized = np.zeros_like(cohort)
    for nn in range(npar):
        for id_ in range(I):
            for to in range(nt):
                source = after_location[:, to, id_, nn, :]
                if np.sum(source) == 0.0 or np.sum(source) < mass_pruning_tolerance:
                    continue
                normalized_probs = None
                if tenure_probs is not None:
                    all_probs = np.asarray(tenure_probs[:, to, id_, j, nn, :, :], dtype=float)
                    prob_sum = np.sum(all_probs, axis=-1)
                    normalized_probs = np.divide(
                        all_probs,
                        prob_sum[:, :, None],
                        out=np.zeros_like(all_probs),
                        where=prob_sum[:, :, None] > 0.0,
                    )
                for tn in range(nt):
                    if tenure_probs is None:
                        selected = tenure_choice[:, to, id_, j, nn, :] == tn
                        mass = source * selected
                    else:
                        mass = source * normalized_probs[:, :, tn]
                    if np.sum(mass) == 0.0 or np.sum(mass) < mass_pruning_tolerance:
                        continue
                    redistributed = np.zeros((Nb, ncs))
                    for cs in range(ncs):
                        idx = tmx_idx[id_, to, tn, nn, cs, :]
                        wt = tmx_wt[id_, to, tn, nn, cs, :]
                        if use_compiled_scatter:
                            redistributed[:, cs] = scatter_vec_kernel(idx, wt, mass[:, cs], Nb)
                        else:
                            redistributed[:, cs] = scatter_redistribute(idx, wt, mass[:, cs], Nb)
                    realized[:, tn, id_, nn, :] += redistributed
    return realized


def realize_current_choices_markov_income(
    cohort: np.ndarray,
    j: int,
    loc_probs: np.ndarray,
    tenure_choice: np.ndarray,
    tenure_probs: np.ndarray | None,
    lmm_idx: np.ndarray,
    lmm_wt: np.ndarray,
    tmx_idx: np.ndarray,
    tmx_wt: np.ndarray,
    *,
    use_compiled_scatter: bool = False,
    mass_pruning_tolerance: float = 1e-15,
) -> np.ndarray:
    """Markov-income counterpart of :func:`realize_current_choices`."""

    Nb, nt, I, Nz, npar, ncs = cohort.shape
    realized = np.zeros_like(cohort)
    for zz in range(Nz):
        realized[:, :, :, zz, :, :] = realize_current_choices(
            cohort[:, :, :, zz, :, :],
            j,
            loc_probs[:, :, :, :, :, zz, :, :],
            tenure_choice[:, :, :, :, zz, :, :],
            None if tenure_probs is None else tenure_probs[:, :, :, :, zz, :, :, :],
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            use_compiled_scatter=use_compiled_scatter,
            mass_pruning_tolerance=mass_pruning_tolerance,
        )
    return realized


def realize_stayer_cross_section(g, loc_probs, tenure_choice, tenure_probs):
    """Owner stayers after location/tenure choices, before saving (7D mass).

    Preserve origin tenure: new buyers with the same destination state are not
    eligible for inherited-debt treatment. Moving location liquidates ownership.
    """
    if g.ndim != 7:
        raise ValueError("stayer mass requires income-resolved seven-dimensional mass")
    out = np.zeros_like(g)
    for i in range(g.shape[2]):
        for ten in range(1, g.shape[1]):
            if tenure_probs is None:
                pr = tenure_choice[:, ten, i] == ten
            else:
                menu = np.asarray(tenure_probs[:, ten, i], dtype=float)
                total = menu.sum(axis=-1)
                pr = np.divide(menu[..., ten], total, out=np.zeros_like(total), where=total > 0)
            out[:, ten, i] = g[:, ten, i] * loc_probs[:, ten, i, i] * pr
    return out


def realize_current_cross_section(
    g: np.ndarray,
    loc_probs: np.ndarray,
    tenure_choice: np.ndarray,
    tenure_probs: np.ndarray | None,
    lmm_idx: np.ndarray,
    lmm_wt: np.ndarray,
    tmx_idx: np.ndarray,
    tmx_wt: np.ndarray,
    *,
    use_compiled_scatter: bool = False,
    mass_pruning_tolerance: float = 1e-15,
) -> np.ndarray:
    """Build the realized current cross-section from beginning-of-period mass."""

    out = np.zeros_like(g)
    if g.ndim == 6:
        for j in range(g.shape[3]):
            out[:, :, :, j, :, :] = realize_current_choices(
                g[:, :, :, j, :, :],
                j,
                loc_probs,
                tenure_choice,
                tenure_probs,
                lmm_idx,
                lmm_wt,
                tmx_idx,
                tmx_wt,
                use_compiled_scatter=use_compiled_scatter,
                mass_pruning_tolerance=mass_pruning_tolerance,
            )
    elif g.ndim == 7:
        for j in range(g.shape[3]):
            out[:, :, :, j, :, :, :] = realize_current_choices_markov_income(
                g[:, :, :, j, :, :, :],
                j,
                loc_probs,
                tenure_choice,
                tenure_probs,
                lmm_idx,
                lmm_wt,
                tmx_idx,
                tmx_wt,
                use_compiled_scatter=use_compiled_scatter,
                mass_pruning_tolerance=mass_pruning_tolerance,
            )
    else:
        raise ValueError(f"unsupported distribution rank: {g.ndim}")
    return out


def _dead_mass_census_at_age(
    g_age: np.ndarray,
    V_age: np.ndarray,
    j: int,
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    *,
    markov_income: bool,
) -> tuple[float, list[dict[str, Any]]]:
    """Return positive mass on Bellman-dead nodes and a compact state census."""

    ga = np.asarray(g_age, dtype=float)
    va = np.asarray(V_age, dtype=float)
    if not markov_income:
        ga = ga[:, :, :, None, :, :]
        va = va[:, :, :, None, :, :]
        z_grid = np.array([1.0])
    else:
        z_grid, _, _ = income_transition_values(P)
    dead_mass_array = np.where(va <= DEAD_VALUE_CUTOFF, ga, 0.0)
    dead_mass = float(np.sum(dead_mass_array))
    if dead_mass <= DEAD_MASS_TOL:
        return dead_mass, []

    census: list[dict[str, Any]] = []
    positive = np.argwhere(dead_mass_array > 0.0)
    phi_state = SD.phi_state
    for b_idx, ten, i, zz, nn, cs in positive[:8]:
        b_now = float(b_grid[b_idx])
        z_value = float(z_grid[zz])
        y_now = penalized_income_at_state(
            P, int(i), int(j), z_value, children_at_home_count(int(nn), int(cs), P)
        )
        resources = float(P.R_gross) * b_now + y_now
        gG = float(SD.g_bar[nn, cs]) if hasattr(SD, "g_bar") else 0.0
        x_test = float(P.R_gross) * max(b_now, 0.0) + y_now
        tr = min(max(gG - x_test, 0.0), gG)
        resources = resources + tr
        if ten == 0:
            floor = float(renter_borrowing_floor(P, b_now, j))
            required_flow = float(SD.c_bar[nn, cs]) + float(r_hat[i]) * float(SD.h_bar[nn, cs])
            unsecured = b_now
        else:
            house = float(P.H_own[ten - 1])
            collateral_floor = -float(phi_state[nn, cs]) * float(p_hat[i]) * house
            effective_collateral = float(effective_owner_collateral_floor(P, collateral_floor, j))
            floor = float(owner_borrowing_floor(P, b_now, collateral_floor, j))
            extra_size_cost = float(getattr(P, "owner_size_cost", 0.0)) * float(p_hat[i]) * max(
                house - float(getattr(P, "owner_size_cost_ref", 6.0)), 0.0
            ) ** float(getattr(P, "owner_size_cost_power", 2.0))
            required_flow = (
                (float(P.delta) + float(P.tau_H)) * float(p_hat[i]) * house
                + extra_size_cost
                + float(SD.c_bar[nn, cs])
            )
            unsecured = b_now - effective_collateral
        census.append(
            {
                "age": float(P.age_start + j * P.da),
                "b": b_now,
                "child_state": int(cs),
                "income": float(y_now),
                "location": int(i),
                "mass": float(dead_mass_array[b_idx, ten, i, zz, nn, cs]),
                "parity": int(nn),
                "slack": float(resources - required_flow - floor),
                "tenure": int(ten),
                "transfer": tr,
                "unsecured_position": float(unsecured),
                "z": z_value,
            }
        )
    return dead_mass, census


def _gate_dead_mass_at_age(
    g_age: np.ndarray,
    V_age: np.ndarray,
    j: int,
    stage: str,
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    *,
    markov_income: bool,
) -> None:
    dead_mass, census = _dead_mass_census_at_age(
        g_age,
        V_age,
        j,
        r_hat,
        p_hat,
        P,
        b_grid,
        SD,
        markov_income=markov_income,
    )
    if dead_mass > DEAD_MASS_TOL:
        raise InfeasibleThetaError(stage, dead_mass, census)


def _censor_entry_dead_mass(g_age: np.ndarray, V_age: np.ndarray) -> float:
    """Relocate entrant mass on Bellman-dead nodes to the feasible frontier.

    For each entry state column, mass on a dead b-node moves to the nearest
    alive node with weakly higher b in the same column, so no mass ever sits
    on a dead node (the July-11 gate invariant is preserved and the gate
    itself stays unchanged as a backstop). Columns with no alive node above
    are left untouched so the gate still rejects them. Operates in place on
    the age-0 view; returns the relocated mass.
    """

    dead = (V_age <= DEAD_VALUE_CUTOFF) & (g_age > 0.0)
    if not np.any(dead):
        return 0.0
    moved = 0.0
    for idx in np.argwhere(dead):
        b_idx, rest = int(idx[0]), tuple(int(k) for k in idx[1:])
        column_alive = V_age[(slice(None),) + rest] > DEAD_VALUE_CUTOFF
        above = np.nonzero(column_alive[b_idx:])[0]
        if above.size == 0:
            continue
        target = b_idx + int(above[0])
        mass = float(g_age[(b_idx,) + rest])
        g_age[(target,) + rest] += mass
        g_age[(b_idx,) + rest] = 0.0
        moved += mass
    return moved


def forward_distribution_markov_income(
    bp_pol: np.ndarray,
    hR_pol: np.ndarray,
    tenure_choice: np.ndarray,
    loc_probs: np.ndarray,
    fert_probs: np.ndarray,
    state_values: np.ndarray,
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    fast_stats: bool = False,
    tenure_probs: np.ndarray | None = None,
    bp_pol_stay: np.ndarray | None = None,
) -> tuple[np.ndarray, SimpleNamespace]:
    fec = get_fecundity_by_age(P)
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    z_grid, z_weights, Pi_z = income_transition_values(P)
    Nz = len(z_grid)
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    bmin = b_grid[0]
    bmax = b_grid[-1]
    g = np.zeros((Nb, nt, I, J, Nz, npar, ncs))
    count_active = birth_count.enabled(P)
    if count_active:
        birth_count.validate_contract(P)
        count_pre = np.zeros_like(g)
        count_post = np.zeros_like(g)
        count_first_tagged = np.zeros_like(g)
        count_risk = np.zeros((3, J))
        count_attempts = np.zeros((3, J))
        count_any_births = np.zeros(J)
    entry_idx, entry_wt = aggregate_entry_wealth_grid_weights(b_grid, P)
    for i in range(I):
        for zz in range(Nz):
            loc_entry_idx, loc_entry_wt = entry_wealth_grid_weights(b_grid, P, i=i, j=0, z_value=float(z_grid[zz]))
            settled_entry_share = (
                readiness_cumulative_probability(P, float(P.age_start))
                if readiness_gate_active(P)
                else 0.0
            )
            for kk, ww in zip(loc_entry_idx, loc_entry_wt):
                entrant_mass = float(ww) * P.entry_by_loc[i] * z_weights[zz]
                g[int(kk), 0, i, 0, zz, 0, 0] += (
                    (1.0 - settled_entry_share) * entrant_mass
                )
                if readiness_gate_active(P):
                    g[int(kk), 0, i, 0, zz, 0, 1] += (
                        settled_entry_share * entrant_mass
                    )
    joint_pre = np.zeros_like(g) if bool(getattr(P, "joint_nested_choice", False)) and not fast_stats else None
    P._entry_censored_mass = 0.0
    P._entry_total_mass = float(np.sum(g[:, :, :, 0, :, :, :]))
    if bool(getattr(P, "entry_wealth_censor_to_frontier", False)):
        P._entry_censored_mass = _censor_entry_dead_mass(
            g[:, :, :, 0, :, :, :], state_values[:, :, :, 0, :, :, :]
        )
    _gate_dead_mass_at_age(
        g[:, :, :, 0, :, :, :],
        state_values[:, :, :, 0, :, :, :],
        0,
        "entry",
        r_hat,
        p_hat,
        P,
        b_grid,
        SD,
        markov_income=True,
    )
    total_births = 0.0
    births_by_loc = np.zeros(I)
    first_births_by_age = np.zeros(J)
    second_births_by_age = np.zeros(J)
    second_attempts_by_age = np.zeros(J)
    second_at_risk_by_age = np.zeros(J)
    third_births_by_age = np.zeros(J)
    third_attempts_by_age = np.zeros(J)
    third_at_risk_by_age = np.zeros(J)
    entrants_mature_by_loc = np.zeros(I)
    entrants_mature_total = 0.0
    # Matured children per new entrant household (0.5 under literal parity:
    # two children pair into one next-generation household). 1.0 nests.
    ecf = float(getattr(P, "entrant_conversion_factor", 1.0))
    # Event-study horizon in four-year model periods.  The generic production
    # default remains zero.  The active E5 transition profile overrides it to
    # one because its PSID calibration row is the four-year change from event
    # time -1 to +3.  The statistic is the birth branch minus an otherwise
    # identical no-birth branch after the declared horizon.
    event_horizon = int(getattr(P, "housing_event_horizon", 0))
    birth_es3_pre_sum = birth_es3_post_sum = birth_es3_mass = 0.0
    birth_es3_control_post_sum = 0.0
    addchild_es3_one_sum = addchild_es3_one_mass = 0.0
    addchild_es3_two_plus_sum = addchild_es3_two_plus_mass = 0.0
    onechild_es3_pre_sum = onechild_es3_post_sum = onechild_es3_mass = 0.0
    twoplus_es3_pre_sum = twoplus_es3_post_sum = twoplus_es3_mass = 0.0

    hc = np.zeros((I, nt))
    he = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hc[i, ten] = p_hat[i] * hs
            he[i, ten] = (1 - P.psi) * p_hat[i] * hs

    phi_choice = SD.phi_choice
    lmm_idx = np.zeros((I, nt, Nb), dtype=np.int64)
    lmm_wt = np.zeros((I, nt, Nb))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + he[io, to], bmin, bmax)
            lmm_idx[io, to, :], lmm_wt[io, to, :] = interp_indices(b_grid, ba)

    tmx_idx, tmx_wt = build_forward_tenure_transition_maps(
        P, b_grid, hc, he, phi_choice, SD.birth_dp, SD.birth_entry_grant
    )

    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2
    ust = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    Pia = P.Pi_child if ust else None

    for j in range(J - 1):
        # Parent-age m-d bookkeeping: post-birth newborn inflow per
        # (parity, at-home) cell at age j.  Stays zero in constant mode.
        newborn_inflow = np.zeros((npar, ncs))
        _gate_dead_mass_at_age(
            g[:, :, :, j, :, :, :],
            state_values[:, :, :, j, :, :, :],
            j,
            f"forward_age_{P.age_start + j * P.da:g}",
            r_hat,
            p_hat,
            P,
            b_grid,
            SD,
            markov_income=True,
        )
        if count_active:
            count_pre[:, :, :, j] = g[:, :, :, j]
        if count_active and (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end):
            flow = birth_count.transition_at_age(g[:, :, :, j].copy(), P, j,
                state_values=state_values[:, :, :, j])
            g[:, :, :, j] = flow["post"]
            count_first_tagged[:, :, :, j] = flow["first_birth_tagged_post"]
            first_births_by_age[j], second_births_by_age[j], third_births_by_age[j] = flow["births_by_order"]
            count_risk[:, j] = flow["at_risk_by_order"]
            count_attempts[:, j] = flow["attempts_by_order"]
            count_any_births[j] = flow["any_birth_mass"]
            second_at_risk_by_age[j], third_at_risk_by_age[j] = flow["at_risk_by_order"][1:]
            second_attempts_by_age[j], third_attempts_by_age[j] = flow["attempts_by_order"][1:]
            total_births += flow["expected_births"]
            births_by_loc += np.sum(flow["expected_births_by_cell"], axis=(0, 1, 3))
            for zz in range(Nz):
                birth_mass = np.sum(count_first_tagged[:, :, :, j, zz], axis=(-2, -1))
                birth_mass_total = float(np.sum(birth_mass))
                if not fast_stats and birth_mass_total > 1e-12 and (j + event_horizon) < J:
                    birth_weight = np.zeros((Nb, nt, I, Nz))
                    birth_weight[:, :, :, zz] = birth_mass
                    pre_h = mean_housing_childless_weighted_markov(birth_weight, j, hR_pol, P)
                    birth_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                    birth_cohort[:, :, :, zz, :, :] = count_first_tagged[:, :, :, j, zz, :, :]
                    birth_cohort = advance_cohort_horizon_markov_income(
                        birth_cohort,
                        j,
                        event_horizon,
                        loc_probs,
                        tenure_choice,
                        tenure_probs,
                        bp_pol,
                        P,
                        b_grid,
                        SD,
                        lmm_idx,
                        lmm_wt,
                        tmx_idx,
                        tmx_wt,
                        ust,
                        Pia,
                        Pi_z,
                    )
                    observed_birth_cohort = (
                        realize_current_choices_markov_income(
                            birth_cohort,
                            j + event_horizon,
                            loc_probs,
                            tenure_choice,
                            tenure_probs,
                            lmm_idx,
                            lmm_wt,
                            tmx_idx,
                            tmx_wt,
                            use_compiled_scatter=use_compiled_scatter,
                        )
                        if bool(getattr(P, "use_postdecision_current_distribution", True))
                        else birth_cohort
                    )
                    post_h = mean_housing_distribution_markov(
                        observed_birth_cohort, j + event_horizon, hR_pol, P
                    )
                    # No-birth control: the same pre-birth households (same
                    # wealth/tenure mass) propagated childless over the same
                    # horizon. Differencing post_h against this nets out the
                    # common lifecycle housing drift, leaving the birth effect.
                    control_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                    control_cohort[
                        :, :, :, zz, 0, readiness_settled_state(P)
                    ] = birth_mass
                    control_cohort = advance_cohort_horizon_markov_income(
                        control_cohort, j, event_horizon, loc_probs, tenure_choice,
                        tenure_probs, bp_pol, P, b_grid, SD, lmm_idx, lmm_wt,
                        tmx_idx, tmx_wt, ust, Pia, Pi_z,
                    )
                    observed_control_cohort = (
                        realize_current_choices_markov_income(
                            control_cohort,
                            j + event_horizon,
                            loc_probs,
                            tenure_choice,
                            tenure_probs,
                            lmm_idx,
                            lmm_wt,
                            tmx_idx,
                            tmx_wt,
                            use_compiled_scatter=use_compiled_scatter,
                        )
                        if bool(getattr(P, "use_postdecision_current_distribution", True))
                        else control_cohort
                    )
                    control_post_h = mean_housing_distribution_markov(
                        observed_control_cohort, j + event_horizon, hR_pol, P
                    )
                    birth_es3_pre_sum += birth_mass_total * pre_h
                    birth_es3_post_sum += birth_mass_total * post_h
                    birth_es3_control_post_sum += birth_mass_total * control_post_h
                    birth_es3_mass += birth_mass_total
        elif bool(getattr(P, "joint_nested_choice", False)):
            if joint_pre is not None:
                joint_pre[:, :, :, j] = g[:, :, :, j]
            post, effective, born, attempts, risk = joint_nested.factor_age(
                g[:, :, :, j], P._joint_choice, P, j
            )
            g[:, :, :, j] = post
            tenure_probs[:, :, :, j] = effective
            first_births_by_age[j] = born[0]
            second_births_by_age[j] = born[1]
            third_births_by_age[j] = born[2]
            second_attempts_by_age[j], second_at_risk_by_age[j] = attempts[1], risk[1]
            third_attempts_by_age[j], third_at_risk_by_age[j] = attempts[2], risk[2]
            total_births += float(born.sum())
            births_by_loc[0] += float(born.sum())
        elif (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end):
            for zz in range(Nz):
                gc = g[:, :, :, j, zz, 0, 0]
                pa = fert_probs[:, :, :, j, zz, :]
                pv = np.arange(npar).reshape(1, 1, 1, npar)
                pi_j = float(fec[j])
                if bool(getattr(P, "sequential_births", False)):
                    # Snapshot every upward at-risk pool BEFORE any birth flow
                    # lands, so no household can chain two births in a period.
                    if independent_child_maturation_active(P):
                        at_risk_up = {
                            (nn, cs): g[:, :, :, j, zz, nn, cs].copy()
                            for nn in range(1, npar - 1)
                            for cs in range(0, nn + 1)
                        }
                    else:
                        at_risk_up = {
                            (nn, 1): g[:, :, :, j, zz, nn, 1].copy()
                            for nn in range(1, npar - 1)
                        }
                    if readiness_gate_active(P):
                        gc = g[:, :, :, j, zz, 0, 1]
                    m1 = gc * pa[:, :, :, 1]
                    realized1 = pi_j * m1
                    if readiness_gate_active(P):
                        g[:, :, :, j, zz, 0, 1] = gc - realized1
                    else:
                        g[:, :, :, j, zz, 0, 0] = gc - realized1
                    g[:, :, :, j, zz, 1, 1] += realized1
                    if parent_age_maturation_active(P):
                        newborn_inflow[1, 1] += float(np.sum(realized1))
                    first_births_by_age[j] += float(np.sum(realized1))
                    total_births += float(np.sum(realized1))
                    for i in range(I):
                        births_by_loc[i] += float(np.sum(realized1[:, :, i]))
                    p2all = getattr(P, "_fert2_probs", None)
                    if p2all is not None:
                        for nn in range(1, npar - 1):
                            child_states = range(0, nn + 1) if independent_child_maturation_active(P) else (1,)
                            for cs in child_states:
                                at_risk = at_risk_up[(nn, cs)]
                                if independent_child_maturation_active(P):
                                    attempt_prob = p2all[:, :, :, j, zz, 1, nn - 1, cs]
                                else:
                                    attempt_prob = p2all[:, :, :, j, zz, 1, nn - 1]
                                m2 = at_risk * attempt_prob
                                realized2 = pi_j * m2
                                g[:, :, :, j, zz, nn, cs] -= realized2
                                destination_cs = birth_destination_child_state(P, cs)
                                g[:, :, :, j, zz, nn + 1, destination_cs] += realized2
                                if parent_age_maturation_active(P):
                                    newborn_inflow[nn + 1, destination_cs] += float(
                                        np.sum(realized2)
                                    )
                                if nn == 1:
                                    second_attempts_by_age[j] += float(np.sum(m2))
                                    second_births_by_age[j] += float(np.sum(realized2))
                                    second_at_risk_by_age[j] += float(np.sum(at_risk))
                                else:
                                    third_attempts_by_age[j] += float(np.sum(m2))
                                    third_births_by_age[j] += float(np.sum(realized2))
                                    third_at_risk_by_age[j] += float(np.sum(at_risk))
                                total_births += float(np.sum(realized2))
                                for i in range(I):
                                    births_by_loc[i] += float(np.sum(realized2[:, :, i]))
                    birth_mass = realized1
                    birth_mass_total = float(np.sum(birth_mass))
                    if not fast_stats and birth_mass_total > 1e-12 and (j + event_horizon) < J:
                        birth_weight = np.zeros((Nb, nt, I, Nz))
                        birth_weight[:, :, :, zz] = birth_mass
                        pre_h = mean_housing_childless_weighted_markov(birth_weight, j, hR_pol, P)
                        birth_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                        birth_cohort[:, :, :, zz, 1, 1] = realized1
                        birth_cohort = advance_cohort_horizon_markov_income(
                            birth_cohort,
                            j,
                            event_horizon,
                            loc_probs,
                            tenure_choice,
                            tenure_probs,
                            bp_pol,
                            P,
                            b_grid,
                            SD,
                            lmm_idx,
                            lmm_wt,
                            tmx_idx,
                            tmx_wt,
                            ust,
                            Pia,
                            Pi_z,
                        )
                        observed_birth_cohort = (
                            realize_current_choices_markov_income(
                                birth_cohort,
                                j + event_horizon,
                                loc_probs,
                                tenure_choice,
                                tenure_probs,
                                lmm_idx,
                                lmm_wt,
                                tmx_idx,
                                tmx_wt,
                                use_compiled_scatter=use_compiled_scatter,
                            )
                            if bool(getattr(P, "use_postdecision_current_distribution", True))
                            else birth_cohort
                        )
                        post_h = mean_housing_distribution_markov(
                            observed_birth_cohort, j + event_horizon, hR_pol, P
                        )
                        # No-birth control: the same pre-birth households (same
                        # wealth/tenure mass) propagated childless over the same
                        # horizon. Differencing post_h against this nets out the
                        # common lifecycle housing drift, leaving the birth effect.
                        control_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                        control_cohort[
                            :, :, :, zz, 0, readiness_settled_state(P)
                        ] = birth_mass
                        control_cohort = advance_cohort_horizon_markov_income(
                            control_cohort, j, event_horizon, loc_probs, tenure_choice,
                            tenure_probs, bp_pol, P, b_grid, SD, lmm_idx, lmm_wt,
                            tmx_idx, tmx_wt, ust, Pia, Pi_z,
                        )
                        observed_control_cohort = (
                            realize_current_choices_markov_income(
                                control_cohort,
                                j + event_horizon,
                                loc_probs,
                                tenure_choice,
                                tenure_probs,
                                lmm_idx,
                                lmm_wt,
                                tmx_idx,
                                tmx_wt,
                                use_compiled_scatter=use_compiled_scatter,
                            )
                            if bool(getattr(P, "use_postdecision_current_distribution", True))
                            else control_cohort
                        )
                        control_post_h = mean_housing_distribution_markov(
                            observed_control_cohort, j + event_horizon, hR_pol, P
                        )
                        birth_es3_pre_sum += birth_mass_total * pre_h
                        birth_es3_post_sum += birth_mass_total * post_h
                        birth_es3_control_post_sum += birth_mass_total * control_post_h
                        birth_es3_mass += birth_mass_total
                    continue
                mbp = gc[:, :, :, None] * pa
                gpf = np.zeros((Nb, nt, I, npar, ncs))
                if pi_j < 1.0:
                    ba_j = pi_j * gc * np.sum(pa * pv, axis=3)
                    gpf[:, :, :, 0, 0] = mbp[:, :, :, 0] + (1.0 - pi_j) * np.sum(mbp[:, :, :, 1:], axis=3)
                    gpf[:, :, :, 1:, 1] = pi_j * mbp[:, :, :, 1:]
                else:
                    ba_j = gc * np.sum(pa * pv, axis=3)
                    gpf[:, :, :, 0, 0] = mbp[:, :, :, 0]
                    gpf[:, :, :, 1:, 1] = mbp[:, :, :, 1:]
                total_births += float(np.sum(ba_j))
                for i in range(I):
                    births_by_loc[i] += float(np.sum(ba_j[:, :, i]))
                g[:, :, :, j, zz, 0, 0] = 0.0
                g[:, :, :, j, zz, :, :] += gpf

                birth_mass = np.sum(mbp[:, :, :, 1:], axis=3)
                birth_mass_total = float(np.sum(birth_mass))
                if not fast_stats and birth_mass_total > 1e-12 and (j + event_horizon) < J:
                    birth_weight = np.zeros((Nb, nt, I, Nz))
                    birth_weight[:, :, :, zz] = birth_mass
                    pre_h = mean_housing_childless_weighted_markov(birth_weight, j, hR_pol, P)
                    birth_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                    birth_cohort[:, :, :, zz, 1:, 1] = mbp[:, :, :, 1:]
                    birth_cohort = advance_cohort_horizon_markov_income(
                        birth_cohort,
                        j,
                        event_horizon,
                        loc_probs,
                        tenure_choice,
                        tenure_probs,
                        bp_pol,
                        P,
                        b_grid,
                        SD,
                        lmm_idx,
                        lmm_wt,
                        tmx_idx,
                        tmx_wt,
                        ust,
                        Pia,
                        Pi_z,
                    )
                    observed_birth_cohort = (
                        realize_current_choices_markov_income(
                            birth_cohort,
                            j + event_horizon,
                            loc_probs,
                            tenure_choice,
                            tenure_probs,
                            lmm_idx,
                            lmm_wt,
                            tmx_idx,
                            tmx_wt,
                            use_compiled_scatter=use_compiled_scatter,
                        )
                        if bool(getattr(P, "use_postdecision_current_distribution", True))
                        else birth_cohort
                    )
                    post_h = mean_housing_distribution_markov(
                        observed_birth_cohort, j + event_horizon, hR_pol, P
                    )
                    # No-birth control: the same pre-birth households (same
                    # wealth/tenure mass) propagated childless over the same
                    # horizon. Differencing post_h against this nets out the
                    # common lifecycle housing drift, leaving the birth effect.
                    control_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                    control_cohort[:, :, :, zz, 0, 0] = birth_mass
                    control_cohort = advance_cohort_horizon_markov_income(
                        control_cohort, j, event_horizon, loc_probs, tenure_choice,
                        tenure_probs, bp_pol, P, b_grid, SD, lmm_idx, lmm_wt,
                        tmx_idx, tmx_wt, ust, Pia, Pi_z,
                    )
                    observed_control_cohort = (
                        realize_current_choices_markov_income(
                            control_cohort,
                            j + event_horizon,
                            loc_probs,
                            tenure_choice,
                            tenure_probs,
                            lmm_idx,
                            lmm_wt,
                            tmx_idx,
                            tmx_wt,
                            use_compiled_scatter=use_compiled_scatter,
                        )
                        if bool(getattr(P, "use_postdecision_current_distribution", True))
                        else control_cohort
                    )
                    control_post_h = mean_housing_distribution_markov(
                        observed_control_cohort, j + event_horizon, hR_pol, P
                    )
                    birth_es3_pre_sum += birth_mass_total * pre_h
                    birth_es3_post_sum += birth_mass_total * post_h
                    birth_es3_control_post_sum += birth_mass_total * control_post_h
                    birth_es3_mass += birth_mass_total

                    one_child = np.zeros_like(observed_birth_cohort)
                    one_child[:, :, :, :, 1, :] = observed_birth_cohort[:, :, :, :, 1, :]
                    one_mass = float(np.sum(one_child))
                    if one_mass > 1e-12:
                        one_child_birth_mass = mbp[:, :, :, 1]
                        one_child_birth_mass_total = float(np.sum(one_child_birth_mass))
                        addchild_es3_one_sum += one_mass * mean_housing_distribution_markov(
                            one_child, j + event_horizon, hR_pol, P
                        )
                        addchild_es3_one_mass += one_mass
                        if one_child_birth_mass_total > 1e-12:
                            one_weight = np.zeros((Nb, nt, I, Nz))
                            one_weight[:, :, :, zz] = one_child_birth_mass
                            one_pre = mean_housing_childless_weighted_markov(one_weight, j, hR_pol, P)
                            onechild_es3_pre_sum += one_child_birth_mass_total * one_pre
                            onechild_es3_post_sum += one_mass * mean_housing_distribution_markov(
                                one_child, j + event_horizon, hR_pol, P
                            )
                            onechild_es3_mass += one_child_birth_mass_total

                    if npar >= 3:
                        two_plus_birth_mass = np.sum(mbp[:, :, :, 2:], axis=3)
                        two_plus_birth_mass_total = float(np.sum(two_plus_birth_mass))
                        two_plus = np.zeros_like(observed_birth_cohort)
                        two_plus[:, :, :, :, 2:, :] = observed_birth_cohort[:, :, :, :, 2:, :]
                        two_plus_mass = float(np.sum(two_plus))
                        if two_plus_mass > 1e-12:
                            addchild_es3_two_plus_sum += two_plus_mass * mean_housing_distribution_markov(
                                two_plus, j + event_horizon, hR_pol, P
                            )
                            addchild_es3_two_plus_mass += two_plus_mass
                            if two_plus_birth_mass_total > 1e-12:
                                two_weight = np.zeros((Nb, nt, I, Nz))
                                two_weight[:, :, :, zz] = two_plus_birth_mass
                                two_pre = mean_housing_childless_weighted_markov(two_weight, j, hR_pol, P)
                                twoplus_es3_pre_sum += two_plus_birth_mass_total * two_pre
                                twoplus_es3_post_sum += two_plus_mass * mean_housing_distribution_markov(
                                    two_plus, j + event_horizon, hR_pol, P
                                )
                                twoplus_es3_mass += two_plus_birth_mass_total

        survival = float(P.survival_probs[j]) if bool(getattr(P, "use_age_survival", False)) else 1.0
        gj = survival * g[:, :, :, j, :, :, :]
        if count_active:
            count_post[:, :, :, j] = g[:, :, :, j]
        gpl = np.zeros((Nb, nt, I, Nz, npar, ncs))
        for zz in range(Nz):
            for io in range(I):
                for to in range(nt):
                    go = flat_nc(gj[:, to, io, zz, :, :], Nb, nc)
                    if np.sum(go) < 1e-15:
                        continue
                    po = np.reshape(loc_probs[:, to, io, :, j, zz, :, :], (Nb, I, nc), order="F")
                    sp = po[:, io, :]
                    gpl[:, to, io, zz, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
                    idx = lmm_idx[io, to, :]
                    wt = lmm_wt[io, to, :]
                    for id_ in range(I):
                        if id_ == io:
                            continue
                        mp = go * po[:, id_, :]
                        if use_compiled_scatter:
                            moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                        else:
                            moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb)
                        gpl[:, 0, id_, zz, :, :] += unflat_nc(moved, Nb, npar, ncs)

        gpt = np.zeros((Nb, nt, I, Nz, npar, ncs))
        gpt_stay = np.zeros((Nb, nt, I, Nz, npar, ncs)) if bp_pol_stay is not None else None
        for zz in range(Nz):
            for nn in range(npar):
                for id_ in range(I):
                    for to in range(nt):
                        gs = gpl[:, to, id_, zz, nn, :]
                        if np.sum(gs) < 1e-15:
                            continue
                        if tenure_probs is not None:
                            all_probs = np.asarray(
                                tenure_probs[:, to, id_, j, zz, nn, :, :], dtype=float
                            )
                            prob_sum = np.sum(all_probs, axis=-1)
                            normalized_probs = np.divide(
                                all_probs,
                                prob_sum[:, :, None],
                                out=np.zeros_like(all_probs),
                                where=prob_sum[:, :, None] > 0,
                            )
                        for tn in range(nt):
                            if tenure_probs is None:
                                tcs = tenure_choice[:, to, id_, j, zz, nn, :]
                                mk = tcs == tn
                                if not np.any(mk):
                                    continue
                                mt = gs * mk
                            else:
                                pr = normalized_probs[:, :, tn]
                                mt = gs * pr
                            if np.sum(mt) < 1e-15:
                                continue
                            rd = np.zeros((Nb, ncs))
                            for cs in range(ncs):
                                idx = tmx_idx[id_, to, tn, nn, cs, :]
                                wt = tmx_wt[id_, to, tn, nn, cs, :]
                                if use_compiled_scatter:
                                    rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                                else:
                                    rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                            if bp_pol_stay is not None and to == tn:
                                gpt_stay[:, tn, id_, zz, nn, :] += rd
                            else:
                                gpt[:, tn, id_, zz, nn, :] += rd

        gps = np.zeros((Nb, nt, I, Nz, npar, ncs))
        for zz in range(Nz):
            for i in range(I):
                for ten in range(nt):
                    gf = flat_nc(gpt[:, ten, i, zz, :, :], Nb, nc)
                    bpv = flat_nc(bp_pol[:, ten, i, j, zz, :, :], Nb, nc)
                    bpc = np.clip(bpv, bmin, bmax)
                    idx, wt = interp_indices(b_grid, bpc)
                    if use_compiled_scatter:
                        g_new = scatter_cols_kernel(idx, wt, gf, Nb)
                    else:
                        g_new = scatter_redistribute_cols(idx, wt, gf, Nb)
                    gps[:, ten, i, zz, :, :] = unflat_nc(g_new, Nb, npar, ncs)
                    if gpt_stay is not None:
                        assert bp_pol_stay is not None
                        gf_s = flat_nc(gpt_stay[:, ten, i, zz, :, :], Nb, nc)
                        bpv_s = flat_nc(bp_pol_stay[:, ten, i, j, zz, :, :], Nb, nc)
                        idx_s, wt_s = interp_indices(b_grid, np.clip(bpv_s, bmin, bmax))
                        if use_compiled_scatter:
                            g_new_s = scatter_cols_kernel(idx_s, wt_s, gf_s, Nb)
                        else:
                            g_new_s = scatter_redistribute_cols(idx_s, wt_s, gf_s, Nb)
                        gps[:, ten, i, zz, :, :] += unflat_nc(g_new_s, Nb, npar, ncs)

        for zz in range(Nz):
            for nn in range(npar):
                for cs in range(ncs):
                    gp = gps[:, :, :, zz, nn, cs]
                    if readiness_gate_active(P) and nn == 0 and cs in (0, 1):
                        current_age = float(P.age_start) + float(j) * float(P.da)
                        next_age = current_age + float(P.da)
                        hazard = readiness_transition_hazard(P, current_age, next_age)
                        readiness_weights = (
                            ((0, 1.0 - hazard), (1, hazard))
                            if cs == 0
                            else ((1, 1.0),)
                        )
                        for csn, readiness_weight in readiness_weights:
                            if readiness_weight <= 0.0:
                                continue
                            for zn in range(Nz):
                                transition_weight = Pi_z[zz, zn]
                                if transition_weight > 0.0:
                                    g[:, :, :, j + 1, zn, nn, csn] += (
                                        transition_weight * readiness_weight * gp
                                    )
                        continue
                    if ust:
                        Pi = _child_Pa_for_age(P, Pia, j)[:, :, nn]
                        if parent_age_maturation_active(P):
                            # m-d exemption: blend the exempt row by the newborn
                            # share of this cell (exact for cell totals/entrants).
                            tot_cell = float(np.sum(gps[:, :, :, :, nn, cs]))
                            if tot_cell > 0.0:
                                surv_j = (
                                    float(P.survival_probs[j])
                                    if bool(getattr(P, "use_age_survival", False))
                                    else 1.0
                                )
                                f_cell = float(
                                    np.clip(
                                        surv_j * newborn_inflow[nn, cs] / tot_cell,
                                        0.0,
                                        1.0,
                                    )
                                )
                                if f_cell > 0.0:
                                    Pi = _blended_child_Pi_for_cell(P, Pia, j, f_cell)[
                                        :, :, nn
                                    ]
                        if not independent_child_maturation_active(P) and cs == K and nn >= 1:
                            pm = Pi[cs, csm1] if nn == 1 else Pi[cs, csm2]
                            if pm > 0:
                                nk = nn
                                for im in range(I):
                                    fi = ecf * nk * pm * float(np.sum(gp[:, :, im]))
                                    entrants_mature_by_loc[im] += fi
                                    entrants_mature_total += fi
                        for csn in range(ncs):
                            wt_child = Pi[cs, csn]
                            if wt_child > 0:
                                if independent_child_maturation_active(P) and csn < cs:
                                    matured = cs - csn
                                    for im in range(I):
                                        fi = ecf * matured * wt_child * float(np.sum(gp[:, :, im]))
                                        entrants_mature_by_loc[im] += fi
                                        entrants_mature_total += fi
                                for zn in range(Nz):
                                    transition_weight = Pi_z[zz, zn]
                                    if transition_weight > 0.0:
                                        g[:, :, :, j + 1, zn, nn, csn] += (
                                            transition_weight * wt_child * gp
                                        )
                    else:
                        if cs == 0:
                            csn = 0
                        elif cs >= csm1:
                            csn = cs
                        elif cs < K:
                            csn = cs + 1
                        else:
                            csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                        if cs == K and csn >= csm1 and nn >= 1:
                            nk = nn
                            for im in range(I):
                                fi = ecf * nk * float(np.sum(gp[:, :, im]))
                                entrants_mature_by_loc[im] += fi
                                entrants_mature_total += fi
                        for zn in range(Nz):
                            transition_weight = Pi_z[zz, zn]
                            if transition_weight > 0.0:
                                g[:, :, :, j + 1, zn, nn, csn] += (
                                    transition_weight * gp
                                )

    _gate_dead_mass_at_age(
        g[:, :, :, J - 1, :, :, :],
        state_values[:, :, :, J - 1, :, :, :],
        J - 1,
        f"forward_age_{P.age_start + (J - 1) * P.da:g}",
        r_hat,
        p_hat,
        P,
        b_grid,
        SD,
        markov_income=True,
    )

    if bool(getattr(P, "joint_nested_choice", False)):
        if joint_pre is not None:
            joint_pre[:, :, :, J - 1] = g[:, :, :, J - 1]
        post, effective, born, attempts, risk = joint_nested.factor_age(
            g[:, :, :, J - 1], P._joint_choice, P, J - 1
        )
        if born.sum() != 0:
            raise NotImplementedError("Fertility in final age is outside this experiment")
        g[:, :, :, J - 1] = post
        tenure_probs[:, :, :, J - 1] = effective

    tm = float(np.sum(g))
    if tm > 1e-12 and normalize_population_mass(P):
        sc = P.N_target / tm
        g *= sc
        if count_active:
            count_pre *= sc
            count_post *= sc
            count_first_tagged *= sc
            count_risk *= sc
            count_attempts *= sc
            count_any_births *= sc
        if joint_pre is not None:
            joint_pre *= sc
        total_births *= sc
        births_by_loc *= sc
        first_births_by_age *= sc
        second_births_by_age *= sc
        second_attempts_by_age *= sc
        second_at_risk_by_age *= sc
        third_births_by_age *= sc
        third_attempts_by_age *= sc
        third_at_risk_by_age *= sc
        entrants_mature_by_loc *= sc
        entrants_mature_total *= sc

    if count_active:
        count_pre[:, :, :, -1] = g[:, :, :, -1]
        count_post[:, :, :, -1] = g[:, :, :, -1]
        P.birth_count_pre_distribution = count_pre
        P.birth_count_post_distribution = count_post
        P.birth_count_first_birth_tagged_distribution = count_first_tagged
    use_postdecision_current = bool(getattr(P, "use_postdecision_current_distribution", True))
    g_current = (
        realize_current_cross_section(
            g,
            loc_probs,
            tenure_choice,
            tenure_probs,
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            use_compiled_scatter=use_compiled_scatter,
        )
        if use_postdecision_current
        else g
    )
    if normalize_population_mass(P):
        assert np.isclose(np.sum(g_current), np.sum(g), rtol=0.0, atol=1e-10)
    if bool(getattr(P, "native_due_stayer_credit", False)):
        P._g_stay_distribution = realize_stayer_cross_section(g, loc_probs, tenure_choice, tenure_probs)
    if fast_stats:
        stats = compute_markov_eq_stats(g_current, P, b_grid, p_hat, hR_pol)
    else:
        stats = compute_markov_statistics(
            g_current,
            fert_probs,
            loc_probs,
            P,
            b_grid,
            p_hat,
            hR_pol,
            asset_g=g,
            bequest_g=g_current,
            bp_pol=bp_pol,
        )
    if bool(getattr(P, "native_due_stayer_credit", False)):
        stats.g_stay_distribution = P._g_stay_distribution.copy()
    if count_active:
        count_pre[:, :, :, -1] = g[:, :, :, -1]
        count_post[:, :, :, -1] = g[:, :, :, -1]
        P.birth_count_pre_distribution = count_pre
        P.birth_count_post_distribution = count_post
        P.birth_count_first_birth_tagged_distribution = count_first_tagged
        stats.birth_count_pre_distribution = count_pre.copy()
        stats.birth_count_post_distribution = count_post.copy()
        stats.birth_count_first_birth_tagged_distribution = count_first_tagged.copy()
        stats.birth_count_births_by_order_by_age = np.stack((first_births_by_age, second_births_by_age, third_births_by_age))
        stats.birth_count_at_risk_by_order_by_age = count_risk.copy()
        stats.birth_count_attempts_by_order_by_age = count_attempts.copy()
        stats.birth_count_any_birth_households_by_age = count_any_births.copy()
        stats.birth_count_expected_children_by_age = np.sum(stats.birth_count_births_by_order_by_age, axis=0)
        stats.birth_count_hazards_by_order_by_age = np.divide(stats.birth_count_births_by_order_by_age, count_risk, out=np.zeros_like(count_risk), where=count_risk > 0)
        stats.birth_count_order_risk_definition = "n<q and min(cap,3-n)>=q-n; crossings of several orders count once per order"
        stats.birth_count_attempt_definition = "P(intended births k >= q-n), conditional on the same order-risk pool"
    P._second_births_by_age = second_births_by_age
    P._second_attempts_by_age = second_attempts_by_age
    P._first_births_by_age = first_births_by_age
    P._second_at_risk_by_age = second_at_risk_by_age
    P._third_births_by_age = third_births_by_age
    P._third_attempts_by_age = third_attempts_by_age
    P._third_at_risk_by_age = third_at_risk_by_age
    if bool(getattr(P, "sequential_births", False)):
        stats.second_attempt_hazard_by_age = second_attempts_by_age / np.maximum(second_at_risk_by_age, 1e-12)
        stats.second_birth_hazard_by_age = second_births_by_age / np.maximum(second_at_risk_by_age, 1e-12)
        stats.parity_progression_1to2_flow = float(np.sum(second_births_by_age) / max(np.sum(first_births_by_age), 1e-12))
        stats.third_attempt_hazard_by_age = third_attempts_by_age / np.maximum(third_at_risk_by_age, 1e-12)
        stats.third_birth_hazard_by_age = third_births_by_age / np.maximum(third_at_risk_by_age, 1e-12)
        stats.parity_progression_2to3_flow = float(np.sum(third_births_by_age) / max(np.sum(second_births_by_age), 1e-12))
    if bool(getattr(P, "joint_nested_choice", False)) and np.any(SD.birth_entry_grant):
        raise NotImplementedError("Joint experimental grant accounting is not implemented")
    grant_recipient_mass, grant_outlays = markov_grant_outlays(
        g,
        tenure_choice,
        tenure_probs,
        P,
        SD,
    )
    property_tax_revenue = property_tax_revenue_from_distribution(
        g_current,
        hR_pol,
        p_hat,
        P,
    )
    transfer_outlays = float(getattr(P, "property_tax_lump_sum_transfer", 0.0)) * float(np.sum(g_current))
    stats.property_tax_revenue = property_tax_revenue
    stats.birth_entry_grant_recipient_mass = grant_recipient_mass
    stats.birth_entry_grant_outlays = grant_outlays
    stats.property_tax_transfer_outlays = transfer_outlays
    stats.property_tax_budget_residual = property_tax_revenue - grant_outlays - transfer_outlays
    stats.entry_censored_mass = float(getattr(P, "_entry_censored_mass", 0.0))
    stats.entry_censored_share = stats.entry_censored_mass / max(
        float(getattr(P, "_entry_total_mass", 0.0)), 1e-300
    )
    stats.total_births_kfe = total_births
    if not fast_stats:
        stats.g_beginning_distribution = g.copy()
        stats.g_cross_sectional_wealth_distribution = g.copy()
    stats.births_by_loc = births_by_loc
    stats.entry_by_loc = np.sum(g[:, :, :, 0, :, :, :], axis=(0, 1, 3, 4, 5))
    stats.entry_rate = float(np.sum(g[:, :, :, 0, :, :, :]))
    stats.total_mass = float(np.sum(g_current))
    stats.entrants_mature_by_loc = entrants_mature_by_loc
    stats.entrants_mature_total = entrants_mature_total
    stats.mature_entry_shares = entrants_mature_by_loc / max(entrants_mature_total, 1e-12)
    if str(getattr(P, "adult_entry_clock", "child_departure")) == "split_birth_vintage":
        if not (bool(getattr(P, "sequential_births", False)) or bool(getattr(P, "joint_nested_choice", False))):
            raise ValueError("Split birth-vintage entry requires observed sequential third-birth flow")
        if float(P.period_years) != 4.0:
            raise ValueError("Split birth-vintage entry requires four-year model periods")
        if int(P.I) != 1 or int(P.n_parity) != 4 or str(getattr(P, "fertility_units", "")) != "literal_topcode":
            raise ValueError("Split birth-vintage stationary entry requires one market and literal 3+ fertility units")
        top_weight = float(getattr(P, "tfr_top_bin_weight"))
        birth_children = adjusted_births(total_births, float(np.sum(third_births_by_age)), top_weight)
        stats.adult_entry_adjusted_birth_children = birth_children
        stats.adult_entry_potential_total = potential_entry_households(birth_children)
        stats.adult_entry_potential_by_loc = np.array([stats.adult_entry_potential_total])
    attach_entry_wealth_stats(stats, b_grid, entry_idx, entry_wt, P)
    stats.housing_increment_0to1_eventstudy_t3 = (
        (birth_es3_post_sum - birth_es3_control_post_sum) / birth_es3_mass if birth_es3_mass > 1e-12 else 0.0
    )
    stats.housing_increment_1to2_proxy_t3 = (
        addchild_es3_two_plus_sum / addchild_es3_two_plus_mass - addchild_es3_one_sum / addchild_es3_one_mass
        if addchild_es3_one_mass > 1e-12 and addchild_es3_two_plus_mass > 1e-12
        else 0.0
    )
    stats.housing_increment_0to1_onechild_eventstudy_t3 = (
        onechild_es3_post_sum / onechild_es3_mass - onechild_es3_pre_sum / onechild_es3_mass
        if onechild_es3_mass > 1e-12
        else 0.0
    )
    stats.housing_increment_0to2plus_eventstudy_t3 = (
        twoplus_es3_post_sum / twoplus_es3_mass - twoplus_es3_pre_sum / twoplus_es3_mass
        if twoplus_es3_mass > 1e-12
        else 0.0
    )
    stats.housing_event_horizon = event_horizon
    stats.current_distribution_timing = (
        "post_housing_choice"
        if bool(getattr(P, "use_postdecision_current_distribution", True))
        else "beginning_of_period_legacy"
    )
    stats.wealth_moment_timing = "beginning_of_period_state"
    if bool(getattr(P, "joint_nested_choice", False)):
        for name in ("housing_increment_0to1_eventstudy_t3", "housing_increment_1to2_proxy_t3",
                     "housing_increment_0to1_onechild_eventstudy_t3", "housing_increment_0to2plus_eventstudy_t3"):
            setattr(stats, name, float("nan"))
        stats.stationary_eventstudy_status = "not_measured_fast_statistics"
        if not fast_stats:
            try:
                stats.housing_increment_0to1_eventstudy_t3 = joint_nested.stationary_first_birth_response(
                    joint_pre, P._joint_choice, P, b_grid, SD, loc_probs, tenure_choice,
                    bp_pol, hR_pol, (lmm_idx, lmm_wt, tmx_idx, tmx_wt))
            except joint_nested.UndefinedStationaryFirstBirthSupport as error:
                # Intermediate fertility-normalization trials may have no births.
                # Keep this conditional moment unavailable; final target and
                # stationary-nesting validation still require every row finite.
                stats.stationary_eventstudy_status = "undefined_first_birth_support"
                stats.stationary_eventstudy_branch_masses = error.masses
            else:
                stats.stationary_eventstudy_status = "joint_matched_one_period_branch"
    return g_current, stats


def collapse_markov_policy(policy: np.ndarray, g: np.ndarray, z_weights: np.ndarray) -> np.ndarray:
    fallback = np.tensordot(policy, z_weights, axes=([4], [0]))
    den = np.sum(g, axis=4)
    num = np.sum(policy * g, axis=4)
    out = fallback.copy()
    mask = den > 1e-15
    out[mask] = num[mask] / den[mask]
    return out


def collapse_markov_fertility_probs(
    fp: np.ndarray,
    g: np.ndarray,
    z_weights: np.ndarray,
    P: SimpleNamespace | None = None,
) -> np.ndarray:
    fallback = np.tensordot(fp, z_weights, axes=([4], [0]))
    settled_cs = readiness_settled_state(P) if P is not None else 0
    mass = g[:, :, :, :, :, 0, settled_cs]
    den = np.sum(mass, axis=4)
    num = np.sum(fp * mass[:, :, :, :, :, None], axis=4)
    out = fallback.copy()
    mask = den > 1e-15
    out[mask, :] = num[mask, :] / den[mask, None]
    return out


def collapse_markov_location_probs(lp: np.ndarray, g: np.ndarray, z_weights: np.ndarray) -> np.ndarray:
    fallback = np.tensordot(lp, z_weights, axes=([5], [0]))
    mass = g[:, :, :, None, :, :, :, :]
    den = np.sum(mass, axis=5)
    num = np.sum(lp * mass, axis=5)
    return np.divide(num, den, out=fallback.copy(), where=den > 1e-15)


def advance_cohort_horizon_markov_income(
    g_in,
    start_age,
    horizon,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
    Pi_z,
    bp_pol_stay=None,
):
    g_out = g_in
    for step in range(1, horizon + 1):
        age_idx = start_age + step - 1
        if age_idx >= P.J - 1:
            break
        g_out = advance_cohort_one_period_markov_income(
            g_out,
            age_idx,
            loc_probs,
            tenure_choice,
            tenure_probs,
            bp_pol,
            P,
            b_grid,
            SD,
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            ust,
            Pia,
            Pi_z,
            bp_pol_stay=bp_pol_stay,
        )
    return g_out


def advance_cohort_one_period_markov_income(
    gj,
    j,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
    Pi_z,
    *,
    mass_pruning_tolerance: float = 1e-15,
    newborn_frac=None,
    bp_pol_stay=None,
):
    """Advance one cohort one period (Markov income).

    ``newborn_frac`` is an optional ``(n_parity, n_child_states)`` array of
    post-birth newborn shares per cell for the parent-age m-d exemption
    (blended standard/exempt rows; ``None`` reproduces the constant path
    bit for bit).
    ``bp_pol_stay`` is an optional stayer savings policy with the same shape
    as ``bp_pol``; when given, mass that stays in its tenure (to == tn) is
    scattered with the stayer policy and all other mass with ``bp_pol``.
    ``None`` reproduces the legacy single-policy scatter bit for bit.
    """
    Nb = len(b_grid)
    nt = 1 + P.n_house
    I = P.I
    Nz = gj.shape[3]
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2

    gpl = np.zeros((Nb, nt, I, Nz, npar, ncs))
    for zz in range(Nz):
        for io in range(I):
            for to in range(nt):
                go = flat_nc(gj[:, to, io, zz, :, :], Nb, nc)
                if np.sum(go) == 0.0 or np.sum(go) < mass_pruning_tolerance:
                    continue
                po = np.reshape(loc_probs[:, to, io, :, j, zz, :, :], (Nb, I, nc), order="F")
                sp = po[:, io, :]
                gpl[:, to, io, zz, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
                idx = lmm_idx[io, to, :]
                wt = lmm_wt[io, to, :]
                for id_ in range(I):
                    if id_ == io:
                        continue
                    mp = go * po[:, id_, :]
                    if use_compiled_scatter:
                        moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                    else:
                        moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb, mass_pruning_tolerance=mass_pruning_tolerance)
                    gpl[:, 0, id_, zz, :, :] += unflat_nc(moved, Nb, npar, ncs)

    gpt = np.zeros((Nb, nt, I, Nz, npar, ncs))
    gpt_stay = np.zeros((Nb, nt, I, Nz, npar, ncs)) if bp_pol_stay is not None else None
    for zz in range(Nz):
        for nn in range(npar):
            for id_ in range(I):
                for to in range(nt):
                    gs = gpl[:, to, id_, zz, nn, :]
                    if np.sum(gs) == 0.0 or np.sum(gs) < mass_pruning_tolerance:
                        continue
                    normalized_probs = None
                    if tenure_probs is not None:
                        all_probs = np.asarray(
                            tenure_probs[:, to, id_, j, zz, nn, :, :], dtype=float
                        )
                        prob_sum = np.sum(all_probs, axis=-1)
                        normalized_probs = np.divide(
                            all_probs,
                            prob_sum[:, :, None],
                            out=np.zeros_like(all_probs),
                            where=prob_sum[:, :, None] > 0.0,
                        )
                    for tn in range(nt):
                        if tenure_probs is None:
                            tcs = tenure_choice[:, to, id_, j, zz, nn, :]
                            mk = tcs == tn
                            if not np.any(mk):
                                continue
                            mt = gs * mk
                        else:
                            pr = normalized_probs[:, :, tn]
                            mt = gs * pr
                        if np.sum(mt) == 0.0 or np.sum(mt) < mass_pruning_tolerance:
                            continue
                        rd = np.zeros((Nb, ncs))
                        for cs in range(ncs):
                            idx = tmx_idx[id_, to, tn, nn, cs, :]
                            wt = tmx_wt[id_, to, tn, nn, cs, :]
                            if use_compiled_scatter:
                                rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                            else:
                                rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                        if gpt_stay is not None and to == tn:
                            gpt_stay[:, tn, id_, zz, nn, :] += rd
                        else:
                            gpt[:, tn, id_, zz, nn, :] += rd

    gps = np.zeros((Nb, nt, I, Nz, npar, ncs))
    for zz in range(Nz):
        for i in range(I):
            for ten in range(nt):
                gf = flat_nc(gpt[:, ten, i, zz, :, :], Nb, nc)
                bpv = flat_nc(bp_pol[:, ten, i, j, zz, :, :], Nb, nc)
                idx, wt = interp_indices(b_grid, np.clip(bpv, b_grid[0], b_grid[-1]))
                if use_compiled_scatter:
                    g_new = scatter_cols_kernel(idx, wt, gf, Nb)
                else:
                    g_new = scatter_redistribute_cols(idx, wt, gf, Nb, mass_pruning_tolerance=mass_pruning_tolerance)
                gps[:, ten, i, zz, :, :] = unflat_nc(g_new, Nb, npar, ncs)
                if gpt_stay is not None:
                    gf_stay = flat_nc(gpt_stay[:, ten, i, zz, :, :], Nb, nc)
                    bpv_stay = flat_nc(bp_pol_stay[:, ten, i, j, zz, :, :], Nb, nc)
                    idx_s, wt_s = interp_indices(b_grid, np.clip(bpv_stay, b_grid[0], b_grid[-1]))
                    if use_compiled_scatter:
                        g_new_stay = scatter_cols_kernel(idx_s, wt_s, gf_stay, Nb)
                    else:
                        g_new_stay = scatter_redistribute_cols(idx_s, wt_s, gf_stay, Nb, mass_pruning_tolerance=mass_pruning_tolerance)
                    gps[:, ten, i, zz, :, :] += unflat_nc(g_new_stay, Nb, npar, ncs)

    g_next = np.zeros_like(gj)
    for zz in range(Nz):
        for nn in range(npar):
            for cs in range(ncs):
                gp = gps[:, :, :, zz, nn, cs]
                if readiness_gate_active(P) and nn == 0 and cs in (0, 1):
                    current_age = float(P.age_start) + float(j) * float(P.da)
                    next_age = current_age + float(P.da)
                    hazard = readiness_transition_hazard(P, current_age, next_age)
                    readiness_weights = (
                        ((0, 1.0 - hazard), (1, hazard))
                        if cs == 0
                        else ((1, 1.0),)
                    )
                    for csn, readiness_weight in readiness_weights:
                        if readiness_weight <= 0.0:
                            continue
                        for zn in range(Nz):
                            transition_weight = Pi_z[zz, zn]
                            if transition_weight > 0.0:
                                g_next[:, :, :, zn, nn, csn] += (
                                    transition_weight * readiness_weight * gp
                                )
                    continue
                if ust:
                    Pi = _child_Pa_for_age(P, Pia, j)[:, :, nn]
                    if parent_age_maturation_active(P) and newborn_frac is not None:
                        f_cell = float(
                            np.clip(
                                float(np.asarray(newborn_frac, dtype=float)[nn, cs]),
                                0.0,
                                1.0,
                            )
                        )
                        Pi_ex_cell = _child_Pa_exempt_for_age(P, j)
                        if f_cell > 0.0 and Pi_ex_cell is not None:
                            Pi = (1.0 - f_cell) * Pi + f_cell * np.asarray(
                                Pi_ex_cell
                            )[:, :, nn]
                    for csn in range(ncs):
                        wt_child = Pi[cs, csn]
                        if wt_child > 0:
                            for zn in range(Nz):
                                transition_weight = Pi_z[zz, zn]
                                if transition_weight > 0.0:
                                    g_next[:, :, :, zn, nn, csn] += (
                                        transition_weight * wt_child * gp
                                    )
                else:
                    if cs == 0:
                        csn = 0
                    elif cs >= csm1:
                        csn = cs
                    elif cs < K:
                        csn = cs + 1
                    else:
                        csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                    for zn in range(Nz):
                        transition_weight = Pi_z[zz, zn]
                        if transition_weight > 0.0:
                            g_next[:, :, :, zn, nn, csn] += transition_weight * gp
    return g_next


def mean_housing_childless_weighted_markov(weight_dist, j, hR_pol, P):
    Nb, nt, I, Nz = weight_dist.shape
    th = mn = 0.0
    childless_cs = readiness_settled_state(P)
    for zz in range(Nz):
        for i in range(I):
            for ten in range(nt):
                gs = weight_dist[:, ten, i, zz]
                mh = float(np.sum(gs))
                if mh < 1e-15:
                    continue
                if ten == 0:
                    th += float(
                        np.sum(gs * hR_pol[:, ten, i, j, zz, 0, childless_cs])
                    )
                else:
                    th += mh * P.H_own[ten - 1]
                mn += mh
    return th / max(mn, 1e-12)


def mean_housing_distribution_markov(g_dist, j, hR_pol, P):
    Nb, nt, I, Nz, npar, ncs = g_dist.shape
    th = mn = 0.0
    for zz in range(Nz):
        for nn in range(npar):
            for cs in range(ncs):
                for i in range(I):
                    for ten in range(nt):
                        gs = g_dist[:, ten, i, zz, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        if ten == 0:
                            th += float(np.sum(gs * hR_pol[:, ten, i, j, zz, nn, cs]))
                        else:
                            th += mh * P.H_own[ten - 1]
                        mn += mh
    return th / max(mn, 1e-12)


def markov_renter_room_moments(g: np.ndarray, hR: np.ndarray, P: SimpleNamespace) -> dict:
    """Renter room threshold/median moments on the FULL income-resolved state.

    `compute_statistics` operates on the income-collapsed renter policy
    `hR_total = E_z[h_R]`. That is correct for linear moments (means, mass
    shares over discrete owner rungs) but WRONG for nonlinear operators on the
    continuous renter policy, because for a policy that varies across income
    states within a (b, j, n, cs) cell,
        E_z[1{h_R >= t}] != 1{E_z[h_R] >= t}  and
        median_z(h_R)    != median(E_z[h_R]).
    These threshold/median moments must be evaluated state-by-income, then
    mass-aggregated. Means are unaffected, so the cross-tenure mean room gap
    (and any mean-based target) is identical under either path.
    """
    Nb, nt, I, J, Nz, npar, ncs = g.shape
    dep_last = P.n_child_stages
    hcut = int(getattr(P, "child_bin_high_cutoff", 2))
    a25s, a45e = age_to_index(P, 25), age_to_index(P, 45)
    a30s, a55e = age_to_index(P, 30), age_to_index(P, 55)
    hRmax = float(P.hR_max)
    ge6_num = ge6_den = 0.0
    cap_all_num = cap_all_den = 0.0
    cap_c0_num = cap_c0_den = 0.0
    cap_c1_num = cap_c1_den = 0.0
    med_vals: list[np.ndarray] = []
    med_wts: list[np.ndarray] = []
    for j in range(J):
        in_3055 = a30s <= j <= a55e
        in_2545 = a25s <= j <= a45e
        if not (in_3055 or in_2545):
            continue
        for i in range(I):
            for zz in range(Nz):
                for nn in range(npar):
                    for cs in range(ncs):
                        cb = current_child_bin_dt(nn, cs, dep_last, hcut, getattr(P, "child_state_mode", "shared_clock"))
                        gr = g[:, 0, i, j, zz, nn, cs]
                        hr = hR[:, 0, i, j, zz, nn, cs]
                        kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                        if not np.any(kr):
                            continue
                        wr = gr[kr]
                        rr = hr[kr]
                        m = float(np.sum(wr))
                        if in_3055 and cb == 2:
                            ge6_den += m
                            ge6_num += float(np.sum(wr[rr >= 6.0 - 1e-8]))
                        if in_2545:
                            cap = float(np.sum(wr[rr >= hRmax - 1e-8]))
                            cap_all_den += m
                            cap_all_num += cap
                            if cb == 2:
                                cap_c0_den += m
                                cap_c0_num += cap
                                med_vals.append(rr)
                                med_wts.append(wr)
                            elif cb == 3:
                                cap_c1_den += m
                                cap_c1_num += cap
    return {
        "prime30_55_childless_renter_share_rooms_ge6": ge6_num / max(ge6_den, 1e-12),
        "prime_childless_renter_median_rooms": weighted_median_from_cells(med_vals, med_wts),
        "renter25_45_all_cap_share": cap_all_num / max(cap_all_den, 1e-12),
        "renter25_45_current0_cap_share": cap_c0_num / max(cap_c0_den, 1e-12),
        "renter25_45_current1_cap_share": cap_c1_num / max(cap_c1_den, 1e-12),
    }


def add_annual_gross_liquid_wealth_moments(stats: SimpleNamespace, g: np.ndarray, P: SimpleNamespace, bg: np.ndarray) -> None:
    """Add wealth/income moments with annual gross-income denominators.

    Core model accounting uses 4-year period after-tax income. The empirical
    PSID entry-wealth summaries are annual-income ratios, so these moments make
    the denominator explicit and avoid comparing annual data to period-income
    model statistics.
    """
    g_arr = np.asarray(g, dtype=float)
    bg_arr = np.asarray(bg, dtype=float).reshape(-1)
    if g_arr.ndim == 6:
        g7 = g_arr[:, :, :, :, None, :, :]
        z_values = np.array([1.0])
    elif g_arr.ndim == 7:
        g7 = g_arr
        z_values = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    else:
        return

    dep_last = int(getattr(P, "n_child_stages", 1))
    hcut = int(getattr(P, "child_bin_high_cutoff", 2))

    def sample_stats(age_lo: float, age_hi: float, *, childless_only: bool, renter_only: bool) -> tuple[float, float, float]:
        vals: list[np.ndarray] = []
        wts: list[np.ndarray] = []
        total_ratio = total_mass = 0.0
        for j in range(age_to_index(P, age_lo), age_to_index(P, age_hi) + 1):
            for i in range(P.I):
                for zz in range(g7.shape[4]):
                    z_value = float(z_values[zz]) if zz < len(z_values) else 1.0
                    for nn in range(P.n_parity):
                        for cs in range(P.n_child_states):
                            if childless_only and current_child_bin_dt(nn, cs, dep_last, hcut, getattr(P, "child_state_mode", "shared_clock")) != 2:
                                continue
                            y = annual_gross_income_at_state(
                                P, i, j, z_value, children_at_home_count(nn, cs, P)
                            )
                            if y <= 0:
                                continue
                            tenures = [0] if renter_only else range(g7.shape[1])
                            for ten in tenures:
                                mass = g7[:, ten, i, j, zz, nn, cs]
                                if not np.any(mass > 0):
                                    continue
                                ratio = bg_arr / y
                                total_ratio += float(np.sum(mass * ratio))
                                total_mass += float(np.sum(mass))
                                positive = mass > 0
                                vals.append(ratio[positive])
                                wts.append(mass[positive])
        mean = total_ratio / max(total_mass, 1e-12)
        median = weighted_median_from_cells(vals, wts)
        return mean, median, total_mass

    samples = {
        "young_all_liquid_wealth_to_annual_gross_income_2530": sample_stats(
            25.0, 30.0, childless_only=False, renter_only=False
        ),
        "young_childless_liquid_wealth_to_annual_gross_income_2535": sample_stats(
            25.0, 35.0, childless_only=True, renter_only=False
        ),
        "young_childless_renter_liquid_wealth_to_annual_gross_income_2535": sample_stats(
            25.0, 35.0, childless_only=True, renter_only=True
        ),
    }
    for name, (mean, median, mass) in samples.items():
        setattr(stats, name, mean)
        setattr(stats, f"{name}_median", median)
        setattr(stats, f"{name}_mass", mass)


def add_annual_gross_old_wealth_moments(
    stats: SimpleNamespace,
    g: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
) -> None:
    """Add cross-sectional old-household wealth moments.

    The empirical sample contains living PSID reference persons aged 76--84;
    it is not a sample of decedents or realized estates. The model therefore
    uses the coherent beginning-of-period living-household balance sheet.
    An owner contributes liquid net worth plus gross housing value, ``b + pH``.
    The transaction wedge ``psi`` is not subtracted: it is neither current
    mortgage debt nor an empirical reduction in home equity. Income is annual
    gross income at the full Markov income state.
    """

    g_arr = np.asarray(g, dtype=float)
    bg_arr = np.asarray(bg, dtype=float).reshape(-1)
    ph_arr = np.asarray(ph, dtype=float).reshape(-1)
    if g_arr.ndim == 6:
        g7 = g_arr[:, :, :, :, None, :, :]
        z_values = np.array([1.0])
    elif g_arr.ndim == 7:
        g7 = g_arr
        z_values = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    else:
        return

    old_vals: list[np.ndarray] = []
    old_wts: list[np.ndarray] = []
    one_vals: list[np.ndarray] = []
    one_wts: list[np.ndarray] = []
    two_plus_vals: list[np.ndarray] = []
    two_plus_wts: list[np.ndarray] = []
    for j in range(int(P.J)):
        age = float(P.age_start) + j * float(P.da)
        in_old_tail = 76.0 <= age <= 84.0
        in_fertility_gap = 65.0 <= age <= 75.0
        if not (in_old_tail or in_fertility_gap):
            continue
        for i in range(int(P.I)):
            price = float(ph_arr[i])
            for zz in range(g7.shape[4]):
                z_value = float(z_values[zz]) if zz < z_values.size else 1.0
                income = annual_gross_income_at_state(P, i, j, z_value)
                if not np.isfinite(income) or income <= 0.0:
                    continue
                for ten in range(g7.shape[1]):
                    housing_value = price * float(P.H_own[ten - 1]) if ten > 0 else 0.0
                    ratio = (bg_arr + housing_value) / income
                    for nn in range(int(P.n_parity)):
                        for cs in range(int(P.n_child_states)):
                            mass = g7[:, ten, i, j, zz, nn, cs]
                            positive = np.isfinite(mass) & (mass > 0.0)
                            if not np.any(positive):
                                continue
                            values = ratio[positive]
                            weights = mass[positive]
                            if in_old_tail:
                                old_vals.append(values)
                                old_wts.append(weights)
                            if in_fertility_gap and nn == 1:
                                one_vals.append(values)
                                one_wts.append(weights)
                            elif in_fertility_gap and nn >= 2:
                                two_plus_vals.append(values)
                                two_plus_wts.append(weights)

    def pooled_quantile(value_cells: list[np.ndarray], weight_cells: list[np.ndarray], prob: float) -> float:
        if not value_cells:
            return float("nan")
        return float(weighted_quantile(np.concatenate(value_cells), np.concatenate(weight_cells), prob))

    old_p50 = pooled_quantile(old_vals, old_wts, 0.5)
    old_p90 = pooled_quantile(old_vals, old_wts, 0.9)
    one_p50 = pooled_quantile(one_vals, one_wts, 0.5)
    two_plus_p50 = pooled_quantile(two_plus_vals, two_plus_wts, 0.5)
    old_p90_p50 = old_p90 / max(old_p50, 1e-12)
    stats.old_total_wealth_to_annual_income_median_7684 = old_p50
    stats.old_total_wealth_to_annual_income_p90_7684 = old_p90
    stats.old_total_wealth_to_annual_income_p90_p50_7684 = old_p90_p50
    # Backward-compatible names for historical target systems.
    stats.old_total_estate_wealth_to_annual_income_median_7684 = old_p50
    stats.old_total_estate_wealth_to_annual_income_p90_7684 = old_p90
    stats.old_total_estate_wealth_to_annual_income_p90_p50_7684 = old_p90_p50
    stats.old_1_total_estate_wealth_to_annual_income_median_6575 = one_p50
    stats.old_2plus_total_estate_wealth_to_annual_income_median_6575 = two_plus_p50
    stats.old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575 = two_plus_p50 - one_p50


def add_old_nonhousing_income_share_moments(
    stats: SimpleNamespace,
    g: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
) -> None:
    """Add the PSID-disciplined old-age nonhousing retention share.

    ``old_nonhousing_ge_1x_income_share_6575`` is the mass share of households
    at ages 65-75 whose liquid wealth ``b`` is at least one year of annual
    gross income at the full Markov income state.  ``b`` enters raw, so
    negative balances count in the denominator but never the numerator.  The
    moment uses the same distribution object, age mapping, and positive-income
    filter as the estate moments, pooled across locations, income states,
    tenures, and children states.
    """

    g_arr = np.asarray(g, dtype=float)
    bg_arr = np.asarray(bg, dtype=float).reshape(-1)
    if g_arr.ndim == 6:
        g7 = g_arr[:, :, :, :, None, :, :]
        z_values = np.array([1.0])
    elif g_arr.ndim == 7:
        g7 = g_arr
        z_values = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    else:
        return

    share_mass = total_mass = 0.0
    for j in range(int(P.J)):
        age = float(P.age_start) + j * float(P.da)
        if not 65.0 <= age <= 75.0:
            continue
        for i in range(int(P.I)):
            for zz in range(g7.shape[4]):
                z_value = float(z_values[zz]) if zz < z_values.size else 1.0
                income = annual_gross_income_at_state(P, i, j, z_value)
                if not np.isfinite(income) or income <= 0.0:
                    continue
                at_least_one_year = (bg_arr / income) >= 1.0
                for ten in range(g7.shape[1]):
                    for nn in range(int(P.n_parity)):
                        for cs in range(int(P.n_child_states)):
                            mass = g7[:, ten, i, j, zz, nn, cs]
                            positive = np.isfinite(mass) & (mass > 0.0)
                            if not np.any(positive):
                                continue
                            share_mass += float(np.sum(mass[positive & at_least_one_year]))
                            total_mass += float(np.sum(mass[positive]))
    stats.old_nonhousing_ge_1x_income_share_6575 = (
        share_mass / total_mass if total_mass > 0.0 else float("nan")
    )


def add_old_wealth_income_moments(
    stats: SimpleNamespace, g: np.ndarray, P: SimpleNamespace, bg: np.ndarray, ph: np.ndarray,
) -> None:
    """Overwrite old-age wealth ratios using the full income-state support.

    The generic statistics routine receives an income-collapsed distribution.
    That is valid only when retirement income is common across income states;
    with ``retirement_income_z_scale != 0`` nonlinear means and medians must
    instead be computed before collapsing the Markov axis.
    """
    g_arr, bg_arr, ph_arr = np.asarray(g, float), np.asarray(bg, float).reshape(-1), np.asarray(ph, float).reshape(-1)
    if g_arr.ndim == 6:
        g7, z_values = g_arr[:, :, :, :, None, :, :], np.array([1.0])
    elif g_arr.ndim == 7:
        g7, z_values = g_arr, np.asarray(getattr(P, "z_grid", [1.0]), float).reshape(-1)
    else:
        raise ValueError("old-age wealth statistics require a six- or seven-dimensional distribution")
    cells: dict[str, tuple[list[np.ndarray], list[np.ndarray]]] = {
        name: ([], []) for name in ("nonhousing", "total", "parent_nonhousing", "childless_nonhousing", "parent_total", "childless_total")
    }
    for j in range(age_to_index(P, 65), age_to_index(P, 75) + 1):
        for i in range(int(P.I)):
            for zz in range(g7.shape[4]):
                income = annual_gross_income_at_state(P, i, j, float(z_values[zz]) if zz < z_values.size else 1.0)
                if not np.isfinite(income) or income <= 0.0:
                    continue
                for tenure in range(g7.shape[1]):
                    equity = (1.0 - float(P.psi)) * ph_arr[i] * float(P.H_own[tenure - 1]) if tenure else 0.0
                    nonhousing, total = bg_arr / income, (bg_arr + equity) / income
                    for parity in range(int(P.n_parity)):
                        for child_state in range(int(P.n_child_states)):
                            mass = g7[:, tenure, i, j, zz, parity, child_state]
                            keep = np.isfinite(mass) & (mass > 0.0)
                            if not np.any(keep):
                                continue
                            weight, nr, tr = mass[keep], nonhousing[keep], total[keep]
                            for key, values in (("nonhousing", nr), ("total", tr)):
                                cells[key][0].append(values); cells[key][1].append(weight)
                            group = (
                                "parent"
                                if parity > 0
                                else "childless"
                                if child_state in readiness_childless_states(P)
                                else None
                            )
                            if group is not None:
                                for key, values in ((f"{group}_nonhousing", nr), (f"{group}_total", tr)):
                                    cells[key][0].append(values); cells[key][1].append(weight)
    def mean(name: str) -> float:
        values, weights = cells[name]
        return sum(float(np.sum(v*w)) for v,w in zip(values,weights)) / max(sum(float(np.sum(w)) for w in weights), 1e-12)
    def median(name: str) -> float:
        values, weights = cells[name]
        return weighted_median_from_cells(values, weights)
    for stem, value in (("old_nonhousing_wealth_to_income_6575", mean("nonhousing")), ("old_total_wealth_to_income_6575", mean("total")), ("old_parent_nonhousing_wealth_to_income_6575", mean("parent_nonhousing")), ("old_childless_nonhousing_wealth_to_income_6575", mean("childless_nonhousing")), ("old_parent_total_wealth_to_income_6575", mean("parent_total")), ("old_childless_total_wealth_to_income_6575", mean("childless_total")), ("old_nonhousing_wealth_to_income_median_6575", median("nonhousing")), ("old_total_wealth_to_income_median_6575", median("total")), ("old_parent_nonhousing_wealth_to_income_median_6575", median("parent_nonhousing")), ("old_childless_nonhousing_wealth_to_income_median_6575", median("childless_nonhousing")), ("old_parent_total_wealth_to_income_median_6575", median("parent_total")), ("old_childless_total_wealth_to_income_median_6575", median("childless_total"))):
        setattr(stats, stem, value)
    stats.old_parent_childless_nonhousing_wealth_to_income_gap_6575 = stats.old_parent_nonhousing_wealth_to_income_6575 - stats.old_childless_nonhousing_wealth_to_income_6575
    stats.old_parent_childless_total_wealth_to_income_gap_6575 = stats.old_parent_total_wealth_to_income_6575 - stats.old_childless_total_wealth_to_income_6575
    stats.old_parent_childless_nonhousing_wealth_to_income_median_gap_6575 = stats.old_parent_nonhousing_wealth_to_income_median_6575 - stats.old_childless_nonhousing_wealth_to_income_median_6575
    stats.old_parent_childless_total_wealth_to_income_median_gap_6575 = stats.old_parent_total_wealth_to_income_median_6575 - stats.old_childless_total_wealth_to_income_median_6575


def compute_markov_statistics(
    g: np.ndarray,
    fp: np.ndarray,
    lp: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
    hR: np.ndarray,
    asset_g: np.ndarray | None = None,
    bequest_g: np.ndarray | None = None,
    bp_pol: np.ndarray | None = None,
) -> SimpleNamespace:
    z_grid, z_weights, Pi_z = income_transition_values(P)
    asset_dist = g if asset_g is None else np.asarray(asset_g, dtype=float)
    if asset_dist.shape != g.shape:
        raise ValueError("asset_g must have the same shape as the realized current distribution")
    g_total = np.sum(g, axis=4)
    asset_total = np.sum(asset_dist, axis=4)
    hR_total = collapse_markov_policy(hR, g, z_weights)
    fp_total = collapse_markov_fertility_probs(fp, g, z_weights, P)
    lp_total = collapse_markov_location_probs(lp, g, z_weights)
    stats = compute_statistics(g_total, fp_total, lp_total, P, bg, ph, hR_total, asset_g=asset_total)
    # Correct the nonlinear renter room moments: compute_statistics applied the
    # >=6 / cap-threshold indicators and the renter median to the income-collapsed
    # renter policy, which understates threshold shares (Jensen on a nonlinear
    # operator). Recompute them on the full income-resolved distribution.
    for _name, _val in markov_renter_room_moments(g, hR, P).items():
        setattr(stats, _name, _val)
    if bequest_g is None or bp_pol is None:
        raise ValueError(
            "Markov statistics require the post-transaction distribution and "
            "post-saving policy for the at-death bequest flow"
        )
    add_aggregate_wealth_bequest_flow_moments(
        stats, asset_dist, bequest_g, bp_pol, P, bg, ph
    )
    stats.bequest_moment_timing = "post_saving_at_death"
    add_annual_gross_liquid_wealth_moments(stats, asset_dist, P, bg)
    add_annual_gross_old_wealth_moments(stats, asset_dist, P, bg, ph)
    add_old_nonhousing_income_share_moments(stats, asset_dist, P, bg)
    if float(getattr(P, "retirement_income_z_scale", 0.0)) != 0.0:
        add_old_wealth_income_moments(stats, asset_dist, P, bg, ph)
    Nz = len(z_grid)
    nt = 1 + P.n_house
    stats.income_state_mass = np.zeros(Nz)
    stats.own_rate_by_income_type = np.zeros(Nz)
    stats.mean_fertility_by_income_type = np.zeros(Nz)
    stats.housing_demand_by_income_type = np.zeros((Nz, P.I))
    for zz in range(Nz):
        gz = g[:, :, :, :, zz, :, :]
        mz = float(np.sum(gz))
        stats.income_state_mass[zz] = mz
        stats.own_rate_by_income_type[zz] = float(np.sum(gz[:, 1:, :, :, :, :]) / max(mz, 1e-12))
        mp = float(np.sum(g[:, :, :, P.A_f_end :, zz, :, :]))
        if mp > 1e-12:
            mean_n = 0.0
            for nn in range(P.n_parity):
                mean_n += nn * float(np.sum(g[:, :, :, P.A_f_end :, zz, nn, :])) / mp
            stats.mean_fertility_by_income_type[zz] = mean_n
        for i in range(P.I):
            Hd = 0.0
            for j in range(P.J):
                for nn in range(P.n_parity):
                    for cs in range(P.n_child_states):
                        Hd += float(np.sum(g[:, 0, i, j, zz, nn, cs] * hR[:, 0, i, j, zz, nn, cs]))
                        for ten in range(1, nt):
                            Hd += float(np.sum(g[:, ten, i, j, zz, nn, cs]) * P.H_own[ten - 1])
            stats.housing_demand_by_income_type[zz, i] = Hd / max(housing_demand_normalizer(P), 1e-12)

    if bool(getattr(P, "permanent_income_levels_enabled", False)):
        group_index = np.asarray(P.permanent_income_group_index, dtype=int).reshape(-1)
        level_values = np.asarray(P.permanent_income_level_values, dtype=float).reshape(-1)
        n_levels = level_values.size
        stats.permanent_income_level_values = level_values.copy()
        stats.permanent_income_completed_mass = np.zeros(n_levels)
        stats.permanent_income_childless_by_level = np.zeros(n_levels)
        stats.permanent_income_completed_fertility_by_level = np.zeros(n_levels)
        stats.permanent_income_own_rate_3055_by_level = np.zeros(n_levels)
        a30s = age_to_index(P, 30)
        a55e = age_to_index(P, 55)
        parity_weights = np.arange(P.n_parity, dtype=float)
        if (
            str(getattr(P, "fertility_units", "parity2x")) == "literal_topcode"
            and parity_weights.size > 0
        ):
            parity_weights[-1] = float(getattr(P, "tfr_top_bin_weight", parity_weights[-1]))
        for level in range(n_levels):
            states = np.flatnonzero(group_index == level)
            completed = np.take(
                g[:, :, :, P.A_f_end :, :, :, :],
                states,
                axis=4,
            )
            completed_mass = float(np.sum(completed))
            stats.permanent_income_completed_mass[level] = completed_mass
            stats.permanent_income_childless_by_level[level] = float(
                np.sum(completed[:, :, :, :, :, 0, :]) / max(completed_mass, 1e-12)
            )
            stats.permanent_income_completed_fertility_by_level[level] = float(
                sum(
                    parity_weights[nn]
                    * np.sum(completed[:, :, :, :, :, nn, :])
                    for nn in range(P.n_parity)
                )
                / max(completed_mass, 1e-12)
            )
            prime = np.take(
                g[:, :, :, a30s : a55e + 1, :, :, :],
                states,
                axis=4,
            )
            prime_mass = float(np.sum(prime))
            stats.permanent_income_own_rate_3055_by_level[level] = float(
                np.sum(prime[:, 1:, :, :, :, :, :]) / max(prime_mass, 1e-12)
            )
        stats.permanent_income_childless_high_minus_low = float(
            stats.permanent_income_childless_by_level[-1]
            - stats.permanent_income_childless_by_level[0]
        )
        stats.permanent_income_completed_fertility_high_minus_low = float(
            stats.permanent_income_completed_fertility_by_level[-1]
            - stats.permanent_income_completed_fertility_by_level[0]
        )
        stats.permanent_income_own_rate_3055_high_minus_low = float(
            stats.permanent_income_own_rate_3055_by_level[-1]
            - stats.permanent_income_own_rate_3055_by_level[0]
        )

    worker_income = worker_mass = 0.0
    young_income = young_mass = young_liquid = 0.0
    payroll_tax_revenue = 0.0
    period_scale = float(getattr(P, "period_years", getattr(P, "da", 1.0))) if bool(getattr(P, "scale_flows_to_period", False)) else 1.0
    a25s = age_to_index(P, 25)
    aye = age_to_index(P, 35)
    for j in range(P.J):
        for i in range(P.I):
            for zz, z_value in enumerate(z_grid):
                yj = income_at_state(P, i, j, float(z_value))
                mass = float(np.sum(g[:, :, i, j, zz, :, :]))
                if j < P.J_R:
                    if child_earnings_penalty_active(P):
                        for nn in range(int(P.n_parity)):
                            for cs in range(int(P.n_child_states)):
                                cell = float(np.sum(g[:, :, i, j, zz, nn, cs]))
                                ycell = penalized_income_at_state(
                                    P, i, j, float(z_value),
                                    children_at_home_count(nn, cs, P),
                                )
                                worker_income += ycell * cell
                                worker_mass += cell
                    else:
                        worker_income += yj * mass
                        worker_mass += mass
                    payroll_tax_revenue += period_scale * P.tau_pay * P.w_hat[i] * P.income_age_profile[j] * float(z_value) * mass
                if a25s <= j <= aye:
                    gm = np.sum(
                        asset_dist[
                            :, 0, i, j, zz, 0, readiness_childless_states(P)
                        ],
                        axis=-1,
                    )
                    mh = float(np.sum(gm))
                    if mh > 1e-15:
                        young_income += penalized_income_at_state(P, i, j, float(z_value), 0) * mh
                        young_mass += mh
                        young_liquid += float(np.sum(gm * bg))
    stats.mean_income = worker_income / max(worker_mass, 1e-12)
    stats.young_childless_renter_income = young_income / max(young_mass, 1e-12)
    stats.young_liquid_wealth = young_liquid / max(young_mass, 1e-12)
    stats.young_liquid_wealth_to_income = young_liquid / max(young_income, 1e-12)
    stats.wealth_to_income = getattr(stats, "mean_wealth_4555", 0.0) / max(stats.mean_income, 1e-12)
    stats.liquid_wealth_to_income = getattr(stats, "liquid_wealth_4555", 0.0) / max(stats.mean_income, 1e-12)
    stats.payroll_tax_revenue = payroll_tax_revenue
    stats.pension_outlays = P.pension * stats.retiree_mass_total
    stats.pension_budget_residual = stats.payroll_tax_revenue - stats.pension_outlays
    stats.implied_balanced_pension = stats.payroll_tax_revenue / max(stats.retiree_mass_total, 1e-12)
    stats.income_transition = Pi_z.copy()
    return stats


def add_aggregate_wealth_bequest_flow_moments(
    stats: SimpleNamespace,
    wealth_g: np.ndarray,
    death_choice_g: np.ndarray,
    bp_pol: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
) -> None:
    """Add repaired wealth and at-death bequest-flow moments.

    The stock numerator is beginning-of-period net worth for all living
    households.  The denominator is annual *gross* labor earnings at ages
    18--65: ``P.income`` is stored after the payroll wedge, so working-age
    earnings are divided by ``(1 - tau_pay)`` to recover the gross object
    that matches the PSID EARNINDRRC construction.  Death estates use
    post-saving ``b'`` and current tenure after the period's transaction.
    """
    wealth_arr = np.asarray(wealth_g, dtype=float)
    death_arr = np.asarray(death_choice_g, dtype=float)
    bp_arr = np.asarray(bp_pol, dtype=float)
    if wealth_arr.ndim != 7:
        return
    if death_arr.shape != wealth_arr.shape or bp_arr.shape != wealth_arr.shape:
        raise ValueError("wealth_g, death_choice_g, and bp_pol must share the full income-resolved state shape")
    bg_arr = np.asarray(bg, dtype=float).reshape(-1)
    ph_arr = np.asarray(ph, dtype=float).reshape(-1)
    z_values = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
    gross_up = 1.0 / max(1.0 - float(getattr(P, "tau_pay", 0.0)), 1e-12)
    aggregate_wealth = aggregate_gross_labor_earnings = annual_bequest_flow = 0.0
    wealth_by_age = np.zeros(int(P.J), dtype=float)
    gross_labor_earnings_by_age = np.zeros(int(P.J), dtype=float)
    for j in range(int(P.J)):
        if bool(getattr(P, "use_age_survival", False)) and j < int(P.J) - 1:
            death_probability = 1.0 - float(P.survival_probs[j])
        elif j == int(P.J) - 1:
            death_probability = 1.0
        else:
            death_probability = 0.0
        for i in range(int(P.I)):
            for zz, z_value in enumerate(z_values):
                state_mass = float(np.sum(wealth_arr[:, :, i, j, zz, :, :]))
                if j < int(P.J_R):
                    if child_earnings_penalty_active(P):
                        for nn in range(int(P.n_parity)):
                            for cs in range(int(P.n_child_states)):
                                cell_mass = float(np.sum(wealth_arr[:, :, i, j, zz, nn, cs]))
                                cell_earn = (
                                    penalized_income_at_state(
                                        P, i, j, float(z_value),
                                        children_at_home_count(nn, cs, P),
                                    )
                                    - float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
                                ) * gross_up / max(period_years, 1e-12)
                                aggregate_gross_labor_earnings += cell_earn * cell_mass
                                gross_labor_earnings_by_age[j] += cell_earn * cell_mass
                    else:
                        gross_earnings = float(P.income[i, j]) * float(z_value) * gross_up / max(period_years, 1e-12)
                        aggregate_gross_labor_earnings += gross_earnings * state_mass
                        gross_labor_earnings_by_age[j] += gross_earnings * state_mass
                for ten in range(wealth_arr.shape[1]):
                    housing_value = float(ph_arr[i]) * float(P.H_own[ten - 1]) if ten > 0 else 0.0
                    mass_by_asset = np.sum(wealth_arr[:, ten, i, j, zz, :, :], axis=(1, 2))
                    total_wealth = float(np.sum(mass_by_asset * (bg_arr + housing_value)))
                    aggregate_wealth += total_wealth
                    wealth_by_age[j] += total_wealth
                    # Death estates are gross by default; the explicit flow
                    # switch or an active estate transfer deducts selling cost.
                    estate_hv = housing_value
                    if ten > 0 and estate_flow_net_active(P):
                        estate_hv = estate_housing_value(
                            P, float(ph_arr[i]), float(P.H_own[ten - 1]), for_accounting=True
                        )
                    estate = bp_arr[:, ten, i, j, zz, :, :] + estate_hv
                    cell_mass = death_arr[:, ten, i, j, zz, :, :]
                    estate_sum = np.sum(cell_mass * np.maximum(estate, 0.0))
                    if bool(getattr(P, "native_due_stayer_credit", False)):
                        stay_mass = P._g_stay_distribution[:, ten, i, j, zz, :, :]
                        stay_estate = P._bp_pol_stay[:, ten, i, j, zz, :, :] + estate_hv
                        estate_sum += np.sum(stay_mass * (np.maximum(stay_estate, 0.0) - np.maximum(estate, 0.0)))
                    annual_bequest_flow += death_probability * float(estate_sum) / max(period_years, 1e-12)
    stats.aggregate_wealth = aggregate_wealth
    stats.aggregate_annual_gross_labor_earnings = aggregate_gross_labor_earnings
    stats.aggregate_wealth_to_annual_gross_labor_earnings = aggregate_wealth / max(aggregate_gross_labor_earnings, 1e-12)
    stats.annual_bequest_flow = annual_bequest_flow
    stats.annual_bequest_flow_to_aggregate_wealth = annual_bequest_flow / max(aggregate_wealth, 1e-12)
    stats.aggregate_wealth_by_age = wealth_by_age
    stats.aggregate_annual_gross_labor_earnings_by_age = gross_labor_earnings_by_age


def compute_markov_eq_stats(g: np.ndarray, P: SimpleNamespace, bg: np.ndarray, ph: np.ndarray, hR: np.ndarray) -> SimpleNamespace:
    """Minimal Markov-income statistics needed for price clearing."""

    _ = bg, ph
    nt = 1 + P.n_house
    npar = P.n_parity
    tm = float(np.sum(g))
    stats = SimpleNamespace()
    stats.own_rate = float(np.sum(g[:, 1:, :, :, :, :, :]) / max(tm, 1e-12))
    stats.pop_share = np.zeros(P.I)
    stats.housing_demand = np.zeros(P.I)
    norm = housing_demand_normalizer(P)
    for i in range(P.I):
        gi = g[:, :, i, :, :, :, :]
        stats.pop_share[i] = float(np.sum(gi)) / max(tm, 1e-12)
        renter_demand = float(np.sum(gi[:, 0, :, :, :, :] * hR[:, 0, i, :, :, :, :]))
        owner_demand = 0.0
        for ten in range(1, nt):
            owner_demand += float(np.sum(gi[:, ten, :, :, :, :])) * float(P.H_own[ten - 1])
        stats.housing_demand[i] = (renter_demand + owner_demand) / max(norm, 1e-12)
    mp = float(np.sum(g[:, :, :, P.A_f_end :, :, :, :]))
    stats.parity_dist = np.zeros(npar)
    for nn in range(npar):
        stats.parity_dist[nn] = np.sum(g[:, :, :, P.A_f_end :, :, nn, :]) / max(mp, 1e-12)
    stats.mean_parity = float(np.sum(np.arange(npar) * stats.parity_dist))
    return stats


def age_to_index(P: SimpleNamespace, age: float) -> int:
    idx = int(round((float(age) - float(P.age_start)) / max(float(P.da), 1e-12)))
    return int(np.clip(idx, 0, P.J - 1))


def compute_statistics(
    g: np.ndarray,
    fp: np.ndarray,
    lp: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
    hR: np.ndarray,
    asset_g: np.ndarray | None = None,
) -> SimpleNamespace:
    J = P.J
    fec = get_fecundity_by_age(P)
    I = P.I
    Nb = len(bg)
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    childless_states = readiness_childless_states(P)
    settled_childless_state = readiness_settled_state(P)
    asset_dist = g if asset_g is None else np.asarray(asset_g, dtype=float)
    if asset_dist.shape != g.shape:
        raise ValueError("asset_g must have the same shape as the realized current distribution")
    tm = float(np.sum(g))
    stats = SimpleNamespace()
    stats.own_rate = float(np.sum(g[:, 1:, :, :, :, :]) / max(tm, 1e-12))
    stats.pop_share = np.zeros(I)
    stats.own_by_loc = np.zeros(I)
    stats.housing_demand = np.zeros(I)
    for i in range(I):
        pi = float(np.sum(g[:, :, i, :, :, :]))
        stats.pop_share[i] = pi / max(tm, 1e-12)
        stats.own_by_loc[i] = np.sum(g[:, 1:, i, :, :, :]) / max(pi, 1e-12)
        Hd = 0.0
        for j in range(J):
            for nn in range(npar):
                for cs in range(ncs):
                    Hd += float(np.sum(g[:, 0, i, j, nn, cs] * hR[:, 0, i, j, nn, cs]))
                    for ten in range(1, nt):
                        Hd += float(np.sum(g[:, ten, i, j, nn, cs]) * P.H_own[ten - 1])
        stats.housing_demand[i] = Hd / housing_demand_normalizer(P)

    stats.worker_mass_by_loc = np.array([np.sum(g[:, :, i, : P.J_R, :, :]) for i in range(I)])
    stats.retiree_mass_by_loc = np.array([np.sum(g[:, :, i, P.J_R :, :, :]) for i in range(I)])
    stats.worker_mass_total = float(np.sum(stats.worker_mass_by_loc))
    stats.retiree_mass_total = float(np.sum(stats.retiree_mass_by_loc))
    stats.parity_dist = np.zeros(npar)
    mp = float(np.sum(g[:, :, :, P.A_f_end :, :, :]))
    for nn in range(npar):
        stats.parity_dist[nn] = np.sum(g[:, :, :, P.A_f_end :, nn, :]) / max(mp, 1e-12)
    stats.mean_parity = float(np.sum(np.arange(npar) * stats.parity_dist))
    stats.own_by_parity = np.zeros(npar)
    for nn in range(npar):
        mn = float(np.sum(g[:, :, :, :, nn, :]))
        stats.own_by_parity[nn] = np.sum(g[:, 1:, :, :, nn, :]) / max(mn, 1e-12)
    stats.mean_parity_by_loc = np.zeros(I)
    stats.frac_childless_by_loc = np.zeros(I)
    for i in range(I):
        mip = float(np.sum(g[:, :, i, P.A_f_end :, :, :]))
        if mip > 1e-12:
            stats.mean_parity_by_loc[i] = sum(
                nn * np.sum(g[:, :, i, P.A_f_end :, nn, :]) / mip for nn in range(npar)
            )
            stats.frac_childless_by_loc[i] = np.sum(g[:, :, i, P.A_f_end :, 0, :]) / mip
    stats.own_by_age = np.zeros(J)
    for jj in range(J):
        gj = g[:, :, :, jj, :, :]
        mj = float(np.sum(gj))
        if mj > 1e-12:
            stats.own_by_age[jj] = np.sum(gj[:, 1:, :, :, :]) / mj
    stats.child_state_dist = np.zeros((J, ncs))
    for jj in range(J):
        mj = float(np.sum(g[:, :, :, jj, :, :]))
        if mj > 1e-12:
            for cs in range(ncs):
                stats.child_state_dist[jj, cs] = np.sum(g[:, :, :, jj, :, cs]) / mj
    stats.fert_by_age = np.zeros(J)
    for j in range(P.A_f_start - 1, P.A_f_end):
        mj = float(np.sum(g[:, :, :, j, 0, childless_states]))
        if mj > 1e-12:
            En = 0.0
            for i in range(I):
                for ten in range(nt):
                    gs = g[:, ten, i, j, 0, settled_childless_state]
                    nz = gs > 1e-15
                    if not np.any(nz):
                        continue
                    pr = fp[:, ten, i, j, :]
                    En += float(np.sum(gs[nz] * (pr[nz, :] @ np.arange(npar))))
            stats.fert_by_age[j] = fec[j] * En / mj

    a22s = age_to_index(P, 22)
    a25s = age_to_index(P, 25)
    a45e = age_to_index(P, 45)
    a30s = age_to_index(P, 30)
    a55e = age_to_index(P, 55)
    a65s = age_to_index(P, 65)
    a75e = age_to_index(P, 75)
    asw = age_to_index(P, 45)
    aew = age_to_index(P, 55)
    newparent_cs = (
        list(range(1, P.n_child_states))
        if independent_child_maturation_active(P)
        else list(range(1, min(P.n_child_stages + 1, 3)))
    )

    if child_earnings_penalty_active(P):
        ti = 0.0
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    cell_mass = float(np.sum(g[:, :, i, : P.J_R, nn, cs]))
                    ti += (
                        float(P.income[i, 0])
                        * child_earnings_multiplier(
                            P, 0, children_at_home_count(nn, cs, P)
                        )
                        * cell_mass
                    )
    else:
        ti = sum(P.income[i, 0] * stats.worker_mass_by_loc[i] for i in range(I))
    tmw = float(np.sum(stats.worker_mass_by_loc))
    mean_income = ti / max(tmw, 1e-12)
    tw = tm4 = 0.0
    for jj in range(asw, aew + 1):
        for i in range(I):
            for ten in range(nt):
                heq = (1 - P.psi) * ph[i] * P.H_own[ten - 1] if ten > 0 else 0.0
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = asset_dist[:, ten, i, jj, nn, cs]
                        tw += float(np.sum(gs * (bg + heq)))
                        tm4 += float(np.sum(gs))
    stats.mean_wealth_4555 = tw / max(tm4, 1e-12)
    stats.mean_income = mean_income
    stats.wealth_to_income = stats.mean_wealth_4555 / max(mean_income, 1e-12)
    stats.entry_mass_by_loc = np.array([np.sum(g[:, :, i, 0, :, :]) for i in range(I)])
    append_pension_budget_stats(stats, g, P)

    tmv = tas = 0.0
    for j in range(J - 1):
        for io in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, io, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        sp = lp[:, ten, io, io, j, nn, cs]
                        tmv += float(np.sum(gs * (1 - sp)))
                        tas += mh
    stats.migration_rate = tmv / max(tas, 1e-12)
    tmv2 = tas2 = 0.0
    for j in range(a22s, min(a45e, J - 1) + 1):
        for io in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, io, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        sp = lp[:, ten, io, io, j, nn, cs]
                        tmv2 += float(np.sum(gs * (1 - sp)))
                        tas2 += mh
    stats.migration_rate_2245 = tmv2 / max(tas2, 1e-12)

    tl = tml = 0.0
    for jj in range(asw, aew + 1):
        for i in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = asset_dist[:, ten, i, jj, nn, cs]
                        tl += float(np.sum(gs * bg))
                        tml += float(np.sum(gs))
    stats.liquid_wealth_4555 = tl / max(tml, 1e-12)
    stats.liquid_wealth_to_income = stats.liquid_wealth_4555 / max(mean_income, 1e-12)

    aye = age_to_index(P, 35)
    young_block = g[:, :, :, a25s : aye + 1, :, :]
    yt = float(np.sum(young_block))
    yo = float(np.sum(young_block[:, 1:, :, :, :, :]))
    stats.young_own_rate = yo / max(yt, 1e-12)
    ylw = yinc = ycm = 0.0
    for jj in range(a25s, aye + 1):
        for i in range(I):
            gs = np.sum(
                asset_dist[:, 0, i, jj, 0, childless_states],
                axis=-1,
            )
            mh = float(np.sum(gs))
            if mh < 1e-15:
                continue
            ylw += float(np.sum(gs * bg))
            if child_earnings_penalty_active(P):
                yinc += penalized_income_at_state(P, i, jj, 1.0, 0) * mh
            else:
                yinc += P.income[i, jj] * mh
            ycm += mh
    stats.young_liquid_wealth = ylw / max(ycm, 1e-12)
    stats.young_childless_renter_income = yinc / max(ycm, 1e-12)
    stats.young_liquid_wealth_to_income = ylw / max(yinc, 1e-12)

    stats.attempt_hazard_by_age = np.zeros(J)
    stats.first_birth_hazard_by_age = np.zeros(J)
    first_birth_flows = np.zeros(J)
    tba = tfb = 0.0
    for j in range(P.A_f_start - 1, P.A_f_end):
        # Birth decisions cover the four-year interval beginning at the state
        # age.  Report timing at the interval midpoint so the model statistic
        # can be measured with the identical binned-age operator in NCHS data.
        ra = P.age_start + (j + 0.5) * P.da
        childless_mass = float(np.sum(g[:, :, :, j, 0, childless_states]))
        attempted = 0.0
        for i in range(I):
            for ten in range(nt):
                gs = g[:, ten, i, j, 0, settled_childless_state]
                nz = gs > 1e-15
                if not np.any(nz):
                    continue
                attempted += float(np.sum(gs[nz] * (1 - fp[nz, ten, i, j, 0])))
                pb = fec[j] * (1 - fp[nz, ten, i, j, 0])
                wb = float(np.sum(gs[nz] * pb))
                tba += ra * wb
                tfb += wb
        if childless_mass > 1e-12:
            stats.attempt_hazard_by_age[j] = attempted / childless_mass
        stats.first_birth_hazard_by_age[j] = fec[j] * stats.attempt_hazard_by_age[j]
        first_birth_flows[j] = childless_mass * stats.first_birth_hazard_by_age[j]
    if birth_count.enabled(P):
        # Full income-resolved pre-birth source distribution is required:
        # inverse binary reconstruction fails with count-valued increments.
        pre = P.birth_count_pre_distribution
        action = P.birth_count_action_probs
        realized = P.birth_count_realized_probs
        tba = tfb = 0.0
        for j in range(P.A_f_start - 1, P.A_f_end):
            source = pre[:, :, :, j, :, 0, 0]
            risk = float(np.sum(source))
            attempted = float(np.sum(source * np.sum(action[:, :, :, j, :, 0, 0, 1:], axis=-1)))
            born = float(np.sum(source * np.sum(realized[:, :, :, j, :, 0, 0, 1:], axis=-1)))
            expected = float(np.sum(source * (realized[:, :, :, j, :, 0, 0, :] @ np.arange(4))))
            stats.attempt_hazard_by_age[j] = attempted / max(risk, 1e-12)
            stats.first_birth_hazard_by_age[j] = born / max(risk, 1e-12)
            stats.fert_by_age[j] = expected / max(risk, 1e-12)
            first_birth_flows[j] = born
            tba += (P.age_start + (j + .5) * P.da) * born
            tfb += born
    stats.mean_age_first_birth = tba / max(tfb, 1e-12)
    first_birth_total = float(np.sum(first_birth_flows))
    stats.first_birth_age_distribution = (
        first_birth_flows / first_birth_total if first_birth_total > 0.0 else np.zeros(J)
    )
    ages = float(P.age_start) + (
        np.arange(J, dtype=float) + 0.5
    ) * float(P.da)
    stats.share_first_births_age30plus = float(np.sum(stats.first_birth_age_distribution[ages >= 30.0]))
    chosen_survival = clock_survival = 1.0
    for j in range(P.A_f_start - 1, P.A_f_end):
        chosen_survival *= 1.0 - stats.attempt_hazard_by_age[j]
        clock_survival *= 1.0 - stats.first_birth_hazard_by_age[j]
    stats.childless_chosen_45 = float(chosen_survival)
    stats.childless_clock_45 = float(max(clock_survival - chosen_survival, 0.0))
    mg1 = float(np.sum(g[:, :, :, P.A_f_end :, 1:, :]))
    stats.parity_progression_1to2 = float(np.sum(g[:, :, :, P.A_f_end :, 2:, :]) / max(mg1, 1e-12))

    ic = 1 if I > 1 else 0
    mcy = float(np.sum(g[:, :, ic, a25s : aye + 1, 0, :]))
    mct = float(np.sum(g[:, :, :, a25s : aye + 1, 0, :]))
    stats.center_share_childless_young = mcy / max(mct, 1e-12)
    stats.center_share_parents = float(np.sum(g[:, :, ic, :, 1:, :]) / max(np.sum(g[:, :, :, :, 1:, :]), 1e-12))
    mpo = float(np.sum(g[:, 1:, :, :, 1:, :]))
    mpa = float(np.sum(g[:, :, :, :, 1:, :]))
    stats.own_rate_parents = mpo / max(mpa, 1e-12)
    mco = float(np.sum(g[:, 1:, :, :, 0, :]))
    mca = float(np.sum(g[:, :, :, :, 0, :]))
    stats.own_rate_childless = mco / max(mca, 1e-12)
    stats.own_family_gap = stats.own_rate_parents - stats.own_rate_childless

    a34e = age_to_index(P, 34)
    a35s = age_to_index(P, 35)
    a44e = age_to_index(P, 44)
    prime_mass = float(np.sum(g[:, :, :, a30s : a55e + 1, :, :]))
    prime_owner = float(np.sum(g[:, 1:, :, a30s : a55e + 1, :, :]))
    stats.own_rate_3055 = prime_owner / max(prime_mass, 1e-12)
    early_mass = float(np.sum(g[:, :, :, a25s : a34e + 1, :, :]))
    early_owner = float(np.sum(g[:, 1:, :, a25s : a34e + 1, :, :]))
    stats.own_rate_2534 = early_owner / max(early_mass, 1e-12)
    mid_mass = float(np.sum(g[:, :, :, a35s : a44e + 1, :, :]))
    mid_owner = float(np.sum(g[:, 1:, :, a35s : a44e + 1, :, :]))
    stats.own_rate_3544 = mid_owner / max(mid_mass, 1e-12)
    if I > 1:
        prime_mass_p = float(np.sum(g[:, :, 0, a30s : a55e + 1, :, :]))
        prime_mass_c = float(np.sum(g[:, :, 1, a30s : a55e + 1, :, :]))
        prime_owner_p = float(np.sum(g[:, 1:, 0, a30s : a55e + 1, :, :]))
        prime_owner_c = float(np.sum(g[:, 1:, 1, a30s : a55e + 1, :, :]))
        stats.own_gradient_3055 = prime_owner_p / max(prime_mass_p, 1e-12) - prime_owner_c / max(prime_mass_c, 1e-12)
    else:
        stats.own_gradient_3055 = 0.0
    nonparent_mass_2245 = float(
        np.sum(g[:, :, :, a22s : a45e + 1, 0, childless_states])
    )
    stats.center_share_nonparents_2245 = (
        float(
            np.sum(g[:, :, 1, a22s : a45e + 1, 0, childless_states])
            / max(nonparent_mass_2245, 1e-12)
        )
        if I > 1
        else 1.0
    )
    nonparent_mass_3055 = float(
        np.sum(g[:, :, :, a30s : a55e + 1, 0, childless_states])
    )
    nonparent_owner_3055 = float(
        np.sum(g[:, 1:, :, a30s : a55e + 1, 0, childless_states])
    )
    stats.own_rate_nonparents_3055 = nonparent_owner_3055 / max(nonparent_mass_3055, 1e-12)
    if not newparent_cs:
        stats.center_share_newparents_2245 = 0.0
        stats.own_rate_newparents_3055 = 0.0
    else:
        newparent_mass_2245 = float(np.sum(g[:, :, :, a22s : a45e + 1, 1:, newparent_cs]))
        stats.center_share_newparents_2245 = (
            float(np.sum(g[:, :, 1, a22s : a45e + 1, 1:, newparent_cs]) / max(newparent_mass_2245, 1e-12))
            if I > 1
            else 1.0
        )
        newparent_mass_3055 = float(np.sum(g[:, :, :, a30s : a55e + 1, 1:, newparent_cs]))
        newparent_owner_3055 = float(np.sum(g[:, 1:, :, a30s : a55e + 1, 1:, newparent_cs]))
        stats.own_rate_newparents_3055 = newparent_owner_3055 / max(newparent_mass_3055, 1e-12)
    stats.own_gap_newparent_nonparent_3055 = stats.own_rate_newparents_3055 - stats.own_rate_nonparents_3055

    old_mass = float(np.sum(g[:, :, :, a65s : a75e + 1, :, :]))
    old_owner = float(np.sum(g[:, 1:, :, a65s : a75e + 1, :, :]))
    stats.old_age_own_rate_6575 = old_owner / max(old_mass, 1e-12)
    old_parent_mass = float(np.sum(g[:, :, :, a65s : a75e + 1, 1:, :]))
    old_parent_owner = float(np.sum(g[:, 1:, :, a65s : a75e + 1, 1:, :]))
    old_childless_mass = float(
        np.sum(g[:, :, :, a65s : a75e + 1, 0, childless_states])
    )
    old_childless_owner = float(
        np.sum(g[:, 1:, :, a65s : a75e + 1, 0, childless_states])
    )
    stats.old_age_own_rate_parents_6575 = old_parent_owner / max(old_parent_mass, 1e-12)
    stats.old_age_own_rate_childless_6575 = old_childless_owner / max(old_childless_mass, 1e-12)
    stats.old_age_parent_childless_gap_6575 = stats.old_age_own_rate_parents_6575 - stats.old_age_own_rate_childless_6575
    old_nonhousing_wealth = old_income = 0.0
    old_parent_nonhousing_wealth = old_parent_income = 0.0
    old_childless_nonhousing_wealth = old_childless_income = 0.0
    old_total_wealth = old_parent_total_wealth = old_childless_total_wealth = 0.0
    old_nonhousing_ratio_vals: list[np.ndarray] = []
    old_nonhousing_ratio_wts: list[np.ndarray] = []
    old_total_ratio_vals: list[np.ndarray] = []
    old_total_ratio_wts: list[np.ndarray] = []
    old_parent_nonhousing_ratio_vals: list[np.ndarray] = []
    old_parent_nonhousing_ratio_wts: list[np.ndarray] = []
    old_childless_nonhousing_ratio_vals: list[np.ndarray] = []
    old_childless_nonhousing_ratio_wts: list[np.ndarray] = []
    old_parent_total_ratio_vals: list[np.ndarray] = []
    old_parent_total_ratio_wts: list[np.ndarray] = []
    old_childless_total_ratio_vals: list[np.ndarray] = []
    old_childless_total_ratio_wts: list[np.ndarray] = []
    for jj in range(a65s, a75e + 1):
        for i in range(I):
            yj = annual_gross_income_at_state(P, i, jj, 1.0)
            for ten in range(nt):
                home_equity = (1 - P.psi) * ph[i] * P.H_own[ten - 1] if ten > 0 else 0.0
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = asset_dist[:, ten, i, jj, nn, cs]
                        mass = float(np.sum(gs))
                        if mass < 1e-15:
                            continue
                        fin = float(np.sum(gs * bg))
                        total = fin + mass * home_equity
                        income = yj * mass
                        old_nonhousing_wealth += fin
                        old_total_wealth += total
                        old_income += income
                        positive = gs > 0
                        if np.any(positive) and yj > 0:
                            wts = gs[positive]
                            nonhousing_ratio = bg[positive] / yj
                            total_ratio = (bg[positive] + home_equity) / yj
                            old_nonhousing_ratio_vals.append(nonhousing_ratio)
                            old_nonhousing_ratio_wts.append(wts)
                            old_total_ratio_vals.append(total_ratio)
                            old_total_ratio_wts.append(wts)
                        if nn > 0:
                            old_parent_nonhousing_wealth += fin
                            old_parent_total_wealth += total
                            old_parent_income += income
                            if np.any(positive) and yj > 0:
                                old_parent_nonhousing_ratio_vals.append(nonhousing_ratio)
                                old_parent_nonhousing_ratio_wts.append(wts)
                                old_parent_total_ratio_vals.append(total_ratio)
                                old_parent_total_ratio_wts.append(wts)
                        elif cs in childless_states:
                            old_childless_nonhousing_wealth += fin
                            old_childless_total_wealth += total
                            old_childless_income += income
                            if np.any(positive) and yj > 0:
                                old_childless_nonhousing_ratio_vals.append(nonhousing_ratio)
                                old_childless_nonhousing_ratio_wts.append(wts)
                                old_childless_total_ratio_vals.append(total_ratio)
                                old_childless_total_ratio_wts.append(wts)
    def weighted_mean_from_cells(value_cells: list[np.ndarray], weight_cells: list[np.ndarray]) -> float:
        numerator = denominator = 0.0
        for values, weights in zip(value_cells, weight_cells):
            if values.size == 0 or weights.size == 0:
                continue
            numerator += float(np.sum(values * weights))
            denominator += float(np.sum(weights))
        return numerator / max(denominator, 1e-12)

    stats.old_nonhousing_wealth_to_income_6575 = weighted_mean_from_cells(
        old_nonhousing_ratio_vals, old_nonhousing_ratio_wts
    )
    stats.old_total_wealth_to_income_6575 = weighted_mean_from_cells(old_total_ratio_vals, old_total_ratio_wts)
    parent_nonhousing_ratio = weighted_mean_from_cells(
        old_parent_nonhousing_ratio_vals, old_parent_nonhousing_ratio_wts
    )
    childless_nonhousing_ratio = weighted_mean_from_cells(
        old_childless_nonhousing_ratio_vals, old_childless_nonhousing_ratio_wts
    )
    parent_total_ratio = weighted_mean_from_cells(old_parent_total_ratio_vals, old_parent_total_ratio_wts)
    childless_total_ratio = weighted_mean_from_cells(old_childless_total_ratio_vals, old_childless_total_ratio_wts)
    stats.old_parent_nonhousing_wealth_to_income_6575 = parent_nonhousing_ratio
    stats.old_childless_nonhousing_wealth_to_income_6575 = childless_nonhousing_ratio
    stats.old_parent_childless_nonhousing_wealth_to_income_gap_6575 = parent_nonhousing_ratio - childless_nonhousing_ratio
    stats.old_parent_total_wealth_to_income_6575 = parent_total_ratio
    stats.old_childless_total_wealth_to_income_6575 = childless_total_ratio
    stats.old_parent_childless_total_wealth_to_income_gap_6575 = parent_total_ratio - childless_total_ratio
    stats.old_nonhousing_wealth_to_income_median_6575 = weighted_median_from_cells(
        old_nonhousing_ratio_vals, old_nonhousing_ratio_wts
    )
    stats.old_total_wealth_to_income_median_6575 = weighted_median_from_cells(old_total_ratio_vals, old_total_ratio_wts)
    old_parent_nonhousing_median = weighted_median_from_cells(
        old_parent_nonhousing_ratio_vals, old_parent_nonhousing_ratio_wts
    )
    old_childless_nonhousing_median = weighted_median_from_cells(
        old_childless_nonhousing_ratio_vals, old_childless_nonhousing_ratio_wts
    )
    old_parent_total_median = weighted_median_from_cells(old_parent_total_ratio_vals, old_parent_total_ratio_wts)
    old_childless_total_median = weighted_median_from_cells(old_childless_total_ratio_vals, old_childless_total_ratio_wts)
    stats.old_parent_nonhousing_wealth_to_income_median_6575 = old_parent_nonhousing_median
    stats.old_childless_nonhousing_wealth_to_income_median_6575 = old_childless_nonhousing_median
    stats.old_parent_childless_nonhousing_wealth_to_income_median_gap_6575 = (
        old_parent_nonhousing_median - old_childless_nonhousing_median
    )
    stats.old_parent_total_wealth_to_income_median_6575 = old_parent_total_median
    stats.old_childless_total_wealth_to_income_median_6575 = old_childless_total_median
    stats.old_parent_childless_total_wealth_to_income_median_gap_6575 = (
        old_parent_total_median - old_childless_total_median
    )
    owner_mass_2545 = float(np.sum(asset_dist[:, 1:, :, a25s : a45e + 1, :, :]))
    stats.owner_neg_liquid_share_2545 = float(
        np.sum(asset_dist[bg < 0, 1:, :, a25s : a45e + 1, :, :]) / max(owner_mass_2545, 1e-12)
    )
    owner_mass_2534 = float(np.sum(asset_dist[:, 1:, :, a25s : a34e + 1, :, :]))
    stats.owner_neg_liquid_share_2534 = float(
        np.sum(asset_dist[bg < 0, 1:, :, a25s : a34e + 1, :, :]) / max(owner_mass_2534, 1e-12)
    )

    dep_last = P.n_child_stages
    hcut = int(getattr(P, "child_bin_high_cutoff", 2))
    renter_vals: list[np.ndarray] = []
    renter_wts: list[np.ndarray] = []
    owner_vals: list[np.ndarray] = []
    owner_wts: list[np.ndarray] = []
    for j in range(a25s, a45e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    if current_child_bin_dt(nn, cs, dep_last, hcut, getattr(P, "child_state_mode", "shared_clock")) != 2:
                        continue
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        renter_vals.append(hr[kr])
                        renter_wts.append(gr[kr])
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        ko = go > 0
                        if np.any(ko):
                            owner_vals.append(P.H_own[ten - 1] * np.ones(np.count_nonzero(ko)))
                            owner_wts.append(go[ko])
    stats.prime_childless_renter_median_rooms = weighted_median_from_cells(renter_vals, renter_wts)
    stats.prime_childless_owner_median_rooms = weighted_median_from_cells(owner_vals, owner_wts)
    renter_mass_3055_childless = 0.0
    renter_rooms_3055_childless = 0.0
    renter_rooms_ge6_3055_childless = 0.0
    owner_mass_3055_childless = 0.0
    owner_rooms_3055_childless = 0.0
    owner_rooms_ge6_3055_childless = 0.0
    owner_mass_3055_parent = 0.0
    owner_rooms_3055_parent = 0.0
    renter_mass_3055_parent = 0.0
    renter_rooms_3055_parent = 0.0
    parent_low_mass_3055 = 0.0
    parent_low_rooms_3055 = 0.0
    parent_high_mass_3055 = 0.0
    parent_high_rooms_3055 = 0.0
    for j in range(a30s, a55e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    child_bin = current_child_bin_dt(nn, cs, dep_last, hcut, getattr(P, "child_state_mode", "shared_clock"))
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        wr = gr[kr]
                        rr = hr[kr]
                        mass = float(np.sum(wr))
                        if child_bin == 2:
                            renter_mass_3055_childless += mass
                            renter_rooms_3055_childless += float(np.sum(wr * rr))
                            renter_rooms_ge6_3055_childless += float(np.sum(wr[rr >= 6.0 - 1e-8]))
                        elif child_bin > 2:
                            renter_mass_3055_parent += mass
                            renter_rooms_3055_parent += float(np.sum(wr * rr))
                            if child_bin == 3:
                                parent_low_mass_3055 += mass
                                parent_low_rooms_3055 += float(np.sum(wr * rr))
                            elif child_bin >= 4:
                                parent_high_mass_3055 += mass
                                parent_high_rooms_3055 += float(np.sum(wr * rr))
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        mo = float(np.sum(go))
                        if mo <= 1e-15:
                            continue
                        rooms = float(P.H_own[ten - 1])
                        if child_bin == 2:
                            owner_mass_3055_childless += mo
                            owner_rooms_3055_childless += mo * rooms
                            if rooms >= 6.0 - 1e-8:
                                owner_rooms_ge6_3055_childless += mo
                        elif child_bin > 2:
                            owner_mass_3055_parent += mo
                            owner_rooms_3055_parent += mo * rooms
                            if child_bin == 3:
                                parent_low_mass_3055 += mo
                                parent_low_rooms_3055 += mo * rooms
                            elif child_bin >= 4:
                                parent_high_mass_3055 += mo
                                parent_high_rooms_3055 += mo * rooms
    stats.prime30_55_childless_renter_mean_rooms = (
        renter_rooms_3055_childless / max(renter_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_mean_rooms = (
        owner_rooms_3055_childless / max(owner_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_minus_renter_mean_rooms = (
        stats.prime30_55_childless_owner_mean_rooms - stats.prime30_55_childless_renter_mean_rooms
    )
    stats.prime30_55_childless_renter_share_rooms_ge6 = (
        renter_rooms_ge6_3055_childless / max(renter_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_share_rooms_ge6 = (
        owner_rooms_ge6_3055_childless / max(owner_mass_3055_childless, 1e-12)
    )
    parent_owner_mean = owner_rooms_3055_parent / max(owner_mass_3055_parent, 1e-12)
    parent_renter_mean = renter_rooms_3055_parent / max(renter_mass_3055_parent, 1e-12)
    stats.prime30_55_parent_owner_minus_renter_mean_rooms = parent_owner_mean - parent_renter_mean
    stats.prime30_55_parent_3plus_minus_1to2_mean_rooms = (
        parent_high_rooms_3055 / max(parent_high_mass_3055, 1e-12)
        - parent_low_rooms_3055 / max(parent_low_mass_3055, 1e-12)
    )
    owner_all_mass_2545 = 0.0
    owner_le6_mass_2545 = 0.0
    owner_7to8_mass_2545 = 0.0
    owner_ge9_mass_2545 = 0.0
    renter_cap_mass_2545_all = 0.0
    renter_cap_mass_2545_child0 = 0.0
    renter_cap_mass_2545_child1 = 0.0
    renter_mass_2545_all = 0.0
    renter_mass_2545_child0 = 0.0
    renter_mass_2545_child1 = 0.0
    renter_rooms_2545_child0 = 0.0
    renter_rooms_2545_child1 = 0.0
    owner_mass_2545_child0 = 0.0
    owner_mass_2545_child1 = 0.0
    owner_rooms_2545_child0 = 0.0
    owner_rooms_2545_child1 = 0.0
    owner_rung_mass_2545 = np.zeros(P.n_house)
    for j in range(a25s, a45e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    child_bin = current_child_bin_dt(nn, cs, dep_last, hcut, getattr(P, "child_state_mode", "shared_clock"))
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        wr = gr[kr]
                        rr = hr[kr]
                        mass = float(np.sum(wr))
                        cap_mass = float(np.sum(wr[rr >= float(P.hR_max) - 1e-8]))
                        renter_mass_2545_all += mass
                        renter_cap_mass_2545_all += cap_mass
                        if child_bin == 2:
                            renter_mass_2545_child0 += mass
                            renter_cap_mass_2545_child0 += cap_mass
                            renter_rooms_2545_child0 += float(np.sum(wr * rr))
                        elif child_bin == 3:
                            renter_mass_2545_child1 += mass
                            renter_cap_mass_2545_child1 += cap_mass
                            renter_rooms_2545_child1 += float(np.sum(wr * rr))
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        mo = float(np.sum(go))
                        if mo <= 1e-15:
                            continue
                        rooms = float(P.H_own[ten - 1])
                        owner_all_mass_2545 += mo
                        owner_rung_mass_2545[ten - 1] += mo
                        if rooms <= 6.0 + 1e-8:
                            owner_le6_mass_2545 += mo
                        elif rooms <= 8.0 + 1e-8:
                            owner_7to8_mass_2545 += mo
                        else:
                            owner_ge9_mass_2545 += mo
                        if child_bin == 2:
                            owner_mass_2545_child0 += mo
                            owner_rooms_2545_child0 += mo * rooms
                        elif child_bin == 3:
                            owner_mass_2545_child1 += mo
                            owner_rooms_2545_child1 += mo * rooms
    for rung in range(P.n_house):
        setattr(
            stats,
            f"owner25_45_rung{rung + 1}_share",
            float(owner_rung_mass_2545[rung] / max(owner_all_mass_2545, 1e-12)),
        )
    stats.owner25_45_rooms_le6_share = float(owner_le6_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.owner25_45_rooms_7to8_share = float(owner_7to8_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.owner25_45_rooms_ge9_share = float(owner_ge9_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.renter25_45_all_cap_share = float(renter_cap_mass_2545_all / max(renter_mass_2545_all, 1e-12))
    stats.renter25_45_current0_cap_share = float(renter_cap_mass_2545_child0 / max(renter_mass_2545_child0, 1e-12))
    stats.renter25_45_current1_cap_share = float(renter_cap_mass_2545_child1 / max(renter_mass_2545_child1, 1e-12))
    stats.renter25_45_current0_mean_rooms = float(renter_rooms_2545_child0 / max(renter_mass_2545_child0, 1e-12))
    stats.renter25_45_current1_mean_rooms = float(renter_rooms_2545_child1 / max(renter_mass_2545_child1, 1e-12))
    stats.owner25_45_current0_mean_rooms = float(owner_rooms_2545_child0 / max(owner_mass_2545_child0, 1e-12))
    stats.owner25_45_current1_mean_rooms = float(owner_rooms_2545_child1 / max(owner_mass_2545_child1, 1e-12))
    stats.mean_housing_by_parity = np.zeros(npar)
    for nn in range(npar):
        th = mn = 0.0
        for i in range(I):
            for ten in range(nt):
                for j in range(J):
                    for cs in range(ncs):
                        gs = g[:, ten, i, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        if ten == 0:
                            th += float(np.sum(gs * hR[:, ten, i, j, nn, cs]))
                        else:
                            th += mh * P.H_own[ten - 1]
                        mn += mh
        stats.mean_housing_by_parity[nn] = th / max(mn, 1e-12)
    stats.housing_increment_1to2 = (
        stats.mean_housing_by_parity[2] - stats.mean_housing_by_parity[1] if npar >= 3 else 0.0
    )
    add_aggregate_wealth_gross_labor_diagnostics(stats, asset_dist, P, bg, ph)
    add_annual_gross_liquid_wealth_moments(stats, asset_dist, P, bg)
    return stats


def append_pension_budget_stats(stats: SimpleNamespace, g: np.ndarray, P: SimpleNamespace) -> None:
    income_profile = P.income_age_profile
    period_scale = float(getattr(P, "period_years", getattr(P, "da", 1.0))) if bool(getattr(P, "scale_flows_to_period", False)) else 1.0
    payroll_tax_revenue = 0.0
    for i in range(P.I):
        for j in range(P.J_R):
            mass_ij = float(np.sum(g[:, :, i, j, :, :]))
            payroll_tax_revenue += period_scale * P.tau_pay * P.w_hat[i] * income_profile[j] * mass_ij
    pension_outlays = P.pension * stats.retiree_mass_total
    stats.payroll_tax_revenue = payroll_tax_revenue
    stats.pension_outlays = pension_outlays
    stats.pension_budget_residual = payroll_tax_revenue - pension_outlays
    stats.implied_balanced_pension = payroll_tax_revenue / max(stats.retiree_mass_total, 1e-12)
    stats.pension = P.pension


def pack_solution_markov_income(
    V,
    c,
    h,
    bp,
    tc,
    tp,
    lp,
    fp,
    fv,
    g,
    st: SimpleNamespace,
    w,
    p,
    P: SimpleNamespace,
) -> SimpleNamespace:
    z_grid, z_weights, Pi_z = income_transition_values(P)
    sol = SimpleNamespace(
        V=V,
        c_pol=c,
        hR_pol=h,
        bp_pol=bp,
        tenure_choice=tc,
        tenure_probs=tp,
        loc_probs=lp,
        fert_probs=fp,
        fert2_probs=getattr(P, "_fert2_probs", None),
        joint_choice=getattr(P, "_joint_choice", None),
        bp_pol_stay=getattr(P, "_bp_pol_stay", None),
        c_pol_stay=getattr(P, "_c_pol_stay", None),
        fert_value=fv,
        g=g,
        g_collapsed=np.sum(g, axis=4),
        w_hat=w,
        p_eq=p,
        type_values=z_grid.copy(),
        type_weights=z_weights.copy(),
        income_transition=Pi_z.copy(),
    )
    for key, value in vars(st).items():
        setattr(sol, key, value)
    if birth_count.enabled(P):
        sol.birth_count_action_probs = P.birth_count_action_probs.copy()
        sol.birth_count_realized_probs = P.birth_count_realized_probs.copy()
        sol.birth_count_policy_axes = P.birth_count_policy_axes
    norm = housing_demand_normalizer(P)
    owner_by_size = np.zeros(P.n_house)
    renter_by_market = np.zeros(P.I)
    owner_by_market = np.zeros(P.I)
    for i in range(P.I):
        renter_by_market[i] = float(np.sum(g[:, 0, i, :, :, :, :] * h[:, 0, i, :, :, :, :])) / max(norm, 1e-12)
        for ten in range(1, 1 + P.n_house):
            demand = float(np.sum(g[:, ten, i, :, :, :, :])) * float(P.H_own[ten - 1]) / max(norm, 1e-12)
            owner_by_size[ten - 1] += demand
            owner_by_market[i] += demand
    user_cost = P.user_cost_rate * np.asarray(p, dtype=float).reshape(-1)
    supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply
    sol.owner_user_cost = user_cost
    sol.owner_asset_price = np.asarray(p, dtype=float).reshape(-1)
    sol.rental_demand_by_market = renter_by_market
    sol.owner_demand_by_market = owner_by_market
    sol.owner_demand_by_size = owner_by_size
    sol.rental_demand_by_size = renter_by_market.copy()
    sol.housing_supply = supply
    sol.aggregate_rental_demand = float(np.sum(renter_by_market))
    sol.aggregate_owner_demand = float(np.sum(owner_by_market))
    sol.aggregate_housing_demand = float(sol.aggregate_rental_demand + sol.aggregate_owner_demand)
    sol.aggregate_housing_supply = float(np.sum(supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((renter_by_market + owner_by_market - supply) / np.maximum(supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.young_owner_rate = float(getattr(sol, "young_own_rate", np.nan))
    sol.old_owner_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(getattr(sol, "parity_dist", np.array([np.nan]))[0])
    if not hasattr(sol, "own_rate_by_income_type"):
        sol.own_rate_by_income_type = np.full(len(z_grid), np.nan)
    if not hasattr(sol, "mean_fertility_by_income_type"):
        sol.mean_fertility_by_income_type = np.full(len(z_grid), np.nan)
    if not hasattr(sol, "housing_demand_by_income_type"):
        sol.housing_demand_by_income_type = np.full((len(z_grid), P.I), np.nan)
    return sol


def pack_fast_solution_markov_income(st: SimpleNamespace, p: np.ndarray, P: SimpleNamespace) -> SimpleNamespace:
    p_vec = np.asarray(p, dtype=float).reshape(-1)
    user_cost = P.user_cost_rate * p_vec
    supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply
    demand = np.asarray(st.housing_demand, dtype=float).reshape(-1)
    sol = SimpleNamespace()
    for key, value in vars(st).items():
        setattr(sol, key, value)
    if birth_count.enabled(P):
        sol.birth_count_action_probs = P.birth_count_action_probs.copy()
        sol.birth_count_realized_probs = P.birth_count_realized_probs.copy()
        sol.birth_count_policy_axes = P.birth_count_policy_axes
    sol.p_eq = p_vec
    sol.owner_user_cost = user_cost
    sol.owner_asset_price = p_vec
    sol.housing_supply = supply
    sol.aggregate_housing_demand = float(np.sum(demand))
    sol.aggregate_housing_supply = float(np.sum(supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((demand - supply) / np.maximum(supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(getattr(sol, "parity_dist", np.array([np.nan]))[0])
    return sol
