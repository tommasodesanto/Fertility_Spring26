"""First-pass solver for the intergenerational housing/fertility model."""

from __future__ import annotations

import time
from types import SimpleNamespace
from typing import Iterable

import numpy as np

from .parameters import apply_overrides, setup_parameters
from .utils import child_goods_cost, housing_need, make_asset_grid, mortgage_payment, interpolate_age_profile


NEG_INF = -1.0e18


def solve_model(
    overrides: SimpleNamespace | dict | None = None,
    *,
    mode: str = "benchmark",
    solve_prices: bool = True,
    verbose: bool = False,
) -> tuple[SimpleNamespace, SimpleNamespace]:
    """Solve the first-pass model.

    Returns `(solution, parameters)`. The first pass solves a stationary
    lifecycle cross-section. It is not a calibrated equilibrium.
    """

    P = apply_overrides(setup_parameters(mode), overrides)
    q = np.asarray(P.owner_user_cost, dtype=float).copy()
    trace: list[dict] = []
    best_q = q.copy()
    best_metric = np.inf
    t0 = time.perf_counter()

    if solve_prices:
        for it in range(int(P.max_iter_eq)):
            sol = solve_fixed_prices(P, q)
            excess = sol.owner_demand_by_size - P.owner_supply
            rel = excess / np.maximum(P.owner_supply, 1e-10)
            max_abs = float(np.max(np.abs(rel)))
            trace.append(
                {
                    "iter": it,
                    "max_abs_rel_excess": max_abs,
                    "owner_user_cost": q.copy(),
                    "owner_demand_by_size": sol.owner_demand_by_size.copy(),
                    "owner_supply": P.owner_supply.copy(),
                }
            )
            if verbose:
                print(f"iter={it:03d} max_rel_excess={max_abs:.4g} q={q}")
            if max_abs < best_metric:
                best_metric = max_abs
                best_q = q.copy()
            if max_abs < float(P.tol_eq):
                break
            q = q * np.exp(float(P.price_damping) * rel)
            q = np.clip(q, float(P.price_min), float(P.price_max))
        q = best_q
    else:
        sol = solve_fixed_prices(P, q)

    sol = solve_fixed_prices(P, q)
    sol.owner_user_cost = q
    sol.owner_asset_price = asset_prices(P, q)
    sol.elapsed_sec = time.perf_counter() - t0
    sol.price_trace = trace
    sol.best_max_abs_rel_excess = float(best_metric) if solve_prices else 0.0
    sol.converged = (not solve_prices) or (sol.best_max_abs_rel_excess < float(P.tol_eq))
    return sol, P


def solve_fixed_prices(P: SimpleNamespace, q_owner: Iterable[float]) -> SimpleNamespace:
    q_owner = np.asarray(q_owner, dtype=float)
    b_grid = make_asset_grid(P)
    age_profile = interpolate_age_profile(P)
    p_asset = asset_prices(P, q_owner)

    shape = (P.J + 1, len(b_grid), len(P.z_grid), P.Nn, P.Nt)
    V = np.zeros(shape, dtype=float)
    policy_shape = (P.J, len(b_grid), len(P.z_grid), P.Nn, P.Nt)
    pol_b = np.zeros(policy_shape, dtype=np.int64)
    pol_ten = np.zeros(policy_shape, dtype=np.int64)
    pol_n = np.zeros(policy_shape, dtype=np.int64)

    for a in range(P.J - 1, -1, -1):
        for ib, b in enumerate(b_grid):
            for iz, z in enumerate(P.z_grid):
                income = income_at(P, a, z, age_profile)
                for in_idx, n_state in enumerate(P.n_child_options):
                    for ten in range(P.Nt):
                        best_v = NEG_INF
                        best_b = 0
                        best_ten = 0
                        best_n = in_idx
                        feasible_ns = fertility_choices(P, a, in_idx)
                        for next_n_idx in feasible_ns:
                            children = int(P.n_child_options[next_n_idx])
                            for next_ten in range(P.Nt):
                                if not tenure_feasible(P, b, income, ten, next_ten, q_owner, p_asset):
                                    continue
                                flow_h = housing_flow_cost(P, next_ten, q_owner)
                                adjust = adjustment_cost(P, a, ten, next_ten, q_owner, p_asset)
                                child_cost = child_goods_cost(P, children, income)
                                resources = income + P.R * b - flow_h - adjust - child_cost
                                c = resources - b_grid
                                u = period_utility(P, c, next_ten, children)
                                cont = expected_continuation(V[a + 1], P, iz, next_n_idx, next_ten)
                                val = u + P.beta * cont
                                idx = int(np.argmax(val))
                                cand_v = float(val[idx])
                                if cand_v > best_v:
                                    best_v = cand_v
                                    best_b = idx
                                    best_ten = next_ten
                                    best_n = next_n_idx
                        V[a, ib, iz, in_idx, ten] = best_v
                        pol_b[a, ib, iz, in_idx, ten] = best_b
                        pol_ten[a, ib, iz, in_idx, ten] = best_ten
                        pol_n[a, ib, iz, in_idx, ten] = best_n

    dist = forward_distribution(P, b_grid, pol_b, pol_ten, pol_n)
    stats = compute_stats(P, b_grid, dist, pol_ten, pol_n, q_owner, p_asset, age_profile)
    stats.V = V
    stats.b_grid = b_grid
    stats.pol_b_idx = pol_b
    stats.pol_ten = pol_ten
    stats.pol_n = pol_n
    stats.dist = dist
    stats.owner_user_cost = q_owner.copy()
    stats.owner_asset_price = p_asset.copy()
    return stats


def asset_prices(P: SimpleNamespace, q_owner: np.ndarray) -> np.ndarray:
    return np.asarray(q_owner, dtype=float) / float(P.rho_property)


def income_at(P: SimpleNamespace, age_index: int, z: float, age_profile: np.ndarray) -> float:
    age = P.age_start + age_index
    if age < P.retire_age:
        return float(P.wage * age_profile[age_index] * z)
    return float(P.pension_replacement * P.wage * z)


def fertility_choices(P: SimpleNamespace, age_index: int, n_idx: int) -> np.ndarray:
    if age_index == P.fertility_choice_index:
        return np.arange(n_idx, P.Nn, dtype=np.int64)
    return np.array([n_idx], dtype=np.int64)


def tenure_feasible(
    P: SimpleNamespace,
    b: float,
    income: float,
    current_ten: int,
    next_ten: int,
    q_owner: np.ndarray,
    p_asset: np.ndarray,
) -> bool:
    if next_ten == 0:
        return True
    if current_ten == next_ten:
        return True
    k = next_ten - 1
    asset_value = float(p_asset[k] * P.owner_h[k])
    down_payment = (1.0 - float(P.phi_ltv)) * asset_value
    if max(float(b), 0.0) + 1e-12 < down_payment:
        return False
    debt = float(P.phi_ltv) * asset_value
    payment = mortgage_payment(debt, float(P.mortgage_rate), int(P.mortgage_maturity))
    return payment <= float(P.psi_pti) * max(float(income), 1e-12) + 1e-12


def housing_flow_cost(P: SimpleNamespace, ten: int, q_owner: np.ndarray) -> float:
    if ten == 0:
        return float(P.rent_user_cost * P.renter_h)
    k = ten - 1
    return float(q_owner[k] * P.owner_h[k])


def adjustment_cost(
    P: SimpleNamespace,
    age_index: int,
    current_ten: int,
    next_ten: int,
    q_owner: np.ndarray,
    p_asset: np.ndarray,
) -> float:
    if current_ten == next_ten:
        return 0.0
    cost = 0.0
    if next_ten > 0 and current_ten != next_ten:
        k_next = next_ten - 1
        cost += float(P.buyer_transaction_cost * p_asset[k_next] * P.owner_h[k_next])
    if current_ten > 0:
        k_cur = current_ten - 1
        cost += float(P.owner_move_cost * p_asset[k_cur] * P.owner_h[k_cur])
        downsizing_or_exit = next_ten == 0 or (next_ten > 0 and next_ten < current_ten)
        if age_index >= P.old_retention_index and downsizing_or_exit:
            cost += float(P.old_retention_wedge * q_owner[k_cur] * P.owner_h[k_cur])
    return cost


def period_utility(P: SimpleNamespace, c: np.ndarray, ten: int, children: int) -> np.ndarray:
    h = float(P.renter_h if ten == 0 else P.owner_h[ten - 1])
    slack = h - housing_need(P, children)
    out = np.full_like(c, NEG_INF, dtype=float)
    feasible = (c > float(P.c_min)) & (slack > 0.0)
    if not np.any(feasible):
        return out
    if abs(float(P.sigma) - 1.0) < 1e-10:
        u_c = np.log(c[feasible])
    else:
        sig = float(P.sigma)
        u_c = (np.power(c[feasible], 1.0 - sig) - 1.0) / (1.0 - sig)
    out[feasible] = u_c + float(P.alpha_h) * np.log(slack) + float(P.beta_n) * np.log1p(children)
    return out


def expected_continuation(
    V_next: np.ndarray,
    P: SimpleNamespace,
    z_idx: int,
    n_idx: int,
    ten: int,
) -> np.ndarray:
    # V_next has shape (Nb, Nz, Nn, Nt). Return expected continuation for each b' grid point.
    return np.tensordot(V_next[:, :, n_idx, ten], P.Pi_z[z_idx], axes=(1, 0))


def forward_distribution(
    P: SimpleNamespace,
    b_grid: np.ndarray,
    pol_b: np.ndarray,
    pol_ten: np.ndarray,
    pol_n: np.ndarray,
) -> np.ndarray:
    dist = np.zeros((P.J, len(b_grid), len(P.z_grid), P.Nn, P.Nt), dtype=float)
    b0 = int(np.argmin(np.abs(b_grid - float(P.b_entry))))
    for iz, mass in enumerate(P.z_dist):
        dist[0, b0, iz, 0, 0] += mass

    for a in range(P.J - 1):
        cur = dist[a]
        nz = np.nonzero(cur > 1e-15)
        for ib, iz, in_idx, ten in zip(*nz):
            mass = cur[ib, iz, in_idx, ten]
            nb = pol_b[a, ib, iz, in_idx, ten]
            nn = pol_n[a, ib, iz, in_idx, ten]
            nt = pol_ten[a, ib, iz, in_idx, ten]
            for izp, prob in enumerate(P.Pi_z[iz]):
                if prob > 0.0:
                    dist[a + 1, nb, izp, nn, nt] += mass * prob
    return dist


def compute_stats(
    P: SimpleNamespace,
    b_grid: np.ndarray,
    dist: np.ndarray,
    pol_ten: np.ndarray,
    pol_n: np.ndarray,
    q_owner: np.ndarray,
    p_asset: np.ndarray,
    age_profile: np.ndarray,
) -> SimpleNamespace:
    total = float(dist.sum())
    owner_mass_by_size_raw = np.zeros(P.K)
    for k in range(P.K):
        owner_mass_by_size_raw[k] = float(dist[:, :, :, :, k + 1].sum())
    owner_mass_by_size = owner_mass_by_size_raw / max(total, 1e-12)
    owner_demand_by_size = owner_mass_by_size * P.owner_h
    owner_mass = float(owner_mass_by_size_raw.sum())

    age_mass = dist.sum(axis=(1, 2, 3, 4))
    own_by_age = np.zeros(P.J)
    children_by_age = np.zeros(P.J)
    for a in range(P.J):
        m = age_mass[a]
        if m <= 0.0:
            continue
        own_by_age[a] = dist[a, :, :, :, 1:].sum() / m
        child_sum = 0.0
        for in_idx, n in enumerate(P.n_child_options):
            child_sum += n * dist[a, :, :, in_idx, :].sum()
        children_by_age[a] = child_sum / m

    post_fert = np.arange(P.fertility_choice_index + 1, P.J)
    post_mass = float(age_mass[post_fert].sum()) if len(post_fert) else 0.0
    completed = 0.0
    childless = 0.0
    if post_mass > 0.0:
        for in_idx, n in enumerate(P.n_child_options):
            m = float(dist[post_fert, :, :, in_idx, :].sum())
            completed += n * m
            if n == 0:
                childless += m
        completed /= post_mass
        childless /= post_mass

    old = np.arange(P.old_retention_index, P.J)
    old_mass = float(age_mass[old].sum()) if len(old) else 0.0
    old_owner_rate = float(dist[old, :, :, :, 1:].sum() / old_mass) if old_mass > 0 else np.nan

    young = np.arange(max(0, P.fertility_choice_index - 5), min(P.J, P.fertility_choice_index + 10))
    young_mass = float(age_mass[young].sum()) if len(young) else 0.0
    young_owner_rate = float(dist[young, :, :, :, 1:].sum() / young_mass) if young_mass > 0 else np.nan

    return SimpleNamespace(
        total_mass=total,
        mean_completed_fertility=float(completed),
        childless_rate=float(childless),
        own_rate=owner_mass / total if total > 0 else np.nan,
        young_owner_rate=young_owner_rate,
        old_owner_rate=old_owner_rate,
        owner_mass_by_size=owner_mass_by_size,
        owner_demand_by_size=owner_demand_by_size,
        owner_supply=np.asarray(P.owner_supply, dtype=float).copy(),
        owner_excess_by_size=owner_demand_by_size - P.owner_supply,
        own_by_age=own_by_age,
        children_by_age=children_by_age,
        q_owner=np.asarray(q_owner, dtype=float).copy(),
        p_asset=np.asarray(p_asset, dtype=float).copy(),
    )
