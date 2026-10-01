"""Shared stage of the extracted solver (bodies byte-identical; see split_receipt.json)."""
from __future__ import annotations
import copy
import math
import time
from types import SimpleNamespace
from typing import Any
import numpy as np
from . import joint_nested
from .adult_entry import adjusted_births, potential_entry_households
from .warm_price import search_warm_price
from .child_preferences import apply_child_preferences
from .parameters import (bequest_utility_net_active, child_earnings_multiplier, child_earnings_penalty_active, children_at_home_count, estate_housing_value, estate_receiver_active, estate_transfer_at_age, get_fecundity_by_age, independent_child_maturation_active, mortgage_stay_floor_active, parent_age_maturation_active, readiness_childless_states, readiness_cumulative_probability, readiness_gate_active, readiness_settled_state, readiness_transition_hazard, rental_wedge_active, unsecured_debt_floor)
from .kernels import (NUMBA_AVAILABLE, full_owner_block_kernel, full_renter_block_kernel, location_logit_kernel, scatter_cols_kernel, scatter_cols_sameidx_kernel, scatter_vec_kernel, tenure_choice_kernel, tenure_logit_kernel)
from .utils import (decode_flat_family_state, flat_nc, interp_indices, interp_on_grid, interp_vector, logsumexp, make_grid, make_value_interp, scatter_redistribute, scatter_redistribute_cols, scatter_redistribute_cols_sameidx, unflat_nc, weighted_median_from_cells, weighted_quantile)


OUTSIDE_OPTION_CLOSURES = {"outside_option", "outside_option_local_births", "local_births_outside", "open_city"}


ACCOUNTING_SCALE_PRICE_CLOSURES = {"accounting_scale_prices", "scaled_housing", "scaled_housing_accounting"}


BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES = {
    "outside_option_benchmark_normalized",
    "benchmark_outside_option_normalized",
}


RENEWAL_VALVE_FIXED_CLOSURES = {"renewal_valve", "renewal_scale", "demographic_valve"}


RENEWAL_VALVE_CALIBRATED_CLOSURES = {
    "renewal_valve_calibrated",
    "renewal_calibrated",
    "demographic_valve_calibrated",
    "benchmark_renewal_valve",
}


RENEWAL_VALVE_CLOSURES = RENEWAL_VALVE_FIXED_CLOSURES | RENEWAL_VALVE_CALIBRATED_CLOSURES


DEAD_VALUE_CUTOFF = -1e9


DEAD_MASS_TOL = 1e-12


class InfeasibleThetaError(RuntimeError):
    """Raised when positive population mass reaches a Bellman-dead state."""

    def __init__(
        self,
        stage: str,
        dead_mass: float,
        census: list[dict[str, Any]] | None = None,
    ) -> None:
        self.stage = str(stage)
        self.dead_mass = float(dead_mass)
        self.census = list(census or [])
        detail = jsonable_feasibility_census(self.census)
        super().__init__(
            f"{self.stage}: dead-node mass {self.dead_mass:.12g} exceeds {DEAD_MASS_TOL:.1e}; "
            f"census={detail}"
        )


def jsonable_feasibility_census(census: list[dict[str, Any]]) -> str:
    """Compact deterministic representation used in exception messages."""

    rows: list[str] = []
    for row in census[:8]:
        fields = ",".join(f"{key}={row[key]}" for key in sorted(row))
        rows.append("{" + fields + "}")
    return "[" + ";".join(rows) + "]"


def uses_outside_option_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in OUTSIDE_OPTION_CLOSURES


def normalize_population_mass(P: SimpleNamespace) -> bool:
    if uses_outside_option_closure(P):
        return False
    return bool(getattr(P, "normalize_population_mass", True))


def income_type_values(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    weights = np.asarray(getattr(P, "z_weights", np.ones(len(z))), dtype=float).reshape(-1)
    if weights.size != z.size:
        weights = np.ones(z.size)
    weights = np.maximum(weights, 0.0)
    weights = weights / weights.sum() if weights.sum() > 0 else np.ones(z.size) / max(z.size, 1)
    return z, weights


def income_transition_values(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z, weights = income_type_values(P)
    Pi = np.asarray(getattr(P, "Pi_z", np.eye(len(z))), dtype=float)
    if Pi.shape != (len(z), len(z)):
        Pi = np.eye(len(z))
    Pi = np.maximum(Pi, 0.0)
    row_sum = Pi.sum(axis=1)
    for row in range(Pi.shape[0]):
        if row_sum[row] > 0:
            Pi[row, :] /= row_sum[row]
        else:
            Pi[row, :] = weights
    return z, weights, Pi


def income_at_state(P: SimpleNamespace, i: int, j: int, z_value: float) -> float:
    y = float(P.income[i, j])
    if j < int(getattr(P, "J_R", P.J)):
        income = y * float(z_value)
    else:
        scale = float(getattr(P, "retirement_income_z_scale", 0.0))
        income = y * (1.0 + scale * (float(z_value) - 1.0))
    # The estate transfer enters exactly like the property-tax rebate: in the
    # budget through this function, never in the down-payment test (which
    # reads dp_arr/bmo, not income). Off, this returns the prior value bit
    # for bit.
    if not estate_receiver_active(P):
        return income + float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
    return (
        income
        + float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
        + estate_transfer_at_age(P, int(j))
    )


def penalized_income_at_state(P: SimpleNamespace, i: int, j: int, z_value: float, m: int) -> float:
    """Working-age after-tax earnings net of the children-at-home time cost.

    Multiplies the earnings part of ``income_at_state`` by
    ``(1 - penalty(m))`` at working ages; the lump-sum fiscal transfer is
    not earnings and is not scaled. Retirement ages are untouched. With the
    switch off this equals ``income_at_state`` bit for bit.
    """
    mult = child_earnings_multiplier(P, int(j), int(m))
    if mult == 1.0:
        return income_at_state(P, i, j, z_value)
    y = float(P.income[i, j])
    if j < int(getattr(P, "J_R", P.J)):
        income = y * float(z_value) * mult
    else:
        scale = float(getattr(P, "retirement_income_z_scale", 0.0))
        income = y * (1.0 + scale * (float(z_value) - 1.0))
    # Fiscal lump sums are not earnings and are not scaled by the time cost;
    # the estate transfer joins the rebate here on the same terms.
    if not estate_receiver_active(P):
        return income + float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
    return (
        income
        + float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
        + estate_transfer_at_age(P, int(j))
    )


def housing_demand_normalizer(P: SimpleNamespace) -> float:
    if normalize_population_mass(P):
        return max(float(getattr(P, "N_target", 1.0)), 1e-12)
    return 1.0


ENTRY_WEALTH_INCOME_RATIO_MODES = {
    "income_ratio",
    "income_ratio_distribution",
    "external_income_ratio",
    "external",
}


def annual_gross_income_at_state(P: SimpleNamespace, i: int, j: int, z_value: float, m: int = 0) -> float:
    """Annual gross-normalized income in the same units as PSID wealth ratios."""
    period_years = float(getattr(P, "period_years", getattr(P, "da", 1.0)))
    tau = float(getattr(P, "tau_pay", 0.0))
    period_income = penalized_income_at_state(P, i, j, z_value, int(m))
    annual_aftertax = period_income / max(period_years, 1e-12)
    if j < int(getattr(P, "J_R", P.J)):
        return annual_aftertax / max(1.0 - tau, 1e-12)
    return annual_aftertax


def _linear_grid_weights_for_points(
    b_grid: np.ndarray, points: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    x = np.asarray(points, dtype=float).reshape(-1)
    w = np.asarray(weights, dtype=float).reshape(-1)
    if x.size != w.size:
        raise ValueError("points and weights must have the same length")
    if bg.size == 0:
        raise ValueError("b_grid must contain at least one node")
    mass = np.zeros(bg.size, dtype=float)
    for val, wt in zip(x, w):
        if wt <= 0.0 or not np.isfinite(val):
            continue
        b = float(np.clip(val, bg[0], bg[-1]))
        hi = int(np.searchsorted(bg, b, side="left"))
        if hi <= 0:
            mass[0] += float(wt)
        elif hi >= bg.size:
            mass[-1] += float(wt)
        elif abs(float(bg[hi]) - b) <= 1e-14:
            mass[hi] += float(wt)
        else:
            lo = hi - 1
            share_hi = (b - float(bg[lo])) / max(float(bg[hi] - bg[lo]), 1e-12)
            mass[lo] += float(wt) * (1.0 - share_hi)
            mass[hi] += float(wt) * share_hi
    total = float(np.sum(mass))
    if total <= 0.0:
        idx = int(np.argmax(bg >= 0.0)) if np.any(bg >= 0.0) else int(bg.size - 1)
        return np.array([idx], dtype=np.int64), np.array([1.0], dtype=float)
    idx = np.flatnonzero(mass > 0.0).astype(np.int64)
    wt = mass[idx] / total
    return idx, wt


def configure_current_household_contract(P: SimpleNamespace) -> SimpleNamespace:
    """Enable the accepted September 27 entry and purchase accounting.

    Caller supplies the authenticated conditional entry arrays and explicit
    wealth grid. Preferences, income, population and fiscal objects are bound
    separately; this helper never fills them from historical defaults.
    """
    for name in ("fixed_reference_entry_conditional", "fixed_reference_entry_grid",
                 "earnings_transaction_grid"):
        if not hasattr(P, name):
            raise ValueError(f"current household contract missing {name}")
    P.native_fixed_reference_entry = True
    P.native_explicit_transaction_grid = True
    P.native_purchase_income = True
    P.native_exact_allocation_output = True
    return P


def precompute_shared(P: SimpleNamespace, b_grid: np.ndarray) -> SimpleNamespace:
    Nb = len(b_grid)
    nc = P.n_parity * P.n_child_states
    K = P.n_child_stages
    csm1 = K + 1
    g0 = float(getattr(P, "transfer_floor_G0", 0.0))
    gn = float(getattr(P, "transfer_floor_Gn", 0.0))
    c_bar = np.zeros((P.n_parity, P.n_child_states))
    h_bar = np.zeros((P.n_parity, P.n_child_states))
    psi_v = np.zeros((P.n_parity, P.n_child_states))
    g_bar = np.zeros((P.n_parity, P.n_child_states))
    alpha_bar = np.full((P.n_parity, P.n_child_states), P.alpha_cons)
    escale = np.ones((P.n_parity, P.n_child_states))
    eqscale_form = str(getattr(P, "eqscale_form", "linear")).lower()
    child_room_floor_active = bool(getattr(P, "child_room_floor", False)) and (
        float(getattr(P, "hbar_child_rooms", 0.0)) > 0.0
        or float(getattr(P, "hbar_first_child_jump", 0.0)) > 0.0
    )
    if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale" and eqscale_form not in {
        "linear", "power", "sqrt"
    }:
        raise ValueError("eqscale_form must be one of: linear, power, sqrt")
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            if independent_child_maturation_active(P):
                nk = cs if cs <= nn else 0
                kp = nk > 0
            else:
                nk = nn
                kp = (cs >= 1) and (cs < csm1)
            if readiness_gate_active(P) and nn == 0 and cs == 1:
                # E6c reuses the otherwise invalid (childless, cs=1) cell for
                # the settled state. It has childless preferences, not the
                # child-at-home consumption and housing adjustments.
                kp = False
            if kp:
                c_bar[nn, cs] = P.c_bar_0 + P.c_bar_n * nk
                if str(P.child_housing_spec).lower() == "linear_only":
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_n * nk
                else:
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * nk
                psi_v[nn, cs] = P.psi_child * nk
                g_bar[nn, cs] = g0 + gn * nk
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    if child_room_floor_active:
                        h_bar[nn, cs] = (
                            float(P.hbar_first_child_jump)
                            + float(P.hbar_child_rooms) * nk
                        )
                        alpha_bar[nn, cs] = P.alpha_cons
                    else:
                        h_bar[nn, cs] = 0.0
                        alpha_bar[nn, cs] = np.clip(
                            P.alpha_cons - (P.delta_alpha_jump + P.delta_alpha * nk), 0.05, 0.95
                        )
                    if eqscale_form == "power":
                        # Imposed Scholz-Seshadri-Khitatrakun (2006, JPE 114(4), p.619;
                        # Citro-Michael 1995) scale relative to a childless couple:
                        #   e(n) = ((2 + 0.7 n)/2)**0.7.
                        # Flow utility is multiplied by escale, u = escale * x**(1-sigma)/(1-sigma),
                        # while per-equivalent CRRA utility is u(x/e) = e**(sigma-1) * x**(1-sigma)/(1-sigma),
                        # so the multiplier is e**(sigma-1); at the baseline sigma = 2 the
                        # multiplier equals the scale itself. n is the literal parity state
                        # (under L4, nn = 3 is the top-coded 3+ bin, scaled at n = 3).
                        escale[nn, cs] = (((2.0 + 0.7 * nk) / 2.0) ** 0.7) ** (float(P.sigma) - 1.0)
                    elif eqscale_form == "sqrt":
                        # Declared robustness alternative: square-root household-size scale,
                        # e(n) = sqrt((2 + n)/2) relative to a childless couple.
                        escale[nn, cs] = (((2.0 + nk) / 2.0) ** 0.5) ** (float(P.sigma) - 1.0)
                    else:
                        escale[nn, cs] = 1.0 + P.gamma_e * nk
            else:
                c_bar[nn, cs] = P.c_bar_0
                h_bar[nn, cs] = P.h_bar_0
                g_bar[nn, cs] = g0
                if str(getattr(P, "preference_spec", "stone_geary")).lower() == "eqscale":
                    c_bar[nn, cs] = 0.0
                    h_bar[nn, cs] = 0.0

    apply_child_preferences(P, alpha_bar, psi_v, escale)
    ces_ec=np.ones_like(escale)
    ces_eh=np.ones_like(escale)
    if bool(getattr(P,"ces_enabled",False)):
        if not independent_child_maturation_active(P):
            raise ValueError("CES requires independent current-child counts")
        if (np.any(c_bar != 0.0) or np.any(h_bar != 0.0)
                or bool(getattr(P,"compensated_child_housing_shares",False))
                or float(P.delta_alpha)!=0.0 or float(P.delta_alpha_jump)!=0.0
                or float(P.sigma)!=2.0 or float(P.alpha_cons)!=.733):
            raise ValueError("CES primitive contract drift: no floors/A/share variation, sigma2, alpha.733")
        escale[:]=1.0
        alpha_bar[:]=P.alpha_cons
        for nn in range(P.n_parity):
            for cs in range(P.n_child_states):
                m=cs if cs<=nn else 0
                ces_ec[nn,cs]=((2.0+.7*m)/2.0)**.7
                ces_eh[nn,cs]=ces_ec[nn,cs]*(1.0+float(P.lambda_housing)*(m>0))
    triples = np.column_stack(
        [
            c_bar.reshape(-1, order="F"),
            h_bar.reshape(-1, order="F"),
            psi_v.reshape(-1, order="F"),
        ]
    )
    unique_triples, type_map = np.unique(triples, axis=0, return_inverse=True)
    birth_dp = np.zeros((P.n_parity, P.n_child_states, 1 + P.n_house, 1 + P.n_house), dtype=bool)
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            for to in range(1 + P.n_house):
                for tn in range(1 + P.n_house):
                    birth_dp[nn, cs, to, tn] = has_birth_dp_grant(P, nn, cs, to, tn)

    return SimpleNamespace(
        c_bar=c_bar,
        h_bar=h_bar,
        psi_v=psi_v,
        g_bar=g_bar,
        cb_flat=c_bar.reshape(1, nc, order="F"),
        hb_flat=h_bar.reshape(1, nc, order="F"),
        psi_flat=psi_v.reshape(1, nc, order="F"),
        gb_flat=g_bar.reshape(1, nc, order="F"),
        alpha_flat=alpha_bar.reshape(1, nc, order="F"),
        escale_flat=escale.reshape(1, nc, order="F"),
        ces_ec_flat=ces_ec.reshape(1,nc,order="F"),
        ces_eh_flat=ces_eh.reshape(1,nc,order="F"),
        nc=nc,
        b=b_grid.reshape(-1, 1),
        bp=b_grid.reshape(1, -1),
        phi_state=get_phi_state_matrix(P),
        phi_choice=get_phi_choice_tensor(P),
        n_types=unique_triples.shape[0],
        type_map=type_map,
        type_cb=unique_triples[:, 0],
        type_hb=unique_triples[:, 1],
        type_psi=unique_triples[:, 2],
        birth_dp=birth_dp,
        birth_entry_grant=get_birth_entry_grant_tensor(P),
    )


def add_aggregate_wealth_gross_labor_diagnostics(
    stats: SimpleNamespace,
    g: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
) -> None:
    """Add the matched aggregate gross/gross ratio and lifecycle diagnostics.

    The numerator is beginning-of-period net worth for every living household.
    The denominator is annual gross labor earnings for working households only.
    Age-binned ratios are robustness diagnostics, never calibration targets.
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

    wealth_by_age = np.zeros(int(P.J), dtype=float)
    gross_labor_earnings_by_age = np.zeros(int(P.J), dtype=float)
    for j in range(int(P.J)):
        for i in range(int(P.I)):
            for zz in range(g7.shape[4]):
                z_value = float(z_values[zz]) if zz < z_values.size else 1.0
                state_mass = float(np.sum(g7[:, :, i, j, zz, :, :]))
                if j < int(P.J_R):
                    if child_earnings_penalty_active(P):
                        for nn in range(int(P.n_parity)):
                            for cs in range(int(P.n_child_states)):
                                cell_mass = float(np.sum(g7[:, :, i, j, zz, nn, cs]))
                                gross_labor_earnings_by_age[j] += (
                                    annual_gross_income_at_state(
                                        P, i, j, z_value, children_at_home_count(nn, cs, P)
                                    )
                                    * cell_mass
                                )
                    else:
                        gross_labor_earnings_by_age[j] += (
                            annual_gross_income_at_state(P, i, j, z_value) * state_mass
                        )
                for ten in range(g7.shape[1]):
                    housing_value = (
                        float(ph_arr[i]) * float(P.H_own[ten - 1])
                        if ten > 0
                        else 0.0
                    )
                    mass_by_asset = np.sum(
                        g7[:, ten, i, j, zz, :, :],
                        axis=(1, 2),
                    )
                    wealth_by_age[j] += float(
                        np.sum(mass_by_asset * (bg_arr + housing_value))
                    )

    aggregate_wealth = float(np.sum(wealth_by_age))
    aggregate_gross_labor_earnings = float(
        np.sum(gross_labor_earnings_by_age)
    )
    stats.aggregate_wealth = aggregate_wealth
    stats.aggregate_annual_gross_labor_earnings = (
        aggregate_gross_labor_earnings
    )
    stats.aggregate_wealth_to_annual_gross_labor_earnings = (
        aggregate_wealth / max(aggregate_gross_labor_earnings, 1e-12)
    )
    stats.aggregate_wealth_by_age = wealth_by_age
    stats.aggregate_annual_gross_labor_earnings_by_age = (
        gross_labor_earnings_by_age
    )
    model_ages = float(P.age_start) + np.arange(int(P.J)) * float(P.da)
    state_years = float(P.da)
    for age_lo, age_hi in ((26, 35), (36, 45), (46, 55), (56, 65)):
        overlap_years = np.maximum(
            0.0,
            np.minimum(model_ages + state_years, float(age_hi + 1))
            - np.maximum(model_ages, float(age_lo)),
        )
        bin_share = overlap_years / state_years
        setattr(
            stats,
            f"aggregate_wealth_to_annual_gross_labor_earnings_{age_lo}_{age_hi}",
            float(np.sum(wealth_by_age * bin_share))
            / max(
                float(np.sum(gross_labor_earnings_by_age * bin_share)),
                1e-12,
            ),
        )


def get_completed_fertility(nn: int, cs: int, P: SimpleNamespace) -> int:
    """Completed births for family state (nn, cs).

    Parity ``nn`` is preserved through maturation (the child transition is
    per-``nn`` over ``cs`` only), so parity-3+ matured states read completed
    fertility from ``nn``.  For nn <= 2 the legacy cs-based mapping is kept
    verbatim — including on unreachable cells such as (nn=2, cs=K+1) — so
    the precomputed bequest table, and hence V, nest bit for bit at
    n_parity=3.
    """
    if independent_child_maturation_active(P):
        return int(nn)
    K = P.n_child_stages
    if cs == 0:
        return 0
    if cs == K + 1:
        return 1 if nn <= 2 else nn
    if cs == K + 2:
        return 2 if nn <= 2 else nn
    return nn


def has_birth_dp_grant(P, nn, cs, to, tn):
    if not bool(getattr(P, "birth_dp_grant", False)):
        return False
    if to != 0 or tn <= 0:
        return False
    if independent_child_maturation_active(P):
        return (nn >= 1) and (1 <= cs <= nn)
    return (nn >= 1) and (cs == 1)


def get_phi_state_matrix(P):
    ps = np.tile(P.phi.reshape(-1, 1), (1, P.n_child_states))
    if not bool(getattr(P, "parent_dp_waiver", False)):
        return ps
    kp = get_parent_target_child_states(P)
    po = getattr(P, "parent_dp_waiver_phi", 1.0)
    if P.n_parity >= 2:
        ps[1:, kp] = np.maximum(ps[1:, kp], po)
    return ps


def get_phi_choice_tensor(P):
    nt = 1 + P.n_house
    pc = np.tile(P.phi.reshape(1, 1, P.n_parity, 1), (P.I, nt, 1, P.n_child_states))
    if not bool(getattr(P, "parent_dp_waiver", False)) or P.n_parity < 2:
        return pc
    po = getattr(P, "parent_dp_waiver_phi", 1.0)
    loc_idx = get_parent_target_locations(P)
    ten_idx = get_parent_target_owner_tenures(P)
    cs_idx = get_parent_target_child_states(P)
    if len(loc_idx) == 0 or len(ten_idx) == 0 or not np.any(cs_idx):
        return pc
    pc[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.where(cs_idx)[0])] = np.maximum(
        pc[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.where(cs_idx)[0])], po
    )
    return pc


def get_parent_target_locations(P):
    loc_idx = np.arange(P.I)
    vals = getattr(P, "parent_dp_waiver_locations", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        loc_idx = np.asarray(vals, dtype=int) - 1
    return loc_idx[(loc_idx >= 0) & (loc_idx < P.I)]


def get_parent_target_owner_tenures(P):
    ten_idx = np.arange(1, 1 + P.n_house)
    vals = getattr(P, "parent_dp_waiver_owner_rungs", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        ten_idx = np.asarray(vals, dtype=int)
    return ten_idx[(ten_idx >= 1) & (ten_idx <= P.n_house)]


def get_parent_target_child_states(P):
    kp = np.zeros(P.n_child_states, dtype=bool)
    if independent_child_maturation_active(P):
        kp[1:] = True
        return kp
    if bool(getattr(P, "parent_dp_waiver_birth_state_only", False)):
        if P.n_child_states >= 2:
            kp[1] = True
    else:
        kp[1 : P.n_child_stages + 1] = True
    return kp


def get_birth_entry_grant_tensor(P):
    nt = 1 + P.n_house
    bg = np.zeros((P.I, nt, P.n_parity, P.n_child_states))
    if not bool(getattr(P, "birth_entry_grant", False)):
        return bg
    grant = getattr(P, "birth_entry_grant_amount", 0.0)
    if not (np.isscalar(grant) and np.isfinite(grant) and grant > 0) or P.n_parity < 2 or P.n_child_states < 2:
        return bg
    loc_idx = np.arange(P.I)
    vals = getattr(P, "birth_entry_grant_locations", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        loc_idx = np.asarray(vals, dtype=int) - 1
    loc_idx = loc_idx[(loc_idx >= 0) & (loc_idx < P.I)]
    ten_idx = np.arange(1, nt)
    vals = getattr(P, "birth_entry_grant_owner_rungs", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        ten_idx = np.asarray(vals, dtype=int)
    ten_idx = ten_idx[(ten_idx >= 1) & (ten_idx < nt)]
    if len(loc_idx) == 0 or len(ten_idx) == 0:
        return bg
    if independent_child_maturation_active(P):
        for parity in range(1, P.n_parity):
            child_states = np.arange(1, min(parity, P.n_child_states - 1) + 1)
            bg[np.ix_(loc_idx, ten_idx, np.array([parity]), child_states)] = grant
    else:
        bg[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.array([1]))] = grant
    return bg


def property_tax_revenue_from_distribution(g, hR_pol, p_hat, P):
    """Stationary property-tax revenue from all occupied housing."""
    tax_rate = max(float(getattr(P, "tau_H", 0.0)), 0.0)
    prices = np.asarray(p_hat, dtype=float).reshape(-1)
    revenue = 0.0
    for market in range(P.I):
        rental_services = float(
            np.sum(g[:, 0, market, :, :, :, :] * hR_pol[:, 0, market, :, :, :, :])
        )
        owner_services = 0.0
        for tenure in range(1, 1 + P.n_house):
            owner_services += float(np.sum(g[:, tenure, market, :, :, :, :])) * float(
                P.H_own[tenure - 1]
            )
        revenue += tax_rate * float(prices[market]) * (rental_services + owner_services)
    return float(revenue)


def markov_grant_outlays(g, tenure_choice, tenure_probs, P, SD):
    """Count actual renter-to-owner grant payments in eligible parent states."""
    grant = np.asarray(SD.birth_entry_grant, dtype=float)
    if not np.any(grant > 0.0):
        return 0.0, 0.0
    if P.I != 1:
        raise NotImplementedError("funded-grant accounting is currently restricted to the one-market model")
    recipient_mass = 0.0
    outlays = 0.0
    for age in range(P.J):
        for income_state in range(g.shape[4]):
            for parity in range(1, P.n_parity):
                child_states = (
                    range(1, min(parity, P.n_child_states - 1) + 1)
                    if independent_child_maturation_active(P)
                    else (1,)
                )
                for child_state in child_states:
                    source = np.asarray(
                        g[:, 0, 0, age, income_state, parity, child_state], dtype=float
                    )
                    if float(np.sum(source)) <= 1e-15:
                        continue
                    eligible = np.where(
                        (grant[0, :, parity, child_state] > 0.0)
                        & (~np.asarray(SD.birth_dp[parity, child_state, 0, :], dtype=bool))
                    )[0]
                    if eligible.size == 0:
                        continue
                    if tenure_probs is None:
                        chosen = np.asarray(
                            tenure_choice[:, 0, 0, age, income_state, parity, child_state],
                            dtype=int,
                        )
                        for tenure in eligible:
                            paid_mass = float(np.sum(source[chosen == tenure]))
                            recipient_mass += paid_mass
                            outlays += paid_mass * float(grant[0, tenure, parity, child_state])
                    else:
                        probabilities = np.asarray(
                            tenure_probs[:, 0, 0, age, income_state, parity, child_state, :],
                            dtype=float,
                        )
                        probability_sum = np.sum(probabilities, axis=1)
                        normalized = np.divide(
                            probabilities,
                            probability_sum[:, None],
                            out=np.zeros_like(probabilities),
                            where=probability_sum[:, None] > 0.0,
                        )
                        eligible_probability = np.sum(normalized[:, eligible], axis=1)
                        recipient_mass += float(np.sum(source * eligible_probability))
                        weighted_grant = normalized[:, eligible] @ grant[
                            0, eligible, parity, child_state
                        ]
                        outlays += float(np.sum(source * weighted_grant))
    return float(recipient_mass), float(outlays)
