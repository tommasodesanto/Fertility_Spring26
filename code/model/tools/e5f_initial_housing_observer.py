"""Default-off initial housing/wealth diagnostics; never a target contract.

The caller supplies a pre-announcement stationary PeriodEvaluation and explicitly
chooses uniform age exposure within each four-year cell. This module does not
solve, infer stationarity, alter targets, or certify the empirical approximations.
"""
from __future__ import annotations

import math
from types import SimpleNamespace
from typing import Any

import numpy as np


AGE_PROJECTION = "uniform_within_age_cell"
MOMENT_NAMES = (
    "aggregate_mean_occupied_rooms_capped9_18_85",
    "own_rate_30_55",
    "own_rate_25_34",
    "prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9",
    "recent_parent_minus_no_resident_child_ownership_30_55",
    "aggregate_wealth_to_annual_gross_labor_earnings",
    "annual_bequest_flow_to_aggregate_wealth",
    "old_total_wealth_to_annual_income_p90_p50_7684",
    "old_total_wealth_to_annual_income_median_7684",
    "housing_increment_0to1",
)


def uniform_age_cell_overlap(parameters, lower: float, upper: float) -> np.ndarray:
    """Fractions of model age intervals overlapping [lower, upper).

    For data ages 76--84 the interval is [76,85), giving 1/2 of the 74 cell,
    all of 78, and 3/4 of 82 on the maintained 18+4*j grid. This is an explicit
    within-cell homogeneity approximation, not reconstructed annual-age states.
    """
    start, width = float(parameters.age_start), float(parameters.da)
    count = int(parameters.J)
    if (not all(math.isfinite(x) for x in (start, width, lower, upper))
            or width <= 0 or count <= 0 or lower >= upper):
        raise ValueError("Age intervals and cell geometry must be finite and ordered")
    left = start + width * np.arange(count, dtype=float)
    return np.maximum(0.0, np.minimum(left + width, upper) - np.maximum(left, lower)) / width


def _row(name: str, reason: str) -> dict[str, Any]:
    return dict(moment=name, model_value=None, available=False,
                status="unavailable", reason=reason, production_eligible=False,
                approximations=[])


def _set_value(row, value, *, approximations=(), **details):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite diagnostic value: {row['moment']}")
    row.update(model_value=value, available=True, status="diagnostic_only",
               reason=None, approximations=list(approximations), **details)


def _positive_ratio(row, numerator, denominator, *, approximations=(), **details):
    numerator, denominator = float(numerator), float(denominator)
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError(f"Nonfinite numerator/denominator: {row['moment']}")
    row.update(numerator=numerator, denominator=denominator, **details)
    if denominator <= 0:
        row["reason"] = "Selected denominator is nonpositive; no denominator floor or invented zero"
        return
    _set_value(row, numerator / denominator, approximations=approximations)


def _validate_current(evaluation, parameters):
    g = np.asarray(evaluation.g_current, dtype=float)
    h_r = np.asarray(evaluation.policy.hR_pol, dtype=float)
    expected = (1 + int(parameters.n_house), int(parameters.I), int(parameters.J),
                len(parameters.z_grid), int(parameters.n_parity), int(parameters.n_child_states))
    if g.ndim != 7 or g.shape[1:] != expected or h_r.shape != g.shape:
        raise ValueError("Current mass and renter policies require matching full income-resolved 7D states")
    if not np.isfinite(g).all() or np.any(g < 0):
        raise ValueError("Current household mass must be finite and nonnegative")
    z_values = np.asarray(parameters.z_grid, dtype=float)
    if z_values.ndim != 1 or not np.isfinite(z_values).all() or np.any(z_values < 0):
        raise ValueError("Income states must be a finite nonnegative vector")
    houses = np.asarray(parameters.H_own, dtype=float)
    if houses.shape != (int(parameters.n_house),) or not np.isfinite(houses).all() or np.any(houses <= 0):
        raise ValueError("Owner products must be finite positive room quantities")
    occupied = g[:, 0] > 0
    if np.any(occupied & (~np.isfinite(h_r[:, 0]) | (h_r[:, 0] <= 0))):
        raise ValueError("Occupied renters require finite positive realized rooms")
    return g, h_r, houses


def _housing_totals(g, h_r, houses, age_weights, family_mask=None):
    # Sum after capping the actual income-state policy. Never cap a policy
    # collapsed over income; unoccupied policies may legitimately be NaN.
    if family_mask is None:
        family_mask = np.ones(g.shape[-2:], dtype=bool)
    mass = rooms = owner_mass = 0.0
    for j, fraction in enumerate(age_weights):
        if fraction <= 0:
            continue
        renters = g[:, 0, :, j, :, :, :] * family_mask
        rented_rooms = np.where(renters > 0,
                               np.minimum(h_r[:, 0, :, j, :, :, :], 9.0), 0.0)
        mass += float(fraction * np.sum(renters))
        rooms += float(fraction * np.sum(renters * rented_rooms))
        for ten, house in enumerate(houses, start=1):
            owners = float(fraction * np.sum(g[:, ten, :, j, :, :, :] * family_mask))
            mass += owners
            owner_mass += owners
            rooms += owners * min(float(house), 9.0)
    return dict(mass=mass, capped_rooms_sum=rooms, owner_mass=owner_mass)


def _wealth_diagnostics(rows, evaluation, parameters, b_grid, age_weights):
    """Reuse the existing balance-sheet primitive, not a second asset formula."""
    from intergen_eqscale_seq_optimized import solver as model
    from intergen_eqscale_seq_optimized.utils import weighted_quantile

    if float(parameters.period_years) != 4.0 or not bool(getattr(parameters, "scale_flows_to_period", False)):
        raise ValueError("Wealth observer requires explicit four-year income-flow units")
    asset_g = np.asarray(evaluation.g_post_fertility, dtype=float)
    death_g = np.asarray(evaluation.g_current, dtype=float)
    bp = np.asarray(evaluation.policy.bp_pol, dtype=float)
    bg = np.asarray(b_grid, dtype=float)
    prices = np.asarray(evaluation.policy.price, dtype=float)
    if asset_g.shape != death_g.shape or bp.shape != death_g.shape:
        raise ValueError("Beginning balance sheets, post-choice deaths and saving require identical 7D shapes")
    if (not np.isfinite(asset_g).all() or np.any(asset_g < 0)
            or bg.shape != (asset_g.shape[0],) or not np.isfinite(bg).all()
            or prices.shape != (int(parameters.I),) or not np.isfinite(prices).all()
            or np.any(prices <= 0)):
        raise ValueError("Invalid asset distribution, financial-wealth grid or prices")
    if not math.isclose(float(np.sum(asset_g)), float(np.sum(death_g)),
                        rel_tol=2e-10, abs_tol=2e-10):
        raise ValueError("Beginning and post-choice distributions must have the same living mass")
    if not np.isfinite(bp).all():
        raise ValueError("Existing at-death primitive requires finite saving policies")
    income = np.asarray(parameters.income, dtype=float)
    tau = float(parameters.tau_pay)
    if (income.shape != (int(parameters.I), int(parameters.J))
            or not np.isfinite(income).all() or np.any(income < 0)
            or not 0 <= tau < 1):
        raise ValueError("Invalid period income or payroll wedge")
    if int(parameters.J_R) != 12:
        raise ValueError("Gross labor earnings require the 12 working cells covering ages 18--65")
    if bool(getattr(parameters, "use_age_survival", False)):
        survival = np.asarray(parameters.survival_probs, dtype=float)
        if (survival.ndim != 1 or len(survival) < int(parameters.J) - 1
                or not np.isfinite(survival).all() or np.any((survival < 0) | (survival > 1))):
            raise ValueError("Death-flow accounting requires valid survival probabilities")
    # This primitive's stock denominator uses ALL ages, so refuse unsupported
    # aggregate age geometry instead of silently observing outside 18--85.
    aggregate_full = np.all(uniform_age_cell_overlap(parameters, 18.0, 86.0) == 1.0)
    aggregate_complete = (float(parameters.age_start) == 18.0
                          and float(parameters.da) * int(parameters.J) == 68.0)
    stats = SimpleNamespace()
    model.add_aggregate_wealth_bequest_flow_moments(
        stats, asset_g, death_g, bp, parameters, bg, prices)
    primitive_note = (
        "Beginning-period net worth uses g_post_fertility before housing transactions; "
        "death flow uses g_current and post-saving policies. Annualizes exactly once."
    )
    if aggregate_full and aggregate_complete:
        _positive_ratio(rows["aggregate_wealth_to_annual_gross_labor_earnings"],
                        stats.aggregate_wealth, stats.aggregate_annual_gross_labor_earnings,
                        approximations=("National PSID household-income normalization; survey-wave stock/flow dates inherited",),
                        accounting=primitive_note)
        _positive_ratio(rows["annual_bequest_flow_to_aggregate_wealth"],
                        stats.annual_bequest_flow, stats.aggregate_wealth,
                        approximations=("External historical restriction; synthetic uncertainty; forced terminal deaths",),
                        accounting=primitive_note)
    else:
        for name in ("aggregate_wealth_to_annual_gross_labor_earnings", "annual_bequest_flow_to_aggregate_wealth"):
            rows[name]["reason"] = "Existing aggregate primitive needs exactly the full 18--85 age coverage"

    # Reuse the model's living-household wealth/income definition and quantile
    # convention, but use the explicit empirical-age overlap [76,85).
    values, weights = [], []
    selected_mass = excluded_income_mass = 0.0
    for j, fraction in enumerate(age_weights):
        if fraction <= 0:
            continue
        for i, price in enumerate(prices):
            for z, z_value in enumerate(parameters.z_grid):
                annual_income = model.annual_gross_income_at_state(parameters, i, j, float(z_value))
                for ten in range(asset_g.shape[1]):
                    w = fraction * np.sum(asset_g[:, ten, i, j, z, :, :], axis=(1, 2))
                    selected_mass += float(np.sum(w))
                    if not math.isfinite(annual_income) or annual_income <= 0:
                        excluded_income_mass += float(np.sum(w))
                        continue
                    housing_value = float(price * parameters.H_own[ten - 1]) if ten > 0 else 0.0
                    positive = w > 0
                    values.append((bg[positive] + housing_value) / annual_income)
                    weights.append(w[positive])
    old_names = ("old_total_wealth_to_annual_income_p90_p50_7684",
                 "old_total_wealth_to_annual_income_median_7684")
    proxy_notes = (
        "Uniform within-age-cell overlap for ages 76--84",
        "Modeled pension plus configured lump-sum income proxy; not PSID total family income",
        "Empirical $1000 real-income and observed-child-history filters are not reproduced",
    )
    details = dict(selected_mass=selected_mass, excluded_nonpositive_income_mass=excluded_income_mass,
                   pension_income_proxy=True,
                   property_tax_lump_sum_transfer=float(getattr(parameters, "property_tax_lump_sum_transfer", 0.0)))
    if selected_mass <= 0 or excluded_income_mass > 0:
        for name in old_names:
            rows[name].update(reason="Old sample is empty or contains occupied nonpositive-income states; no silent exclusion", **details)
        return
    v, w = np.concatenate(values), np.concatenate(weights)
    p50, p90 = float(weighted_quantile(v, w, .5)), float(weighted_quantile(v, w, .9))
    _set_value(rows[old_names[1]], p50, approximations=proxy_notes, **details)
    _positive_ratio(rows[old_names[0]], p90, p50, approximations=proxy_notes, **details)


def _stationary_birth_diagnostic(evaluation, parameters, b_grid, shared):
    import run_e5f_transition_calibration as measurement
    branch = measurement.begin_dated_first_birth_housing_branch(
        evaluation, parameters, b_grid, shared, origin_period=0)
    result = measurement.finish_dated_first_birth_housing_branch(
        branch, evaluation, parameters, b_grid, shared, destination_period=1)
    if result.get("census_age_bridge_applied") is not False:
        raise RuntimeError("Stationary matched birth branches must exclude the Census age bridge")
    for key in ("housing_response", "treated_mean_housing", "control_mean_housing",
                "treated_continuation_births", "origin_mass", "destination_mass"):
        if key not in result or not math.isfinite(float(result[key])):
            raise RuntimeError(f"Incomplete stationary matched-birth diagnostic: {key}")
    if float(result["origin_mass"]) <= 0 or float(result["destination_mass"]) <= 0:
        raise RuntimeError("Matched-birth diagnostics require positive branch mass")
    if float(result["treated_continuation_births"]) < 0:
        raise RuntimeError("Matched-birth continuation mass cannot be negative")
    if not math.isclose(float(result["housing_response"]),
                        float(result["treated_mean_housing"]) - float(result["control_mean_housing"]),
                        rel_tol=0, abs_tol=1e-10):
        raise RuntimeError("Matched-birth response must equal treated minus control rooms")
    return result


def observe_initial_housing_wealth(
    evaluation, parameters, b_grid, shared, *, diagnostic_enabled: bool = False,
    age_projection: str | None = None, diagnostic_allow_family_proxies: bool = False,
    include_wealth: bool = True, include_birth_response: bool = False,
) -> dict[str, Any]:
    """Observe supplied stationary choices with visible measurement limits.

    Disabled by default, with no model imports. Enabling requires the exact
    AGE_PROJECTION string. Family proxy and matched birth branching are separate
    opt-ins. The caller must establish that the supplied evaluation is stationary;
    this diagnostic neither solves nor certifies equilibrium/fiscal conditions.
    """
    for flag in (diagnostic_enabled, diagnostic_allow_family_proxies, include_wealth, include_birth_response):
        if type(flag) is not bool:
            raise TypeError("Diagnostic switches must be explicit booleans")
    rows = {name: _row(name, "Diagnostic observer is disabled") for name in MOMENT_NAMES}
    result = dict(observer_id="e5f_initial_housing_wealth_diagnostic_v1",
                  status="disabled", production_eligible=False, target_contract_activated=False,
                  stationary_input_certified=False, rows=list(rows.values()),
                  moments={name: None for name in MOMENT_NAMES}, age_projection=None)
    if not diagnostic_enabled:
        return result
    if age_projection != AGE_PROJECTION:
        raise ValueError("Explicit age_projection='uniform_within_age_cell' is required")
    if (float(parameters.period_years) != 4.0 or float(parameters.da) != 4.0
            or float(parameters.age_start) != 18.0 or int(parameters.J) != 17):
        raise ValueError("Initial observer requires 17 four-year cells spanning ages 18--85")
    g, h_r, houses = _validate_current(evaluation, parameters)
    ages = {key: uniform_age_cell_overlap(parameters, lo, hi) for key, lo, hi in (
        ("18_85", 18.0, 86.0), ("30_55", 30.0, 56.0),
        ("25_34", 25.0, 35.0), ("76_84", 76.0, 85.0))}
    for row in rows.values():
        row["reason"] = "Not implemented for this empirical definition"
    common = ("Uniform annual-age exposure and constant policies/distribution within each four-year cell",
              "National model versus 42-MET2013 housing sample")
    all_housing = _housing_totals(g, h_r, houses, ages["18_85"])
    _positive_ratio(rows[MOMENT_NAMES[0]], all_housing["capped_rooms_sum"], all_housing["mass"],
                    approximations=common, rooms_capped_at=9.0,
                    cap_before_income_aggregation=True)
    for label in ("30_55", "25_34"):
        totals = _housing_totals(g, h_r, houses, ages[label])
        _positive_ratio(rows[f"own_rate_{label}"], totals["owner_mass"], totals["mass"],
                        approximations=(*common, "ACS DUE structure restriction is not represented by a model state"))
    family = rows[MOMENT_NAMES[3]]
    if not diagnostic_allow_family_proxies:
        family["reason"] = "Requires diagnostic_allow_family_proxies=True; dependents are not exact ACS resident-child groups"
    elif str(getattr(parameters, "child_state_mode", "")) != "independent_count" or int(getattr(parameters, "child_bin_high_cutoff", -1)) != 3:
        family["reason"] = "Dependent-count proxy requires independent_count and child_bin_high_cutoff=3"
    else:
        n, m = np.indices(g.shape[-2:])
        # Never use lifetime parity to assign the current family-size bin.
        # Readiness states may exceed n for childless households; exclude them.
        low = _housing_totals(g, h_r, houses, ages["30_55"], (m >= 1) & (m <= 2) & (m <= n))
        high = _housing_totals(g, h_r, houses, ages["30_55"], (m >= 3) & (m <= n))
        family.update(low_mass=low["mass"], high_mass=high["mass"], model_dependent_proxy=True)
        if low["mass"] <= 0 or high["mass"] <= 0:
            family["reason"] = "One dependent-count group has zero mass; no denominator floor"
        else:
            _set_value(family, high["capped_rooms_sum"] / high["mass"] - low["capped_rooms_sum"] / low["mass"],
                       approximations=(*common, "Model dependent counts replace ACS resident own children with YNGCH<18"),
                       rooms_capped_at=9.0, cap_before_income_aggregation=True)
    rows[MOMENT_NAMES[4]]["reason"] = (
        "Unavailable: ACS oldest resident child under four and no-resident-child controls "
        "cannot be reconstructed from dependent counts and lifetime parity")
    if include_wealth:
        _wealth_diagnostics(rows, evaluation, parameters, b_grid, ages["76_84"])
    else:
        for name in MOMENT_NAMES[5:9]:
            rows[name]["reason"] = "Wealth diagnostics were not requested"
    if include_birth_response:
        birth = _stationary_birth_diagnostic(evaluation, parameters, b_grid, shared)
        _set_value(rows[MOMENT_NAMES[9]], birth["housing_response"],
                   approximations=("Pooled PSID response applied to initial stationary economy under stability assumption",
                                   "Matched model risk-set weights differ from Sun-Abraham cohort/event weights; non-flat prepath caveat retained"),
                   uncapped_rooms=True, destination_continuation_births_allowed=True,
                   stationary_policy_used_at_both_dates=True, branch=birth)
    else:
        rows[MOMENT_NAMES[9]]["reason"] = "Stationary matched-birth branching was not requested"
    result.update(status="diagnostic_only_measurement_approximations_unresolved",
                  age_projection=AGE_PROJECTION,
                  age_cell_labels=(float(parameters.age_start) + float(parameters.da) * np.arange(int(parameters.J))).tolist(),
                  age_overlap_weights={k: v.tolist() for k, v in ages.items()},
                  diagnostic_allow_family_proxies=diagnostic_allow_family_proxies,
                  moments={name: row["model_value"] for name, row in rows.items()})
    return result
