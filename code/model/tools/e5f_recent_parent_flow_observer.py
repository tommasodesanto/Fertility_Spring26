"""Default-off, passive selected-birth ownership diagnostic; no model solve.

This is a synchronized post-fertility snapshot, not an exact annual ACS
oldest-resident-child-under-four observer and not a causal birth response.
The caller must already have configured the sequential calendar model and
must supply a matching, price-gated PeriodEvaluation. Nothing is configured,
re-solved, aged, reweighted, or promoted to a production target here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from e5f_initial_housing_observer import AGE_PROJECTION, uniform_age_cell_overlap


SNAPSHOT = "synchronized_post_fertility_snapshot"
MOMENT = "ownership_current_birth_from_empty_dependent_home_minus_current_empty_home_30_55"
MASS_ATOL = 2e-10
PRUNING_TOLERANCE = 1e-15


def _array(value, name, shape, *, probability=False):
    result = np.asarray(value)
    if result.shape != shape or not np.isfinite(result).all() or np.any(result < 0):
        raise ValueError(f"{name} must have shape {shape} and finite nonnegative values")
    if probability and np.any(result > 1):
        raise ValueError(f"{name} contains a probability above one")
    return result


def _error(left, right):
    return float(np.max(np.abs(left - right)))


def _check(error, name):
    if not np.isfinite(error) or error > MASS_ATOL:
        raise ValueError(f"{name} failed: {error:.6g} exceeds {MASS_ATOL}")


def _group(g, age_weights, family_mask=None):
    selected = g if family_mask is None else g * family_mask
    mass_by_age = selected.sum(axis=(0, 1, 2, 4, 5, 6))
    owners_by_age = selected[:, 1:].sum(axis=(0, 1, 2, 4, 5, 6))
    denominator = float(age_weights @ mass_by_age)
    numerator = float(age_weights @ owners_by_age)
    if not np.isfinite([numerator, denominator]).all():
        raise ValueError("Nonfinite group numerator or denominator")
    return dict(owner_numerator=numerator, denominator=denominator,
                ownership_rate=numerator / denominator if denominator > 0 else None,
                mass_by_age=mass_by_age.tolist(), owner_mass_by_age=owners_by_age.tolist(),
                unweighted_mass=float(mass_by_age.sum()))


def observe_recent_parent_flow(
    evaluation: Any,
    parameters: Any,
    *,
    diagnostic_enabled: bool = False,
    snapshot: str | None = None,
    age_projection: str | None = None,
    diagnostic_allow_residence_proxy: bool = False,
    input_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Observe actual births into empty-dependent homes using fixed policies.

All three interpretation choices are explicit opt-ins. Enabled observation
raises on invalid inputs or nonpositive main-group denominators; it never
invents zeros or floors denominators. First/continuation subgroup rates may
be None if that subgroup is empty. Attribution uses destination parity after
the real one-birth kernel, which current housing transport preserves.
"""
    for flag in (diagnostic_enabled, diagnostic_allow_residence_proxy):
        if type(flag) is not bool:
            raise TypeError("Diagnostic opt-ins must be explicit booleans")
    packet = dict(observer_id="e5f_recent_parent_flow_diagnostic_v1", status="disabled",
                  moment=MOMENT, model_value=None, available=False,
                  production_eligible=False, target_contract_activated=False,
                  actual_weight=None, loss_contribution=None)
    if not diagnostic_enabled:
        return packet
    if snapshot != SNAPSHOT:
        raise ValueError(f"Explicit snapshot={SNAPSHOT!r} is required")
    if age_projection != AGE_PROJECTION:
        raise ValueError(f"Explicit age_projection={AGE_PROJECTION!r} is required")
    if not diagnostic_allow_residence_proxy:
        raise ValueError("Requires diagnostic_allow_residence_proxy=True")
    if input_provenance is not None:
        if not isinstance(input_provenance, dict):
            raise TypeError("input_provenance must be a dictionary or None")
        # Make a detached JSON receipt; reject NaNs or unserializable claims.
        input_provenance = json.loads(json.dumps(input_provenance, allow_nan=False))

    P = parameters
    geometry = (P.age_start, P.da, P.period_years, P.J, P.n_parity, P.n_child_states)
    if tuple(float(x) for x in geometry) != (18., 4., 4., 17., 4., 4.):
        raise ValueError("Requires 17 four-year age cells from 18, four parity/count states")
    if (not bool(getattr(P, "sequential_births", False))
            or getattr(P, "child_state_mode", None) != "independent_count"
            or getattr(P, "fertility_units", None) != "literal_topcode"
            or any(bool(getattr(P, key, False)) for key in
                   ("joint_nested_choice", "fertility_nest_choice", "two_shock_choice"))):
        raise ValueError("Unsupported architecture: requires maintained sequential independent_count")
    if (float(P.A_f_start), float(P.A_f_end)) != (1., 7.):
        raise ValueError("Requires maintained fertility cells beginning at ages 18 through 42")

    # Import only after opt-in. Never configure or mutate calendar.model/P.
    import run_e5f_open_population_transition as transition
    from intergen_eqscale_seq_optimized import solver as model

    if transition.calendar.model is not model:
        raise ValueError("Caller must configure the sequential calendar model before observation")
    policy = evaluation.policy
    if getattr(policy, "joint_choice", None) is not None:
        raise ValueError("Unsupported architecture: supplied policy owns a joint-choice object")
    pre = np.asarray(evaluation.g_pre)
    if pre.ndim != 7 or pre.shape[0] < 2 or pre.shape[4] < 1:
        raise ValueError("Mass requires seven axes, at least two wealth nodes and one income state")
    nb, nt, locations, ages, incomes, npar, ncs = pre.shape
    if (float(nt) != float(P.n_house) + 1 or float(locations) != float(P.I)
            or ages != 17 or npar != 4 or ncs != 4 or nt < 2 or locations < 1):
        raise ValueError("Mass dimensions do not match model geometry")
    pre = _array(pre, "g_pre", pre.shape)
    post = _array(evaluation.g_post_fertility, "g_post_fertility", pre.shape)
    current = _array(evaluation.g_current, "g_current", pre.shape)
    if float(pre.sum()) <= 0:
        raise ValueError("Input population mass must be positive")
    parity, child_state = np.indices((npar, ncs))
    childless_states = tuple(model.readiness_childless_states(P))
    never = (parity == 0) & np.isin(child_state, childless_states)
    former = (parity > 0) & (child_state == 0)
    empty = never | former
    valid = never | ((parity > 0) & (child_state <= parity))
    for name, mass in (("g_pre", pre), ("g_post_fertility", post), ("g_current", current)):
        if np.any(mass[..., ~valid] > 0):
            raise ValueError(f"{name} carries invalid family-state mass")

    fertility_shape = (nb, nt, locations, ages, incomes, npar)
    first = _array(policy.fert_probs, "fert_probs", fertility_shape, probability=True)
    continuation = transition.calendar.policy_continuation_birth_probs(policy, P)
    _array(continuation, "fert2_probs", (nb, nt, locations, ages, incomes, 2, 2, ncs),
           probability=True)
    _array(model.get_fecundity_by_age(P), "fecundity", (ages,), probability=True)
    location_probs = _array(
        policy.loc_probs, "loc_probs", (nb, nt, locations, locations, ages, incomes, npar, ncs),
        probability=True)
    # Occupied post-fertility origins require a complete location lottery.
    # Tenure lotteries, in contrast, are normalized by the existing transport.
    location_sum_error = float(np.max(np.abs(location_probs.sum(axis=3)[post > 0] - 1.),
                                      initial=0.))
    _check(location_sum_error, "occupied location probability sum")
    choices = _array(policy.tenure_choice, "tenure_choice", pre.shape)
    if np.any(choices != np.floor(choices)) or np.any(choices >= nt):
        raise ValueError("tenure_choice contains an invalid tenure index")
    if policy.tenure_probs is not None:
        _array(policy.tenure_probs, "tenure_probs", pre.shape + (nt,), probability=True)
    for prefix, shape in (("lmm", (locations, nt, nb)),
                          ("tmx", (locations, nt, nt, npar, ncs, nb))):
        indices = _array(getattr(policy.maps, prefix + "_idx"), prefix + "_idx", shape)
        if not np.issubdtype(indices.dtype, np.integer) or np.any(indices >= nb - 1):
            raise ValueError(f"{prefix}_idx must contain integer lower-bracket wealth indices")
        _array(getattr(policy.maps, prefix + "_wt"), prefix + "_wt", shape, probability=True)
    price = _array(policy.price, "price", (locations,))
    if np.any(price <= 0):
        raise ValueError("Policy prices must be positive")
    projection = float(evaluation.feasibility_projection_mass)
    if not np.isfinite(projection) or projection < 0:
        raise ValueError("Invalid feasibility projection mass")

    def fertility(mass):
        return transition.apply_sequential_fertility(mass, first, P, continuation)

    def transport(mass):
        return model.realize_current_cross_section(
            mass, policy.loc_probs, policy.tenure_choice, policy.tenure_probs,
            policy.maps.lmm_idx, policy.maps.lmm_wt, policy.maps.tmx_idx, policy.maps.tmx_wt,
            use_compiled_scatter=bool(getattr(P, "use_numba_scatter", False)),
            mass_pruning_tolerance=PRUNING_TOLERANCE)

    # Full replay checks that the selected flow uses this evaluation's policies.
    replay_post, all_births, _ = fertility(pre)
    accounting = dict(full_fertility_replay_max_error=_error(replay_post, post),
                      all_birth_flow_error=abs(float(evaluation.births) - all_births))
    del replay_post
    for key, value in accounting.items():
        _check(value, key)
    accounting["full_current_replay_max_error"] = _error(transport(post), current)
    _check(accounting["full_current_replay_max_error"], "full current replay")

    # The official kernel is linear at fixed probabilities. Its original risk
    # pools prevent a new first birth receiving a second birth in this call.
    empty_pre = pre * empty
    empty_post, selected_births, _ = fertility(empty_pre)
    birth_post = empty_post * ((parity > 0) & (child_state == 1))
    control_post = post * empty
    accounting.update(
        empty_control_fertility_max_error=_error(empty_post * empty, control_post),
        selected_birth_flow_error=abs(float(birth_post.sum()) - selected_births),
        empty_partition_max_error=_error(empty_post, birth_post + control_post),
        empty_fertility_mass_error=abs(float(empty_pre.sum() - empty_post.sum())),
    )
    for key, value in accounting.items():
        _check(value, key)
    birth_current = transport(birth_post)
    control_current = current * empty
    accounting.update(
        selected_birth_transport_mass_error=abs(float(birth_current.sum() - birth_post.sum())),
        empty_control_transport_max_error=_error(transport(control_post), control_current),
        selected_birth_subset_max_excess=float(np.maximum(birth_current - current, 0).max()),
        full_fertility_mass_error=abs(float(pre.sum() - post.sum())),
        full_current_mass_error=abs(float(post.sum() - current.sum())),
    )
    for key, value in accounting.items():
        _check(value, key)

    weights = uniform_age_cell_overlap(P, 30., 56.)
    groups = dict(
        selected_birth=_group(birth_current, weights),
        current_empty=_group(control_current, weights),
        first_birth=_group(birth_current, weights, parity == 1),
        continuation_birth=_group(birth_current, weights, parity >= 2),
        empty_never_parent=_group(control_current, weights, never),
        empty_former_parent=_group(control_current, weights, former),
    )
    for key in ("selected_birth", "current_empty"):
        if groups[key]["denominator"] <= 0:
            raise ValueError(f"{key} denominator is nonpositive; no denominator floor")
    accounting.update(occupied_location_probability_sum_error=location_sum_error,
                      all_births=float(all_births), empty_home_births=float(selected_births),
                      pre_empty_mass=float(empty_pre.sum()), feasibility_projection_mass=projection,
                      absolute_tolerance=MASS_ATOL, transport_pruning_tolerance=PRUNING_TOLERANCE)
    source_files = (Path(__file__), Path(transition.__file__), Path(model.__file__),
                    Path(transition.calendar.__file__),
                    Path(model.readiness_childless_states.__code__.co_filename),
                    Path(uniform_age_cell_overlap.__code__.co_filename))
    packet.update(
        status="diagnostic_only_measurement_approximations_unresolved", available=True,
        model_value=groups["selected_birth"]["ownership_rate"] - groups["current_empty"]["ownership_rate"],
        groups=groups, accounting=accounting,
        metadata=dict(
            snapshot=snapshot, age_projection=age_projection,
            diagnostic_allow_residence_proxy=True,
            age_interval=[30., 56.], age_overlap_weights=weights.tolist(),
            age_cell_starts=(18. + 4. * np.arange(17)).tolist(),
            readiness_enabled=bool(model.readiness_gate_active(P)),
            childless_readiness_states=list(childless_states),
            settled_readiness_state=int(model.readiness_settled_state(P)),
            birth_probability_source="supplied policy.fert_probs and owned policy.fert2_probs",
            population_source="evaluation.g_pre after its existing price feasibility gate",
            ownership_definition="realized current tenure index > 0, after location and tenure transactions",
            price=price.tolist(), policy_input_provenance=input_provenance,
            provenance_independently_certified=False, equilibrium_or_stationarity_certified=False,
            source_pins_are_complete_science_manifest=False,
            source_sha256={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
            warnings=[
                "Selected realized-birth ownership contrast; not a causal forced-birth response",
                "Synchronized post-fertility snapshot; annual bridge spreads births over t+1 through t+4",
                "Uniform annual-age exposure and constant within-cell outcomes are explicit approximations",
                "Model dependents proxy ACS resident own children; resident adult children and returns are absent",
                "Current empty control includes former parents and every valid childless readiness state",
                "ACS oldest resident own child under four is not reconstructed at annual interview dates",
                "Independent stochastic child departure is not a deterministic age-18 residence cutoff",
                "At most one explicit birth per four-year period; top-code representative does not scale household mass",
                "National model versus 42-MET2013 ACS sample; DUE UNITSSTR 3:10 restriction lacks a model state",
                "Caller-supplied date, policy and equilibrium provenance is recorded, not independently certified",
                "Empirical target, uncertainty, weights and production objective are unchanged",
            ],
        ),
    )
    return packet
