from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np


TOOLS = Path(__file__).resolve().parent
MODEL = TOOLS.parent
for path in (MODEL, TOOLS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as calibration
import collect_e5f_transition_calibration as collector
from intergen_eqscale_seq_optimized import parameters


def test_replacement_conversion_and_old_stationary_identity() -> None:
    conversion = transition.effective_birth_to_household_conversion(2.1)
    assert conversion == 1.0 / 2.1
    accounting = transition.stationary_renewal_from_births(
        entry_flow=0.05,
        topcode_adjusted_births=2.1 * 0.05,
        outside_origin_entry_share=0.169,
        birth_to_household_conversion=conversion,
    )
    assert math.isclose(accounting["queue_B_over_E"], 1.0, abs_tol=1e-15)
    assert math.isclose(accounting["retention_rho"], 0.831, abs_tol=1e-15)
    assert math.isclose(accounting["identity_residual"], 0.0, abs_tol=1e-15)


def test_historical_bridge_is_invariant_to_open_share() -> None:
    E = 0.0617334556564513
    B = E
    conversion = 1.0 / 2.1
    for share in (1e-6, 0.05, 0.169, 0.6, 0.999999):
        accounting = transition.stationary_renewal_from_births(
            entry_flow=E,
            topcode_adjusted_births=B / conversion,
            outside_origin_entry_share=share,
            birth_to_household_conversion=conversion,
        )
        entrants = (
            accounting["outside_flow_M"]
            + accounting["retention_rho"] * B
        )
        assert math.isclose(entrants, E, rel_tol=0.0, abs_tol=1e-16)

    # The observed-age bridge changes only age-cell scale. Two provisional
    # distributions with the same within-age composition must therefore map
    # to the same dated state even when their pre-bridge age masses differ.
    conditional = np.arange(1.0, 9.0).reshape(1, 2, 1, 1, 2, 2, 1)
    ages = 18.0 + 4.0 * np.arange(17)
    first = conditional * np.linspace(0.5, 1.5, 17).reshape(1, 1, 1, 17, 1, 1, 1)
    second = conditional * np.linspace(1.7, 0.3, 17).reshape(1, 1, 1, 17, 1, 1, 1)
    first_bridged, _ = transition.reweight_distribution_to_observed_age_path(
        first, ages, year=2011, initial_mass=1.0
    )
    second_bridged, _ = transition.reweight_distribution_to_observed_age_path(
        second, ages, year=2011, initial_mass=1.0
    )
    assert np.allclose(first_bridged, second_bridged, rtol=0.0, atol=2e-16)


def test_collector_rejects_wrong_replacement_contract() -> None:
    E = 0.05
    B = E
    share = 0.169
    summary = {
        "old_model_completed_fertility": 2.1,
        "renewal_accounting_contract": {
            "replacement_fertility": 2.1,
            "effective_birth_to_household_conversion": 1.0 / 2.1,
            "birth_vintage_queue_waiting_slots": 4,
            "birth_to_entry_effect_lag_dates": 5,
            "birth_to_entry_effect_lag_years": 20.0,
        },
        "renewal_accounting_old_state": {
            "old_entry_flow_E": E,
            "old_queue_mature_flow_B": B,
            "old_queue_B_over_E": 1.0,
            "old_renewal_residual": 0.0,
            "outside_flow_M": share * E,
            "outside_origin_entry_share": share,
            "retention_rho": 1.0 - share,
        },
    }
    collector.validate_renewal_accounting(summary)
    summary["renewal_accounting_contract"]["replacement_fertility"] = 2.12
    try:
        collector.validate_renewal_accounting(summary)
    except RuntimeError:
        pass
    else:
        raise AssertionError("collector accepted a mixed replacement contract")


def test_collector_uses_driver_gate_for_independent_fertility_measurement() -> None:
    E = 0.05
    share = 0.169
    summary = {
        "old_model_completed_fertility": 2.1 + 2.1e-8,
        "renewal_accounting_contract": {
            "replacement_fertility": 2.1,
            "effective_birth_to_household_conversion": 1.0 / 2.1,
            "birth_vintage_queue_waiting_slots": 4,
            "birth_to_entry_effect_lag_dates": 5,
            "birth_to_entry_effect_lag_years": 20.0,
        },
        "renewal_accounting_old_state": {
            "old_entry_flow_E": E,
            "old_queue_mature_flow_B": E,
            "old_queue_B_over_E": 1.0,
            "old_renewal_residual": 0.0,
            "outside_flow_M": share * E,
            "outside_origin_entry_share": share,
            "retention_rho": 1.0 - share,
        },
    }
    collector.validate_renewal_accounting(summary)
    summary["old_model_completed_fertility"] = 2.1 + 5.0e-8
    try:
        collector.validate_renewal_accounting(summary)
    except RuntimeError:
        pass
    else:
        raise AssertionError("collector accepted a fertility identity outside 2e-8")


def test_collector_rejects_post_target_or_policy_calibration() -> None:
    collector.validate_calibration_scope(
        {"post_2023_periods": 0, "policy_case": "none"}
    )
    for bad in (
        {"post_2023_periods": 1, "policy_case": "none"},
        {"post_2023_periods": 0, "policy_case": "dependent-child-ltv95"},
    ):
        try:
            collector.validate_calibration_scope(bad)
        except RuntimeError:
            pass
        else:
            raise AssertionError("collector accepted an out-of-scope calibration task")


def test_collector_reconstructs_target_loss_and_fingerprint() -> None:
    rows = [
        {
            "candidate": "task_001",
            "moment": "a",
            "target": "1.0",
            "model": "1.5",
            "gap": "0.5",
            "weight": "4.0",
            "loss_contribution": "1.0",
            "standardized_gap": "1.0",
        },
        {
            "candidate": "task_001",
            "moment": "b",
            "target": "2.0",
            "model": "1.0",
            "gap": "-1.0",
            "weight": "9.0",
            "loss_contribution": "9.0",
            "standardized_gap": "-3.0",
        },
    ]
    payload = [["a", 1.0, 4.0], ["b", 2.0, 9.0]]
    fingerprint = hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    summary = {"target_count": 2, "target_fingerprint": fingerprint}
    candidate = {"candidate": "task_001", "transition_loss": 10.0}
    assert collector.validate_target_fit(rows, summary, candidate) == {
        "a": 1.5,
        "b": 1.0,
    }
    rows[0]["loss_contribution"] = "0.9"
    try:
        collector.validate_target_fit(rows, summary, candidate)
    except RuntimeError:
        pass
    else:
        raise AssertionError("collector accepted corrupted target-loss algebra")


def test_collector_gates_dated_housing_ledger() -> None:
    fields = [
        "candidate",
        "period",
        "status",
        "origin_period",
        "destination_period",
        "origin_price",
        "destination_price",
        "origin_mass",
        "destination_mass",
        "survival_loss",
        "treated_feasibility_projection_mass",
        "control_feasibility_projection_mass",
        "treated_continuation_births",
        "treated_mean_housing",
        "control_mean_housing",
        "housing_response",
        "census_age_bridge_applied",
    ]
    rows = []
    for period in range(5):
        rows.append(
            {
                "candidate": "task_001",
                "period": period,
                "status": (
                    "same_policy_one_period_diagnostic_no_prior_date"
                    if period == 0
                    else "complete_one_period_dated_matched_branch_did"
                ),
                "origin_period": "" if period == 0 else period - 1,
                "destination_period": period,
                "origin_price": "" if period == 0 else 1.0,
                "destination_price": "" if period == 0 else 1.1,
                "origin_mass": "" if period == 0 else 0.2,
                "destination_mass": "" if period == 0 else 0.19,
                "survival_loss": "" if period == 0 else 0.01,
                "treated_feasibility_projection_mass": "" if period == 0 else 0.0,
                "control_feasibility_projection_mass": "" if period == 0 else 0.0,
                "treated_continuation_births": "" if period == 0 else 0.02,
                "treated_mean_housing": "" if period == 0 else 6.0,
                "control_mean_housing": "" if period == 0 else 5.7 - 0.1 * period,
                "housing_response": 0.3 + 0.1 * period,
                "census_age_bridge_applied": "False",
            }
        )
    summary = {
        "dated_measurement_contract": {
            "first_birth_housing_horizon_periods": 1,
            "first_birth_housing_horizon_years": 4.0,
            "first_birth_housing_terminal_origin_year": 2019,
            "first_birth_housing_terminal_destination_year": 2023,
        }
    }
    candidate = {
        "candidate": "task_001",
        "terminal_first_birth_housing_response": 0.7,
    }
    with tempfile.TemporaryDirectory() as temporary:
        task_dir = Path(temporary)
        case_dir = task_dir / "cases" / "task_001"
        case_dir.mkdir(parents=True)
        ledger = case_dir / "dated_first_birth_housing_did.csv"

        def save() -> None:
            with ledger.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)

        save()
        collector.validate_dated_housing_ledger(
            task_dir,
            summary,
            candidate,
            {"housing_increment_0to1": 0.7},
        )
        rows[-1]["census_age_bridge_applied"] = "True"
        save()
        try:
            collector.validate_dated_housing_ledger(
                task_dir,
                summary,
                candidate,
                {"housing_increment_0to1": 0.7},
            )
        except RuntimeError:
            pass
        else:
            raise AssertionError("collector accepted a bridged housing contrast")


def test_population_bridge_contract_ignores_candidate_diagnostics_only() -> None:
    initial = {
        "mapping": "four-year bins",
        "model_age_grid": [18.0, 22.0],
        "sample": "heads",
        "source": "ACS",
        "source_path": "/data/acs.dta",
        "source_receipt": {"sample": 200701},
        "source_sha256": "a" * 64,
        "source_year": 2007,
        "units": "households",
        "period_years": 4.0,
        "horizon_periods": 4,
        "horizon_calendar_year": 2023.0,
        "future_entry_rule": "fixed",
        "stationary_age_mass": [0.4, 0.6],
        "reweighted_initial_age_mass": [0.3, 0.7],
    }
    bridge = {
        "status": "externally_normalized_not_estimated",
        "target_indices": {"0": 1.0, "4": 1.1},
        "total_households": "Census HH-3",
        "age_distribution": "ACS head shares",
        "initial_age_reweight": initial,
    }
    first = collector.population_bridge_contract(bridge)
    changed_diagnostics = {
        **bridge,
        "initial_age_reweight": {
            **initial,
            "stationary_age_mass": [0.2, 0.8],
            "reweighted_initial_age_mass": [0.1, 0.9],
        },
    }
    assert collector.population_bridge_contract(changed_diagnostics) == first
    changed_source = {
        **bridge,
        "initial_age_reweight": {**initial, "source_sha256": "b" * 64},
    }
    assert collector.population_bridge_contract(changed_source) != first


def test_birth_vintage_impulse_reaches_entry_after_twenty_years() -> None:
    conversion = 1.0 / 2.1
    baseline_births = 2.1
    queue = [conversion * baseline_births] * 4
    due_flows = []
    for period in range(5):
        births = baseline_births + (0.21 if period == 0 else 0.0)
        due, queue = transition.advance_birth_vintage_queue(
            queue, births, conversion
        )
        due_flows.append(due)
        assert len(queue) == 4
    assert due_flows[:4] == [1.0] * 4
    assert math.isclose(due_flows[4], 1.1, abs_tol=1e-15)
    # The flow popped during t=4 is injected into the next dated distribution,
    # so a birth during t=0 first changes entrants at t=5 (20 years later).


def test_fixed_cohort_timing_uses_hazards_not_source_cell_masses() -> None:
    original_model = calibration.calendar.model
    calibration.calendar.model = SimpleNamespace(
        age_to_index=lambda _P, age: int((age - 18.0) / 4.0)
    )
    P = SimpleNamespace(age_start=18.0, da=4.0, J=7)
    old = {
        "hazard": np.asarray([0.10, 0.20, 0.0, 0.0, 0.0, 0.0, 0.0]),
        "at_risk": np.asarray([1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        "flow": np.asarray([0.10, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    }
    dated_hazards = [0.30, 0.40, 0.50, 0.60, 0.70]
    records = {}
    for period, hazard in enumerate(dated_hazards):
        age_index = 2 + period
        hazard_vector = np.zeros(7)
        risk_vector = np.zeros(7)
        flow_vector = np.zeros(7)
        hazard_vector[age_index] = hazard
        risk_vector[age_index] = 100.0 * (period + 1)
        flow_vector[age_index] = hazard * risk_vector[age_index]
        records[period] = {
            "first_birth_accounting_by_age": {
                "hazard": hazard_vector,
                "at_risk": risk_vector,
                "flow": flow_vector,
            }
        }
    try:
        result = calibration.cohort_timing_moments(old, records, P)
    finally:
        calibration.calendar.model = original_model
    hazards = [0.10, 0.20, *dated_hazards]
    survivor = 1.0
    expected = []
    for hazard in hazards:
        expected.append(survivor * hazard)
        survivor *= 1.0 - hazard
    ages = np.arange(20.0, 45.0, 4.0)
    expected_mean = float(np.dot(ages, expected) / np.sum(expected))
    expected_late = float(
        np.sum(np.asarray(expected)[ages >= 30.0]) / np.sum(expected)
    )
    assert math.isclose(result["mean_age_first_birth"], expected_mean, abs_tol=1e-14)
    assert math.isclose(
        result["share_first_births_age30plus"], expected_late, abs_tol=1e-14
    )
    assert result["first_birth_age_measurement"] == (
        "four-year model cells labeled by interval midpoints"
    )
    assert math.isclose(result["synthetic_childless_probability"], survivor, abs_tol=1e-14)

    normalized = calibration.cohort_timing_moments(
        old,
        records,
        P,
        terminal_childless_probability=survivor + 5.0e-10,
    )
    assert normalized["synthetic_childless_probability"] == survivor + 5.0e-10
    assert math.isclose(
        normalized["synthetic_ever_first_birth_probability"]
        + normalized["synthetic_childless_probability"],
        1.0,
        rel_tol=0.0,
        abs_tol=1e-16,
    )
    assert math.isclose(normalized["mean_age_first_birth"], expected_mean, abs_tol=1e-14)
    gate = normalized["terminal_childless_identity_normalization"]
    assert gate["status"] == "terminal_stock_identity_after_2e-9_roundoff_gate"
    try:
        calibration.cohort_timing_moments(
            old,
            records,
            P,
            terminal_childless_probability=survivor + 3.0e-9,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("timing identity normalization accepted a material gap")


def test_period_tfr_reconstructs_parity_flows() -> None:
    top_weight = 3.602359422009
    dummy_model = SimpleNamespace(
        get_fecundity_by_age=lambda P: np.ones(P.J),
        readiness_settled_state=lambda P: 0,
    )
    original_model = calibration.calendar.model
    calibration.calendar.model = dummy_model
    P = SimpleNamespace(
        J=3,
        A_f_start=1,
        A_f_end=3,
        n_parity=4,
        age_start=18.0,
        da=4.0,
        tfr_top_bin_weight=top_weight,
        fertility_units="literal_topcode",
    )
    P._fert2_probs = np.zeros((1, 1, 1, 3, 1, 2, 2, 4))
    P._fert2_probs[:, :, :, :, :, 1, 0, :] = 0.10
    P._fert2_probs[:, :, :, :, :, 1, 1, :] = 0.05
    g_pre = np.zeros((1, 1, 1, 3, 1, 4, 4))
    g_pre[:, :, :, :, :, 0, 0] = 0.50
    g_pre[:, :, :, :, :, 1, 0] = 0.30
    g_pre[:, :, :, :, :, 2, 0] = 0.20
    fert_probs = np.zeros((1, 1, 1, 3, 1, 2))
    fert_probs[:, :, :, :, :, 1] = 0.20
    explicit_by_age = 0.10 + 0.03 + 0.01
    g_post = g_pre.copy()
    g_post[:, :, :, :, :, 3, 0] += 0.01
    evaluation = SimpleNamespace(
        g_pre=g_pre,
        g_post_fertility=g_post,
        births=3.0 * explicit_by_age,
        policy=SimpleNamespace(fert_probs=fert_probs, fert2_probs=P._fert2_probs.copy()),
    )
    try:
        result = calibration.period_fertility_diagnostics(evaluation, P)
    finally:
        calibration.calendar.model = original_model
    expected_adjusted_by_age = 0.10 + 0.03 + (top_weight - 2.0) * 0.01
    assert math.isclose(
        result["period_tfr_explicit_states"],
        3.0 * explicit_by_age,
        abs_tol=1e-14,
    )
    assert math.isclose(
        result["period_tfr_topcode_adjusted"],
        3.0 * expected_adjusted_by_age,
        abs_tol=1e-14,
    )


def test_first_birth_housing_response_uses_one_four_year_branch() -> None:
    calls: list[tuple[int, int]] = []

    def advance(cohort, start_age, horizon, *_args):
        calls.append((int(start_age), int(horizon)))
        return cohort

    def mean_housing(cohort, age, *_args):
        assert age == 1
        parent_mass = float(np.sum(cohort[..., 1, :]))
        return 5.0 if parent_mass > 0.0 else 2.0

    dummy_model = SimpleNamespace(
        get_fecundity_by_age=lambda _P: np.asarray([0.5, 0.0]),
        readiness_settled_state=lambda _P: 0,
        income_transition_values=lambda _P: (
            np.asarray([1.0]),
            np.asarray([1.0]),
            np.eye(1),
        ),
        advance_cohort_horizon_markov_income=advance,
        realize_current_choices_markov_income=lambda cohort, *_args, **_kwargs: cohort,
        mean_housing_distribution_markov=mean_housing,
    )
    original_model = calibration.calendar.model
    calibration.calendar.model = dummy_model
    P = SimpleNamespace(
        J=2,
        A_f_start=1,
        A_f_end=1,
        n_house=0,
        I=1,
        n_parity=2,
        n_child_states=2,
        use_stochastic_aging=False,
        use_numba_scatter=False,
        housing_event_horizon=1,
    )
    g_pre = np.zeros((1, 1, 1, 2, 1, 2, 2))
    g_pre[0, 0, 0, 0, 0, 0, 0] = 2.0
    fert_probs = np.zeros((1, 1, 1, 2, 1, 2))
    fert_probs[0, 0, 0, 0, 0, 1] = 1.0
    maps = SimpleNamespace(
        lmm_idx=np.zeros((1, 1, 1), dtype=int),
        lmm_wt=np.zeros((1, 1, 1)),
        tmx_idx=np.zeros((1, 1, 1, 2, 2, 1), dtype=int),
        tmx_wt=np.zeros((1, 1, 1, 2, 2, 1)),
    )
    policy = SimpleNamespace(
        bp_pol=np.zeros((1, 1, 1, 2, 1, 2, 2)),
        fert_probs=fert_probs,
        loc_probs=np.zeros(1),
        tenure_choice=np.zeros(1),
        tenure_probs=None,
        hR_pol=np.zeros(1),
        maps=maps,
    )
    evaluation = SimpleNamespace(g_pre=g_pre, policy=policy)
    try:
        response = calibration.first_birth_housing_response(
            evaluation,
            P,
            np.asarray([0.0]),
            SimpleNamespace(nc=4),
        )
        assert math.isclose(response, 3.0, abs_tol=1e-15)
        assert calls == [(0, 1), (0, 1)]
        P.housing_event_horizon = 0
        try:
            calibration.first_birth_housing_response(
                evaluation,
                P,
                np.asarray([0.0]),
                SimpleNamespace(nc=4),
            )
        except ValueError:
            pass
        else:
            raise AssertionError("housing response accepted the contemporaneous clock")
    finally:
        calibration.calendar.model = original_model


def test_dated_first_birth_branch_uses_origin_and_destination_kernels() -> None:
    advance_calls: list[int] = []
    advance_policy_sentinels: list[float] = []
    destination_policy_sentinels: list[float] = []
    fertility_calls: list[float] = []

    def advance(cohort, age, loc_probs, *_args):
        advance_calls.append(int(age))
        advance_policy_sentinels.append(float(np.asarray(loc_probs).reshape(-1)[0]))
        return cohort

    def apply_continuation(cohort, fert_probs, parameters, fert2_probs):
        # The control branch is never passed through this fertility operator.
        assert float(np.sum(cohort[..., 0, :])) == 0.0
        assert float(np.asarray(fert_probs).reshape(-1)[0]) == 33.0
        assert float(np.asarray(fert2_probs).reshape(-1)[0]) == 77.0
        assert float(np.asarray(parameters._fert2_probs).reshape(-1)[0]) == -99.0
        fertility_calls.append(float(np.sum(cohort)))
        return cohort.copy(), 0.1, np.asarray([0.1])

    def realize_destination(cohort, loc_probs, *_args, **_kwargs):
        destination_policy_sentinels.append(
            float(np.asarray(loc_probs).reshape(-1)[0])
        )
        return cohort

    def demand(cohort, _policy, _P):
        parent_mass = float(np.sum(cohort[..., 1:, :]))
        mass = float(np.sum(cohort))
        return np.asarray([(5.0 if parent_mass > 0.0 else 2.0) * mass])

    dummy_model = SimpleNamespace(
        get_fecundity_by_age=lambda _P: np.asarray([0.5, 0.0]),
        readiness_settled_state=lambda _P: 0,
        income_transition_values=lambda _P: (
            np.asarray([1.0]),
            np.asarray([1.0]),
            np.eye(1),
        ),
        advance_cohort_one_period_markov_income=advance,
        realize_current_cross_section=realize_destination,
    )
    original_model = calibration.calendar.model
    original_gate = calibration.calendar.gate_pre_fertility_distribution
    original_demand = calibration.calendar.housing_demand_by_location
    original_apply = calibration.transition.apply_sequential_fertility
    calibration.calendar.model = dummy_model
    calibration.calendar.gate_pre_fertility_distribution = (
        lambda cohort, *_args: (np.asarray(cohort).copy(), 0.0)
    )
    calibration.calendar.housing_demand_by_location = demand
    calibration.transition.apply_sequential_fertility = apply_continuation
    P = SimpleNamespace(
        J=2,
        A_f_start=1,
        A_f_end=1,
        n_house=0,
        I=1,
        n_parity=2,
        n_child_states=2,
        use_stochastic_aging=False,
        use_numba_scatter=False,
        use_age_survival=True,
        survival_probs=np.asarray([0.8]),
        _fert2_probs=np.asarray([-99.0]),
    )
    g_pre = np.zeros((1, 1, 1, 2, 1, 2, 2))
    g_pre[0, 0, 0, 0, 0, 0, 0] = 2.0
    fert_probs = np.zeros((1, 1, 1, 2, 1, 2))
    fert_probs[0, 0, 0, 0, 0, 1] = 1.0
    maps = SimpleNamespace(
        lmm_idx=np.zeros((1, 1, 1), dtype=int),
        lmm_wt=np.zeros((1, 1, 1)),
        tmx_idx=np.zeros((1, 1, 1, 2, 2, 1), dtype=int),
        tmx_wt=np.zeros((1, 1, 1, 2, 2, 1)),
    )
    policy = SimpleNamespace(
        price=np.asarray([1.0]),
        bp_pol=np.zeros((1, 1, 1, 2, 1, 2, 2)),
        fert_probs=fert_probs,
        loc_probs=np.asarray([11.0]),
        tenure_choice=np.zeros(1),
        tenure_probs=None,
        hR_pol=np.zeros(1),
        maps=maps,
    )
    origin = SimpleNamespace(g_pre=g_pre, policy=policy)
    destination_policy = SimpleNamespace(**vars(policy))
    destination_policy.price = np.asarray([2.0])
    destination_policy.loc_probs = np.asarray([22.0])
    destination_policy.fert_probs = np.full_like(fert_probs, 33.0)
    destination_policy.fert2_probs = np.asarray([77.0])
    destination = SimpleNamespace(policy=destination_policy)
    try:
        branch = calibration.begin_dated_first_birth_housing_branch(
            origin,
            P,
            np.asarray([0.0]),
            SimpleNamespace(nc=4),
            origin_period=3,
        )
        assert math.isclose(branch["origin_mass"], 1.0, abs_tol=1e-15)
        assert math.isclose(
            branch["survivor_mass_before_destination_gating"],
            0.8,
            abs_tol=1e-15,
        )
        result = calibration.finish_dated_first_birth_housing_branch(
            branch,
            destination,
            P,
            np.asarray([0.0]),
            SimpleNamespace(nc=4),
            destination_period=4,
        )
        assert advance_calls == [0, 0]
        assert advance_policy_sentinels == [11.0, 11.0]
        assert destination_policy_sentinels == [22.0, 22.0]
        assert fertility_calls == [0.8]
        assert math.isclose(result["housing_response"], 3.0, abs_tol=1e-15)
        assert math.isclose(result["treated_continuation_births"], 0.1)
        assert result["origin_period"] == 3
        assert result["destination_period"] == 4
        assert result["census_age_bridge_applied"] is False
    finally:
        calibration.calendar.model = original_model
        calibration.calendar.gate_pre_fertility_distribution = original_gate
        calibration.calendar.housing_demand_by_location = original_demand
        calibration.transition.apply_sequential_fertility = original_apply


def test_parameter_table_uses_actual_transition_domain() -> None:
    original_domain = calibration.TRANSITION_SEARCH_DOMAIN
    calibration.TRANSITION_SEARCH_DOMAIN = (
        ("beta_annual", 0.94, 0.9995, "discount"),
        ("first_birth_fixed_cost", 0.0, 8.0, "softzero"),
        ("psi_child_change_2023", -1.5, 0.2, "asinh"),
    )
    try:
        rows = calibration.parameter_rows(
            {"beta": 0.98**4, "first_birth_fixed_cost": 2.5},
            {"new_psi_child": -0.4},
            old_psi_child=0.1,
            joint_panel=True,
        )
    finally:
        calibration.TRANSITION_SEARCH_DOMAIN = original_domain
    free = [row for row in rows if row["is_free_parameter"]]
    assert [row["parameter"] for row in free] == [
        "beta_annual",
        "first_birth_fixed_cost",
        "psi_child_change_2023",
    ]
    assert math.isclose(free[0]["value"], 0.98, abs_tol=1e-15)
    assert math.isclose(free[2]["value"], -0.5, abs_tol=1e-15)
    assert free[2]["lower_bound"] == -1.5
    assert free[2]["upper_bound"] == 0.2
    assert all(row["status"] == "estimated_free_transition_parameter" for row in free)
    derived = [row for row in rows if not row["is_free_parameter"]]
    assert [row["parameter"] for row in derived] == [
        "psi_child_2007",
        "psi_child_2023",
    ]


def test_first_child_jump_can_be_fixed_or_jointly_estimated_but_not_both() -> None:
    base_domain = (("first_birth_fixed_cost", 0.0, 8.0, "softzero"),)
    theta: dict[str, float] = {}
    domain, overrides, metadata = calibration.configure_first_child_room_jump(
        model_profile_name=calibration.REPAIRED_MODEL_PROFILE,
        theta=theta,
        active_domain=base_domain,
        profile_overrides={},
        model_profile={"name": calibration.REPAIRED_MODEL_PROFILE},
        fixed_jump=None,
        estimate_jump=True,
    )
    assert domain[-1] == ("hbar_first_child_jump", 0.0, 0.5, "softzero")
    assert theta["hbar_first_child_jump"] == 0.0
    assert overrides["hbar_first_child_jump"] == 0.0
    assert metadata["first_child_room_jump_status"] == (
        "jointly estimated transition parameter"
    )

    fixed_theta: dict[str, float] = {}
    fixed_domain, fixed_overrides, fixed_metadata = (
        calibration.configure_first_child_room_jump(
            model_profile_name=calibration.REPAIRED_MODEL_PROFILE,
            theta=fixed_theta,
            active_domain=base_domain,
            profile_overrides={},
            model_profile={"name": calibration.REPAIRED_MODEL_PROFILE},
            fixed_jump=0.2,
            estimate_jump=False,
        )
    )
    assert fixed_domain == base_domain
    assert fixed_theta["hbar_first_child_jump"] == 0.2
    assert fixed_overrides["hbar_first_child_jump"] == 0.2
    assert fixed_metadata["first_child_room_jump"] == 0.2

    try:
        calibration.configure_first_child_room_jump(
            model_profile_name=calibration.REPAIRED_MODEL_PROFILE,
            theta={},
            active_domain=base_domain,
            profile_overrides={},
            model_profile={},
            fixed_jump=0.2,
            estimate_jump=True,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("fixed and estimated jump must be mutually exclusive")


def test_structural_age_recursion_removes_only_solver_roundoff() -> None:
    entry = 0.06
    survival = np.asarray([1.0, 1.0, 0.9, 0.8])
    exact = np.asarray([entry, entry, entry, 0.9 * entry, 0.72 * entry])
    noisy = exact * np.asarray([1.0, 1.0 + 4e-10, 1.0 - 3e-10, 1.0, 1.0])
    reconstructed, gate = calibration.validated_structural_stationary_age_mass(
        noisy,
        entry_flow=entry,
        structural_survival=survival,
    )
    assert np.allclose(reconstructed, exact, rtol=0.0, atol=1e-17)
    assert gate["max_relative_mass_gap_to_entry"] < 5e-9
    materially_wrong = noisy.copy()
    materially_wrong[2] *= 1.0 + 1e-6
    try:
        calibration.validated_structural_stationary_age_mass(
            materially_wrong,
            entry_flow=entry,
            structural_survival=survival,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("structural age recursion accepted a material KFE mismatch")

    ages = 18.0 + 4.0 * np.arange(17)
    observed = np.full(17, entry)
    observed[4] *= 1.0 + 4.0e-10
    diagnostic = transition.acs_2007_age_reweight_diagnostic(
        observed,
        ages,
        entry,
        structural_survival=np.ones(16),
    )
    assert diagnostic["age_survival"] == [1.0] * 16
    assert diagnostic["survival_source"] == (
        "structural_survival_after_1e-8_age_mass_roundoff_gate"
    )
    assert math.isclose(
        sum(diagnostic["reweighted_initial_age_mass"]),
        float(np.sum(observed)),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_branch_mass_normalization_is_strict_and_mass_preserving() -> None:
    raw = np.asarray([0.01, 0.02, 0.010000000024])
    normalized, gate = calibration.normalize_distribution_mass_roundoff(
        raw,
        expected_mass=0.04,
        stage="unit_test",
    )
    assert math.isclose(float(np.sum(normalized)), 0.04, rel_tol=0.0, abs_tol=1e-17)
    assert gate["relative_gap"] < 5e-9
    try:
        calibration.normalize_distribution_mass_roundoff(
            np.asarray([0.01, 0.02, 0.010001]),
            expected_mass=0.04,
            stage="unit_test_bad",
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("branch normalization accepted a material mass error")


def test_calendar_mass_normalization_is_default_off_and_fail_closed() -> None:
    P = parameters.setup_parameters()
    assert P.normalize_transition_mass_roundoff is False
    raw = np.asarray([0.01, 0.02, 0.0100000003])
    normalized = transition.normalize_pure_transition_mass_roundoff(
        raw,
        expected_mass=0.04,
        stage="unit_test_calendar",
    )
    assert math.isclose(float(np.sum(normalized)), 0.04, rel_tol=0.0, abs_tol=1e-17)
    try:
        transition.normalize_pure_transition_mass_roundoff(
            np.asarray([0.01, 0.02, 0.0100000005]),
            expected_mass=0.04,
            stage="unit_test_calendar_bad",
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("calendar normalization accepted a material mass error")


if __name__ == "__main__":
    test_replacement_conversion_and_old_stationary_identity()
    test_historical_bridge_is_invariant_to_open_share()
    test_collector_rejects_wrong_replacement_contract()
    test_collector_uses_driver_gate_for_independent_fertility_measurement()
    test_collector_rejects_post_target_or_policy_calibration()
    test_collector_reconstructs_target_loss_and_fingerprint()
    test_collector_gates_dated_housing_ledger()
    test_population_bridge_contract_ignores_candidate_diagnostics_only()
    test_birth_vintage_impulse_reaches_entry_after_twenty_years()
    test_fixed_cohort_timing_uses_hazards_not_source_cell_masses()
    test_period_tfr_reconstructs_parity_flows()
    test_first_birth_housing_response_uses_one_four_year_branch()
    test_dated_first_birth_branch_uses_origin_and_destination_kernels()
    test_parameter_table_uses_actual_transition_domain()
    test_first_child_jump_can_be_fixed_or_jointly_estimated_but_not_both()
    test_structural_age_recursion_removes_only_solver_roundoff()
    test_branch_mass_normalization_is_strict_and_mass_preserving()
    test_calendar_mass_normalization_is_default_off_and_fail_closed()
    print("E5F_TRANSITION_ACCOUNTING_TESTS_PASS tests=18")
