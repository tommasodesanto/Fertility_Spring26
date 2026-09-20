from __future__ import annotations

import math

import numpy as np

from diagnose_e5f_income_time_aggregation import (
    annual_level_covariance,
    block_average_level_covariance,
    continuous_endpoint_moments,
    final_status,
    main,
    markov_payload_moments,
)


def test_one_year_block_is_the_annual_covariance():
    vp, ve, rho = 0.35, 0.22, 0.81
    for lag in (0, 1, 2, 4):
        np.testing.assert_allclose(
            block_average_level_covariance(vp, ve, rho, 1, lag),
            annual_level_covariance(vp, ve, rho, lag),
            atol=1e-13,
        )


def test_four_year_iid_limit_has_variance_divided_by_four():
    ve = 0.48
    expected = (math.exp(ve) - 1.0) / 4.0
    np.testing.assert_allclose(block_average_level_covariance(0.0, ve, 0.7, 4, 0), expected)
    for lag in (1, 2, 4):
        np.testing.assert_allclose(block_average_level_covariance(0.0, ve, 0.7, 4, lag), 0.0)


def test_zero_risk_process_has_zero_level_covariance_and_unit_mean():
    assert annual_level_covariance(0.0, 0.0, 0.9, 0) == 0.0
    assert block_average_level_covariance(0.0, 0.0, 0.9, 4, 4) == 0.0
    continuous = continuous_endpoint_moments(0.0, 0.0, 0.9, (0, 1, 2, 4))
    assert continuous["mean"]["value"] == 1.0
    for row in continuous["moments"].values():
        assert row["log_covariance"] == 0.0
        assert row["level_covariance"] == 0.0


def test_markov_moments_use_stationary_transition_powers():
    # A two-state symmetric Markov chain gives a hand-checkable covariance.
    z = np.array([0.5, 1.5])
    weights = np.array([0.5, 0.5])
    pi = np.array([[0.75, 0.25], [0.25, 0.75]])
    result = markov_payload_moments(z, weights, pi, (0, 1, 2))
    centered = z - 1.0
    np.testing.assert_allclose(result["mean"]["value"], 1.0)
    np.testing.assert_allclose(result["moments"]["0"]["level_covariance"], 0.25)
    np.testing.assert_allclose(result["moments"]["1"]["level_covariance"], 0.125)
    np.testing.assert_allclose(result["moments"]["2"]["level_covariance"], 0.0625)


def test_failed_gate_returns_failed_status_without_running_simulation():
    assert final_status(
        requested_batches=2,
        completed_batches=2,
        exact_checks_pass=False,
        mean_checks_pass=True,
        plots_present=True,
        dependencies_unchanged=True,
        stop_reason=None,
    ) == "failed"
    assert final_status(
        requested_batches=2,
        completed_batches=2,
        exact_checks_pass=True,
        mean_checks_pass=False,
        plots_present=True,
        dependencies_unchanged=True,
        stop_reason=None,
    ) == "failed"
    assert final_status(
        requested_batches=2,
        completed_batches=2,
        exact_checks_pass=True,
        mean_checks_pass=True,
        plots_present=False,
        dependencies_unchanged=True,
        stop_reason=None,
    ) == "failed"


def test_main_returns_nonzero_for_injected_failed_receipt(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "diagnose_e5f_income_time_aggregation.run_diagnostic",
        lambda *args, **kwargs: {"status": "failed"},
    )
    assert main(["--mode", "smoke", "--output", str(tmp_path / "unused")]) == 1
