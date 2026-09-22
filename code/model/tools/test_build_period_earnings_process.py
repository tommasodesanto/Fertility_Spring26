from __future__ import annotations

import numpy as np
import pytest

from build_period_earnings_process import build_period_earnings_process


def _candidate(**kwargs):
    values = dict(
        rho_period=0.81,
        persistent_innovation_sd_period=0.38,
        transitory_sd_period=0.085,
        n_persistent=5,
        n_iid=3,
    )
    values.update(kwargs)
    return build_period_earnings_process(**values)


def test_four_year_constructor_is_stochastic_stationary_and_mean_one():
    overrides, meta = _candidate()
    z, w, pi = overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"]
    assert z.size == 15
    np.testing.assert_allclose(w.sum(), 1.0, atol=1e-14)
    np.testing.assert_allclose(w @ z, 1.0, atol=1e-14)
    np.testing.assert_allclose(pi.sum(axis=1), 1.0, atol=1e-14)
    np.testing.assert_allclose(w @ pi, w, atol=1e-12)
    assert meta["frequency"] == "four_year"
    assert meta["period_years"] == 4
    assert meta["stationarity"]
    assert meta["mean"] == pytest.approx(1.0, abs=1e-14)


def test_log_covariances_match_rho_power_and_iid_rows_are_independent():
    rho, persistent_sd, transitory_sd = 0.81, 0.38, 0.085
    overrides, meta = _candidate()
    vp = persistent_sd**2 / (1.0 - rho**2)
    expected = np.array([vp + transitory_sd**2] + [vp * rho**lag for lag in range(1, 5)])
    np.testing.assert_allclose(meta["continuous_period_log_covariances"], expected, atol=1e-14)
    np.testing.assert_allclose(meta["discrete_period_log_covariances"], expected, atol=1e-12)
    reshaped = overrides["Pi_z"].reshape(5, 3, 5, 3)
    for previous_persistent in range(5):
        for next_persistent in range(5):
            np.testing.assert_allclose(
                reshaped[previous_persistent, 0, next_persistent, :],
                reshaped[previous_persistent, 2, next_persistent, :],
                atol=1e-14,
            )


def test_component_order_is_persistent_then_iid():
    overrides, meta = _candidate(n_persistent=5, n_iid=5)
    components = meta["components"]
    np.testing.assert_allclose(
        overrides["z_grid"].reshape(5, 5),
        np.multiply.outer(components["persistent_levels"], components["iid_levels"]),
    )
    np.testing.assert_allclose(
        overrides["z_weights"].reshape(5, 5),
        np.multiply.outer(components["persistent_weights"], components["iid_weights"]),
    )


def test_covariance_errors_are_recorded_for_requested_resolution_grid():
    _, meta = _candidate(n_persistent=9, n_iid=5)
    rows = meta["resolution_table"]
    assert {(row["n_persistent"], row["n_iid"]) for row in rows} == {
        (5, 3), (5, 5), (9, 3), (9, 5), (15, 3), (15, 5)
    }
    for row in rows:
        assert len(row["errors_vs_continuous"]) == 5
        assert np.isfinite(row["maximum_absolute_error"])
    errors = {
        (row["n_persistent"], row["n_iid"]): row["maximum_absolute_error"]
        for row in rows
    }
    assert errors[(15, 5)] <= errors[(5, 3)] + 1e-14


def test_zero_iid_supports_a_pure_persistent_ar1_with_one_iid_state():
    rho = 0.7345934905942886
    stationary_variance = 0.5084845767213341
    innovation_sd = np.sqrt(stationary_variance * (1.0 - rho**2))
    overrides, meta = _candidate(
        rho_period=rho,
        persistent_innovation_sd_period=innovation_sd,
        transitory_sd_period=0.0,
        n_persistent=7,
        n_iid=1,
    )
    assert overrides["z_grid"].size == 7
    np.testing.assert_array_equal(meta["components"]["iid_levels"], [1.0])
    np.testing.assert_array_equal(meta["components"]["iid_weights"], [1.0])
    np.testing.assert_array_equal(overrides["Pi_z"].shape, (7, 7))
    expected = stationary_variance * rho ** np.arange(5)
    np.testing.assert_allclose(meta["continuous_period_log_covariances"], expected, atol=1e-14)
    np.testing.assert_allclose(meta["discrete_period_log_covariances"], expected, atol=1e-12)
    assert meta["transitory_sd_period"] == 0.0
    assert {row["joint_states"] for row in meta["resolution_table"]} == {5, 9, 15}


@pytest.mark.parametrize(
    "kwargs",
    [
        {"transitory_sd_period": 0.0, "n_iid": 3},
        {"transitory_sd_period": 0.085, "n_iid": 1},
    ],
)
def test_iid_variance_and_iid_dimension_must_collapse_together(kwargs):
    with pytest.raises(ValueError):
        _candidate(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"rho_period": 4.0},
        {"rho_period": np.nan},
        {"persistent_innovation_sd_period": 0.0},
        {"persistent_innovation_sd_period": np.inf},
        {"transitory_sd_period": -0.1},
        {"n_persistent": 1},
        {"n_iid": 3.5},
    ],
)
def test_invalid_units_or_collapsed_variance_are_rejected(kwargs):
    with pytest.raises(ValueError):
        _candidate(**kwargs)
