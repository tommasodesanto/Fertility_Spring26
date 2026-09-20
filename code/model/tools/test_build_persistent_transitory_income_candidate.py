from __future__ import annotations

import numpy as np

from build_persistent_transitory_income_candidate import (
    build_persistent_transitory_income_candidate,
)


def _candidate(sd: float = 0.16):
    return build_persistent_transitory_income_candidate(
        rho_annual=0.964,
        persistent_innovation_sd_annual=0.150,
        transitory_log_sd_period=sd,
        period_years=4.0,
    )


def test_joint_rule_has_valid_rows_stationarity_and_mean_one():
    overrides, meta = _candidate()
    z, w, pi = overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"]
    np.testing.assert_allclose(w.sum(), 1.0)
    np.testing.assert_allclose(w @ z, 1.0)
    np.testing.assert_allclose(pi.sum(axis=1), 1.0)
    np.testing.assert_allclose(w @ pi, w, atol=1e-12)
    assert z.size == 15 and meta["joint_states"] == 15
    assert overrides["permanent_income_levels_enabled"] is False


def test_iid_transition_is_independent_of_current_transitory_state():
    overrides, _ = _candidate()
    pi = overrides["Pi_z"].reshape(5, 3, 5, 3)
    for p0 in range(5):
        for p1 in range(5):
            np.testing.assert_allclose(pi[p0, 0, p1, :], pi[p0, 2, p1, :])


def test_zero_transitory_risk_preserves_persistent_expectations():
    overrides, _ = _candidate(0.0)
    z = overrides["z_grid"].reshape(5, 3)
    w = overrides["z_weights"].reshape(5, 3)
    np.testing.assert_allclose(z[:, 0], z[:, 1])
    np.testing.assert_allclose(z[:, 1], z[:, 2])
    np.testing.assert_allclose(w.sum(axis=1), np.array([.0625, .25, .375, .25, .0625]))
    np.testing.assert_allclose((w * z).sum(axis=1), z[:, 0] * np.array([.0625, .25, .375, .25, .0625]))


def test_log_variance_is_explicit_and_period_is_frozen():
    from build_persistent_transitory_income_candidate import _iid_lognormal_rule
    z, w = _iid_lognormal_rule(0.2)
    logz = np.log(z)
    np.testing.assert_allclose(w @ ((logz - w @ logz) ** 2), 0.04, atol=1e-14)
    _, meta = _candidate(0.2)
    assert np.isclose(meta["transitory_log_variance_period"], 0.04)
    import pytest
    with pytest.raises(ValueError):
        build_persistent_transitory_income_candidate(
            rho_annual=.964, persistent_innovation_sd_annual=.15,
            transitory_log_sd_period=.2, period_years=1.0,
        )


def test_economic_log_covariances_match_persistent_plus_iid_process():
    rho, vp, ve = .9703033301891563, .692826827789733, .09785297241466638
    overrides, _ = build_persistent_transitory_income_candidate(
        rho_annual=rho, persistent_innovation_sd_annual=np.sqrt((1-rho**2)*vp),
        transitory_log_sd_period=np.sqrt(ve), period_years=4.0)
    w, pi = overrides['z_weights'], overrides['Pi_z']
    y = np.log(overrides['z_grid']); y -= w @ y
    np.testing.assert_allclose(w @ y**2, vp + ve, atol=1e-12)
    np.testing.assert_allclose((w*y) @ pi @ y, vp*rho**4, atol=1e-12)
    np.testing.assert_allclose((w*y) @ pi @ pi @ y, vp*rho**8, atol=1e-12)
