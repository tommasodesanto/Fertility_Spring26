from __future__ import annotations

import numpy as np
import pytest

from build_literature_period_income import (
    annual_block_level_covariances,
    build_literature_period_income_adapter,
)


def test_annual_iid_only_block_has_analytic_variance_and_zero_lags():
    sd = 0.17
    got = annual_block_level_covariances(0.95, 1e-14, sd, period_years=4, max_lag=4)
    expected = np.expm1(sd**2) / 4.0
    np.testing.assert_allclose(got[0], expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(got[1:], 0.0, atol=1e-12)


def test_block_covariance_sum_identity():
    rho, eta, eps = 0.95, 0.21, 0.17
    got = annual_block_level_covariances(rho, eta, eps)
    vp = eta**2 / (1 - rho**2)
    direct = sum(np.expm1(vp * rho**abs(j - i) + (eps**2 if i == j else 0.0)) for i in range(4) for j in range(4)) / 16
    assert np.isclose(got[0], direct)


def test_proxy_matches_c0_c1_c2_and_records_finite_lag4():
    overrides, meta = build_literature_period_income_adapter(transitory_sd_annual=0.8)
    target = np.asarray(meta["continuous_covariance_targets"])
    proxy = np.asarray(meta["proxy_fit"]["proxy_covariances"])
    np.testing.assert_allclose(proxy[:3], target[:3], atol=2e-14)
    assert np.isfinite(proxy[4] - target[4])
    assert 0.0 < meta["proxy_fit"]["rho_period"] < 1.0
    # log covariance metadata is the autocovariance of log nodes, not log1p(level covariance).
    assert not np.allclose(meta["actual_chain_log_covariances"], np.log1p(meta["actual_chain_level_covariances"]))
    assert overrides["permanent_income_levels_enabled"] is False


def test_chain_is_stochastic_stationary_and_iid_rows_are_independent():
    overrides, meta = build_literature_period_income_adapter(transitory_sd_annual=0.8)
    z, w, pi = overrides["z_grid"], overrides["z_weights"], overrides["Pi_z"]
    assert z.size == 15
    np.testing.assert_allclose(pi.sum(axis=1), 1.0, atol=1e-14)
    np.testing.assert_allclose(w @ pi, w, atol=1e-12)
    np.testing.assert_allclose(w @ z, 1.0, atol=1e-14)
    reshaped = pi.reshape(5, 3, 5, 3)
    for p0 in range(5):
        for p1 in range(5):
            np.testing.assert_allclose(reshaped[p0, 0, p1, :], reshaped[p0, 2, p1, :])
    assert len(meta["grid_resolution_table"]) == 3
    assert [row["n_persistent"] for row in meta["grid_resolution_table"]] == [5, 9, 15]


def test_inadmissible_proxy_is_rejected_without_clipping(monkeypatch):
    import build_literature_period_income as module

    monkeypatch.setattr(module, "annual_block_level_covariances", lambda *args, **kwargs: np.array([0.1, 0.2, 0.01, 0.0, 0.0]))
    with pytest.raises(ValueError, match="inadmissible period proxy"):
        module.build_literature_period_income_adapter()


def test_sommer_pins_report_explicit_negative_transitory_variance_blocker():
    import build_literature_period_income as module
    diagnostic = module.diagnose_literature_period_income_adapter()
    assert diagnostic["status"] == "blocked_inadmissible_proxy"
    assert diagnostic["proxy_fit"]["Ve"] < 0.0
    with pytest.raises(ValueError, match="inadmissible period proxy"):
        module.build_literature_period_income_adapter()
    endpoint = diagnostic["conventional_endpoint"]
    assert endpoint["label"].endswith("diagnostic_only")
    assert len(endpoint["exact_annual_block_covariances"]) == 5
    assert [x["n_persistent"] for x in endpoint["grid_resolution_table"]] == [5, 9, 15]
    assert all(len(x["errors_vs_exact_block"]) == 5 for x in endpoint["grid_resolution_table"])


def test_endpoint_constructor_is_deterministic_and_separates_error_benchmarks():
    import build_literature_period_income as module
    a, ma = module.build_conventional_endpoint_income_adapter(n_persistent=5)
    b, mb = module.build_conventional_endpoint_income_adapter(n_persistent=5)
    np.testing.assert_array_equal(a["z_grid"], b["z_grid"])
    np.testing.assert_array_equal(a["Pi_z"], b["Pi_z"])
    assert len(ma["discrete_level_errors_vs_continuous_endpoint"]) == 5
    assert len(ma["discrete_level_errors_vs_exact_block"]) == 5
    assert ma["iid_transition_independent"] and ma["stationary"]
    assert ma == mb
    with pytest.raises(ValueError, match="period_years"):
        module.build_conventional_endpoint_income_adapter(period_years=4.5)
