from __future__ import annotations

import numpy as np

from diagnose_e5f_income_grid_resolution import (
    _iid_lognormal_rule,
    _weighted_quantile,
    block_average_level_covariance,
)


def test_iid_rule_has_mean_one_and_exact_log_variance():
    levels, weights = _iid_lognormal_rule(5, 0.31)
    np.testing.assert_allclose(weights @ levels, 1.0, atol=1e-14)
    logs = np.log(levels)
    np.testing.assert_allclose(weights @ (logs - weights @ logs) ** 2, 0.31**2, atol=1e-14)


def test_weighted_quantile_is_upper_step_quantile():
    levels = np.array([1.0, 2.0, 5.0])
    weights = np.array([0.2, 0.5, 0.3])
    assert _weighted_quantile(levels, weights, 0.2) == 1.0
    assert _weighted_quantile(levels, weights, 0.21) == 2.0
    assert _weighted_quantile(levels, weights, 0.99) == 5.0


def test_block_formula_iid_limit_is_level_variance_over_four():
    # This checks the block reference used by the prior MC receipt without
    # rerunning any annual panel.
    np.testing.assert_allclose(
        block_average_level_covariance(0.0, 0.4, 0.8, 4, 0),
        (np.exp(0.4) - 1.0) / 4.0,
    )
