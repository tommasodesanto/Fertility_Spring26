from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.e6c_profile import (
    E6C_DOMAIN,
    E6C_SEED,
    e6c_metadata,
    e6c_overrides,
)
from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    readiness_cumulative_probability,
    readiness_transition_hazard,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import apply_child_aging, precompute_shared


def readiness_parameters() -> SimpleNamespace:
    return apply_overrides(
        setup_parameters(),
        {
            "sequential_births": True,
            **e6c_overrides(),
            **E6C_SEED,
        },
    )


def test_readiness_profile_has_exactly_two_free_coordinates() -> None:
    assert [name for name, *_ in E6C_DOMAIN] == [
        "readiness_location_age",
        "readiness_spread_years",
    ]
    assert e6c_overrides() == {"readiness_gate_enabled": True}
    assert e6c_metadata()["identifying_moments"] == [
        "mean_age_first_birth",
        "share_first_births_age30plus",
    ]


def test_logistic_readiness_arrival_is_irreversible_and_consistent() -> None:
    parameters = readiness_parameters()
    ages = np.array([18.0, 22.0, 26.0, 30.0, 34.0])
    cumulative = np.array(
        [readiness_cumulative_probability(parameters, age) for age in ages]
    )

    assert np.all(np.diff(cumulative) > 0.0)
    for current, nxt, f_current, f_next in zip(
        ages[:-1], ages[1:], cumulative[:-1], cumulative[1:]
    ):
        hazard = readiness_transition_hazard(parameters, current, nxt)
        assert 0.0 < hazard < 1.0
        assert f_current + (1.0 - f_current) * hazard == pytest.approx(f_next)


def test_readiness_gate_requires_sequential_births() -> None:
    with pytest.raises(ValueError, match="sequential_births"):
        apply_overrides(setup_parameters(), e6c_overrides())


def test_settled_childless_state_has_childless_flow_preferences() -> None:
    parameters = readiness_parameters()
    b_grid = np.array([-1.0, 0.0, 1.0])
    shared = precompute_shared(parameters, b_grid)

    assert shared.cb_flat[0, parameters.n_parity] == shared.cb_flat[0, 0]
    assert shared.hb_flat[0, parameters.n_parity] == shared.hb_flat[0, 0]
    assert shared.psi_flat[0, parameters.n_parity] == shared.psi_flat[0, 0]


def test_backward_readiness_transition_uses_unsettled_and_settled_values() -> None:
    parameters = readiness_parameters()
    shape = (1, 1, 1, parameters.n_parity, parameters.n_child_states)
    next_value = np.zeros(shape)
    next_value[..., 0, 0] = 2.0
    next_value[..., 0, 1] = 10.0
    current = apply_child_aging(
        next_value,
        parameters,
        1,
        1,
        1,
        parameters.n_parity,
        parameters.n_child_states,
        age_index=0,
    )
    hazard = readiness_transition_hazard(parameters, 18.0, 22.0)

    assert current[..., 0, 0].item() == pytest.approx(
        (1.0 - hazard) * 2.0 + hazard * 10.0
    )
    assert current[..., 0, 1].item() == 10.0


def test_default_off_precompute_is_bitwise_identical() -> None:
    baseline = setup_parameters()
    explicit_off = apply_overrides(
        setup_parameters(),
        {"readiness_gate_enabled": False},
    )
    b_grid = np.array([-1.0, 0.0, 1.0])
    first = precompute_shared(baseline, b_grid)
    second = precompute_shared(explicit_off, b_grid)

    for name in ("cb_flat", "hb_flat", "psi_flat", "gb_flat", "alpha_flat", "escale_flat"):
        assert np.array_equal(getattr(first, name), getattr(second, name))


def test_winner_diagnostic_configures_combined_e6c_arm() -> None:
    from intergen_eqscale_seq_optimized import build_e6_winner_diagnostics

    with patch.dict("os.environ", {}, clear=True):
        build_e6_winner_diagnostics.configure_environment("e6abc")
        assert build_e6_winner_diagnostics.expected_arm("e6abc") == "E6ABC"
        assert all(
            os.environ[name] == "1"
            for name in ("E6A", "E6B", "E6C")
        )


def test_strict_collector_accepts_twelve_parameter_e6c_contract() -> None:
    from intergen_eqscale_seq_optimized import collect_e1
    from intergen_eqscale_seq_optimized.e5_profile import E5_DOMAIN

    domain = E5_DOMAIN + E6C_DOMAIN
    theta = {
        ("beta" if name == "beta_annual" else name): (
            0.96**4 if name == "beta_annual" else (lower + upper) / 2.0
        )
        for name, lower, upper, _transform in domain
    }
    tight = {
        "strict_converged": True,
        "rank_loss": 1.0,
        "market_residual": 1e-6,
        "target_fit": [],
        "theta": theta,
    }
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory) / "chains"
        chain = root / "chain_1"
        chain.mkdir(parents=True)
        summary = {
            "metadata": {
                "arm": "E6ABC",
                "free_parameter_count": 12,
                "target_count": 12,
                "seed": 1,
                "active_domain": [
                    {
                        "name": name,
                        "lower": lower,
                        "upper": upper,
                        "transform": transform,
                    }
                    for name, lower, upper, transform in domain
                ],
            },
            "best_tight": tight,
            "tight_repeat_check": {
                "both_strict": True,
                "loss_abs_difference": 0.0,
                "max_abs_moment_difference": 0.0,
            },
            "n_cases_completed": 1,
            "n_strict": 1,
        }
        (chain / "summary.json").write_text(json.dumps(summary))
        outdir = Path(temporary_directory) / "report"
        with patch(
            "sys.argv",
            [
                "collect_e1.py",
                "--results-root",
                str(root),
                "--outdir",
                str(outdir),
            ],
        ):
            collect_e1.main()

        parameter_rows = (outdir / "parameter_table_full.csv").read_text()
        assert parameter_rows.count("\n") == 13
        assert "readiness_location_age" in parameter_rows
        assert "readiness_spread_years" in parameter_rows
