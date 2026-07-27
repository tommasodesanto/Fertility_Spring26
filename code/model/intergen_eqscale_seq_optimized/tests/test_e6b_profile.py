from __future__ import annotations

import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.e6b_profile import (
    E6B_FIXED_LOG_VARIANCE,
    e6b_metadata,
    e6b_overrides,
    permanent_income_levels,
)
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters


def test_permanent_level_rule_matches_mean_and_log_variance() -> None:
    levels, weights, log_nodes = permanent_income_levels()

    assert levels.shape == weights.shape == log_nodes.shape == (3,)
    assert weights == pytest.approx([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0])
    assert np.all(np.diff(levels) > 0.0)
    assert float(weights @ levels) == pytest.approx(1.0, abs=1e-14)
    log_mean = float(weights @ log_nodes)
    log_variance = float(weights @ (log_nodes - log_mean) ** 2)
    assert log_variance == pytest.approx(E6B_FIXED_LOG_VARIANCE, abs=1e-14)


def test_combined_process_keeps_permanent_levels_fixed() -> None:
    overrides = e6b_overrides()
    grid = np.asarray(overrides["z_grid"])
    weights = np.asarray(overrides["z_weights"])
    transition = np.asarray(overrides["Pi_z"])
    groups = np.asarray(overrides["permanent_income_group_index"])

    assert grid.shape == weights.shape == groups.shape == (15,)
    assert transition.shape == (15, 15)
    assert float(weights @ grid) == pytest.approx(1.0, abs=1e-14)
    assert weights @ transition == pytest.approx(weights, abs=1e-14)
    assert transition.sum(axis=1) == pytest.approx(np.ones(15), abs=1e-14)
    for row in range(15):
        outside_group = groups != groups[row]
        assert np.all(transition[row, outside_group] == 0.0)


def test_e6b_overrides_pass_parameter_validation() -> None:
    parameters = apply_overrides(setup_parameters(), e6b_overrides())

    assert parameters.permanent_income_levels_enabled is True
    assert parameters.Nz == 15
    assert parameters.permanent_income_group_index.shape == (15,)
    assert parameters.permanent_income_base_state_index.shape == (15,)


def test_e6b_metadata_discloses_external_identification() -> None:
    metadata = e6b_metadata()

    assert metadata["combined_income_state_count"] == 15
    assert metadata["fixed_log_variance"] == pytest.approx(E6B_FIXED_LOG_VARIANCE)
    assert "no added free" in metadata["identification"]


def test_default_parameters_leave_permanent_levels_disabled() -> None:
    parameters = setup_parameters()

    assert parameters.permanent_income_levels_enabled is False
    assert parameters.permanent_income_log_variance == 0.0


@pytest.mark.parametrize("arm", ["E6A", "E6B", "E6AB"])
def test_strict_collector_accepts_e6_contracts(arm: str) -> None:
    from intergen_eqscale_seq_optimized import collect_e1

    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory) / "chains"
        chain = root / "chain_1"
        chain.mkdir(parents=True)
        tight = {
            "strict_converged": True,
            "rank_loss": 1.0,
            "target_fit": [],
        }
        summary = {
            "metadata": {
                "arm": arm,
                "free_parameter_count": 10,
                "target_count": 12,
                "seed": 1,
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
        assert (outdir / "results.json").exists()
