import numpy as np
import pytest
from pathlib import Path

from intergen_eqscale_seq_optimized.build_e2_packet import (
    DEFAULT_OUTDIR,
    E4_DEFAULT_OUTDIR,
    E4_SOURCE,
    SOURCE,
    arm_externals,
    load_winner_theta,
    parse_args,
)
from intergen_eqscale_seq_optimized.externals import flhsv_income_overrides


E4_RESULTS = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/"
    "output/model/eqscale_seq_e4_split_recalibration_20260723/report/results.json"
)


def test_e4_winner_keys_are_accepted_only_by_e4_arm():
    theta = load_winner_theta(E4_RESULTS, "e4")
    assert "kappa_fert_continuation" in theta
    assert "gamma_e" not in theta
    with pytest.raises(ValueError, match="unexpected E2 winner keys"):
        load_winner_theta(E4_RESULTS, "e2")


def test_e4_arm_externals_use_power_scale_and_flhsv_income():
    assert arm_externals("e2") == {}
    overrides = arm_externals("e4")
    assert overrides["eqscale_form"] == "power"
    np.testing.assert_array_equal(overrides["Pi_z"], flhsv_income_overrides()["Pi_z"])


def test_e4_postprocesses_defaults_without_changing_e2_defaults(monkeypatch):
    monkeypatch.setattr("sys.argv", ["build_e2_packet.py"])
    e2_args = parse_args()
    assert e2_args.arm == "e2"
    assert e2_args.source == SOURCE
    assert e2_args.outdir == DEFAULT_OUTDIR

    monkeypatch.setattr("sys.argv", ["build_e2_packet.py", "--arm", "e4"])
    e4_args = parse_args()
    assert e4_args.source == E4_SOURCE
    assert e4_args.outdir == E4_DEFAULT_OUTDIR
