import numpy as np

from intergen_eqscale_seq_optimized.externals import (
    FLHSV_ANNUAL_INNOVATION_SD,
    FL_INNOVATION_VARIANCE_ANNUAL,
    FL_RHO_ANNUAL,
    HSV_TAU_US,
    flhsv_income_overrides,
)
from intergen_eqscale_seq_optimized.local_panel import income_process_overrides
from intergen_eqscale_seq_optimized.run_fertility_frontier_scan import cell_overrides


def test_flhsv_external_literals_and_income_overrides():
    assert FL_RHO_ANNUAL == 0.9136
    assert FL_INNOVATION_VARIANCE_ANNUAL == 0.0426
    assert HSV_TAU_US == 0.181
    assert FLHSV_ANNUAL_INNOVATION_SD == (0.0426 ** 0.5) * 0.819

    actual = flhsv_income_overrides(5)
    expected = income_process_overrides(5, "rouwenhorst", (0.0426 ** 0.5) * 0.819, 0.9136)
    assert actual.keys() == expected.keys()
    for key in actual:
        if isinstance(actual[key], np.ndarray):
            np.testing.assert_array_equal(actual[key], expected[key])
        else:
            assert actual[key] == expected[key]


def test_v2_cell_one_uses_power_scale_and_flhsv_income():
    overrides = cell_overrides(1, "v2_imposed_scale")
    assert overrides["eqscale_form"] == "power"
    assert overrides["gamma_e"] == 0.0
    np.testing.assert_array_equal(overrides["Pi_z"], flhsv_income_overrides()["Pi_z"])


def test_v3_entry_noise_cell_overrides_and_prior_designs_are_unchanged():
    cell1 = cell_overrides(1, "v3_entry_noise")
    assert cell1["psi_child"] == 0.5
    assert cell1["delta_alpha_jump"] == 0.087
    assert cell1["kappa_fert"] == 15.9109
    assert "kappa_fert_continuation" not in cell1 or cell1["kappa_fert_continuation"] is None

    cell50 = cell_overrides(50, "v3_entry_noise")
    assert cell50["kappa_fert"] == 3.0
    assert cell50["kappa_fert_continuation"] == 0.3

    cell99 = cell_overrides(99, "v3_entry_noise")
    assert cell99["kappa_fert"] == 10.0
    assert cell99["kappa_fert_continuation"] == 0.3

    v1 = cell_overrides(1, "v1")
    assert v1["psi_child"] == -0.5
    assert v1["kappa_fert"] == 0.5
    assert v1["gamma_e"] == 0.0085
    v2 = cell_overrides(1, "v2_imposed_scale")
    assert v2["eqscale_form"] == "power"
    assert v2["gamma_e"] == 0.0
