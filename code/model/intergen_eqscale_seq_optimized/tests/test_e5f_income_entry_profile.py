from __future__ import annotations

import numpy as np

from intergen_eqscale_seq_optimized.e5f_floor_profile import E5F_DOMAIN
from intergen_eqscale_seq_optimized.e5f_income_entry_profile import (
    E5F_INCOME_ENTRY_DOMAIN,
    e5f_income_entry_metadata,
    e5f_income_entry_overrides,
)
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters


def test_profile_adds_only_one_free_parameter() -> None:
    assert len(E5F_INCOME_ENTRY_DOMAIN) == len(E5F_DOMAIN) + 1
    assert E5F_INCOME_ENTRY_DOMAIN[-1] == (
        "first_birth_fixed_cost",
        0.0,
        8.0,
        "softzero",
    )


def test_profile_reuses_measured_permanent_income_process() -> None:
    parameters = apply_overrides(setup_parameters(), e5f_income_entry_overrides())
    assert parameters.permanent_income_levels_enabled is True
    assert parameters.Nz == 15
    assert parameters.first_birth_fixed_cost == 0.0
    np.testing.assert_allclose(
        parameters.z_weights @ parameters.Pi_z,
        parameters.z_weights,
        atol=1.0e-12,
        rtol=0.0,
    )


def test_profile_metadata_separates_external_and_free_objects() -> None:
    metadata = e5f_income_entry_metadata()
    assert metadata["added_free_parameter"] == "first_birth_fixed_cost"
    assert "no added free" in metadata["permanent_income"]["identification"]
