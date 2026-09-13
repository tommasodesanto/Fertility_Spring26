from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code/model"))

from intergen_eqscale_seq_optimized.solver import (  # noqa: E402
    advance_cohort_one_period_markov_income,
)


def fixture() -> dict:
    nb, nt, locations, ages, incomes, families, child_states = 2, 3, 1, 2, 1, 1, 1
    b_grid = np.array([0.0, 1.0])
    P = SimpleNamespace(
        n_house=nt - 1,
        I=locations,
        J=ages,
        n_parity=families,
        n_child_states=child_states,
        n_child_stages=1,
        use_numba_scatter=False,
    )
    SD = SimpleNamespace(nc=families * child_states)
    gj = np.zeros((nb, nt, locations, incomes, families, child_states))
    gj[0, 0, 0, 0, 0, 0] = 0.25
    gj[1, 0, 0, 0, 0, 0] = 0.75
    loc_probs = np.ones((nb, nt, locations, locations, ages, incomes, families, child_states))
    tenure_choice = np.zeros((nb, nt, locations, ages, incomes, families, child_states), dtype=np.int16)
    tenure_probs = np.zeros(
        (nb, nt, locations, ages, incomes, families, child_states, nt), dtype=np.float32
    )
    tenure_probs[...] = np.array([0.2, 0.3, 0.5], dtype=np.float32)
    bp_pol = np.broadcast_to(
        b_grid.reshape(nb, 1, 1, 1, 1, 1, 1),
        (nb, nt, locations, ages, incomes, families, child_states),
    ).copy()
    lmm_idx = np.zeros((locations, nt, nb), dtype=np.int64)
    lmm_wt = np.broadcast_to(np.array([0.0, 1.0]), (locations, nt, nb)).copy()
    tmx_idx = np.zeros((locations, nt, nt, families, child_states, nb), dtype=np.int64)
    tmx_wt = np.broadcast_to(
        np.array([0.0, 1.0]), (locations, nt, nt, families, child_states, nb)
    ).copy()
    return dict(
        gj=gj,
        j=0,
        loc_probs=loc_probs,
        tenure_choice=tenure_choice,
        tenure_probs=tenure_probs,
        bp_pol=bp_pol,
        P=P,
        b_grid=b_grid,
        SD=SD,
        lmm_idx=lmm_idx,
        lmm_wt=lmm_wt,
        tmx_idx=tmx_idx,
        tmx_wt=tmx_wt,
        ust=False,
        Pia=None,
        Pi_z=np.ones((incomes, incomes)),
    )


class TenureProbabilityMassConservationTests(unittest.TestCase):
    def test_float32_tenure_rows_are_normalized_once_before_scatter(self) -> None:
        inputs = fixture()
        raw = np.array([0.2, 0.3, 0.5], dtype=np.float32).astype(float)
        self.assertEqual(float(np.sum(raw)) - 1.0, 1.4901161193847656e-08)

        result = advance_cohort_one_period_markov_income(**inputs)
        normalized = raw / np.sum(raw)
        manual = np.zeros_like(result)
        for b, source_mass in enumerate((0.25, 0.75)):
            manual[b, :, 0, 0, 0, 0] = source_mass * normalized

        self.assertEqual(float(np.sum(result)), float(np.sum(inputs["gj"])))
        np.testing.assert_array_equal(result, manual)

    def test_zero_probability_row_remains_zero_for_existing_mass_gate(self) -> None:
        inputs = fixture()
        inputs["tenure_probs"].fill(0.0)
        result = advance_cohort_one_period_markov_income(**inputs)
        self.assertEqual(float(np.sum(result)), 0.0)

    def test_deterministic_tenure_path_is_unchanged(self) -> None:
        inputs = fixture()
        inputs["tenure_probs"] = None
        inputs["tenure_choice"][..., 0] = 2
        result = advance_cohort_one_period_markov_income(**inputs)
        manual = np.zeros_like(result)
        manual[0, 2, 0, 0, 0, 0] = 0.25
        manual[1, 2, 0, 0, 0, 0] = 0.75
        self.assertEqual(float(np.sum(result)), 1.0)
        np.testing.assert_array_equal(result, manual)


if __name__ == "__main__":
    unittest.main()
