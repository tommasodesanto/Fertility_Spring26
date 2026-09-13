"""Small pure test for the diagnostic hook and gate classification."""
from pathlib import Path

import numpy as np

import diagnose_e5f_transition_mass as diagnostic


def advance_cohort_one_period_markov_income():
    gj = np.ones((2, 1, 1, 1, 1, 1))
    gpl = gj.copy()
    gpt = gj.copy()
    gps = gj.copy()
    g_next = gj.copy()
    tenure_probs = np.ones((2, 1, 1, 20, 1, 1, 1, 1), dtype=np.float32)
    j = 14
    return g_next


def test_profile_captures_stages(tmp_path: Path):
    hook, state = diagnostic._profile_factory(tmp_path, 1.0e-8)
    import sys
    previous = sys.getprofile()
    sys.setprofile(hook)
    try:
        advance_cohort_one_period_markov_income()
    finally:
        sys.setprofile(previous)
    assert state["calls"] == 1
    assert state["first_failure"] is None
    record = (tmp_path / "transition_calls.jsonl").read_text()
    assert '"gpl"' in record and '"g_next"' in record
    assert '"j"' not in record


def test_gate_exception_classification():
    assert diagnostic._is_expected_gate_error(
        {"type": "RuntimeError", "message": "x mass gate failed: y"}
    )
    assert not diagnostic._is_expected_gate_error(
        {"type": "RuntimeError", "message": "unrelated failure"}
    )
