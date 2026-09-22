from __future__ import annotations

import json
import random
import copy
from pathlib import Path

import pytest

import run_e5f_earnings_wealth_search as search


def _plan(tmp_path: Path) -> dict:
    names = search.PARAMETERS
    start = {n: 1.0 for n in names}
    start.update(beta_annual=.99, h_P=2.3)
    bounds = {n: [.1, 2.] for n in names}
    bounds["beta_annual"] = [.94, .99]
    bounds["h_P"] = [.1, 2.3]
    return {
        "objective_canonical_sha256": search.OBJECTIVE,
        "source_root": str(tmp_path),
        "files": {"adapter": {"path": str(tmp_path / "adapter.py"), "sha256": "x"}, "initial_contract": {"path": str(tmp_path / "initial.json"), "sha256": "pinned"}},
        "parameter_bounds": bounds,
        "starting_structural_parameters": start,
        "structural_parameters": start,
        "initial_psi": .15,
        "search_config": {"arm": "reference", "max_proposals": 4, "workers": 2, "search_seconds": 3, "verification_seconds": 1, "case_timeout": 1, "seed": 5},
    }


def _receipt(point: dict[str, float], loss: float, reps: int = 1) -> dict:
    params = [{"parameter": n, "estimate": v, "structural_coordinate": True} for n, v in point.items()]
    params += [{"parameter": f"fixed_{i}", "estimate": 0.0} for i in range(17 - len(params))]
    return {"schema": "e5f_initial_minimum_distance_result_v1", "loss": loss, "target_fit": [{"restriction_id": str(i), "target": 0., "model": 0., "gap": 0., "loss_contribution": 0.} for i in range(13)], "parameters": params, "_summary": {"status": "verified_scored_candidate", "objective_canonical_sha256": search.OBJECTIVE, "repetitions": reps, "exact_loss_equality": reps == 2, "original_graphs": [{} for _ in range(17)]}}


def test_plan_and_proposals_are_nine_coordinate_and_bounded(tmp_path):
    plan = _plan(tmp_path)
    search.validate_plan(plan)
    points = search.proposal_batch(plan["starting_structural_parameters"], plan, random.Random(2), 4, {search.fingerprint(plan["starting_structural_parameters"])}, 1.)
    assert len(points) == 4
    assert all(set(p) == set(search.PARAMETERS) for p in points)
    assert all(plan["parameter_bounds"][n][0] <= p[n] <= plan["parameter_bounds"][n][1] for p in points for n in search.PARAMETERS)


def test_gate_rejection_is_reviewable_and_contract_is_fatal():
    assert search.classify_failure(RuntimeError("market clearing gate failed")) == "numerical_gate_rejection"
    assert search.classify_failure(RuntimeError("source fingerprint mismatch")) == "fatal_contract_error"


def test_bounded_search_keeps_cases_and_exact_verification(tmp_path):
    plan = _plan(tmp_path)
    start = plan["starting_structural_parameters"]
    anchor = _receipt(start, 10., reps=2)
    calls = []

    def fake(point, case, **kwargs):
        calls.append((point, kwargs["repetitions"]))
        case.mkdir(parents=True, exist_ok=True)
        if kwargs["repetitions"] == 2:
            return _receipt(point, 1., reps=2)
        return _receipt(point, 1., reps=1)

    result = search.run_search(plan, tmp_path / "search", fake, anchor=anchor, search_seconds=4)
    assert result["status"] == "verified_selection"
    assert any(reps == 2 for _, reps in calls)
    assert (tmp_path / "search" / "cases.jsonl").exists()
    assert (tmp_path / "search" / "heartbeat.json").exists()


def test_receipt_requires_complete_native_tables(tmp_path):
    point = _plan(tmp_path)["starting_structural_parameters"]
    bad = _receipt(point, 1.)
    bad["parameters"] = bad["parameters"][:-1]
    with pytest.raises(search.ContractError, match="13 targets and 17 parameters"):
        search._validate_receipt(bad, point, repetitions=1)


def test_numeric_replay_excludes_only_wall_clock_time(tmp_path):
    a = _receipt(_plan(tmp_path)["starting_structural_parameters"], 1.)
    a["normalization"] = dict(psi_child=.2, completed_fertility=2.1,
        target=2.1, absolute_gap=0., stationary_solves=6, stationary_solve_seconds=100.)
    b = copy.deepcopy(a)
    b["normalization"]["stationary_solve_seconds"] = 200.
    assert search._numeric_fit(a) == search._numeric_fit(b)
    for key in ("psi_child", "completed_fertility", "stationary_solves"):
        changed = copy.deepcopy(b)
        changed["normalization"][key] += .001
        assert search._numeric_fit(a) != search._numeric_fit(changed)


def test_unexpected_failure_is_fatal_and_all_gate_failures_stop(tmp_path):
    plan = _plan(tmp_path)
    anchor = _receipt(plan["starting_structural_parameters"], 10., reps=2)

    def fatal(point, case, **kwargs):
        raise RuntimeError("unexpected coding failure")

    with pytest.raises(search.ContractError, match="unexpected coding failure"):
        search.run_search(plan, tmp_path / "fatal", fatal, anchor=anchor, search_seconds=4)

    def gates(point, case, **kwargs):
        raise RuntimeError("market clearing gate failed")

    with pytest.raises(search.ContractError, match="all dispatched cases"):
        search.run_search(plan, tmp_path / "gates", gates, anchor=anchor, search_seconds=4)
