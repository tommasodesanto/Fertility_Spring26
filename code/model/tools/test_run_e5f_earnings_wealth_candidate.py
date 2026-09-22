from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pytest

import run_e5f_earnings_wealth_candidate as candidate


ROOT = Path(candidate.__file__).resolve().parents[3]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plan(tmp_path: Path, *, arm: str = "reference", repetitions: int = 1) -> dict:
    wrapper = tmp_path / "wrapper.py"
    marker = tmp_path / "preflight.marker"
    wrapper.write_text(
        f"from pathlib import Path\n"
        f"def preflight(*args): Path({str(marker)!r}).write_text('called')\n"
        "def run(*args): raise AssertionError('run called')\n"
    )
    initial = tmp_path / "initial.json"
    initial.write_text('{"source_sha256": {}}')
    run_contract = tmp_path / "run.json"
    run_contract.write_text("{}")
    files = {
        "adapter": {"path": str(Path(candidate.__file__).resolve()), "sha256": _sha(Path(candidate.__file__))},
        "accounting": {"path": str(Path(candidate.accounting.__file__).resolve()), "sha256": _sha(Path(candidate.accounting.__file__))},
        "income": {"path": str(Path(candidate.income.__file__).resolve()), "sha256": _sha(Path(candidate.income.__file__))},
        "wrapper": {"path": str(wrapper), "sha256": _sha(wrapper)},
        "initial_contract": {"path": str(initial), "sha256": _sha(initial)},
        "run_contract": {"path": str(run_contract), "sha256": _sha(run_contract)},
    }
    cases = [{"id": "reference_case", "arm": arm, "repetitions": repetitions,
              "native_seconds": 5, "wrapper_seconds": 5}]
    return {
        "objective_canonical_sha256": candidate.OBJECTIVE,
        "files": files,
        "cases": cases,
        "source_root": str(ROOT),
        "structural_parameters": {"beta_annual": 0.99},
        "initial_psi": 0.1,
        "reference_loss": 0.0,
        "income_specification": {
            "author_decision": "approved_diagnostic",
            "mapping": "conventional_endpoint",
            "constructor_arguments": {},
        },
    }


def test_candidate_refuses_pending_decision_and_unknown_mapping(monkeypatch, tmp_path):
    plan = _plan(tmp_path)
    plan["income_specification"]["author_decision"] = "pending"
    with pytest.raises(ValueError, match="decision is pending"):
        candidate.candidate(plan)
    plan["income_specification"]["author_decision"] = "approved_diagnostic"
    plan["income_specification"]["mapping"] = "unknown"
    with pytest.raises(ValueError, match="unsupported income-period mapping"):
        candidate.candidate(plan)


def test_verify_plan_rejects_file_mismatch_and_loaded_identity(tmp_path):
    plan = _plan(tmp_path)
    plan["files"]["wrapper"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed contracted file"):
        candidate.verify_plan(plan)
    plan = _plan(tmp_path)
    copied = tmp_path / "accounting_copy.py"
    shutil.copyfile(plan["files"]["accounting"]["path"], copied)
    plan["files"]["accounting"] = {"path": str(copied), "sha256": _sha(copied)}
    with pytest.raises(ValueError, match="loaded a different accounting"):
        candidate.verify_plan(plan)


def test_run_case_reference_bypasses_candidate_and_preflight_skips_run(monkeypatch, tmp_path):
    plan = _plan(tmp_path)
    marker = tmp_path / "preflight.marker"
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr(candidate, "candidate", lambda _: (_ for _ in ()).throw(AssertionError("candidate called")))
    out = tmp_path / "case"
    result = candidate.run_case(plan_path, plan, "reference", out, 1, preflight=True)
    assert result["household_solves"] == 0
    assert marker.read_text() == "called"


def test_run_case_pins_repetitions_and_calls_candidate(monkeypatch, tmp_path):
    plan = _plan(tmp_path, arm="literature_income", repetitions=1)
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    calls = []
    monkeypatch.setattr(candidate, "candidate", lambda _: calls.append(True) or ({}, {"components": {}}))
    with pytest.raises(ValueError, match="contracted case"):
        candidate.run_case(plan_path, plan, "literature_income", tmp_path / "case", 2, preflight=True)
    assert calls == []
    candidate.run_case(plan_path, plan, "literature_income", tmp_path / "case", 1, preflight=True)
    assert calls == [True]


def test_grid_gate_measures_error_against_intended_endpoint_process(tmp_path):
    plan = _plan(tmp_path)
    spec = plan['income_specification']
    spec['max_relative_discrete_level_covariance_error'] = .05
    spec['constructor_arguments'] = {'n_persistent': 15}
    overrides, metadata = candidate.candidate(plan)
    assert len(overrides['z_grid']) == 45
    assert metadata['maximum_relative_discrete_level_covariance_error'] < .05
    spec['constructor_arguments'] = {'n_persistent': 5}
    with pytest.raises(ValueError, match='distribution-approximation gate'):
        candidate.candidate(plan)
