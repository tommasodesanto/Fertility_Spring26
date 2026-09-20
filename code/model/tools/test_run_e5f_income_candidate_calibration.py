from pathlib import Path
import json
import pytest

from run_e5f_income_candidate_calibration import OBJECTIVE, assemble_pilot_contract, validate_plan


ROOT = Path(__file__).resolve().parents[3]
PLAN = ROOT / "output/model/native_financing_diagnostic_20260919/earnings_candidate/calibration_plan.json"


def local_plan():
    plan = json.loads(PLAN.read_text())
    plan['source_manifest_path'] = str(ROOT / 'output/model/paper_baseline_sep14/native/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_template_v6/inputs/working_contract.json')
    return plan


def test_plan_preserves_objective_and_nine_coordinates_without_solving():
    plan = local_plan()
    assert plan["objective_canonical_sha256"] == OBJECTIVE
    assert plan["free_parameter_count"] == 9
    manifest = json.loads(Path(plan["source_manifest_path"]).read_text())
    assert len(manifest["source_provenance"]["observation_snapshot"]["source_sha256"]) == 641
    candidate = validate_plan(plan, require_source=False)
    assert candidate["state_count"] == 15


def test_source_gate_is_fail_closed(tmp_path):
    plan = local_plan()
    plan["source_root"] = str(tmp_path)
    with pytest.raises(RuntimeError, match="BLOCKED"):
        validate_plan(plan, require_source=True)


def test_pilot_contract_is_cold_and_keeps_frozen_objective():
    plan = local_plan()
    candidate = validate_plan(plan, require_source=False)
    pilot = assemble_pilot_contract(plan, candidate)
    assert pilot["status"] == "assembled_review_required"
    assert pilot["objective_canonical_sha256"] == OBJECTIVE
    assert pilot["cold_solve_required"] is True
    assert pilot["cached_checkpoint_policy"].startswith("reference-only")
