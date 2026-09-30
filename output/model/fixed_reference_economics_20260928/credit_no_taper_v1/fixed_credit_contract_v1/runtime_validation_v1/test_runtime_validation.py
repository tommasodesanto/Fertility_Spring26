import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def test_plan_preserves_bounded_contract():
    p = json.loads((ROOT / "plan.json").read_text())
    assert p["planned_lifecycle_solves"] == 1 and p["maximum_lifecycle_solves"] == 2
    assert p["case_seconds"] == 300 and p["total_seconds"] == 900
    assert p["unsecured_credit_limit"] == 0.0 and p["cases"] == ["control"]

def test_cash_diagnosis_blocks_strict_solve():
    d = json.loads((ROOT / "affordability_diagnosis.json").read_text())
    assert d["status"] == "strict_zero_lifecycle_blocked_before_submission"
    assert all(x < 0.0 for x in d["cash_before_consumption_rent_saving"].values())
    assert all(x < 0.0 for x in d["purchase_cash_before_downpayment"].values())
