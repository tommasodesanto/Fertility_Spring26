"""Validate the saved lifecycle CSV and suppress non-stayer-aware raw fields."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = list(csv.DictReader((HERE / "lifecycle.csv").open(newline="")))
RECEIPT = json.loads((HERE / "receipt.json").read_text())
assert RECEIPT["status"] == "passed_saved_policy_lifecycle" and RECEIPT["model_solves"] == 0
A = "reference_policy_reference_PRE"
B = "credit_policy_reference_PRE"
C = "credit_policy_credit_PRE"
INVALID = (
    "mean_chosen_next_financial_assets",
    "mean_next_minus_post_transaction_financial_assets",
    "mean_nonhousing_consumption",
    "owner_next_negative_asset_mass",
    "owner_next_negative_asset_rate",
    "owner_next_at_grid_min_mass",
)
safe = [dict(row) for row in RAW]
for row in safe:
    if row["scenario"] == A and row["inherited_tenure"] == "owner":
        for key in INVALID:
            row[key] = ""
with (HERE / "lifecycle_validated.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(safe[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(safe)

summary = []
for age in range(18, 43, 4):
    d = {}
    for scenario in (A, B, C):
        rows = [r for r in safe if r["scenario"] == scenario and int(r["age"]) == age]
        d[scenario] = {
            "mass": sum(float(r["origin_never_parent_mass"]) for r in rows),
            "births": sum(float(r["first_birth_flow"]) for r in rows),
            "inherited_renter_mass": sum(float(r["origin_never_parent_mass"]) for r in rows if r["inherited_tenure"] == "renter"),
            "inherited_owner_mass": sum(float(r["origin_never_parent_mass"]) for r in rows if r["inherited_tenure"] == "owner"),
        }
    assert abs(d[A]["mass"] - d[B]["mass"]) < 1e-12
    summary.append({
        "age": age,
        "reference_PRE_never_parent_mass": d[A]["mass"],
        "credit_PRE_never_parent_mass": d[C]["mass"],
        "reference_PRE_reference_policy_first_births": d[A]["births"],
        "reference_PRE_credit_policy_first_births": d[B]["births"],
        "credit_PRE_credit_policy_first_births": d[C]["births"],
        "fixed_PRE_policy_effect": d[B]["births"] - d[A]["births"],
        "fixed_credit_policy_distribution_effect": d[C]["births"] - d[B]["births"],
        "reference_PRE_inherited_owner_mass": d[B]["inherited_owner_mass"],
        "credit_PRE_inherited_owner_mass": d[C]["inherited_owner_mass"],
    })
assert len(safe) == 39
for scenario in (A, B, C):
    actual = sum(float(r["first_birth_flow"]) for r in safe if r["scenario"] == scenario)
    assert abs(actual - RECEIPT["expected_first_births"][scenario]) < 2e-10
assert abs(sum(r["fixed_credit_policy_distribution_effect"] for r in summary) + 0.00540370586708255) < 2e-10
with (HERE / "age_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary)
(HERE / "postprocess_checks.json").write_text(json.dumps({
    "status": "passed",
    "source_rows": len(RAW),
    "validated_rows": len(safe),
    "suppressed_reference_owner_columns": list(INVALID),
    "first_birth_distribution_effect": sum(r["fixed_credit_policy_distribution_effect"] for r in summary),
    "model_solves": 0,
}, indent=2) + "\n")
