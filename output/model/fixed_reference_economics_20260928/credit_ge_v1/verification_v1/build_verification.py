import csv
import hashlib
import json
from pathlib import Path


ROOT = Path("output/model/fixed_reference_economics_20260928")
GE = ROOT / "credit_ge_v1" / "solve_v1"
VERIFY = ROOT / "credit_ge_v1" / "verification_v1"


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observer_values(path):
    data = load(path)
    fertility = data["fertility"]["constant_post_cell"]["moments"]
    housing = data["housing_wealth"]["moments"]
    return {
        "mean_age_first_birth": fertility["period_mean_age_first_birth"],
        "mean_occupied_rooms_uncapped": housing["aggregate_mean_occupied_rooms_ahs_uncapped_18_85"],
        "mean_occupied_rooms_capped9": housing["aggregate_mean_occupied_rooms_capped9_18_85"],
    }


def table_count(path):
    with Path(path).open(newline="") as handle:
        return len(list(csv.DictReader(handle)))


def main():
    remote_manifest = VERIFY / "remote_hash_manifest.sha256"
    manifests = []
    mismatches = []
    for line in remote_manifest.read_text().splitlines():
        expected, remote = line.split(maxsplit=1)
        remote = remote.strip()
        if remote.startswith("sources_solve_v1/"):
            local = VERIFY / "source_receipts" / Path(remote).name
        else:
            local = GE / remote
        record = {"remote_path": remote, "local_path": str(local), "remote_sha256": expected}
        if local.exists():
            record["local_sha256"] = sha(local)
            record["match"] = record["local_sha256"] == expected
        else:
            record["local_sha256"] = None
            record["match"] = False
        manifests.append(record)
        if not record["match"]:
            mismatches.append(record)

    completed = load(GE / "completed.json")
    selected = load(GE / "selected.json")
    repeat_arrays = load(GE / "selected_repeat" / "selected_repeat_arrays.json")
    root_receipt = load(GE / "root_03" / "receipt.json")
    repeat_receipt = load(GE / "selected_repeat" / "receipt.json")
    q0_receipt = load(GE / "q0_smoke" / "receipt.json")
    root_gates = load(GE / "root_03" / "gates.json")
    repeat_gates = load(GE / "selected_repeat" / "gates.json")
    source_manifest = (VERIFY / "source_receipts" / "source.sha256").read_text().splitlines()
    local = ROOT / "credit_v1" / "solve_v1"
    comparisons = {
        "original_credit_v1_control": observer_values(local / "control" / "observers.json"),
        "matched_grid_baseline": observer_values(local / "grid_control" / "observers.json"),
        "partial_equilibrium_credit": observer_values(local / "credit" / "observers.json"),
        "general_equilibrium_selected": observer_values(GE / "root_03" / "observers.json"),
    }
    output = {
        "status": "pass" if not mismatches else "hash_mismatch",
        "reference_label": completed["reference_label"],
        "plan_sha256": completed["plan_sha256"],
        "reference_checkpoint_sha256": root_receipt["reference_checkpoint"]["sha256"],
        "reference_manifest_sha256": root_receipt["reference_manifest_sha256"],
        "selected": selected,
        "repeat": completed["repeat"],
        "completion": {"lifecycle_solves": completed["total_lifecycle_solves"], "elapsed_seconds": completed["elapsed_seconds"]},
        "repeat_verification": {
            "array_count": repeat_arrays["array_count"],
            "every_array_exact": all(value.get("exact") is True and value.get("max_abs") == 0.0 for value in repeat_arrays["arrays"].values()),
            "selected_repeat_target_fit_bytes_equal": (GE / "root_03" / "target_fit.csv").read_bytes() == (GE / "selected_repeat" / "target_fit.csv").read_bytes(),
            "selected_repeat_parameters_bytes_equal": (GE / "root_03" / "parameters.csv").read_bytes() == (GE / "selected_repeat" / "parameters.csv").read_bytes(),
            "target_fit_rows": table_count(GE / "root_03" / "target_fit.csv"),
            "parameter_rows": table_count(GE / "root_03" / "parameters.csv"),
        },
        "selected_gate_excerpt": {
            "renewal_residual": root_receipt["renewal_residual"],
            "relative_market_residual": root_receipt["ge_housing_certificate"]["relative_residual"],
            "queue_raw_l1": root_receipt["native_scaled_step"]["raw_distribution_l1"],
            "queue_adjusted_l1": root_receipt["native_scaled_step"]["distribution_l1_after_accounted_renewal_error"],
            "paygo_scaled_residual": root_receipt["paygo_residual"],
            "credit_floor_occupied_mass": root_gates["credit_solvency"]["occupied_mass_at_conservative_grid_floor"],
            "repeat_gates_bytes_equal": (GE / "root_03" / "gates.json").read_bytes() == (GE / "selected_repeat" / "gates.json").read_bytes(),
            "repeat_receipt_sha256": completed["repeat_receipt_sha256"],
        },
        "source_manifest_lines": source_manifest,
        "comparison_values": comparisons,
        "saved_supply_and_rooms_decomposition": {
            "q0_absolute_supply": q0_receipt["ge_housing_certificate"]["absolute_supply"],
            "qGE_absolute_supply": root_receipt["ge_housing_certificate"]["absolute_supply"],
            "Hs_qGE_over_Hs_q0_minus_1": root_receipt["ge_housing_certificate"]["absolute_supply"] / q0_receipt["ge_housing_certificate"]["absolute_supply"] - 1.0,
            "GE_rooms_over_original160_rooms_minus_1": root_receipt["cohort_summary"]["rooms_per_household"] / comparisons["original_credit_v1_control"]["mean_occupied_rooms_uncapped"] - 1.0,
        },
        "received_file_hashes": manifests,
        "hash_mismatches": mismatches,
        "notes": ["No conditional_cohort_state.pkl.gz checkpoint was retained locally.", "Contact-sheet Slurm submission was rejected by Torch for missing account allocation; original selected and repeat PNG packets were collected."],
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
