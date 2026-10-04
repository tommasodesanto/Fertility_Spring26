"""One fixed-point Estate-A cap-two GE at the verified binary winner."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection/binary"
OUTPUT = ROOT / "output/model/experiments/birth_count_choice/cap2_at_binary_winner_v1"
SOURCE_SHA256 = "2e104f260050b802370de3d8bdd0d0792c8a7af9e13be849fadfdd6171106f71"
SOURCE_LOSS = 21.275413361071312
BUDGET_SECONDS = 1500


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def inputs_and_receipt():
    from model import calibration, estate_contract
    from model.inputs import load_inputs

    source_path = SOURCE / "provenance/search_completed.json"
    if sha(source_path) != SOURCE_SHA256:
        raise RuntimeError("Binary-winner search receipt SHA-256 drift")
    source = json.loads(source_path.read_text())
    completed = json.loads((SOURCE / "provenance/completed.json").read_text())
    native = json.loads((SOURCE / "provenance/native_postcheck/completed.json").read_text())
    selected = source["selected"]
    if not (source["arm"] == "binary" and source["birth_cap"] == 1 and source["chain"] == 1
            and completed["status"] == "selected_numerically_verified"
            and native["status"] == "full_native_postcheck_passed"
            and native["search_receipt_sha256"] == SOURCE_SHA256
            and completed["repeat"]["status"] == "exact_full_ge_repeat_passed"
            and selected["status"] == "passed" and selected["loss"] == SOURCE_LOSS
            and selected["parameters"] == completed["selected"]["parameters"]):
        raise RuntimeError("Verified binary-winner identity drift")
    point = selected["parameters"]
    contract, target_pin, weight_pin, bounds = estate_contract.contract()
    if (source["target_fingerprint"], source["weight_fingerprint"]) != (target_pin, weight_pin):
        raise RuntimeError("Estate-A target/weight fingerprint drift")
    if set(point) != set(bounds) or any(not lo <= point[k] <= hi for k, (lo, hi) in bounds.items()):
        raise RuntimeError("Winner's ten coordinates or bounds drift")
    source_rows = rows(SOURCE / "selected_root/target_fit_new_contract.csv")
    source_params = rows(SOURCE / "selected_root/parameters_estate_a.csv")
    if len(source_rows) != 14 or len(source_params) != 31 or source_rows != selected["target_fit"]:
        raise RuntimeError("Downloaded binary-winner full tables drift")
    if source_params != native["parameters"] or any(
            float(next(r for r in source_params if r["parameter"] == k)["estimate"]) != v
            for k, v in point.items()):
        raise RuntimeError("Downloaded winner parameter table drift")
    if [{k: r[k] for k in ("moment", "target", "weight", "role")} for r in source_rows] != contract:
        raise RuntimeError("Binary-winner target rows differ from live contract")
    P, grid = load_inputs(parameters=point)
    estate_contract.apply_experiment_flags(P, estate_contract.experiment_flags(1))
    if calibration.effective_input_fingerprint(P, grid) != selected["input_fingerprint"]:
        raise RuntimeError("Rebuilt winner caller inputs/grid differ from verified source")
    derived_h0 = float(selected["H0_derived"])
    if derived_h0 != float(selected["closure"]["H0_derived"]):
        raise RuntimeError("Binary-winner derived H0 drift")
    Q, other_grid = load_inputs(parameters=point, external_inputs={"H0": [derived_h0]})
    estate_contract.apply_experiment_flags(Q, estate_contract.experiment_flags(2))
    if not np.array_equal(grid, other_grid):
        raise RuntimeError("Wealth grid drift")
    changed = sorted(k for k in vars(P) if not np.array_equal(np.asarray(getattr(P, k)), np.asarray(getattr(Q, k))))
    if changed != ["H0", "birth_count_choice_cap"]:
        raise RuntimeError(f"Unexpected economic input difference: {changed}")
    return source, point, derived_h0, Q, other_grid, target_pin, weight_pin


def preflight():
    from model.reporting import build_context, production_model_facade
    from model.engine import solver
    source, point, h0, P, grid, target_pin, weight_pin = inputs_and_receipt()
    with tempfile.TemporaryDirectory(prefix="cap2_winner_preflight_") as directory:
        context = build_context(P, grid, directory, price_start=source["selected"]["price"],
                                deadline=time.time() + 120, max_lifecycle=32, closure="fixed_h0")
        if len(context["manifest"]["full_parameter_table"]) != 31:
            raise RuntimeError("Native parameter table drift")
        if production_model_facade().solve_markov_income_at_prices is not solver.solve_markov_income_at_prices:
            raise RuntimeError("Wrong birth-count engine bound")
        adapters = [r["function"] for r in context["birth_count_observer_adapters"]]
    return dict(status="passed_zero_solves", source_search_sha256=SOURCE_SHA256,
                source_loss=SOURCE_LOSS, source_chain=1, source_point=point,
                source_input_fingerprint=source["selected"]["input_fingerprint"],
                changed_executed_fields=["H0", "birth_count_choice_cap"],
                H0_source="binary winner's derived population-one coefficient",
                fixed_H0=h0, birth_cap=2, closure="fixed_h0", price_start=source["selected"]["price"],
                target_fingerprint=target_pin, weight_fingerprint=weight_pin,
                target_rows=14, parameter_rows=31, max_lifecycle=32,
                native_budget_seconds=BUDGET_SECONDS, external_limit_seconds=1800,
                observer_functions_adapted=adapters, output_root=str(OUTPUT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    receipt = preflight()
    if args.preflight:
        print(json.dumps(receipt, indent=2, sort_keys=True))
        return
    from model import estate_contract
    from model.workflow import run_stationary
    source, point, h0, _, _, target_pin, weight_pin = inputs_and_receipt()
    result, case = run_stationary(point, {"H0": [h0]}, {},
        price_guess=source["selected"]["price"], budget_seconds=BUDGET_SECONDS,
        max_lifecycle=32, closure="fixed_h0", output_root=OUTPUT,
        experiment_flags=estate_contract.experiment_flags(2),
        parameter_file_metadata=dict(source_search_sha256=SOURCE_SHA256,
            source_selected_loss=SOURCE_LOSS, source_target_fingerprint=target_pin,
            source_weight_fingerprint=weight_pin, source_closure="population_one",
            H0_source="derived binary-winner coefficient; held fixed in diagnostic"))
    fit = rows(case / "target_fit_new_contract.csv")
    params = rows(case / "parameters.csv")
    figures = list((case / "standard_diagnostics").glob("*.png"))
    if len(fit) != 14 or len(params) != 31 or len(figures) != 17:
        raise RuntimeError("Incomplete final fit, parameter table, or standard figures")
    print(json.dumps(dict(status="complete", case=str(case), price=result.price,
                          closure=result.closure, target_rows=len(fit), parameter_rows=len(params),
                          standard_figures=len(figures), new_contract_loss=json.loads(
                              (case / "estate_a_rescore_receipt.json").read_text())["loss"]), indent=2))


if __name__ == "__main__":
    main()
