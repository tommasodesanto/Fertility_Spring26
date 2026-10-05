#!/usr/bin/env python3
"""Bridge the saved chain-11 Estate-A search point to a dated reference.

This performs exactly one fixed-price lifecycle solve. It never searches over
prices or parameters. The search point remains a provisional calibration.
Run in the pinned overnight source environment, with one Numba/BLAS thread.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np


LABEL = "0149_nm"
PRICE = 0.7551066122345627
H0 = 6.299624572680468
LOSS = 15.021541735257825
TARGET_PIN = "c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70"
WEIGHT_PIN = "f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4"
SOURCE_FILES = (
    "inputs.py", "reporting.py", "equilibrium.py", "native_price.py",
    "native_phase_b.py", "estate_contract.py", "storage.py",
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def near(a, b, tolerance=1e-10) -> bool:
    return np.isfinite(float(a)) and np.isfinite(float(b)) and abs(float(a)-float(b)) <= tolerance


def check_rows(actual: Path, expected: list[dict[str, str]], keys: tuple[str, ...], *, tol=1e-10) -> None:
    got = rows(actual)
    if len(got) != len(expected):
        raise RuntimeError(f"Report row count differs: {actual}")
    for index, (a, b) in enumerate(zip(got, expected)):
        if a.get(keys[0]) != b.get(keys[0]):
            raise RuntimeError(f"Report row identity differs: {actual}:{index}")
        for key in keys[1:]:
            x, y = a.get(key, ""), b.get(key, "")
            if x == y:
                continue
            if not x or not y or not near(x, y, tol):
                raise RuntimeError(f"Report differs: {actual}:{index}:{key}: {x} vs {y}")


def preflight(run_root: Path, inventory: Path, source_root: Path) -> dict:
    checkpoint = json.loads((run_root / "best_so_far.json").read_text())
    best = checkpoint["best"]
    if checkpoint["status"] != "provisional_until_fresh_native_postcheck" or best["label"] != LABEL:
        raise RuntimeError("Wrong provisional candidate")
    if not (near(best["loss"], LOSS, 1e-12) and near(best["price"], PRICE, 1e-14)
            and near(best["H0_derived"], H0, 1e-12)):
        raise RuntimeError("Chain-11 candidate numeric identity differs")
    if (best["target_fingerprint"], best["weight_fingerprint"]) != (TARGET_PIN, WEIGHT_PIN):
        raise RuntimeError("Lower-wealth target or weight contract differs")
    if best["closure"]["closure_mode"] != "population_one" or best["closure"]["population_scale"] != 1.:
        raise RuntimeError("Source reference is not normalized at one household")
    if best["experiment_flags"] != dict(birth_count_choice_enabled=True,
            birth_count_choice_cap=1, bequest_net_of_selling_cost=True,
            estate_flow_net_of_selling_cost=True):
        raise RuntimeError("Estate-A one-birth contract differs")
    case = run_root / LABEL
    phase = case / "phase_b_ge"
    selected = json.loads((phase / "selected.json").read_text())
    if selected["status"] != "passed" or not near(selected["selected_price"], PRICE, 1e-14):
        raise RuntimeError("Saved native price/repeat was not selected")
    report = phase / "selected_root"
    repeat = phase / "selected_repeat_final"
    archive = phase / "selected_repeat" / "stage" / "solution_arrays.npz"
    if not archive.is_file():
        raise RuntimeError("Saved exact-repeat native solution arrays missing")
    if len(list((report / "standard_diagnostics").glob("*.png"))) != 17:
        raise RuntimeError("Selected source diagnostic packet incomplete")
    if len(rows(report / "parameters.csv")) != 31 or len(rows(report / "target_fit.csv")) != 14:
        raise RuntimeError("Selected source tables incomplete")
    if len(rows(repeat / "parameters.csv")) != 31 or len(rows(repeat / "target_fit.csv")) != 14:
        raise RuntimeError("Saved native repeat tables incomplete")
    inv = json.loads(inventory.read_text())
    files = inv["files"]
    root = source_root.resolve()
    for name in SOURCE_FILES:
        rel = f"code/model/experiments/birth_count_choice/model/{name}"
        if rel not in files or digest(root / rel) != files[rel]:
            raise RuntimeError(f"Executing model source differs from overnight stage: {rel}")
    return dict(best=best, report=report, repeat=repeat, archive=archive,
                source_inventory_sha256=digest(inventory), bridge_sha256=digest(Path(__file__).resolve()),
                source_case=str(case), selected_sha256=digest(phase / "selected.json"),
                array_archive_sha256=digest(archive))


def check_arrays(solution, shared, archive: Path) -> dict:
    current = {key: value for key, value in vars(solution).items()
               if isinstance(value, np.ndarray) and value.dtype != object}
    shared_arrays = {"shared." + key: value for key, value in vars(shared).items()
                     if isinstance(value, np.ndarray) and value.dtype != object}
    with np.load(archive, allow_pickle=False) as saved:
        # H0 is deliberately replaced, so precomputed shared arrays may change.
        # Every solution array must still replay at the identical selected price.
        expected = set(current)
        source_solution = {key for key in saved.files if not key.startswith("shared.")
                           and key not in {"parameters.H0", "normalization.population"}}
        extra = source_solution - expected
        missing = expected - source_solution
        if extra or missing:
            raise RuntimeError(f"Saved array inventory differs: missing={sorted(missing)}, extra={sorted(extra)}")
        if len(current) != 78 or not set(shared_arrays).issubset(saved.files):
            raise RuntimeError("Complete 78-array selected solution or shared-array inventory missing")
        max_gap = 0.
        for key, value in current.items():
            prior = saved[key]
            if value.shape != prior.shape or value.dtype != prior.dtype:
                raise RuntimeError(f"Saved array shape/dtype differs: {key}")
            if value.dtype.kind in "biu":
                gap = 0. if np.array_equal(value, prior) else float("inf")
            else:
                if not (np.array_equal(np.isnan(value), np.isnan(prior))
                        and np.array_equal(np.isposinf(value), np.isposinf(prior))
                        and np.array_equal(np.isneginf(value), np.isneginf(prior))):
                    raise RuntimeError(f"Saved array nonfinite mask differs: {key}")
                finite = np.isfinite(value)
                gap = float(np.max(np.abs(value[finite]-prior[finite]))) if finite.any() else 0.
            if not np.isfinite(gap) or gap > 1e-10:
                raise RuntimeError(f"Saved array differs: {key}: {gap}")
            max_gap = max(max_gap, gap)
        if not (np.array_equal(saved["parameters.H0"], np.asarray([H0]))
                and np.array_equal(saved["normalization.population"], np.asarray([1.]))):
            raise RuntimeError("Saved normalized H0/population differs")
    return dict(solution_arrays=len(current), shared_arrays=len(shared_arrays), max_absolute_gap=max_gap)


def build(args: argparse.Namespace, proof: dict) -> Path:
    sys.path.insert(0, str(args.source_root.resolve() / "code/model/experiments/birth_count_choice"))
    from model.inputs import load_inputs
    from model.estate_contract import apply_experiment_flags, experiment_flags, rescore_report
    from model.reporting import build_context
    from model.native_price import solve_fixed_price
    from model.native_phase_b import observe_price
    from model.equilibrium import Budget
    from model.storage import StoredResult, save_case, load_case

    case = args.output_case.resolve()
    if case.exists():
        raise RuntimeError("Refusing to overwrite output case")
    case.mkdir(parents=True)
    try:
        native = case / "native"
        native.mkdir()
        best = proof["best"]
        P, grid = load_inputs(parameters=best["parameters"], external_inputs={"H0": [H0]})
        apply_experiment_flags(P, experiment_flags(1))
        deadline = time.time() + args.budget_seconds
        context = build_context(P, grid, native, price_start=PRICE,
                                deadline=deadline, max_lifecycle=1, closure="fixed_h0")
        budget = Budget(native, deadline, 1)
        live = solve_fixed_price(context, 0., PRICE, budget, "selected_bridge",
                                 native / "phase_b_ge" / "selected_bridge" / "stage")
        # This is the only Bellman/KFE solve. The observer reconstructs and
        # checks the stationary distribution using its saved policies.
        closure = observe_price(context, live, "selected_root", final=True)
        report = native / "phase_b_ge" / "selected_root"
        if not near(closure["population_scale"], 1., 1e-10) or not near(closure["H0_derived"] if "H0_derived" in closure else H0, H0, 1e-10):
            raise RuntimeError("Fixed-H0 bridge changed physical population or H0")
        for key in ("renewal_residual", "absolute_housing_residual", "actual_paygo_residual"):
            if abs(float(closure[key])) > 1e-6:
                raise RuntimeError(f"Stationary accounting fails: {key}")
        arrays = check_arrays(live["sol"], live["sd"], proof["archive"])
        check_rows(report / "target_fit.csv", rows(proof["report"] / "target_fit.csv"),
                   ("moment", "target", "model", "gap", "weight", "loss_contribution"))
        check_rows(report / "parameters.csv", rows(proof["report"] / "parameters.csv"),
                   ("parameter", "estimate"), tol=2e-12)
        revised, residual, _, _ = rescore_report(report, case)
        if len(revised) != 14 or not near(float(residual @ residual), LOSS, 1e-8):
            raise RuntimeError("New wealth-target score differs")
        check_rows(case / "target_fit_new_contract.csv", best["target_fit"],
                   ("moment", "target", "model", "gap", "weight", "loss_contribution"))
        shutil.copy2(report / "target_fit.csv", case / "target_fit.csv")
        shutil.copy2(report / "parameters.csv", case / "parameters.csv")
        shutil.copytree(report / "standard_diagnostics", case / "standard_diagnostics")
        contract = dict(closure="fixed_h0", experiment_flags=experiment_flags(1),
                        parameters=best["parameters"], external_inputs={"H0": [H0]},
                        target_fingerprint=TARGET_PIN, weight_fingerprint=WEIGHT_PIN,
                        experimental_economic_change=(
                            "One-intended-birth Estate A with post-interest soft financing; "
                            "ten coordinates estimated against lower aggregate wealth/earnings "
                            "target 4.45838713455674. Other targets and weights retained; "
                            "adult estate receiver remains none."),
                        source_status="provisional_search_checkpoint", native_bridge="one_fixed_price_solve",
                        population_scale=1., price=PRICE)
        (case / "input_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
        result = StoredResult(live["sol"], live["P"], grid, PRICE,
                              parameters=best["parameters"], label="provisional lower-wealth Estate-A reference")
        result.closure = dict(closure, closure_mode="fixed_h0", status="converged")
        result.report_directory = str(report)
        metadata = dict(closure=result.closure, report_directory=str(report),
                        input_contract_file="input_contract.json",
                        input_contract_sha256=digest(case / "input_contract.json"),
                        provenance=dict(source_case=proof["source_case"],
                                        source_selected_sha256=proof["selected_sha256"],
                                        source_array_archive_sha256=proof["array_archive_sha256"],
                                        source_inventory_sha256=proof["source_inventory_sha256"],
                                        bridge_sha256=proof["bridge_sha256"],
                                        provisional_calibration=True, optimizer_postcheck=False,
                                        lifecycle_solves=budget.used_lifecycle, array_comparison=arrays))
        save_case(result, case, metadata=metadata)
        load_case(case)
        (case / "bridge_receipt.json").write_text(json.dumps(metadata["provenance"], indent=2, sort_keys=True) + "\n")
        return case
    except BaseException:
        (case / "bridge_failed.json").write_text(json.dumps(dict(status="failed", provisional_source=True)) + "\n")
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True,
                        help=".../results/production_binary_chain_11/run")
    parser.add_argument("--inventory", type=Path, required=True,
                        help="pinned overnight A stage inventory.json")
    parser.add_argument("--source-root", type=Path, required=True,
                        help="mounted canonical repository root containing overnight pinned source")
    parser.add_argument("--output-case", type=Path, required=True,
                        help="new canonical case path, under a separate output root/cases")
    parser.add_argument("--budget-seconds", type=int, default=1800)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if args.budget_seconds < 300 or args.budget_seconds > 1800:
        parser.error("budget must be 300..1800 seconds")
    proof = preflight(args.run_root.resolve(), args.inventory.resolve(), args.source_root.resolve())
    if args.preflight_only:
        print(json.dumps(dict(status="preflight_passed_zero_solves", label=LABEL,
                              selected_sha256=proof["selected_sha256"],
                              array_archive_sha256=proof["array_archive_sha256"])))
        return
    path = build(args, proof)
    print(json.dumps(dict(status="provisional_reference_bridge_passed", case=str(path),
                          lifecycle_solves=1, optimizer_postcheck=False)))


if __name__ == "__main__":
    main()
