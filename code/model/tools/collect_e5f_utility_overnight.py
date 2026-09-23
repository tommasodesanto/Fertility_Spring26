#!/usr/bin/env python3
"""Read-only collector for the adaptive E5F utility comparison."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

CELLS = ("B_floor", "B_shares", "D_floor", "D_shares")
TARGET_COUNT = 13
GRAPH_COUNT = 17
SOURCE_COUNT = 641
BUNDLE_ROOT = None


class CollectionError(ValueError):
    pass


def read(path):
    return json.loads(Path(path).read_text())


def optional(path):
    return read(path) if Path(path).is_file() else {}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_path(path):
    """Resolve frozen remote bundle paths to the local mirror when present."""
    path = Path(path)
    if path.exists():
        return path
    if BUNDLE_ROOT is not None:
        marker = f"/{BUNDLE_ROOT.name}/"
        raw = str(path)
        if marker in raw:
            local = BUNDLE_ROOT / raw.split(marker, 1)[1]
            if local.exists():
                return local
    return path


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def expected_parameters(cell):
    return 17 if cell.endswith("floor") else 19


def validate_row_counts(targets, parameters, cell):
    if len(targets) != TARGET_COUNT or len(parameters) != expected_parameters(cell):
        raise CollectionError(f"full fit row counts differ: targets={len(targets)}, parameters={len(parameters)}")


def validate_science_rows(targets, objective_rows, loss):
    """Reconcile target values, weights, residuals and every SMM contribution."""
    if len(targets) != TARGET_COUNT or len(objective_rows) != TARGET_COUNT:
        raise CollectionError("objective and score must both contain 13 target rows")
    by_id = {row.get("restriction_id"): row for row in targets}
    frozen = {row.get("restriction_id"): row for row in objective_rows}
    if len(by_id) != TARGET_COUNT or set(by_id) != set(frozen):
        raise CollectionError("scored target identities differ from the frozen objective")
    scored, normalization = [], []
    for restriction_id, row in by_id.items():
        source = frozen[restriction_id]
        if float(row["target"]) != float(source["target"]):
            raise CollectionError(f"scored target changed: {restriction_id}")
        if row.get("role") != source.get("role"):
            raise CollectionError(f"scored target role changed: {restriction_id}")
        if restriction_id == "initial_normalization":
            normalization.append(row)
            if row.get("scored") is not False or row.get("actual_weight") is not None or row.get("loss_contribution") is not None:
                raise CollectionError("normalization was included in the scored objective")
            if float(row["target"]) != 2.1:
                raise CollectionError("separate initial normalization target is not 2.1")
            continue
        scored.append(row)
        if row.get("scored") is not True:
            raise CollectionError(f"empirical target is not scored: {restriction_id}")
        if float(row["actual_weight"]) != float(source["actual_weight"]):
            raise CollectionError(f"target weight changed: {restriction_id}")
        target, model = float(row["target"]), float(row["model"])
        gap, weight = float(row["gap"]), float(row["actual_weight"])
        contribution = float(row["loss_contribution"])
        if not all(math.isfinite(x) for x in (target, model, gap, weight, contribution)) or weight <= 0:
            raise CollectionError(f"nonfinite/invalid scored target fields: {restriction_id}")
        if not math.isclose(model-target, gap, rel_tol=1e-12, abs_tol=1e-12):
            raise CollectionError(f"target residual does not reconcile: {restriction_id}")
        if not math.isclose(weight*gap**2, contribution, rel_tol=1e-12, abs_tol=1e-9):
            raise CollectionError(f"weighted squared-gap contribution does not reconcile: {restriction_id}")
    if len(scored) != 12 or len(normalization) != 1:
        raise CollectionError("score must contain 12 scored targets and one separate normalization")
    if not math.isclose(math.fsum(float(row["loss_contribution"]) for row in scored),
                        float(loss), rel_tol=1e-12, abs_tol=1e-9):
        raise CollectionError("target loss contributions do not sum to objective")


def assert_one_target_system(target_hashes, expected):
    if set(target_hashes) != {expected}:
        raise CollectionError("mixed target/weight fingerprints across utility cells")


def selection_candidates(inventory, verified):
    """Return verified smoke/production cases; exact repeats are never selected."""
    return [(stage, cell, case_id, evidence)
            for (stage, cell, case_id), evidence in verified.items()
            if stage in ("smoke", "production")]


def count_statuses(inventory, cell):
    rows = [row for row in inventory if row["cell"] == cell]
    return {state: sum(row["status"] == state for row in rows)
            for state in ("unrun", "verified_scored", "rejected", "failed", "running", "incomplete", "collection_rejected")}


def copy_selected_support(evidence, destination):
    """Copy compact native receipts needed for independent interpretation."""
    destination = Path(destination)
    result_root, evaluation = evidence["evaluation"].parent, evidence["evaluation"]
    copied = []
    for name in ("initial_contract.json", "run_contract.json", "runtime_contract.json"):
        source = result_root / name
        if source.is_file():
            target = destination / name
            shutil.copy2(source, target)
            copied.append(str(target))
    for source_name, target_name in (("summary.json", "native_summary.json"),
                                     ("lifecycle_2023.csv", "lifecycle_2023.csv")):
        source = evaluation / "raw" / "repetition_01" / source_name
        if source.is_file():
            target = destination / target_name
            shutil.copy2(source, target)
            copied.append(str(target))
    for name in ("entry_wealth.json", "income_process.json", "wealth_grid.json"):
        source = evaluation / name
        if source.is_file():
            target = destination / name
            shutil.copy2(source, target)
            copied.append(str(target))
    return copied


def objective_and_contract(plan):
    """Resolve the pinned objective through the case plan's run contract."""
    contracts = plan.get("files", {})
    entry = contracts.get("run_contract")
    if not entry:
        raise CollectionError("plan lacks files.run_contract")
    path = resolve_path(entry["path"])
    if not path.is_file() or sha(path) != entry.get("sha256"):
        raise CollectionError("run contract missing or file pin mismatch")
    run = read(path)
    objective_info = run.get("working_objective")
    if not objective_info:
        raise CollectionError("run contract lacks working_objective")
    obj_path = resolve_path(objective_info["path"])
    if not obj_path.is_file() or sha(obj_path) != objective_info.get("sha256"):
        raise CollectionError("working objective missing or file pin mismatch")
    objective = read(obj_path)
    if canonical(objective) != plan.get("objective_canonical_sha256"):
        raise CollectionError("objective canonical fingerprint differs from plan")
    target = {key: value for key, value in objective.items() if key != "parameter_restrictions"}
    target_hash = canonical(target)
    if target_hash != plan.get("target_system_sha256"):
        raise CollectionError("target/weight identity differs from plan")
    return run, objective, target_hash


def validate_source_inventory(plan, manifest):
    inventory = plan.get("source_manifest", {}).get("files")
    if not isinstance(inventory, dict) or len(inventory) != SOURCE_COUNT:
        raise CollectionError("source inventory missing or not 641 files")
    fingerprint = canonical(inventory)
    expected = manifest.get("source_manifest_sha256")
    # The preparation contract fingerprints the full source manifest, while
    # the inherited battery source_inventory_sha256 fingerprints its file rows.
    if expected and fingerprint != expected:
        raise CollectionError("source inventory fingerprint differs from manifest")
    source_root = resolve_path(plan.get("source_root", ""))
    for name, expected_hash in inventory.items():
        source = source_root / name
        if not source.is_file() or sha(source) != expected_hash:
            raise CollectionError(f"source file missing or hash mismatch: {source}")
    return fingerprint


def validate_template_assets(plan, template, manifest):
    """Verify frozen source, runtime file, and generated-adapter pins once/cell."""
    validate_source_inventory(plan, manifest)
    files = plan.get("files", {})
    for name, item in files.items():
        path = resolve_path(item["path"])
        if not path.is_file() or sha(path) != item.get("sha256"):
            raise CollectionError(f"frozen runtime file pin mismatch: {name}")
    adapter = files.get("adapter", {})
    if (plan.get("adapter_path") != adapter.get("path")
            or plan.get("adapter_sha256") != adapter.get("sha256")
            or adapter.get("sha256") != manifest.get("adapter_source_change", {}).get("bundled_sha256")):
        raise CollectionError("generated utility adapter reference differs from the manifest and plan")
    preference_pin = files.get("preference_contract")
    if not preference_pin:
        raise CollectionError("cell template lacks generated preference contract")
    preference = read(resolve_path(preference_pin["path"]))
    runtime = preference.get("runtime_extension", {})
    if runtime != adapter:
        raise CollectionError("preference contract runtime adapter differs from plan adapter")
    for name, ref in preference.get("generated_files", {}).items():
        if name not in files or ref != files[name]:
            raise CollectionError(f"generated utility reference differs from plan: {name}")
    if preference.get("target_system_sha256") != manifest.get("target_system_sha256"):
        raise CollectionError("preference contract target-system fingerprint differs")
    if preference.get("parameter_rows") != expected_parameters(Path(template["path"]).parent.name):
        raise CollectionError("preference contract parameter row count differs")


def verify_case(case_dir, cell, manifest, template, pin=None):
    case_dir = Path(case_dir)
    evaluation = case_dir / "result" / "evaluation"
    status = read(case_dir / "status.json")
    plan_path = case_dir / "plan.json"
    plan_hash = sha(plan_path)
    status_plan_hash = status.get("plan_sha256", status.get("statusplan_sha256"))
    if not status_plan_hash or plan_hash != status_plan_hash:
        raise CollectionError("dynamic plan hash differs from status pin")
    plan = read(plan_path)
    if pin and pin.get("sha256") and plan_hash != pin["sha256"]:
        raise CollectionError("dynamic plan hash differs from manifest pin")
    template_path = resolve_path(template["path"])
    if template.get("sha256") and sha(template_path) != template["sha256"]:
        raise CollectionError("template file hash differs from manifest pin")
    if plan.get("objective_canonical_sha256") != template.get("objective_canonical_sha256"):
        raise CollectionError("case objective differs from cell template")
    if plan.get("target_system_sha256") != manifest.get("target_system_sha256"):
        raise CollectionError("case target-system fingerprint differs from manifest")
    if plan.get("structural_parameters") != status.get("parameters"):
        raise CollectionError("case parameters differ from status record")
    template_plan = read(template_path)
    if canonical(plan.get("parameter_bounds")) != canonical(template_plan.get("parameter_bounds")):
        raise CollectionError("case parameter bounds differ from frozen cell template")
    if canonical(plan.get("parameter_bounds")) != canonical(template.get("bounds")):
        raise CollectionError("manifest bounds differ from frozen cell template")
    if canonical(plan.get("files")) != canonical(template_plan.get("files")):
        raise CollectionError("case runtime/input file pins differ from frozen cell template")
    for field in ("source_root", "income_specification", "entry_specification",
                  "preference_specification", "wealth_grid"):
        if plan.get(field) != template_plan.get(field):
            raise CollectionError(f"case changed frozen economic/runtime field: {field}")
    for name, value in plan["structural_parameters"].items():
        low, high = map(float, plan["parameter_bounds"][name])
        if not math.isfinite(float(value)) or not low <= float(value) <= high:
            raise CollectionError(f"planned value outside actual bounds: {name}")

    if canonical(plan.get("source_manifest")) != canonical(template_plan.get("source_manifest")):
        raise CollectionError("dynamic case source inventory differs from cell template")
    run, objective, target_hash = objective_and_contract(plan)
    score_dir = evaluation / "scored_repetition_01"
    score = read(score_dir / "score.json")
    summary = read(evaluation / "summary.json")
    receipt = read(score_dir / "verified_evaluation_receipt.json")
    case_run_path = evaluation.parent / "run_contract.json"
    initial_path = evaluation.parent / "initial_contract.json"
    case_run, initial = read(case_run_path), read(initial_path)
    if sha(case_run_path) != summary.get("run_contract_sha256"):
        raise CollectionError("case run_contract.json hash differs from evaluation summary")
    if sha(initial_path) != summary.get("initial_solve_contract_sha256"):
        raise CollectionError("case initial_contract.json hash differs from evaluation summary")
    if initial.get("earnings_wealth_plan_sha256") != plan_hash:
        raise CollectionError("initial contract does not pin this dynamic case plan")
    if initial.get("structural_candidate") != plan.get("structural_parameters"):
        raise CollectionError("initial structural candidate differs from dynamic plan")
    source_hashes = initial.get("source_sha256")
    if not isinstance(source_hashes, dict) or len(source_hashes) != SOURCE_COUNT:
        raise CollectionError("initial contract source inventory is incomplete")
    if canonical(source_hashes) != canonical(template_plan["source_manifest"]["files"]):
        raise CollectionError("initial contract source inventory differs from frozen template")
    if canonical(source_hashes) != manifest.get("source_manifest_sha256"):
        raise CollectionError("initial contract source inventory differs from frozen manifest")
    if status.get("status") == "completed":
        if status.get("score_sha256") != sha(score_dir / "score.json"):
            raise CollectionError("completed status score hash differs from score.json")
        if status.get("summary_sha256") != sha(evaluation / "summary.json"):
            raise CollectionError("completed status summary hash differs from summary.json")
    for field in ("working_objective", "objective_source_files", "expected_graph_filenames",
                  "scorer", "validator", "wrapper_sha256", "source_root"):
        if case_run.get(field) != run.get(field):
            raise CollectionError(f"case run contract changed frozen field: {field}")
    runtime_path = evaluation.parent / "runtime_contract.json"
    runtime_contract = read(runtime_path)
    if runtime_contract.get("plan_sha256") != plan_hash:
        raise CollectionError("runtime contract does not pin the dynamic plan")
    if runtime_contract.get("frozen_source_manifest_sha256") != canonical(source_hashes):
        raise CollectionError("runtime contract source inventory fingerprint differs")
    expected_runtime_files = {name: files for name, files in plan["files"].items()
                              if name in {"adapter", "accounting", "income", "period_income", "earnings_adapter"}
                              or name.startswith("preference_")}
    if runtime_contract.get("additional_runtime_files") != expected_runtime_files:
        raise CollectionError("runtime adapter file references differ from dynamic plan")
    if summary.get("status") != "verified_scored_candidate" or summary.get("repetitions") != 1:
        raise CollectionError("evaluation summary is not a single scored candidate")
    if summary.get("objective_canonical_sha256") != template.get("objective_canonical_sha256"):
        raise CollectionError("summary objective fingerprint differs from template")
    if score.get("contract_sha256") != template.get("objective_canonical_sha256"):
        raise CollectionError("score objective fingerprint differs from template")
    if score.get("target_system_sha256", target_hash) != target_hash:
        raise CollectionError("score target-system fingerprint differs")

    # Verify all objective source fingerprints in both independent receipts.
    expected_sources = {}
    for name, item in run.get("objective_source_files", {}).items():
        source = resolve_path(item["path"])
        if not source.is_file():
            raise CollectionError(f"objective source missing: {source}")
        expected_sources[name] = canonical(read(source)) if item.get("hash_kind") == "canonical_json" else sha(source)
    if expected_sources != objective.get("source_fingerprints"):
        raise CollectionError("objective source fingerprint values differ from pinned source files")
    if score.get("source_fingerprints") != expected_sources or receipt.get("source_fingerprints") != expected_sources:
        raise CollectionError("score/receipt source fingerprints differ from pinned objective sources")
    if receipt.get("status") != "verified" or receipt.get("numerical_gates_verified") is not True:
        raise CollectionError("evaluation receipt did not verify numerical gates")
    receipt_hash = canonical(receipt)
    receipt_list = summary.get("evaluation_receipt_sha256", [])
    if receipt_hash not in receipt_list or score.get("evaluation_receipt_sha256") != receipt_hash:
        raise CollectionError("score/summary evaluation-receipt fingerprint mismatch")
    if receipt.get("checkpoint_sha256") != score.get("checkpoint_sha256"):
        raise CollectionError("score and receipt checkpoint hashes differ")

    targets, parameters = score.get("target_fit", []), score.get("parameters", [])
    expected_count = expected_parameters(cell)
    validate_row_counts(targets, parameters, cell)
    validate_science_rows(targets, objective.get("target_rows", []), score.get("loss"))
    if score.get("restriction_count") != TARGET_COUNT:
        raise CollectionError("score restriction_count differs from the 13-row target system")
    if score.get("free_parameter_count") != len(plan["parameter_bounds"]):
        raise CollectionError("score free-parameter dimension differs from the searched plan")
    if len({row.get("restriction_id") for row in targets}) != TARGET_COUNT:
        raise CollectionError("target rows contain duplicate/missing restriction IDs")
    if score.get("scored_moment_count") != 12:
        raise CollectionError("scored_moment_count differs from the 12 empirical restrictions")
    if len({row.get("parameter") for row in parameters}) != expected_count:
        raise CollectionError("parameter rows contain duplicate/missing parameter names")
    csv_targets = {row["restriction_id"]: row for row in csv.DictReader((score_dir / "target_fit.csv").open())}
    csv_parameters = {row["parameter"]: row for row in csv.DictReader((score_dir / "parameters.csv").open())}
    if set(csv_targets) != {row["restriction_id"] for row in targets}:
        raise CollectionError("target CSV row identity differs from score JSON")
    if set(csv_parameters) != {row["parameter"] for row in parameters}:
        raise CollectionError("parameter CSV row identity differs from score JSON")
    for row in targets:
        csvrow = csv_targets[row["restriction_id"]]
        for field in ("target", "model", "gap", "actual_weight", "loss_contribution"):
            if field in row and row[field] is not None and float(csvrow[field]) != float(row[field]):
                raise CollectionError(f"target CSV differs from JSON: {row['restriction_id']} {field}")
        if not math.isclose(float(row["model"]) - float(row["target"]), float(row["gap"]), rel_tol=1e-12, abs_tol=1e-12):
            raise CollectionError(f"target gap does not reconcile: {row['restriction_id']}")
    for row in parameters:
        csvrow = csv_parameters[row["parameter"]]
        if float(csvrow["estimate"]) != float(row["estimate"]):
            raise CollectionError(f"parameter CSV differs from JSON: {row['parameter']}")
        name = row["parameter"]
        if name in plan["structural_parameters"] and not math.isclose(float(row["estimate"]), float(plan["structural_parameters"][name]), rel_tol=1e-12, abs_tol=1e-12):
            raise CollectionError(f"parameter estimate differs from proposed value: {name}")
    loss = float(score["loss"])
    if not math.isfinite(loss) or not math.isclose(loss, float(summary["loss"]), rel_tol=0, abs_tol=1e-12):
        raise CollectionError("score loss differs from summary")

    graph_rows = summary.get("original_graphs", [])
    expected_names = set(run.get("expected_graph_filenames", []))
    if len(graph_rows) != GRAPH_COUNT or len(expected_names) != GRAPH_COUNT:
        raise CollectionError("standard diagnostic graph count changed")
    if {Path(item["path"]).name for item in graph_rows} != expected_names:
        raise CollectionError("standard diagnostic graph identity changed")
    for item in graph_rows:
        graph = evaluation / "raw" / item["path"]
        if not graph.is_file() or sha(graph) != item.get("sha256"):
            raise CollectionError(f"standard diagnostic graph missing/hash mismatch: {item['path']}")

    checkpoint = evaluation / "raw" / "repetition_01" / "initial_state.pkl.gz"
    if receipt.get("checkpoint_sha256") != score.get("checkpoint_sha256"):
        raise CollectionError("score and receipt checkpoint hashes differ")
    # Nonselected checkpoint bytes are not hashed to avoid repeated large-file
    # I/O. The selected checkpoint is hashed once after selection.
    return {"plan": plan, "score": score, "summary": summary, "receipt": receipt,
            "target_hash": target_hash, "case_dir": case_dir, "evaluation": evaluation,
            "loss": loss, "checkpoint_sha256": score["checkpoint_sha256"],
            "checkpoint_path": checkpoint}


def _count_solves(evaluation):
    records = []
    if evaluation.is_dir():
        for path in evaluation.glob("raw/repetition_*/stationary_solves.json"):
            try:
                records.extend(read(path))
            except (OSError, ValueError, TypeError):
                pass
    heartbeat = optional(evaluation / "raw/heartbeat.json") if evaluation.is_dir() else {}
    started = max(len(records), int(heartbeat.get("stationary_solves", 0)))
    completed = sum(item.get("status") == "completed" for item in records)
    return {"started_solves": started, "completed_solves": completed,
            "unfinished_solves": max(0, started - completed)}


def _attempt_status(case_dir, evaluation, verified, error=None):
    status = optional(case_dir / "status.json")
    native_failure = optional(evaluation / "raw/failure.json")
    if verified:
        return "verified_scored", ""
    if error:
        return "collection_rejected", error
    if status.get("status") == "incomplete_timeout":
        return "incomplete", status.get("error", "timeout")
    if status.get("status") == "numerical_gate_rejection":
        return "rejected", native_failure.get("error", status.get("error", "numerical gate rejection"))
    if status.get("status") == "fatal_contract_or_code_error":
        return "failed", status.get("error", "fatal contract/code error")
    if native_failure:
        message = str(native_failure.get("error", "native failure"))
        kind = str(native_failure.get("error_type", ""))
        if "timeout" in (kind + " " + message).lower() or "time limit" in message.lower():
            return "incomplete", message
        if kind == "InfeasibleThetaError":
            return "rejected", message
        return "failed", message
    if status.get("status") == "completed":
        return "collection_rejected", "runner marked case completed but score/summary is missing"
    if status.get("status") == "started" or optional(case_dir / "heartbeat.json"):
        return "running", "case has a live runner status/heartbeat and no score yet"
    return "incomplete", "case started without a verified score/status"


def collect(manifest_path, results_root, output):
    global BUNDLE_ROOT
    manifest_path, results_root, output = map(Path, (manifest_path, results_root, output))
    BUNDLE_ROOT = manifest_path.parent.resolve()
    manifest = read(manifest_path)
    if manifest.get("schema") != "e5f_utility_overnight_manifest_v1" or manifest.get("status") != "ready":
        raise CollectionError("utility manifest schema/status mismatch")
    if set(manifest.get("templates", {})) != set(CELLS):
        raise CollectionError("manifest must contain all four utility cells")
    if int(manifest.get("source_file_count", -1)) != SOURCE_COUNT:
        raise CollectionError("manifest source inventory count changed")
    if (manifest.get("workers") != {"B_floor": 18, "B_shares": 18, "D_floor": 2, "D_shares": 2}
            or manifest.get("total_workers") != 40
            or manifest.get("maximum_proposals_per_worker") != 18
            or manifest.get("maximum_objectives") != {"production": 720, "smoke": 8, "total": 736, "verification": 8}):
        raise CollectionError("worker allocation or maximum planned-case cap changed")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite collection directory: {output}")
    output.mkdir(parents=True)

    templates = {cell: read(resolve_path(manifest["templates"][cell]["path"])) for cell in CELLS}
    target_hashes = {plan.get("target_system_sha256") for plan in templates.values()}
    assert_one_target_system(target_hashes, manifest.get("target_system_sha256"))
    for cell in CELLS:
        validate_template_assets(templates[cell], manifest["templates"][cell], manifest)
        _, _, target_hash = objective_and_contract(templates[cell])
        if target_hash != manifest["target_system_sha256"]:
            raise CollectionError(f"{cell}: template target/weight identity differs from manifest")

    inventory, verified = [], {}
    case_specs = []
    # The runner's hard cap is 18 proposals for each of 40 worker slots.
    task_id = 0
    for cell in CELLS:
        workers = manifest["workers"][cell]
        for local_worker in range(1, int(workers) + 1):
            task_id += 1
            # Paths follow the runner's task allocation, not a shared parameter
            # anchor: each adaptive worker can move its own center.
            for proposal in range(1, 19):
                case_id = f"worker{task_id:02d}_proposal{proposal:02d}"
                case_specs.append(("production", cell, case_id, None,
                                   {"task_id": task_id, "proposal_index": proposal}))
    # Eight smokes: two cases in each cell. The runner's smoke stage uses these
    # exact IDs and makes smoke_01 the common parameter anchor within each cell.
    for cell in CELLS:
        for number in (1, 2):
            case_specs.append(("smoke", cell, f"smoke_{number:02d}", None,
                               {"smoke_index": number, "common_parameter_anchor": number == 1}))
    # Repetitions are verification evidence, never members of the selection pool.
    for cell in CELLS:
        for number in (1, 2):
            case_specs.append(("verification", cell, f"selected_repeat_{number:02d}", None,
                               {"repeat_index": number}))

    target_hash_by_cell = {}
    for stage, cell, case_id, pin, metadata in case_specs:
        case_dir = results_root / stage / cell / case_id
        evaluation = case_dir / "result" / "evaluation"
        row = {"stage": stage, "cell": cell, "case_id": case_id,
               "case_directory": str(case_dir), "status": "unrun", "loss": None, **metadata}
        if not case_dir.exists():
            inventory.append(row)
            continue
        row.update(_count_solves(evaluation))
        score_path = evaluation / "scored_repetition_01" / "score.json"
        summary_path = evaluation / "summary.json"
        if score_path.is_file() and summary_path.is_file():
            try:
                evidence = verify_case(case_dir, cell, manifest, manifest["templates"][cell], pin)
                target_hash_by_cell[cell] = evidence["target_hash"]
                row.update(status="verified_scored", loss=evidence["loss"],
                           score_path=str(score_path), checkpoint_sha256=evidence["checkpoint_sha256"])
                verified[(stage, cell, case_id)] = evidence
            except (AssertionError, CollectionError, KeyError, OSError, TypeError, ValueError) as exc:
                row["status"], row["validity_exception"] = _attempt_status(
                    case_dir, evaluation, False, f"{type(exc).__name__}: {exc}")
        else:
            row["status"], row["validity_exception"] = _attempt_status(case_dir, evaluation, False)
        inventory.append(row)

    # Keep all collected score rows, including smoke and adaptive proposals.
    all_targets, all_parameters = [], []
    for (stage, cell, case_id), evidence in verified.items():
        for item in evidence["score"]["target_fit"]:
            all_targets.append({"stage": stage, "cell": cell, "case_id": case_id, **item})
        for item in evidence["score"]["parameters"]:
            item = dict(item)
            bounds = evidence["plan"].get("parameter_bounds", {})
            name = item["parameter"]
            if name in bounds:
                low, high = map(float, bounds[name])
                estimate = float(item["estimate"])
                item.update(actual_lower=low, actual_upper=high,
                            near_actual_bound=min(estimate-low, high-estimate) <= .01*(high-low),
                            restriction_type="searched")
            else:
                item.update(actual_lower=None, actual_upper=None, near_actual_bound=None,
                            restriction_type="externally_fixed_or_derived")
            all_parameters.append({"stage": stage, "cell": cell, "case_id": case_id, **item})
    write_csv(output / "all_target_fits.csv", all_targets)
    write_csv(output / "all_parameters.csv", all_parameters)

    selected, receipts, copied_graphs, selected_support = [], [], [], []
    for cell in CELLS:
        candidate_rows = [(stage, c, case_id, evidence)
                          for stage, c, case_id, evidence in selection_candidates(inventory, verified)
                          if c == cell]
        candidate_rows.sort(key=lambda row: row[3]["loss"])
        candidates = []
        for stage, _, case_id, evidence in candidate_rows:
            checkpoint = evidence["checkpoint_path"]
            if checkpoint.is_file() and sha(checkpoint) == evidence["checkpoint_sha256"]:
                candidates.append((stage, case_id, evidence))
                break
            for row in inventory:
                if row["stage"] == stage and row["cell"] == cell and row["case_id"] == case_id:
                    row["status"] = "collection_rejected"
                    row["validity_exception"] = "selected checkpoint is missing or its bytes differ from the score/receipt hash"
                    break
            verified.pop((stage, cell, case_id), None)
        if not candidates:
            continue
        stage, case_id, best = candidates[0]
        selected_record = {"cell": cell, "stage": stage, "case_id": case_id,
                           "loss": best["loss"], "checkpoint_sha256": best["checkpoint_sha256"],
                           "exact_repeats_are_selection": False}
        selected.append(selected_record)
        dest = output / cell / "selected"
        dest.mkdir(parents=True, exist_ok=True)
        score_dir = best["evaluation"] / "scored_repetition_01"
        for name in ("score.json", "target_fit.csv", "parameters.csv", "verified_evaluation_receipt.json"):
            shutil.copy2(score_dir / name, dest / name)
        shutil.copy2(best["evaluation"] / "summary.json", dest / "wrapper_summary.json")
        shutil.copy2(best["case_dir"] / "plan.json", dest / "plan.json")
        selected_support.extend({"cell": cell, "path": str(Path(path).relative_to(output))}
                               for path in copy_selected_support(best, dest))
        # Actual plan bounds are recorded beside every estimated parameter.
        annotated = []
        for row in best["score"]["parameters"]:
            item = dict(row)
            name = item["parameter"]
            if name in best["plan"]["parameter_bounds"]:
                low, high = map(float, best["plan"]["parameter_bounds"][name])
                item.update(actual_lower=low, actual_upper=high,
                            near_actual_bound=min(float(item["estimate"]) - low,
                                                  high - float(item["estimate"])) <= .01 * (high-low),
                            restriction_type="searched")
            else:
                item.update(actual_lower=None, actual_upper=None, near_actual_bound=None,
                            restriction_type="externally_fixed_or_derived")
            annotated.append(item)
        write_csv(dest / "parameters_actual_bounds.csv", annotated)
        graph_dir = dest / "standard_diagnostics"
        graph_dir.mkdir()
        for item in best["summary"]["original_graphs"]:
            source = best["evaluation"] / "raw" / item["path"]
            target = graph_dir / Path(item["path"]).name
            shutil.copy2(source, target)
            copied_graphs.append({"cell": cell, "path": str(target.relative_to(output)),
                                  "sha256": sha(target), "original_sha256": item["sha256"]})

        smoke_evidence = verified.get(("smoke", cell, "smoke_01"))
        if smoke_evidence:
            smoke_dest = output / cell / "smoke_anchor"
            smoke_dest.mkdir(parents=True, exist_ok=True)
            smoke_score_dir = smoke_evidence["evaluation"] / "scored_repetition_01"
            for name in ("score.json", "target_fit.csv", "parameters.csv", "verified_evaluation_receipt.json"):
                shutil.copy2(smoke_score_dir / name, smoke_dest / name)
            shutil.copy2(smoke_evidence["evaluation"] / "summary.json", smoke_dest / "wrapper_summary.json")
            shutil.copy2(smoke_evidence["case_dir"] / "plan.json", smoke_dest / "plan.json")
            selected_support.extend({"cell": cell, "stage": "smoke_anchor",
                                     "path": str(Path(path).relative_to(output))}
                                    for path in copy_selected_support(smoke_evidence, smoke_dest))

        repeats = []
        for number in (1, 2):
            item = verified.get(("verification", cell, f"selected_repeat_{number:02d}"))
            status_path = results_root / "verification" / cell / f"selected_repeat_{number:02d}" / "status.json"
            if status_path.is_file():
                repeat_row = read(status_path)
                if item:
                    receipt_path = item["evaluation"] / "scored_repetition_01" / "verified_evaluation_receipt.json"
                    receipt_dest = output / cell / "verification" / f"repeat_{number:02d}_receipt.json"
                    receipt_dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(receipt_path, receipt_dest)
                    repeat_row["receipt_copy"] = str(receipt_dest.relative_to(output))
                repeats.append(repeat_row)
        receipts.append({**selected_record, "repeats": repeats,
                         "repeat_count": len([item for item in repeats if item.get("status") == "completed"]),
                         "selection_includes_smoke": True})

    # Detect manifest target inconsistencies even if no case happened to run.
    if set(target_hash_by_cell.values()) - {manifest["target_system_sha256"]}:
        raise CollectionError("collected case has a mixed target/weight fingerprint")
    by_cell = {cell: count_statuses(inventory, cell) for cell in CELLS}
    solve_counts = {name: sum(row.get(name, 0) for row in inventory)
                    for name in ("started_solves", "completed_solves", "unfinished_solves")}
    attempted = sum(row["status"] != "unrun" for row in inventory)
    expected_total = 40 * 18 + 8 + 8
    if len(inventory) != expected_total:
        raise CollectionError(f"planned-case count mismatch: {len(inventory)} != {expected_total}")
    write(output / "inventory.json", inventory)
    write_csv(output / "inventory.csv", inventory)
    write(output / "selection.json", selected)
    write(output / "verification_receipts.json", receipts)
    write(output / "figure_manifest.json", copied_graphs)
    write(output / "selected_support_manifest.json", selected_support)
    write(output / "collection_summary.json", {
        "manifest_sha256": sha(manifest_path), "source_manifest_sha256": manifest["source_manifest_sha256"],
        "target_system_sha256": manifest["target_system_sha256"], "source_file_count": SOURCE_COUNT,
        "planned_cases": expected_total, "attempted_cases": attempted,
        "case_counts_by_cell": by_cell, "stationary_solves": solve_counts,
        "selected_cells": len(selected), "selection_stages": ["smoke", "production"],
        "verification_is_selection_pool": False,
    })
    (output / "README.md").write_text(
        "# Utility overnight collection\n\n"
        "Read-only collection of the adaptive 40-worker comparison. Selection is the lowest verified loss among smoke and production cases; exact repeats are verification only. The inventory includes every planned case, with unrun, incomplete, rejected and collection-rejected cases shown separately. Each selected cell has its full 13-target fit, all 17 or 19 estimated parameters with actual bounds, and the unchanged 17 standard diagnostic graphs. Checkpoint hashes are recorded; checkpoint files are not copied or loaded.\n")
    return {"attempted": attempted, "planned": expected_total, "verified": len(verified),
            "selected_cells": len(selected), "output": str(output)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(collect(args.manifest, args.results_root, args.output), sort_keys=True))


if __name__ == "__main__":
    main()
