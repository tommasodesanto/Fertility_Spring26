#!/usr/bin/env python3
"""Collect bounded direct-search cases and write the next immutable E5F plan.

No solver or optimizer is run here. Joint patterns combine successful coordinate
moves or rebalance the first-child and per-child room floors. All original
moments and weights determine candidate ranking.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
import shutil
from pathlib import Path

import run_e5f_bounded_calibration_refinement as adapter


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def collect(plan_path, expected, require_complete=True):
    plan = adapter.load_plan(plan_path, expected)
    rows, fits, parameters, provenance = [], [], [], set()
    missing = []
    for case in plan["cases"]:
        out = plan_path.parent / case["output"]
        receipt_path = out / "case_receipt.json"
        if not receipt_path.exists():
            missing.append(case["id"])
            continue
        receipt = adapter.read_json(receipt_path)
        if receipt["status"] != "complete" or receipt["plan_sha256"] != expected or receipt["case_id"] != case["id"]:
            raise RuntimeError("Case receipt belongs to another plan or is incomplete")
        for name, sha in receipt["artifact_sha256"].items():
            adapter.verify(out / name, sha)
        s, models, gates = adapter.validate_result(out, plan, case)
        import collect_e5f_transition_calibration as collector
        shared = {name: s[name] for name in ("target_measurements", "dated_measurement_contract", "renewal_accounting_contract", "external_closure_contract", "model_profile")}
        shared["population_bridge"] = collector.population_bridge_contract(s["population_bridge"])
        provenance.add(json.dumps(shared, sort_keys=True))
        rows.append({"id": case["id"], "label": case["label"], "loss": float(s["best_candidate"]["transition_loss"]),
            "first_birth_rooms": models["housing_increment_0to1"], "mean_rooms": models["aggregate_mean_occupied_rooms_18_85"],
            "ownership": models["own_rate"], "elapsed_seconds": receipt["elapsed_seconds"],
            "occupied_value_drops": receipt["policy_array_diagnostic"]["occupied_negative_steps"],
            "budget_excess_share": receipt["budget_diagnostic"]["budget_excess_share"],
            "plan": str(plan_path.resolve()), "plan_sha256": expected,
            "summary": str((out / "summary.json").resolve()), "summary_sha256": adapter.digest(out / "summary.json")})
        for row in adapter.read_csv(out / "target_fit_long.csv"):
            fits.append({"case_id": case["id"], "label": case["label"], **row})
        for row in adapter.read_csv(out / "parameter_table.csv"):
            parameters.append({"case_id": case["id"], "label": case["label"], **row})
    if len(provenance) > 1:
        raise RuntimeError("Mixed scientific/population/measurement contracts")
    report = plan_path.parent / "report"
    report.mkdir(exist_ok=True)
    rows.sort(key=lambda r: (r["loss"], r["id"]))
    write_csv(report / "all_candidates.csv", rows)
    write_csv(report / "all_target_fits.csv", fits)
    write_csv(report / "all_parameters.csv", parameters)
    status = {"status": "complete" if not missing else "partial", "plan_sha256": expected,
              "completed_cases": len(rows), "missing_cases": missing, "best": rows[0] if rows else None,
              "production_promoted": False}
    adapter.write_json(report / "summary.json", status)
    if rows:
        adapter.write_json(plan_path.parent.parent / "latest_completed_case.json", {"stage": plan["stage"], "latest": max(rows, key=lambda r: (Path(r["summary"]).stat().st_mtime, r["id"]))})
        best_path = plan_path.parent.parent / "best_so_far.json"
        prior = adapter.read_json(best_path) if best_path.exists() else {"loss": math.inf}
        prior_loss = float(prior.get("loss", prior.get("best", {}).get("loss", math.inf)))
        if rows[0]["loss"] < prior_loss:
            adapter.write_json(best_path, {"stage": plan["stage"], "loss": rows[0]["loss"], "best": rows[0], "production_promoted": False})
    if require_complete and missing:
        raise RuntimeError(f"Incomplete stage: {missing}")
    return plan, status, rows


def fresh_plan(prior, out, stage):
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"Refusing nonempty plan directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    plan = copy.deepcopy(prior)
    plan.update(stage=stage, cases=[], input_sha256={}, production_promoted=False)
    plan["planner_sha256"] = adapter.digest(__file__)
    return plan


def add_case(plan, out, center, label, *, panel_id=1, panel_size=1, panel_design="mixed", radius=.005, panel_seed=2026090506):
    i = len(plan["cases"]) + 1
    path = out / f"center_{i:03d}.json"
    adapter.write_json(path, center)
    plan["cases"].append({"id": i, "label": label, "center": path.name,
        "center_sha256": adapter.digest(path), "panel_task_id": panel_id,
        "panel_size": panel_size, "panel_design": panel_design,
        "panel_seed": panel_seed, "radius": radius, "output": f"task_{i:03d}"})


def finish_plan(plan, out):
    path = out / "plan.json"
    adapter.write_json(path, plan)
    print(json.dumps({"path": str(path), "sha256": adapter.digest(path), "cases": len(plan["cases"])}))


def coordinate(prior, status, rows, out):
    # Both smoke paths must reproduce the same original candidate exactly.
    if len(rows) != 2 or rows[0]["loss"] != rows[1]["loss"]:
        raise RuntimeError("Two exact smoke repeats required")
    for row in rows:
        receipt = adapter.read_json(Path(row["summary"]).parent / "case_receipt.json")
        if not receipt.get("reference", {}).get("exact_twelve_row_fit"):
            raise RuntimeError("Missing exact inherited-candidate reference check")
    center = adapter.read_json(rows[0]["summary"])
    plan = fresh_plan(prior, out, "round1_small_coordinate")
    plan["starting_summary_sha256"] = rows[0]["summary_sha256"]
    for i in range(1, 24):
        add_case(plan, out, center, "anchor" if i == 1 else f"coordinate_{(i-2)//2}_{'minus' if i%2==0 else 'plus'}",
                 panel_id=i, panel_size=23, panel_design="coordinate")
    finish_plan(plan, out)


def joint(prior, status, rows, out):
    if prior["stage"] != "round1_small_coordinate" or len(rows) != 23:
        raise RuntimeError("Joint proposals require the complete small coordinate panel")
    import run_e5f_transition_calibration as calibration
    by_id = {r["id"]: r for r in rows}
    anchor = adapter.read_json(by_id[1]["summary"])
    domain = anchor["panel_design"]["domain"]
    units = list(anchor["panel_design"]["unit_vector"])
    center = copy.deepcopy(adapter.read_json(rows[0]["summary"]))
    best_units = list(center["panel_design"]["unit_vector"])
    gains = []
    direction = [0.0] * 11
    for j in range(11):
        minus, plus = by_id[2+2*j], by_id[3+2*j]
        winner = min((minus, plus), key=lambda r: (r["loss"], r["id"]))
        gain = by_id[1]["loss"] - winner["loss"]
        if gain > 1e-8:
            direction[j] = -.005 if winner["id"] % 2 == 0 else .005
            gains.append((gain, j))
    top = {j for _, j in sorted(gains, reverse=True)[:3]}
    top_direction = [d if j in top else 0.0 for j,d in enumerate(direction)]
    plan = fresh_plan(prior, out, "round2_joint_patterns")
    plan.update(starting_summary_sha256=rows[0]["summary_sha256"],
        search_rule="full-objective direct patterns; no inverse Jacobian; all existing parameters remain free",
        coordinate_gain_direction=direction, top_three_indices=sorted(top))
    names = [d["name"] for d in domain]
    index_floor, index_jump = names.index("hbar_child_rooms"), names.index("hbar_first_child_jump")
    seen = {tuple(best_units)}
    patterns = []
    for label, d in (("all_improving_coordinates", direction), ("top_three_improving_coordinates", top_direction)):
        for scale in (.5, 1., 2.):
            patterns.append((f"{label}_scale_{scale:g}", [u+scale*v for u,v in zip(best_units,d)], None))
    # For m>0, hbar(m)=jump+m*floor. Increasing jump by d and reducing
    # floor by d/3 preserves m=3 needs and raises m=1 needs by 2d/3.
    for compensated in (False, True):
        for delta in (.025, .05, .10):
            u = [x+(d if compensated else 0) for x,d in zip(best_units,direction)]
            theta = center["best_candidate"]["theta"]
            proposed = {index_jump: theta["hbar_first_child_jump"]+delta,
                        index_floor: theta["hbar_child_rooms"]-delta/3}
            for j, value in proposed.items():
                spec = domain[j]
                if not spec["lower"] <= value <= spec["upper"]:
                    raise RuntimeError("Predeclared child-space pattern leaves existing bounds")
                u[j] = calibration.inverse_unit(value, spec["lower"], spec["upper"], spec["transform"])
            patterns.append((f"child_space_rebalance_{delta:g}" + ("_with_successful_moves" if compensated else ""), u, delta))
    for label, unit, delta in patterns:
        unit = [min(1., max(0., float(u))) for u in unit]
        key = tuple(unit)
        if key in seen:
            continue
        seen.add(key)
        payload = {"old_psi_child": center["old_psi_child"],
                   "best_candidate": {"theta": copy.deepcopy(center["best_candidate"]["theta"]),
                                      "old_psi_child": center["old_psi_child"],
                                      "new_psi_child": center["best_candidate"]["new_psi_child"]}}
        theta = payload["best_candidate"]["theta"]
        for u, spec in zip(unit, domain):
            name = spec["name"]
            value = calibration.transform_unit(u, spec["lower"], spec["upper"], spec["transform"])
            if name == "beta_annual":
                theta["beta"] = value**4
            elif name == "psi_child_change_2023":
                payload["best_candidate"]["new_psi_child"] = payload["old_psi_child"] + value
            else:
                theta[name] = value
        payload["proposal"] = {"label": label, "unit_vector": unit, "delta_first_child_rooms": delta,
            "status": "unevaluated joint candidate; no fitted moments or loss are supplied"}
        add_case(plan, out, payload, label)
    if len(plan["cases"]) > 12:
        raise RuntimeError("Joint case budget exceeded")
    finish_plan(plan, out)


def repeats(prior, status, rows, out):
    # Carry forward the best across *all* completed stages, even when the
    # final joint stage contains no improvement.
    best = adapter.read_json(Path(rows[0]["plan"]).parent.parent / "best_so_far.json")["best"]
    best_path = Path(best["summary"])
    chosen_plan_path = Path(best["plan"])
    chosen_plan = adapter.load_plan(chosen_plan_path, best["plan_sha256"])
    chosen_case = next(c for c in chosen_plan["cases"] if c["id"] == best["id"])
    center = adapter.read_json(chosen_plan_path.parent / chosen_case["center"])
    plan = fresh_plan(prior, out, "final_exact_repeats")
    plan["selected_summary_sha256"] = adapter.digest(best_path)
    reference = out / "selected_reference"
    reference.mkdir()
    summary = adapter.read_json(best_path)
    reference_case = reference / "cases" / summary["best_candidate"]["candidate"]
    reference_case.mkdir(parents=True)
    for name in ("summary.json", "target_fit_long.csv"):
        shutil.copy2(best_path.parent / name, reference / name)
    shutil.copy2(best_path.parent / "cases" / summary["best_candidate"]["candidate"] / "transition_path.csv", reference_case / "transition_path.csv")
    plan["input_sha256"] = {str(p.relative_to(out)): adapter.digest(p) for p in reference.rglob("*") if p.is_file()}
    for i in range(2):
        add_case(plan, out, center, f"final_exact_repeat_{i+1}", panel_id=chosen_case["panel_task_id"],
                 panel_size=chosen_case["panel_size"], panel_design=chosen_case["panel_design"],
                 radius=chosen_case["radius"], panel_seed=chosen_case["panel_seed"])
        plan["cases"][-1].update(reference="selected_reference/summary.json", reference_sha256=adapter.digest(best_path))
    finish_plan(plan, out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("collect", "coordinate", "joint", "repeat"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--expected-planner-sha256", help="Required predeclared source hash for every candidate-generation action")
    args = parser.parse_args()
    if args.action != "collect":
        if not args.expected_planner_sha256:
            raise RuntimeError("Candidate generation requires an explicit planner source hash")
        adapter.verify(__file__, args.expected_planner_sha256)
    prior, status, rows = collect(args.plan.resolve(), args.plan_sha256, require_complete=not args.allow_partial)
    if args.action == "collect":
        print(json.dumps(status))
        return
    if args.allow_partial or not args.outdir:
        raise RuntimeError("Planning requires a complete stage and a new output directory")
    prior = copy.deepcopy(prior)
    prior["parent_plan_sha256"] = args.plan_sha256
    prior["parent_planner_sha256"] = prior.get("planner_sha256")
    {"coordinate": coordinate, "joint": joint, "repeat": repeats}[args.action](prior, status, rows, args.outdir.resolve())


if __name__ == "__main__":
    main()
