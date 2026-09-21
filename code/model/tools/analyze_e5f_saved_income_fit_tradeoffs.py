#!/usr/bin/env python3
"""Audit the saved 96-case income search without running any model solves.

This is a descriptive reader for the terminal overnight search.  It refuses to
continue if completed receipts do not share the frozen target/source contract,
or if their stored objective cannot be reconstructed from target rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SEARCH = ROOT / "output/model/native_financing_diagnostic_20260919/overnight/final_search"
DEFAULT_PLAN = ROOT / "output/model/native_financing_diagnostic_20260919/overnight/plan.remote.json"
DEFAULT_INCUMBENT = ROOT / "output/model/native_financing_diagnostic_20260919/income_search_18047156"
DEFAULT_OUT = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/quantification_v1/saved_fit"

FREE = [
    "beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
    "theta0", "theta1", "first_birth_fixed_cost", "h_P",
]

BLOCKS = {
    "Childless women, ages 40–44": "fertility",
    "Exactly one child among mothers, ages 40–44": "fertility",
    "Period mean first-birth age": "fertility",
    "First births at age 30+": "fertility",
    "Wealth / annual gross labor earnings": "wealth",
    "Annual bequests / aggregate wealth": "wealth",
    "Old wealth/income p90 / median, ages 76–84": "wealth",
    "Mean occupied rooms, capped at 9": "housing",
    "First-birth room response, −1 to +3": "housing",
    "Rooms: 3+ versus 1–2 resident children (model dependent proxy)": "housing",
    "Ownership, heads 30–55": "ownership",
    "Recent-parent ownership gap": "ownership",
}

TARGET_SIGNATURE_FIELDS = (
    "label", "target", "actual_weight", "scored", "role", "restriction_id",
    "empirical_record_id", "evaluation_input", "definition", "sample",
    "empirical_provenance_contract_id", "empirical_source_path", "weight_status",
    "weight_rationale", "working_scale",
)


class ContractError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def finite(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{name} is not numeric: {value!r}") from exc
    if not math.isfinite(out):
        raise ContractError(f"{name} is non-finite")
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def target_signature(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{k: row.get(k) for k in TARGET_SIGNATURE_FIELDS} for row in rows]


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ContractError(f"cannot read JSON {path}: {exc}") from exc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--search-dir", type=Path, default=DEFAULT_SEARCH)
    ap.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    ap.add_argument("--incumbent-dir", type=Path, default=DEFAULT_INCUMBENT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    search = args.search_dir.resolve()
    plan_path = args.plan.resolve()
    incumbent = args.incumbent_dir.resolve()
    out = args.output_dir.resolve()
    cases_path, summary_path, receipt_path = (search / x for x in ("cases.json", "summary.json", "receipt.json"))
    selected_fit_path = search / "selected_target_fit.csv"
    selected_params_path = search / "selected_parameters.csv"
    readout_path = search / "readout.md"
    for p in (cases_path, summary_path, receipt_path, plan_path, selected_fit_path, selected_params_path, readout_path):
        if not p.is_file():
            raise ContractError(f"missing required input: {p}")
    incumbent_receipt_path = incumbent / "receipt.json"
    incumbent_fit_path = incumbent / "selected_target_fit.csv"
    for p in (incumbent_receipt_path, incumbent_fit_path):
        if not p.is_file():
            raise ContractError(f"missing required incumbent comparison input: {p}")
    cases = read_json(cases_path)
    summary = read_json(summary_path)
    receipt = read_json(receipt_path)
    plan = read_json(plan_path)
    if not isinstance(cases, list) or len(cases) != 96:
        raise ContractError(f"expected exactly 96 case entries, got {type(cases).__name__}/{len(cases) if isinstance(cases, list) else 'n/a'}")
    if summary.get("proposal_count") != 96 or receipt.get("proposal_count") != 96:
        raise ContractError("summary/receipt proposal count is not 96")
    actual_bounds = plan.get("parameter_bounds")
    if not isinstance(actual_bounds, dict) or set(actual_bounds) != set(FREE):
        raise ContractError("plan parameter_bounds does not exactly cover the nine free parameters")
    actual_bounds = {k: [finite(v[0], f"lower bound {k}"), finite(v[1], f"upper bound {k}")] for k, v in actual_bounds.items()}
    if actual_bounds["beta_annual"][1] != 0.99:
        raise ContractError("plan beta upper bound is not the actual 0.99 search restriction")

    valid = [x for x in cases if x.get("status") == "completed"]
    invalid = [x for x in cases if x.get("status") != "completed"]
    if len(valid) != 89 or len(invalid) != 7:
        raise ContractError(f"expected 89 completed and 7 rejected cases, got {len(valid)} and {len(invalid)}")
    fingerprints = None
    sig = None
    target_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    case_rows: list[dict[str, Any]] = []
    by_case: dict[int, dict[str, Any]] = {}
    for entry in valid:
        case = int(entry.get("case"))
        r = entry.get("receipt")
        if not isinstance(r, dict) or not isinstance(r.get("target_fit"), list) or len(r["target_fit"]) != 13:
            raise ContractError(f"case {case}: missing 13-row target_fit receipt")
        if not isinstance(r.get("parameters"), list) or len(r["parameters"]) != 17:
            raise ContractError(f"case {case}: missing 17-row parameter receipt")
        if not isinstance(r.get("source_fingerprints"), dict):
            raise ContractError(f"case {case}: missing source_fingerprints")
        this_sig = target_signature(r["target_fit"])
        this_fp = r["source_fingerprints"]
        if sig is None:
            sig, fingerprints = this_sig, this_fp
        elif canonical(this_sig) != canonical(sig):
            raise ContractError(f"case {case}: target definitions/values/weights differ from first completed receipt")
        elif canonical(this_fp) != canonical(fingerprints):
            raise ContractError(f"case {case}: source fingerprints differ from first completed receipt")
        score = 0.0
        block_loss = {b: 0.0 for b in ("fertility", "housing", "ownership", "wealth")}
        abs_raw = {}
        for j, row in enumerate(r["target_fit"]):
            label = row.get("label")
            if label == "Initial model completed fertility":
                if row.get("scored") is not False or row.get("actual_weight") is not None or row.get("loss_contribution") is not None:
                    raise ContractError(f"case {case}: normalization row is incorrectly scored")
                block = "normalization"
                contrib = None
            else:
                if label not in BLOCKS or row.get("scored") is not True:
                    raise ContractError(f"case {case}: unknown or unscored target row {label!r}")
                block = BLOCKS[label]
                gap = finite(row.get("gap"), f"case {case} row {j} gap")
                weight = finite(row.get("actual_weight"), f"case {case} row {j} weight")
                recorded = finite(row.get("loss_contribution"), f"case {case} row {j} contribution")
                contrib = weight * gap * gap
                if not math.isclose(contrib, recorded, rel_tol=2e-11, abs_tol=2e-10):
                    raise ContractError(f"case {case} row {j}: stored loss contribution does not equal weight*gap^2")
                score += contrib
                block_loss[block] += contrib
                abs_raw[label] = abs(gap)
            target_rows.append({"case": case, "case_status": entry.get("status"), "seed_label": entry.get("seed_label"), "direction": entry.get("direction"), "scale": entry.get("scale"), "objective": entry.get("objective"), "target_index": j, "block": block, **{k: row.get(k) for k in ("label", "target", "model", "gap", "actual_weight", "loss_contribution", "scored", "role", "restriction_id", "empirical_record_id", "evaluation_input", "definition", "sample", "weight_status")}})
        reported_loss = finite(r.get("loss"), f"case {case} loss")
        objective = finite(entry.get("objective"), f"case {case} objective")
        if not math.isclose(score, reported_loss, rel_tol=2e-11, abs_tol=2e-8) or not math.isclose(score, objective, rel_tol=2e-11, abs_tol=2e-8):
            raise ContractError(f"case {case}: recomputed total {score} disagrees with receipt/objective ({reported_loss}, {objective})")
        params = {}
        for row in r["parameters"]:
            name = row.get("parameter")
            if name in params:
                raise ContractError(f"case {case}: duplicate parameter {name}")
            est = finite(row.get("estimate"), f"case {case} parameter {name}")
            if name in actual_bounds:
                lo, hi = actual_bounds[name]
                if est < lo - 1e-12 or est > hi + 1e-12:
                    raise ContractError(f"case {case}: {name} outside actual plan bounds")
                near = min(est - lo, hi - est) <= 0.01 * (hi - lo)
                lower, upper = lo, hi
            else:
                near = False
                lower = upper = ""
            params[name] = est
            parameter_rows.append({"case": case, "case_status": entry.get("status"), "seed_label": entry.get("seed_label"), "direction": entry.get("direction"), "scale": entry.get("scale"), "parameter": name, "estimate": est, "actual_lower": lower, "actual_upper": upper, "near_bound_1pct_actual_span": near, "receipt_near_bound": row.get("near_bound"), "structural_coordinate": row.get("structural_coordinate"), "status": row.get("status"), "interpretation": row.get("interpretation"), "transform": row.get("transform")})
        if set(params) != {x["parameter"] for x in r["parameters"]} or not set(FREE).issubset(params):
            raise ContractError(f"case {case}: malformed parameter rows")
        case_row = {"case": case, "seed_label": entry.get("seed_label"), "direction": entry.get("direction"), "scale": entry.get("scale"), "objective": objective, "recomputed_loss": score, **{f"{b}_loss": block_loss[b] for b in block_loss}, "rooms_ownership_joint_loss": block_loss["housing"] + block_loss["ownership"], "fertility_wealth_joint_loss": block_loss["fertility"] + block_loss["wealth"], "rooms_abs_gap": abs_raw["Mean occupied rooms, capped at 9"], "ownership_abs_gap": abs_raw["Ownership, heads 30–55"], "first_birth_rooms_abs_gap": abs_raw["First-birth room response, −1 to +3"], "all_four_block_loss": sum(block_loss.values())}
        case_rows.append(case_row)
        by_case[case] = {"entry": entry, "row": case_row}

    if fingerprints != receipt.get("source_checkpoint_sha256") and not receipt.get("source_checkpoint_sha256"):
        raise ContractError("top-level receipt is missing source checkpoint fingerprint")
    selected_case = int(summary.get("selected", {}).get("case"))
    if selected_case != int(receipt.get("selected_case")) or selected_case not in by_case:
        raise ContractError("selected case is inconsistent or not completed")
    if not math.isclose(by_case[selected_case]["row"]["objective"], finite(receipt.get("selected_loss"), "selected loss"), rel_tol=2e-11, abs_tol=2e-8):
        raise ContractError("selected case loss does not reproduce top-level receipt")
    with selected_fit_path.open(newline="", encoding="utf-8") as f:
        selected_fit_csv = list(csv.DictReader(f))
    with selected_params_path.open(newline="", encoding="utf-8") as f:
        selected_params_csv = list(csv.DictReader(f))
    if len(selected_fit_csv) != 13 or len(selected_params_csv) != 17:
        raise ContractError("saved selected tables do not contain 13 target and 17 parameter rows")
    selected_receipt_rows = by_case[selected_case]["entry"]["receipt"]["target_fit"]
    for i, (a, b) in enumerate(zip(selected_receipt_rows, selected_fit_csv)):
        for field in ("label", "restriction_id", "role", "scored"):
            if str(a.get(field)) != str(b.get(field)):
                raise ContractError(f"selected target row {i}: {field} differs from saved selected table")
        for field in ("target", "model", "gap", "actual_weight", "loss_contribution"):
            av, bv = a.get(field), b.get(field)
            if av is None or bv in (None, ""):
                if av != (None if bv in (None, "") else bv):
                    raise ContractError(f"selected target row {i}: {field} null mismatch")
            elif not math.isclose(finite(av, f"selected receipt {field}"), finite(bv, f"selected table {field}"), rel_tol=2e-10, abs_tol=2e-9):
                raise ContractError(f"selected target row {i}: {field} differs from saved selected table")
    for i, (a, b) in enumerate(zip(by_case[selected_case]["entry"]["receipt"]["parameters"], selected_params_csv)):
        if a.get("parameter") != b.get("parameter") or not math.isclose(finite(a.get("estimate"), "selected estimate"), finite(b.get("estimate"), "selected estimate csv"), rel_tol=2e-10, abs_tol=2e-10):
            raise ContractError(f"selected parameter row {i} differs from saved selected table")
    readout_text = readout_path.read_text(encoding="utf-8")
    if "selected case 60" not in readout_text or "353.6588729140903" not in readout_text:
        raise ContractError("saved readout does not identify selected case 60 and its loss")

    case_rows.sort(key=lambda x: (x["objective"], x["case"]))
    for rank, row in enumerate(case_rows, 1):
        row["total_loss_rank"] = rank
    invalid_rows = [{"case": x.get("case"), "status": x.get("status"), "seed_label": x.get("seed_label"), "direction": x.get("direction"), "scale": x.get("scale"), "objective": x.get("objective"), "reason": x.get("error") or "no verified scored receipt"} for x in invalid]

    # Observed sample nondominance in two raw target-unit deviations.
    frontier = []
    for row in case_rows:
        dominated = any((other["rooms_abs_gap"] <= row["rooms_abs_gap"] and other["ownership_abs_gap"] <= row["ownership_abs_gap"] and (other["rooms_abs_gap"] < row["rooms_abs_gap"] or other["ownership_abs_gap"] < row["ownership_abs_gap"])) for other in case_rows)
        if not dominated:
            frontier.append({**row, "frontier_label": "observed evaluated-search frontier; not attainable frontier"})
    frontier.sort(key=lambda x: (x["rooms_abs_gap"], x["ownership_abs_gap"], x["case"]))
    for rank, row in enumerate(frontier, 1): row["frontier_rank"] = rank

    out.mkdir(parents=True, exist_ok=True)
    target_fields = ["case", "case_status", "seed_label", "direction", "scale", "objective", "target_index", "block", "label", "target", "model", "gap", "actual_weight", "loss_contribution", "scored", "role", "restriction_id", "empirical_record_id", "evaluation_input", "definition", "sample", "weight_status"]
    parameter_fields = ["case", "case_status", "seed_label", "direction", "scale", "parameter", "estimate", "actual_lower", "actual_upper", "near_bound_1pct_actual_span", "receipt_near_bound", "structural_coordinate", "status", "interpretation", "transform"]
    case_fields = list(case_rows[0].keys())
    write_csv(out / "all_valid_target_fit.csv", target_rows, target_fields)
    write_csv(out / "all_valid_parameters.csv", parameter_rows, parameter_fields)
    write_csv(out / "ranked_fit_tradeoffs.csv", case_rows, case_fields)
    write_csv(out / "rejected_cases.csv", invalid_rows, ["case", "status", "seed_label", "direction", "scale", "objective", "reason"])
    write_csv(out / "observed_room_ownership_frontier.csv", frontier, list(frontier[0].keys()) if frontier else case_fields)
    selected_reproduction = {
        "schema": "e5f_saved_selected_case_reproduction_v1",
        "case": selected_case,
        "receipt_loss": by_case[selected_case]["row"]["objective"],
        "selected_table_target_rows": len(selected_fit_csv),
        "selected_table_parameter_rows": len(selected_params_csv),
        "selected_readout_contains_case_and_loss": True,
        "target_table_sha256": sha256(selected_fit_path),
        "parameter_table_sha256": sha256(selected_params_path),
        "readout_sha256": sha256(readout_path),
    }
    (out / "selected_case60_reproduction.json").write_text(json.dumps(selected_reproduction, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    selected_target = [x for x in target_rows if x["case"] == selected_case and x["block"] != "normalization"]
    with incumbent_fit_path.open(newline="", encoding="utf-8") as f:
        incumbent_target = list(csv.DictReader(f))
    incumbent_receipt = read_json(incumbent_receipt_path)
    incumbent_case = int(incumbent_receipt.get("selected_case"))
    incumbent_loss = finite(incumbent_receipt.get("selected_loss"), "incumbent selected loss")
    incumbent_target_scored = [r for r in incumbent_target if str(r.get("scored")).lower() == "true"]
    if len(selected_target) != 12 or len(incumbent_target_scored) != 12:
        raise ContractError("selected/incumbent target comparison does not have 12 scored rows")
    for a, b in zip(selected_target, incumbent_target_scored):
        if a["label"] != b.get("label") or not math.isclose(finite(a["target"], "target"), finite(b.get("target"), "incumbent target"), rel_tol=1e-10, abs_tol=1e-12) or not math.isclose(finite(a["actual_weight"], "weight"), finite(b.get("actual_weight"), "incumbent weight"), rel_tol=1e-10, abs_tol=1e-9):
            raise ContractError("incumbent target definitions/values/weights do not match overnight search")
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ContractError("matplotlib is required to make the supplemental figure") from exc
    labels = [r["label"] for r in selected_target]
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    y = np.arange(12)
    axes[0].barh(y + 0.18, [float(r["gap"]) for r in incumbent_target_scored], height=0.34, label=f"retained incumbent (prior case {incumbent_case})", color="#9ecae1")
    axes[0].barh(y - 0.18, [float(r["gap"]) for r in selected_target], height=0.34, label=f"overnight selected case {selected_case}", color="#2171b5")
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_yticks(y, labels, fontsize=7)
    axes[0].set_xlabel("Model minus target (target units)")
    axes[0].set_title("Scored target deviations")
    axes[0].legend(fontsize=8)
    xs = [r["rooms_abs_gap"] for r in case_rows]; ys = [r["ownership_abs_gap"] for r in case_rows]
    axes[1].scatter(xs, ys, s=20, alpha=0.35, color="#636363", label="89 evaluated cases")
    axes[1].scatter([by_case[selected_case]["row"]["rooms_abs_gap"]], [by_case[selected_case]["row"]["ownership_abs_gap"]], s=65, color="#d62728", label=f"selected {selected_case}")
    axes[1].scatter([r["rooms_abs_gap"] for r in frontier], [r["ownership_abs_gap"] for r in frontier], s=32, facecolors="none", edgecolors="#ff7f0e", label="observed evaluated-search frontier")
    axes[1].set_xlabel("Absolute rooms gap (rooms)")
    axes[1].set_ylabel("Absolute ownership gap (share)")
    axes[1].set_title("Observed room/ownership deviations")
    axes[1].legend(fontsize=7)
    fig.savefig(out / "saved_fit_tradeoffs.png", dpi=190)
    plt.close(fig)

    readme = f"""# Saved income-fit tradeoffs

This deterministic audit reads the 96 saved proposals from `overnight/final_search`: {len(valid)} completed receipts and {len(invalid)} rejected proposals. It runs no household or equilibrium solves. The completed receipts agree exactly on the 13 target rows (12 scored plus the separate 2.1 normalization), the nine free-parameter bounds from `plan.remote.json` (with the actual beta upper bound 0.99), and the source-fingerprint dictionary. Every scored contribution is independently recomputed as `actual_weight * gap^2`; each total reproduces the saved objective.

Reproduce with: `code/model/.venv/bin/python code/model/tools/analyze_e5f_saved_income_fit_tradeoffs.py`.

The ranked table sums losses by the explicit blocks `fertility` (four rows), `housing` (three), `ownership` (two), and `wealth` (three). Near-bound flags use a distance of at most 1% of the actual plan-bound span. The room/ownership file is a nondominated set under absolute target-unit gaps among these evaluated cases. It is an observed evaluated-search frontier only; it does not establish an attainable frontier, identification, derivatives, or causal effects.

Five factual findings from the saved search:

1. The overnight selected case {selected_case} has total loss {by_case[selected_case]['row']['objective']:.12g}; the retained incumbent from the prior search (case {incumbent_case}) has loss {incumbent_loss:.12g}.
2. The selected case ranks {by_case[selected_case]['row']['total_loss_rank']} of {len(case_rows)} by the unchanged weighted objective.
3. The selected case's largest block is {max((b for b in ('fertility','housing','ownership','wealth')), key=lambda b: by_case[selected_case]['row'][f'{b}_loss'])}, with loss {max(by_case[selected_case]['row'][f'{b}_loss'] for b in ('fertility','housing','ownership','wealth')):.6g}; this is a descriptive decomposition, not a new score.
4. The observed room/ownership nondominated set contains {len(frontier)} evaluated cases; its coordinates are raw absolute target-unit deviations.
5. Seven proposals are retained separately as rejected cases and are excluded from fit rankings and frontier construction.

Limits: this is a finite, adaptive proposal sample. It cannot establish global optimization, local identification, an attainable tradeoff curve, or causal comparative statics. The target rows retain the project's existing model-observer and empirical-sample qualifications; no target, weight, parameter bound, or objective was changed.
"""
    (out / "README.md").write_text(readme, encoding="utf-8")

    outputs = [p for p in out.iterdir() if p.is_file() and p.name != "analysis_receipt.json"]
    analysis_receipt = {
        "schema": "e5f_saved_income_fit_tradeoffs_v1",
        "status": "verified",
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (cases_path, summary_path, receipt_path, plan_path, incumbent_receipt_path, incumbent_fit_path)},
        "script_sha256": sha256(Path(__file__).resolve()),
        "case_counts": {"total": 96, "valid": len(valid), "rejected": len(invalid)},
        "selected_case": selected_case,
        "incumbent_case": incumbent_case,
        "incumbent_loss": incumbent_loss,
        "target_signature_sha256": hashlib.sha256(canonical(sig).encode()).hexdigest(),
        "source_fingerprints_sha256": hashlib.sha256(canonical(fingerprints).encode()).hexdigest(),
        "source_fingerprints": fingerprints,
        "actual_bounds": actual_bounds,
        "near_bound_rule": "min(estimate-lower, upper-estimate) <= 0.01*(upper-lower)",
        "block_assignment": BLOCKS,
        "recomputed_total_loss_max_abs_error": max(abs(r["objective"] - r["recomputed_loss"]) for r in case_rows),
        "output_sha256": {p.name: sha256(p) for p in sorted(outputs)},
    }
    (out / "analysis_receipt.json").write_text(json.dumps(analysis_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "verified", "valid": len(valid), "rejected": len(invalid), "selected_case": selected_case, "incumbent_case": incumbent_case, "frontier_cases": len(frontier), "output_dir": str(out)}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except ContractError as exc:
        raise SystemExit(f"CONTRACT ERROR: {exc}")
