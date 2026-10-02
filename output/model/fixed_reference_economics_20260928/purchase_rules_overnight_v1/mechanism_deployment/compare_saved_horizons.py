"""Read-only overlap of two accepted dated paths and their saved period-4 states.

The pinned macro helper is invoked unchanged. The saved-state portion repeats
only the gate-field formulas in one_shock_floor.compare_2023_states: saved
checkpoints lack the value arrays and forecast paths needed for its full call.
"""
from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np


EXPECTED = {
    "pinned_tools/run_e5f_preference_transition.py": "26802381cec56be33794a526eeec80eb6066c09ae59ab6e5957d0cee2e26f1ac",
    "one_shock_floor.py": "a6499da273ff3a607dd741a283e6ba8415e5c9a4c161a18d1b1eab07cce900a0",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-root", type=Path, required=True)
    ap.add_argument("--frozen-root", type=Path, required=True)
    ap.add_argument("--short-case", type=Path, required=True)
    ap.add_argument("--long-case", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    require(not args.out.exists(), "Refusing to replace an existing overlap report")

    transition = args.source_root / "source/code/model/experiments/transition_readiness"
    inventory = json.loads((args.source_root / "inventory.json").read_text())
    source_shas = {}
    for rel, expected in EXPECTED.items():
        path = transition / rel
        digest = sha(path)
        inventory_rel = "code/model/experiments/transition_readiness/" + rel
        require(digest == expected == inventory["files"][inventory_rel], "Pinned helper identity changed: " + rel)
        source_shas[rel] = digest

    native_source_checks = {}
    for rel in (
        "code/model/tools/run_e5f_perfect_foresight_transition.py",
        "code/model/tools/run_e5f_open_population_transition.py",
        "code/model/tools/run_dynamic_population_transition.py",
        "code/model/intergen_eqscale_seq_optimized/adult_entry.py",
    ):
        frozen_path = args.frozen_root / rel
        digest = sha(frozen_path)
        for overlay in (Path("/scratch/td2248/projects/grid_resolution_credit053_v2"),
                        Path("/scratch/td2248/projects/normalized_floor_calibration_v1")):
            require(not (overlay / "source" / rel).exists(), "Mounted overlay replaces frozen helper: " + rel)
        require(not (args.source_root / "source" / rel).exists(), "Extension overlay replaces frozen helper: " + rel)
        receipts = {}
        for case in (args.short_case, args.long_case):
            pin_file = case / "run/runtime/runtime_auth/runtime_preparation/native_preparation/preparation.json"
            pinned = json.loads(pin_file.read_text())["current_source_files"]
            require(pinned.get(rel) == digest, "Saved native preparation does not pin frozen helper: " + str(case) + ": " + rel)
            receipts[case.name] = pinned[rel]
        native_source_checks[rel] = dict(path=str(frozen_path), sha256=digest, production_receipts=receipts)

    state_tree = ast.parse((transition / "one_shock_floor.py").read_text())
    gates_node = next(n for n in state_tree.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "GATES" for t in n.targets))
    gates = eval(compile(ast.Expression(gates_node.value), str(transition / "one_shock_floor.py"), "eval"), {"dict": dict})
    tolerance = gates["horizon_relative_tolerance"]
    require(tolerance == 1e-3, "Original horizon tolerance changed")

    helper = transition / "pinned_tools/run_e5f_preference_transition.py"
    helper_tree = ast.parse(helper.read_text())
    node = next(n for n in helper_tree.body if isinstance(n, ast.FunctionDef) and n.name == "compare_horizons")
    scope = {"math": math, "require": require}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(helper), "exec"), scope)
    compare_horizons = scope["compare_horizons"]

    cases = (args.short_case, args.long_case)
    completed = [json.loads((p / "run/dated_path/completed.json").read_text()) for p in cases]
    contracts = [json.loads((p / "run/run_contract.json").read_text()) for p in cases]
    for p, d in zip(cases, completed):
        require(json.loads((p / "launcher_terminal.json").read_text())["exit_code"] == 0,
                "Case launcher did not exit successfully")
        require(d["status"] == "passed" and d["root"]["converged"]
                and all(d["root"]["gates"].values()) and d["terminal"]["all_checks_pass"]
                and d["terminal"]["raw_queue_pass"], "Case is not accepted")
    require(completed[0]["reference_identity"] == completed[1]["reference_identity"], "Reference identity changed")
    require(completed[0]["phi_path"] == completed[1]["phi_path"][:completed[0]["horizon"]], "Early financing path changed")
    require(contracts[0]["selected_sha256"] == contracts[1]["selected_sha256"]
            and contracts[0]["fixed_H0"] is contracts[1]["fixed_H0"] is True, "Housing contract changed")
    require(contracts[0]["arm"] == contracts[1]["arm"]
            and contracts[0]["kind"] == contracts[1]["kind"], "Policy contract changed")
    identity = completed[0]["reference_identity"]
    wrappers = []
    for d, c in zip(completed, contracts):
        wrappers.append(dict(
            horizon=d["horizon"], reference_manifest_sha256=identity["reference_sha256"],
            source_pins=identity["source_pins"],
            housing=(c["selected_sha256"], c["fixed_H0"]),
            shock_contract=(c["arm"], c["kind"], tuple(d["phi_path"][:completed[0]["horizon"]])),
            root_and_terminal_pass=True, rows=d["rows"],
        ))
    macro = compare_horizons(*wrappers, 4, tolerance)

    sys.path.insert(0, str(args.frozen_root / "code/model/tools"))
    import run_e5f_perfect_foresight_transition as pf

    states, checkpoint_hashes = [], []
    for p, d in zip(cases, completed):
        path = p / "run/dated_path/accepted_period4_state.pkl.gz"
        digest = sha(path)
        require(digest == d["horizon_overlap_period4_state"]["sha256"], "Saved period-4 state hash changed")
        with gzip.open(path, "rb") as stream:
            states.append(pickle.load(stream)["state"])
        checkpoint_hashes.append(digest)
    g = [np.asarray(s.g_pre, float) for s in states]
    require(g[0].shape == g[1].shape and all(np.isfinite(v).all() and (v >= 0).all() for v in g),
            "Saved period-4 distributions invalid")
    population = [float(v.sum()) for v in g]
    require(min(population) > 0, "Saved period-4 population invalid")
    state_gaps = dict(
        population_absolute_gap=abs(population[0] - population[1]),
        population_relative_gap=abs(population[0] - population[1]) / max(population),
        normalized_distribution_l1=float(np.abs(g[0] / population[0] - g[1] / population[1]).sum()),
    )
    for name in ("scheduled_entries", "scheduled_raw_entries"):
        x, y = [np.asarray(pf.birth_queue_values(getattr(s, name)), float) for s in states]
        require(x.shape == y.shape and np.isfinite(x).all() and np.isfinite(y).all(), "Saved queue invalid")
        state_gaps[name + "_relative_l1"] = float(np.abs(x - y).sum() /
                                                max(float(np.abs(x).sum()), float(np.abs(y).sum()), 1e-12))
    state_pass = all(v <= 1e-3 for k, v in state_gaps.items() if k != "population_absolute_gap")
    flow_gaps = [abs(a["flow"] - b["flow"]) for a, b in
                 zip(completed[0]["first_births"][:4], completed[1]["first_births"][:4])]
    report = dict(
        schema="saved_accepted_horizon_overlap_v1", source_shas=source_shas,
        native_source_checks=native_source_checks,
        completed_paths=[str(p / "run/dated_path/completed.json") for p in cases],
        verified_period4_state_hashes=checkpoint_hashes,
        macro_original_helper="pinned_tools/run_e5f_preference_transition.py:compare_horizons",
        original_usage_periods=4, original_horizon_tolerance=tolerance, macro=macro,
        saved_state_original_comparator_gate_fields=state_gaps,
        saved_state_gate_fields_pass=state_pass, saved_state_gate_tolerance=1e-3,
        full_original_state_comparator_available=False,
        unavailable_fields=["current_2023_V", "continuation_2027_V", "forecast_prices",
                            "forecast_pensions", "forecast_psi"],
        first_birth_flow_first4_absolute_gaps=flow_gaps,
        first_birth_flow_first4_max_absolute_gap=max(flow_gaps),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        stream.write(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
