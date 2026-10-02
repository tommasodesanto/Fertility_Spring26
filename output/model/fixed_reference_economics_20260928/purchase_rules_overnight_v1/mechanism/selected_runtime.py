"""Authenticate a selected purchase-rule fit before using dated model callbacks."""
from __future__ import annotations

import hashlib
import gzip
import importlib.util
import json
import pickle
import re
import sys
from pathlib import Path
from types import SimpleNamespace, MethodType
import numpy as np

PACKET = Path(__file__).resolve().parent.parent
ROOT = PACKET.parents[3]
OLD = PACKET.parent / "utility_floor_round2_v1"
TRANSITION = ROOT / "code/model/experiments/transition_readiness"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def validate_compact_repeat(report, repeat):
    """Authenticate the saved compact native repeat without demanding full reports."""
    require((repeat / "stage/solution_arrays.npz").is_file(),
            "Selected native state arrays missing")
    compact_closure = read(repeat / "closure.json")
    full_closure = read(report / "closure.json")
    require(bool(compact_closure) and all(full_closure.get(key) == value
                                          for key, value in compact_closure.items()),
            "Saved selected repeat closure differs from selected root")


def renderer_identity(report):
    """Read the renderer written inside every standard PNG artifact."""
    from PIL import Image
    plots = sorted((Path(report) / "standard_diagnostics").glob("*.png"))
    require(len(plots) == 17, "Standard diagnostic plot count differs")
    versions = set()
    for plot in plots:
        with Image.open(plot) as image:
            marker = image.info.get("Software")
        require(isinstance(marker, str) and marker.startswith("Matplotlib version"),
                "Standard plot renderer metadata missing: " + plot.name)
        versions.add(marker)
    require(len(versions) == 1, "Mixed renderer identities within standard plots")
    return next(iter(versions)), [plot.name for plot in plots]


def full_native_arrays(packet):
    """Recreate the 95 saved postcheck arrays from the fresh native solution."""
    sol, shared, P, ev = (packet[key] for key in ("solution", "shared", "parameters", "evaluation"))
    values = {key: value for key, value in vars(sol).items()
              if isinstance(value, np.ndarray) and value.dtype != object}
    values.update({"shared." + key: value for key, value in vars(shared).items()
                   if isinstance(value, np.ndarray) and value.dtype != object})
    for key in ("g_pre", "g_post_fertility", "g_current", "g_stay_distribution"):
        values["distribution." + key] = getattr(ev, key)
    demand = np.asarray(ev.demand_by_loc)
    supply = np.asarray(ev.supply_by_loc)
    values["market_housing_demand"] = demand
    values["market_housing_excess"] = demand - supply
    values["parameters.H0"] = np.asarray(P.H0)
    values["normalization.population"] = np.array([1.])
    return values


def compare_full_native_arrays(saved_path, packet, tolerance=1e-10):
    fresh = full_native_arrays(packet)
    with np.load(saved_path, allow_pickle=False) as saved:
        require(len(saved.files) == len(fresh) == 95 and set(saved.files) == set(fresh),
                "Saved native array key set differs")
        maximum = 0.
        for key in sorted(saved.files):
            old, new = saved[key], np.asarray(fresh[key])
            require(old.shape == new.shape and old.dtype == new.dtype,
                    "Saved native array shape/dtype differs: " + key)
            require(old.dtype.kind in "biufc", "Unsupported saved native array dtype: " + key)
            if old.dtype.kind in "biu":
                require(np.array_equal(old, new), "Saved discrete native array differs: " + key)
                continue
            for mask in (np.isnan, np.isposinf, np.isneginf):
                require(np.array_equal(mask(old), mask(new)),
                        "Saved native array nonfinite mask differs: " + key)
            finite = np.isfinite(old)
            delta = float(np.max(np.abs(old[finite] - new[finite]))) if finite.any() else 0.
            require(delta <= tolerance, "Saved native numeric array differs: " + key)
            maximum = max(maximum, delta)
    return {"fields": len(fresh), "max_abs_delta": maximum, "numeric_tolerance": tolerance}


def compare_renderer_aware_reference(runtime, reference, current, packet, saved_arrays):
    """Preserve the full existing CSV/closure check; isolate renderer-only PNG drift."""
    from floor_runtime import compare_normalized_reference_reports
    old_renderer, old_names = renderer_identity(reference)
    new_renderer, new_names = renderer_identity(current)
    require(old_names == new_names, "Standard plot names differ")
    arrays = compare_full_native_arrays(saved_arrays, packet)
    if old_renderer == new_renderer:
        comparison = compare_normalized_reference_reports(reference, current)
        reason = "same_renderer_exact_pngs"
    else:
        try:
            comparison = compare_normalized_reference_reports(reference, current)
            reason = "different_renderer_exact_pngs"
        except RuntimeError as exc:
            require(str(exc) == "Normalized reference standard plots differ",
                    "Full native reference report differs: " + str(exc))
            comparison = {"status": "renderer_only_plot_bytes_differ",
                          "target_rows": 14, "parameter_rows": 31,
                          "standard_plot_names": old_names}
            reason = "cross_renderer_full_arrays_and_reports"
    return dict(status="passed", reason=reason, saved_renderer=old_renderer,
                fresh_renderer=new_renderer, arrays=arrays, report=comparison,
                standard_plot_names=old_names)


def reconstruct_purchase_reference(self, folder):
    """Exact native double replay with the isolated renderer-aware comparison."""
    from floor_runtime import normalized_housing_contract, write
    folder = Path(folder)
    results = []
    receipts = []
    saved_arrays = self.selected_repeat / "stage/solution_arrays.npz"
    for index in range(2):
        packet, record = self.stationary(self.P.psi_child, self.reference_price,
                                         folder / f"repeat_{index}")
        require(abs(record["renewal_residual"]) <= 1e-6, "Selected renewal differs")
        require(abs(record["population_scale"] / self.population_scale - 1) <= 1e-10,
                "Selected population differs")
        live = dict(P=packet["parameters"], b_grid=self.grid, sd=packet["shared"],
                    sol=packet["solution"], price=np.array([self.reference_price]))
        context = dict(self.ctx, out=folder / f"repeat_{index}", deadline_epoch=float("inf"))
        with self.native_bindings():
            self.ge.observe_price(context, live, "selected_root", final=True)
        report = folder / f"repeat_{index}/phase_b_ge/selected_root"
        require(normalized_housing_contract(self.handoff), "Selected normalized housing contract missing")
        receipts.append(compare_renderer_aware_reference(self, self.report, report, packet, saved_arrays))
        results.append((packet, record))
    self.runner.compare_repeated(folder / "repeat_0/phase_b_ge/selected_root",
                                 folder / "repeat_1/phase_b_ge/selected_root")
    self.packet = results[-1][0]
    self.initial_state = self.stationary_state(self.packet, self.population_scale)
    self.reference_verified = True
    checkpoint = folder / "selected_native_packet.pkl.gz"
    with gzip.open(checkpoint, "wb") as stream:
        pickle.dump(self.packet, stream)
    write(folder / "reference_reconstruction.json", dict(
        schema="current_floor_reference_reconstruction_v1", status="passed",
        policy_calls=sum(row[1]["policy_calls"] for row in results), lifecycle_calls=2,
        checkpoint=dict(path=str(checkpoint), sha256=sha(checkpoint)),
        checkpoint_sha256=sha(checkpoint), identity=self.identity(),
        renderer_aware_replays=receipts,
        reports={str(report): {str(path.relative_to(report)): sha(path)
                              for path in report.rglob("*") if path.is_file() and
                              (path.name in ("parameters.csv", "target_fit.csv", "closure.json")
                               or path.suffix == ".png")}
                 for report in (folder / "repeat_0/phase_b_ge/selected_root",
                                folder / "repeat_1/phase_b_ge/selected_root")}))
    return self.packet


def authenticate_selected(arm, completed):
    """Require an actual fresh selected postcheck and matching 14/31/17 reports."""
    require(arm in {"hard", "quarter"}, "Unknown purchase rule")
    completed = Path(completed).resolve()
    receipt = read(completed)
    require(receipt.get("status") == "selected_numerically_verified", "No passed selected postcheck")
    postcheck = receipt.get("selected_postcheck", {})
    require(postcheck.get("status") == "passed", "Selected native postcheck failed")
    require(receipt.get("selected", {}).get("parameters") == postcheck.get("parameters", receipt.get("selected", {}).get("parameters")),
            "Selected and postchecked coordinates differ")
    chain = completed.parent.parent
    require(re.fullmatch(r"chain_?\d+", chain.name) is not None,
            "Selected source is not a named calibration chain")
    root = completed.parent / "selected_postcheck/phase_b_ge/selected_root"
    repeat = completed.parent / "selected_postcheck/phase_b_ge/selected_repeat"
    require(root.is_dir() and repeat.is_dir(), "Selected native ROOT/REPEAT missing")
    validate_compact_repeat(root, repeat)
    input_contract = read(completed.parent / "input_contract.json")
    require(input_contract.get("purchase_rule") == arm and input_contract.get("normalized_population") == 1.0
            and input_contract.get("owner_financed_share") == 0.8
            and input_contract.get("entry", {}).get("arm") == "nonnegative_mean",
            "Selected purchase/entry/population contract differs")
    require(receipt["selected"].get("weight_fingerprint") == input_contract.get("weight_contract_sha256")
            and postcheck.get("weight_fingerprint") == input_contract.get("weight_contract_sha256")
            and input_contract.get("parameter_contract_sha256"),
            "Selected target-and-weight fingerprint differs")
    return receipt, input_contract, root, repeat


def construct(arm, completed, output):
    """Build the existing native floor runtime with this arm's isolated engine."""
    receipt, contract, report, repeat = authenticate_selected(arm, completed)
    for rel, digest in read(PACKET / "engine_pins.json").items():
        require(sha(PACKET / rel) == digest, "Isolated engine source drift: " + rel)
    for rel, digest in read(PACKET / "source_pins.json").items():
        require(sha(ROOT / rel) == digest, "Pinned integration source drift: " + rel)
    engine_root = PACKET / "engines" / arm
    sys.path.insert(0, str(engine_root))
    from small_credit_lab.engine import solver
    from small_credit_lab import credit
    from refactor_lab.engine import solver as checked_solver
    require(Path(solver.__file__).resolve().is_relative_to(engine_root)
            and Path(checked_solver.__file__).resolve().is_relative_to(engine_root),
            "Purchase-rule engine import escaped selected arm")
    sys.path[:0] = [str(OLD), str(ROOT / "code/model/tools"), str(ROOT / "code/model"), str(TRANSITION)]
    import inputs
    import runner
    from floor_runtime import FloorRuntime, normalized_housing_contract, bind_normalized_reference_housing
    # Saved selected_repeat is compact; reconstruct_reference compares two fresh
    # full 14/31/17 reports against selected_root and against each other.
    rows = runner.readtable(report / "target_fit.csv")
    runner.residual(rows)
    parameter_rows = runner.readtable(report / "parameters.csv")
    require(len(parameter_rows) == 31, "Selected complete parameter table missing")
    actual = {row["parameter"]: float(row["estimate"]) for row in parameter_rows}
    closure = read(report / "closure.json")
    require(float(closure["normalized_population"]) == 1.0 and float(closure["population_scale"]) == 1.0,
            "Selected initial population differs")
    require(float(closure["H0_derived"]) == actual["H0"], "Selected H0 differs")
    bounds = contract["bounds"]
    point = receipt["selected"]["parameters"]
    require(set(point) == set(contract["parameter_contract"]["free_coordinates"]), "Selected coordinate set differs")
    require(all(abs(actual[k] - float(v)) <= 1e-10 for k, v in point.items()),
            "Selected parameter report differs")
    lane = contract["lane"]
    P, grid = inputs.proposal(lane)
    P, entry = inputs.entry(P, grid, "nonnegative_mean")
    require(np.asarray(P.phi).shape == (4,), "Financed-share family dimensions differ")
    P.phi = np.full_like(np.asarray(P.phi, dtype=float), 0.8)
    P.experimental_purchase_saving_fraction = 0.25 if arm == "quarter" else 1.0
    P = inputs.bind(P, point, bounds, "floor")
    handoff = {"normalized_housing_contract": dict(schema="normalized_n0_fixed_h0_v1", N0=1.0,
              H0_source="authenticated_parameter_table", counterfactual_H0_fixed=True)}
    require(normalized_housing_contract(handoff), "Normalized fixed-H0 contract missing")
    P = bind_normalized_reference_housing(P, handoff, closure, actual)
    credit.bind_engine_credit(P, "corrected", 0.0)
    sys.path.insert(0, str(runner.BASE))
    base = load("purchase_base_workflow", runner.BASE / "run_comparison.py")
    ge = load("purchase_native_observer", OLD.parent / "utility_calibration_round1_v1" / "phase_b_pilot.py")
    runner.install_reporter_on_authored(base.authored)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    ctx = base.authored.context_from_bundle(SimpleNamespace(
        bundle=ROOT / "output/model/publication_refactor_20260929/local_export_v1/inputs",
        reference_root=ROOT, out=output))
    base.authored.authenticate_frozen(ctx)
    ctx.update(P=P, b_grid=grid, selected_d_bar=0.0,
               reference_psi=float(P.psi_child), expected_parameters=actual,
               free_coordinates=list(point), out=output, deadline_epoch=float("inf"))
    runner.install_observer_metadata(ge, "floor", bounds)
    for row in ctx["manifest"]["full_parameter_table"]:
        if row["parameter"] in bounds:
            row["lower"], row["upper"] = map(str, bounds[row["parameter"]])
    rt = FloorRuntime()
    rt.handoff = handoff
    rt.handoff_pin = {"sha256": sha(completed)}
    rt.sourcepins = {**read(PACKET / "source_pins.json"), **{str(PACKET.relative_to(ROOT) / k): v
                    for k, v in read(PACKET / "engine_pins.json").items() if k.startswith("engines/" + arm + "/")}}
    rt.report = report
    rt.selected_repeat = repeat
    rt.reconstruct_reference = MethodType(reconstruct_purchase_reference, rt)
    rt.runner = runner
    rt.ge = ge
    rt.ctx = ctx
    rt.P = P
    rt.grid = grid
    rt.model = solver
    rt.rt = ctx["prepared"].rt
    rt.pf = rt.rt["primitive"].pf
    rt.parameter_rows = parameter_rows
    rt.reference_price = float(closure["price"])
    rt.population_scale = 1.0
    rt.packet = None
    rt.initial_state = None
    rt.reference_verified = False
    rt.total_native_calls = 0
    rt.folder = output
    rt.scaffold = load("purchase_dated_scaffold", TRANSITION / "pinned_tools/run_e5f_preference_transition.py")
    rt.install_observer_adapters(output / "observer_adapters")
    with rt.native_bindings():
        got = ctx["fp"].actual_parameters(ctx["prepared"], P, grid)
        ge.validate_parameter_estimates(dict(expected_parameters=actual), parameter_rows, got)
        cohort = rt.pf.calendar.entrant_cohort(np.array([1.0]), P, grid)
        np.testing.assert_allclose(cohort.sum(axis=(1, 2, 4, 5)),
                                   P.fixed_reference_entry_conditional * P.z_weights[None, :],
                                   rtol=0, atol=2e-16)
    return rt
