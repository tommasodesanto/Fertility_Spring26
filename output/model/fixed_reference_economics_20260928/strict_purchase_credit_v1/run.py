"""Three fresh fixed-price strict-purchase lifecycle solves: 80, 90, 90 repeat."""
from __future__ import annotations

import csv
import difflib
import hashlib
import inspect
import json
import signal
import sys
import time
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STRICT = HERE / "input/strict80"
SANDBOX = ROOT / "code/model/experiments/strict_purchase_sandbox/source"
sys.path.insert(0, str(SANDBOX))
from small_credit_lab import credit  # noqa: E402
from small_credit_lab.engine import solver  # noqa: E402
V2 = HERE.parent / "normalized_calibration_v2"
sys.path.insert(0, str(V2))
import run_psi as v2  # noqa: E402

PRICE = 0.7152515073815459
H0 = 6.778473404808042
LIMIT = 1200
LANE = "floor_s0"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_once(source, before, after):
    assert source.count(before) == 1, before[:100]
    return source.replace(before, after, 1)


def variant_observer(source):
    source = replace_once(source,
        "    population = physical_supply / demand\n",
        "    population = 1.0  # Fixed physical population, not housing-clearing scale.\n")
    source = replace_once(source,
        "        _require(abs(renewal) <= RENEWAL_TOL, \"Birth renewal root fails\")\n"
        "        result[\"native_population_step\"] = native_population_step(\n"
        "            context, packet, population, actual_births, entry)\n",
        "        result[\"native_population_step\"] = native_population_step(\n"
        "            context, packet, population, actual_births, entry)\n"
        "        result[\"renewal_status\"] = \"measured_fixed_price_diagnostic\"\n"
        "        result[\"housing_market_status\"] = \"measured_fixed_H0_diagnostic\"\n")
    source = replace_once(source,
        '            if key in context["free_coordinates"]:\n',
        '            if key in context["fixed_coordinates"]:\n'
        '                row["status"] = "fixed at strict-80 reference estimate for diagnostic"\n'
        '            elif key in context["free_coordinates"]:\n')
    source = replace_once(source,
        '                row["status"] = "fixed population scale; not identified by per-household targets"\n',
        '                row["status"] = "fixed strict-80 physical supply coefficient for diagnostic"\n')
    source = replace_once(source,
        '            elif key == "psi_child":\n',
        '            elif key == "financed_share":\n'
        '                row["status"] = "experimental fixed policy input"\n'
        '            elif key == "psi_child":\n')
    source = replace_once(source,
        '        result["standard_plot_supply_units"] = "physical supply divided by endogenous household population"\n',
        '        result["standard_plot_supply_units"] = "physical supply and demand at fixed N0=1"\n')
    return source


FIXED_PHASE = '''def run_phase_b(context, phase_a_result, budget):
    """One fresh lifecycle solve at the pinned price; no renewal or housing root."""
    _observer_context(context)
    d_bar = float(phase_a_result["selected_d_bar"])
    _require(d_bar >= 0 and math.isfinite(d_bar), "Invalid credit input")
    context["selected_d_bar"] = d_bar
    context["reference_psi"] = float(context["P"].psi_child)
    context["deadline_epoch"] = float(budget.deadline_epoch)
    q = float(context["fixed_price"])
    _require(q > 0 and math.isfinite(q), "Invalid fixed price")
    live = solve_fixed_price(context, d_bar, q, budget, "selected_root",
        Path(context["out"]) / "phase_b_ge" / "selected_root" / "stage")
    observed = _observe_with_deadline(context, live, "selected_root", final=True)
    return dict(status="passed_fixed_price_diagnostic", selected_price=q,
                selected=observed, lifecycle_solves=int(budget.used_lifecycle))
'''


def make_ge(out):
    native = v2.native
    native.verify_sources()
    for rel, digest in json.loads((V2 / "source_pins.json").read_text()).items():
        assert sha(ROOT / rel) == digest, rel
    sys.path.insert(0, str(native.BASE))
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("strict_credit_initializer", native.BASE / "run_comparison.py")
    initializer = module_from_spec(spec)
    spec.loader.exec_module(initializer)
    ge_path = native.HERE.parent / "utility_calibration_round1_v1/phase_b_pilot.py"
    sys.path.insert(0, str(ge_path.parent))
    spec = spec_from_file_location("strict_credit_original_ge", ge_path)
    original = module_from_spec(spec)
    spec.loader.exec_module(original)
    original_observer = inspect.getsource(original.observe_price)
    original_phase = inspect.getsource(original.run_phase_b)
    observer = variant_observer(original_observer)
    review = out / "source_review"
    review.mkdir()
    for name, old, new in (("observe_price", original_observer, observer),
                           ("run_phase_b", original_phase, FIXED_PHASE)):
        (review / (name + ".py")).write_text(new)
        (review / (name + ".diff")).write_text("".join(difflib.unified_diff(
            old.splitlines(True), new.splitlines(True), fromfile="native/" + name,
            tofile="fixed_price/" + name)))
    ge = ModuleType("strict_credit_fixed_price_ge")
    ge.__dict__.update(original.__dict__)
    exec(compile(observer, str(review / "observe_price.py"), "exec"), ge.__dict__)
    exec(compile(FIXED_PHASE, str(review / "run_phase_b.py"), "exec"), ge.__dict__)
    wrapper = inspect.getsource(original._observe_with_deadline)
    (review / "_observe_with_deadline.py").write_text(wrapper)
    exec(compile(wrapper, str(review / "_observe_with_deadline.py"), "exec"), ge.__dict__)
    v2.write(review / "receipt.json", dict(
        original_observer_sha256=hashlib.sha256(original_observer.encode()).hexdigest(),
        variant_observer_sha256=hashlib.sha256(observer.encode()).hexdigest(),
        original_phase_sha256=hashlib.sha256(original_phase.encode()).hexdigest(),
        variant_phase_sha256=hashlib.sha256(FIXED_PHASE.encode()).hexdigest(),
        fixed_price=PRICE, fixed_H0=H0, fixed_population=1.0,
        no_renewal_or_housing_root=True, unchanged_household_and_kernels=True))
    return initializer, ge


def read_table(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), "Refusing existing output"
    args.out.mkdir(parents=True)
    start = time.time()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("20-minute cap")))
    signal.setitimer(signal.ITIMER_REAL, LIMIT)
    try:
        strict = json.loads((STRICT / "collected/run/completed.json").read_text())
        incumbent = json.loads((STRICT / "incumbent.json").read_text())
        assert strict["status"] == "fixed_coordinate_strict_purchase_experiment_passed"
        assert strict["target_fingerprint"] == incumbent["target_fingerprint"]
        assert strict["weight_fingerprint"] == incumbent["weight_fingerprint"]
        assert abs(strict["price"] - PRICE) < 1e-12
        assert abs(strict["H0_derived"] - H0) < 1e-12
        manifest = json.loads((STRICT / "manifest.json").read_text())
        for pair, pins in manifest["source_pairs"].items():
            original, sandbox = (ROOT / rel for rel in pair.split("|"))
            assert sha(original) == pins["original"] and sha(sandbox) == pins["sandbox"]
        assert sha(STRICT / "incumbent.json") == manifest["incumbent_sha256"]
        point = incumbent["parameters"]
        assert len(point) == 10 and strict["all_ten_coordinates_fixed"]
        seed, bounds, _ = v2.inputs.seed_and_bounds(LANE)
        bounds = {k: tuple(x) for k, x in bounds.items()}
        bounds["psi_child"] = tuple(v2.CONFIG["psi_bounds"])
        v2.inputs.LANES[LANE].update(seed=dict(point), bounds=bounds,
                                     free_coordinates=list(point))
        P, grid = v2.inputs.proposal(LANE)
        P, _ = v2.inputs.entry(P, grid, "nonnegative_mean")
        assert P.native_purchase_income is True and P.native_due_stayer_credit is True
        assert P.N_target == 1.0 and np.all(np.asarray(P.phi) == 0.8)
        initializer, ge = make_ge(args.out)
        native = v2.native
        native.install_reporter_on_authored(initializer.authored)
        base_ctx = initializer.authored.context_from_bundle(SimpleNamespace(
            bundle=ROOT / "output/model/publication_refactor_20260929/local_export_v1/inputs",
            reference_root=ROOT, out=args.out))
        initializer.authored.authenticate_frozen(base_ctx)
        fixed_rows = {row["parameter"]: row for row in strict["parameters"]}
        assert len(fixed_rows) == 31
        for row in base_ctx["manifest"]["full_parameter_table"]:
            if row["parameter"] in point:
                fixed = fixed_rows[row["parameter"]]
                row["lower"], row["upper"] = fixed["lower"], fixed["upper"]
                row["near_bound"] = fixed["near_bound"]
        cases = []
        for name, phi in (("strict80_replay", 0.8), ("strict90", 0.9),
                          ("strict90_repeat", 0.9)):
            case_out = args.out / name
            case_out.mkdir()
            P_case = v2.inputs.bind(P, point, bounds, "floor")
            P_case.phi = np.full_like(np.asarray(P_case.phi), phi, dtype=float)
            P_case.H0 = np.full_like(np.asarray(P_case.H0), H0, dtype=float)
            credit.bind_engine_credit(P_case, "corrected", 0.0)
            ctx = dict(base_ctx)
            expected = native.expected_parameters(point, (120, 9), "floor")
            expected.update(H0=H0, financed_share=phi)
            ctx.update(P=P_case, b_grid=grid, out=case_out,
                       selected_d_bar=0.0, reference_psi=float(P_case.psi_child),
                       expected_parameters=expected, expected_dimensions={
                           "wealth_grid_nodes": len(grid), "income_states": len(P_case.z_grid)},
                       free_coordinates=[], fixed_coordinates=list(point),
                       fixed_price=PRICE, deadline_epoch=start + LIMIT,
                       price_start=PRICE)
            actual = ctx["fp"].actual_parameters(ctx["prepared"], P_case, grid)
            ge.validate_parameter_estimates(ctx, ctx["manifest"]["full_parameter_table"], actual)
            budget = initializer.ArmBudget(case_out, start + LIMIT)
            budget.max_lifecycle = 1
            result = ge.run_phase_b(ctx, dict(selected_d_bar=0.0), budget)
            assert result["status"] == "passed_fixed_price_diagnostic" and budget.used_lifecycle == 1
            report = case_out / "phase_b_ge/selected_root"
            target, params = read_table(report / "target_fit.csv"), read_table(report / "parameters.csv")
            assert len(target) == 14 and len(params) == 31
            assert native.target_identity(target) == v2.CONFIG["base_target_contract"]
            assert len(list((report / "standard_diagnostics").glob("*.png"))) == 17
            assert all(abs(float(r["estimate"]) - point[r["parameter"]]) < 1e-11
                       for r in params if r["parameter"] in point)
            if name == "strict80_replay":
                for left, right in zip(target, strict["target_fit"]):
                    assert left["moment"] == right["moment"]
                    for key in ("target", "model", "gap", "weight", "loss_contribution"):
                        if left[key] or right[key]:
                            assert abs(float(left[key]) - float(right[key])) <= 1e-10, (key, left["moment"])
            case = dict(label=name, phi=phi, status=result["status"], price=PRICE, H0=H0,
                        target_fit=target, parameters=params,
                        loss=sum(float(r["loss_contribution"] or 0) for r in target),
                        renewal_residual=result["selected"]["renewal_residual"],
                        housing_residual=result["selected"]["absolute_housing_residual"],
                        report=str(report), lifecycle_solves=budget.used_lifecycle)
            cases.append(case)
            v2.write(args.out / "latest_completed.json", case)
        repeat = native.compare_repeated(Path(cases[1]["report"]), Path(cases[2]["report"]))
        v2.write(args.out / "completed.json", dict(status="fixed_price_diagnostic_passed",
            cases=cases, repeat=repeat, total_lifecycle_solves=3,
            fixed_price=PRICE, fixed_H0=H0, fixed_population=1.0,
            target_fingerprint=manifest["target_fingerprint"],
            weight_fingerprint=manifest["weight_fingerprint"],
            elapsed_seconds=time.time() - start))
    except BaseException as exc:
        v2.write(args.out / "failure.json", dict(type=type(exc).__name__, message=str(exc),
                 completed_cases=len(list(args.out.glob("*/phase_b_ge/selected_root/closure.json")))))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
