"""Exploratory adapter: native root/observers without plots or price repeat.

Generated, reviewable source variants preserve the native numerical statements.
No native module or acceptance gate is monkeypatched. Final selected verification
must use the original ``runner.native_evaluator`` with all repeats and plots.
"""
from __future__ import annotations

import difflib
import hashlib
import importlib.util
import inspect
import json
import sys
from pathlib import Path
from types import ModuleType

HERE = Path(__file__).resolve().parent
NATIVE = HERE.parent / "utility_calibration_round1_v1"


def _replace_once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError("Native exploratory source boundary changed: " + old[:80])
    return source.replace(old, new, 1)


def variants(observer_source, root_source, evaluator_source):
    """Return literal-source variants; reject ambiguous native source boundaries."""
    plotting_start = "        # Native plotting expects supply and demand in the same units."
    plotting_end = "        result[\"target_fit_rows\"] = len(fits)"
    if observer_source.count(plotting_start) != 1 or observer_source.count(plotting_end) != 1:
        raise RuntimeError("Native plotting boundary changed")
    start, end = observer_source.index(plotting_start), observer_source.index(plotting_end)
    observer = observer_source[:start] + observer_source[end:]
    observer = _replace_once(observer,
        '        result["standard_plot_supply_units"] = "physical supply divided by endogenous household population"\n', '')
    observer = _replace_once(observer, '        fp.write(out / "closure.json", result)\n',
        '        result["exploration_unverified"] = True\n'
        '        result["selected_price_repeat_performed"] = False\n'
        '        result["standard_plots_deferred_to_final_selection"] = True\n'
        '        fp.write(out / "closure.json", result)\n')
    repeat_start = '    repeat, repeat_live = trial(root["price"], "selected_repeat", repeat=True, reason="fresh exact selected-price repeat")'
    if root_source.count(repeat_start) != 1:
        raise RuntimeError("Native selected-repeat boundary changed")
    root = root_source[:root_source.index(repeat_start)] + '''    result = dict(status="passed", selected_d_bar=d_bar, selected_price=root["price"],
                  selected=certified, points=points, price_search=search,
                  remaining_lifecycle=int(budget.remaining_lifecycle),
                  exploration_unverified=True, selected_price_repeat_performed=False)
    context["fp"].write(Path(context["out"]) / "phase_b_ge" / "selected.json", result)
    return result
'''
    evaluator = _replace_once(evaluator_source, '    import phase_b_pilot as ge\n',
                              '    ge = _get_fast_ge(out)\n')
    repeat_block = '''        # The native selected-price repeat verifies arrays/tables; verify actual PNGs too.
        repeat=directory/'phase_b_ge/selected_repeat_final'
        write(directory/'native_selected_repeat.json',compare_repeated(report,repeat))
'''
    evaluator = _replace_once(evaluator, repeat_block, '')
    evaluator = _replace_once(evaluator, "return dict(status='passed',residual=rr.tolist(),",
        "return dict(status='passed',exploration_unverified=True,experimental_not_adopted=True,residual=rr.tolist(),")
    return observer, root, evaluator


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def make_evaluator(out, lane, P, grid, deadline, price_start=None, *, native_runner=None):
    """Initialize once and return the native evaluator's three-argument closure.

    The caller supplies its authenticated local/cluster runner when applicable.
    The price anchor is fixed for this closure; no solution or price is recycled
    across proposals. Native 300-second cases/32-call caps and reserve stay intact.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if native_runner is None:
        sys.path.insert(0, str(NATIVE))
        native_runner = _load(NATIVE / "runner.py", "floor_fast_native_runner")
    native_runner.verify_sources()
    original_evaluator = inspect.getsource(native_runner.native_evaluator)
    # Loading the original integration establishes the exact single_price path.
    sys.path.insert(0, str(native_runner.BASE))
    _load(native_runner.BASE / "run_comparison.py", "floor_fast_path_initializer")
    sys.path.insert(0, str(NATIVE))
    native_ge = _load(NATIVE / "phase_b_pilot.py", "floor_fast_original_ge")
    originals = (inspect.getsource(native_ge.observe_price),
                 inspect.getsource(native_ge.run_phase_b), original_evaluator)
    sources = variants(*originals)
    receipt = dict(experimental_not_adopted=True, exploration_unverified=True,
        native_root_tolerances_unchanged=True, source_sha256={},
        removed_work=["17 diagnostic plots", "fresh selected-price repeat"],
        runtime_initialized_once=True, cross_candidate_solution_reuse=False,
        cross_candidate_price_anchor_reuse=False)
    review = out / "fast_source_review"
    review.mkdir(exist_ok=True)
    names = ("observe_price", "run_phase_b", "native_evaluator")
    for name, original, variant in zip(names, originals, sources):
        (review / (name + ".py")).write_text(variant)
        (review / (name + ".diff")).write_text("".join(difflib.unified_diff(
            original.splitlines(True), variant.splitlines(True),
            fromfile="native/" + name, tofile="exploratory/" + name)))
        receipt["source_sha256"][name] = dict(
            native=hashlib.sha256(original.encode()).hexdigest(),
            exploratory=hashlib.sha256(variant.encode()).hexdigest())
    (review / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    # Distinct namespace: the original GE module and all its gates are intact.
    fast_ge = ModuleType("floor_fast_ge")
    fast_ge.__dict__.update(native_ge.__dict__)
    exec(compile(sources[0], str(review / "observe_price.py"), "exec"), fast_ge.__dict__)
    exec(compile(sources[1], str(review / "run_phase_b.py"), "exec"), fast_ge.__dict__)
    # A copied function retains its original __globals__; recompile this exact
    # unchanged deadline wrapper so its observe_price lookup uses this namespace.
    wrapper_source = inspect.getsource(native_ge._observe_with_deadline)
    (review / "_observe_with_deadline.py").write_text(wrapper_source)
    exec(compile(wrapper_source, str(review / "_observe_with_deadline.py"), "exec"),
         fast_ge.__dict__)
    namespace = dict(native_runner.__dict__)
    namespace["_get_fast_ge"] = lambda _out: fast_ge
    exec(compile(sources[2], str(review / "native_evaluator.py"), "exec"), namespace)
    return namespace["native_evaluator"](out, lane, P, grid, deadline, price_start)


def compare_saved_baseline(result, verified_report, *, atol=1e-10):
    """Authenticate numerical equivalence before enabling exploratory search."""
    import csv
    import math
    def rows(path):
        with Path(path).open(newline="") as handle:
            return list(csv.DictReader(handle))
    fast, full = Path(result["report"]), Path(verified_report)
    checks = {}
    for name, size in (("target_fit.csv", 14), ("parameters.csv", 31)):
        left, right = rows(fast / name), rows(full / name)
        if len(left) != size or len(right) != size:
            raise RuntimeError("Baseline complete row count differs: " + name)
        key = "moment" if name == "target_fit.csv" else "parameter"
        if [r[key] for r in left] != [r[key] for r in right]:
            raise RuntimeError("Baseline row identity differs: " + name)
        fields = ("target", "model", "gap", "weight", "loss_contribution") if size == 14 else ("estimate",)
        maximum = 0.
        for a, b in zip(left, right):
            for field in fields:
                if a[field] == b[field] == "":
                    continue
                x, y = float(a[field]), float(b[field])
                if not math.isfinite(x) or not math.isfinite(y) or abs(x-y) > atol:
                    raise RuntimeError("Fast/full baseline differs: " + name + ":" + a[key] + ":" + field)
                maximum = max(maximum, abs(x-y))
        checks[name] = dict(rows=size,maximum_absolute_error=maximum)
    return dict(status="matched_saved_full_baseline",absolute_tolerance=atol,checks=checks)
