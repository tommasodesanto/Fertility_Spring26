"""Bounded Phase A controller. Phase B is a separate reviewed entry point."""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from phase_a import entrant_budget_bound, run_phase_a

BUNDLE_SHA = "427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7"
AUTH_SHA = "96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44"


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temp.replace(path)


class Budget:
    def __init__(self, out, deadline_epoch, *, smoke=False):
        self.out = Path(out)
        self.deadline_epoch = float(deadline_epoch)
        self.smoke = smoke
        self.used_lifecycle = 0
        self.phase_a_lifecycle = 0
        self.max_lifecycle = 10
        self.stage_deadline_seconds = 300
        self.phase_a_started_epoch = time.time()
        self.progress("initialized")

    @property
    def remaining_lifecycle(self):
        return self.max_lifecycle - self.used_lifecycle

    def claim_lifecycle(self, label):
        if self.smoke:
            raise RuntimeError("Smoke must not claim a lifecycle evaluation")
        if self.used_lifecycle >= self.max_lifecycle:
            raise RuntimeError("Total 10-lifecycle cap reached")
        if label.startswith("phase_a"):
            if self.phase_a_lifecycle >= 3 or time.time() - self.phase_a_started_epoch >= 900:
                raise RuntimeError("Phase A 3-lifecycle or 900-second cap reached")
            self.phase_a_lifecycle += 1
        if time.time() + 1 >= self.deadline_epoch:
            raise RuntimeError("2400-second global deadline reached")
        self.used_lifecycle += 1
        self.progress("lifecycle_claimed", label=label)

    def progress(self, status, **extra):
        receipt = dict(status=status, time_epoch=time.time(), deadline_epoch=self.deadline_epoch,
                       lifecycle_used=self.used_lifecycle, lifecycle_remaining=self.remaining_lifecycle,
                       phase_a_lifecycle=self.phase_a_lifecycle, **extra)
        write_json(self.out / "latest.json", receipt)
        if status == "completed_case":
            write_json(self.out / "latest_completed_case.json", receipt)


def context_from_bundle(args):
    from small_credit_lab import inputs
    from small_credit_lab.engine.utils import make_grid
    loaded = inputs.load_inputs(args.bundle, args.reference_root, BUNDLE_SHA)
    P = copy.deepcopy(loaded.parameters)
    required = ("native_due_stayer_credit", "native_exact_inherited_distribution",
                "native_fixed_reference_entry", "native_explicit_transaction_grid",
                "native_purchase_income", "native_exact_allocation_output")
    if any(getattr(P, flag, None) is not True for flag in required):
        raise RuntimeError("Required frozen entrant/accounting flag absent")
    np.testing.assert_array_equal(make_grid(P), loaded.b_grid)
    q_ref = float(np.asarray(loaded.reference_price).reshape(-1)[0])
    return dict(loaded=loaded, P=P, b_grid=loaded.b_grid.copy(), q_ref=q_ref,
                reference_root=args.reference_root, bundle=args.bundle, bundle_sha=BUNDLE_SHA,
                source_root=Path(__file__).parent / "source", out=args.out, write_json=write_json)


def authenticate_frozen(context):
    path = (Path(context["reference_root"]) /
            "output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != AUTH_SHA:
        raise RuntimeError("Frozen observer authenticator differs from pin")
    spec = importlib.util.spec_from_file_location("small_credit_frozen_observer_v1", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Frozen observer cannot load")
    fp = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = fp
    spec.loader.exec_module(fp)
    auth_dir = Path(context["out"]) / "runtime_auth"
    auth_dir.mkdir()
    manifest, _unused, objective, runtime, prepared, reference = fp.authenticate(auth_dir)
    if manifest["checkpoint"]["sha256"] != "b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d":
        raise RuntimeError("Frozen observer checkpoint differs")
    if not np.array_equal(context["b_grid"], np.asarray(reference["b_grid"])):
        raise RuntimeError("Bundle and frozen observer wealth grids differ")
    if not np.array_equal(np.array([context["q_ref"]]), np.asarray(reference["solution"].p_eq).reshape(-1)):
        raise RuntimeError("Bundle and frozen observer prices differ")
    context.update(fp=fp, prepared=prepared, manifest=manifest, objective=objective,
                   runtime=runtime, reference=reference)
    write_json(Path(context["out"]) / "frozen_observer_identity.json",
               dict(authenticator_sha256=digest, checkpoint_sha256=manifest["checkpoint"]["sha256"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("smoke", "full"))
    ap.add_argument("--reference-root", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--deadline-epoch", type=float, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit("Refusing to overwrite existing output")
    args.out.mkdir(parents=True)
    budget = Budget(args.out, args.deadline_epoch, smoke=args.mode == "smoke")
    try:
        context = context_from_bundle(args)
        context["deadline_epoch"] = budget.deadline_epoch
        write_json(args.out / "input_identity.json", context["loaded"].identity)
        authenticate_frozen(context)
        if args.mode == "smoke":
            # Traverse the exact Phase A controller with authenticated inputs,
            # replacing only the lifecycle call with a zero-solve sentinel.
            import single_price
            original = single_price.solve_fixed_price
            calls = []
            def mock_solve(context, d_bar, q, budget, label, stage_dir):
                calls.append(dict(d_bar=d_bar, q=q, label=label))
                if len(calls) != 1 or budget.used_lifecycle:
                    raise AssertionError("Smoke lifecycle loop violated")
                return dict(stage_dir=str(stage_dir), summary={"mock": True})
            single_price.solve_fixed_price = mock_solve
            try:
                result = run_phase_a(context, budget)
            finally:
                single_price.solve_fixed_price = original
            bound = result["entry_budget_bound"]
            assert bound["selected_d_bar"] > bound["unrounded_requirement"] >= 0
            assert budget.used_lifecycle == 0 and len(calls) == 1
            from phase_b_ge import smoke_phase_b
            phase_b_mock = smoke_phase_b(context, budget)
            assert phase_b_mock["lifecycle_solves"] == 0 and budget.used_lifecycle == 0
            receipt = dict(status="passed", lifecycle_solves=0, mock_calls=calls,
                           selected_d_bar=result["selected_d_bar"], phase_b_mock=phase_b_mock,
                           input_identity=context["loaded"].identity)
        else:
            result = run_phase_a(context, budget)
            # Phase B is a separate module with its own economic gates. It
            # reuses the selected q_ref live solution and the same 10-call/
            # 2400-second budget; no default lab GE closure is invoked.
            from phase_b_ge import run_phase_b
            phase_b = run_phase_b(context, result, budget)
            receipt = dict(status="full_passed" if phase_b.get("status") == "passed" else "uncomputed_bounded_budget",
                           lifecycle_solves=budget.used_lifecycle,
                           selected_d_bar=result["selected_d_bar"], summary=result["summary"],
                           phase_b=phase_b,
                           input_identity=context["loaded"].identity,
                           note="Revised positive-credit experiment; frozen reference unchanged")
        write_json(args.out / "completed.json", receipt)
        budget.progress(receipt["status"])
        print(json.dumps(receipt, indent=2, default=str))
    except BaseException as exc:
        failure = dict(status="failed", type=type(exc).__name__, message=str(exc),
                       lifecycle_used=budget.used_lifecycle, no_auto_retry=True)
        write_json(args.out / "failure.json", failure)
        budget.progress("failed", failure=failure)
        raise


if __name__ == "__main__":
    main()
