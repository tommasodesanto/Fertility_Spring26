#!/usr/bin/env python3
"""Three bounded fixed-q0 lifecycle cells for renter-only finite-grid credit."""
import argparse
import copy
import csv
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"

HERE = Path(__file__).resolve().parent
WINNER = HERE.parents[1]
ROOT = HERE.parents[7]
sys.path.insert(0, str(WINNER))
import fixed_price_responses as d
from support_audit import audit_renter_support

SOURCE = ROOT / "code/model/refactor_lab/engine/household.py"
OVERRIDE = HERE / "household_override.py"
EXPECTED_SOURCE = "2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082"
COMPARATOR = WINNER / "purchase_ltv_v1/local_run/retry5/results/baseline_80_80"
CASES = (("baseline_d0", False), ("renter_natural", True), ("renter_natural_repeat", True))


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_override():
    spec = importlib.util.spec_from_file_location("refactor_lab.engine.renter_only_household_override", OVERRIDE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def alarm(*_):
    raise TimeoutError("120-second case cap or ten-minute total deadline")


def isolated_cell(original, support):
    source = inspect.getsource(original)
    save = 'if regime == "reference" and float(factor) == 1.0:'
    replay = 'if regime=="reference" and float(factor)==1.0:'
    closure = '    closure = {"status":'
    require(source.count(save) == 1 and source.count(replay) == 1 and source.count(closure) == 1,
            "Fixed-price driver hook identity changed")
    source = source.replace(save, save[:-1] + ' and out.name == "baseline_d0":')
    source = source.replace(replay, replay[:-1] + ' and out.name=="baseline_d0":')
    source = source.replace(closure,
        '    if bool(getattr(P, "renter_only_natural_credit", False)):\n'
        '        support = audit_renter_support(P, ev, sol, grid, support_trace)\n'
        + closure)
    namespace = dict(d.__dict__, audit_renter_support=support, support_trace=[])
    exec(compile(source, str(HERE / "run_renter_only.py") + "::isolated_cell", "exec"), namespace)
    return namespace


def compare_csv(observed, reference, fields, numeric):
    with open(observed, newline="") as f:
        a = list(csv.DictReader(f))
    with open(reference, newline="") as f:
        b = list(csv.DictReader(f))
    require(len(a) == len(b), "Reference row count changed: " + observed.name)
    for i, (left, right) in enumerate(zip(a, b)):
        for field in fields:
            require(left[field] == right[field], f"Reference identity mismatch {field} row {i}")
        for field in numeric:
            lv, rv = left[field], right[field]
            if lv or rv:
                require(bool(lv) and bool(rv) and abs(float(lv) - float(rv)) <= 2e-12,
                        f"Reference numeric mismatch {field} row {i}")


def check_native_contract(auth):
    P = auth["P"]
    from refactor_lab.engine.parameters import parent_age_maturation_active
    require(len(auth["params_rows"]) == 31 and len(auth["grid"]) == 120 and int(P.Nz) == 9,
            "Frozen winner dimensions changed")
    require(abs(float(d.BINDING["candidate_price"]) - 0.719168368828958) < 1e-15,
            "Frozen winner price changed")
    require(not bool(getattr(P, "native_solvency_credit", False))
            and bool(getattr(P, "native_due_stayer_credit", False))
            and bool(getattr(P, "native_purchase_income", False))
            and getattr(P, "unsecured_credit_limit", None) == 0
            and not bool(getattr(P, "use_pti_constraint", False))
            and not bool(getattr(P, "joint_nested_choice", False))
            and not parent_age_maturation_active(P) and int(P.I) == 1
            and all(float(x) == .8 for x in P.phi), "Native credit or maturation contract changed")


def main(out, total_deadline):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    require(time.time() < total_deadline <= time.time() + 600,
            "Ten-minute total execution deadline missing or expired")
    require(sha(SOURCE) == EXPECTED_SOURCE, "Pinned native household source changed")
    require(COMPARATOR.joinpath("target_fit.csv").is_file(), "Source-consistent local baseline unavailable")
    d.write(out / "source_receipt.json", dict(original_household_sha256=sha(SOURCE),
        override_sha256=sha(OVERRIDE), support_audit_sha256=sha(HERE / "support_audit.py"),
        fixed_price_driver_sha256=sha(d.__file__), driver_sha256=sha(__file__),
        comparator=str(COMPARATOR), absolute_deadline_epoch=total_deadline,
        candidate="verified winner chain7/0173_nm", fixed_price=0.719168368828958))
    cell_namespace = isolated_cell(d.make_price_cell, audit_renter_support)
    price_cell = cell_namespace["make_price_cell"]
    auth = d.authenticate_candidate(out / "runtime_preparation")
    override = load_override()
    P = auth["P"]
    check_native_contract(auth)
    import refactor_lab.engine.equilibrium as eq
    original_bellman = eq.solve_bellman_full_markov_income
    completed = []
    d.write(out / "latest_completed.json", dict(completed=completed, lifecycle_claimed=0))
    try:
        for label, treatment in CASES:
            require(time.time() < total_deadline, "Total deadline reached before " + label)
            case_out = out / label
            case_out.mkdir()
            if treatment:
                P.unsecured_credit_limit = None
                P.renter_only_natural_credit = True
            else:
                P.unsecured_credit_limit = 0.0
                if hasattr(P, "renter_only_natural_credit"):
                    del P.renter_only_natural_credit
            eq.solve_bellman_full_markov_income = override.solve_bellman_full_markov_income
            require(eq.solve_bellman_full_markov_income is override.solve_bellman_full_markov_income,
                    "Isolated Bellman binding failed")
            cell_namespace["support_trace"] = override.RENTER_SUPPORT_TRACE
            case_deadline = min(total_deadline, time.time() + 120)
            d.write(case_out / "attempt.json", dict(label=label, treatment=treatment,
                fixed_price=0.719168368828958, lifecycle_claimed=len(completed) + 1,
                case_deadline_epoch=case_deadline, total_deadline_epoch=total_deadline))
            old = signal.signal(signal.SIGALRM, alarm)
            signal.setitimer(signal.ITIMER_REAL, max(.001, case_deadline - time.time()))
            try:
                closure = price_cell(auth, "reference", 1., case_out, case_deadline)
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
                signal.signal(signal.SIGALRM, old)
            require(time.time() < total_deadline, "Total deadline reached after " + label)
            if not treatment:
                compare_csv(case_out / "target_fit.csv", COMPARATOR / "target_fit.csv",
                    ("moment", "role", "target", "weight"), ("model", "gap", "loss_contribution"))
                compare_csv(case_out / "parameters.csv", COMPARATOR / "parameters.csv",
                    ("parameter", "lower", "upper", "near_bound"), ("estimate",))
            if label == "renter_natural_repeat":
                first = out / "renter_natural"
                for name in ("target_fit.csv", "parameters.csv"):
                    require((case_out / name).read_bytes() == (first / name).read_bytes(),
                            "Selected repeat table mismatch: " + name)
                a = {p.name: sha(p) for p in (first / "standard_diagnostics").glob("*.png")}
                b = {p.name: sha(p) for p in (case_out / "standard_diagnostics").glob("*.png")}
                require(len(a) == 17 and a == b, "Selected repeat standard plots differ")
            closure["renter_only_credit_experiment"] = dict(
                fixed_price=0.719168368828958, renter_floor="finite-grid natural support" if treatment else "zero",
                buyer_origination_financed_share=.8, owner_stayer_DUE_financed_share=.8,
                same_control_PRE=True, not_GE=True, not_calibration=True, production_adoption=False,
                natural_support_certified=False,
                economic_changes=["Replace only renter fixed unsecured floor with finite-grid natural continuation support"] if treatment else [])
            closure["economic_changes"] = closure["renter_only_credit_experiment"]["economic_changes"]
            d.write(case_out / "closure.json", closure)
            require(auth.get("lifecycle_solves", 0) == len(completed) + 1 <= 3,
                    "Three-lifecycle cap or count violated")
            completed.append(dict(label=label, status="completed", target_fit_sha256=sha(case_out / "target_fit.csv"),
                parameters_sha256=sha(case_out / "parameters.csv"),
                standard_plot_sha256={p.name: sha(p) for p in (case_out / "standard_diagnostics").glob("*.png")}))
            d.write(out / "latest_completed.json", dict(completed=completed, lifecycle_claimed=len(completed)))
        d.write(out / "completed.json", dict(status="three_cells_repeat_verified_support_limited",
            completed=completed, lifecycle_claimed=len(completed), natural_support_certified=False,
            production_adoption=False))
    except BaseException as exc:
        d.write(out / "failure.json", dict(status="failed_no_retry", completed=completed,
            lifecycle_claimed=min(3, len(completed) + 1), error=str(exc), traceback=traceback.format_exc()))
        raise
    finally:
        eq.solve_bellman_full_markov_income = original_bellman


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--deadline-epoch", required=True, type=float)
    args = parser.parse_args()
    main(args.out, args.deadline_epoch)
