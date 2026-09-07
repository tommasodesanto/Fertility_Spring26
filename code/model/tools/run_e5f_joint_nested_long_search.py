#!/usr/bin/env python3
"""Bounded receipt-driven long search for the experimental joint-nested E5F model.

This controller never manufactures an objective value: an incumbent is updated only
from a complete ``case_receipt.json`` written by the full historical adapter.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import signal
import subprocess
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

import run_e5f_joint_overnight_case as adapter
import build_e5f_bounded_refinement_plan as planner
planner.adapter = adapter

N, RNG_SEED = 11, 20260906
ROOT = Path(__file__).resolve().parents[3]
FINAL_VERIFICATION_HISTORIES = 24  # 22 Jacobian probes plus two exact repeats.
RUN_PROFILES = {
    "v1": {"max_workers": 12, "case_timeout_seconds": 3600, "max_histories": 360,
           "max_search_seconds": 32400, "max_total_seconds": 43200, "population_size": 32,
           "max_generations": 8, "polish_rounds": 2, "smoke_histories": 4},
    "wide32": {"max_workers": 32, "case_timeout_seconds": 3600, "max_histories": 640,
               "max_search_seconds": 32400, "max_total_seconds": 43200, "population_size": 64,
               "max_generations": 8, "polish_rounds": 2, "smoke_histories": 4,
               "final_reserve_seconds": 16200},
    "parallel32": {"max_workers": 32, "case_timeout_seconds": 3600, "max_histories": 640,
                   "max_search_seconds": 32400, "max_total_seconds": 43200, "population_size": 32,
                   "max_generations": 8, "polish_rounds": 2, "smoke_histories": 4,
                   "final_reserve_seconds": 12600},
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def transform(u, spec):
    lo, hi, kind = spec["lower"], spec["upper"], spec["transform"]
    if kind == "log": return lo * (hi / lo) ** u
    if kind == "discount": return lo + (hi - lo) * (1 - (1 - u) ** 2)
    if kind == "softzero": return lo + (hi - lo) * u ** 2
    if kind == "asinh": return math.sinh((1-u)*math.asinh(lo)+u*math.asinh(hi))
    raise ValueError(f"unknown transform {kind}")


def unit_to_payload(seed, domain, unit, label):
    if len(domain) != N or len(unit) != N or any(not math.isfinite(x) or not 0 <= x <= 1 for x in unit):
        raise ValueError("Invalid eleven-dimensional proposal")
    out = {"old_psi_child": seed["old_psi_child"], "best_candidate": {
        "theta": copy.deepcopy(seed["best_candidate"]["theta"]),
        "old_psi_child": seed["old_psi_child"], "new_psi_child": seed["best_candidate"]["new_psi_child"]}}
    for spec, u in zip(domain, unit):
        value, name = transform(float(u), spec), spec["name"]
        if name == "beta_annual": out["best_candidate"]["theta"]["beta"] = value ** 4
        elif name == "psi_child_change_2023": out["best_candidate"]["new_psi_child"] = out["old_psi_child"] + value
        else: out["best_candidate"]["theta"][name] = value
    out["proposal"] = {"label": label, "unit_vector": list(unit), "status": "unevaluated full-objective proposal"}
    return out


def classify_failure(error, error_type=""):
    """Only declared feasibility/numerical failures may be skipped."""
    text = str(error)
    if error_type == "InfeasibleThetaError" or "InfeasibleThetaError" in text:
        return "infeasible_theta"
    if "Housing market did not clear" in text or "market gate failed" in text:
        return "market_nonconvergence"
    if text.startswith(("Old-steady-state fertility normalization is not bracketed", "Old-steady-state fertility normalization missed tolerance")):
        return "fertility_normalization"
    if text.startswith(("The dated first-birth branch has zero treated mass", "The dated first-birth branch has zero destination mass")):
        return "undefined_first_birth_support"
    return None


def bounded_perturbation(center, radius, rng):
    return [min(1., max(0., x+rng.uniform(-radius, radius))) for x in center]


def best_completed(rows):
    """Select only receipt-backed rows; rejection records intentionally lack loss."""
    return min(rows, key=lambda row: row["loss"]) if rows else None


def validate_policy_receipt(receipt, contract, *, smoke, selected_hashes):
    expected_cases = {"baseline", "supply-plus-20", "dependent-child-ltv95", "property-tax-2pct-no-rebate"}
    if (receipt.get("status") != "complete" or receipt.get("failures")
            or receipt.get("smoke") is not smoke or set(receipt.get("cases", {})) != expected_cases):
        raise RuntimeError("Incomplete policy-loop receipt")
    if receipt.get("policy_workers", 1) != contract.get("policy_workers", 1):
        raise RuntimeError("Policy receipt has different contracted concurrency")
    if (receipt.get("scientific_bundle") != contract["code_bundle_sha256"]
            or receipt.get("target_fingerprint") != contract["target_fingerprint"]
            or receipt.get("selected_summary_sha256") not in selected_hashes):
        raise RuntimeError("Policy receipt source or selected history changed")
    for case in receipt["cases"].values():
        if (case.get("status") != "complete" or case.get("dates") != (2 if smoke else 11)
                or case.get("source_summary_sha256") != receipt["selected_summary_sha256"]):
            raise RuntimeError("Incomplete policy branch or mixed selected history")
        gates = case["gates"]
        for key, limit in (("maximum_market_residual", 2e-4), ("maximum_mass_residual", 2e-10)):
            if not math.isfinite(gates[key]) or gates[key] > limit:
                raise RuntimeError("Policy receipt violates unchanged market or mass gates")



def verify_parallel_policy_budget(contract):
    """Require measured complete four-process smoke before using its shorter reserve."""
    if contract.get("policy_workers") != 4:
        raise RuntimeError("Parallel search requires four contracted policy workers")
    evidence = contract.get("parallel_policy_timing", {})
    path = Path(evidence.get("receipt", ""))
    adapter.verify(path, evidence.get("receipt_sha256", ""))
    receipt = adapter.read_json(path)
    validate_policy_receipt(receipt, contract, smoke=True,
                            selected_hashes={receipt.get("selected_summary_sha256")})
    elapsed = receipt.get("elapsed_seconds")
    if not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed <= 0:
        raise RuntimeError("Parallel policy smoke has no positive finite elapsed time")
    projected = elapsed * 44 / 8
    # Two one-hour historical waves, a measured 44-date policy forecast, and
    # twenty minutes of buffer must fit the fixed three-and-a-half-hour reserve.
    if projected > 4200 or evidence.get("projected_full_policy_seconds") != projected:
        raise RuntimeError("Measured policies do not fit the declared final reserve")
    imported = contract.get("imported_smoke") or {}
    proof_path = Path(imported.get("root", "")) / "policy_loop_verification.json"
    adapter.verify(proof_path, imported.get("policy_loop_verification_sha256", ""))
    proof = adapter.read_json(proof_path)
    if (Path(proof.get("receipt", "")).resolve() != path.resolve()
            or proof.get("sha256") != evidence["receipt_sha256"]):
        raise RuntimeError("Policy timing must come from the imported complete smoke")
    return projected


def initial_population(center, domain, rng, profile="v1"):
    if profile not in RUN_PROFILES:
        raise ValueError(f"Unsupported run profile: {profile}")
    jump = next(i for i, d in enumerate(domain) if d["name"] == "hbar_first_child_jump")
    kappa = next(i for i, d in enumerate(domain) if d["name"] == "tenure_choice_kappa")
    lam = next(i for i, d in enumerate(domain) if d["name"] == "joint_nest_lambda")
    # Convert diagnostic physical scales by a small monotone grid search.
    inverse = lambda value, s: math.log(value/s["lower"])/math.log(s["upper"]/s["lower"])
    result = [list(center)]
    if profile == "wide32":
        kappa_grid = (.01, .03, .1, .3, 1., 2., 4., 8.)
        lambda_grid = (.02, .05, .2, .5, .8, 1.)
        omitted = (2., .8)  # The exact anchor replaces this nearest grid row.
        for kappa_value in kappa_grid:
            for lambda_value in lambda_grid:
                if (kappa_value, lambda_value) == omitted:
                    continue
                u = bounded_perturbation(center, .08, rng)
                u[kappa] = inverse(kappa_value, domain[kappa])
                u[lam] = inverse(lambda_value, domain[lam])
                result.append(u)
        # Retain the original grid and pair its small-inner-scale rows with
        # smaller preference declines. This changes starting points only:
        # all eleven coordinates remain unrestricted in subsequent DE.
        psi = next(i for i, d in enumerate(domain) if d["name"] == "psi_child_change_2023")
        psi_spec = domain[psi]
        if psi_spec["transform"] != "asinh" or not psi_spec["lower"] < 0 < psi_spec["upper"]:
            raise ValueError("Preference-change proposals require the declared signed asinh domain")
        anchor_inner = transform(center[kappa], domain[kappa]) * transform(center[lam], domain[lam])
        for row in result[1:]:
            inner = transform(row[kappa], domain[kappa]) * transform(row[lam], domain[lam])
            if not any(math.isclose(transform(row[lam], domain[lam]), value) for value in (.02, .2)):
                continue
            u = list(row)
            delta = transform(row[psi], psi_spec) * min(.5, inner / anchor_inner)
            u[psi] = (math.asinh(delta) - math.asinh(psi_spec["lower"])) / (math.asinh(psi_spec["upper"]) - math.asinh(psi_spec["lower"]))
            result.append(u)
        if len(result) != RUN_PROFILES[profile]["population_size"]:
            raise RuntimeError("wide32 initial population does not match its profile")
        return result
    for i in range(RUN_PROFILES[profile]["population_size"] - 1):
        u = bounded_perturbation(center, .05, rng)
        u[kappa] = inverse((.5, 1., 2., 4.)[i % 4], domain[kappa])
        u[lam] = inverse((.2, .5, .8, 1.)[(i // 4) % 4], domain[lam])
        # Preserve an anchor-near jump as well as meaningful bounded variants.
        u[jump] = min(1., max(0., center[jump] + (-.04, -.015, .015, .04)[i % 4]))
        result.append(u)
    return result


def de_trials(pop, losses, rng):
    best = min(range(len(pop)), key=lambda i: (losses[i], i)); trials = []
    for i, parent in enumerate(pop):
        a, b, c = rng.sample([j for j in range(len(pop)) if j != i], 3); f = rng.uniform(.5, .9)
        mutant = [parent[k]+f*(pop[best][k]-parent[k])+f*(pop[a][k]-pop[b][k]) for k in range(N)]
        mutant = [min(1., max(0., x)) for x in mutant]; forced = rng.randrange(N)
        trials.append([mutant[k] if k == forced or rng.random() < .9 else parent[k] for k in range(N)])
    return trials


def verify_contract(path, expected_sha):
    adapter.verify(path, expected_sha); c = adapter.read_json(path)
    if c.get("schema") != "e5f_joint_nested_long_v1": raise RuntimeError("Unknown long-search contract")
    for key, source in (("controller_sha256", __file__), ("adapter_sha256", adapter.__file__),
                        ("planner_sha256", planner.__file__)):
        if key not in c: raise RuntimeError(f"Missing {key}")
        adapter.verify(source, c[key])
    for name, sha in c.get("reused_helper_sha256", {}).items(): adapter.verify(ROOT / "code/model/tools" / name, sha)
    for key in ("seed_center", "seed_reference"):
        adapter.verify(c[key], c[f"{key}_sha"] if f"{key}_sha" in c else c[f"{key}_sha256"])
    profile = c.get("run_profile", "v1")
    if profile not in RUN_PROFILES:
        raise RuntimeError("Unsupported bounded controller profile")
    expected = {**RUN_PROFILES[profile], "random_seed": RNG_SEED}
    if any(c.get(k) != v for k, v in expected.items()): raise RuntimeError("Changed bounded controller design")
    if c["search_domain"] != adapter.SEARCH_DOMAIN: raise RuntimeError("Changed parameter domain")
    base = c["base_plan"]
    if "base_plan_sha256" not in c:
        raise RuntimeError("Missing immutable base_plan_sha256")
    canonical = json.dumps(base, indent=2, sort_keys=True, allow_nan=False).encode() + b"\n"
    if hashlib.sha256(canonical).hexdigest() != c["base_plan_sha256"]:
        raise RuntimeError("base_plan JSON does not match base_plan_sha256")
    if base.get("schema") != "e5f_joint_nested_overnight_v1": raise RuntimeError("Wrong base plan")
    for key in ("source_sha256", "code_bundle_sha256", "target_fingerprint", "search_domain"):
        if c.get(key) != base.get(key): raise RuntimeError(f"Mixed {key}")
    if c["target_fingerprint"] != adapter.TARGET or c["source_sha256"] != adapter.SOURCE: raise RuntimeError("Adapter/source target mismatch")
    if profile == "parallel32": verify_parallel_policy_budget(c)
    c["run_profile"] = profile
    c["contract_path"] = str(Path(path).resolve())
    return c


class Search:
    def __init__(self, contract, mode):
        self.c, self.mode = contract, mode; self.root = Path(contract["output_root"]) / mode
        if self.root.exists(): raise RuntimeError(f"Refusing existing run directory: {self.root}")
        self.root.mkdir(parents=True); self.started = time.monotonic(); self.wall = time.time()
        self.finish = min(self.wall + contract["max_total_seconds"], contract["absolute_finish_epoch"])
        self.search_finish = min(self.wall + contract["max_search_seconds"],
                                 self.finish - contract.get("final_reserve_seconds", 10800))
        self.seed = adapter.read_json(contract["seed_center"]); self.ledger, self.rejects, self.best = [], [], None
        self.stop_event = threading.Event(); self.lock = threading.Lock()
        self.completed = 0; self.consecutive_timeouts = 0; self.active = {}; self.phase = "initializing"
        self.search_stop_reason = None
        self.state("running")

    def state(self, status, **more):
        write_json(self.root / "search_state.json", {"status": status, "phase": self.phase,
            "elapsed_seconds": time.monotonic()-self.started, "completed_histories": self.completed,
            "max_histories": self.c["max_histories"], "best_loss": self.best["loss"] if self.best else None,
            "search_stop_reason": self.search_stop_reason,
            "absolute_finish_epoch": self.finish, "epoch": time.time(), "production_promoted": False, **more})

    def can_fit(self, n, final=False):
        if not final and getattr(self, "search_stop_reason", None): return False
        if self.completed + n + (0 if final else FINAL_VERIFICATION_HISTORIES) > self.c["max_histories"]: return False
        # A conservative full timeout is required per wave; final steps may use reserved time.
        cutoff = self.finish if final else self.search_finish
        return time.time() + math.ceil(n/self.c["max_workers"]) * self.c["case_timeout_seconds"] <= cutoff

    def new_plan(self, stage, vectors, labels, case_offset=0):
        if len(vectors) != len(labels):
            raise ValueError("Every proposal requires exactly one label")
        dest = self.root / stage; dest.mkdir(parents=True, exist_ok=False); plan = copy.deepcopy(self.c["base_plan"])
        plan.update(stage=stage, cases=[], input_sha256={}, controller_sha256=self.c["controller_sha256"],
                    adapter_sha256=self.c["adapter_sha256"], planner_sha256=self.c["planner_sha256"],
                    launch_deadline_epoch=(self.finish-60 if stage.startswith("final") else min(self.search_finish, self.finish-60)))
        for i, (u, label) in enumerate(zip(vectors, labels), case_offset + 1):
            center = dest / f"center_{i:03d}.json"; adapter.write_json(center, unit_to_payload(self.seed, self.c["search_domain"], u, label))
            plan["cases"].append({"id": i, "label": label, "center": center.name, "center_sha256": digest(center),
                "panel_task_id": 1, "panel_size": 1, "panel_design": "mixed", "panel_seed": RNG_SEED,
                "radius": .05, "output": f"task_{i:03d}"})
        if stage == "final_repeats":
            origin = self.repeat_origin
            source_plan_path = Path(origin["plan"])
            source_plan = adapter.load_plan(source_plan_path, origin["plan_sha256"])
            original = next(c for c in source_plan["cases"] if c["id"] == origin["id"])
            for case in plan["cases"]:
                center = dest/case["center"]
                center.write_bytes((source_plan_path.parent/original["center"]).read_bytes())
                for name in ("panel_task_id", "panel_size", "panel_design", "panel_seed", "radius"):
                    case[name] = original[name]
                case["center_sha256"] = digest(center)
        path = dest / "plan.json"; adapter.write_json(path, plan); return path, digest(path)

    def _child(self, plan, sha, case):
        out = plan.parent / case["output"]; log = plan.parent / f"case_{case['id']:03d}.log"
        if self.stop_event.is_set():
            raise RuntimeError("Controller stopped before candidate launch")
        if time.time() > adapter.read_json(plan)["launch_deadline_epoch"]:
            raise RuntimeError("Candidate launch deadline passed")
        with log.open("w") as stream:
            proc = subprocess.Popen([sys.executable, adapter.__file__, "--plan", str(plan), "--plan-sha256", sha, "--case-id", str(case["id"])], stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            with self.lock:
                self.active[proc.pid] = proc
                if self.stop_event.is_set():
                    os.killpg(proc.pid, signal.SIGTERM)
            try: code = proc.wait(timeout=min(self.c["case_timeout_seconds"], max(1, self.finish-time.time())))
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGTERM)
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: os.killpg(proc.pid, signal.SIGKILL); proc.wait()
                return case, out, "timeout", "case timeout"
            finally:
                with self.lock: self.active.pop(proc.pid, None)
        if code == 0: return case, out, "complete", ""
        failure = adapter.read_json(out / "adapter_failure.json") if (out / "adapter_failure.json").exists() else {"error": f"case exited {code}", "type": "ProcessExit"}
        return case, out, "failed", failure

    def stop_active(self):
        self.stop_event.set()
        with self.lock:
            active = list(self.active.values())
        for proc in active:
            if proc.poll() is None:
                try: os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError: pass
        for proc in active:
            try: proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try: os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                proc.wait(timeout=10)

    def _record_completed(self, plan, sha, case, out):
        receipt = adapter.read_json(out / "case_receipt.json")
        if receipt.get("status") != "complete" or receipt.get("plan_sha256") != sha: raise RuntimeError("Invalid complete receipt")
        for name, value in receipt["artifact_sha256"].items(): adapter.verify(out / name, value)
        summary, _, _ = adapter.validate_result(out, adapter.load_plan(plan, sha), case)
        if receipt.get("case_id") != case["id"]:
            raise RuntimeError("Receipt belongs to another case")
        loss = float(receipt["loss"])
        if loss != float(summary["best_candidate"]["transition_loss"]):
            raise RuntimeError("Receipt objective differs from complete target fit")
        if not math.isfinite(loss): raise RuntimeError("Nonfinite actual calibrated loss")
        row = {"stage": plan.parent.name, "id": case["id"], "label": case["label"], "loss": loss,
               "summary": str((out/"summary.json").resolve()), "plan": str(plan.resolve()), "plan_sha256": sha,
               "unit_vector": summary["panel_design"]["unit_vector"], "elapsed_seconds": receipt["elapsed_seconds"]}
        self.ledger.append(row); self.completed += 1; self.consecutive_timeouts = 0
        if self.best is None or loss < self.best["loss"]: self.best = row

    def _reject(self, case, out, kind, error):
        self.completed += 1
        if kind == "timeout": self.consecutive_timeouts += 1
        else: self.consecutive_timeouts = 0
        unit = adapter.read_json(out.parent / case["center"]).get("proposal", {}).get("unit_vector")
        self.rejects.append({"stage": self.phase, "case_id": case["id"], "label": case["label"], "rejection_type": kind,
            "error": str(error), "unit_vector": unit, "source_sha256": self.c["source_sha256"],
            "target_fingerprint": self.c["target_fingerprint"], "elapsed_seconds": time.monotonic()-self.started})
        if self.consecutive_timeouts >= 3 and not getattr(self, "search_stop_reason", None):
            self.search_stop_reason = {"reason": "three_consecutive_timeouts", "stage": self.phase,
                "action": "stop new search proposals; finish active cases under existing caps; verify the valid incumbent",
                "elapsed_seconds": time.monotonic()-self.started, "consecutive_timeouts": self.consecutive_timeouts}
            write_json(self.root/"search_stop.json", self.search_stop_reason)

    def reports(self):
        if self.ledger:
            write_json(self.root/"latest_completed_case.json", {"latest": self.ledger[-1], "production_promoted": False})
            write_json(self.root/"best_so_far.json", {"best": self.best, "loss": self.best["loss"], "production_promoted": False})
        write_csv(self.root/"all_cases.csv", self.ledger); write_csv(self.root/"rejects_ledger.csv", self.rejects)
        fits=[]; parameters=[]
        for row in self.ledger:
            base=Path(row["summary"]).parent
            for x in adapter.read_csv(base/"target_fit_long.csv"): fits.append({"stage":row["stage"],"case_id":row["id"],**x})
            for x in adapter.read_csv(base/"parameter_table.csv"): parameters.append({"stage":row["stage"],"case_id":row["id"],**x})
        write_csv(self.root/"all_target_fits.csv", fits); write_csv(self.root/"all_parameters.csv", parameters)

    def batch(self, stage, vectors, labels, *, smoke=False):
        if not self.can_fit(len(vectors), final=stage.startswith("final")):
            raise RuntimeError("Budget/deadline cannot fit declared stage")
        self.phase = stage
        # The unchanged adapter accepts at most forty cases in a plan.
        # Preserve global population IDs across bounded plans in one worker queue.
        plans = []
        if len(vectors) > 40:
            for start in range(0, len(vectors), 32):
                plans.append(self.new_plan(f"{stage}_part_{start // 32 + 1:02d}",
                    vectors[start:start+32], labels[start:start+32], case_offset=start))
        else:
            plans.append(self.new_plan(stage, vectors, labels))
        futures = {}; rejected_before = len(self.rejects)
        pool = ThreadPoolExecutor(max_workers=self.c["max_workers"])
        try:
            cases = iter((plan, sha, case) for plan, sha in plans
                         for case in adapter.load_plan(plan, sha)["cases"])
            def submit():
                if getattr(self, "search_stop_reason", None) and not stage.startswith("final"):
                    return
                entry = next(cases, None)
                if entry is not None:
                    plan, sha, case = entry
                    futures[pool.submit(self._child, plan, sha, case)] = (plan, sha)
            for _ in range(min(len(vectors), self.c["max_workers"])): submit()
            while futures:
                done, _ = wait(futures, timeout=60, return_when=FIRST_COMPLETED)
                if not done:
                    self.reports(); self.state("running")
                    if time.time() >= self.finish: raise RuntimeError("Absolute finish deadline reached")
                    continue
                for future in done:
                    plan, sha = futures.pop(future)
                    case, out, status, detail = future.result()
                    if status == "complete": self._record_completed(plan, sha, case, out)
                    else:
                        typ = "timeout" if status == "timeout" else classify_failure(detail.get("error", ""), detail.get("type", ""))
                        if typ is None: raise RuntimeError(f"Fatal candidate failure: {detail}")
                        self._reject(case, out, typ, detail)
                        if smoke: raise RuntimeError(f"Required verification case rejected: {typ}: {detail}")
                self.reports(); self.state("running")
                for _ in done: submit()
            unstarted = [{"plan": str(plan), "plan_sha256": sha, **case} for plan, sha, case in cases]
            if unstarted:
                write_json(self.root/f"{stage}_unstarted.json", {"status": "not_run_after_timeout_stop",
                    "cases": unstarted, "counted_as_completed_or_rejected": False, "search_stop_reason": self.search_stop_reason})
            for plan, sha in plans:
                planner.collect(plan, sha, require_complete=len(self.rejects) == rejected_before)
        except BaseException:
            for future in futures: future.cancel()
            self.stop_active()
            raise
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
        hashes = {sha for _, sha in plans}
        return [r for r in self.ledger if r["plan_sha256"] in hashes]

    def smoke(self):
        u=list(self.seed["panel_design"]["unit_vector"])
        paired=[]
        for sign in (-1,1): paired.append([min(1.,max(0.,x+sign*.00125)) for x in u])
        rows=self.batch("smoke_histories", [u,u,*paired],
            ["anchor_1","anchor_2","all_minus","all_plus"], smoke=True)
        anchors=[r for r in rows if r["label"] in ("anchor_1","anchor_2")]
        probes=[r for r in rows if r["label"] in ("all_minus","all_plus")]
        if len(anchors)!=2 or len(probes)!=2:
            raise RuntimeError("Four completed smoke histories required")
        paths=[Path(r["summary"]).parent for r in sorted(anchors,key=lambda r:r["id"])]
        exact=adapter.compare_reference(paths[1], paths[0]/"summary.json")
        graphs=sorted((paths[0]/"standard_diagnostics").glob("*.png"))
        if len(graphs)!=17: raise RuntimeError("Missing standard diagnostic packet")
        for p in graphs: adapter.verify(paths[1]/"standard_diagnostics"/p.name, digest(p))
        if any(sum(a!=b for a,b in zip(r["unit_vector"],u)) != N for r in probes): raise RuntimeError("Smoke probes did not move all eleven dimensions")
        write_json(self.root/"smoke_verification.json", {"status":"pass", "exact":exact, "exact_standard_pngs":17,
            "anchor_receipts":[str(Path(r["summary"]).parent/"case_receipt.json") for r in anchors],
            "probe_receipts":[str(Path(r["summary"]).parent/"case_receipt.json") for r in probes]})
        if self.c.get("finalizer_driver"):
            driver = self.c["finalizer_driver"]
            adapter.verify(driver, self.c["finalizer_sha256"])
            output = Path(self.c["output_root"])/"policy_loop_smoke"
            with (self.root/"policy_loop_smoke.log").open("w") as log:
                result = subprocess.run([sys.executable, driver, "--selected-summary", str(paths[0]/"summary.json"),
                    "--outdir", str(output), "--contract", str(self.c["contract_path"]), "--smoke"],
                    stdout=log, stderr=subprocess.STDOUT, timeout=3600)
            if result.returncode: raise RuntimeError("Policy-loop smoke failed")
            receipt = adapter.read_json(output/"equilibrium_receipt.json")
            if receipt.get("status") != "complete" or not receipt.get("smoke") or len(receipt.get("cases", {})) != 4:
                raise RuntimeError("Incomplete policy-loop smoke receipt")
            write_json(self.root/"policy_loop_verification.json", {"status": "pass",
                "receipt": str(output/"equilibrium_receipt.json"), "sha256": digest(output/"equilibrium_receipt.json"),
                "driver_sha256": self.c["finalizer_sha256"]})
        self.summary("smoke_passed")

    def require_smoke(self):
        imported = self.c.get("imported_smoke")
        if imported:
            required = ("root", "original_contract", "original_contract_sha256", "smoke_verification_sha256")
            if any(key not in imported for key in required):
                raise RuntimeError("Incomplete imported smoke provenance")
            root = Path(imported["root"])
            original_contract = Path(imported["original_contract"])
            adapter.verify(original_contract, imported["original_contract_sha256"])
            original = adapter.read_json(original_contract)
            for key in ("source_sha256", "code_bundle_sha256", "target_fingerprint", "search_domain",
                        "adapter_sha256", "planner_sha256", "finalizer_sha256", "numerical_gates", "closure"):
                if original.get(key) != self.c.get(key):
                    raise RuntimeError(f"Imported smoke original contract has mixed {key}")
            if original.get("policy_workers", 1) != self.c.get("policy_workers", 1):
                raise RuntimeError("Imported smoke has different contracted policy concurrency")
            proof_path = root/"smoke_verification.json"
            adapter.verify(proof_path, imported["smoke_verification_sha256"])
        else:
            if self.c["run_profile"] != "v1":
                raise RuntimeError("Hash-pinned imported full-loop smoke is required")
            root = Path(self.c["output_root"])/"smoke"
            original_contract = None
            proof_path = root/"smoke_verification.json"
        proof = adapter.read_json(proof_path)
        if proof.get("status") != "pass" or proof.get("exact_standard_pngs") != 17:
            raise RuntimeError("Unchanged full-loop smoke required")
        if len(proof.get("anchor_receipts", [])) != 2 or len(proof.get("probe_receipts", [])) != 2:
            raise RuntimeError("Four completed smoke cases required")
        for value in proof["anchor_receipts"] + proof["probe_receipts"]:
            receipt = adapter.read_json(value); out = Path(value).parent; plan_path = out.parent/"plan.json"
            plan = adapter.load_plan(plan_path, receipt["plan_sha256"])
            for key in ("source_sha256", "code_bundle_sha256", "target_fingerprint", "search_domain"):
                if plan.get(key) != self.c.get(key):
                    raise RuntimeError(f"Smoke plan belongs to different {key}")
            if imported:
                for key in ("controller_sha256", "adapter_sha256", "planner_sha256"):
                    if plan.get(key) != original.get(key):
                        raise RuntimeError(f"Smoke plan differs from original contract {key}")
                for key in ("helper_sha256", "target_set", "first_child_jump_upper"):
                    if plan.get(key) != original["base_plan"].get(key):
                        raise RuntimeError(f"Smoke plan differs from original base plan {key}")
            case = next(c for c in plan["cases"] if c["id"] == receipt["case_id"])
            self._record_completed(plan_path, receipt["plan_sha256"], case, out)
        paths = [Path(p).parent for p in proof["anchor_receipts"]]
        adapter.compare_reference(paths[1], paths[0]/"summary.json")
        graphs = sorted((paths[0]/"standard_diagnostics").glob("*.png"))
        if len(graphs) != 17: raise RuntimeError("Incomplete smoke graphs")
        for graph in graphs: adapter.verify(paths[1]/"standard_diagnostics"/graph.name, digest(graph))
        if self.c.get("finalizer_driver"):
            policy_path = root/"policy_loop_verification.json"
            if imported:
                if "policy_loop_verification_sha256" not in imported:
                    raise RuntimeError("Missing imported policy smoke provenance")
                adapter.verify(policy_path, imported["policy_loop_verification_sha256"])
            policy_proof = adapter.read_json(policy_path)
            adapter.verify(policy_proof["receipt"], policy_proof["sha256"])
            if policy_proof["status"] != "pass" or policy_proof["driver_sha256"] != self.c["finalizer_sha256"]:
                raise RuntimeError("Policy-loop smoke source changed")
            receipt = adapter.read_json(policy_proof["receipt"])
            validate_policy_receipt(receipt, self.c, smoke=True,
                selected_hashes={digest(path/"summary.json") for path in paths})
            selected_path = Path(receipt["selected_summary"])
            if not selected_path.is_absolute():
                if "receipt_working_directory" not in policy_proof:
                    raise RuntimeError("Relative policy source requires its original working directory")
                selected_path = Path(policy_proof["receipt_working_directory"])/selected_path
            if selected_path.resolve() not in {(path/"summary.json").resolve() for path in paths}:
                raise RuntimeError("Policy receipt does not refer to an original smoke anchor")
            adapter.verify(selected_path, receipt["selected_summary_sha256"])
            adapter.verify(Path(policy_proof["receipt"]).parent/"inherited_state_verification.json",
                           receipt["inherited_state_verification_sha256"])
        record = {"path": str(proof_path), "sha256": digest(proof_path), "case_count": self.completed}
        if imported:
            record.update(original_contract=str(original_contract), original_contract_sha256=imported["original_contract_sha256"])
        write_json(self.root/"inherited_smoke.json", record)

    def search(self):
        self.require_smoke(); rng=random.Random(RNG_SEED); center=list(self.seed["panel_design"]["unit_vector"])
        pop=initial_population(center,self.c["search_domain"],rng,self.c["run_profile"])
        initial=self.batch("initial_population",pop,["seed"]+[f"initial_{i}" for i in range(1,len(pop))])
        # Every slot starts from an actually evaluated feasible point. A
        # rejected initial proposal inherits a valid point, never a fake loss.
        available = sorted(initial, key=lambda row: row["loss"]) or [self.best]
        by_id = {row["id"]: row for row in initial}
        seeded = [by_id.get(i + 1, available[i % len(available)]) for i in range(self.c["population_size"])]
        pop = [list(row["unit_vector"]) for row in seeded]
        losses = [row["loss"] for row in seeded]
        for gen in range(1, self.c["max_generations"] + 1):
            if not self.can_fit(self.c["population_size"]): break
            trials=de_trials(pop,losses,rng); results=self.batch(f"de_generation_{gen:02d}",trials,[f"de_{gen}_{i}" for i in range(self.c["population_size"])])
            lookup = {row["id"]: row for row in results}
            for i in range(len(trials)):
                trial = lookup.get(i + 1)
                if trial and trial["loss"] < losses[i]:
                    pop[i], losses[i] = list(trial["unit_vector"]), trial["loss"]
            write_json(self.root/"population.json", {"generation":gen,"unit_vectors":pop,"evaluated_actual_losses":losses})
        for round_ in range(1, self.c["polish_rounds"] + 1):
            if not self.can_fit(2*N): break
            anchor_loss = self.best["loss"]
            base=list(self.best["unit_vector"]); radius=.00125/(2**(round_-1)); vectors=[]; labels=[]
            for j in range(N):
                for sign,name in ((-1,"minus"),(1,"plus")):
                    u=list(base); u[j]=min(1.,max(0.,u[j]+sign*radius)); vectors.append(u); labels.append(f"coordinate_{j}_{name}")
            probes=self.batch(f"polish_{round_}_coordinates",vectors,labels)
            direction=[0.]*N
            for j in range(N):
                options=[r for r in probes if r["label"].startswith(f"coordinate_{j}_")]
                if options:
                    winner=min(options,key=lambda r:r["loss"])
                    if winner["loss"]<anchor_loss: direction[j]=winner["unit_vector"][j]-base[j]
            joint=[]; jl=[]
            for scale in (.5,1.,2.,4.):
                u=[min(1.,max(0.,x+scale*d)) for x,d in zip(base,direction)]
                if u!=base and u not in joint: joint.append(u); jl.append(f"joint_all_{scale:g}")
            if joint and self.can_fit(len(joint)): self.batch(f"polish_{round_}_joint",joint,jl)
        self.final_assessment()

    def final_assessment(self):
        if self.best is None: return self.summary("no_valid_completed_case")
        anchor = copy.deepcopy(self.best); base = list(anchor["unit_vector"])
        vectors, labels = [], []
        for j in range(N):
            for sign, name in ((-1, "minus"), (1, "plus")):
                u = list(base); u[j] = min(1., max(0., u[j] + sign * .00125))
                vectors.append(u); labels.append(f"jacobian_{j}_{name}")
        if self.completed + FINAL_VERIFICATION_HISTORIES <= self.c["max_histories"] and self.can_fit(2*N, final=True):
            rows = self.batch("final_jacobian", vectors, labels)
            self.write_jacobian(rows, anchor)
        else:
            write_json(self.root/"jacobian_diagnostics.json", {"status": "not_run_budget", "anchor": anchor})
        selected = copy.deepcopy(self.best); base = list(selected["unit_vector"])
        if self.can_fit(2, final=True):
            self.repeat_origin = selected
            repeats = self.batch("final_repeats", [base, base], ["selected_repeat_1", "selected_repeat_2"], smoke=True)
            if len(repeats) != 2: raise RuntimeError("Two final repetitions required")
            origin = Path(selected["summary"]).parent
            graphs = sorted((origin/"standard_diagnostics").glob("*.png"))
            if len(graphs) != 17: raise RuntimeError("Selected graph packet is incomplete")
            for row in repeats:
                dest = Path(row["summary"]).parent
                adapter.compare_reference(dest, origin/"summary.json")
                for graph in graphs: adapter.verify(dest/"standard_diagnostics"/graph.name, digest(graph))
            write_json(self.root/"final_verification.json", {"status": "pass", "exact_repeats": 2,
                "exact_standard_pngs": 17, "selected": selected, "production_promoted": False})
        self.reports()
        self.run_finalizer_if_pinned()
        self.summary(self.completion_status())

    def completion_status(self):
        if not (self.root/"final_verification.json").exists():
            return "best_valid_without_final_repeats"
        if self.c.get("finalizer_driver"):
            path = self.root/"finalizer_status.json"
            if not path.exists() or adapter.read_json(path).get("status") != "complete":
                return "calibration_verified_policy_incomplete"
        return "complete_verified"

    def write_jacobian(self, rows, anchor):
        import numpy as np
        fit0 = adapter.read_csv(Path(anchor["summary"]).parent/"target_fit_long.csv")
        moments = [x["moment"] for x in fit0]
        weights = {x["moment"]: float(x["weight"]) for x in fit0}
        center_values = {x["moment"]: float(x["model"]) for x in fit0}
        center = anchor["unit_vector"]; table = {}
        for row in rows:
            j = int(row["label"].split("_")[1])
            values = {x["moment"]: float(x["model"]) for x in adapter.read_csv(Path(row["summary"]).parent/"target_fit_long.csv")}
            table.setdefault(j, []).append((row["unit_vector"][j], values))
        matrix = np.full((len(moments), N), np.nan); methods = []; widths = []
        for j in range(N):
            pairs = [(center[j], center_values), *table.get(j, [])]
            pairs.sort(key=lambda x: x[0]); lo, hi = pairs[0], pairs[-1]
            width = hi[0] - lo[0]; widths.append(width)
            if width <= 1e-12:
                methods.append("missing_or_collapsed"); continue
            methods.append("central" if lo[0] < center[j] < hi[0] else "one_sided")
            for i, moment in enumerate(moments):
                matrix[i, j] = math.sqrt(weights[moment])*(hi[1][moment]-lo[1][moment])/width
        names = [x["name"] for x in self.c["search_domain"]]
        write_csv(self.root/"jacobian_weighted_moments.csv", [{"moment": m, **{names[j]: float(matrix[i,j]) for j in range(N)}} for i,m in enumerate(moments)])
        meta = {"status": "complete" if np.isfinite(matrix).all() else "incomplete", "anchor": anchor,
                "methods": dict(zip(names, methods)), "unit_widths": dict(zip(names, widths)),
                "selected_may_include_better_jacobian_probe": True,
                "derivative_units": "sqrt(weight) times model moment per transformed unit coordinate"}
        if np.isfinite(matrix).all():
            singular = np.linalg.svd(matrix, compute_uv=False)
            meta.update(singular_values=singular.tolist(), numerical_rank=int(np.linalg.matrix_rank(matrix)),
                relative_rank_1e6=int(np.sum(singular > singular[0]*1e-6)),
                relative_rank_1e3=int(np.sum(singular > singular[0]*1e-3)),
                condition_number=float(singular[0]/singular[-1]) if singular[-1] > 0 else None)
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(10,7), constrained_layout=True)
            normalized = matrix / np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
            im = ax.imshow(normalized, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
            ax.set_xticks(range(N), names, rotation=50, ha="right", fontsize=7)
            ax.set_yticks(range(len(moments)), moments, fontsize=7)
            ax.set_title("Supplemental: local target sensitivity (columns normalized)")
            fig.colorbar(im, ax=ax); fig.savefig(self.root/"jacobian_supplemental.png", dpi=180); plt.close(fig)
        write_json(self.root/"jacobian_diagnostics.json", meta)

    def run_finalizer_if_pinned(self):
        driver = self.c.get("finalizer_driver")
        if not driver: return
        adapter.verify(driver, self.c["finalizer_sha256"])
        remaining = self.finish - time.time()
        if remaining <= 60:
            write_json(self.root/"finalizer_status.json", {"status": "skipped_insufficient_time"}); return
        self.phase = "equilibrium_paths_and_pdf"
        out = Path(self.c["output_root"])/"equilibrium_path"
        with (self.root/"finalizer.log").open("w") as log:
            proc = subprocess.Popen([sys.executable, driver, "--selected-summary", self.best["summary"],
                "--outdir", str(out), "--contract", self.c["contract_path"]],
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            with self.lock: self.active[proc.pid] = proc
            try:
                while proc.poll() is None:
                    self.state("running")
                    if time.time() > self.finish - 30: raise RuntimeError("Finalizer exhausted remaining time")
                    try: proc.wait(timeout=min(60, max(1, self.finish-time.time()-30)))
                    except subprocess.TimeoutExpired: pass
            except BaseException:
                self.stop_active(); raise
            finally:
                with self.lock: self.active.pop(proc.pid, None)
        receipt_path = out/"equilibrium_receipt.json"
        receipt = adapter.read_json(receipt_path) if receipt_path.exists() else {}
        write_json(self.root/"finalizer_status.json", {"status": receipt.get("status", "failed"),
            "returncode": proc.returncode, "receipt": str(receipt_path),
            "receipt_sha256": digest(receipt_path) if receipt_path.exists() else None})
        if proc.returncode: raise RuntimeError("Pinned finalizer failed")
        validate_policy_receipt(receipt, self.c, smoke=False,
            selected_hashes={digest(self.best["summary"])})

    def summary(self,status):
        self.reports(); lines=["# Morning summary","",f"Status: {status}. No estimate is promoted.",""]
        if self.best: lines += [f"Best actual completed calibrated loss: {self.best['loss']:.10g}.",f"Selected receipt source: `{Path(self.best['summary']).parent / 'case_receipt.json'}`.","Complete target-fit and parameter tables: `all_target_fits.csv`, `all_parameters.csv`."]
        if self.rejects: lines += ["",f"Outstanding rejected numerical/economic proposals: {len(self.rejects)}; see `rejects_ledger.csv`."]
        if self.search_stop_reason:
            lines += ["", "Search proposals stopped after three consecutive timeouts; see `search_stop.json`.",
                "Completed cases remain eligible for exact repetition and policy verification. Unstarted proposals have no loss."]
        (self.root/"MORNING_SUMMARY.md").write_text("\n".join(lines)+"\n")
        self.state(status)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--contract",type=Path,required=True); p.add_argument("--contract-sha256",required=True); p.add_argument("--mode",choices=("smoke","search"),required=True)
    a=p.parse_args(); run=Search(verify_contract(a.contract,a.contract_sha256),a.mode)
    try: getattr(run,a.mode)()
    except BaseException as e:
        run.stop_active(); run.summary("failed"); run.state("failed", error=str(e)); raise

if __name__ == "__main__": main()
