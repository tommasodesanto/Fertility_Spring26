#!/usr/bin/env python3
"""Bounded overnight calibration controller for the simple fertility nest.

This is orchestration only.  It invokes the frozen refinement adapter, keeps
failed proposals out of ranking, and never changes the scientific model.
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import random
import signal
import shutil
import subprocess
import sys
import threading
import time
import traceback

import build_e5f_bounded_refinement_plan as planner
import run_e5f_bounded_calibration_refinement as adapter


SCHEMA = "e5f_simple_fertility_overnight_v1"
EXPECTED_IMPORTED_LOSS = 26.249682727266702
SEQUENTIAL_LOSS = 30.408527701170645


def _now():
    return time.monotonic()


def _digest(path):
    return adapter.digest(Path(path))


def _absolute(path, name):
    path = Path(path)
    if not path.is_absolute():
        raise RuntimeError(f"{name} must be an absolute path")
    return path


def classify_failure(returncode, timed_out, log_path, outdir):
    """Return an explicitly declared rejection kind, else raise later as fatal."""
    if timed_out:
        return "timeout"
    failure = Path(outdir) / "adapter_failure.json"
    if not failure.exists():
        return None
    payload = adapter.read_json(failure)
    text, error_type = str(payload.get("error", "")), str(payload.get("type", ""))
    if error_type == "InfeasibleThetaError":
        return "infeasible_theta"
    if text.startswith(("Housing market did not clear:", "Could not bracket the contemporaneous housing-market root",
                        "market gate failed")):
        return "market_nonconvergence"
    if text.startswith(("Old-steady-state fertility normalization is not bracketed",
                        "Old-steady-state fertility normalization missed tolerance")):
        return "old_fertility_no_bracket"
    if text.startswith(("The dated first-birth branch has zero treated mass",
                        "The dated first-birth branch has zero destination mass")):
        return "undefined_firstbirth_support"
    return None


def repair_unit(value):
    """Deterministic reflection/clip repair to the closed normalized cube."""
    if not math.isfinite(value):
        raise ValueError("non-finite normalized coordinate")
    value = abs(float(value))
    value = value % 2.0
    return value if value <= 1.0 else 2.0 - value


def inward_probe(unit, step=.002):
    """Move every coordinate toward the cube interior by exactly ``step``."""
    return [repair_unit(x + (step if x <= .5 else -step)) for x in unit]


def de_proposal(target, best, a, b, rng, *, rule, c=None):
    """Pure bounded current-to-best/1 or rand/1 proposal used by this controller."""
    c = a if c is None else c
    if not (len(target) == len(best) == len(a) == len(b) == len(c)):
        raise ValueError("DE vectors have different dimensions")
    scale = .55 + .25 * rng.random()
    cross = .65 + .25 * rng.random()
    base = target if rule == "current_to_best_1" else c
    raw = [(target[j] + scale * (best[j] - target[j]) + scale * (a[j] - b[j])
           if rule == "current_to_best_1" else base[j] + scale * (a[j] - b[j]))
           for j in range(len(target))]
    forced = rng.randrange(len(target))
    return [repair_unit(raw[j]) if (j == forced or rng.random() < cross) else target[j]
            for j in range(len(target))]


def initial_population(anchor, seed, size=23):
    """Deterministic 23-point design around the verified anchor, in all 11 coordinates."""
    if len(anchor) != 11 or size != 23:
        raise ValueError("overnight controller requires 23 11-dimensional population slots")
    rng = random.Random(seed)
    points = [("seed_anchor", list(anchor))]
    radii = (.005, .015, .04, .10)
    for radius in radii:
        for sign in (-1., 1.):
            direction = [sign if (j + int(radius * 1000)) % 2 else -sign for j in range(11)]
            points.append((f"nearby_r{radius:g}_{'plus' if sign > 0 else 'minus'}",
                           [repair_unit(x + radius * d) for x, d in zip(anchor, direction)]))
    # Mixed starts deliberately move the first-child jump/floor coordinates,
    # while the remaining points use broad deterministic directions.
    for n in range(size - len(points)):
        radius = .20 if n < 5 else (.10 if n < 9 else .04)
        direction = [rng.choice((-1., 1.)) for _ in anchor]
        direction[-2] = 1. if n % 2 else -1.
        direction[-1] = -1. if n % 2 else 1.
        points.append((f"mixed_{n + 1:02d}_r{radius:g}",
                       [repair_unit(x + radius * d) for x, d in zip(anchor, direction)]))
    return points


def _write_json(path, value):
    adapter.write_json(Path(path), value)


def _theta_center(summary, units, label):
    """Make the unchanged driver's panel-center payload for a normalized proposal."""
    import run_e5f_transition_calibration as calibration
    domain = summary["panel_design"]["domain"]
    if len(domain) != 11 or len(units) != 11:
        raise RuntimeError("The original 11-coordinate domain has changed")
    payload = {"old_psi_child": summary["old_psi_child"],
               "best_candidate": {"theta": copy.deepcopy(summary["best_candidate"]["theta"]),
                                  "old_psi_child": summary["old_psi_child"],
                                  "new_psi_child": summary["best_candidate"]["new_psi_child"]},
               "proposal": {"label": label, "unit_vector": list(units),
                            "status": "unevaluated proposal; no fitted moments or loss are supplied"}}
    theta = payload["best_candidate"]["theta"]
    for unit, spec in zip(units, domain):
        value = calibration.transform_unit(unit, spec["lower"], spec["upper"], spec["transform"])
        if spec["name"] == "beta_annual":
            theta["beta"] = value ** 4
        elif spec["name"] == "psi_child_change_2023":
            payload["best_candidate"]["new_psi_child"] = payload["old_psi_child"] + value
        else:
            theta[spec["name"]] = value
    return payload


def _new_plan(template, outdir, stage, cases, seed):
    outdir.mkdir(parents=True, exist_ok=False)
    plan = copy.deepcopy(template)
    plan.update(schema="e5f_bounded_refinement_v1", stage=stage, cases=[], input_sha256={},
                planner_sha256=_digest(planner.__file__), production_promoted=False)
    for number, (label, center, design) in enumerate(cases, 1):
        center_path = outdir / f"center_{number:03d}.json"
        _write_json(center_path, center)
        plan["cases"].append({"id": number, "label": label, "center": center_path.name,
                              "center_sha256": _digest(center_path), "panel_task_id": 1,
                              "panel_size": 1, "panel_design": "mixed", "radius": .0025,
                              "panel_seed": seed, "output": f"task_{number:03d}"})
    if len(plan["cases"]) > 23:
        raise RuntimeError("adapter plan exceeds 23 cases")
    path = outdir / "plan.json"
    _write_json(path, plan)
    return path, _digest(path)


def _run_process(command, log, timeout, active, lock):
    with Path(log).open("w") as stream:
        proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        with lock:
            active.add(proc)
        try:
            try:
                code = proc.wait(timeout=timeout)
                return {"returncode": code, "timeout": False}
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
                return {"returncode": 124, "timeout": True}
        finally:
            with lock:
                active.discard(proc)


class Controller:
    def __init__(self, contract_path, contract_sha):
        adapter.verify(contract_path, contract_sha)
        self.contract_path, self.contract_sha = Path(contract_path), contract_sha
        self.c = adapter.read_json(contract_path)
        self._validate_contract()
        self.root = _absolute(self.c["output"], "output")
        self.start = _now()
        wall = min(float(self.c["total_seconds"]), float(self.c.get("hardfinish_epoch", time.time() + self.c["total_seconds"])) - time.time())
        self.hard_deadline = self.start + max(0., wall)
        self.search_deadline = self.start + min(float(self.c["search_seconds"]), wall - float(self.c["final_reserve_seconds"]))
        self.active, self.lock, self.rows = set(), threading.Lock(), []
        self.launched = 0
        self.valid = []
        self.rejections = []
        self.search_rejection_waves = []
        self.search_stopped = False
        self.cancelled = False

    def _validate_contract(self):
        c = self.c
        required = ("output", "imported_repeat_plan", "comparison_reference", "source",
                    "source_sha256", "code_sha256", "imported_repeat_plan_sha256", "reference_sha256")
        if c.get("schema") != SCHEMA or any(key not in c for key in required):
            raise RuntimeError("Incomplete overnight controller contract")
        # ``smoke_plan`` is the current contract spelling.  Keep the former
        # internal spelling as a compatibility alias for prewritten launchers.
        smoke_path = c.get("smoke_plan", c.get("smoke_template_plan"))
        smoke_sha = c.get("smoke_plan_sha256", c.get("smoke_template_plan_sha256"))
        if not smoke_path or not smoke_sha:
            raise RuntimeError("Contract omits immutable orchestration smoke plan")
        if (c.get("workers"), c.get("total_seconds"), c.get("search_seconds"), c.get("final_reserve_seconds"),
            c.get("case_timeout_seconds"), c.get("max_evaluations"), c.get("seed")) != (23, 43200, 32400, 10800, 5400, 300, 20260908):
            raise RuntimeError("Contract does not pin the authorized overnight budget")
        for path, sha in c["code_sha256"].items():
            adapter.verify(_absolute(path, "code hash path"), sha)
        adapter.verify(_absolute(smoke_path, "smoke_plan"), smoke_sha)
        adapter.verify(_absolute(c["imported_repeat_plan"], "imported_repeat_plan"), c["imported_repeat_plan_sha256"])
        for path, sha in c["reference_sha256"].items():
            adapter.verify(_absolute(path, "reference hash path"), sha)
        adapter.verify(_absolute(c["source"], "source"), c["source_sha256"])
        if c["source_sha256"] != adapter.SOURCE:
            raise RuntimeError("Contract source fingerprint differs from frozen adapter source")
        if time.time() > c.get("launch_deadline_epoch", float("inf")):
            raise RuntimeError("Predeclared launch deadline exceeded")
        template = adapter.load_plan(Path(smoke_path), smoke_sha)
        if template["source"] != c["source"] or template["source_sha256"] != c["source_sha256"]:
            raise RuntimeError("Template/source contract mismatch")
        if c.get("sequential_reference_loss", SEQUENTIAL_LOSS) != SEQUENTIAL_LOSS:
            raise RuntimeError("Sequential comparison reference loss changed")

    def _heartbeat(self, stage, **more):
        _write_json(self.root / "heartbeat.json", {"stage": stage, "epoch": time.time(),
                    "elapsed_seconds": _now() - self.start, "launched": self.launched,
                    "valid": len(self.valid), "rejections": len(self.rejections), **more})

    def _record(self, row, *, update_best=True):
        self.rows.append(row)
        _write_json(self.root / "case_state.json", self.rows)
        _write_json(self.root / "latest_completed_case.json", row)
        if row.get("status") == "valid":
            self.valid.append(row)
            current = adapter.read_json(self.root / "best_so_far.json") if (self.root / "best_so_far.json").exists() else None
            if update_best and (current is None or row["loss"] < current["loss"]):
                _write_json(self.root / "best_so_far.json", row)

    def _import_best(self):
        plan_path = _absolute(self.c["imported_repeat_plan"], "imported_repeat_plan")
        plan, status, rows = planner.collect(plan_path, self.c["imported_repeat_plan_sha256"], require_complete=True)
        if len(rows) != 2 or any(abs(row["loss"] - EXPECTED_IMPORTED_LOSS) > 0 for row in rows):
            raise RuntimeError("Imported job 17152974 must contain two exact 26.249682727266702 repeats")
        for row in rows:
            receipt = adapter.read_json(Path(row["summary"]).parent / "case_receipt.json")
            if not receipt.get("reference", {}).get("exact_twelve_row_fit"):
                raise RuntimeError("Imported repeat lacks exact-reference receipt")
        return adapter.read_json(rows[0]["summary"]), rows

    def _stage(self, plan_path, plan_sha, stage, *, deadline, allow_rejections=True, update_best=True, search=False):
        plan = adapter.load_plan(plan_path, plan_sha)
        if search and self.search_stopped:
            return []
        ceiling = self.c["max_evaluations"] - (24 if search else 0)
        if self.launched + len(plan["cases"]) > ceiling:
            raise RuntimeError("Consumed launch-slot budget exceeded")
        if _now() + self.c["case_timeout_seconds"] > deadline:
            return []
        self.launched += len(plan["cases"])
        self._heartbeat(stage, running_cases=len(plan["cases"]))
        outputs = []
        stage_rejections = []
        stage_epoch = time.time()
        with futures.ThreadPoolExecutor(max_workers=23) as pool:
            tasks = {}
            for case in plan["cases"]:
                log = plan_path.parent / f"case_{case['id']:03d}.log"
                cmd = [sys.executable, str(Path(adapter.__file__).resolve()), "--plan", str(plan_path),
                       "--plan-sha256", plan_sha, "--case-id", str(case["id"])]
                tasks[pool.submit(_run_process, cmd, log, self.c["case_timeout_seconds"], self.active, self.lock)] = (case, log)
            try:
                pending = set(tasks)
                while pending:
                    done, pending = futures.wait(pending, timeout=60, return_when=futures.FIRST_COMPLETED)
                    for future in done:
                        case, log = tasks[future]
                        result = future.result()
                        out = plan_path.parent / case["output"]
                        row = {"stage": stage, "case_id": case["id"], "label": case["label"],
                               "plan": str(plan_path), "plan_sha256": plan_sha, **result}
                        if result["returncode"] == 0:
                            receipt = adapter.read_json(out / "case_receipt.json")
                            if receipt.get("status") != "complete" or receipt.get("plan_sha256") != plan_sha:
                                raise RuntimeError("missing or mismatched completed-case receipt")
                            for name, sha in receipt["artifact_sha256"].items():
                                adapter.verify(out / name, sha)
                            adapter.validate_result(out, plan, case)
                            row.update(status="valid", loss=receipt["loss"], summary=str(out / "summary.json"),
                                       receipt=str(out / "case_receipt.json"))
                            outputs.append(row)
                        else:
                            reason = classify_failure(result["returncode"], result["timeout"], log, out)
                            if not allow_rejections or reason is None:
                                self._reap()
                                raise RuntimeError(f"Fatal case failure {stage}/{case['id']}; see {log}")
                            row.update(status="rejected", rejection=reason, log=str(log))
                            self.rejections.append(row)
                            stage_rejections.append(row)
                        self._record(row, update_best=update_best)
                    for future in pending:
                        case, _ = tasks[future]
                        heartbeat = plan_path.parent / case["output"] / "heartbeat.json"
                        if ((heartbeat.exists() and time.time() - heartbeat.stat().st_mtime > 1800)
                                or (not heartbeat.exists() and time.time() - stage_epoch > 1800)):
                            self._reap()
                            raise RuntimeError("stale case heartbeat exceeds 30 minutes")
                    self._heartbeat(stage, running_cases=len(pending))
            except BaseException:
                self._reap()
                for task in tasks:
                    task.cancel()
                raise
        if search:
            self.search_rejection_waves.append(dict(count=len(plan["cases"]), rejected=stage_rejections))
            recent = self.search_rejection_waves[-2:]
            all_timeout = len(self.search_rejection_waves) >= 3 and all(
                len(wave["rejected"]) == wave["count"] and all(row["rejection"] == "timeout" for row in wave["rejected"])
                for wave in self.search_rejection_waves[-3:])
            half_rejected = len(recent) == 2 and all(len(wave["rejected"]) * 2 >= wave["count"] for wave in recent)
            if all_timeout or half_rejected:
                self.search_stopped = True
                _write_json(self.root / "search_stop.json", {"reason": "three_all_timeout_waves" if all_timeout else "two_half_rejection_waves",
                            "stage": stage, "epoch": time.time()})
        return outputs

    def _reap(self):
        with self.lock:
            processes = list(self.active)
        for proc in processes:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

    def _smoke(self, template, imported):
        plan_path = _absolute(self.c.get("smoke_plan", self.c.get("smoke_template_plan")), "smoke_plan")
        sha = self.c.get("smoke_plan_sha256", self.c.get("smoke_template_plan_sha256"))
        rows = self._stage(plan_path, sha, "smoke", deadline=self.search_deadline, allow_rejections=False)
        if len(rows) != 2 or not any(row["loss"] == EXPECTED_IMPORTED_LOSS for row in rows):
            raise RuntimeError("Smoke baseline did not exactly reproduce imported best")
        for row in rows:
            if row["loss"] == EXPECTED_IMPORTED_LOSS:
                receipt = adapter.read_json(Path(row["receipt"]))
                if not receipt.get("reference", {}).get("exact_twelve_row_fit"):
                    raise RuntimeError("Smoke baseline lacks exact-reference receipt")
        return rows

    def _search(self, template, imported):
        anchor = imported["panel_design"]["unit_vector"]
        population = initial_population(anchor, self.c["seed"])
        names = [spec["name"] for spec in imported["panel_design"]["domain"]]
        floor, jump = names.index("hbar_child_rooms"), names.index("hbar_first_child_jump")
        adjusted = []
        for label, unit in population:
            if label.startswith("mixed_"):
                number = int(label.split("_")[1])
                unit = list(unit)
                unit[jump] = repair_unit(anchor[jump] + (.20 if number % 2 else -.20))
                unit[floor] = repair_unit(anchor[floor] + (-.10 if number % 2 else .10))
            adjusted.append((label, unit))
        population = adjusted
        case_data = [(label, _theta_center(imported, units, label), "initial_population") for label, units in population]
        path, sha = _new_plan(template, self.root / "search_initial", "initial_population", case_data, self.c["seed"])
        initial = self._stage(path, sha, "initial_population", deadline=self.search_deadline, search=True)
        archive = [(r, adapter.read_json(r["summary"])["panel_design"]["unit_vector"]) for r in self.valid]
        if not archive:
            return
        no_improve = 0
        rng = random.Random(self.c["seed"])
        for generation in range(1, 7):
            if self.search_stopped:
                break
            if _now() + self.c["case_timeout_seconds"] > self.search_deadline:
                break
            before = min(r["loss"] for r, _ in archive)
            unique = {}
            for row, vector in sorted(archive, key=lambda x: (x[0]["loss"], x[0]["stage"], x[0]["case_id"])):
                unique.setdefault(tuple(vector), (row, vector))
            ranked = list(unique.values())
            if len(ranked) < 4:
                # Rejected slots deliberately have no objective and cannot be
                # donors; stop DE and retain the valid archived incumbent.
                break
            best = ranked[0][1]
            proposals = []
            for slot, (_, target) in enumerate(ranked[:23]):
                choices = [u for _, u in ranked if u != target]
                if len(choices) < 3:
                    break
                a, b, donor_c = rng.sample(choices, 3)
                rule = "current_to_best_1" if generation % 2 else "rand_1"
                proposal = de_proposal(target, best, a, b, rng, rule=rule, c=donor_c)
                label = f"de_g{generation:02d}_{slot + 1:02d}_{rule}"
                proposals.append((label, _theta_center(imported, proposal, label), "differential_evolution"))
            if not proposals:
                break
            path, sha = _new_plan(template, self.root / f"de_generation_{generation:02d}",
                                  f"de_generation_{generation:02d}", proposals, self.c["seed"] + generation)
            accepted = self._stage(path, sha, f"de_generation_{generation:02d}", deadline=self.search_deadline, search=True)
            for row in accepted:
                archive.append((row, adapter.read_json(row["summary"])["panel_design"]["unit_vector"]))
            after = min(r["loss"] for r, _ in archive)
            no_improve = no_improve + 1 if before - after < 1e-4 * max(1., abs(before)) else 0
            if no_improve >= 2:
                break
        self._polish(template, imported, archive)

    def _polish(self, template, imported, archive):
        for round_id in range(1, 3):
            if self.search_stopped or _now() + self.c["case_timeout_seconds"] > self.search_deadline:
                return
            best_row, best_unit = min(archive, key=lambda x: (x[0]["loss"], x[0]["stage"], x[0]["case_id"]))
            candidates = []
            for j in range(11):
                for sign in (-1., 1.):
                    unit = list(best_unit); unit[j] = min(1., max(0., unit[j] + sign * .0025))
                    if unit == best_unit:
                        continue
                    label = f"polish_coordinate_{j}_{'minus' if sign < 0 else 'plus'}"
                    candidates.append((label, _theta_center(imported, unit, label), "coordinate_polish"))
            stage = f"polish_{round_id}_coordinate"
            path, sha = _new_plan(template, self.root / stage, stage, candidates, self.c["seed"])
            accepted = self._stage(path, sha, stage, deadline=self.search_deadline, search=True)
            winners = {}
            for row in accepted:
                actual = adapter.read_json(row["summary"])["panel_design"]["unit_vector"]
                archive.append((row, actual))
                index = int(row["label"].split("_")[2])
                gain = best_row["loss"] - row["loss"]
                if gain > 0 and (index not in winners or gain > winners[index][0]):
                    winners[index] = (gain, actual[index] - best_unit[index])
            ranked = sorted(winners, key=lambda j: (-winners[j][0], j))
            patterns, seen = [], {tuple(best_unit)}
            for count in sorted(set((len(ranked), min(3, len(ranked)), min(6, len(ranked)))), reverse=True):
                if not count:
                    continue
                direction = [winners[j][1] if j in ranked[:count] else 0. for j in range(11)]
                for scale in (.5, 1., 2., 4.):
                    unit = [min(1., max(0., x + scale*d)) for x, d in zip(best_unit, direction)]
                    if tuple(unit) in seen:
                        continue
                    seen.add(tuple(unit))
                    label = f"polish_combined_top{count}_scale{scale:g}"
                    patterns.append((label, _theta_center(imported, unit, label), "combined_successful_coordinate_pattern"))
            if patterns and not self.search_stopped and _now() + self.c["case_timeout_seconds"] <= self.search_deadline:
                stage = f"polish_{round_id}_combined"
                path, sha = _new_plan(template, self.root / stage, stage, patterns[:12], self.c["seed"])
                accepted = self._stage(path, sha, stage, deadline=self.search_deadline, search=True)
                archive.extend((row, adapter.read_json(row["summary"])["panel_design"]["unit_vector"]) for row in accepted)

    def _final_checks(self, template):
        best = adapter.read_json(self.root / "best_so_far.json")
        selected = adapter.read_json(best["summary"])
        unit = selected["panel_design"]["unit_vector"]
        cases = []
        selected_plan = adapter.load_plan(Path(best["plan"]), best["plan_sha256"])
        selected_case = next(case for case in selected_plan["cases"] if case["id"] == best["case_id"])
        selected_center = adapter.read_json(Path(best["plan"]).parent / selected_case["center"])
        for i in range(2):
            cases.append((f"final_exact_repeat_{i+1}", copy.deepcopy(selected_center), "final_exact_repeat"))
        for j in range(11):
            for sign in (-1., 1.):
                probe = list(unit); probe[j] = min(1., max(0., probe[j] + sign * .0025))
                label = f"final_jacobian_{j}_{'minus' if sign < 0 else 'plus'}"
                cases.append((label, _theta_center(selected, probe, label), "final_jacobian_probe"))
        # The adapter caps a plan at 23; reserve two waves (23 + 1) within the
        # contract's three-hour final reserve.  Selection is already frozen.
        first, second = cases[:23], cases[23:]
        all_rows = []
        reference_root = self.root / "final_selected_reference"
        reference_root.mkdir()
        selected_dir = Path(best["summary"]).parent
        selected_case_dir = selected_dir / "cases" / selected["best_candidate"]["candidate"]
        (reference_root / "cases" / selected["best_candidate"]["candidate"]).mkdir(parents=True)
        for name in ("summary.json", "target_fit_long.csv"):
            shutil.copy2(selected_dir / name, reference_root / name)
        shutil.copy2(selected_case_dir / "transition_path.csv", reference_root / "cases" /
                     selected["best_candidate"]["candidate"] / "transition_path.csv")
        reference_summary = reference_root / "summary.json"
        for wave, subset in enumerate((first, second), 1):
            path, sha = _new_plan(template, self.root / f"final_wave_{wave}", f"final_wave_{wave}", subset, self.c["seed"])
            if wave == 1:
                plan = adapter.read_json(path)
                for case in plan["cases"][:2]:
                    # An exact repeat retains the selected candidate's original
                    # generator metadata as well as its unrounded center file.
                    for name in ("panel_task_id", "panel_size", "panel_design", "radius", "panel_seed"):
                        case[name] = selected_case[name]
                    case.update(reference=str(reference_summary), reference_sha256=_digest(reference_summary))
                _write_json(path, plan); sha = _digest(path)
            old_timeout = self.c["case_timeout_seconds"]
            self.c["case_timeout_seconds"] = 5100  # 85 min; reserve receipt/write margin.
            try:
                all_rows += self._stage(path, sha, f"final_wave_{wave}", deadline=self.hard_deadline, update_best=False)
            finally:
                self.c["case_timeout_seconds"] = old_timeout
        repeats = [r for r in all_rows if r["label"].startswith("final_exact_repeat")]
        verified = len(repeats) == 2 and repeats[0]["loss"] == repeats[1]["loss"] == best["loss"]
        jacobian = {f"final_jacobian_{j}_{sign}": {"status": "not_completed"}
                    for j in range(11) for sign in ("minus", "plus")}
        for row in self.rows:
            if row["label"] in jacobian:
                jacobian[row["label"]] = {k: row[k] for k in ("status", "loss", "summary", "rejection") if k in row}
        _write_json(self.root / "identification_probe_status.json", jacobian)
        self._write_jacobian(selected, jacobian)
        return best, verified, jacobian

    def _write_jacobian(self, selected, probes):
        import numpy as np
        base_units = selected["panel_design"]["unit_vector"]
        best = adapter.read_json(self.root/"best_so_far.json")
        fit = adapter.read_csv(Path(best["summary"]).parent/"target_fit_long.csv")
        if len(fit) != 12 or len(selected["panel_design"]["domain"]) != 11:
            raise RuntimeError("Jacobian requires the complete 12 by 11 system")
        base = np.array([float(r["standardized_gap"]) for r in fit])
        matrix = np.full((12, 11), np.nan)
        methods, rows = [], []
        for j, spec in enumerate(selected["panel_design"]["domain"]):
            valid = []
            for sign in ("minus", "plus"):
                info = probes[f"final_jacobian_{j}_{sign}"]
                if info["status"] != "valid":
                    continue
                summary = adapter.read_json(info["summary"])
                probe_fit = adapter.read_csv(Path(info["summary"]).parent/"target_fit_long.csv")
                if len(probe_fit) != 12 or any((a["moment"], a["target"], a["weight"]) != (b["moment"], b["target"], b["weight"]) for a,b in zip(fit,probe_fit)):
                    raise RuntimeError("Jacobian target/weight provenance mismatch")
                valid.append((summary["panel_design"]["unit_vector"][j], np.array([float(r["standardized_gap"]) for r in probe_fit])))
            if len(valid)==2 and abs(valid[1][0]-valid[0][0])>1e-12:
                matrix[:,j]=(valid[1][1]-valid[0][1])/(valid[1][0]-valid[0][0]); method="two_sided_or_boundary_pair"
            else:
                distinct=[v for v in valid if abs(v[0]-base_units[j])>1e-12]
                if distinct:
                    matrix[:,j]=(distinct[0][1]-base)/(distinct[0][0]-base_units[j]);method="one_sided"
                else:
                    method="unavailable"
            methods.append(dict(parameter=spec["name"],method=method))
            for k,r in enumerate(fit):
                rows.append(dict(moment=r["moment"],parameter=spec["name"],derivative=float(matrix[k,j]) if np.isfinite(matrix[k,j]) else None,method=method))
        planner.write_csv(self.root/"weighted_moment_jacobian.csv",rows)
        complete=bool(np.isfinite(matrix).all())
        singular=np.linalg.svd(matrix,compute_uv=False) if complete else np.array([])
        _write_json(self.root/"jacobian_diagnostic.json",dict(complete_matrix=complete,columns=methods,
            rank=int(np.linalg.matrix_rank(matrix)) if complete else None,
            singular_values=singular.tolist(),condition_number=float(singular[0]/singular[-1]) if complete and singular[-1]>0 else None,
            interpretation="Weighted residual derivatives in normalized parameter coordinates; local sensitivity only, not proof of global identification."))

    def run(self):
        self.root.mkdir(parents=True, exist_ok=False)
        smoke_path = Path(self.c.get("smoke_plan", self.c.get("smoke_template_plan")))
        smoke_sha = self.c.get("smoke_plan_sha256", self.c.get("smoke_template_plan_sha256"))
        template = adapter.load_plan(smoke_path, smoke_sha)
        if template.get("choice_model") != "fertility_nest" or not template.get("suppress_plots"):
            raise RuntimeError("Template must pin simple fertility nest and no figures")
        imported, imported_rows = self._import_best()
        _write_json(self.root / "contract_receipt.json", {"contract_sha256": self.contract_sha,
                    "imported_repeats": imported_rows, "selection_seed_loss": EXPECTED_IMPORTED_LOSS,
                    "no_policy": True, "no_figures": True})
        self._smoke(template, imported)
        self._search(template, imported)
        best, repeats_verified, jacobian = self._final_checks(template)
        reference = adapter.read_json(Path(self.c["comparison_reference"]) / "summary.json")
        if reference["best_candidate"]["transition_loss"] != SEQUENTIAL_LOSS:
            raise RuntimeError("Sequential reference summary loss changed")
        selected_dir = Path(best["summary"]).parent
        selected_fit = adapter.read_csv(selected_dir / "target_fit_long.csv")
        reference_fit = adapter.read_csv(Path(self.c["comparison_reference"]) / "target_fit_long.csv")
        if len(selected_fit) != 12 or len(reference_fit) != 12:
            raise RuntimeError("Final comparison lacks the full 12-row target table")
        comparison = []
        for nested, sequential in zip(selected_fit, reference_fit):
            if any(nested[k] != sequential[k] for k in ("moment", "target", "weight")):
                raise RuntimeError("Final comparison changes target/weight contract")
            comparison.append({"moment": nested["moment"], "target": nested["target"], "weight": nested["weight"],
                               "sequential_model": sequential["model"], "nested_model": nested["model"],
                               "sequential_gap": sequential["gap"], "nested_gap": nested["gap"],
                               "sequential_loss": sequential["loss_contribution"], "nested_loss": nested["loss_contribution"]})
        planner.write_csv(self.root / "comparison_target_fits.csv", comparison)
        shutil.copy2(selected_dir / "parameter_table.csv", self.root / "selected_parameter_table.csv")
        shutil.copy2(Path(self.c["comparison_reference"]) / "parameter_table.csv", self.root / "sequential_parameter_table.csv")
        report = ["# Simple fertility-nest overnight controller", "", f"Selected loss: {best['loss']}",
                  f"Sequential fixed reference loss: {SEQUENTIAL_LOSS}", f"Valid evaluations: {len(self.valid)}",
                  f"Consumed launch slots: {self.launched}", f"Declared rejections: {len(self.rejections)}",
                  f"Exact final repeats verified: {repeats_verified}", "",
                  "Full target comparison: `comparison_target_fits.csv`.",
                  "Selected parameters and bounds: `selected_parameter_table.csv`.",
                  "Sequential parameters and bounds: `sequential_parameter_table.csv`.",
                  "Local sensitivity: `weighted_moment_jacobian.csv` and `jacobian_diagnostic.json`.",
                  "This bounded search does not establish a global optimum. No policy simulations were run.", "",
                  "Selection was frozen before final Jacobian probes."]
        (self.root / "FINAL_READOUT.md").write_text("\n".join(report) + "\n")
        _write_json(self.root / "final_receipt.json", {"status": ("complete" if repeats_verified and all(v["status"] == "valid" for v in jacobian.values())
                               else "selected_verified_partial_identification" if repeats_verified else "partial_final_verification"), "contract_sha256": self.contract_sha,
                    "selected": best, "selected_verified_by_two_exact_repeats": repeats_verified,
                    "identification_probe": jacobian, "coverage": {"launch_slots": self.launched,
                    "valid": len(self.valid), "rejections": self.rejections}, "sequential_reference_loss": SEQUENTIAL_LOSS,
                    "selected_artifact_sha256": {name: _digest(selected_dir / name) for name in
                    ("summary.json", "target_fit_long.csv", "parameter_table.csv")},
                    "no_policy": True, "no_figures": True, "production_promoted": False})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--contract-sha256", required=True)
    args = parser.parse_args()
    controller = None
    try:
        controller = Controller(args.contract.resolve(), args.contract_sha256)
        def interrupt(signum, frame):
            controller._reap()
            raise KeyboardInterrupt(f"controller received signal {signum}")
        signal.signal(signal.SIGTERM, interrupt)
        signal.signal(signal.SIGINT, interrupt)
        controller.run()
    except BaseException as error:
        if controller is not None:
            controller._reap()
            _write_json(controller.root / "failure.json", {"error": repr(error), "traceback": traceback.format_exc(),
                        "elapsed_seconds": _now() - controller.start, "launched": controller.launched})
        raise


if __name__ == "__main__":
    main()
