#!/usr/bin/env python3
"""Isolated rental-size wedge: five fixed-price household solves on the original checkpoint.

Design (fixed-price partial equilibrium; no GE, no stationary recalibration):

    cap6zero   hR_max=6,  slope=0     exact replay of the original control
    cap10zero  hR_max=10, slope=0     replay of the old factorial cap10/phi=.8/lambda=0 arm
    cap10s005  hR_max=10, slope=0.05  finding
    cap10s02   hR_max=10, slope=0.2   finding (solved in the smoke stage, reused)
    cap10s1    hR_max=10, slope=1.0   finding

Every renter pays the full outside-landlord cost

    C(h) = rent * h + slope * h * max(h - 6, 0),

with zero intercept and the knee fixed at six rooms.  Only the cap and the slope
differ across arms; chi, preferences, prices, phi=.8, lambda=0, the raw saved
``stationary_g_pre`` population, entry and fiscal contracts are untouched.

The lifecycle modes run only from a staged, hash-pinned experiment root that
holds (i) the reviewed rental-wedge source snapshot (frozen September core plus
the explicit isolated patch) under ``source/code/model`` and (ii) the nine
runtime helpers copied from the finance_dose_v1 runtime under
``code/model/tools``.  The active checkout is never importable.  Wedge arms
refuse to solve until ``source_manifest.rental_wedge_port_reviewed`` is true
and a passed original control receipt exists; production arms additionally
require the verified smoke receipts (cap6zero and cap10s02).

The experiment-owned budget audit (``audit_budget_rows``) is kept as a third,
independent budget check on realized renter rows.  The ported snapshot audits
(``run_e5f_matched_pf_smoke.dated_budget`` and
``run_e5f_independent_numerical_audit.budget_audit``) already charge the full
cost and are both applied to every solved arm, together with the independent
segment-wise saving oracle (``saving_audit``), the standard policy-array audit,
and the native mass/probability gates.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import json
import math
import os
import signal
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[3]
CHECKPOINT_SHA256 = "3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993"
SCHEMA = "e5f_isolated_rental_wedge_v1"
KNEE = 6.0
ZERO_INTERCEPT = 0.0
MAX_LIFECYCLE_SOLVES = 8
CORE = "intergen_eqscale_seq_optimized"
TOL = 1e-10
SAVING_DRAWS = 32
# Experiment gate on the independent saving oracle.  The retained audit itself
# has no hard gate: it raises only when the oracle is worse than the saved
# choice by more than 1e-8 and otherwise reports maximum_value_gain plus the
# sample shares above 1e-6, 1e-4 and 1e-3.  Both facts are written to receipts.
SAVING_GAIN_TOLERANCE = 1e-7
RETAINED_SAVING_AUDIT_DEFINITION = (
    "run_e5f_independent_numerical_audit.saving_audit: weighted random sample plus deterministic "
    "mass-bearing age/tenure boundary states (first, last, modal occupied); per-state exact "
    "per-segment maximization against fixed production continuation; hard error only if the "
    "oracle is worse than the saved choice by more than 1e-8; reports maximum_value_gain and the "
    "draw-weighted shares of value gain above 1e-6, 1e-4 and 1e-3 (boundary states carry zero "
    "sample multiplicity and enter only the maximum)."
)
BUDGET_MASS_TOLERANCE = 2e-10
BUDGET_EXCESS_TOLERANCE = 1e-9

REQUIRED_SOURCE_FILES = (
    "code/model/intergen_eqscale_seq_optimized/__init__.py",
    "code/model/intergen_eqscale_seq_optimized/parameters.py",
    "code/model/intergen_eqscale_seq_optimized/kernels.py",
    "code/model/intergen_eqscale_seq_optimized/solver.py",
    "code/model/intergen_eqscale_seq_optimized/utils.py",
    "code/model/intergen_eqscale_seq_optimized/joint_nested.py",
    "code/model/tools/run_e5f_matched_pf_smoke.py",
    "code/model/tools/run_dynamic_population_transition.py",
    "code/model/tools/run_e5f_open_population_transition.py",
    "code/model/tools/run_e5f_perfect_foresight_transition.py",
    "code/model/tools/run_e5f_independent_numerical_audit.py",
)

PORT_BLOCKERS = (
    "Frozen parameters.py has no rental_wedge_intercept/slope/knee fields or helper.",
    "Frozen kernels.py and solver.py have no wedge-aware renter optimizer.",
    "The active wedge is restricted to the Markov-income factored path and rejects exhaustive saving.",
    "Existing dated_budget, budget_audit, and boundary gates omit the wedge cost.",
    "The active path cannot be imported from this driver; a snapshot-local runner and module-origin checks are required.",
)

# Runtime helpers staged from the finance_dose_v1 runtime (nine pinned files).
RUNTIME_HELPERS = (
    "build_persistent_transitory_income_candidate",
    "run_e5f_income_candidate_calibration",
    "run_e5f_income_candidate_search",
    "run_e5f_income_overnight_search",
    "run_e5f_financing_factorial",
    "build_e5f_native_financing_report",
    "run_e5f_native_financing_diagnostic",
    "run_e5f_native_income_cohort_diagnostic",
    "run_e5f_native_rental_access_diagnostic",
)
# Helpers this driver actually imports (the factorial runner is staged for provenance only;
# its run_case compares baselines on phi/lambda/cap and its allowlist does not know the wedge).
RUNTIME_HELPERS_USED = (
    "run_e5f_native_rental_access_diagnostic",
    "run_e5f_native_income_cohort_diagnostic",
    "build_e5f_native_financing_report",
)
POLICY_NAMES = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price")
EXTRA_ARRAYS = ("g_pre", "g_post_fertility", "g_current", "births")
METRIC_KEYS = ("birth_flow", "first_birth_flow", "mean_rooms", "ownership", "renter_mass", "owner_mass", "market_residual", "pre_mass")
COHORT_KEYS = ("cumulative_explicit_births_per_initial_household", "first_births_per_initial_household", "first_birth_mean_age", "initial_mass", "rows")
AUDIT_KEYS = ("maximum_occupied_value_drop", "occupied_negative_steps", "lower_node_mass_at_negative_steps", "share_pre_choice_mass_at_negative_steps")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMBA_DISABLE_JIT", "NUMBA_CACHE_DIR")
OLD_POPULATION_LABEL = "common saved stationary pre-choice mass"
OLD_PREFERENCES_LABEL = "family checkpoint unchanged"
OLD_ENTRY_RULE = "native conditional entrant for each process"


class ContractError(ValueError):
    """Raised when an immutable input contract is absent or inconsistent."""


class ReviewRequiredError(RuntimeError):
    """Raised before any unreviewed wedge lifecycle solve."""


@dataclass(frozen=True)
class Case:
    name: str
    cap: float
    slope: float


CASES = (
    Case("cap6zero", 6.0, 0.0),
    Case("cap10zero", 10.0, 0.0),
    Case("cap10s005", 10.0, 0.05),
    Case("cap10s02", 10.0, 0.2),
    Case("cap10s1", 10.0, 1.0),
)
CASE_BY_NAME = {case.name: case for case in CASES}
SMOKE_CASES = ("cap6zero", "cap10s02")
PRODUCTION_CASES = ("cap10zero", "cap10s005", "cap10s1")
STAGE_CASES = {"smoke": SMOKE_CASES, "production": PRODUCTION_CASES}
CASE_TIMEOUT_SECONDS = {"cap6zero": 600, "cap10zero": 600, "cap10s005": 900, "cap10s02": 900, "cap10s1": 900}
STAGE_SECONDS = {"smoke": 2700, "production": 3600}
HOUSEHOLD_SOLVES_PLANNED = len(CASES)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def write_json(path: Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)


def hash_tree(root: Path, *, skip_suffixes: tuple[str, ...] = ()) -> dict[str, str]:
    root = Path(root)
    return {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*")) if p.is_file() and not p.name.endswith(skip_suffixes)}


def under(path: Path, root: Path) -> bool:
    try:
        Path(path).relative_to(root)
        return True
    except ValueError:
        return False


def stage_of(case_name: str) -> str:
    for stage, names in STAGE_CASES.items():
        if case_name in names:
            return stage
    raise ContractError(f"case is not assigned to a stage: {case_name}")


def _safe_relative(path: str) -> Path:
    relative = Path(path)
    if relative.is_absolute() or ".." in relative.parts:
        raise ContractError(f"source manifest path is not repository-relative: {path}")
    return relative


def load_source_manifest(source_root: Path, manifest_path: Path) -> dict[str, Any]:
    """Verify every hash in an explicitly supplied immutable source snapshot."""
    root = Path(source_root).resolve()
    if root in (ROOT.resolve(), (ROOT / "code/model").resolve()):
        raise ContractError("refusing the active checkout; supply an immutable source snapshot")
    if not root.is_dir():
        raise ContractError(f"source snapshot is not a directory: {root}")
    try:
        manifest = json.loads(Path(manifest_path).resolve().read_text())
    except FileNotFoundError as exc:
        raise ContractError(f"source manifest missing: {manifest_path}") from exc
    pins = manifest.get("source_files") or manifest.get("source_sha256")
    if not isinstance(pins, dict) or not pins:
        raise ContractError("source manifest must contain a source_files hash map")
    missing_required = sorted(set(REQUIRED_SOURCE_FILES) - set(pins))
    if missing_required:
        raise ContractError(f"source manifest omits required runtime files: {missing_required}")
    for relative_text, expected in pins.items():
        relative = _safe_relative(str(relative_text))
        path = root / relative
        if not path.is_file():
            raise ContractError(f"pinned source file is missing: {relative}")
        actual = sha256(path)
        if actual != str(expected):
            raise ContractError(f"source hash mismatch for {relative}: {actual} != {expected}")
    # Manifest keys are snapshot-root relative (code/model/...); the root is the snapshot root.
    present = {str(p.relative_to(root)) for p in (root / "code/model").rglob("*.py") if "__pycache__" not in p.parts}
    extra = sorted(present - set(pins))
    if extra:
        raise ContractError(f"unpinned Python files inside the source snapshot: {extra}")
    if manifest.get("checkpoint_sha256", CHECKPOINT_SHA256) != CHECKPOINT_SHA256:
        raise ContractError("source manifest names a different checkpoint")
    return {
        "source_root": str(root),
        "manifest": str(Path(manifest_path).resolve()),
        "manifest_sha256": sha256(Path(manifest_path)),
        "source_files": len(pins),
        "wedge_port_reviewed": bool(manifest.get("rental_wedge_port_reviewed", False)),
        "control_required_before_wedge": bool(manifest.get("control_required_before_wedge", True)),
        "schema": manifest.get("schema"),
        "changed_source_vs_frozen": changed_source_hashes(manifest),
    }


def changed_source_hashes(manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Record the explicit isolated patch: frozen base hash versus ported hash, never rewritten."""
    base = manifest.get("base_hashes_frozen") or {}
    ported = manifest.get("ported_hashes") or {}
    return {
        rel: {"frozen": base.get(rel), "ported": ported.get(rel), "changed": base.get(rel) != ported.get(rel)}
        for rel in sorted(manifest.get("relevant_files", set(base) | set(ported)))
    }


def validate_checkpoint(checkpoint: Path, selected_checkpoint: Path | None = None) -> dict[str, str]:
    """Validate the retained original checkpoint without loading model code."""
    checkpoint = Path(checkpoint).resolve()
    if not checkpoint.is_file():
        raise ContractError(f"checkpoint missing: {checkpoint}")
    actual = sha256(checkpoint)
    if actual != CHECKPOINT_SHA256:
        raise ContractError(f"checkpoint hash mismatch: {actual} != {CHECKPOINT_SHA256}")
    if selected_checkpoint is not None:
        try:
            selected = json.loads(Path(selected_checkpoint).resolve().read_text())
        except FileNotFoundError as exc:
            raise ContractError(f"selected checkpoint record missing: {selected_checkpoint}") from exc
        if selected.get("checkpoint_sha256") != CHECKPOINT_SHA256:
            raise ContractError("selected checkpoint record does not identify the retained checkpoint")
    return {"checkpoint": str(checkpoint), "checkpoint_sha256": actual}


def rental_cost(rent: float, rooms: float, slope: float, *, knee: float = KNEE, intercept: float = ZERO_INTERCEPT) -> float:
    """Return the renter's full outside-landlord cost at the realized rooms."""
    values = (rent, rooms, slope, knee, intercept)
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("rental cost inputs must be finite")
    if rooms < 0 or slope < 0 or knee < 0 or intercept < 0:
        raise ValueError("rental cost inputs must be weakly nonnegative")
    return float(rooms * (rent + intercept + slope * max(rooms - knee, 0.0)))


def audit_budget_rows(
    rows: Iterable[Mapping[str, Any]],
    *,
    slope: float,
    knee: float = KNEE,
    intercept: float = ZERO_INTERCEPT,
    tolerance: float = BUDGET_EXCESS_TOLERANCE,
) -> dict[str, float]:
    """Audit budget rows with the experiment-owned full cost.

    Renter rows require ``rent`` and realized ``rooms``.  Owner rows require
    ``housing_cost``.  A missing cost input is an error, so a gate cannot pass
    by silently dropping the wedge.  ``mass`` is used for the mass-weighted
    excess gate and defaults to one only for synthetic/unit tests.
    """
    if slope < 0 or intercept != ZERO_INTERCEPT or knee != KNEE:
        raise ValueError("the experiment fixes zero intercept and a six-room knee")
    bad_mass = 0.0
    weighted_positive_gap = 0.0
    maximum_occupied_excess = 0.0
    row_count = 0
    for row in rows:
        row_count += 1
        tenure = row.get("tenure")
        mass = float(row.get("mass", 1.0))
        if not math.isfinite(mass) or mass < 0:
            raise ValueError("budget row mass must be finite and nonnegative")
        if tenure in (0, "renter"):
            if "rent" not in row or "rooms" not in row:
                raise ContractError("renter budget row must include rent and realized rooms")
            cost = rental_cost(float(row["rent"]), float(row["rooms"]), slope, knee=knee, intercept=intercept)
        elif tenure in (1, "owner", "owner_cost"):
            if "housing_cost" not in row:
                raise ContractError("owner budget row must include housing_cost")
            cost = float(row["housing_cost"])
        else:
            raise ContractError(f"unknown tenure in budget row: {tenure!r}")
        resources = float(row["resources"])
        gap = float(row["consumption"]) + float(row["saving"]) + cost - resources
        if not all(math.isfinite(value) for value in (cost, resources, gap)):
            raise ValueError("budget row contains a nonfinite value")
        if gap > tolerance:
            bad_mass += mass
        weighted_positive_gap += mass * max(gap, 0.0)
        if mass > 1e-12:
            maximum_occupied_excess = max(maximum_occupied_excess, gap)
    summary = {
        "rows": float(row_count),
        "budget_excess_mass": bad_mass,
        "weighted_positive_budget_gap": weighted_positive_gap,
        "maximum_occupied_excess": maximum_occupied_excess,
        "budget_tolerance": float(tolerance),
        "slope": float(slope),
        "knee": float(knee),
        "intercept": float(intercept),
    }
    if bad_mass > BUDGET_MASS_TOLERANCE or maximum_occupied_excess > tolerance:
        raise RuntimeError(f"rental budget gate failed: {summary}")
    return summary


def dated_budget_gate(rows: Iterable[Mapping[str, Any]], *, slope: float) -> dict[str, float]:
    """Named entry point for the dated-period gate."""
    return audit_budget_rows(rows, slope=slope)


def independent_budget_gate(rows: Iterable[Mapping[str, Any]], *, slope: float) -> dict[str, float]:
    """Named entry point for the independent audit gate."""
    return audit_budget_rows(rows, slope=slope)


def control_spec(source: Mapping[str, Any], checkpoint: Mapping[str, str]) -> dict[str, Any]:
    """Describe the exact-baseline control that must pass before any wedge arm."""
    return {
        "status": "control_required_before_wedge",
        "case": "cap6zero",
        "cap": 6.0,
        "slope": 0.0,
        "intercept": 0.0,
        "knee": KNEE,
        "source": dict(source),
        "checkpoint": dict(checkpoint),
        "comparisons": [
            "all mandatory policy arrays versus the saved evaluation policy at 1e-10",
            "g_current and births versus the saved evaluation at 1e-10",
            "scalar metrics and cohort summaries versus every old cap6/phi=.8/lambda=0 receipt at 1e-10",
            "raw stationary_g_pre bitwise unchanged; mass, budget, value, probability, and saving-oracle gates",
            "module origins inside the experiment root only",
        ],
        "claim": "control reproduction is a prerequisite; no active off-equivalence is asserted here",
    }


def require_reviewed_wedge_port(case: Case, source: Mapping[str, Any], control_receipt: Mapping[str, Any] | None) -> None:
    """Refuse every non-control arm until the source port and original control are reviewed."""
    if case.name == "cap6zero":
        return
    if not source.get("wedge_port_reviewed", False):
        raise ReviewRequiredError(
            f"refusing {case.name}: rental wedge port is not marked reviewed; blockers: "
            + "; ".join(PORT_BLOCKERS)
        )
    if not control_receipt or control_receipt.get("status") != "passed_control":
        raise ReviewRequiredError("refusing wedge case until a passed original-checkpoint control receipt is supplied")


def arm_parameters(base: Mapping[str, Any], case_name: str, *, source: Mapping[str, Any], control_receipt: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Return the declared arm overrides on a plain mapping (contract-level; used by unit tests)."""
    try:
        case = CASE_BY_NAME[case_name]
    except KeyError as exc:
        raise ContractError(f"unknown experiment case: {case_name}") from exc
    require_reviewed_wedge_port(case, source, control_receipt)
    if "chi" not in base:
        raise ContractError("checkpoint parameters must expose chi for the fixed-preference contract")
    if case.name == "cap6zero":
        if float(base.get("hR_max", case.cap)) != case.cap:
            raise ContractError("original control must retain hR_max=6")
        return dict(base)
    out = dict(base)
    out.update({"hR_max": case.cap})
    if case.slope > 0:
        out.update({"rental_wedge_intercept": ZERO_INTERCEPT, "rental_wedge_slope": case.slope, "rental_wedge_knee": KNEE})
    return out


def planned_cli(experiment_root: str = "<EXPERIMENT_ROOT>") -> dict[str, Any]:
    driver = f"{experiment_root}/code/model/tools/run_e5f_isolated_rental_wedge.py"
    manifest = f"{experiment_root}/launch_manifest.json"
    results = f"{experiment_root}/results"
    common = f"--launch-manifest {manifest} --results {results}"
    case = lambda name: f"timeout --signal=TERM --kill-after=30s {CASE_TIMEOUT_SECONDS[name]}s python -B {driver} --mode case --case {name} {common} --deadline-epoch $DEADLINE"
    return {
        "smoke": [f"python -B {driver} --mode verify --stage smoke --phase before {common}",
                  f"DEADLINE=$(( $(date +%s) + {STAGE_SECONDS['smoke']} ))", *(case(n) for n in SMOKE_CASES),
                  f"python -B {driver} --mode check-smoke {common}",
                  f"python -B {driver} --mode verify --stage smoke --phase after {common}"],
        "production": [f"python -B {driver} --mode verify --stage production --phase before {common}",
                       f"python -B {driver} --mode check-smoke {common}",
                       f"DEADLINE=$(( $(date +%s) + {STAGE_SECONDS['production']} ))", *(case(n) for n in PRODUCTION_CASES),
                       f"python -B {driver} --mode combine {common}",
                       f"python -B {driver} --mode verify --stage production --phase after {common}"],
    }


def experiment_plan(source_root: Path, manifest: Path, checkpoint: Path | None, selected_checkpoint: Path | None = None) -> dict[str, Any]:
    source = load_source_manifest(source_root, manifest)
    selected = validate_checkpoint(checkpoint, selected_checkpoint) if checkpoint else {"checkpoint": None, "checkpoint_sha256": CHECKPOINT_SHA256, "note": "not verified locally; the remote stage verifies the retained checkpoint hash"}
    return {
        "schema": SCHEMA,
        "status": "ready_for_staging" if source["wedge_port_reviewed"] else "dryrun_only_port_unreviewed",
        "source": source,
        "checkpoint": selected,
        "cases": [case.__dict__ for case in CASES],
        "stages": {stage: list(names) for stage, names in STAGE_CASES.items()},
        "household_solves_planned": HOUSEHOLD_SOLVES_PLANNED,
        "max_lifecycle_solves": MAX_LIFECYCLE_SOLVES,
        "case_timeout_seconds": dict(CASE_TIMEOUT_SECONDS),
        "stage_seconds": dict(STAGE_SECONDS),
        "saving_audit": {"draws": SAVING_DRAWS, "gain_tolerance": SAVING_GAIN_TOLERANCE, "retained_definition": RETAINED_SAVING_AUDIT_DEFINITION},
        "port_blockers": list(PORT_BLOCKERS),
        "runtime_pin_strategy": "nine runtime helpers copied from the finance_dose_v1 runtime tools directory and pinned by the plan.remote.json hashes; the ported snapshot supplies the frozen core and the five patched/frozen helpers; the two sets live in separate directories and never overwrite each other",
        "planned_cli": planned_cli(),
        "control": control_spec(source, selected),
        "scope": "fixed-price partial equilibrium on the original checkpoint; explicit lifetime births are cohort accounting, not a stationary normalization or a frictionless benchmark; positive slopes are findings",
    }


# ----------------------------------------------------------------------------------------------
# Remote (staged experiment root) modes.  Everything below imports numpy and the pinned runtime.
# ----------------------------------------------------------------------------------------------

def mod(name: str) -> Any:
    return importlib.import_module(name)


def manifest_checks(m: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Fail-closed SHA256 checks of every pinned input, in both copied and original locations."""
    exp = Path(m["experiment_root"])
    rows: dict[str, dict[str, Any]] = {}

    def check(key: str, path: Any, expected: str) -> None:
        p = Path(path)
        actual = sha256(p) if p.is_file() else None
        rows[key] = {"path": str(p), "expected": expected, "actual": actual, "ok": actual == expected}

    check("checkpoint", m["checkpoint"]["path"], m["checkpoint"]["sha256"])
    check("summary_runtime", m["summary"]["remote_path"], m["summary"]["sha256"])
    check("plan_runtime", m["plan"]["remote_path"], m["plan"]["sha256"])
    check("plan_copy", m["plan"]["copy"], m["plan"]["sha256"])
    for name, pin in m["runtime_helpers"].items():
        check(f"helper_copy/{name}", exp / "code/model/tools" / name, pin)
        check(f"helper_runtime/{name}", Path(m["runtime_root"]) / "code/model/tools" / name, pin)
    check("driver", m["driver"]["path"], m["driver"]["sha256"])
    check("port_manifest_copy", m["port_manifest"]["copy"], m["port_manifest"]["sha256"])
    port = read_json(m["port_manifest"]["copy"])
    for rel, pin in port["source_files"].items():
        check(f"source_copy/{rel}", Path(m["source_copy"]) / rel, pin)
    for rel, entry in m["changed_source_vs_frozen"].items():
        if entry["frozen"] is not None:
            check(f"frozen_original/{rel}", Path(m["frozen_root"]) / rel, entry["frozen"])
    for name in m.get("generated_files", {}):
        check(f"generated/{name}", exp / name, m["generated_files"][name])
    if Path(__file__).resolve() != Path(m["driver"]["path"]).resolve():
        raise RuntimeError(f"running driver {Path(__file__).resolve()} is not the pinned driver {m['driver']['path']}")
    copied = sorted(p.name for p in (exp / "code/model/tools").glob("*.py"))
    expected = sorted(set(m["runtime_helpers"]) | {Path(m["driver"]["path"]).name})
    if copied != expected:
        raise RuntimeError(f"copied helper set {copied} differs from allowed {expected}")
    if sorted(m["runtime_helpers"]) != sorted(f"{n}.py" for n in RUNTIME_HELPERS):
        raise RuntimeError("launch manifest runtime helper set differs from the driver's pinned nine")
    bad = [k for k, r in rows.items() if not r["ok"]]
    if bad:
        raise RuntimeError("fail-closed hash check failed: " + ", ".join(bad))
    return rows


def check_sys_path(exp: Path) -> None:
    for entry in sys.path:
        p = (Path(entry) if entry else Path.cwd()).resolve()
        if ("Fertility_Spring26" in str(p) or entry == "") and not under(p, exp):
            raise RuntimeError(f"foreign sys.path entry: {entry!r}")


def install(m: Mapping[str, Any]) -> dict[str, Any]:
    """Put the copied ported source first on sys.path and import the pinned runtime helpers."""
    exp = Path(m["experiment_root"]).resolve()
    frozen = Path(m["source_copy"]).resolve() / "code/model"
    if ROOT.resolve() != exp:
        raise RuntimeError(f"driver ROOT {ROOT} is not the experiment root {exp}")
    if not (frozen / CORE / "solver.py").is_file():
        raise RuntimeError(f"copied ported core missing under {frozen}")
    check_sys_path(exp)
    for p in (frozen, frozen / "tools"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    check_sys_path(exp)
    handles = {
        "params": mod(CORE + ".parameters"),
        "model": mod(CORE + ".solver"),
        "primitive": mod("run_e5f_matched_pf_smoke"),
        "audit": mod("run_e5f_independent_numerical_audit"),
        "rental": mod("run_e5f_native_rental_access_diagnostic"),
        "cohort": mod("run_e5f_native_income_cohort_diagnostic"),
        "report": mod("build_e5f_native_financing_report"),
    }
    for key, module in handles.items():
        file = Path(getattr(module, "__file__", "")).resolve()
        expected_root = exp / "code/model/tools" if key in ("rental", "cohort", "report") else frozen
        if not under(file, expected_root):
            raise RuntimeError(f"{key} imported from {file}, expected under {expected_root}")
    for name in ("rental_wedge_active", "rental_wedge_total_cost", "build_debt_caps"):
        if not hasattr(handles["params"], name):
            raise RuntimeError(f"ported parameters module lacks {name}")
    for name in ("budget_audit", "saving_audit", "policy_array_audit", "standard_diagnostics"):
        if not hasattr(handles["audit"], name):
            raise RuntimeError(f"ported audit module lacks {name}")
    check_sys_path(exp)
    return handles


def import_origins(m: Mapping[str, Any]) -> dict[str, str]:
    """Every project module must come from the copied ported source or the copied runtime helpers."""
    exp = Path(m["experiment_root"]).resolve()
    frozen = Path(m["source_copy"]).resolve() / "code/model"
    helpers = exp / "code/model/tools"
    driver = Path(m["driver"]["path"]).resolve()
    rows: dict[str, str] = {}
    for name, module in list(sys.modules.items()):
        file = getattr(module, "__file__", None)
        if not file:
            continue
        p = Path(file).resolve()
        if "Fertility_Spring26" not in str(p) and not under(p, exp):
            continue
        rows[name] = str(p)
        top = name.split(".")[0]
        if name in ("__main__", "__mp_main__") or top == driver.stem:
            ok = p == driver
        elif top == CORE:
            ok = under(p, frozen / CORE)
        elif top in RUNTIME_HELPERS:
            ok = p == helpers / (top + ".py")
        else:
            ok = under(p, frozen)
        if not ok:
            raise RuntimeError(f"import origin violation: {name} loaded from {p}")
    for req in (CORE + ".solver", CORE + ".parameters", CORE + ".kernels", "run_e5f_matched_pf_smoke", "run_e5f_independent_numerical_audit", *RUNTIME_HELPERS_USED):
        if req not in rows:
            raise RuntimeError(f"required module not loaded: {req}")
    check_sys_path(exp)
    return rows


def packet(path: Path) -> Any:
    import gzip
    import pickle
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def changed(a: Any, b: Any) -> list[str]:
    import numpy as np

    def eq(x: Any, y: Any) -> bool:
        try:
            return np.array_equal(x, y) if isinstance(x, np.ndarray) or isinstance(y, np.ndarray) else bool(x == y)
        except Exception:
            return False
    return sorted(k for k in set(vars(a)) | set(vars(b)) if not eq(getattr(a, k, None), getattr(b, k, None)))


def policy_arrays(policy: Any) -> dict[str, Any]:
    import numpy as np
    missing = [n for n in POLICY_NAMES if getattr(policy, n, None) is None]
    if missing:
        raise ValueError(f"missing mandatory policy arrays: {missing}")
    return {n: np.asarray(getattr(policy, n)) for n in POLICY_NAMES}


def compare_exact(saved: Any, observed: Any) -> dict[str, float]:
    import numpy as np
    out = {}
    for n, x in policy_arrays(saved).items():
        y = policy_arrays(observed)[n]
        np.testing.assert_allclose(y, x, atol=TOL, rtol=0, err_msg=n)
        out[n] = float(np.max(np.abs(y.astype(float) - x.astype(float)), initial=0.0))
    return out


def population_audit(x: Mapping[str, Any]) -> dict[str, Any]:
    """Raw saved stationary_g_pre must equal the saved evaluation g_pre; the raw array is used."""
    import numpy as np
    raw = np.asarray(x["stationary_g_pre"])
    evaluated = np.asarray(x["evaluation"].g_pre)
    if raw.shape != evaluated.shape:
        raise ContractError(f"population shape mismatch: {raw.shape} vs {evaluated.shape}")
    if not np.isfinite(raw).all() or not np.isfinite(evaluated).all():
        raise ContractError("population contains nonfinite values")
    delta = evaluated - raw
    out = {"raw_vs_evaluation_l1": float(np.abs(delta).sum()), "raw_vs_evaluation_linf": float(np.abs(delta).max(initial=0.0)),
           "raw_vs_evaluation_changed_count": int(np.count_nonzero(delta)), "raw_vs_evaluation_total_gap": float(abs(evaluated.sum() - raw.sum())),
           "source": "raw stationary_g_pre", "shape": list(raw.shape), "mass": float(raw.sum())}
    if out["raw_vs_evaluation_l1"] > 1e-12 or out["raw_vs_evaluation_linf"] > 1e-12 or out["raw_vs_evaluation_total_gap"] > 1e-12:
        raise ContractError(f"population mismatch exceeds tolerance: {out}")
    return out


def build_arm(base: Any, case: Case, params: Any) -> tuple[Any, dict[str, Any]]:
    """Copy the checkpoint parameters and apply only the cap and (for positive slopes) the wedge fields."""
    import numpy as np
    if float(base.hR_max) != 6.0:
        raise ContractError(f"original checkpoint hR_max is {base.hR_max}, expected 6")
    if not np.allclose(np.asarray(base.phi, dtype=float), 0.8, atol=0, rtol=0):
        raise ContractError("original checkpoint phi is not 0.8 everywhere")
    if float(getattr(base, "lambda_d", 0.0)) != 0.0:
        raise ContractError("original checkpoint lambda_d is not zero")
    if getattr(base, "income_candidate_fingerprint", None) is not None:
        raise ContractError("original checkpoint carries a candidate income fingerprint")
    for name in ("rental_wedge_intercept", "rental_wedge_slope", "rental_wedge_knee"):
        if hasattr(base, name):
            raise ContractError(f"original checkpoint already carries {name}; the control would not be a frozen replay")
    P = copy.deepcopy(base)
    if case.cap != float(base.hR_max):
        P.hR_max = float(case.cap)
    if case.slope > 0:
        P.rental_wedge_intercept = ZERO_INTERCEPT
        P.rental_wedge_slope = float(case.slope)
        P.rental_wedge_knee = KNEE
    params.build_debt_caps(P)
    altered = changed(base, P)
    allowed = {"hR_max"} | ({"rental_wedge_intercept", "rental_wedge_slope", "rental_wedge_knee"} if case.slope > 0 else set())
    if set(altered) - allowed:
        raise ContractError(f"unexpected parameter changes for {case.name}: {sorted(set(altered) - allowed)}")
    if case.name == "cap6zero" and altered:
        raise ContractError(f"control changed parameters: {altered}")
    active = bool(params.rental_wedge_active(P))
    if active != (case.slope > 0):
        raise ContractError(f"rental_wedge_active={active} disagrees with slope={case.slope}")
    record = {"case": case.name, "hR_max": float(P.hR_max), "slope": float(case.slope), "intercept": ZERO_INTERCEPT, "knee": KNEE,
              "effective_wedge": {"intercept": float(getattr(P, "rental_wedge_intercept", 0.0)), "slope": float(getattr(P, "rental_wedge_slope", 0.0)),
                                  "knee": float(getattr(P, "rental_wedge_knee", 6.0)), "active": active},
              "changed_fields": altered, "chi": float(P.chi), "phi": np.asarray(P.phi, dtype=float).tolist(), "lambda_d": float(getattr(P, "lambda_d", 0.0)),
              "cost_function": "C(h) = rent*h + slope*h*max(h-6,0); intercept 0; knee 6"}
    return P, record


def old_rows(summary: Mapping[str, Any], cap: float) -> list[dict[str, Any]]:
    if summary.get("status") != "complete" or summary.get("design") != "dose":
        raise ContractError("old summary is not a complete dose summary")
    rows = [r for r in summary["cases"] if (float(r["phi"]), float(r["lambda"]), float(r["rental_cap"])) == (0.8, 0.0, float(cap))]
    if not rows:
        raise ContractError(f"no old receipt for phi=0.8 lambda=0 cap={cap}")
    return rows


def reconcile_contract(row: Mapping[str, Any], x: Any, m: Mapping[str, Any], plan: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    import numpy as np
    c = row["contract"]
    checks = {
        "checkpoint_path": (c["checkpoint"], m["checkpoint"]["path"]),
        "checkpoint_sha256": (c["checkpoint_sha256"], m["checkpoint"]["sha256"]),
        "source_root": (c["source_root"], m["frozen_root"]),
        "plan_source_root": (plan["source_root"], m["frozen_root"]),
        "source_manifest": (c["source_manifest"], plan["source_manifest_path"]),
        "family": (c["family"], "original"),
        "price": (c["price"], np.asarray(x["evaluation"].policy.price).tolist()),
        "initial_population_shape": (c["initial_population_shape"], list(np.asarray(x["stationary_g_pre"]).shape)),
        "candidate_fingerprint": (c.get("candidate_fingerprint"), getattr(x["parameters"], "income_candidate_fingerprint", None)),
        "population": (row["population"], OLD_POPULATION_LABEL),
        "preferences": (row["preferences"], OLD_PREFERENCES_LABEL),
        "entry_rule": (row["cohort"]["entry_rule"], OLD_ENTRY_RULE),
        "phi": (float(row["phi"]), 0.8),
        "lambda": (float(row["lambda"]), 0.0),
    }
    out = {k: {"old": a, "expected": b, "ok": a == b} for k, (a, b) in checks.items()}
    bad = [k for k, v in out.items() if not v["ok"]]
    if bad:
        raise ContractError(f"contract reconciliation failed for {row['label']}: " + ", ".join(bad))
    return out


def compare_scalars(new: Mapping[str, Any], old: Mapping[str, Any], keys: Iterable[str]) -> dict[str, dict[str, Any]]:
    rows = {}
    for k in keys:
        a, b = new.get(k), old.get(k)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
            d = abs(float(a) - float(b))
            rows[k] = {"new": a, "old": b, "abs_diff": d, "within": bool(d <= TOL)}
        else:
            rows[k] = {"new": a, "old": b, "abs_diff": None, "within": a == b}
    return rows


def reproduce(new: Mapping[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    per, passed = [], True
    for r in rows:
        metrics_cmp = compare_scalars(new["metrics"], r["metrics"], METRIC_KEYS)
        cohort_cmp = compare_scalars(new["cohort"], r["cohort"], COHORT_KEYS)
        ok = all(v["within"] for v in {**metrics_cmp, **cohort_cmp}.values())
        passed &= ok
        per.append({"old_label": r["label"], "metrics": metrics_cmp, "cohort": cohort_cmp, "metrics_and_cohort_within_tolerance": ok,
                    "budget_informational": compare_scalars(new["dated_budget"], r["budget"], tuple(r["budget"])),
                    "policy_audit_informational": compare_scalars(new["policy_audit"], r["policy_audit"], AUDIT_KEYS)})
    return {"tolerance": TOL, "gate": "scalar metrics and cohort summaries versus every matching old receipt", "passed": bool(passed), "rows": per}


def metrics(ev: Any, P: Any, report: Any) -> dict[str, float]:
    """Mirror run_e5f_financing_factorial.metrics (the old receipts' scalar definitions)."""
    import numpy as np
    gpre, gpost, gc = (np.asarray(a) for a in (ev.g_pre, ev.g_post_fertility, ev.g_current))
    first = float(gpre[:, :, :, :, :, 0, :].sum() - gpost[:, :, :, :, :, 0, :].sum())
    owner, rooms, _ = report._physical_housing({"g_current": gc, "hR_pol": np.asarray(ev.policy.hR_pol)}, {"H_own": np.asarray(P.H_own)})
    return {"birth_flow": float(np.asarray(ev.births).sum()), "first_birth_flow": first, "mean_rooms": float(rooms), "ownership": float(owner),
            "renter_mass": float(gc[:, 0, ...].sum()), "owner_mass": float(gc[:, 1:, ...].sum()),
            "market_residual": float(ev.relative_market_residual), "pre_mass": float(gpre.sum())}


def renter_budget_rows(ev: Any, P: Any, shared: Any, grid: Any, rent: float, model: Any) -> list[dict[str, Any]]:
    """Occupied realized renter rows for the experiment-owned budget audit (same resources as the ported audits)."""
    import numpy as np
    p, g = ev.policy, np.asarray(ev.g_current)
    grid = np.asarray(grid, dtype=float)
    rows: list[dict[str, Any]] = []
    for j in range(int(P.J)):
        for zz, z in enumerate(P.z_grid):
            income = model.income_at_state(P, 0, j, float(z))
            base_resources = P.R_gross * grid + income
            for nn in range(int(P.n_parity)):
                for cs in range(int(P.n_child_states)):
                    idx = (slice(None), 0, 0, j, zz, nn, cs)
                    gm = g[idx]
                    occupied = np.flatnonzero(gm > 1e-12)
                    if not len(occupied):
                        continue
                    grant = float(np.asarray(shared.gb_flat).reshape(-1)[nn + int(P.n_parity) * cs])
                    resources = base_resources + np.clip(grant - (P.R_gross * np.maximum(grid, 0) + income), 0, grant)
                    rooms, cons, save = np.asarray(p.hR_pol[idx]), np.asarray(p.c_pol[idx]), np.asarray(p.bp_pol[idx])
                    for b in occupied:
                        rows.append({"tenure": "renter", "rent": float(rent), "rooms": float(rooms[b]), "consumption": float(cons[b]),
                                     "saving": float(save[b]), "resources": float(resources[b]), "mass": float(gm[b])})
    if not rows:
        raise ContractError("no occupied renter states; the experiment-owned budget audit has nothing to check")
    return rows


def plots_manifest(out: Path) -> dict[str, Any]:
    pngs = sorted((out / "standard_diagnostics").glob("*.png"))
    if len(pngs) != 17:
        raise RuntimeError(f"expected 17 standard PNGs, found {len(pngs)} under {out / 'standard_diagnostics'}")
    return {"count": 17, "files": {p.name: {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for p in pngs}}


class Heartbeat:
    def __init__(self, path: Path, payload: dict[str, Any], interval: float = 60.0) -> None:
        self.path, self.payload, self.interval = path, dict(payload), interval
        self.started = time.time()
        self.phase = "starting"
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._loop, daemon=True)

    def _write(self) -> None:
        write_json(self.path, {**self.payload, "phase": self.phase, "elapsed_seconds": time.time() - self.started, "updated": time.time(),
                               "utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())})

    def _loop(self) -> None:
        while not self.stop.wait(self.interval):
            self._write()

    def __enter__(self) -> "Heartbeat":
        self._write()
        self.thread.start()
        return self

    def set_phase(self, phase: str) -> None:
        self.phase = phase
        self._write()

    def __exit__(self, *exc: Any) -> None:
        self.stop.set()
        self.phase = "stopped" if exc[0] is None else f"failed:{exc[0].__name__}"
        self._write()


def stage_partial_summary(results: Path, stage: str) -> dict[str, Any]:
    rows = {}
    for name in STAGE_CASES[stage]:
        receipt = results / name / "receipt.json"
        if receipt.is_file():
            r = read_json(receipt)
            rows[name] = {"status": r.get("status"), "metrics": r.get("metrics"), "cohort": {k: r.get("cohort", {}).get(k) for k in COHORT_KEYS},
                          "household_solves": r.get("household_solves"), "solve_wall_seconds": r.get("solve_wall_seconds"), "receipt_sha256": sha256(receipt)}
    summary = {"status": "complete" if len(rows) == len(STAGE_CASES[stage]) else "partial", "stage": stage, "completed": sorted(rows), "pending": [n for n in STAGE_CASES[stage] if n not in rows],
               "cases": rows, "job": os.environ.get("SLURM_JOB_ID"), "updated": time.time()}
    write_json(results / f"{stage}_partial_summary.json", summary)
    return summary


def load_control_receipt(results: Path) -> dict[str, Any] | None:
    path = results / "cap6zero" / "receipt.json"
    return read_json(path) if path.is_file() else None


def run_case(a: argparse.Namespace, m: Mapping[str, Any]) -> None:
    import numpy as np
    case = CASE_BY_NAME[a.case]
    stage = stage_of(case.name)
    results = Path(a.results)
    out = results / case.name
    if out.exists():
        raise RuntimeError(f"case output already exists (fresh directory required; no automatic retries): {out}")
    if a.deadline_epoch and time.time() > a.deadline_epoch:
        raise RuntimeError("stage time budget exhausted before case start")
    out.mkdir(parents=True)
    hashes = manifest_checks(m)
    source = load_source_manifest(Path(m["source_copy"]), Path(m["port_manifest"]["copy"]))
    selected = validate_checkpoint(Path(m["checkpoint"]["path"]))
    control_receipt = load_control_receipt(results)
    require_reviewed_wedge_port(case, source, control_receipt)
    if stage == "production":
        smoke = results / "smoke_check.json"
        if not smoke.is_file() or read_json(smoke).get("status") != "verified":
            raise ReviewRequiredError("production arm refused: smoke_check.json is missing or not verified")
    handles = install(m)
    origins_before = import_origins(m)
    plan = read_json(m["plan"]["copy"])
    old_summary = read_json(m["summary"]["remote_path"])
    rows = old_rows(old_summary, case.cap) if case.slope == 0 else []
    x = packet(Path(m["checkpoint"]["path"]))
    population = population_audit(x)
    P, arm = build_arm(x["parameters"], case, handles["params"])
    contract = {r["label"]: reconcile_contract(r, x, m, plan) for r in rows}
    grid = np.asarray(x["b_grid"])
    price = np.asarray(x["evaluation"].policy.price)
    rent = float(P.user_cost_rate * price[0])
    model_mod = handles["model"]
    write_json(out / "arm.json", {**arm, "stage": stage, "job": os.environ.get("SLURM_JOB_ID"), "started": time.time()})
    calls: list[float] = []
    original_solve = model_mod.solve_markov_income_at_prices

    def counted(*args: Any, **kwargs: Any) -> Any:
        calls.append(time.time())
        return original_solve(*args, **kwargs)

    # The counter stays installed through graphs and the cohort so any hidden extra solve is caught.
    model_mod.solve_markov_income_at_prices = counted
    try:
        status = _solve_and_gate(a, m, case, stage, results, out, x, P, arm, rows, contract, population, grid, rent, handles, hashes, source, selected, origins_before, calls)
    finally:
        model_mod.solve_markov_income_at_prices = original_solve
    write_json(results / "latest_completed.json", {"case": case.name, "stage": stage, "status": status, "receipt_sha256": sha256(out / "receipt.json"), "updated": time.time()})
    stage_partial_summary(results, stage)
    if status == "failed_reproduction":
        raise RuntimeError(f"{case.name} does not reproduce the old receipt within {TOL}; see {out / 'receipt.json'}")


def _solve_and_gate(a: argparse.Namespace, m: Mapping[str, Any], case: Case, stage: str, results: Path, out: Path, x: Any, P: Any, arm: dict[str, Any],
                    rows: list[dict[str, Any]], contract: dict[str, Any], population: dict[str, Any], grid: Any, rent: float, handles: dict[str, Any],
                    hashes: dict[str, Any], source: dict[str, Any], selected: dict[str, str], origins_before: dict[str, str], calls: list[float]) -> str:
    import numpy as np
    rental, cohort_mod, audit, report, model_mod = handles["rental"], handles["cohort"], handles["audit"], handles["report"], handles["model"]
    with Heartbeat(results / "heartbeat.json", {"case": case.name, "stage": stage, "job": os.environ.get("SLURM_JOB_ID")}) as hb:
        hb.set_phase("household_solve")
        started = time.time()
        ev, dated_budget, shared, model = rental.native_solve(x, P)
        solve_wall = time.time() - started
        if len(calls) != 1:
            raise RuntimeError(f"expected exactly one household solve, observed {len(calls)}")
        if model is not model_mod:
            raise RuntimeError("native_solve returned a solver module other than the ported one")
        captured = {n: np.array(getattr(ev.policy, n), copy=True) for n in POLICY_NAMES}
        captured.update({n: np.array(getattr(ev, n), copy=True) for n in EXTRA_ARRAYS})
        arrays = out / "policy_arrays.npz"
        np.savez_compressed(arrays, **captured)
        write_json(out / "solve_receipt.json", {"status": "solved_awaiting_gates", "household_solves": len(calls),
                   "solve_wall_seconds": solve_wall, "policy_arrays_sha256": sha256(arrays)})
        hb.set_phase("gates")
        if not np.array_equal(np.asarray(ev.g_pre), np.asarray(x["stationary_g_pre"])):
            raise RuntimeError("period mapper changed the raw saved stationary_g_pre")
        mass = rental.gates(ev)
        if float(dated_budget.get("budget_excess_mass", np.inf)) > BUDGET_MASS_TOLERANCE or float(dated_budget.get("maximum_occupied_excess", np.inf)) > BUDGET_EXCESS_TOLERANCE:
            raise RuntimeError(f"dated budget gate failed: {dated_budget}")
        audit_packet = {"evaluation": ev, "parameters": P, "b_grid": grid, "shared": shared}
        (out / "budget_audit").mkdir()
        independent_budget = audit.budget_audit(audit_packet, out / "budget_audit")
        if float(independent_budget["budget_excess_mass"]) > BUDGET_MASS_TOLERANCE or float(independent_budget["maximum_occupied_excess"]) > BUDGET_EXCESS_TOLERANCE:
            raise RuntimeError(f"independent budget gate failed: {independent_budget}")
        experiment_budget = audit_budget_rows(renter_budget_rows(ev, P, shared, grid, rent, model), slope=case.slope)
        pa = audit.policy_array_audit(audit_packet, out)
        if pa["occupied_negative_steps"]:
            raise RuntimeError(f"occupied value monotonicity gate failed: {pa}")
        hb.set_phase("saving_oracle")
        (out / "saving_audit").mkdir()
        saving = audit.saving_audit(audit_packet, out / "saving_audit", SAVING_DRAWS)
        saving_gate = {"maximum_value_gain": float(saving["maximum_value_gain"]), "tolerance": SAVING_GAIN_TOLERANCE,
                       "passed": bool(float(saving["maximum_value_gain"]) <= SAVING_GAIN_TOLERANCE), "draws": SAVING_DRAWS,
                       "unique_states": saving["unique_states"], "retained_definition": RETAINED_SAVING_AUDIT_DEFINITION,
                       "wedge_branch": "segment_oracle_wedge for renters" if case.slope > 0 else "segment_oracle (legacy renter flow)"}
        if not saving_gate["passed"]:
            raise RuntimeError(f"independent saving oracle gate failed: {saving_gate}")
        control_comparison = None
        if case.name == "cap6zero":
            hb.set_phase("control_identity")
            control_comparison = {"policy_max_abs_diff": compare_exact(x["evaluation"].policy, ev.policy)}
            for n in ("g_current", "births"):
                np.testing.assert_allclose(np.asarray(getattr(ev, n)), np.asarray(getattr(x["evaluation"], n)), atol=TOL, rtol=0, err_msg=n)
                control_comparison[n + "_max_abs_diff"] = float(np.max(np.abs(np.asarray(getattr(ev, n)) - np.asarray(getattr(x["evaluation"], n))), initial=0.0))
            control_comparison["tolerance"] = TOL
        hb.set_phase("standard_graphs")
        graphs = rental.standard_graphs(x, P, ev, shared, model, out)
        if graphs.get("status") != "completed" or graphs.get("count") != 17:
            raise RuntimeError(f"stable 17-plot packet failed: {graphs}")
        plots = plots_manifest(out)
        hb.set_phase("cohort")
        cohort = cohort_mod.run_cohort(x, P, ev.policy, case.name, out / "cohort")
        if cohort.get("status") != "completed" or int(cohort.get("rows", -1)) != int(P.J):
            raise RuntimeError(f"cohort receipt incomplete: {cohort}")
        if len(calls) != 1:
            raise RuntimeError(f"cohort or graphs triggered extra household solves: {len(calls)}")
        hb.set_phase("arrays")
        new = {"metrics": metrics(ev, P, report), "cohort": cohort, "dated_budget": dated_budget, "policy_audit": pa}
        repro = reproduce(new, rows) if rows else {"gate": "none: positive-slope arm is a finding, not a reproduction", "passed": None, "rows": []}
        status = "passed_control" if case.name == "cap6zero" else "completed"
        if rows and not repro["passed"]:
            status = "failed_reproduction"
        receipt = {"schema": SCHEMA, "status": status, "case": case.name, "stage": stage, "arm": arm, "old_labels": [r["label"] for r in rows],
                   "reproduction": repro, "contract_reconciliation": contract, "population": population, "mass_gates": mass,
                   "dated_budget": dated_budget, "independent_budget": independent_budget, "experiment_budget": experiment_budget,
                   "policy_audit": pa, "saving_oracle": saving, "saving_gate": saving_gate, "control_comparison": control_comparison,
                   "standard_graphs": graphs, "plots": plots, "cohort": cohort, "metrics": new["metrics"], "hashes": hashes, "hashes_after": manifest_checks(m), "source": source,
                   "checkpoint": selected, "import_origins_before_solve": origins_before, "import_origins": import_origins(m),
                   "household_solves": len(calls), "solve_wall_seconds": solve_wall, "case_wall_seconds": time.time() - started,
                   "policy_arrays": {"path": str(arrays), "sha256": sha256(arrays), "names": sorted(captured), "shapes": {n: list(v.shape) for n, v in captured.items()}},
                   "threads": {k: os.environ.get(k) for k in THREAD_VARS}, "job": os.environ.get("SLURM_JOB_ID"), "launch_manifest_sha256": sha256(a.launch_manifest),
                   "scope": "fixed-price partial equilibrium on the original checkpoint; no GE, entry, fiscal, preference, price or population change; explicit lifetime births are cohort accounting, not a stationary normalization or frictionless benchmark"}
        receipt["file_hashes"] = hash_tree(out)
        write_json(out / "receipt.json", receipt)
    return status


def check_receipt(results: Path, name: str, a: argparse.Namespace) -> dict[str, Any]:
    import numpy as np
    out = results / name
    r = read_json(out / "receipt.json")
    problems = []
    expected_status = "passed_control" if name == "cap6zero" else "completed"
    if r.get("status") != expected_status or r.get("case") != name:
        problems.append("status")
    g = r.get("standard_graphs", {})
    paths = [Path(p) for p in g.get("paths", [])]
    if g.get("count") != 17 or len(paths) != 17 or not all(p.is_file() for p in paths):
        problems.append("standard17")
    plots = r.get("plots", {}).get("files", {})
    if len(plots) != 17 or any(not Path(v["path"]).is_file() or sha256(Path(v["path"])) != v["sha256"] for v in plots.values()):
        problems.append("plots_manifest")
    arrays = out / "policy_arrays.npz"
    pa = r.get("policy_arrays", {})
    if not arrays.is_file() or sha256(arrays) != pa.get("sha256"):
        problems.append("policy_arrays")
    else:
        names = sorted(np.load(arrays).files)
        if names != pa.get("names") or not set(EXTRA_ARRAYS) <= set(names) or not set(POLICY_NAMES) <= set(names):
            problems.append("policy_array_names")
    if not r.get("hashes") or not all(v.get("ok") for v in r["hashes"].values()):
        problems.append("hashes")
    if not r.get("hashes_after") or not all(v.get("ok") for v in r["hashes_after"].values()):
        problems.append("hashes_after")
    if name in ("cap6zero", "cap10zero") and not r.get("reproduction", {}).get("passed"):
        problems.append("reproduction")
    if name == "cap6zero" and not r.get("control_comparison"):
        problems.append("control_comparison")
    if r.get("household_solves") != 1 or r.get("cohort", {}).get("status") != "completed":
        problems.append("cohort_or_solve_count")
    if not r.get("saving_gate", {}).get("passed"):
        problems.append("saving_gate")
    if r.get("policy_audit", {}).get("occupied_negative_steps"):
        problems.append("value_monotonicity")
    for key in ("dated_budget", "independent_budget", "experiment_budget"):
        b = r.get(key, {})
        if float(b.get("budget_excess_mass", np.inf)) > BUDGET_MASS_TOLERANCE or float(b.get("maximum_occupied_excess", np.inf)) > BUDGET_EXCESS_TOLERANCE:
            problems.append(key)
    if r.get("launch_manifest_sha256") != sha256(a.launch_manifest):
        problems.append("manifest")
    if not (out / "cohort" / "cohort_arrays.npz").is_file() or not (out / "cohort" / "initial_native_entry.npz").is_file():
        problems.append("cohort_arrays")
    stale = [rel for rel, h in r.get("file_hashes", {}).items() if not (out / rel).is_file() or sha256(out / rel) != h]
    if stale:
        problems.append("file_hashes:" + ",".join(stale[:5]))
    if problems:
        raise RuntimeError(f"receipt check failed for {name}: " + ", ".join(problems))
    return {"case": name, "status": r["status"], "receipt_sha256": sha256(out / "receipt.json")}


def check_smoke(a: argparse.Namespace, m: Mapping[str, Any]) -> None:
    results = Path(a.results)
    rows = {name: check_receipt(results, name, a) for name in SMOKE_CASES}
    write_json(results / "smoke_check.json", {"status": "verified", "cases": rows, "job": os.environ.get("SLURM_JOB_ID"), "checked": time.time(),
                                              "meaning": "exact baseline control and positive-wedge exact-loop smoke both passed every gate"})


def combine(a: argparse.Namespace, m: Mapping[str, Any]) -> None:
    import csv
    import numpy as np
    results = Path(a.results)
    labels = [c.name for c in CASES]
    checks = {name: check_receipt(results, name, a) for name in labels}
    receipts = {name: read_json(results / name / "receipt.json") for name in labels}
    loaded = {name: dict(np.load(results / name / "policy_arrays.npz")) for name in labels}
    entry = {name: np.asarray(np.load(results / name / "cohort" / "initial_native_entry.npz")["g_pre"]) for name in labels}
    names = sorted(set(POLICY_NAMES) | set(EXTRA_ARRAYS))
    for name, arrays in loaded.items():
        if sorted(arrays) != names:
            raise RuntimeError(f"mandatory array set mismatch: {name}")
    reference = labels[0]
    population_identity = {name: bool(np.array_equal(loaded[name]["g_pre"], loaded[reference]["g_pre"])) for name in labels}
    entry_identity = {name: bool(np.array_equal(entry[name], entry[reference])) for name in labels}
    if not all(population_identity.values()):
        raise RuntimeError(f"common population is not bitwise identical across arms: {population_identity}")
    if not all(entry_identity.values()):
        raise RuntimeError(f"entry cohort is not bitwise identical across arms: {entry_identity}")
    pairs = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            A, B, per = loaded[labels[i]], loaded[labels[j]], {}
            for n in names:
                if A[n].shape != B[n].shape:
                    per[n] = {"shape_mismatch": [list(A[n].shape), list(B[n].shape)], "within_tolerance": False}
                    continue
                finite = np.isfinite(A[n]) & np.isfinite(B[n])
                d = np.abs(A[n][finite].astype(float) - B[n][finite].astype(float))
                mx = float(d.max(initial=0.0))
                if not np.array_equal(A[n][~finite], B[n][~finite]):
                    mx = float("inf")
                per[n] = {"max_abs_diff": mx, "count_above_tolerance": int(np.count_nonzero(d > TOL)), "identical": bool(np.array_equal(A[n], B[n])), "within_tolerance": bool(mx <= TOL)}
            pairs[f"{labels[i]}__vs__{labels[j]}"] = {"arrays": per, "policies_coincide_within_tolerance": all(per[n]["within_tolerance"] for n in POLICY_NAMES),
                                                      "distributions_coincide_within_tolerance": all(per[n]["within_tolerance"] for n in EXTRA_ARRAYS)}
    fields = ["case", "cap", "slope", "status", *METRIC_KEYS, *COHORT_KEYS, "saving_max_value_gain", "dated_budget_max_excess", "independent_budget_max_excess", "solve_wall_seconds"]
    rows = []
    for name in labels:
        r = receipts[name]
        rows.append({"case": name, "cap": r["arm"]["hR_max"], "slope": r["arm"]["slope"], "status": r["status"], **{k: r["metrics"].get(k) for k in METRIC_KEYS},
                     **{k: r["cohort"].get(k) for k in COHORT_KEYS}, "saving_max_value_gain": r["saving_gate"]["maximum_value_gain"],
                     "dated_budget_max_excess": r["dated_budget"].get("maximum_occupied_excess"), "independent_budget_max_excess": r["independent_budget"].get("maximum_occupied_excess"),
                     "solve_wall_seconds": r.get("solve_wall_seconds")})
    with (results / "comparisons.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    total_solves = sum(int(receipts[n]["household_solves"]) for n in labels)
    if total_solves != HOUSEHOLD_SOLVES_PLANNED or total_solves > MAX_LIFECYCLE_SOLVES:
        raise RuntimeError(f"household solve count {total_solves} differs from the planned {HOUSEHOLD_SOLVES_PLANNED}")
    write_json(results / "combined_summary.json", {
        "schema": SCHEMA, "status": "complete", "cases": rows, "receipts": checks, "household_solves_total": total_solves,
        "population_identity_bitwise": population_identity, "entry_cohort_identity_bitwise": entry_identity, "array_names": names,
        "pairwise_policy_comparison": pairs, "identity_is_gate": {"population": True, "entry_cohort": True, "policies": False},
        "tolerance": TOL, "job": os.environ.get("SLURM_JOB_ID"), "launch_manifest_sha256": sha256(a.launch_manifest),
        "scope": "fixed-price partial equilibrium on the original checkpoint; slopes are findings; explicit lifetime births are cohort accounting, not a stationary normalization or a frictionless benchmark",
    })


def verify(a: argparse.Namespace, m: Mapping[str, Any]) -> None:
    exp = Path(m["experiment_root"]).resolve()
    results = Path(a.results)
    hashes = manifest_checks(m)
    source = load_source_manifest(Path(m["source_copy"]), Path(m["port_manifest"]["copy"]))
    selected = validate_checkpoint(Path(m["checkpoint"]["path"]))
    tree = {"source": hash_tree(exp / "source"), "code": hash_tree(exp / "code"),
            "inputs": {name: sha256(exp / name) for name in sorted(p.name for p in exp.iterdir() if p.is_file())}}
    receipt = {"status": "verified", "stage": a.stage, "phase": a.phase, "hashes": hashes, "hashes_after": manifest_checks(m), "source": source, "checkpoint": selected, "tree": tree,
               "tree_sha256": hashlib.sha256(json.dumps(tree, sort_keys=True).encode()).hexdigest(), "sys_path": list(sys.path), "job": os.environ.get("SLURM_JOB_ID"), "checked": time.time()}
    if a.phase == "before":
        install(m)
        # Serialization dependencies are not visible to a static import scan.
        # Resolve the actual checkpoint before allowing a household case.
        checkpoint_packet = packet(Path(m["checkpoint"]["path"]))
        receipt["checkpoint_deserialization"] = {"passed": True, "household_solves": 0}
        del checkpoint_packet
        receipt["import_origins"] = import_origins(m)
    else:
        before_path = results / f"verify_{a.stage}_before.json"
        if not before_path.is_file():
            raise RuntimeError(f"missing before-phase verification receipt: {before_path}")
        before = read_json(before_path)
        diffs = {}
        for section in ("source", "code", "inputs"):
            b, n = before["tree"][section], tree[section]
            changed_paths = sorted(set(b) ^ set(n)) + sorted(k for k in set(b) & set(n) if b[k] != n[k])
            if changed_paths:
                diffs[section] = changed_paths
        receipt["immutability"] = {"before_tree_sha256": before["tree_sha256"], "after_tree_sha256": receipt["tree_sha256"], "identical": not diffs, "differences": diffs}
        write_json(results / f"source_immutability_{a.stage}.json", receipt["immutability"])
        if diffs:
            raise RuntimeError(f"staged inputs changed during the {a.stage} stage: {diffs}")
    write_json(a.receipt or results / f"verify_{a.stage}_{a.phase}.json", receipt)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=("plan", "verify", "case", "check-smoke", "combine"), default="plan")
    parser.add_argument("--source-root", type=Path, help="plan mode: immutable ported source snapshot code/model directory")
    parser.add_argument("--manifest", type=Path, help="plan mode: port source_manifest.json")
    parser.add_argument("--checkpoint", type=Path, help="plan mode: optional local copy of the retained original checkpoint")
    parser.add_argument("--selected-checkpoint", type=Path)
    parser.add_argument("--output", type=Path, help="plan mode: write the plan JSON here instead of stdout")
    parser.add_argument("--launch-manifest", type=Path, help="remote modes: staged launch_manifest.json")
    parser.add_argument("--results", type=Path, help="remote modes: results directory inside the experiment root")
    parser.add_argument("--case", choices=tuple(CASE_BY_NAME))
    parser.add_argument("--stage", choices=tuple(STAGE_CASES))
    parser.add_argument("--phase", choices=("before", "after"), default="before")
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--receipt", type=Path)
    a = parser.parse_args(argv)
    if a.mode == "plan":
        if not a.source_root or not a.manifest:
            parser.error("plan mode requires --source-root and --manifest")
        plan = experiment_plan(a.source_root, a.manifest, a.checkpoint, a.selected_checkpoint)
        rendered = json.dumps(plan, indent=2, sort_keys=True) + "\n"
        if a.output:
            a.output.parent.mkdir(parents=True, exist_ok=True)
            a.output.write_text(rendered)
        else:
            print(rendered, end="")
        return 0
    if not a.launch_manifest or not a.results:
        parser.error(f"{a.mode} mode requires --launch-manifest and --results")
    m = read_json(a.launch_manifest)
    if m.get("schema") != "e5f_isolated_rental_wedge_launch_v1":
        raise ContractError(f"unexpected launch manifest schema: {m.get('schema')}")

    def stop(signum: int, frame: Any) -> None:
        raise KeyboardInterrupt(f"terminated by signal {signum}")
    signal.signal(signal.SIGTERM, stop)
    a.results.mkdir(parents=True, exist_ok=True)
    try:
        if a.mode == "verify":
            if not a.stage:
                raise ContractError("--stage required for verify")
            verify(a, m)
        elif a.mode == "case":
            if not a.case:
                raise ContractError("--case required")
            run_case(a, m)
        elif a.mode == "check-smoke":
            check_smoke(a, m)
        else:
            combine(a, m)
    except BaseException as exc:
        failure = {"status": "failed", "mode": a.mode, "case": a.case, "stage": a.stage, "phase": a.phase, "error": f"{type(exc).__name__}: {exc}",
                   "job": os.environ.get("SLURM_JOB_ID"), "failed_at": time.time()}
        write_json(a.results / f"failure_{a.mode}_{a.case or a.stage or ''}_{int(time.time())}.json", failure)
        if a.mode == "case" and a.case and (a.results / a.case).is_dir():
            write_json(a.results / a.case / "failed_receipt.json", failure)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
