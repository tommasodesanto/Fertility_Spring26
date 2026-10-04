"""Runtime adapter for normalized Cobb--Douglas housing shares.

This experiment uses
``Q=c**a*s**(1-a)/(a**a*(1-a)**(1-a))``.  It deliberately has no
reference-rent or ``alpha0`` numerator.  The production engine is never
edited: ``install`` temporarily replaces its imported preference callback and
restores every import when the context manager exits.
"""
from __future__ import annotations

import contextlib
import copy
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np

from production import calibration as canonical
from production import inputs as production_inputs
from production import reporting as production_reporting

ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT = "ces_normalized_shares_v1"
ALPHA0 = 0.733
DELTA_ALPHA = 0.0
SHARE_BOUNDS = (0.0, 0.25)
_CANONICAL_PARAMETERS = production_inputs.DEFAULT_PARAMETERS
PARAMETERS = {**{k: v for k, v in _CANONICAL_PARAMETERS.items() if k != "h_P"},
              "delta_alpha_jump": 0.07780442689806688,
              "delta_alpha": 0.03897536437154123}


def _canonical(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def contract() -> dict[str, Any]:
    """Return the isolated 14-row contract with family_rooms promoted to scored."""
    targets, target_pin, weight_pin, bounds = canonical._contract()
    targets = copy.deepcopy(targets)
    family = [row for row in targets if row["moment"] == "family_rooms"]
    if len(family) != 1 or family[0] != {"moment": "family_rooms", "target": "0.38509964969278165", "weight": "0.0", "role": "validation"}:
        raise RuntimeError("canonical family_rooms baseline row drift")
    family[0].update(weight="280.52808370152104", role="scored")
    result_bounds = dict(bounds)
    if "h_P" not in result_bounds:
        raise RuntimeError("canonical inherited-coordinate mapping no longer has h_P")
    result_bounds.pop("h_P")
    result_bounds["delta_alpha_jump"] = SHARE_BOUNDS
    result_bounds["delta_alpha"] = SHARE_BOUNDS
    if len(result_bounds) != 11:
        raise RuntimeError("experimental mapping must have exactly eleven coordinates")
    return dict(experiment=EXPERIMENT, target_fit=copy.deepcopy(targets),
                baseline_target_fingerprint=target_pin, baseline_weight_fingerprint=weight_pin,
                target_fingerprint=_canonical(targets),
                weight_fingerprint=_canonical(dict(base_contract=targets, multipliers={})),
                bounds=result_bounds, scored_count=sum(r["role"] == "scored" for r in targets),
                total_target_rows=len(targets), alpha0=ALPHA0, delta_alpha="estimated; bounds [0, .25]",
                contract_id="ces_normalized_jump_slope_family_rooms_v1",
                material_multiplier="e(m)^(sigma-1)*[a(m)^a(m)*(1-a(m))^(1-a(m))]^(sigma-1)")


def bind_parameters(P, grid, parameters: Mapping[str, float]):
    """Bind the experimental eleven-coordinate vector without a housing floor."""
    spec = contract(); bounds = spec["bounds"]
    if set(parameters) != set(bounds):
        raise ValueError("wrong normalized-share eleven-coordinate mapping")
    for key, value in parameters.items():
        if not math.isfinite(float(value)) or not bounds[key][0] <= float(value) <= bounds[key][1]:
            raise ValueError("out of experimental bounds: " + key)
    base = {k: float(v) for k, v in parameters.items() if k not in ("delta_alpha_jump", "delta_alpha")}
    # Preserve the production input mapper for all nine unchanged coordinates.
    base["h_P"] = 0.0
    Q = production_inputs.bind_parameters(P, grid, base)
    Q.alpha_cons = ALPHA0
    Q.delta_alpha_jump = float(parameters["delta_alpha_jump"])
    Q.delta_alpha = float(parameters["delta_alpha"])
    Q.hbar_first_child_jump = 0.0
    Q.hbar_child_rooms = 0.0
    Q.child_room_floor = False
    Q.preference_spec = "eqscale"
    Q.eqscale_form = "power"
    Q.compensated_child_housing_shares = False
    Q.normalized_ces_limit_shares = True
    Q.normalized_ces_limit_identity = EXPERIMENT
    production_inputs.validate_inputs(Q, grid)
    return Q


def load_inputs(parameters: Mapping[str, float] | None = None):
    """Load original fixed primitives/grid and only bind this experiment's vector."""
    P, grid = production_inputs.load_inputs()
    point = dict(PARAMETERS if parameters is None else parameters)
    return bind_parameters(P, grid, point), grid.copy()


def _children_at_home(P, n: int, cs: int) -> int:
    """Return the independent-count child state, rejecting shared-clock inputs."""
    from production.engine.parameters import independent_child_maturation_active

    if not independent_child_maturation_active(P):
        raise ValueError("normalized CES shares require child_state_mode='independent_count'")
    n, cs = int(n), int(cs)
    return cs if 0 <= cs <= n else 0


def normalized_share_callback(P, alpha, benefit, material_multiplier):
    """Apply curved benefits and the normalized-Cobb--Douglas multiplier.

    The multiplier is set for *every* family-state cell, including childless
    cells.  The canonical renter/owner ``chi`` mapping remains in the native
    kernels; this callback changes neither rent nor housing-service mapping.
    """
    if not bool(getattr(P, "normalized_ces_limit_shares", False)):
        raise RuntimeError("normalized-share callback installed without its explicit switch")
    curvature = float(P.child_benefit_curvature)
    if not 0.0 <= curvature < 1.0:
        raise ValueError("child-benefit curvature must be in [0,1)")
    jump = float(P.delta_alpha_jump)
    slope = float(P.delta_alpha)
    if not SHARE_BOUNDS[0] <= jump <= SHARE_BOUNDS[1] or not SHARE_BOUNDS[0] <= slope <= SHARE_BOUNDS[1]:
        raise ValueError("normalized-share jump/slope outside [0, .25]")
    if float(P.alpha_cons) != ALPHA0:
        raise ValueError("normalized-share contract fixes alpha0=.733")
    for n in range(int(P.n_parity)):
        for cs in range(int(P.n_child_states)):
            m = _children_at_home(P, n, cs)
            a = ALPHA0 if m == 0 else float(np.clip(ALPHA0 - jump - slope * m, .05, .95))
            alpha[n, cs] = a
            benefit[n, cs] = 0.0 if m == 0 else float(P.psi_child) * m ** (1.0 - curvature)
            e = ((2.0 + .7 * m) / 2.0) ** .7
            k = a ** a * (1.0 - a) ** (1.0 - a)
            material_multiplier[n, cs] = (e * k) ** (float(P.sigma) - 1.0)


def _utility_contract(P, grid, *, residual=None, target_fit=None) -> dict[str, Any]:
    """Durable receipt for the immutable normalized-share utility contract."""
    jump = float(P.delta_alpha_jump)
    receipt = dict(
        experiment=EXPERIMENT,
        contract_id=contract()["contract_id"],
        immutable_contract_fingerprint=_canonical(contract()),
        target_fingerprint=contract()["target_fingerprint"],
        weight_fingerprint=contract()["weight_fingerprint"],
        baseline_target_fingerprint=contract()["baseline_target_fingerprint"],
        baseline_weight_fingerprint=contract()["baseline_weight_fingerprint"],
        effective_input_fingerprint=canonical.effective_input_fingerprint(P, grid),
        alpha_childless=ALPHA0,
        alpha_parent=float(np.clip(ALPHA0 - jump - float(P.delta_alpha), .05, .95)),
        delta_alpha_jump=jump,
        delta_alpha=float(P.delta_alpha),
        sigma=float(P.sigma),
        equivalence_scale="((2 + 0.7*m) / 2)**0.7",
        material_multiplier=contract()["material_multiplier"],
        child_state_mode=str(getattr(P, "child_state_mode", "")),
        housing_floor=False,
        reference_rent=None,
    )
    if target_fit is not None:
        receipt["target_fit"] = target_fit
    if residual is not None:
        receipt["residual"] = [float(x) for x in np.asarray(residual)]
        receipt["residual_count"] = len(receipt["residual"])
        receipt["loss"] = float(np.asarray(residual) @ np.asarray(residual))
    return receipt


def _write_utility_contract(directory, P, grid, *, residual=None, target_fit=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    receipt = _utility_contract(P, grid, residual=residual, target_fit=target_fit)
    (directory / "utility_contract.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return receipt


class _Installation(contextlib.AbstractContextManager):
    def __init__(self):
        from production.engine import shared
        from production import equilibrium
        self.patches = ((shared, "apply_child_preferences", normalized_share_callback),
                        (equilibrium, "build_context", build_reporting_context))
        self.original = []
    def __enter__(self):
        for module, name, replacement in self.patches:
            existed = hasattr(module, name)
            self.original.append((module, name, existed, getattr(module, name, None)))
            setattr(module, name, replacement)
        return self
    def __exit__(self, *exc):
        for module, name, existed, original in reversed(self.original):
            if existed:
                setattr(module, name, original)
            else:
                delattr(module, name)
        return False


def install():
    """Return a reversible runtime-hook context manager; it makes no source edit."""
    return _Installation()


def _experimental_actual_parameters(native_actual):
    def actual(prepared, P, grid):
        values = dict(native_actual(prepared, P, grid))
        values["h_P"] = 0.0
        values["delta_alpha_jump"] = float(P.delta_alpha_jump)
        values["delta_alpha"] = float(P.delta_alpha)
        return values
    return actual


def build_reporting_context(P, grid, out, *, price_start, deadline, max_lifecycle, closure):
    """Build one actual native reporting context with accurate experiment rows."""
    Path(out).mkdir(parents=True, exist_ok=True)
    context = production_reporting.build_context(P, grid, out, price_start=price_start,
                                                  deadline=deadline, max_lifecycle=max_lifecycle,
                                                  closure=closure)
    native_actual = context["fp"].actual_parameters
    context["fp"].actual_parameters = _experimental_actual_parameters(native_actual)
    rows = copy.deepcopy(context["manifest"]["full_parameter_table"])
    if sum(r["parameter"] == "delta_alpha_jump" for r in rows) != 1:
        raise RuntimeError("native parameter report has no unique share-jump row")
    for row in rows:
        if row["parameter"] == "delta_alpha_jump":
            row.update(lower=str(SHARE_BOUNDS[0]), upper=str(SHARE_BOUNDS[1]),
                       status="free normalized-CES-limit consumption-share jump; alpha0 fixed .733",
                       reference_estimate=str(DELTA_ALPHA))
        elif row["parameter"] == "delta_alpha":
            row.update(lower=str(SHARE_BOUNDS[0]), upper=str(SHARE_BOUNDS[1]),
                       status="free normalized-CES-limit post-first-child consumption-share slope; alpha0 fixed .733",
                       reference_estimate=str(DELTA_ALPHA))
        elif row["parameter"] == "h_P":
            row.update(lower="0.0", upper="0.0", estimate="0.0",
                       reference_estimate="0.0",
                       status="fixed zero; housing floor removed in this experiment")
        elif row["parameter"] == "utility_reference_rent":
            row["status"] = "inactive legacy input; unused by normalized CES-limit utility"
        elif row["parameter"] == "alpha_cons":
            row["status"] = "fixed alpha0=.733 in normalized CES-limit share experiment"
    context["manifest"]["full_parameter_table"] = rows
    context["free_coordinates"] = ["delta_alpha_jump" if key == "h_P" else key
                                  for key in context["free_coordinates"]] + ["delta_alpha"]
    context["expected_parameters"] = context["fp"].actual_parameters(context["prepared"], P, grid)
    if {r["parameter"] for r in rows} != set(context["expected_parameters"]):
        raise RuntimeError("experimental parameter table/actual report mapping mismatch")
    context["experiment"] = dict(identity=EXPERIMENT, alpha0=ALPHA0, delta_alpha=float(P.delta_alpha),
        housing_floor=False, reference_rent=None,
        material_multiplier=contract()["material_multiplier"],
        housing_services="canonical production renter/owner chi mapping")
    _write_utility_contract(out, P, grid)
    return context


def residual_from_report(report_directory):
    """Authenticate raw native reports, then independently rescore all 14 rows."""
    _, raw_fits, rows = canonical.residual_from_report(report_directory)
    if [{k: row[k] for k in ("moment", "target", "weight", "role")} for row in raw_fits] != canonical._contract()[0]:
        raise RuntimeError("native raw target contract must remain the baseline 14-row CSV")
    fits = copy.deepcopy(raw_fits)
    family = [row for row in fits if row["moment"] == "family_rooms"]
    if len(family) != 1:
        raise RuntimeError("missing unique family_rooms report row")
    family[0].update(weight="280.52808370152104", role="scored")
    for row in fits:
        row["loss_contribution"] = str(float(row["weight"]) * float(row["gap"]) ** 2 if row["role"] == "scored" else 0.0)
    scored = [row for row in fits if row["role"] == "scored"]
    residual = np.asarray([math.sqrt(float(row["weight"])) * float(row["gap"]) for row in scored])
    loss = sum(float(row["loss_contribution"]) for row in scored)
    if residual.shape != (11,) or not np.isfinite(residual).all() or abs(float(residual @ residual) - loss) > 1e-8:
        raise RuntimeError("eleven finite scored residuals with exact experimental loss arithmetic required")
    return residual, fits, rows


def _write_experimental_fit(report_directory, fits):
    path = Path(report_directory) / "target_fit_experimental.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("moment", "target", "model", "gap", "weight", "loss_contribution", "role"))
        writer.writeheader()
        writer.writerows(fits)
    return path


def make_evaluator(out, lane, P, grid, deadline, price_start=None, *, solver=None,
                   target_fingerprint=None, weight_fingerprint=None, **_ignored):
    """Eleven-coordinate evaluator with native raw and experimental rescoring."""
    del lane
    spec = contract()
    if target_fingerprint not in (None, spec["target_fingerprint"]):
        raise RuntimeError("caller target fingerprint differs from unchanged contract")
    if weight_fingerprint not in (None, spec["weight_fingerprint"]):
        raise RuntimeError("caller weight fingerprint differs from unchanged contract")
    base_P, base_grid = (load_inputs() if P is None else (copy.deepcopy(P), np.asarray(grid).copy()))
    if solver is None:
        from production.equilibrium import solve_stationary_ge
        solver = solve_stationary_ge
    start = production_inputs.DEFAULT_PRICE if price_start is None else float(price_start)
    destination = Path(out)
    def evaluate(label, point, end):
        effective_end = min(float(end), float(deadline))
        if time.time() >= effective_end:
            return dict(status="budget_exhausted", lifecycle_solves=0)
        Q = bind_parameters(base_P, base_grid.copy(), point)
        effective_input_fingerprint = canonical.effective_input_fingerprint(Q, base_grid)
        with install():
            result = solver(Q, base_grid.copy(), out=destination / str(label), price_start=start,
                            budget_seconds=effective_end-time.time(), max_lifecycle=32,
                            closure="population_one")
        if result.get("status", "passed") != "passed":
            return {**{k: result[k] for k in ("status", "reason", "lifecycle_solves", "price_search") if k in result},
                    "effective_input_fingerprint": effective_input_fingerprint}
        residual, fits, rows = residual_from_report(result["report_directory"])
        report_directory = Path(result["report_directory"])
        _write_experimental_fit(report_directory, fits)
        utility = _write_utility_contract(report_directory, Q, base_grid, residual=residual, target_fit=fits)
        repeat_directory = report_directory.parent / "selected_repeat_final"
        if repeat_directory.is_dir():
            repeat_residual, repeat_fits, _ = residual_from_report(repeat_directory)
            _write_experimental_fit(repeat_directory, repeat_fits)
            if repeat_fits != fits or not np.array_equal(repeat_residual, residual):
                raise RuntimeError("independent repeat experimental rescore differs from root")
            _write_utility_contract(repeat_directory, Q, base_grid, residual=repeat_residual, target_fit=repeat_fits)
        h0 = float(result["closure"]["H0_derived"])
        if not canonical.H0_BOUNDS[0] <= h0 <= canonical.H0_BOUNDS[1]:
            return dict(status="inadmissible_numerical", lifecycle_solves=int(result["lifecycle_solves"]),
                        effective_input_fingerprint=effective_input_fingerprint)
        return dict(status="passed", report=str(result["report_directory"]), residual=residual.tolist(),
            loss=float(residual @ residual), objective=float(residual @ residual), target_fit=fits,
            parameter_table=rows, closure=result["closure"], H0_derived=h0, price=float(result["price"]),
            lifecycle_solves=int(result["lifecycle_solves"]), target_fingerprint=spec["target_fingerprint"],
            weight_fingerprint=spec["weight_fingerprint"], experiment=EXPERIMENT,
            effective_input_fingerprint=effective_input_fingerprint, utility_contract=utility)
    return evaluate


def preflight_contexts(out):
    """Build three successive real reporting contexts, zero lifecycle solves."""
    out = Path(out)
    if os.environ.get("CES_NORMALIZED_SHARES_STAGED_CONTEXT") != "1":
        result = dict(status="deferred_unavailable_dependency", actual_lifecycle_solves=0,
                      reason="actual frozen reporting dependencies are unavailable locally; run only in the staged dependency snapshot",
                      required_environment="CES_NORMALIZED_SHARES_STAGED_CONTEXT=1")
        out.mkdir(parents=True, exist_ok=True)
        (out / "context_preflight_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
        return result
    receipts = []
    from production import equilibrium
    for i in range(3):
        P, grid = load_inputs()
        with install():
            # Call the same global lookup that solve_stationary_ge uses.
            context = equilibrium.build_context(P, grid, out / f"context_{i}", price_start=production_inputs.DEFAULT_PRICE,
                                                deadline=time.time()+120., max_lifecycle=32, closure="population_one")
        receipt = dict(index=i, target_rows=len(contract()["target_fit"]), parameter_rows=len(context["manifest"]["full_parameter_table"]),
                       expected_parameter_rows=len(context["expected_parameters"]), experiment=context["experiment"],
                       lifecycle_solves=0)
        if receipt["target_rows"] != 14 or receipt["parameter_rows"] != 31 or receipt["expected_parameter_rows"] != 31:
            raise RuntimeError("actual reporting context does not retain 14/31 report shape")
        receipts.append(receipt)
    result = dict(status="passed_zero_solves", actual_lifecycle_solves=0, contexts=receipts)
    out.mkdir(parents=True, exist_ok=True)
    (out / "context_preflight_receipt.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
