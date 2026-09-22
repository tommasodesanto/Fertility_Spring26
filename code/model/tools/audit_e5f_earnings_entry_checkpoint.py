#!/usr/bin/env python3
"""Read-only audit of entrant distributions and feasibility in a native E5F checkpoint.

This loads one trusted native checkpoint and inspects saved arrays only. It never
calls the household optimizer, Bellman solver, or equilibrium routine.
"""
from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
import pickle
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _install_pickle_compatibility() -> None:
    """Allow known NumPy 2 / pathlib 3.12 pickle names in the project venv."""
    import pathlib

    if "numpy._core" not in sys.modules:
        sys.modules["numpy._core"] = np.core
        sys.modules["numpy._core.numeric"] = np.core.numeric
        sys.modules["numpy._core.multiarray"] = np.core.multiarray
    local = types.ModuleType("pathlib._local")
    for name in ("Path", "PosixPath", "WindowsPath"):
        if hasattr(pathlib, name):
            setattr(local, name, getattr(pathlib, name))
    sys.modules.setdefault("pathlib._local", local)


class _NativeCheckpointUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        if module == "pathlib._local":
            import pathlib
            return getattr(pathlib, name)
        return super().find_class(module, name)


def _load_checkpoint(path: Path, source_root: Path) -> Any:
    _install_pickle_compatibility()
    model_root = source_root / "code" / "model"
    tools_root = model_root / "tools"
    if not model_root.is_dir() or not tools_root.is_dir():
        raise FileNotFoundError(f"source root must contain code/model and tools: {source_root}")
    sys.path[:0] = [str(tools_root), str(model_root)]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as stream:
        return _NativeCheckpointUnpickler(stream).load()


def _get(obj: Any, name: str, default: Any = None) -> Any:
    return obj.get(name, default) if isinstance(obj, dict) else getattr(obj, name, default)


def _array_summary(values: Any) -> dict[str, Any] | None:
    if values is None:
        return None
    a = np.asarray(values)
    return {"shape": list(a.shape), "dtype": str(a.dtype),
            "finite": bool(np.all(np.isfinite(a))),
            "min": float(np.nanmin(a)) if a.size else None,
            "max": float(np.nanmax(a)) if a.size else None}


def _find_native_state(payload: Any) -> Any:
    if _get(payload, "parameters") is not None and _get(payload, "solution") is not None:
        return payload
    for key in ("old", "state", "checkpoint", "normalized_old"):
        candidate = _get(payload, key)
        if candidate is not None:
            try:
                return _find_native_state(candidate)
            except ValueError:
                pass
    raise ValueError("checkpoint does not contain a native {parameters, solution} state")


def _collapse_wealth(distribution: Any, age_axis: int | None = None) -> tuple[np.ndarray, float]:
    a = np.asarray(distribution, dtype=float)
    if a.ndim < 1:
        raise ValueError("saved distribution must have a wealth axis")
    if age_axis is not None:
        if a.ndim != 7 or age_axis != 3:
            raise ValueError(f"expected native 7D distribution with age axis 3, got {a.shape}")
        a = a[:, :, :, 0, :, :, :]
    if a.ndim not in (6, 7):
        raise ValueError(f"expected native 6D age slice or 7D distribution, got {a.shape}")
    raw_mass = float(np.sum(a))
    marginal = np.sum(a, axis=tuple(range(1, a.ndim)))
    total = float(np.sum(marginal))
    if total > 0:
        marginal = marginal / total
    return marginal, raw_mass


def _normalized_gaps(reference: np.ndarray | None, realized: np.ndarray | None) -> dict[str, Any]:
    if reference is None or realized is None:
        return {"available": False}
    if reference.shape != realized.shape:
        return {"available": False, "reason": "wealth-grid dimensions differ",
                "reference_shape": list(reference.shape), "realized_shape": list(realized.shape)}
    ref_total, real_total = float(reference.sum()), float(realized.sum())
    if ref_total <= 0 or real_total <= 0:
        return {"available": False, "reason": "zero marginal mass"}
    p, q = reference / ref_total, realized / real_total
    diff = q - p
    return {"available": True, "reference_total_before_normalization": ref_total,
            "realized_total_before_normalization": real_total,
            "normalized_l1": float(np.sum(np.abs(diff))),
            "total_variation_distance": float(0.5 * np.sum(np.abs(diff))),
            "max_abs_gap": float(np.max(np.abs(diff))),
            "mean_wealth_gap_model_units": None}


def _function_ast(path: Path, names: list[str]) -> dict[str, str | None]:
    tree = ast.parse(path.read_text())
    result: dict[str, str | None] = {}
    for name in names:
        node = next((x for x in tree.body if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)) and x.name == name), None)
        if node is None:
            result[name] = None
            continue
        node = ast.fix_missing_locations(node)
        if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:]
        result[name] = ast.dump(node, annotate_fields=True, include_attributes=False)
    return result


def audit(checkpoint: Path, source_root: Path) -> dict[str, Any]:
    checkpoint = checkpoint.resolve()
    source_root = source_root.resolve()
    payload = _load_checkpoint(checkpoint, source_root)
    state = _find_native_state(payload)
    P, solution = _get(state, "parameters"), _get(state, "solution")
    if P is None or solution is None:
        raise ValueError("native state is missing parameters or solution")
    model_root = source_root / "code" / "model"
    seq_solver = model_root / "intergen_eqscale_seq" / "solver.py"
    opt_solver = model_root / "intergen_eqscale_seq_optimized" / "solver.py"
    from intergen_eqscale_seq import solver as reference_solver

    # Read arrays directly; avoid copies of full value/distribution tensors.
    V = np.asarray(_get(solution, "V", _get(_get(state, "policy"), "V")))
    g = _get(solution, "g")
    beginning = _get(solution, "g_beginning_distribution")
    entry_censored_mass = _get(solution, "entry_censored_mass", _get(P, "_entry_censored_mass"))
    entry_total = _get(P, "_entry_total_mass", _get(solution, "entry_rate"))
    cutoff = float(getattr(reference_solver, "DEAD_VALUE_CUTOFF", -1e9))
    death_thresholds = (-1e6, -1e8, cutoff)

    realized_marginal = None
    realized_total = None
    distribution_timing = str(_get(solution, "current_distribution_timing", "unknown"))
    prechoice_array = None
    if beginning is not None:
        b = np.asarray(beginning)
        if b.ndim == 7:
            prechoice_array = b[:, :, :, 0, :, :, :]
            distribution_timing = "g_beginning_distribution age index 0 (beginning-of-period/pre-choice)"
    if prechoice_array is None and g is not None:
        b = np.asarray(g)
        if b.ndim == 7:
            prechoice_array = b[:, :, :, 0, :, :, :]
            distribution_timing = "g age index 0; saved current_distribution_timing=" + distribution_timing
    dead_occupancy: dict[str, Any] = {"available": False}
    if prechoice_array is not None:
        realized_marginal, realized_total = _collapse_wealth(prechoice_array)
        if V.ndim == 7 and V.shape[3] > 0:
            V0 = V[:, :, :, 0, :, :, :]
            if V0.shape == prechoice_array.shape:
                mass_total = float(np.sum(prechoice_array))
                rows: dict[str, Any] = {}
                for threshold in death_thresholds:
                    mask = V0 <= threshold
                    occupied = mask & (prechoice_array > 0)
                    rows[str(threshold)] = {
                        "value_state_count_all_choices": int(np.count_nonzero(mask)),
                        "value_state_share_all_choices": float(np.mean(mask)),
                        "occupied_state_count": int(np.count_nonzero(occupied)),
                        "occupied_probability_mass_conditional_on_age18": float(np.sum(prechoice_array[occupied]) / mass_total) if mass_total else None,
                        "occupied_raw_population_mass": float(np.sum(prechoice_array[occupied])),
                    }
                dead_occupancy = {"available": True, "value_shape_age18": list(V0.shape),
                                  "age18_prechoice_mass": mass_total,
                                  "thresholds": rows,
                                  "interpretation": "diagnostic occupancy only; additional thresholds are not new feasibility gates"}

    all_age_low_values = {"available": False}
    if beginning is not None and np.shape(beginning) == V.shape:
        mass = np.asarray(beginning)
        total = float(np.sum(mass))
        rows = {}
        for threshold in death_thresholds:
            occupied = (V <= threshold) & (mass > 0)
            affected = float(np.sum(mass[occupied]))
            rows[str(threshold)] = {
                "occupied_state_count": int(np.count_nonzero(occupied)),
                "occupied_raw_population_mass": affected,
                "occupied_population_share": affected / total if total else None,
            }
        all_age_low_values = {"available": True, "timing": "saved beginning-of-period distribution",
                              "total_mass": total, "thresholds": rows,
                              "interpretation": "Diagnostic of very low saved values; does not by itself establish choice contamination."}

    ref_marginal = None
    ref_source = None
    conditional = _get(P, "fixed_reference_entry_conditional")
    fixed_grid = _get(P, "fixed_reference_entry_grid", _get(state, "b_grid"))
    income_weights = None
    try:
        z_grid, income_weights, _ = reference_solver.income_transition_values(P)
    except Exception as exc:  # checkpoint may hold only an already-marginal receipt
        z_grid = _get(P, "z_grid")
        income_weights = _get(P, "z_weights")
        income_weight_error = str(exc)
    else:
        income_weight_error = None
    if conditional is not None and income_weights is not None:
        C = np.asarray(conditional, dtype=float)
        wz = np.asarray(income_weights, dtype=float).reshape(-1)
        if C.ndim == 2 and C.shape[1] == wz.size:
            ref_marginal = C @ wz
            ref_source = "P.fixed_reference_entry_conditional weighted by current income-state probabilities"
    if ref_marginal is None:
        idx = _get(solution, "entry_wealth_grid_indices")
        wt = _get(solution, "entry_wealth_grid_weights")
        bg = np.asarray(fixed_grid if fixed_grid is not None else _get(state, "b_grid", []), dtype=float)
        if idx is not None and wt is not None and bg.size:
            marginal = np.zeros(bg.size, dtype=float)
            ii, ww = np.asarray(idx, dtype=int), np.asarray(wt, dtype=float)
            if ii.size == ww.size and np.all((ii >= 0) & (ii < bg.size)):
                marginal[ii] += ww
                ref_marginal = marginal
                ref_source = "saved unconditional solution.entry_wealth_grid_indices/weights (conditional matrix absent)"
    gaps = _normalized_gaps(ref_marginal, realized_marginal)
    if gaps.get("available") and fixed_grid is not None:
        bg = np.asarray(fixed_grid, dtype=float)
        if bg.size == ref_marginal.size == realized_marginal.size:
            gaps["mean_wealth_gap_model_units"] = float(bg @ realized_marginal - bg @ ref_marginal)

    helper_names = ["income_at_state", "annual_gross_income_at_state", "entry_wealth_ratio_distribution",
                    "_linear_grid_weights_for_points", "entry_wealth_grid_weights", "_censor_entry_dead_mass"]
    seq_ast = _function_ast(seq_solver, helper_names)
    opt_ast = _function_ast(opt_solver, helper_names)
    equivalence = {name: (seq_ast[name] is not None and seq_ast[name] == opt_ast[name]) for name in helper_names}
    from intergen_eqscale_seq_optimized import solver as optimized_solver
    runtime_entry_income_gaps = []
    if z_grid is not None:
        for z in np.asarray(z_grid, dtype=float).reshape(-1):
            y_seq = reference_solver.annual_gross_income_at_state(P, 0, 0, float(z))
            y_opt = optimized_solver.annual_gross_income_at_state(P, 0, 0, float(z))
            runtime_entry_income_gaps.append(float(y_opt - y_seq))
    z_summary = _array_summary(z_grid)
    wealth_grid = _get(state, "b_grid", _get(P, "b_grid"))
    return {
        "schema": "e5f_earnings_entry_checkpoint_audit_v1",
        "checkpoint": str(checkpoint), "checkpoint_sha256": sha256(checkpoint),
        "source_root": str(source_root), "native_state_type": f"{type(state).__module__}.{type(state).__name__}",
        "entry": {
            "age": float(_get(P, "age_start", 18.0)),
            "mode": str(_get(P, "entry_wealth_mode", "scalar")),
            "active_entry_rule_metadata": _get(P, "entry_rule", _get(P, "entry_specification", None)),
            "censor_to_frontier": bool(_get(P, "entry_wealth_censor_to_frontier", False)),
            "saved_censored_mass": None if entry_censored_mass is None else float(entry_censored_mass),
            "saved_censored_share": None if _get(solution, "entry_censored_share") is None else float(_get(solution, "entry_censored_share")),
            "saved_total_entry_mass": None if entry_total is None else float(entry_total),
            "ratio_nodes": np.asarray(_get(P, "entry_wealth_ratio_nodes", []), dtype=float).tolist(),
            "ratio_weights": np.asarray(_get(P, "entry_wealth_ratio_weights", []), dtype=float).tolist(),
            "ratio_source": _get(P, "entry_wealth_ratio_source"),
            "fixed_reference_conditional_shape": None if conditional is None else list(np.shape(conditional)),
            "fixed_reference_conditional_meaning": ("present; weighted marginal is computed from this matrix" if conditional is not None else "absent"),
            "reference_marginal_source": ref_source,
            "reference_marginal_sum": None if ref_marginal is None else float(np.sum(ref_marginal)),
            "reference_marginal_mean_wealth": None if ref_marginal is None or wealth_grid is None else float(np.asarray(wealth_grid) @ ref_marginal),
            "realized_age18_prechoice_marginal_sum_before_normalization": realized_total,
            "realized_age18_prechoice_marginal_source": distribution_timing if realized_marginal is not None else None,
            "realized_age18_prechoice_mean_wealth": None if realized_marginal is None or wealth_grid is None else float(np.asarray(wealth_grid) @ realized_marginal),
            "marginal_comparison": gaps,
        },
        "saved_arrays": {
            "V": _array_summary(V), "g": _array_summary(g),
            "g_beginning_distribution": _array_summary(beginning),
            "income_grid": z_summary, "income_state_count": None if z_grid is None else int(np.asarray(z_grid).size),
            "wealth_grid": _array_summary(wealth_grid),
            "wealth_grid_min": None if wealth_grid is None else float(np.min(wealth_grid)),
            "wealth_grid_max": None if wealth_grid is None else float(np.max(wealth_grid)),
            "age18_prechoice_timing": distribution_timing if realized_marginal is not None else "unavailable",
        },
        "occupied_low_value_states": dead_occupancy,
        "all_age_occupied_low_value_states": all_age_low_values,
        "source_equivalence": {
            "seq_solver_path": str(seq_solver), "seq_solver_sha256": sha256(seq_solver),
            "optimized_solver_path": str(opt_solver), "optimized_solver_sha256": sha256(opt_solver),
            "helper_ast_equivalence_docstrings_ignored": equivalence,
            "all_entry_helpers_equivalent": bool(all(equivalence.values())),
            "annual_gross_income_runtime_entry_gap_by_income_state": runtime_entry_income_gaps,
            "annual_gross_income_runtime_entry_max_abs_gap": max(map(abs, runtime_entry_income_gaps), default=None),
            "runtime_formula_note": ("The source AST differs for income_at_state because optimized code adds property_tax_lump_sum_transfer; "
                                     "the entry conversion is equal for this checkpoint if that transfer is zero. The numeric gap is reported above."),
        },
        "audit_limits": ["Read-only checkpoint inspection; no household optimizer or equilibrium solver was called.",
                         "V<=-1e6 and V<=-1e8 masses are diagnostic summaries, not new gates.",
                         "A missing saved conditional entry matrix is reported; the unconditional marginal is used only when stored grid indices/weights exist."],
        "income_weight_error": income_weight_error,
        "source_hashes": {"checkpoint": sha256(checkpoint), "seq_solver": sha256(seq_solver),
                          "optimized_solver": sha256(opt_solver), "audit_cli": sha256(Path(__file__).resolve())},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.checkpoint, args.source_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output.resolve()), "checkpoint_sha256": result["checkpoint_sha256"],
                      "entry_censored_mass": result["entry"]["saved_censored_mass"],
                      "marginal_comparison": result["entry"]["marginal_comparison"],
                      "all_entry_helpers_equivalent": result["source_equivalence"]["all_entry_helpers_equivalent"]}, indent=2))


if __name__ == "__main__":
    main()
