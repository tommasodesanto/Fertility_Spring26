#!/usr/bin/env python3
"""Analyze saved E5F financing arrays without Bellman or equilibrium solves.

The script consumes the completed factorial summary for each checkpoint family.
It deliberately opens one arm at a time: the saved policy arrays are the
scientific object, while the frozen source is used only for the debt-cap
formula.  Missing remote-only artifacts are reported in the receipt and cause
the command to exit nonzero in strict mode.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import pickle
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BASE = ROOT / "output/model/native_financing_diagnostic_20260919/overnight/final_mechanisms"
FAMILIES = ("original", "stationary_new_income", "refit_new_income")
DOSES = (0.0, 0.25, 1.0, 5.0)
ATOL = 1.0e-12
FLOW_ATOL = 1.0e-10
PINNED_SOLVER_SHA = "2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da"
PINNED_PARAMETERS_SHA = "c0c1c18500fba069152659eaf588c3c895993cdcecb104cee7d6edca6bfae6a5"
POLICY_KEYS = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs",
               "loc_probs", "fert_probs", "fert_value", "fert2_probs", "price")
CONTRACT_KEYS = ("checkpoint", "checkpoint_sha256", "source_root", "source_manifest")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    tmp.replace(path)


def finite_float(value: Any) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def row_contract(row: Mapping[str, Any]) -> Mapping[str, Any]:
    value = row.get("contract", {})
    return value if isinstance(value, Mapping) else {}


def arm_label(row: Mapping[str, Any]) -> str:
    return str(row.get("label") or row.get("arm") or "unknown")


def infer_case_dir(row: Mapping[str, Any]) -> Path | None:
    direct = row.get("case_dir") or row.get("array_dir")
    if direct:
        return Path(str(direct))
    for key in ("arrays_path", "array_path", "cohort_arrays_path", "cohort_csv_path"):
        if row.get(key):
            return Path(str(row[key])).parent
    graphs = row.get("standard_graphs", {})
    paths = graphs.get("paths", []) if isinstance(graphs, Mapping) else []
    for raw in paths:
        path = Path(str(raw))
        for parent in path.parents:
            if parent.name.startswith("arm_"):
                return parent
    return None


def artifact_paths(row: Mapping[str, Any]) -> dict[str, Path | None]:
    case = infer_case_dir(row)
    def pick(*keys: str) -> Path | None:
        for key in keys:
            if row.get(key):
                return Path(str(row[key]))
        return None
    cohort_arrays = pick("cohort_arrays_path")
    entry = pick("initial_native_entry_path", "entry_arrays_path")
    cohort_csv = pick("cohort_csv_path", "cohort_by_age_path")
    if case is not None:
        cohort_dir = case / "cohort"
        cohort_arrays = cohort_arrays or cohort_dir / "cohort_arrays.npz"
        entry = entry or cohort_dir / "initial_native_entry.npz"
        cohort_csv = cohort_csv or cohort_dir / "cohort_by_age.csv"
    return {"case": case, "cohort_arrays": cohort_arrays,
            "entry": entry, "cohort_csv": cohort_csv}


def parse_summary_specs(values: Iterable[str], named: Mapping[str, Path | None]) -> dict[str, Path]:
    result: dict[str, Path] = {k: v for k, v in named.items() if v is not None}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"--summary requires family=PATH, got {raw!r}")
        family, path = raw.split("=", 1)
        if family not in FAMILIES:
            raise ValueError(f"unknown family {family!r}")
        result[family] = Path(path)
    missing = [family for family in FAMILIES if family not in result]
    if missing:
        raise ValueError("missing summaries: " + ", ".join(missing))
    return result


def summary_rows(summary: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    if isinstance(summary.get("cases"), list):
        return [row for row in summary["cases"] if isinstance(row, Mapping)]
    if isinstance(summary.get("rows"), list):
        return [row for row in summary["rows"] if isinstance(row, Mapping)]
    raise ValueError("summary has no cases/rows list")


def select_rows(family: str, summary_path: Path, smoke: bool) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    summary = read_json(summary_path)
    if summary.get("status") != "complete":
        raise ValueError(f"{family} summary is not completed: {summary.get('status')!r}")
    rows = summary_rows(summary)
    selected: dict[float, dict[str, Any]] = {}
    validation: dict[str, Any] = {"summary": str(summary_path), "rows_seen": len(rows), "errors": []}
    contracts: list[Mapping[str, Any]] = []
    for row0 in rows:
        row = dict(row0)
        phi = finite_float(row.get("phi"))
        cap = finite_float(row.get("rental_cap"))
        lam = finite_float(row.get("lambda"))
        if phi is None or cap is None or lam is None:
            continue
        if abs(phi - 0.8) > ATOL or abs(cap - 6.0) > ATOL or not any(abs(lam - dose) <= ATOL for dose in DOSES):
            continue
        if row.get("status") != "completed":
            validation["errors"].append(f"{arm_label(row)} status={row.get('status')!r}")
        contract = row_contract(row)
        missing_contract = [key for key in CONTRACT_KEYS if not contract.get(key)]
        if missing_contract:
            validation["errors"].append(f"{arm_label(row)} missing contract keys: {missing_contract}")
        contracts.append(contract)
        if lam not in selected:
            selected[lam] = row
    # Smoke is deliberately numerically small only in plotting/output work; it
    # still loads all four doses so identity and entry checks are exercised.
    required = DOSES
    missing = [dose for dose in required if dose not in selected]
    if missing:
        validation["errors"].append("missing doses: " + ", ".join(map(str, missing)))
    validation["selected_labels"] = {str(k): arm_label(v) for k, v in selected.items()}
    validation["selected_doses"] = sorted(selected)
    for key in ("source_root", "source_manifest", "checkpoint_sha256"):
        vals = {str(c.get(key)) for c in contracts if c.get(key) is not None}
        validation[f"{key}_values"] = sorted(vals)
        if len(vals) > 1:
            validation["errors"].append(f"within-family {key} differs")
    return [selected[dose] for dose in required if dose in selected], {"family": family, **validation, "selected": selected}


def validate_source_checkpoint(row: Mapping[str, Any]) -> dict[str, Any]:
    contract = row_contract(row)
    out: dict[str, Any] = {
        "checkpoint": contract.get("checkpoint"),
        "checkpoint_sha256": contract.get("checkpoint_sha256"),
        "source_root": contract.get("source_root"),
        "source_manifest": contract.get("source_manifest"),
        "checkpoint_status": "missing",
        "source_manifest_status": "missing",
        "source_status": "missing",
        "parameters_sha256": None,
    }
    checkpoint = Path(str(contract["checkpoint"])) if contract.get("checkpoint") else None
    if checkpoint is not None and checkpoint.exists():
        actual = sha256(checkpoint)
        out["checkpoint_actual_sha256"] = actual
        out["checkpoint_status"] = "verified" if actual == contract.get("checkpoint_sha256") else "hash_mismatch"
    elif checkpoint is not None:
        out["checkpoint_status"] = "remote_or_missing"
    manifest = Path(str(contract["source_manifest"])) if contract.get("source_manifest") else None
    if manifest is not None and manifest.exists():
        out["source_manifest_status"] = "exists"
    elif manifest is not None:
        out["source_manifest_status"] = "remote_or_missing"
    source = Path(str(contract["source_root"])) if contract.get("source_root") else None
    if source is not None and source.exists():
        solver = source / "intergen_eqscale_seq_optimized" / "solver.py"
        if not solver.exists():
            solver = source / "code" / "model" / "intergen_eqscale_seq_optimized" / "solver.py"
        if solver.exists():
            actual = sha256(solver)
            out["solver_sha256"] = actual
            parameters = source / "intergen_eqscale_seq_optimized" / "parameters.py"
            if not parameters.exists():
                parameters = source / "code" / "model" / "intergen_eqscale_seq_optimized" / "parameters.py"
            if parameters.exists():
                out["parameters_sha256"] = sha256(parameters)
            out["source_status"] = ("verified" if actual == PINNED_SOLVER_SHA and
                                     out["parameters_sha256"] == PINNED_PARAMETERS_SHA else "source_hash_mismatch")
        else:
            out["source_status"] = "solver_missing"
    elif source is not None:
        out["source_status"] = "remote_or_missing"
    return out


def load_npz(path: Path, keys: Iterable[str] | None = None) -> dict[str, np.ndarray]:
    if path is None or not path.exists():
        raise FileNotFoundError(str(path) if path else "array path unavailable")
    with np.load(path, allow_pickle=False) as archive:
        wanted = list(keys) if keys is not None else list(archive.files)
        return {key: np.asarray(archive[key]) for key in wanted if key in archive.files}


def frozen_model_root(source_root: str | None) -> Path | None:
    if not source_root:
        return None
    root = Path(source_root)
    for candidate in (root / "code" / "model", root):
        if (candidate / "intergen_eqscale_seq_optimized").exists():
            return candidate
    return None


def frozen_import_audit(model_root: Path) -> dict[str, Any]:
    bad = {}
    for name, module in sys.modules.items():
        if name.startswith("intergen_eqscale_seq_optimized") and getattr(module, "__file__", None):
            path = Path(module.__file__).resolve()
            try:
                path.relative_to(model_root.resolve())
            except ValueError:
                bad[name] = str(path)
    return {"status": "verified" if not bad else "failed", "bad_modules": bad,
            "model_root": str(model_root)}


def load_checkpoint(path: Path | None, source_root: str | None) -> tuple[np.ndarray | None, Any | None, str, dict[str, Any]]:
    if path is None or not path.exists():
        return None, None, "remote_or_missing", {}
    model_root = frozen_model_root(source_root)
    if model_root is None:
        return None, None, "frozen_source_missing", {}
    for candidate in (model_root, model_root / "tools"):
        if str(candidate) not in sys.path:
            sys.path.insert(0, str(candidate))
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as stream:
            packet = pickle.load(stream)
        grid = np.asarray(packet.get("b_grid"), dtype=float) if packet.get("b_grid") is not None else None
        audit = frozen_import_audit(model_root)
        if audit["status"] != "verified":
            return None, None, "active_module_origin", audit
        return grid, packet.get("parameters"), "loaded", audit
    except Exception as exc:  # the receipt remains useful when remote pickles need cluster packages
        return None, None, f"unloadable:{type(exc).__name__}:{exc}", {"model_root": str(model_root)}


def source_file(source_root: str | None, relative: str) -> Path | None:
    if not source_root:
        return None
    root = Path(source_root)
    for candidate in (root / relative, root / "code" / "model" / relative):
        if candidate.exists():
            return candidate
    return root / relative


def frozen_debt_caps(P: Any, lam: float, source_root: str | None) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Call the verified frozen parameters.build_debt_caps; never reimplement it."""
    module_path = source_file(source_root, "intergen_eqscale_seq_optimized/parameters.py")
    if module_path is None or not module_path.exists():
        raise FileNotFoundError("frozen parameters.py unavailable")
    module_sha = sha256(module_path)
    if module_sha != PINNED_PARAMETERS_SHA:
        raise ValueError(f"parameters.py hash mismatch: {module_sha}")
    spec = importlib.util.spec_from_file_location("frozen_e5f_parameters", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    Q = copy.deepcopy(P)
    Q.lambda_d = float(lam)
    if lam > 0:
        Q.debt_taper_start_age, Q.debt_taper_end_age = 82.0, 86.0
    module.build_debt_caps(Q)
    J = int(Q.J)
    ages = float(Q.age_start) + np.arange(J, dtype=float) * float(Q.da)
    meta = {"module_file": str(module_path), "module_sha256": module_sha,
            "build_function": "intergen_eqscale_seq_optimized.parameters.build_debt_caps",
            "taper_override": [82.0, 86.0] if lam > 0 else None}
    return np.asarray(Q.debt_caps, dtype=float), np.asarray(Q.debt_taper_weights, dtype=float), ages, meta


def weighted_mass_stats(mass: np.ndarray, grid: np.ndarray) -> tuple[float, float | None, float | None]:
    mass = np.asarray(mass, dtype=float)
    grid = np.asarray(grid, dtype=float)
    total = float(mass.sum())
    if total <= 0:
        return 0.0, None, None
    support = mass > 0.0
    return total, float(np.sum(mass * grid) / total), (float(grid[support][0]), float(grid[support][-1]))


def group_rows(g_pre: np.ndarray, grid: np.ndarray | None, caps: np.ndarray | None,
               taper: np.ndarray | None, ages: np.ndarray | None, lam: float,
               age_index: int | None = None) -> list[dict[str, Any]]:
    if grid is None or caps is None or taper is None or ages is None:
        return [{"lambda": lam, "status": "debt_support_unavailable", "reason": "checkpoint/grid/parameters unavailable"}]
    if g_pre.ndim != 7:
        raise ValueError(f"g_pre must be 7D, got {g_pre.shape}")
    nb, nt, I, J, Nz, npar, ncs = g_pre.shape
    if I != 1:
        raise ValueError(f"location axis I={I}; this diagnostic requires I=1 and does not sum locations")
    if len(grid) != nb:
        raise ValueError(f"wealth grid length {len(grid)} does not match g_pre axis {nb}")
    rows: list[dict[str, Any]] = []
    age_indices = range(J) if age_index is None else (age_index,)
    for j in age_indices:
        s_next = float(taper[min(j + 1, len(taper) - 1)])
        d_next = float(caps[min(j + 1, len(caps) - 1)])
        for z in range(Nz):
            for ten in range(nt):
                for n in range(npar):
                    for cs in range(ncs):
                        mass_b = np.asarray(g_pre[:, ten, 0, j, z, n, cs], dtype=float)
                        total, mean_b, endpoints = weighted_mass_stats(mass_b, grid)
                        if total <= 0:
                            continue
                        debt_mass = float(mass_b[grid < 0].sum())
                        endpoint_support = endpoints or (None, None)
                        row: dict[str, Any] = {
                            "lambda": lam, "age_index": j, "age_years": float(ages[j]),
                            "income_state": z, "tenure": ten, "number_children_state": n,
                            "child_state": cs, "mass": total, "debt_mass": debt_mass,
                            "debt_mass_share": debt_mass / total, "wealth_mean": mean_b,
                            "wealth_min": endpoint_support[0], "wealth_max": endpoint_support[1],
                            "grid_floor": float(grid[0]), "statutory_debt_cap_next": d_next,
                            "taper_next": s_next,
                            "mass_at_grid_floor": float(mass_b[0]),
                            "mass_at_grid_ceiling": float(mass_b[-1]),
                        }
                        rows.append(row)
    return rows


def cohort_group_rows(g_by_age: np.ndarray, grid: np.ndarray, caps: np.ndarray,
                      taper: np.ndarray, ages: np.ndarray, lam: float) -> list[dict[str, Any]]:
    """Weight support by the actual saved age-specific cohort distributions."""
    if g_by_age.ndim != 8:
        raise ValueError(f"g_pre_by_age must be 8D, got {g_by_age.shape}")
    rows: list[dict[str, Any]] = []
    for cohort_age in range(g_by_age.shape[0]):
        rows.extend(group_rows(g_by_age[cohort_age], grid, caps, taper, ages, lam, age_index=cohort_age))
    return rows


def flow_report(csv_path: Path | None, receipt: Mapping[str, Any]) -> dict[str, Any]:
    if csv_path is None or not csv_path.exists():
        return {"status": "failed", "path": str(csv_path) if csv_path else None}
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    explicit = sum(float(row["explicit_birth_flow"]) for row in rows)
    first = sum(float(row["exact_first_birth_flow"]) for row in rows)
    continuation = explicit - first
    target = receipt.get("cohort", {}) if isinstance(receipt.get("cohort"), Mapping) else receipt
    target_explicit = finite_float(target.get("cumulative_explicit_births_per_initial_household"))
    target_first = finite_float(target.get("first_births_per_initial_household"))
    gaps = [None if target_explicit is None else explicit - target_explicit,
            None if target_first is None else first - target_first]
    status = "complete" if rows and all(gap is not None and abs(gap) <= FLOW_ATOL for gap in gaps) else "failed"
    return {"status": status, "path": str(csv_path), "rows": len(rows),
            "explicit_sum": explicit, "first_sum": first, "continuation_sum": continuation,
            "explicit_gap": gaps[0], "first_gap": gaps[1],
            "age_rows": [{"age_years": float(r["age_years"]), "explicit": float(r["explicit_birth_flow"]),
                          "first": float(r["exact_first_birth_flow"]),
                          "continuation": float(r["explicit_birth_flow"]) - float(r["exact_first_birth_flow"])} for r in rows]}


def identity_report(path1: Path | None, path5: Path | None) -> dict[str, Any]:
    if path1 is None or path5 is None or not path1.exists() or not path5.exists():
        return {"status": "failed", "lambda1": str(path1) if path1 else None, "lambda5": str(path5) if path5 else None}
    with np.load(path1, allow_pickle=False) as a1, np.load(path5, allow_pickle=False) as a5:
        missing = {"lambda1": sorted(set(POLICY_KEYS) - set(a1.files)),
                   "lambda5": sorted(set(POLICY_KEYS) - set(a5.files))}
        if any(missing.values()):
            return {"status": "failed", "missing_keys": missing}
        rows = []
        for key in POLICY_KEYS:
            x, y = np.asarray(a1[key]), np.asarray(a5[key])
            if x.shape != y.shape:
                rows.append({"array": key, "shape_equal": False, "dtype_equal": False,
                             "max_abs_diff": None, "different_entries": None})
                continue
            if x.dtype != y.dtype:
                rows.append({"array": key, "shape_equal": True, "dtype_equal": False,
                             "max_abs_diff": None, "different_entries": None})
                continue
            if x.dtype.kind == "b" or y.dtype.kind == "b" or x.dtype.kind in "iu" or y.dtype.kind in "iu":
                exact = bool(np.array_equal(x, y))
                tol = exact
                different = int(np.count_nonzero(x != y))
                maxdiff = None
            else:
                exact = bool(np.array_equal(x, y, equal_nan=True))
                tol = bool(np.allclose(x, y, atol=ATOL, rtol=0.0, equal_nan=True))
                xf, yf = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
                finite_mask = np.isfinite(xf) & np.isfinite(yf)
                finite_delta = np.abs(yf[finite_mask] - xf[finite_mask])
                maxdiff = float(finite_delta.max()) if finite_delta.size else 0.0
                different = int(np.count_nonzero(~np.isclose(x, y, atol=ATOL, rtol=0.0, equal_nan=True)))
            rows.append({"array": key, "shape_equal": True, "dtype_equal": True, "exact_bitwise": exact,
                         "max_abs_diff": maxdiff, "different_entries": different,
                         "identical_atol": tol})
    structural_ok = all(r.get("shape_equal") and r.get("dtype_equal") for r in rows)
    return {"status": "complete" if structural_ok else "failed", "identical_all_atol": bool(structural_ok and all(r.get("identical_atol") for r in rows)),
            "required_keys": list(POLICY_KEYS), "rows": rows}


def entry_identity(paths: list[Path | None]) -> dict[str, Any]:
    if len(paths) != len(DOSES) or any(path is None or not path.exists() for path in paths):
        return {"status": "failed", "paths": [str(path) if path else None for path in paths]}
    existing = [path for path in paths if path is not None]
    baseline = load_npz(existing[0]).get("g_pre")
    if baseline is None:
        return {"status": "failed", "reason": "entry file lacks g_pre"}
    out = []
    for path in existing[1:]:
        arr = load_npz(path).get("g_pre")
        if arr is None or arr.shape != baseline.shape:
            out.append({"path": str(path), "identical": False, "reason": "missing_or_shape_mismatch"})
            continue
        delta = arr - baseline
        out.append({"path": str(path), "identical": bool(np.array_equal(arr, baseline)),
                    "l1": float(np.abs(delta).sum()), "linf": float(np.abs(delta).max(initial=0.0))})
    return {"status": "complete" if all(item["identical"] for item in out) else "failed",
            "baseline": str(existing[0]), "comparisons": out}


def cohort_arrays_check(path: Path | None, arrays: Mapping[str, np.ndarray] | None = None) -> dict[str, Any]:
    """Validate the one retained cohort archive and its diagonal age support."""
    if path is None or not path.exists():
        return {"status": "unavailable", "path": str(path) if path else None}
    try:
        arrays = arrays or load_npz(path, keys=("g_final", "g_pre_by_age", "lifetime_g_pre"))
        missing = sorted({"g_final", "g_pre_by_age", "lifetime_g_pre"} - set(arrays))
        g = arrays.get("g_pre_by_age")
        if missing or g is None:
            return {"status": "failed", "path": str(path), "missing_keys": missing}
        finite = bool(np.isfinite(g).all())
        nonnegative = bool(np.all(g >= -1e-12))
        mass_matrix = g.sum(axis=(1, 2, 3, 5, 6, 7)) if g.ndim == 8 else np.zeros((1, 1))
        offdiag = mass_matrix.copy(); np.fill_diagonal(offdiag, 0.0)
        diagonal = float(np.max(np.abs(offdiag)))
        return {"status": "complete" if finite and nonnegative and g.ndim == 8 and g.shape[0] == g.shape[4] and diagonal <= FLOW_ATOL else "failed",
                "path": str(path), "missing_keys": missing, "finite": finite, "nonnegative": nonnegative,
                "shape": list(g.shape), "diagonal_support_offdiag_max": diagonal,
                "shapes": {key: list(value.shape) for key, value in arrays.items()},
                "masses": {key: float(value.sum()) for key, value in arrays.items()}}
    except Exception as exc:
        return {"status": "failed", "path": str(path), "reason": f"{type(exc).__name__}: {exc}"}


def standard_plot_check(row: Mapping[str, Any]) -> dict[str, Any]:
    graphs = row.get("standard_graphs", {})
    paths = graphs.get("paths", []) if isinstance(graphs, Mapping) else []
    paths = [Path(str(path)) for path in paths]
    return {"status": "complete" if len(paths) == 17 and all(path.exists() for path in paths) else "failed",
            "count": len(paths), "paths": [str(path) for path in paths]}


def cohort_csv_mass_check(g_by_age: np.ndarray, csv_path: Path | None) -> dict[str, Any]:
    if csv_path is None or not csv_path.exists():
        return {"status": "failed", "path": str(csv_path) if csv_path else None}
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != g_by_age.shape[0] or "surviving_mass" not in rows[0]:
        return {"status": "failed", "rows": len(rows)}
    mass = g_by_age.sum(axis=(1, 2, 3, 4, 5, 6, 7))
    gaps = [float(mass[j]) - float(rows[j]["surviving_mass"]) for j in range(len(rows))]
    return {"status": "complete" if all(abs(gap) <= FLOW_ATOL for gap in gaps) else "failed",
            "max_abs_gap": max(map(abs, gaps), default=np.inf), "rows": len(rows)}


def write_supplemental_plots(case_dir: Path, group: list[Mapping[str, Any]], lam: float) -> dict[str, Any]:
    """Write only the two new debt/support views; the standard 17 are untouched."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return {"status": "unavailable", "reason": f"{type(exc).__name__}: {exc}"}
    if not group or "age_index" not in group[0]:
        return {"status": "unavailable", "reason": "debt/support rows unavailable"}
    by_age: dict[int, dict[str, float]] = {}
    for row in group:
        age = int(row["age_index"])
        rec = by_age.setdefault(age, {"mass": 0.0, "debt": 0.0, "cap": 0.0})
        rec["mass"] += float(row.get("mass", 0.0) or 0.0)
        rec["debt"] += float(row.get("debt_mass", 0.0) or 0.0)
        rec["cap"] = float(row.get("statutory_debt_cap_next", 0.0) or 0.0)
    ages = sorted(by_age)
    age_years = {int(r["age_index"]): float(r["age_years"]) for r in group}
    x = [age_years[a] for a in ages]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5))
    axes[0].plot(x, [by_age[a]["debt"] for a in ages], marker="o")
    axes[0].set(xlabel="Age", ylabel="Cohort mass", title="Liquid wealth below zero")
    axes[1].plot(x, [by_age[a]["cap"] for a in ages], marker="s")
    axes[1].set(xlabel="Age", ylabel="Model wealth units", title="Next-period statutory credit cap")
    fig.suptitle(f"Saved cohort support, lambda={lam:g}"); fig.tight_layout()
    debt_path = case_dir / "debt_support_by_age.png"
    fig.savefig(debt_path, dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.0, 3.5))
    ax.plot(x, [by_age[a]["debt"] / by_age[a]["mass"] if by_age[a]["mass"] else 0.0 for a in ages], marker="o")
    ax.set(xlabel="Age", ylabel="Share of surviving cohort", title=f"Liquid wealth below zero, lambda={lam:g}", ylim=(0.0, 1.0))
    fig.tight_layout()
    share_path = case_dir / "debt_support_share_by_age.png"
    fig.savefig(share_path, dpi=150); plt.close(fig)
    return {"status": "complete", "paths": [str(debt_path), str(share_path)]}


def analyze_family(family: str, rows: list[dict[str, Any]], info: Mapping[str, Any], output: Path,
                   started: float, wall_seconds: float, case_seconds: float) -> dict[str, Any]:
    family_dir = output / family
    family_dir.mkdir(parents=True, exist_ok=True)
    case_results: list[dict[str, Any]] = []
    entry_paths: list[Path | None] = []
    for row in rows:
        if time.monotonic() - started >= wall_seconds:
            break
        case_started = time.monotonic()
        lam = float(row["lambda"])
        paths = artifact_paths(row)
        entry_paths.append(paths["entry"])
        contract = validate_source_checkpoint(row)
        grid, P, checkpoint_status, import_audit = load_checkpoint(Path(str(contract["checkpoint"])) if contract.get("checkpoint") else None, contract.get("source_root"))
        case: dict[str, Any] = {"family": family, "lambda": lam, "label": arm_label(row),
                                "paths": {key: str(value) if value else None for key, value in paths.items()},
                                "contract": contract, "checkpoint_load": checkpoint_status,
                                "receipt_cohort": row.get("cohort", {}), "frozen_import_audit": import_audit,
                                "standard_plots": standard_plot_check(row)}
        try:
            cohort = load_npz(paths["cohort_arrays"], keys=("g_pre_by_age", "lifetime_g_pre", "g_final"))
            g_by_age = cohort.get("g_pre_by_age")
            if g_by_age is None:
                raise ValueError("cohort_arrays.npz lacks g_pre_by_age")
            if P is None or grid is None:
                raise ValueError("verified checkpoint parameters and wealth grid are required")
            caps, taper, ages, cap_meta = frozen_debt_caps(P, lam, contract.get("source_root"))
            if g_by_age.ndim != 8 or g_by_age.shape[0] != int(P.J):
                raise ValueError(f"g_pre_by_age shape {g_by_age.shape} is not (J,...), J={P.J}")
            case["status"] = "complete"
            case["cohort_array_shapes"] = {key: list(value.shape) for key, value in cohort.items()}
            case["debt_support_status"] = "complete"
            case["debt_cap_source"] = cap_meta
            case["debt_cap_by_age"] = caps.tolist()
            case["debt_taper_by_age"] = taper.tolist()
            case["cohort_group_rows"] = cohort_group_rows(g_by_age, grid, caps, taper, ages, lam)
            case["flow"] = flow_report(paths["cohort_csv"], row)
            case["cohort_arrays"] = cohort_arrays_check(paths["cohort_arrays"], cohort)
            case["cohort_csv_mass"] = cohort_csv_mass_check(g_by_age, paths["cohort_csv"])
            case_dir = family_dir / f"lambda_{lam:g}"
            write_json(case_dir / "debt_support.json", case)
            group = case["cohort_group_rows"]
            if group and "age_index" in group[0]:
                with (case_dir / "debt_support.csv").open("w", newline="", encoding="utf-8") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(group[0]))
                    writer.writeheader(); writer.writerows(group)
            case["supplemental_plots"] = write_supplemental_plots(case_dir, group, lam)
            case["elapsed_seconds"] = time.monotonic() - case_started
            if case["elapsed_seconds"] > case_seconds:
                raise TimeoutError(f"case exceeded {case_seconds:g}s")
            case.pop("cohort_group_rows", None)
            case.pop("debt_cap_by_age", None)
            case.pop("debt_taper_by_age", None)
            case_results.append(case)
            write_json(family_dir / f"lambda_{lam:g}" / "latest_completed.json", case)
            write_json(family_dir / "heartbeat.json", {"family": family, "last_lambda": lam, "status": "running"})
        except Exception as exc:
            case["status"] = "failed"
            case["error"] = f"{type(exc).__name__}: {exc}"
            case["elapsed_seconds"] = time.monotonic() - case_started
            case_results.append(case)
            write_json(family_dir / f"lambda_{lam:g}" / "latest_completed.json", case)
    identity = {"status": "unavailable", "reason": "factorial retains no per-lambda policy arrays"}
    entries = entry_identity(entry_paths)
    result = {"family": family, "status": "complete", "validation": dict(info),
              "cases": case_results, "entry_identity": entries,
              "lambda1_vs_lambda5": identity,
              "decomposition": {"status": "unavailable", "reason": "ordered saved-policy replay was not attempted; no symmetric cross-evaluation"}}
    write_json(family_dir / "summary.json", result)
    return result


def write_compact_tables(output: Path, results: Mapping[str, Any]) -> None:
    rows = []
    for family, result in results.items():
        for case in result.get("cases", []):
            flow = case.get("flow", {})
            rows.append({"family": family, "lambda": case.get("lambda"), "status": case.get("status"),
                         "explicit_sum": flow.get("explicit_sum"), "first_sum": flow.get("first_sum"),
                         "continuation_sum": flow.get("continuation_sum"), "cohort_array_path": case.get("paths", {}).get("cohort_arrays")})
    with (output / "cohort_flows.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else ["family", "lambda", "status"])
        writer.writeheader(); writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", action="append", default=[], help="family=PATH; repeat once per family")
    parser.add_argument("--summary-original", type=Path)
    parser.add_argument("--summary-pilot", type=Path)
    parser.add_argument("--summary-refit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="run the same four-dose numerical checks with compact output")
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=list(FAMILIES))
    parser.add_argument("--wall-time-seconds", type=float, default=1200.0)
    parser.add_argument("--per-case-seconds", type=float, default=180.0)
    parser.add_argument("--strict", action="store_true", help="fail on any incomplete case or invalid receipt/source/artifact gate")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    named = {"original": args.summary_original, "stationary_new_income": args.summary_pilot, "refit_new_income": args.summary_refit}
    if not args.summary and all(value is None for value in named.values()):
        named = {family: DEFAULT_BASE / ("stationary_new_income" if family == "stationary_new_income" else family) / "summary.json" for family in FAMILIES}
    summaries = parse_summary_specs(args.summary, named)
    if args.wall_time_seconds <= 0 or args.per_case_seconds <= 0:
        raise ValueError("time limits must be positive")
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"output must be fresh and empty: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    write_json(args.output / "plan.json", {"status": "planned", "summaries": {key: str(value) for key, value in summaries.items()},
                                            "families": args.families, "doses": list(DOSES),
                                            "phi": 0.8, "rental_cap": 6.0,
                                            "wall_time_seconds": args.wall_time_seconds,
                                            "per_case_seconds": args.per_case_seconds,
                                            "scope": "saved cohort arrays only; no solves"})
    results: dict[str, Any] = {}
    errors: list[str] = []
    for family in args.families:
        rows, info = select_rows(family, summaries[family], args.smoke)
        result = analyze_family(family, rows, info, args.output, started, args.wall_time_seconds, args.per_case_seconds)
        results[family] = result
        errors.extend(f"{family}: {error}" for error in info.get("errors", []))
        if args.smoke:
            errors.extend(f"{family}/{case.get('label')}: cohort smoke case incomplete"
                          for case in result.get("cases", []) if case.get("status") != "complete")
        if args.strict:
            if len(result.get("cases", [])) != len(rows):
                errors.append(f"{family}: wall-time cap left incomplete cases")
            if result.get("entry_identity", {}).get("status") != "complete":
                errors.append(f"{family}: entry identity failed")
            for case in result.get("cases", []):
                contract = case.get("contract", {})
                if case.get("status") != "complete":
                    errors.append(f"{family}/{case.get('label')}: status={case.get('status')}")
                if case.get("checkpoint_load") != "loaded" or case.get("frozen_import_audit", {}).get("status") != "verified":
                    errors.append(f"{family}/{case.get('label')}: frozen checkpoint import failed")
                for artifact in ("cohort_arrays", "cohort_csv_mass", "flow", "standard_plots", "supplemental_plots"):
                    if case.get(artifact, {}).get("status") != "complete":
                        errors.append(f"{family}/{case.get('label')}: {artifact} failed")
                for key in ("checkpoint_status", "source_manifest_status", "source_status"):
                    if contract.get(key) not in ("verified", "exists"):
                        errors.append(f"{family}/{case.get('label')}: {key}={contract.get(key)}")
    write_compact_tables(args.output, results)
    receipt = {"status": "complete" if not errors else "validation_failed", "strict": args.strict,
               "smoke": args.smoke, "summaries": {key: str(value) for key, value in summaries.items()},
               "families": results, "errors": errors,
               "elapsed_seconds": time.monotonic() - started,
               "scope": "saved cohort arrays only; no Bellman, equilibrium, or current-model fallback"}
    write_json(args.output / "receipt.json", receipt)
    if (args.strict or args.smoke) and errors:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
