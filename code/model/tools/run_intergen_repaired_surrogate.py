#!/usr/bin/env python3
"""Generate ExtraTrees proposals for the repaired 15-moment calibration.

This is a proposal generator, not a calibration evaluator. Every emitted point
must subsequently be solved by ``run_intergen_combined_recalibration.py``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import time
from pathlib import Path
from typing import Any

import numpy as np

try:
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.model_selection import GroupKFold
except ImportError as exc:  # pragma: no cover - depends on runtime environment
    raise SystemExit("scikit-learn is required for repaired-surrogate generation") from exc

from intergen_housing_fertility.calibration import get_target_set
from intergen_housing_fertility.production_profile import (
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_TARGET_SET,
)


ROOMS_NAME = "aggregate_mean_occupied_rooms_18_85"
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0
H0_BOUND = ("H0", 1.0, 20.0)
LOSS_TOL = 1e-8
RESIDUAL_GATE = 1e-4
CANDIDATE_COUNT = 250_000
LOCAL_SHARE = 0.875
LOCAL_SCALES = np.array([0.005, 0.010, 0.020, 0.040], dtype=float)
CV_OVERALL_RHO_GATE = 0.60
CV_LOCAL_RHO_GATE = 0.40
CV_LOCAL_LOSS_CUTOFF = 12.0
CV_PREDICTED_BEST_GATE = 9.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--records-root",
        type=Path,
        action="append",
        required=True,
        help="Repeat for each task tree containing task_*/cases.jsonl.",
    )
    parser.add_argument("--incumbent", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--trees", type=int, default=600)
    parser.add_argument("--proposals", type=int, default=64)
    parser.add_argument("--seed", type=int, default=2026071301)
    return parser.parse_args()


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def active_spec() -> tuple[list[tuple[str, str, float, float]], dict[str, float], dict[str, float]]:
    bounds = []
    for report_name, lower, upper in (*PRODUCTION_SEARCH_BOUNDS, H0_BOUND):
        model_name = "beta" if report_name == "beta_annual" else report_name
        bounds.append((report_name, model_name, float(lower), float(upper)))
    if len(bounds) != 14:
        raise RuntimeError(f"expected 14 search coordinates, found {len(bounds)}")

    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    targets = dict(targets)
    weights = dict(weights)
    targets[ROOMS_NAME] = ROOMS_TARGET
    weights[ROOMS_NAME] = ROOMS_WEIGHT
    if len(targets) != 15 or set(targets) != set(weights):
        raise RuntimeError("expected the repaired 15-moment target/weight system")
    return bounds, targets, weights


def theta_to_unit(theta: dict[str, Any], bounds: list[tuple[str, str, float, float]]) -> np.ndarray:
    values = []
    for report_name, model_name, lower, upper in bounds:
        if not finite(theta.get(model_name)):
            raise ValueError(f"missing/nonfinite theta coordinate {model_name}")
        value = float(theta[model_name])
        if report_name == "beta_annual":
            if value <= 0.0:
                raise ValueError("beta must be positive")
            value = value ** 0.25
        unit = (value - lower) / (upper - lower)
        if unit < -1e-10 or unit > 1.0 + 1e-10:
            raise ValueError(f"theta coordinate {report_name} is outside its active bound")
        values.append(float(np.clip(unit, 0.0, 1.0)))
    return np.asarray(values, dtype=float)


def unit_to_theta(unit: np.ndarray, bounds: list[tuple[str, str, float, float]]) -> dict[str, float]:
    theta: dict[str, float] = {}
    for u, (report_name, model_name, lower, upper) in zip(unit, bounds):
        value = lower + float(np.clip(u, 0.0, 1.0)) * (upper - lower)
        theta[model_name] = value**4 if report_name == "beta_annual" else value
    return theta


def recomputed_loss(moments: dict[str, Any], targets: dict[str, float], weights: dict[str, float]) -> float:
    total = 0.0
    for name, target in targets.items():
        if not finite(moments.get(name)):
            raise ValueError(f"missing/nonfinite active moment {name}")
        gap = float(moments[name]) - float(target)
        total += float(weights[name]) * gap * gap
    return float(total)


def validate_record(
    record: dict[str, Any],
    bounds: list[tuple[str, str, float, float]],
    targets: dict[str, float],
    weights: dict[str, float],
) -> tuple[np.ndarray, float]:
    if record.get("status") != "ok" or record.get("strict_converged") is not True:
        raise ValueError("record is not status-ok and strict-converged")
    residual = record.get("market_residual", dict(record.get("moments", {})).get("market_residual"))
    if not finite(residual) or float(residual) > RESIDUAL_GATE:
        raise ValueError("record fails the strict market-residual gate")
    if not finite(record.get("rank_loss")):
        raise ValueError("record has no finite rank_loss")
    theta = record.get("theta")
    moments = record.get("moments")
    if not isinstance(theta, dict) or not isinstance(moments, dict):
        raise ValueError("record lacks theta or moments object")
    unit = theta_to_unit(theta, bounds)
    loss = recomputed_loss(moments, targets, weights)
    reported = float(record["rank_loss"])
    if abs(reported - loss) > LOSS_TOL:
        raise ValueError(f"stored/recomputed loss mismatch: {reported} versus {loss}")
    return unit, loss


def read_candidate_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"JSON object required: {path}")
    candidate = payload.get("best", payload)
    if not isinstance(candidate, dict):
        raise ValueError(f"candidate object required: {path}")
    return candidate


def load_dataset(
    roots: list[Path],
    bounds: list[tuple[str, str, float, float]],
    targets: dict[str, float],
    weights: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    unique: dict[tuple[float, ...], tuple[np.ndarray, float, str]] = {}
    counts = {
        "files": 0,
        "records_seen": 0,
        "records_noneligible": 0,
        "eligible_contract_records": 0,
        "duplicate_theta_records": 0,
    }
    root_rows = []
    for root_index, root in enumerate(roots):
        root = root.resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"record root is not a directory: {root}")
        files = sorted(root.rglob("cases.jsonl"))
        if not files:
            raise FileNotFoundError(f"no cases.jsonl beneath record root: {root}")
        root_count = 0
        for path in files:
            counts["files"] += 1
            try:
                relative_group = str(path.parent.relative_to(root))
            except ValueError:
                relative_group = str(path.parent)
            group = f"root{root_index}:{relative_group}"
            with path.open() as handle:
                for line_number, line in enumerate(handle, start=1):
                    if not line.strip():
                        continue
                    counts["records_seen"] += 1
                    try:
                        record = json.loads(line)
                    except json.JSONDecodeError as exc:
                        raise ValueError(f"malformed JSON at {path}:{line_number}") from exc
                    if not isinstance(record, dict):
                        raise ValueError(f"non-object record at {path}:{line_number}")
                    if record.get("status") != "ok" or record.get("strict_converged") is not True:
                        counts["records_noneligible"] += 1
                        continue
                    try:
                        unit, loss = validate_record(record, bounds, targets, weights)
                    except ValueError as exc:
                        raise ValueError(f"eligible-record contract failure at {path}:{line_number}: {exc}") from exc
                    counts["eligible_contract_records"] += 1
                    root_count += 1
                    key = tuple(np.round(unit, 12))
                    previous = unique.get(key)
                    if previous is not None:
                        counts["duplicate_theta_records"] += 1
                    if previous is None or loss < previous[1]:
                        unique[key] = (unit, loss, group)
        root_rows.append({"root": str(root), "files": len(files), "eligible_rows": root_count})

    if not unique:
        raise RuntimeError("no records satisfy the repaired-objective data contract")
    values = list(unique.values())
    X = np.vstack([row[0] for row in values])
    y = np.asarray([row[1] for row in values], dtype=float)
    groups = np.asarray([row[2] for row in values], dtype=object)
    counts["unique_theta"] = int(len(values))
    counts["group_count"] = int(len(set(groups.tolist())))
    counts["roots"] = root_rows
    return X, y, groups, counts


def rank_correlation(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if actual.size < 3:
        return math.nan
    ra = np.argsort(np.argsort(actual, kind="mergesort"), kind="mergesort").astype(float)
    rp = np.argsort(np.argsort(predicted, kind="mergesort"), kind="mergesort").astype(float)
    ra -= ra.mean()
    rp -= rp.mean()
    denom = math.sqrt(float(ra @ ra) * float(rp @ rp))
    return float((ra @ rp) / denom) if denom > 0.0 else math.nan


def make_model(trees: int, seed: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=max(10, int(trees)),
        min_samples_leaf=2,
        max_features=0.8,
        n_jobs=1,
        random_state=int(seed),
    )


def grouped_cv(X: np.ndarray, y: np.ndarray, groups: np.ndarray, trees: int, seed: int) -> dict[str, Any]:
    n_groups = len(set(groups.tolist()))
    folds = min(4, n_groups)
    if folds < 3:
        raise RuntimeError("at least three task groups are required for grouped CV")
    fold_rows = []
    passed = True
    for fold, (train, test) in enumerate(GroupKFold(folds).split(X, y, groups)):
        model = make_model(min(int(trees), 300), seed + fold)
        model.fit(X[train], y[train])
        predicted = model.predict(X[test])
        overall_rho = rank_correlation(y[test], predicted)
        local = y[test] < CV_LOCAL_LOSS_CUTOFF
        local_rho = rank_correlation(y[test][local], predicted[local]) if int(local.sum()) >= 3 else math.nan
        top_n = max(1, int(math.ceil(0.01 * len(test))))
        predicted_top = np.argsort(predicted)[:top_n]
        predicted_top_actual_min = float(np.min(y[test][predicted_top]))
        fold_pass = bool(
            finite(overall_rho)
            and overall_rho >= CV_OVERALL_RHO_GATE
            and finite(local_rho)
            and local_rho >= CV_LOCAL_RHO_GATE
            and predicted_top_actual_min < CV_PREDICTED_BEST_GATE
        )
        passed = passed and fold_pass
        fold_rows.append(
            {
                "fold": fold,
                "train_rows": int(len(train)),
                "test_rows": int(len(test)),
                "test_local_rows_loss_lt_12": int(local.sum()),
                "overall_spearman": overall_rho,
                "local_spearman_loss_lt_12": local_rho,
                "median_absolute_error": float(np.median(np.abs(y[test] - predicted))),
                "predicted_best_1pct_actual_min": predicted_top_actual_min,
                "actual_test_min": float(np.min(y[test])),
                "pass": fold_pass,
            }
        )
    return {
        "passed": passed,
        "fold_count": folds,
        "gates": {
            "overall_spearman_min": CV_OVERALL_RHO_GATE,
            "local_spearman_loss_lt_12_min": CV_LOCAL_RHO_GATE,
            "predicted_best_1pct_actual_loss_below": CV_PREDICTED_BEST_GATE,
        },
        "folds": fold_rows,
    }


def global_probe(n: int, dim: int, seed: int) -> np.ndarray:
    try:
        from scipy.stats import qmc

        power = int(math.ceil(math.log2(max(1, n))))
        return qmc.Sobol(dim, scramble=True, seed=seed).random_base2(power)[:n]
    except ImportError:
        return np.random.default_rng(seed).random((n, dim))


def generate_candidates(X: np.ndarray, y: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    local_n = int(round(CANDIDATE_COUNT * LOCAL_SHARE))
    global_n = CANDIDATE_COUNT - local_n
    center_n = min(200, len(y))
    centers = X[np.argsort(y)[:center_n]]
    center_index = rng.integers(0, center_n, size=local_n)
    scale_index = rng.integers(0, len(LOCAL_SCALES), size=local_n)
    local = centers[center_index] + rng.normal(size=(local_n, X.shape[1])) * LOCAL_SCALES[scale_index, None]
    local = np.clip(local, 0.0, 1.0)
    global_points = global_probe(global_n, X.shape[1], seed + 17)
    candidates = np.vstack((local, global_points))
    source = np.concatenate((np.full(local_n, "local", dtype=object), np.full(global_n, "global", dtype=object)))
    scale = np.concatenate((LOCAL_SCALES[scale_index], np.full(global_n, np.nan)))
    return candidates, source, scale


def ensemble_prediction(model: ExtraTreesRegressor, X: np.ndarray, chunk: int = 5000) -> tuple[np.ndarray, np.ndarray]:
    mean = np.empty(len(X), dtype=float)
    std = np.empty(len(X), dtype=float)
    for start in range(0, len(X), chunk):
        end = min(start + chunk, len(X))
        draws = np.vstack([tree.predict(X[start:end]) for tree in model.estimators_])
        mean[start:end] = draws.mean(axis=0)
        std[start:end] = draws.std(axis=0)
    return mean, std


def select_diverse(
    candidates: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    count: int,
    incumbent_unit: np.ndarray,
) -> tuple[list[int], list[str]]:
    modes = {
        "mean": mean,
        "lcb_0p5": mean - 0.5 * std,
        "lcb_1": mean - std,
        "lcb_2": mean - 2.0 * std,
    }
    pools = {name: np.argsort(score)[: max(2000, count * 50)] for name, score in modes.items()}
    cursors = {name: 0 for name in modes}
    selected: list[int] = []
    labels: list[str] = []
    names = list(modes)
    attempts = 0
    while len(selected) < count and attempts < sum(len(x) for x in pools.values()):
        name = names[attempts % len(names)]
        pool = pools[name]
        cursor = cursors[name]
        cursors[name] += 1
        attempts += 1
        if cursor >= len(pool):
            continue
        idx = int(pool[cursor])
        point = candidates[idx]
        if np.linalg.norm(point - incumbent_unit) < 1e-6:
            continue
        if any(np.linalg.norm(point - candidates[old]) < 0.02 for old in selected):
            continue
        selected.append(idx)
        labels.append(name)
    if len(selected) < count:
        for idx in np.argsort(mean):
            idx = int(idx)
            if idx in selected or np.linalg.norm(candidates[idx] - incumbent_unit) < 1e-6:
                continue
            selected.append(idx)
            labels.append("mean_topup")
            if len(selected) == count:
                break
    if len(selected) != count:
        raise RuntimeError(f"could select only {len(selected)} of {count} proposals")
    return selected, labels


def nearest_training_distance(point: np.ndarray, X: np.ndarray) -> float:
    return float(np.sqrt(np.min(np.sum((X - point) ** 2, axis=1))))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def json_dump(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main() -> None:
    args = parse_args()
    if args.trees < 10 or args.proposals < 1:
        raise ValueError("trees must be at least 10 and proposals positive")
    started = time.time()
    bounds, targets, weights = active_spec()
    X, y, groups, counts = load_dataset(args.records_root, bounds, targets, weights)

    incumbent = read_candidate_json(args.incumbent.resolve())
    incumbent_unit, incumbent_loss = validate_record(incumbent, bounds, targets, weights)
    if incumbent_loss > 8.925945490698217 + 1e-6:
        raise ValueError(
            f"incumbent loss {incumbent_loss} is worse than the repaired final 8.9259454907"
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    cv = grouped_cv(X, y, groups, args.trees, args.seed)
    json_dump(args.outdir / "cv.json", cv)
    if not cv["passed"]:
        raise RuntimeError("task-grouped surrogate CV gates failed; no proposals emitted")

    model = make_model(args.trees, args.seed)
    model.fit(X, y)
    with (args.outdir / "surrogate.pkl").open("wb") as handle:
        pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    candidate_X, candidate_source, candidate_scale = generate_candidates(X, y, args.seed + 101)
    predicted_mean, predicted_std = ensemble_prediction(model, candidate_X)
    selected, selection_modes = select_diverse(
        candidate_X, predicted_mean, predicted_std, args.proposals, incumbent_unit
    )
    selected_rank = {idx: rank for rank, idx in enumerate(selected, start=1)}

    proposal_rows = []
    proposal_root = args.outdir / "proposals"
    for rank, (idx, mode) in enumerate(zip(selected, selection_modes), start=1):
        point = candidate_X[idx]
        theta = unit_to_theta(point, bounds)
        row = {
            "proposal_rank": rank,
            "candidate_index": idx,
            "selection_mode": mode,
            "candidate_source": str(candidate_source[idx]),
            "local_scale": None if not finite(candidate_scale[idx]) else float(candidate_scale[idx]),
            "predicted_mean_loss": float(predicted_mean[idx]),
            "predicted_tree_std": float(predicted_std[idx]),
            "nearest_training_distance": nearest_training_distance(point, X),
            "theta": theta,
            "unit_vector": point.tolist(),
        }
        proposal_rows.append(row)
        task_dir = proposal_root / f"task_{rank}"
        task_dir.mkdir(parents=True, exist_ok=True)
        json_dump(task_dir / "best.json", row)

    score_candidates: set[int] = set(selected)
    acquisitions = [predicted_mean, predicted_mean - 0.5 * predicted_std, predicted_mean - predicted_std, predicted_mean - 2.0 * predicted_std]
    for acquisition in acquisitions:
        score_candidates.update(int(x) for x in np.argsort(acquisition)[:1250])
    score_rows = []
    for idx in sorted(score_candidates, key=lambda i: float(predicted_mean[i])):
        score_rows.append(
            {
                "candidate_index": idx,
                "candidate_source": candidate_source[idx],
                "local_scale": "" if not finite(candidate_scale[idx]) else float(candidate_scale[idx]),
                "predicted_mean_loss": float(predicted_mean[idx]),
                "predicted_tree_std": float(predicted_std[idx]),
                "lcb_0p5": float(predicted_mean[idx] - 0.5 * predicted_std[idx]),
                "lcb_1": float(predicted_mean[idx] - predicted_std[idx]),
                "lcb_2": float(predicted_mean[idx] - 2.0 * predicted_std[idx]),
                "selected_rank": selected_rank.get(idx, ""),
            }
        )
    write_csv(args.outdir / "candidate_scores.csv", score_rows)
    write_csv(
        args.outdir / "proposals.csv",
        [
            {
                **{key: value for key, value in row.items() if key not in {"theta", "unit_vector"}},
                **row["theta"],
            }
            for row in proposal_rows
        ],
    )

    manifest = {
        "status": "proposals_only_real_solver_validation_required",
        "algorithm": "extra_trees_trust_region_lcb",
        "seed": int(args.seed),
        "trees": int(args.trees),
        "proposal_count": int(args.proposals),
        "candidate_count": CANDIDATE_COUNT,
        "local_candidate_share": LOCAL_SHARE,
        "local_scales_normalized": LOCAL_SCALES.tolist(),
        "record_roots": [str(path.resolve()) for path in args.records_root],
        "incumbent_path": str(args.incumbent.resolve()),
        "incumbent_loss": incumbent_loss,
        "target_set": PRODUCTION_TARGET_SET,
        "targets": targets,
        "weights": weights,
        "bounds": [
            {"name": report_name, "model_key": model_name, "lower": lower, "upper": upper}
            for report_name, model_name, lower, upper in bounds
        ],
        "data_contract": {
            "status": "ok",
            "strict_converged": True,
            "market_residual_max": RESIDUAL_GATE,
            "stored_recomputed_loss_tolerance": LOSS_TOL,
            "parameter_count": 14,
            "moment_count": 15,
        },
        "dataset_counts": counts,
        "dataset_loss_min": float(np.min(y)),
        "dataset_loss_median": float(np.median(y)),
        "cv_passed": bool(cv["passed"]),
        "elapsed_sec": time.time() - started,
    }
    json_dump(args.outdir / "dataset_manifest.json", manifest)
    json_dump(args.outdir / "proposals.json", proposal_rows)
    print(
        f"wrote {args.proposals} proposals from {len(y)} unique repaired records; "
        f"CV passed; incumbent={incumbent_loss:.10f}; elapsed={time.time()-started:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
