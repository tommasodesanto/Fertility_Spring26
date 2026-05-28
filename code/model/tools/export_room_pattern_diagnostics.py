"""Export targeted room-pattern diagnostics for direct-calibration records.

The goal is to diagnose the housing-size failure without changing the SMM
objective.  The script solves one or more saved ``best.json`` records and writes
CSV tables for realized renter rooms, owner rung use, ownership profiles, and
the support behind the housing event moments.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np

MODEL_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.objective import extract_moments  # noqa: E402
from dt_cp_model.parameters import asdict  # noqa: E402
from dt_cp_model.solver import make_grid, run_model_cp_dt  # noqa: E402
from dt_cp_model.theta import apply_theta  # noqa: E402


DEFAULT_OUTDIR = REPO_ROOT / "output" / "model" / "room_pattern_diagnostics_20260527"
ROOM_BINS = [
    ("le4", lambda x: x <= 4.0),
    ("5", lambda x: (x > 4.0) & (x <= 5.0)),
    ("6", lambda x: (x > 5.0) & (x <= 6.0)),
    ("7_8", lambda x: (x > 6.0) & (x <= 8.0)),
    ("9_10", lambda x: (x > 8.0) & (x <= 10.0)),
    ("11plus", lambda x: x > 10.0),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Export room-pattern diagnostics for calibration records.")
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help="Case spec LABEL=path/to/best.json,hR_max[,H1;H2;...]",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_case_rows: list[dict] = []
    all_room_rows: list[dict] = []
    all_cap_rows: list[dict] = []
    all_rung_rows: list[dict] = []
    all_own_rows: list[dict] = []
    all_age_rows: list[dict] = []
    all_event_rows: list[dict] = []

    for spec in args.case:
        label, record_path, hR_max, H_own = parse_case(spec)
        print(f"Solving {label}: {record_path} with hR_max={hR_max:g}", flush=True)
        record, sol, P, p_eq, b_grid = solve_record(record_path, hR_max, H_own)
        fresh = extract_moments(sol, P, p_eq, 0.0, 0.0, 0.0, True)

        case_out = args.outdir / label
        case_out.mkdir(parents=True, exist_ok=True)

        case_rows = case_summary_rows(label, record_path, record, fresh, P, p_eq)
        room_rows = room_distribution_rows(label, sol, P)
        cap_rows = renter_cap_rows(label, sol, P)
        rung_rows = owner_rung_rows(label, sol, P)
        own_rows = ownership_profile_rows(label, sol, P)
        age_rows = ownership_by_age_rows(label, sol, P)
        event_rows = event_context_rows(label, record, fresh, sol, P)

        write_csv(case_out / "case_summary.csv", case_rows)
        write_csv(case_out / "room_distribution.csv", room_rows)
        write_csv(case_out / "renter_cap_shares.csv", cap_rows)
        write_csv(case_out / "owner_rung_shares.csv", rung_rows)
        write_csv(case_out / "ownership_profiles.csv", own_rows)
        write_csv(case_out / "ownership_by_age.csv", age_rows)
        write_csv(case_out / "event_context.csv", event_rows)

        all_case_rows.extend(case_rows)
        all_room_rows.extend(room_rows)
        all_cap_rows.extend(cap_rows)
        all_rung_rows.extend(rung_rows)
        all_own_rows.extend(own_rows)
        all_age_rows.extend(age_rows)
        all_event_rows.extend(event_rows)

    write_csv(args.outdir / "case_summary.csv", all_case_rows)
    write_csv(args.outdir / "room_distribution.csv", all_room_rows)
    write_csv(args.outdir / "renter_cap_shares.csv", all_cap_rows)
    write_csv(args.outdir / "owner_rung_shares.csv", all_rung_rows)
    write_csv(args.outdir / "ownership_profiles.csv", all_own_rows)
    write_csv(args.outdir / "ownership_by_age.csv", all_age_rows)
    write_csv(args.outdir / "event_context.csv", all_event_rows)
    print(f"Wrote diagnostics to {args.outdir}", flush=True)


def parse_case(spec: str) -> tuple[str, Path, float, list[float] | None]:
    if "=" not in spec:
        raise ValueError(f"Case must be LABEL=path,hR_max[,H1;H2;...], got {spec!r}")
    label, rhs = spec.split("=", 1)
    parts = rhs.split(",", 2)
    if len(parts) not in (2, 3):
        raise ValueError(f"Case must be LABEL=path,hR_max[,H1;H2;...], got {spec!r}")
    H_own = None
    if len(parts) == 3 and parts[2].strip():
        H_own = [float(x) for x in parts[2].split(";") if x.strip()]
    return label, Path(parts[0]), float(parts[1]), H_own


def solve_record(
    record_path: Path,
    hR_max: float,
    H_own: list[float] | None,
) -> tuple[dict, SimpleNamespace, SimpleNamespace, np.ndarray, np.ndarray]:
    record = json.loads(record_path.read_text())
    owner_h_bar_scale = record.get("owner_h_bar_scale")
    owner_size_cost = record.get("owner_size_cost")
    owner_size_cost_ref = record.get("owner_size_cost_ref")
    owner_size_cost_power = record.get("owner_size_cost_power")
    tenure_choice_kappa = record.get("tenure_choice_kappa")
    alpha_cons = record.get("alpha_cons")
    setup = build_direct_calibration_setup(
        "benchmark",
        geo_weight=100.0,
        population_closure="outside_option_benchmark_normalized",
        scale_target=1.0,
        scale_weight=100.0,
        hR_max=hR_max,
        alpha_cons=alpha_cons,
        owner_h_bar_scale=owner_h_bar_scale,
        owner_size_cost=owner_size_cost,
        owner_size_cost_ref=owner_size_cost_ref,
        owner_size_cost_power=owner_size_cost_power,
        tenure_choice_kappa=tenure_choice_kappa,
        H_own=H_own,
    )
    theta = np.asarray(record["theta"], dtype=float)
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    P = SimpleNamespace(**copy.deepcopy(asdict(setup.P_base)))
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    P = apply_theta(P, [theta_dict[name] for name in structural_names], structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    sol, P, p_eq = run_model_cp_dt(P, verbose=False)
    return record, sol, P, p_eq, make_grid(P)


def case_summary_rows(label: str, record_path: Path, record: dict, fresh, P, p_eq: np.ndarray) -> list[dict]:
    rows = []
    rows.extend(
        [
            metric_row(label, "record_path", str(record_path)),
            metric_row(label, "run_tag", record.get("run_tag", "")),
            metric_row(label, "job_id", record.get("job_id", "")),
            metric_row(label, "eval_id", record.get("eval_id", "")),
            metric_row(label, "record_loss", record.get("loss", math.nan)),
            metric_row(label, "record_solve_ok", record.get("solve_ok", "")),
            metric_row(label, "record_converged", record.get("converged", "")),
            metric_row(label, "hR_max", P.hR_max),
            metric_row(label, "H_own", " ".join(f"{x:.6g}" for x in np.asarray(P.H_own))),
            metric_row(label, "p_P", float(p_eq[0])),
            metric_row(label, "p_C", float(p_eq[1])),
        ]
    )
    for key in [
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "own_rate",
        "own_gradient",
        "own_family_gap",
        "own_lifecycle_slope",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_0to1_onechild",
        "housing_increment_0to2plus",
        "housing_increment_1to2",
        "parity_progression_1to2",
    ]:
        rows.append(metric_row(label, key, getattr(fresh, key, math.nan)))
    return rows


def metric_row(case: str, metric: str, value) -> dict:
    return {"case": case, "metric": metric, "value": value}


def room_distribution_rows(case: str, sol, P) -> list[dict]:
    rows: list[dict] = []
    for age_label, ages in age_windows(P).items():
        for tenure in ["renter", "owner", "all"]:
            for concept in ["current_children", "completed_fertility"]:
                for child_bin in ["0", "1", "2plus", "all"]:
                    cells = collect_selected_cells(sol, P, ages=ages, tenure=tenure, location="all", concept=concept, child_bin=child_bin)
                    if cells["weight"].size == 0 or float(np.sum(cells["weight"])) <= 0:
                        continue
                    rows.append(room_distribution_row(case, age_label, tenure, concept, child_bin, cells, P))
    return rows


def room_distribution_row(case: str, age_label: str, tenure: str, concept: str, child_bin: str, cells: dict, P) -> dict:
    w = cells["weight"]
    rooms = cells["rooms"]
    total = float(np.sum(w))
    row = {
        "case": case,
        "age_window": age_label,
        "tenure": tenure,
        "child_concept": concept,
        "child_bin": child_bin,
        "mass": total,
        "mean_rooms": weighted_mean(rooms, w),
        "p10_rooms": weighted_quantile(rooms, w, 0.10),
        "p25_rooms": weighted_quantile(rooms, w, 0.25),
        "p50_rooms": weighted_quantile(rooms, w, 0.50),
        "p75_rooms": weighted_quantile(rooms, w, 0.75),
        "p90_rooms": weighted_quantile(rooms, w, 0.90),
    }
    for name, mask_fn in ROOM_BINS:
        row[f"share_rooms_{name}"] = weighted_share(mask_fn(rooms), w)
    if tenure == "renter":
        row["share_at_renter_cap"] = weighted_share(rooms >= float(P.hR_max) - 1e-8, w)
    else:
        row["share_at_renter_cap"] = math.nan
    return row


def renter_cap_rows(case: str, sol, P) -> list[dict]:
    rows: list[dict] = []
    for age_label, ages in age_windows(P).items():
        for concept in ["current_children", "completed_fertility"]:
            for child_bin in ["0", "1", "2plus", "all"]:
                cells = collect_selected_cells(sol, P, ages=ages, tenure="renter", location="all", concept=concept, child_bin=child_bin)
                w = cells["weight"]
                if w.size == 0 or float(np.sum(w)) <= 0:
                    continue
                rooms = cells["rooms"]
                rows.append(
                    {
                        "case": case,
                        "age_window": age_label,
                        "child_concept": concept,
                        "child_bin": child_bin,
                        "hR_max": P.hR_max,
                        "renter_mass": float(np.sum(w)),
                        "share_at_cap": weighted_share(rooms >= float(P.hR_max) - 1e-8, w),
                        "share_le4": weighted_share(rooms <= 4.0, w),
                        "mean_rooms": weighted_mean(rooms, w),
                        "p50_rooms": weighted_quantile(rooms, w, 0.50),
                    }
                )
    return rows


def owner_rung_rows(case: str, sol, P) -> list[dict]:
    rows: list[dict] = []
    for age_label, ages in age_windows(P).items():
        for concept in ["current_children", "completed_fertility"]:
            for child_bin in ["0", "1", "2plus", "all"]:
                cells = collect_selected_cells(sol, P, ages=ages, tenure="owner", location="all", concept=concept, child_bin=child_bin)
                total = float(np.sum(cells["weight"]))
                if total <= 0:
                    continue
                for rung, rooms in enumerate(P.H_own, start=1):
                    mass = owner_rung_mass(sol, P, ages, concept, child_bin, rung)
                    rows.append(
                        {
                            "case": case,
                            "age_window": age_label,
                            "child_concept": concept,
                            "child_bin": child_bin,
                            "owner_rung": rung,
                            "rooms": float(rooms),
                            "owner_mass_total": total,
                            "rung_mass": mass,
                            "rung_share": mass / max(total, 1e-14),
                        }
                    )
    return rows


def ownership_profile_rows(case: str, sol, P) -> list[dict]:
    rows: list[dict] = []
    for age_label, ages in age_windows(P).items():
        for location in ["all", "periphery", "center"]:
            for concept in ["current_children", "completed_fertility"]:
                for child_bin in ["0", "1", "2plus", "all"]:
                    cells = collect_selected_cells(sol, P, ages=ages, tenure="all", location=location, concept=concept, child_bin=child_bin)
                    w = cells["weight"]
                    if w.size == 0 or float(np.sum(w)) <= 0:
                        continue
                    owner = cells["owner"]
                    rooms = cells["rooms"]
                    row = {
                        "case": case,
                        "age_window": age_label,
                        "location": location,
                        "child_concept": concept,
                        "child_bin": child_bin,
                        "mass": float(np.sum(w)),
                        "own_rate": weighted_mean(owner, w),
                        "mean_rooms": weighted_mean(rooms, w),
                        "p50_rooms": weighted_quantile(rooms, w, 0.50),
                    }
                    for tenure in ["renter", "owner"]:
                        tcells = collect_selected_cells(sol, P, ages=ages, tenure=tenure, location=location, concept=concept, child_bin=child_bin)
                        tw = tcells["weight"]
                        prefix = "renter" if tenure == "renter" else "owner"
                        if tw.size == 0 or float(np.sum(tw)) <= 0:
                            row[f"{prefix}_mass"] = 0.0
                            row[f"{prefix}_mean_rooms"] = math.nan
                            row[f"{prefix}_p50_rooms"] = math.nan
                        else:
                            row[f"{prefix}_mass"] = float(np.sum(tw))
                            row[f"{prefix}_mean_rooms"] = weighted_mean(tcells["rooms"], tw)
                            row[f"{prefix}_p50_rooms"] = weighted_quantile(tcells["rooms"], tw, 0.50)
                    rows.append(row)
    return rows


def ownership_by_age_rows(case: str, sol, P) -> list[dict]:
    rows: list[dict] = []
    for age in range(P.age_start, P.age_start + P.J):
        ages = [age_index(P, age)]
        for location in ["all", "periphery", "center"]:
            for child_bin in ["0", "1", "2plus", "all"]:
                cells = collect_selected_cells(
                    sol,
                    P,
                    ages=ages,
                    tenure="all",
                    location=location,
                    concept="current_children",
                    child_bin=child_bin,
                )
                w = cells["weight"]
                if w.size == 0 or float(np.sum(w)) <= 0:
                    continue
                rows.append(
                    {
                        "case": case,
                        "age": age,
                        "location": location,
                        "current_child_bin": child_bin,
                        "mass": float(np.sum(w)),
                        "own_rate": weighted_mean(cells["owner"], w),
                        "mean_rooms": weighted_mean(cells["rooms"], w),
                    }
                )
    return rows


def event_context_rows(case: str, record: dict, fresh, sol, P) -> list[dict]:
    rows = []
    for metric in [
        "housing_increment_0to1",
        "housing_increment_0to1_onechild",
        "housing_increment_0to2plus",
        "housing_increment_1to2",
        "parity_progression_1to2",
    ]:
        rows.append(event_row(case, "fresh_moment", metric, getattr(fresh, metric, math.nan)))
        if metric in record.get("moments", {}):
            rows.append(event_row(case, "record_moment", metric, record["moments"][metric]))
    rows.append(event_row(case, "fresh_moment", "housing_event_horizon", getattr(sol, "housing_event_horizon", math.nan)))

    for age_label in ["25_45", "fertility_complete", "65_75"]:
        ages = age_windows(P)[age_label]
        for concept in ["current_children", "completed_fertility"]:
            for child_bin in ["0", "1", "2plus", "all"]:
                cells = collect_selected_cells(sol, P, ages=ages, tenure="all", location="all", concept=concept, child_bin=child_bin)
                w = cells["weight"]
                if w.size == 0 or float(np.sum(w)) <= 0:
                    continue
                stem = f"{age_label}_{concept}_{child_bin}"
                rows.append(event_row(case, "support_mass", f"{stem}_mass", float(np.sum(w))))
                rows.append(event_row(case, "support_mean", f"{stem}_mean_rooms", weighted_mean(cells["rooms"], w)))
                rows.append(event_row(case, "support_mean", f"{stem}_own_rate", weighted_mean(cells["owner"], w)))
    return rows


def event_row(case: str, source: str, metric: str, value) -> dict:
    return {"case": case, "source": source, "metric": metric, "value": value}


def collect_selected_cells(sol, P, *, ages: Iterable[int], tenure: str, location: str, concept: str, child_bin: str) -> dict:
    cols = {k: [] for k in ["weight", "rooms", "owner"]}
    locs = location_indices(P, location)
    for j in ages:
        if j < 0 or j >= P.J:
            continue
        for loc in locs:
            for ten in range(1 + P.n_house):
                if not tenure_match(ten, tenure):
                    continue
                for nn in range(P.n_parity):
                    for cs in range(P.n_child_states):
                        if not child_match(P, nn, cs, concept, child_bin):
                            continue
                        mass = sol.g[:, ten, loc, j, nn, cs]
                        nz = mass > 0
                        if not np.any(nz):
                            continue
                        w = mass[nz]
                        if ten == 0:
                            rooms = sol.hR_pol[nz, ten, loc, j, nn, cs]
                            owner = np.zeros(np.count_nonzero(nz))
                        else:
                            rooms = np.full(np.count_nonzero(nz), float(P.H_own[ten - 1]))
                            owner = np.ones(np.count_nonzero(nz))
                        cols["weight"].append(w)
                        cols["rooms"].append(np.asarray(rooms, dtype=float))
                        cols["owner"].append(owner)
    return {name: concat(parts) for name, parts in cols.items()}


def owner_rung_mass(sol, P, ages: Iterable[int], concept: str, child_bin: str, rung: int) -> float:
    total = 0.0
    for j in ages:
        for loc in range(P.I):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    if child_match(P, nn, cs, concept, child_bin):
                        total += float(np.sum(sol.g[:, rung, loc, j, nn, cs]))
    return total


def age_windows(P) -> dict[str, range]:
    return {
        "25_34": range(age_index(P, 25), age_index(P, 34) + 1),
        "35_44": range(age_index(P, 35), age_index(P, 44) + 1),
        "25_45": range(age_index(P, 25), age_index(P, 45) + 1),
        "45_55": range(age_index(P, 45), age_index(P, 55) + 1),
        "fertility_complete": range(P.A_f_end, P.J),
        "65_75": range(age_index(P, 65), age_index(P, 75) + 1),
        "all": range(P.J),
    }


def age_index(P, real_age: int) -> int:
    return int(np.clip(round(real_age - P.age_start), 0, P.J - 1))


def location_indices(P, location: str) -> list[int]:
    if location == "all":
        return list(range(P.I))
    if location == "periphery":
        return [0]
    if location == "center":
        return [1]
    raise ValueError(f"Unknown location {location!r}")


def tenure_match(ten: int, tenure: str) -> bool:
    if tenure == "all":
        return True
    if tenure == "renter":
        return ten == 0
    if tenure == "owner":
        return ten > 0
    raise ValueError(f"Unknown tenure {tenure!r}")


def child_match(P, nn: int, cs: int, concept: str, child_bin: str) -> bool:
    if child_bin == "all":
        return True
    if concept == "current_children":
        value = current_children(P, nn, cs)
    elif concept == "completed_fertility":
        value = nn
    else:
        raise ValueError(f"Unknown child concept {concept!r}")
    if child_bin == "0":
        return value == 0
    if child_bin == "1":
        return value == 1
    if child_bin == "2plus":
        return value >= 2
    raise ValueError(f"Unknown child bin {child_bin!r}")


def current_children(P, nn: int, cs: int) -> int:
    if nn <= 0 or cs == 0 or cs > P.n_child_stages:
        return 0
    return int(nn)


def concat(parts: list[np.ndarray]) -> np.ndarray:
    if not parts:
        return np.array([], dtype=float)
    return np.concatenate(parts).astype(float, copy=False)


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    denom = float(np.sum(w))
    return float(np.sum(x * w) / denom) if denom > 0 else math.nan


def weighted_share(mask: np.ndarray, w: np.ndarray) -> float:
    denom = float(np.sum(w))
    return float(np.sum(w[mask]) / denom) if denom > 0 else math.nan


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    if x.size == 0 or float(np.sum(w)) <= 0:
        return math.nan
    order = np.argsort(x)
    xs = x[order]
    ws = w[order]
    cw = np.cumsum(ws)
    cutoff = q * float(cw[-1])
    idx = int(np.searchsorted(cw, cutoff, side="left"))
    return float(xs[min(idx, xs.size - 1)])


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    for row in rows[1:]:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
