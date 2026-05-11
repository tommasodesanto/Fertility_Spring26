"""Export distributional model moments from the active Python benchmark.

This is deliberately a diagnostic/export tool, not a calibration objective. It
solves the current direct-geometry benchmark and writes model analogs for the
income, wealth, room, tenure, location, and mobility objects used in the
distributional discipline report.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np

MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from dt_cp_model.direct_calibration import (
    DIRECT_GEOMETRY_NAMES,
    OUTSIDE_VALUE_NAME,
    RENEWAL_FLOW_NAME,
    build_direct_calibration_setup,
)
from dt_cp_model.parameters import asdict
from dt_cp_model.solver import make_grid, run_model_cp_dt
from dt_cp_model.theta import apply_theta


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BEST = (
    REPO_ROOT
    / "code"
    / "cluster"
    / "results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506"
    / "direct_geometry_best.json"
)
DEFAULT_OUT = REPO_ROOT / "output" / "model" / "distributional_discipline_current"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--best-json", type=Path, default=DEFAULT_BEST)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--rank", type=int, default=0)
    parser.add_argument("--max-iter-eq", type=int, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    record = load_best_record(args.best_json, args.rank)
    sol, P, p_eq = solve_record(record, max_iter_eq=args.max_iter_eq, verbose=not args.quiet)
    b_grid = make_grid(P)

    export_metadata(args.outdir, record, sol, P, p_eq)
    export_wealth_bins(args.outdir, sol, P, b_grid)
    export_income_bins(args.outdir, sol, P, b_grid)
    export_room_distribution(args.outdir, sol, P, b_grid)
    export_fertility_slices(args.outdir, sol, P, b_grid)
    export_location_tenure(args.outdir, sol, P, b_grid)
    export_location_transitions(args.outdir, sol, P)
    export_mobility(args.outdir, sol, P)


def load_best_record(path: Path, rank: int) -> dict:
    records = json.loads(path.read_text())
    if not isinstance(records, list):
        records = [records]
    if rank < 0 or rank >= len(records):
        raise IndexError(f"rank {rank} outside available records 0..{len(records)-1}")
    return records[rank]


def solve_record(record: dict, *, max_iter_eq: int | None, verbose: bool):
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        geo_weight=100.0,
        population_closure="renewal_valve_calibrated",
        scale_target=1.0,
        renewal_retention=float(record.get("moments", {}).get("renewal_retention", 1.0)),
    )
    theta = np.asarray(record["theta"], dtype=float)
    theta_dict = {name: float(value) for name, value in zip(setup.names, theta)}
    P = SimpleNamespace(**asdict(setup.P_base))
    structural_names = [
        name
        for name in setup.names
        if name not in DIRECT_GEOMETRY_NAMES and name not in (OUTSIDE_VALUE_NAME, RENEWAL_FLOW_NAME)
    ]
    structural_theta = [theta_dict[name] for name in structural_names]
    P = apply_theta(P, structural_theta, structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    if max_iter_eq is not None:
        P.max_iter_eq = int(max_iter_eq)
    return run_model_cp_dt(P, verbose=verbose)


def export_metadata(outdir: Path, record: dict, sol, P, p_eq) -> None:
    rows = [
        ("source_best_json", str(DEFAULT_BEST)),
        ("job_id", record.get("job_id")),
        ("eval_id", record.get("eval_id")),
        ("loss", record.get("loss")),
        ("solve_ok_record", record.get("solve_ok")),
        ("p_P", float(p_eq[0])),
        ("p_C", float(p_eq[1])),
        ("tfr", float(2 * sol.mean_parity)),
        ("childless_rate", float(sol.parity_dist[0])),
        ("own_rate_all_ages", float(sol.own_rate)),
        ("own_rate_3055", float(getattr(sol, "own_rate_3055", np.nan))),
        ("pop_share_C", float(sol.pop_share[1])),
        ("rent_ratio_C_over_P", float(p_eq[1] / p_eq[0])),
        ("population_closure", getattr(P, "population_closure", "")),
        ("equilibrium_reason", getattr(sol, "timings", {}).get("convergence_reason", "")),
        ("final_eq_error", getattr(sol, "timings", {}).get("final_eq_error", np.nan)),
    ]
    write_rows(outdir / "model_metadata.csv", ["key", "value"], rows)
    (outdir / "theta.json").write_text(json.dumps(record.get("parameters", {}), indent=2, sort_keys=True))


def export_wealth_bins(outdir: Path, sol, P, b_grid: np.ndarray) -> None:
    rows: list[list] = []
    for age_label, ages in age_windows(P).items():
        for wealth_concept in ["liquid", "sale_net", "accounting_networth"]:
            for group_name, mask_fn in groups(P).items():
                cells = collect_cells(sol, P, b_grid, ages=ages, mask_fn=mask_fn, wealth_concept=wealth_concept)
                rows.extend(bin_rows("wealth", age_label, wealth_concept, group_name, cells, nbins=5))
    write_rows(
        outdir / "model_wealth_bins.csv",
        [
            "domain",
            "age_bin",
            "concept",
            "group",
            "bin",
            "mass",
            "mass_share",
            "x_mean",
            "x_median",
            "income_mean",
            "parent_share",
            "childless_share",
            "mean_completed_fertility",
            "owner_rate",
            "center_share",
            "mean_rooms",
        ],
        rows,
    )


def export_income_bins(outdir: Path, sol, P, b_grid: np.ndarray) -> None:
    rows: list[list] = []
    for age_label, ages in {
        "ages25_30": range(age_index(P, 25), age_index(P, 30) + 1),
        "ages25_45": range(age_index(P, 25), age_index(P, 45) + 1),
        "ages35_45": range(age_index(P, 35), age_index(P, 45) + 1),
        "age45": range(age_index(P, 45), age_index(P, 45) + 1),
    }.items():
        cells = collect_cells(sol, P, b_grid, ages=ages, mask_fn=lambda *_: True, wealth_concept="income")
        rows.extend(bin_rows("income", age_label, "current_income", "all", cells, nbins=5))
    write_rows(
        outdir / "model_income_bins.csv",
        [
            "domain",
            "age_bin",
            "concept",
            "group",
            "bin",
            "mass",
            "mass_share",
            "x_mean",
            "x_median",
            "income_mean",
            "parent_share",
            "childless_share",
            "mean_completed_fertility",
            "owner_rate",
            "center_share",
            "mean_rooms",
        ],
        rows,
    )


def export_room_distribution(outdir: Path, sol, P, b_grid: np.ndarray) -> None:
    rows_summary: list[list] = []
    rows_bins: list[list] = []
    for age_label, ages in {"all": range(P.J), "25_45": range(age_index(P, 25), age_index(P, 45) + 1)}.items():
        for tenure_name, tenures in {"Renter": [0], "Owner": list(range(1, 1 + P.n_house))}.items():
            for child_label, child_filter in child_bins().items():
                vals, wts = room_cells(sol, P, b_grid, ages, tenures, child_filter)
                if wts.size == 0 or float(np.sum(wts)) <= 0:
                    continue
                rows_summary.append(
                    [
                        age_label,
                        tenure_name,
                        child_label,
                        float(np.sum(wts)),
                        weighted_mean(vals, wts),
                        weighted_quantile(vals, wts, 0.25),
                        weighted_quantile(vals, wts, 0.50),
                        weighted_quantile(vals, wts, 0.75),
                        weighted_share(vals >= 7, wts),
                        weighted_share(vals >= 8, wts),
                        weighted_share(vals >= 11, wts),
                    ]
                )
                for bin_name, bin_mask in room_bin_masks(vals).items():
                    mass = float(np.sum(wts[bin_mask]))
                    rows_bins.append([age_label, tenure_name, child_label, bin_name, mass, mass / max(float(np.sum(wts)), 1e-14)])
    write_rows(
        outdir / "model_rooms_summary.csv",
        ["age_window", "tenure", "child_bin", "mass", "mean", "p25", "p50", "p75", "share_ge_7", "share_ge_8", "share_ge_11"],
        rows_summary,
    )
    write_rows(outdir / "model_rooms_bins.csv", ["age_window", "tenure", "child_bin", "bin", "mass", "share"], rows_bins)


def export_location_tenure(outdir: Path, sol, P, b_grid: np.ndarray) -> None:
    rows = []
    for label, child_filter in {
        "Non-Parents": lambda nn, cs: current_children(nn, cs, P) == 0,
        "New Parents": lambda nn, cs: nn >= 1 and 1 <= cs <= min(P.n_child_stages, 2),
        "Older Parents": lambda nn, cs: nn >= 1 and not (1 <= cs <= min(P.n_child_stages, 2)),
        "All Parents": lambda nn, cs: nn >= 1,
    }.items():
        vals, wts = cell_arrays(sol, P, b_grid, ages=range(age_index(P, 22), age_index(P, 45) + 1), child_filter=child_filter)
        if wts.size == 0:
            continue
        rows.append(
            [
                label,
                float(np.sum(wts)),
                weighted_mean(vals["center"], wts),
                weighted_mean(vals["owner"], wts),
                weighted_mean(vals["rooms"], wts),
                weighted_mean(vals["completed_fertility"], wts),
            ]
        )
    write_rows(
        outdir / "model_location_tenure_by_parent.csv",
        ["parent_status", "mass", "center_share", "owner_rate", "mean_rooms", "mean_completed_fertility"],
        rows,
    )


def export_fertility_slices(outdir: Path, sol, P, b_grid: np.ndarray) -> None:
    rows = []
    header = [
        "slice",
        "age_window",
        "age",
        "location",
        "tenure",
        "mass",
        "current_parent_share",
        "completed_parent_share",
        "childless_completed_share",
        "mean_completed_fertility",
        "tfr",
        "owner_rate",
        "center_share",
        "mean_rooms",
    ]
    locations = {
        "All": lambda ten, loc, nn, cs: True,
        "Periphery": lambda ten, loc, nn, cs: loc == 0,
        "Center": lambda ten, loc, nn, cs: loc == 1,
    }
    tenures = {
        "All": lambda ten, loc, nn, cs: True,
        "Renter": lambda ten, loc, nn, cs: ten == 0,
        "Owner": lambda ten, loc, nn, cs: ten > 0,
    }

    for age in range(22, 46):
        j = age_index(P, age)
        for location, loc_filter in locations.items():
            vals, wts = cell_arrays(sol, P, b_grid, ages=[j], mask_fn=loc_filter)
            rows.append(["age", "single_age", age, location, "All"] + fertility_summary(vals, wts))

    for window, ages in {
        "ages22_45": range(age_index(P, 22), age_index(P, 45) + 1),
        "ages25_45": range(age_index(P, 25), age_index(P, 45) + 1),
        "ages30_45": range(age_index(P, 30), age_index(P, 45) + 1),
        "age45": range(age_index(P, 45), age_index(P, 45) + 1),
    }.items():
        for location, loc_filter in locations.items():
            vals, wts = cell_arrays(sol, P, b_grid, ages=ages, mask_fn=loc_filter)
            rows.append(["location", window, "", location, "All"] + fertility_summary(vals, wts))

        for tenure, tenure_filter in tenures.items():
            vals, wts = cell_arrays(sol, P, b_grid, ages=ages, mask_fn=tenure_filter)
            rows.append(["tenure", window, "", "All", tenure] + fertility_summary(vals, wts))

        for location, loc_filter in locations.items():
            for tenure, tenure_filter in tenures.items():
                if location == "All" and tenure == "All":
                    continue
                vals, wts = cell_arrays(
                    sol,
                    P,
                    b_grid,
                    ages=ages,
                    mask_fn=lambda ten, loc, nn, cs, lf=loc_filter, tf=tenure_filter: lf(ten, loc, nn, cs) and tf(ten, loc, nn, cs),
                )
                rows.append(["location_tenure", window, "", location, tenure] + fertility_summary(vals, wts))

    write_rows(outdir / "model_fertility_slices.csv", header, rows)


def fertility_summary(vals: dict[str, np.ndarray], wts: np.ndarray) -> list[float]:
    if wts.size == 0 or float(np.sum(wts)) <= 0:
        return [np.nan] * 9
    return [
        float(np.sum(wts)),
        weighted_mean(vals["current_parent"], wts),
        weighted_mean(vals["parent"], wts),
        weighted_mean(vals["childless"], wts),
        weighted_mean(vals["completed_fertility"], wts),
        2.0 * weighted_mean(vals["completed_fertility"], wts),
        weighted_mean(vals["owner"], wts),
        weighted_mean(vals["center"], wts),
        weighted_mean(vals["rooms"], wts),
    ]


def export_location_transitions(outdir: Path, sol, P) -> None:
    rows = []
    labels = {0: "Periphery", 1: "Center"}
    for group, child_filter in {
        "All Parents": lambda nn, cs: nn >= 1,
        "Non-Parents": lambda nn, cs: current_children(nn, cs, P) == 0,
    }.items():
        for origin in range(P.I):
            flows = np.zeros(P.I)
            total = 0.0
            for j in range(min(P.J - 1, age_index(P, 22)), min(P.J - 1, age_index(P, 45)) + 1):
                for ten in range(1 + P.n_house):
                    for nn in range(P.n_parity):
                        for cs in range(P.n_child_states):
                            if not child_filter(nn, cs):
                                continue
                            mass_vec = sol.g[:, ten, origin, j, nn, cs]
                            mass = float(np.sum(mass_vec))
                            if mass <= 0:
                                continue
                            total += mass
                            for dest in range(P.I):
                                probs = sol.loc_probs[:, ten, origin, dest, j, nn, cs]
                                flows[dest] += float(np.sum(mass_vec * probs))
            for dest in range(P.I):
                rows.append([group, labels[origin], labels[dest], flows[dest], flows[dest] / max(total, 1e-14)])
    write_rows(
        outdir / "model_origin_transition_by_parent.csv",
        ["parent_compare_all", "origin_label", "dest_label", "model_flow", "share"],
        rows,
    )


def export_mobility(outdir: Path, sol, P) -> None:
    rows = []
    for label, child_filter in {
        "Non-Parents": lambda nn, cs: current_children(nn, cs, P) == 0,
        "Current Parents": lambda nn, cs: current_children(nn, cs, P) > 0,
        "Completed Parents": lambda nn, cs: nn >= 1,
    }.items():
        num = den = 0.0
        for j in range(min(P.J - 1, age_index(P, 22)), min(P.J - 1, age_index(P, 45)) + 1):
            for i in range(P.I):
                for ten in range(1 + P.n_house):
                    for nn in range(P.n_parity):
                        for cs in range(P.n_child_states):
                            if not child_filter(nn, cs):
                                continue
                            mass_vec = sol.g[:, ten, i, j, nn, cs]
                            mass = float(np.sum(mass_vec))
                            if mass <= 0:
                                continue
                            stay_prob = sol.loc_probs[:, ten, i, i, j, nn, cs]
                            num += float(np.sum(mass_vec * (1 - stay_prob)))
                            den += mass
        rows.append([label, den, num / max(den, 1e-14), "location-change probability; model has no move-for-size reason code"])
    write_rows(outdir / "model_mobility_by_parent.csv", ["group", "mass", "model_move_probability", "interpretation"], rows)


def age_index(P, real_age: int) -> int:
    return int(np.clip(round(real_age - P.age_start), 0, P.J - 1))


def age_windows(P):
    return {
        "ages25_30": range(age_index(P, 25), age_index(P, 30) + 1),
        "age35": range(age_index(P, 35), age_index(P, 35) + 1),
        "age45": range(age_index(P, 45), age_index(P, 45) + 1),
        "ages25_45": range(age_index(P, 25), age_index(P, 45) + 1),
        "ages35_45": range(age_index(P, 35), age_index(P, 45) + 1),
        "ages45_55": range(age_index(P, 45), age_index(P, 55) + 1),
    }


def groups(P):
    return {
        "all": lambda ten, loc, nn, cs: True,
        "renters": lambda ten, loc, nn, cs: ten == 0,
        "owners": lambda ten, loc, nn, cs: ten > 0,
        "periphery": lambda ten, loc, nn, cs: loc == 0,
        "center": lambda ten, loc, nn, cs: loc == 1,
        "current_parents": lambda ten, loc, nn, cs: current_children(nn, cs, P) > 0,
        "nonparents": lambda ten, loc, nn, cs: current_children(nn, cs, P) == 0,
    }


def collect_cells(sol, P, b_grid, *, ages: Iterable[int], mask_fn, wealth_concept: str) -> dict[str, np.ndarray]:
    vals, wts = cell_arrays(sol, P, b_grid, ages=ages, mask_fn=mask_fn)
    if vals["weight"].size == 0:
        return vals
    if wealth_concept == "liquid":
        x = vals["liquid"]
    elif wealth_concept == "sale_net":
        x = vals["sale_net"]
    elif wealth_concept == "accounting_networth":
        x = vals["accounting_networth"]
    elif wealth_concept == "income":
        x = vals["income"]
    else:
        raise ValueError(wealth_concept)
    vals["x"] = x
    vals["weight"] = wts
    return vals


def cell_arrays(sol, P, b_grid, *, ages: Iterable[int], mask_fn=None, child_filter=None):
    cols = {k: [] for k in ["weight", "liquid", "sale_net", "accounting_networth", "income", "parent", "current_parent", "childless", "completed_fertility", "owner", "center", "rooms"]}
    for j in ages:
        if j < 0 or j >= P.J:
            continue
        for loc in range(P.I):
            for ten in range(1 + P.n_house):
                heq = (1 - P.psi) * sol.p_eq[loc] * P.H_own[ten - 1] if ten > 0 else 0.0
                owner = 1.0 if ten > 0 else 0.0
                rooms_owner = P.H_own[ten - 1] if ten > 0 else None
                income = float(P.income[loc, j])
                for nn in range(P.n_parity):
                    for cs in range(P.n_child_states):
                        if mask_fn is not None and not mask_fn(ten, loc, nn, cs):
                            continue
                        if child_filter is not None and not child_filter(nn, cs):
                            continue
                        w = sol.g[:, ten, loc, j, nn, cs]
                        nz = w > 0
                        if not np.any(nz):
                            continue
                        liq = b_grid[nz]
                        wt = w[nz]
                        rooms = sol.hR_pol[nz, ten, loc, j, nn, cs] if ten == 0 else np.full(np.count_nonzero(nz), rooms_owner)
                        fert = np.full(np.count_nonzero(nz), float(nn))
                        parent = fert >= 1
                        current_parent = current_children(nn, cs, P) > 0
                        cols["weight"].append(wt)
                        cols["liquid"].append(liq)
                        cols["sale_net"].append(liq + heq)
                        cols["accounting_networth"].append(liq + sol.p_eq[loc] * P.H_own[ten - 1] if ten > 0 else liq)
                        cols["income"].append(np.full(np.count_nonzero(nz), income))
                        cols["parent"].append(parent.astype(float))
                        cols["current_parent"].append(np.full(np.count_nonzero(nz), float(current_parent)))
                        cols["childless"].append((fert == 0).astype(float))
                        cols["completed_fertility"].append(fert)
                        cols["owner"].append(np.full(np.count_nonzero(nz), owner))
                        cols["center"].append(np.full(np.count_nonzero(nz), 1.0 if loc == 1 else 0.0))
                        cols["rooms"].append(np.asarray(rooms, dtype=float))
    vals = {k: concat(v) for k, v in cols.items()}
    return vals, vals["weight"]


def bin_rows(domain: str, age_label: str, concept: str, group: str, cells: dict[str, np.ndarray], nbins: int):
    w = cells.get("weight", np.array([]))
    x = cells.get("x", np.array([]))
    if w.size == 0 or float(np.sum(w)) <= 0 or x.size == 0:
        return []
    qs = [weighted_quantile(x, w, q / nbins) for q in range(1, nbins)]
    bin_id = np.searchsorted(qs, x, side="right") + 1
    rows = []
    total = float(np.sum(w))
    for b in range(1, nbins + 1):
        m = bin_id == b
        wb = w[m]
        if wb.size == 0 or float(np.sum(wb)) <= 0:
            rows.append([domain, age_label, concept, group, b] + [np.nan] * 11)
            continue
        rows.append(
            [
                domain,
                age_label,
                concept,
                group,
                b,
                float(np.sum(wb)),
                float(np.sum(wb)) / total,
                weighted_mean(x[m], wb),
                weighted_quantile(x[m], wb, 0.50),
                weighted_mean(cells["income"][m], wb),
                weighted_mean(cells["parent"][m], wb),
                weighted_mean(cells["childless"][m], wb),
                weighted_mean(cells["completed_fertility"][m], wb),
                weighted_mean(cells["owner"][m], wb),
                weighted_mean(cells["center"][m], wb),
                weighted_mean(cells["rooms"][m], wb),
            ]
        )
    return rows


def room_cells(sol, P, b_grid, ages, tenures, child_filter):
    vals, wts = cell_arrays(sol, P, b_grid, ages=ages, mask_fn=lambda ten, loc, nn, cs: ten in tenures, child_filter=child_filter)
    return vals["rooms"], wts


def child_bins():
    return {
        "0": lambda nn, cs: current_children(nn, cs, None) == 0,
        "1": lambda nn, cs: current_children(nn, cs, None) == 1,
        "2+": lambda nn, cs: current_children(nn, cs, None) >= 2,
        "all": lambda nn, cs: True,
    }


def current_children(nn: int, cs: int, P) -> int:
    dep_last = 4 if P is None else P.n_child_stages
    if nn <= 0 or cs == 0 or cs > dep_last:
        return 0
    return int(nn)


def room_bin_masks(vals):
    return {
        "<=4": vals <= 4,
        "5": (vals > 4) & (vals <= 5),
        "6": (vals > 5) & (vals <= 6),
        "7-8": (vals > 6) & (vals <= 8),
        "9-10": (vals > 8) & (vals <= 10),
        "11+": vals > 10,
    }


def weighted_mean(x, w) -> float:
    return float(np.sum(np.asarray(x) * np.asarray(w)) / max(float(np.sum(w)), 1e-14))


def weighted_share(mask, w) -> float:
    return float(np.sum(np.asarray(w)[np.asarray(mask)]) / max(float(np.sum(w)), 1e-14))


def weighted_quantile(x, w, q: float) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return float("nan")
    x = x[ok]
    w = w[ok]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cw = np.cumsum(w)
    return float(x[np.searchsorted(cw, q * cw[-1], side="left")])


def concat(parts):
    return np.concatenate(parts) if parts else np.array([])


def write_rows(path: Path, header: list[str], rows: Iterable[Iterable]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
