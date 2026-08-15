#!/usr/bin/env python3
"""Build and verify current one-shot calendar-transition initial states.

The stationary packet is tied to the July-23 timing-repaired Stone--Geary
estimate.  A separate diagnostic packet reweights only adult-age cells to the
2023 ACS household-head distribution.  No shock or empirical trend is added.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from intergen_housing_fertility_optimized.calibration import diagnostic_loss, extract_moments
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.solver import (
    advance_cohort_one_period_markov_income,
    build_forward_tenure_transition_maps,
    interp_indices,
    precompute_shared,
    realize_current_cross_section,
    run_model_cp_dt,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/current_one_shot_transition_initial_state"
DEFAULT_ACS_EXTRACT = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
)
DEFAULT_MMS_LOOKUP = ROOT / "code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv"

# Source-specific regression anchors from the exact retained solve.  These
# are verification receipts, not empirical targets or transition assumptions.
EXPECTED_CURRENT_EFFECTIVE_ENTRY = 0.061733455656451254
EXPECTED_CURRENT_MATURE_FLOW = 0.060701144771170716


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--acs-extract", type=Path, default=DEFAULT_ACS_EXTRACT)
    parser.add_argument("--mms-lookup", type=Path, default=DEFAULT_MMS_LOOKUP)
    parser.add_argument("--acs-year", type=int, default=2023)
    parser.add_argument("--atol", type=float, default=2e-9)
    parser.add_argument("--verbose-solver", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_acs_head_age_distribution(
    extract_path: Path,
    lookup_path: Path,
    *,
    year: int,
    age_start: int,
    period_years: int,
    n_periods: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Read exact four-year head-age bins from the raw IPUMS ACS extract.

    The Stata file is sorted by year.  We use its fixed-width record metadata
    to seek directly to the requested year, then aggregate only the required
    numeric fields.  This avoids materializing the 9.2 GB extract in memory.
    """

    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - repository venv has pandas
        raise RuntimeError("The repository Python environment must provide pandas") from exc

    age_end = int(age_start + period_years * n_periods - 1)
    lookup: dict[tuple[int, int, int], str] = {}
    original_location: dict[tuple[int, int, int], str] = {}
    with lookup_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if not (int(row["year_start"]) <= year <= int(row["year_end"])):
                continue
            key = (int(row["statefip"]), int(row["puma"]), int(row["cbsacode"]))
            if key in lookup:
                raise ValueError(f"Duplicate MMS lookup key: {key}")
            raw_location = str(row["mms_location"]).strip().lower()
            if raw_location not in {"center", "middle", "periphery"}:
                continue
            original_location[key] = raw_location
            lookup[key] = "center" if raw_location == "middle" else raw_location
    if not lookup:
        raise ValueError(f"No lookup rows cover ACS year {year}: {lookup_path}")

    reader = pd.read_stata(
        extract_path,
        iterator=True,
        convert_categoricals=False,
        convert_dates=False,
    )
    try:
        dtype = getattr(reader, "_dtype", None)
        data_location = getattr(reader, "data_location", None)
        raw_names = getattr(dtype, "names", None)
        if dtype is None or data_location is None or raw_names is None:
            raise RuntimeError("Unsupported pandas StataReader: fixed-width metadata unavailable")
        field = dict(zip(reader.varlist, raw_names))
        required = {
            "year",
            "sample",
            "statefip",
            "puma",
            "met2013",
            "gq",
            "pernum",
            "relate",
            "hhwt",
            "age",
        }
        missing = sorted(required - set(field))
        if missing:
            raise ValueError(f"ACS extract is missing required variables: {missing}")

        record_size = int(dtype.itemsize)
        stream = reader.path_or_buf

        def read_records(start: int, count: int) -> np.ndarray:
            stream.seek(int(data_location) + int(start) * record_size)
            raw = np.frombuffer(
                stream.read(int(count) * record_size), dtype=dtype, count=int(count)
            )
            if reader.byteorder != reader._native_byteorder:
                raw = raw.byteswap().newbyteorder()
            return raw

        def lower_year_bound(target_year: int) -> int:
            lower, upper = 0, int(reader.nobs)
            while lower < upper:
                midpoint = (lower + upper) // 2
                observed_year = int(read_records(midpoint, 1)[field["year"]][0])
                if observed_year < target_year:
                    lower = midpoint + 1
                else:
                    upper = midpoint
            return lower

        first = lower_year_bound(int(year))
        last = lower_year_bound(int(year) + 1)
        if first >= last:
            raise ValueError(f"ACS extract has no rows for year {year}")

        weights = np.zeros(n_periods, dtype=float)
        records = np.zeros(n_periods, dtype=np.int64)
        location_weights = {
            "center": np.zeros(n_periods, dtype=float),
            "periphery": np.zeros(n_periods, dtype=float),
        }
        location_records = {
            "center": np.zeros(n_periods, dtype=np.int64),
            "periphery": np.zeros(n_periods, dtype=np.int64),
        }
        year_sample_codes: set[int] = set()
        eligible_head_records = 0
        metro_head_records = 0
        matched_records = 0
        absorbed_middle_records = 0
        absorbed_middle_weight = 0.0
        chunk_size = 250_000
        f = {name: field[name] for name in required}
        for start in range(first, last, chunk_size):
            raw = read_records(start, min(chunk_size, last - start))
            year_sample_codes.update(int(value) for value in np.unique(raw[f["sample"]]))
            base = (
                (raw[f["year"]] == year)
                & np.isin(raw[f["gq"]], [1, 2])
                & (raw[f["pernum"]] == 1)
                & (raw[f["relate"]] == 1)
                & np.isfinite(raw[f["hhwt"]])
                & (raw[f["hhwt"]] > 0.0)
                & (raw[f["age"]] >= age_start)
                & (raw[f["age"]] <= age_end)
            )
            eligible_head_records += int(np.sum(base))
            metro = base & (raw[f["met2013"]] > 0)
            indices = np.flatnonzero(metro)
            metro_head_records += int(indices.size)
            for row_index in indices:
                key = (
                    int(raw[f["statefip"]][row_index]),
                    int(raw[f["puma"]][row_index]),
                    int(raw[f["met2013"]][row_index]),
                )
                location = lookup.get(key)
                if location is None:
                    continue
                age = int(raw[f["age"]][row_index])
                age_index = (age - age_start) // period_years
                weight = float(raw[f["hhwt"]][row_index])
                weights[age_index] += weight
                records[age_index] += 1
                location_weights[location][age_index] += weight
                location_records[location][age_index] += 1
                matched_records += 1
                if original_location[key] == "middle":
                    absorbed_middle_records += 1
                    absorbed_middle_weight += weight
    finally:
        reader.close()

    total_weight = float(np.sum(weights))
    if total_weight <= 0.0 or np.any(weights <= 0.0):
        raise ValueError(
            f"ACS head-age bins must all have positive weight; observed {weights.tolist()}"
        )
    shares = weights / total_weight
    receipt = {
        "year": int(year),
        "ipums_sample_codes": sorted(year_sample_codes),
        "raw_year_records": int(last - first),
        "eligible_head_records_before_metro_filter": int(eligible_head_records),
        "metro_head_records_before_mms_match": int(metro_head_records),
        "matched_mms_head_records": int(matched_records),
        "unmatched_metro_head_records": int(metro_head_records - matched_records),
        "weighted_head_mass": total_weight,
        "age_start": int(age_start),
        "age_end": int(age_end),
        "period_years": int(period_years),
        "age_bin_rule": "closed integer bins [18+4j, 21+4j], j=0,...,16",
        "sample_filter": (
            f"YEAR={year}; GQ in {{1,2}}; PERNUM=1; RELATE=1; HHWT>0; "
            f"AGE={age_start}..{age_end}; MET2013>0; matched 2020-PUMA MMS geography"
        ),
        "weight": "HHWT",
        "middle_geography_rule": "middle MMS PUMAs absorbed into center, matching the ownership-target audit convention",
        "absorbed_middle_records": int(absorbed_middle_records),
        "absorbed_middle_weight": float(absorbed_middle_weight),
        "age_bin_records": records,
        "age_bin_weights": weights,
        "age_bin_shares": shares,
        "center_age_bin_records": location_records["center"],
        "center_age_bin_weights": location_weights["center"],
        "periphery_age_bin_records": location_records["periphery"],
        "periphery_age_bin_weights": location_weights["periphery"],
        "pandas_version": str(pd.__version__),
        "extract_path": str(extract_path),
        "extract_sha256": sha256(extract_path),
        "lookup_path": str(lookup_path),
        "lookup_sha256": sha256(lookup_path),
    }
    return shares, receipt


def reweight_age_cells(distribution: np.ndarray, target_age_shares: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Replace age-cell masses while preserving every within-age conditional state."""

    age_mass = np.sum(distribution, axis=(0, 1, 2, 4, 5, 6))
    target = np.asarray(target_age_shares, dtype=float).reshape(-1)
    if target.size != distribution.shape[3] or np.any(target < 0.0):
        raise ValueError("Observed age shares must be nonnegative and match the model age axis")
    if np.any(age_mass <= 0.0) or not np.isclose(np.sum(target), 1.0):
        raise ValueError("Cannot reweight zero-mass ages or unnormalized observed shares")
    factors = target / age_mass
    return distribution * factors.reshape(1, 1, 1, -1, 1, 1, 1), factors


def max_within_age_conditional_l1(original: np.ndarray, reweighted: np.ndarray) -> float:
    gaps: list[float] = []
    for age_index in range(original.shape[3]):
        source = original[:, :, :, age_index, :, :, :]
        target = reweighted[:, :, :, age_index, :, :, :]
        gaps.append(
            float(
                np.sum(
                    np.abs(
                        source / float(np.sum(source))
                        - target / float(np.sum(target))
                    )
                )
            )
        )
    return max(gaps)


def housing_demand_by_market(distribution: np.ndarray, solution: Any, parameters: Any) -> np.ndarray:
    """Measure occupied housing services using the retained reporting convention."""

    demand = np.zeros(int(parameters.I), dtype=float)
    renter_policy = np.asarray(solution.hR_pol, dtype=float)
    for location in range(int(parameters.I)):
        demand[location] += float(
            np.sum(
                distribution[:, 0, location, :, :, :, :]
                * renter_policy[:, 0, location, :, :, :, :]
            )
        )
        for tenure in range(1, 1 + int(parameters.n_house)):
            demand[location] += float(
                np.sum(distribution[:, tenure, location, :, :, :, :])
            ) * float(parameters.H_own[tenure - 1])
    return demand


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def load_selected_source(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("profile") != NEW_MOMENT_PROFILE_NAME:
        raise ValueError(
            f"Expected profile {NEW_MOMENT_PROFILE_NAME!r}, found {payload.get('profile')!r}."
        )
    selected = payload.get("selected")
    if not isinstance(selected, dict) or selected.get("status") != "ok":
        raise ValueError(f"Source has no successful selected result: {path}")
    theta = selected.get("theta")
    if not isinstance(theta, dict) or not theta:
        raise ValueError(f"Selected source has no theta: {path}")
    return payload, selected


def build_transition_maps(solution: Any, parameters: Any, price: np.ndarray) -> tuple[Any, ...]:
    """Rebuild the exact housing-transaction maps used by the retained KFE."""

    b_grid = np.asarray(solution.b_grid, dtype=float)
    nb = b_grid.size
    nt = 1 + int(parameters.n_house)
    nloc = int(parameters.I)
    shared = precompute_shared(parameters, b_grid)

    purchase_cost = np.zeros((nloc, nt))
    sale_proceeds = np.zeros((nloc, nt))
    for location in range(nloc):
        for tenure in range(1, nt):
            housing = float(parameters.H_own[tenure - 1])
            purchase_cost[location, tenure] = float(price[location]) * housing
            sale_proceeds[location, tenure] = (
                (1.0 - float(parameters.psi)) * float(price[location]) * housing
            )

    move_idx = np.zeros((nloc, nt, nb), dtype=np.int64)
    move_weight = np.zeros((nloc, nt, nb))
    for origin in range(nloc):
        for tenure in range(nt):
            moved_assets = np.clip(
                b_grid + sale_proceeds[origin, tenure], b_grid[0], b_grid[-1]
            )
            move_idx[origin, tenure], move_weight[origin, tenure] = interp_indices(
                b_grid, moved_assets
            )

    transaction_idx, transaction_weight = build_forward_tenure_transition_maps(
        parameters,
        b_grid,
        purchase_cost,
        sale_proceeds,
        shared.phi_choice,
        shared.birth_dp,
        shared.birth_entry_grant,
    )
    return shared, move_idx, move_weight, transaction_idx, transaction_weight


def apply_fertility_choice(cohort: np.ndarray, age_index: int, solution: Any, parameters: Any) -> np.ndarray:
    """Apply the current-age one-shot family-size draw without other choices."""

    out = np.array(cohort, dtype=float, copy=True)
    if not (parameters.A_f_start <= age_index + 1 <= parameters.A_f_end):
        return out

    for income_state in range(out.shape[3]):
        childless = out[:, :, :, income_state, 0, 0].copy()
        probs = np.asarray(
            solution.fert_probs[:, :, :, age_index, income_state, :], dtype=float
        )
        split = childless[:, :, :, None] * probs
        out[:, :, :, income_state, 0, 0] = split[:, :, :, 0]
        out[:, :, :, income_state, 1:, 1] += split[:, :, :, 1:]
    return out


def advance_age(
    cohort: np.ndarray,
    age_index: int,
    solution: Any,
    parameters: Any,
    maps: tuple[Any, ...],
) -> np.ndarray:
    """Advance a post-fertility cohort to next date's pre-fertility state."""

    shared, move_idx, move_weight, transaction_idx, transaction_weight = maps
    survival = (
        float(parameters.survival_probs[age_index])
        if bool(parameters.use_age_survival)
        else 1.0
    )
    return advance_cohort_one_period_markov_income(
        survival * cohort,
        age_index,
        solution.loc_probs,
        solution.tenure_choice,
        solution.tenure_probs,
        solution.bp_pol,
        parameters,
        np.asarray(solution.b_grid, dtype=float),
        shared,
        move_idx,
        move_weight,
        transaction_idx,
        transaction_weight,
        bool(parameters.use_stochastic_aging and hasattr(parameters, "Pi_child")),
        np.asarray(parameters.Pi_child, dtype=float),
        np.asarray(solution.income_transition, dtype=float),
    )


def reconstruct_pre_fertility_state(
    post_fertility: np.ndarray,
    solution: Any,
    parameters: Any,
    maps: tuple[Any, ...],
) -> np.ndarray:
    """Recover the start-of-date state implied by the stationary retained KFE."""

    pre_fertility = np.zeros_like(post_fertility)

    # Entrants are childless before the age-18 family-size draw.  Collapsing
    # the solved age-18 split exactly preserves the normalized entrant wealth,
    # location, tenure, and income distribution, including frontier censoring.
    entry_by_cell = np.sum(post_fertility[:, :, :, 0, :, :, :], axis=(-2, -1))
    pre_fertility[:, :, :, 0, :, 0, 0] = entry_by_cell

    # Every older start-of-date cohort is last period's post-fertility cohort
    # after housing, saving, income, survival, and child-stage transitions.
    for age_index in range(1, int(parameters.J)):
        pre_fertility[:, :, :, age_index, :, :, :] = advance_age(
            post_fertility[:, :, :, age_index - 1, :, :, :],
            age_index - 1,
            solution,
            parameters,
            maps,
        )
    return pre_fertility


def apply_fertility_cross_section(pre_fertility: np.ndarray, solution: Any, parameters: Any) -> np.ndarray:
    out = np.zeros_like(pre_fertility)
    for age_index in range(int(parameters.J)):
        out[:, :, :, age_index, :, :, :] = apply_fertility_choice(
            pre_fertility[:, :, :, age_index, :, :, :],
            age_index,
            solution,
            parameters,
        )
    return out


def one_calendar_step(
    post_fertility: np.ndarray,
    stationary_entry_pre_fertility: np.ndarray,
    solution: Any,
    parameters: Any,
    maps: tuple[Any, ...],
) -> np.ndarray:
    """Advance all incumbent ages once and insert the stationary age-18 flow."""

    next_pre = np.zeros_like(post_fertility)
    next_pre[:, :, :, 0, :, :, :] = stationary_entry_pre_fertility
    for age_index in range(int(parameters.J) - 1):
        next_pre[:, :, :, age_index + 1, :, :, :] = advance_age(
            post_fertility[:, :, :, age_index, :, :, :],
            age_index,
            solution,
            parameters,
            maps,
        )
    return next_pre


def child_clock_receipt(post_fertility: np.ndarray, solution: Any, parameters: Any) -> dict[str, Any]:
    """Return exact household and model reproduction-unit maturation flows."""

    child_transition = np.asarray(parameters.Pi_child, dtype=float)
    dependent_state = int(parameters.n_child_stages)
    household_flow = 0.0
    reproduction_flow = 0.0
    by_family_size: list[dict[str, float]] = []
    for family_size in range(int(parameters.n_parity)):
        matured_state = (
            0
            if family_size == 0
            else int(parameters.n_child_stages) + 1
            if family_size == 1
            else int(parameters.n_child_stages) + 2
        )
        maturation_probability = float(
            child_transition[dependent_state, matured_state, family_size]
        )
        surviving_dependent_mass = 0.0
        for age_index in range(int(parameters.J) - 1):
            survival = (
                float(parameters.survival_probs[age_index])
                if bool(parameters.use_age_survival)
                else 1.0
            )
            surviving_dependent_mass += survival * float(
                np.sum(
                    post_fertility[
                        :, :, :, age_index, :, family_size, dependent_state
                    ]
                )
            )
        matured_households = maturation_probability * surviving_dependent_mass
        matured_reproduction_units = float(family_size) * matured_households
        household_flow += matured_households
        reproduction_flow += matured_reproduction_units
        by_family_size.append(
            {
                "family_size_bin": family_size,
                "surviving_dependent_household_mass": surviving_dependent_mass,
                "maturation_probability": maturation_probability,
                "matured_household_flow": matured_households,
                "matured_model_reproduction_units": matured_reproduction_units,
            }
        )
    return {
        "matured_household_flow": household_flow,
        "matured_model_reproduction_units": reproduction_flow,
        "solver_matured_model_reproduction_units": float(solution.entrants_mature_total),
        "by_family_size": by_family_size,
    }


def make_check(value: float, target: float, tolerance: float) -> dict[str, Any]:
    gap = float(value) - float(target)
    return {
        "value": float(value),
        "target": float(target),
        "gap": gap,
        "tolerance": float(tolerance),
        "passed": bool(abs(gap) <= tolerance),
    }


def write_age_summary(path: Path, pre: np.ndarray, post: np.ndarray, current: np.ndarray, parameters: Any) -> None:
    rows: list[dict[str, Any]] = []
    dependent_state = int(parameters.n_child_stages)
    matured_one = dependent_state + 1
    matured_two = dependent_state + 2
    for age_index in range(int(parameters.J)):
        post_age = post[:, :, :, age_index, :, :, :]
        family_mass = np.sum(post_age, axis=(0, 1, 2, 3, 5))
        dependent_mass = float(np.sum(post_age[..., dependent_state]))
        dependent_reproduction_units = sum(
            family_size
            * float(np.sum(post_age[:, :, :, :, family_size, dependent_state]))
            for family_size in range(int(parameters.n_parity))
        )
        matured_reproduction_units = float(
            np.sum(post_age[:, :, :, :, 1, matured_one])
        )
        if int(parameters.n_parity) >= 3:
            matured_reproduction_units += 2.0 * float(
                np.sum(post_age[:, :, :, :, 2, matured_two])
            )
        rows.append(
            {
                "age_index": age_index,
                "age": float(parameters.age_start + age_index * parameters.period_years),
                "pre_fertility_mass": float(np.sum(pre[:, :, :, age_index, :, :, :])),
                "post_fertility_pre_housing_mass": float(np.sum(post_age)),
                "realized_post_housing_mass": float(
                    np.sum(current[:, :, :, age_index, :, :, :])
                ),
                "family_size_0_mass": float(family_mass[0]),
                "family_size_bin1_mass": float(family_mass[1]),
                "family_size_bin2_mass": float(family_mass[2]),
                "dependent_household_mass": dependent_mass,
                "dependent_model_reproduction_units": dependent_reproduction_units,
                "matured_model_reproduction_units_stock": matured_reproduction_units,
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_acs_age_bins(path: Path, receipt: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for age_index, share in enumerate(receipt["age_bin_shares"]):
        age_low = int(receipt["age_start"] + receipt["period_years"] * age_index)
        rows.append(
            {
                "age_index": age_index,
                "age_low": age_low,
                "age_high": age_low + int(receipt["period_years"]) - 1,
                "records": int(receipt["age_bin_records"][age_index]),
                "hhwt_sum": float(receipt["age_bin_weights"][age_index]),
                "share": float(share),
                "center_records_middle_absorbed": int(
                    receipt["center_age_bin_records"][age_index]
                ),
                "center_hhwt_sum_middle_absorbed": float(
                    receipt["center_age_bin_weights"][age_index]
                ),
                "periphery_records": int(
                    receipt["periphery_age_bin_records"][age_index]
                ),
                "periphery_hhwt_sum": float(
                    receipt["periphery_age_bin_weights"][age_index]
                ),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_source_to_state(path: Path) -> None:
    rows = [
        {
            "packet_key": "g0_stationary_pre_fertility",
            "source": "reconstructed from the stationary KFE age predecessor and stationary entrant slice",
            "timing": "start of date, before the current-age one-shot family-size draw",
        },
        {
            "packet_key": "g0_stationary_post_fertility_pre_housing",
            "source": "solution.g_beginning_distribution",
            "timing": "after current-age family-size draw; before location, tenure, and housing transaction",
        },
        {
            "packet_key": "g0_stationary_realized_post_housing",
            "source": "solution.g",
            "timing": "after current location, tenure, and housing transaction; before saving transition",
        },
        {
            "packet_key": "entry_stationary_pre_fertility",
            "source": "age-18 solved KFE collapsed back to childless status",
            "timing": "effective normalized entrant flow before the age-18 family-size draw",
        },
        {
            "packet_key": "child_transition",
            "source": "parameters.Pi_child",
            "timing": "applied once per four-year adult transition after saving",
        },
        {
            "packet_key": "g0_observed_age_pre_fertility (and later timing siblings)",
            "source": "stationary timing layers rescaled by 2023 ACS HHWT head-age shares",
            "timing": "same three within-date timings; within-age conditional states unchanged",
        },
        {
            "packet_key": "entry_observed_age_pre_fertility",
            "source": "2023 ACS age-18--21 share times the stationary conditional entrant state",
            "timing": "constant-youngest-cohort diagnostic continuation, not an estimated trend",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    source_path = args.source.resolve()
    outdir = args.outdir.resolve()
    tolerance = float(args.atol)
    if tolerance <= 0.0:
        raise ValueError("--atol must be positive")

    payload, selected = load_selected_source(source_path)
    theta = {str(key): float(value) for key, value in selected["theta"].items()}
    overrides = new_moment_overrides(tight=True, optimized=True)
    overrides.update(theta)
    solution, parameters, price = run_model_cp_dt(
        overrides, verbose=bool(args.verbose_solver)
    )
    price = np.asarray(price, dtype=float).reshape(-1)
    moments = extract_moments(solution, parameters)
    target_system = new_moment_target_system()
    reproduced_loss = diagnostic_loss(
        moments,
        targets=target_system.targets_dict(),
        weights=target_system.weights_dict(),
    )

    post_fertility = np.asarray(solution.g_beginning_distribution, dtype=float)
    current = np.asarray(solution.g, dtype=float)
    maps = build_transition_maps(solution, parameters, price)
    pre_fertility = reconstruct_pre_fertility_state(
        post_fertility, solution, parameters, maps
    )
    recovered_post = apply_fertility_cross_section(
        pre_fertility, solution, parameters
    )
    entry_pre_fertility = pre_fertility[:, :, :, 0, :, :, :].copy()
    next_pre = one_calendar_step(
        recovered_post, entry_pre_fertility, solution, parameters, maps
    )
    next_post = apply_fertility_cross_section(next_pre, solution, parameters)

    _, move_idx, move_weight, transaction_idx, transaction_weight = maps
    recovered_current = realize_current_cross_section(
        post_fertility,
        solution.loc_probs,
        solution.tenure_choice,
        solution.tenure_probs,
        move_idx,
        move_weight,
        transaction_idx,
        transaction_weight,
        use_compiled_scatter=bool(getattr(parameters, "use_numba_scatter", False)),
    )

    acs_age_shares, acs_receipt = read_acs_head_age_distribution(
        args.acs_extract.resolve(),
        args.mms_lookup.resolve(),
        year=int(args.acs_year),
        age_start=int(parameters.age_start),
        period_years=int(parameters.period_years),
        n_periods=int(parameters.J),
    )
    observed_pre, observed_pre_factors = reweight_age_cells(
        pre_fertility, acs_age_shares
    )
    observed_post, observed_post_factors = reweight_age_cells(
        post_fertility, acs_age_shares
    )
    observed_current, observed_current_factors = reweight_age_cells(
        current, acs_age_shares
    )
    observed_recovered_post = apply_fertility_cross_section(
        observed_pre, solution, parameters
    )
    observed_recovered_current = realize_current_cross_section(
        observed_post,
        solution.loc_probs,
        solution.tenure_choice,
        solution.tenure_probs,
        move_idx,
        move_weight,
        transaction_idx,
        transaction_weight,
        use_compiled_scatter=bool(getattr(parameters, "use_numba_scatter", False)),
    )
    observed_entry_pre_fertility = observed_pre[:, :, :, 0, :, :, :].copy()
    observed_next_pre = one_calendar_step(
        observed_recovered_post,
        observed_entry_pre_fertility,
        solution,
        parameters,
        maps,
    )
    observed_next_age_mass = np.sum(
        observed_next_pre, axis=(0, 1, 2, 4, 5, 6)
    )

    child_transition = np.asarray(parameters.Pi_child, dtype=float)
    row_sum_gap = float(np.max(np.abs(np.sum(child_transition, axis=1) - 1.0)))
    child_receipt = child_clock_receipt(post_fertility, solution, parameters)
    observed_child_receipt = child_clock_receipt(
        observed_post, solution, parameters
    )
    effective_entry = float(np.sum(entry_pre_fertility))
    raw_parameter_entry = float(np.sum(parameters.entry_by_loc))
    terminal_exit = float(np.sum(post_fertility[:, :, :, -1, :, :, :]))
    death_exit = 0.0
    for age_index in range(int(parameters.J) - 1):
        death_exit += (1.0 - float(parameters.survival_probs[age_index])) * float(
            np.sum(post_fertility[:, :, :, age_index, :, :, :])
        )

    family_child_mass = np.sum(post_fertility, axis=(0, 1, 2, 3, 4))
    valid = np.zeros_like(family_child_mass, dtype=bool)
    dependent_state = int(parameters.n_child_stages)
    valid[0, 0] = True
    valid[1, dependent_state] = True
    valid[1, dependent_state + 1] = True
    valid[2, dependent_state] = True
    valid[2, dependent_state + 2] = True
    invalid_family_child_mass = float(np.sum(family_child_mass[~valid]))

    stationary_housing_demand = housing_demand_by_market(
        current, solution, parameters
    )
    observed_housing_demand = housing_demand_by_market(
        observed_current, solution, parameters
    )
    supply_scale_at_p0 = (
        parameters.user_cost_rate * price / np.asarray(parameters.r_bar, dtype=float)
    ) ** np.asarray(parameters.xi_supply, dtype=float)
    observed_p0_clearing_supply_intercept = observed_housing_demand / supply_scale_at_p0

    checks = {
        "source_price": make_check(float(price[0]), float(selected["price"]), 1e-10),
        "source_loss": make_check(float(reproduced_loss), float(selected["rank_loss"]), 1e-10),
        "market_residual": {
            "value": float(solution.best_max_abs_rel_excess),
            "target": 0.0,
            "tolerance": float(parameters.tol_eq),
            "passed": bool(solution.best_max_abs_rel_excess <= parameters.tol_eq),
        },
        "pre_fertility_total_mass": make_check(float(np.sum(pre_fertility)), 1.0, tolerance),
        "post_fertility_total_mass": make_check(float(np.sum(post_fertility)), 1.0, tolerance),
        "current_total_mass": make_check(float(np.sum(current)), 1.0, tolerance),
        "fertility_timing_recovery_l1": make_check(
            float(np.sum(np.abs(recovered_post - post_fertility))), 0.0, tolerance
        ),
        "calendar_one_step_pre_fertility_l1": make_check(
            float(np.sum(np.abs(next_pre - pre_fertility))), 0.0, tolerance
        ),
        "calendar_one_step_post_fertility_l1": make_check(
            float(np.sum(np.abs(next_post - post_fertility))), 0.0, tolerance
        ),
        "current_choice_timing_recovery_l1": make_check(
            float(np.sum(np.abs(recovered_current - current))), 0.0, tolerance
        ),
        "entry_replaces_terminal_and_death_exit": make_check(
            effective_entry, terminal_exit + death_exit, tolerance
        ),
        "effective_entry_regression_anchor": make_check(
            effective_entry, EXPECTED_CURRENT_EFFECTIVE_ENTRY, tolerance
        ),
        "mature_child_flow_solver_identity": make_check(
            float(child_receipt["matured_model_reproduction_units"]),
            float(solution.entrants_mature_total),
            tolerance,
        ),
        "mature_child_flow_regression_anchor": make_check(
            float(solution.entrants_mature_total),
            EXPECTED_CURRENT_MATURE_FLOW,
            tolerance,
        ),
        "child_transition_row_sum": make_check(row_sum_gap, 0.0, tolerance),
        "invalid_family_child_state_mass": make_check(
            invalid_family_child_mass, 0.0, tolerance
        ),
        "acs_age_shares_sum": make_check(float(np.sum(acs_age_shares)), 1.0, tolerance),
        "observed_pre_total_mass": make_check(float(np.sum(observed_pre)), 1.0, tolerance),
        "observed_post_total_mass": make_check(float(np.sum(observed_post)), 1.0, tolerance),
        "observed_current_total_mass": make_check(float(np.sum(observed_current)), 1.0, tolerance),
        "observed_age_share_match": make_check(
            float(
                np.max(
                    np.abs(
                        np.sum(observed_pre, axis=(0, 1, 2, 4, 5, 6))
                        - acs_age_shares
                    )
                )
            ),
            0.0,
            tolerance,
        ),
        "observed_within_age_conditionals_preserved": make_check(
            max_within_age_conditional_l1(pre_fertility, observed_pre),
            0.0,
            tolerance,
        ),
        "observed_fertility_timing_recovery_l1": make_check(
            float(np.sum(np.abs(observed_recovered_post - observed_post))),
            0.0,
            tolerance,
        ),
        "observed_current_choice_timing_recovery_l1": make_check(
            float(np.sum(np.abs(observed_recovered_current - observed_current))),
            0.0,
            tolerance,
        ),
        "observed_p0_supply_normalization": make_check(
            float(
                np.max(
                    np.abs(
                        observed_p0_clearing_supply_intercept
                        * supply_scale_at_p0
                        - observed_housing_demand
                    )
                )
            ),
            0.0,
            tolerance,
        ),
    }
    failed = [name for name, receipt in checks.items() if not receipt["passed"]]
    if failed:
        detail = {name: checks[name] for name in failed}
        raise RuntimeError(f"Initial-state validation failed: {json.dumps(detail, indent=2)}")

    outdir.mkdir(parents=True, exist_ok=True)
    age_grid = parameters.age_start + parameters.period_years * np.arange(parameters.J)
    effective_entry_by_location = np.sum(
        entry_pre_fertility, axis=(0, 1, 3, 4, 5)
    )
    common_arrays = {
        "liquid_wealth_grid": np.asarray(solution.b_grid, dtype=float),
        "adult_age_grid": np.asarray(age_grid, dtype=float),
        "income_state_grid": np.asarray(solution.type_values, dtype=float),
        "income_transition": np.asarray(solution.income_transition, dtype=float),
        "child_transition": child_transition,
        "survival_probability_to_next_age": np.asarray(
            parameters.survival_probs, dtype=float
        ),
        "owner_housing_rungs": np.asarray(parameters.H_own, dtype=float),
        "date0_equilibrium_price": price,
        "raw_parameter_entry_by_location": np.asarray(
            parameters.entry_by_loc, dtype=float
        ),
        "effective_stationary_entry_by_location": np.asarray(
            effective_entry_by_location, dtype=float
        ),
        "tenure_codes": np.arange(1 + int(parameters.n_house), dtype=np.int64),
        "family_size_codes": np.arange(int(parameters.n_parity), dtype=np.int64),
        "child_status_codes": np.arange(int(parameters.n_child_states), dtype=np.int64),
    }
    np.savez_compressed(
        outdir / "initial_state_packet.npz",
        g0_stationary_pre_fertility=pre_fertility,
        g0_stationary_post_fertility_pre_housing=post_fertility,
        g0_stationary_realized_post_housing=current,
        entry_stationary_pre_fertility=entry_pre_fertility,
        g0_observed_age_pre_fertility=observed_pre,
        g0_observed_age_post_fertility_pre_housing=observed_post,
        g0_observed_age_realized_post_housing=observed_current,
        entry_observed_age_pre_fertility=observed_entry_pre_fertility,
        observed_age_shares=acs_age_shares,
        stationary_age_shares=np.sum(
            pre_fertility, axis=(0, 1, 2, 4, 5, 6)
        ),
        age_reweight_factors=observed_pre_factors,
        date0_stationary_housing_demand=stationary_housing_demand,
        date0_stationary_supply_intercept=np.asarray(parameters.H0, dtype=float),
        date0_observed_age_housing_demand=observed_housing_demand,
        date0_housing_stock_K0=observed_housing_demand,
        date0_observed_age_supply_intercept=observed_p0_clearing_supply_intercept,
        closure_reference_E0=np.array([effective_entry], dtype=float),
        closure_reference_B0=np.array(
            [float(solution.entrants_mature_total)], dtype=float
        ),
        closure_reference_M0=np.array(
            [effective_entry - float(solution.entrants_mature_total)], dtype=float
        ),
        baseline_one_step_next_pre_fertility_age_mass=observed_next_age_mass,
        **common_arrays,
    )
    write_age_summary(
        outdir / "stationary_age_state_summary.csv",
        pre_fertility,
        post_fertility,
        current,
        parameters,
    )
    write_age_summary(
        outdir / f"observed_age_{int(args.acs_year)}_state_summary.csv",
        observed_pre,
        observed_post,
        observed_current,
        parameters,
    )
    write_acs_age_bins(
        outdir / f"acs_{int(args.acs_year)}_head_age_bins.csv", acs_receipt
    )
    write_source_to_state(outdir / "source_to_state.csv")

    maturation_probability = float(
        child_transition[dependent_state, dependent_state + 1, 1]
    )
    metadata = {
        "status": "verified",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "purpose": "stationary and observed-age initial states for a calendar-time transition",
        "model": "current Stone--Geary one-shot completed-family-size model",
        "source": {
            "results_json": str(source_path),
            "results_sha256": sha256(source_path),
            "internal_profile": NEW_MOMENT_PROFILE_NAME,
            "selected_status": selected["status"],
            "theta": theta,
            "source_price": float(selected["price"]),
            "source_loss": float(selected["rank_loss"]),
        },
        "reproduction": {
            "price": float(price[0]),
            "loss": float(reproduced_loss),
            "market_residual": float(solution.best_max_abs_rel_excess),
            "strict_converged": bool(
                getattr(solution, "timings", {}).get("strict_converged", False)
            ),
            "overrides": overrides,
            "solver_sha256": sha256(
                ROOT / "code/model/intergen_housing_fertility_optimized/solver.py"
            ),
            "profile_sha256": sha256(
                ROOT
                / "code/model/intergen_housing_fertility_optimized/new_moment_profile.py"
            ),
        },
        "state_axes": [
            "liquid_wealth",
            "tenure_or_owner_rung",
            "location",
            "adult_age",
            "income_state",
            "one_shot_family_size_bin",
            "child_status",
        ],
        "state_shape": list(post_fertility.shape),
        "packet": {
            "path": str(outdir / "initial_state_packet.npz"),
            "primary_observed_key": "g0_observed_age_pre_fertility",
            "stationary_nesting_key": "g0_stationary_pre_fertility",
            "future_fixed_entry_reference_key": "entry_stationary_pre_fertility",
            "initial_youngest_head_stock_key": "entry_observed_age_pre_fertility",
            "date0_fixed_housing_stock_key": "date0_housing_stock_K0",
        },
        "timing": {
            "g0_pre_fertility": "start of calendar date, before the current-age family-size choice",
            "g0_post_fertility_pre_housing": "after family-size choice, before current location/tenure/housing transaction",
            "g0_realized_post_housing": "after current location/tenure/housing transaction, before saving and aging",
            "first_counterfactual_date": "apply date-0 policies to g0_pre_fertility; no four-year implementation delay is imposed",
        },
        "state_codes": {
            "tenure_or_owner_rung": {
                "0": "renter",
                **{
                    str(index + 1): f"owner housing rung {float(size):g}"
                    for index, size in enumerate(parameters.H_own)
                },
            },
            "one_shot_family_size_bin": {
                "0": "0 children",
                "1": "1--2 children empirical bin",
                "2": "3+ children empirical bin",
            },
            "child_status": {
                "0": "no dependent child / childless",
                "1": "dependent children on the shared child clock",
                "2": "family-size bin 1 matured / no dependent children",
                "3": "family-size bin 2 matured / no dependent children",
            },
        },
        "model_clock": {
            "period_years": float(parameters.period_years),
            "adult_ages": age_grid,
            "fertility_choice_ages": age_grid[
                int(parameters.A_f_start) - 1 : int(parameters.A_f_end)
            ],
            "child_stage_duration_periods": np.asarray(
                parameters.stage_durations, dtype=float
            ),
            "child_maturation_probability_per_period": maturation_probability,
            "implied_expected_child_stage_years": float(
                parameters.period_years / maturation_probability
            ),
            "child_clock_interpretation": (
                "one shared geometric clock for all children selected in the one-shot family-size draw; "
                "the maintained model has no separate child ages or birth-vintage states"
            ),
        },
        "stationary_flows": {
            "effective_entry_E": effective_entry,
            "raw_parameter_entry": raw_parameter_entry,
            "normalization_scale_effective_over_raw_entry": effective_entry
            / raw_parameter_entry,
            "terminal_exit": terminal_exit,
            "death_exit_before_terminal": death_exit,
            "mature_child_flow_B": float(solution.entrants_mature_total),
            "B_over_E": float(solution.entrants_mature_total) / effective_entry,
            "note": (
                "Use effective_entry_E, not raw_parameter_entry, to nest the maintained normalized "
                "stationary KFE. The gap is the solver's ex-post population normalization."
            ),
        },
        "child_clock_receipt": child_receipt,
        "observed_age_initializer": {
            "status": "diagnostic_empirical_age_stock",
            "acs": acs_receipt,
            "mapping": (
                "Each stationary age cell is multiplied by ACS_share_j/model_share_j. "
                "All conditional wealth, tenure, income, family-size, and child-status mass within age is unchanged."
            ),
            "age_reweight_factors_pre_fertility": observed_pre_factors,
            "age_reweight_factors_post_fertility": observed_post_factors,
            "age_reweight_factors_realized_current": observed_current_factors,
            "observed_youngest_age_cell_mass": float(
                np.sum(observed_entry_pre_fertility)
            ),
            "initial_youngest_stock_warning": (
                "entry_observed_age_pre_fertility is the date-0 age-18--21 head stock, not the future entry rule. "
                "Repeating it would conflate demographic entry with household formation because co-resident nonheads are outside the head sample. "
                "Use entry_stationary_pre_fertility for the agreed E0 fixed-entry reference."
            ),
            "baseline_one_step_next_total_mass": float(np.sum(observed_next_pre)),
            "baseline_one_step_next_age_mass": observed_next_age_mass,
            "implied_date0_mature_reproduction_flow": float(
                observed_child_receipt["matured_model_reproduction_units"]
            ),
        },
        "date0_market_normalization": {
            "common_price_p0": price,
            "stationary_housing_demand_at_p0": stationary_housing_demand,
            "observed_age_housing_demand_at_p0": observed_housing_demand,
            "fixed_housing_stock_K0": observed_housing_demand,
            "calibrated_stationary_supply_intercept": np.asarray(
                parameters.H0, dtype=float
            ),
            "observed_age_p0_clearing_supply_intercept": observed_p0_clearing_supply_intercept,
            "supply_elasticity": np.asarray(parameters.xi_supply, dtype=float),
            "convention": (
                "For every observed-G0 scenario, hold p0 fixed at the retained current equilibrium and use "
                "observed_age_p0_clearing_supply_intercept so date-0 housing stock K0 and market clearing are identical. "
                "Subsequent main runs keep that stock fixed; the static supply curve is diagnostic only."
            ),
        },
        "scenario_reference_flows": {
            "E0": effective_entry,
            "B0": float(solution.entrants_mature_total),
            "open_closure_M0_equals_E0_minus_B0": effective_entry
            - float(solution.entrants_mature_total),
            "contract": (
                "same observed-age G0 and date-0 housing stock in all scenarios; no fertility-preference shock; "
                "endogenous_closed M=0; endogenous_open M=E0-B0; fixed-entry benchmark entrants=E0"
            ),
        },
        "validation": checks,
        "transition_inputs_not_added": {
            "empirical_trends": "none",
            "shock_path": "none",
            "fertility_preference_shift": "none",
        },
    }
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "verified current one-shot initial states: "
        f"E={effective_entry:.15f}, B={float(solution.entrants_mature_total):.15f}, "
        f"B/E={float(solution.entrants_mature_total) / effective_entry:.12f}"
    )
    print(
        f"ACS {int(args.acs_year)} heads: n={int(acs_receipt['matched_mms_head_records'])}, "
        f"HHWT={float(acs_receipt['weighted_head_mass']):.0f}"
    )
    print(f"wrote {outdir / 'initial_state_packet.npz'}")
    print(f"wrote {outdir / 'metadata.json'}")


if __name__ == "__main__":
    main()
