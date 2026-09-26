"""Pure, finite search design for the experimental four-arm utility comparison.

This module imports only the standard library and neither reads files nor runs a
model.  Generate/validate plans on Torch.  A plan is preparation, not permission
to launch or evidence that the proposed share-utility normalization is adopted.
The runtime must separately pin source, observation, entry, and utility-unit
contracts and enforce deadlines, gates, and original-versus-both repeat checks.

The initial population includes the supplied seeds.  Each later DE generation
has exactly that many objective attempts; a rejected attempt consumes its slot.
All four arms get the same counts, workers, bounds on common coordinates, and
named random streams.  Arm-specific survivor selection may subsequently diverge.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import random
from typing import Any, Mapping, Sequence


COMMON_NAMES = ("H0", "beta_annual", "chi", "first_birth_fixed_cost",
                "kappa_fert", "kappa_fert_continuation", "theta0")
FLOOR_NAMES = COMMON_NAMES + ("h_P",)
SHARE_NAMES = COMMON_NAMES + ("delta_alpha_jump", "delta_alpha")
ARM_NAMES = ("floor_linear", "floor_concave", "shares_linear", "shares_concave")
TARGET_IDS = (
    "initial_normalization", "cps_childlessness", "cps_exactly_one",
    "nchs_mean_age", "nchs_share30", "wealth_earnings", "bequest_wealth",
    "old_dispersion", "mean_rooms", "ownership_30_55", "first_birth_rooms",
    "family_rooms", "recent_parent_ownership",
)
# These are the approved full intervals in the September 25 frozen contract,
# not the narrower support of yesterday's selected/local proposal bank.
_FLOOR_RESTRICTIONS = {
    "H0": (.2, 80., "log"), "beta_annual": (.94, .99, "discount"),
    "chi": (.1, 5., "log"), "first_birth_fixed_cost": (0., 8., "softzero"),
    "kappa_fert": (.02, 50., "log"),
    "kappa_fert_continuation": (.02, 50., "log"),
    "theta0": (0., 8., "softzero"), "h_P": (.1, 2.3, "log"),
}
# One initial solve + at most eight bracketing + fourteen refinement calls.
NORMALIZATION_MAX_STATIONARY_CALLS = 23


def _number(value: Any, name: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise ValueError(f"{name} must be finite and >= {minimum}")
    return result


def _count(value: Any, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _positive(value: Any, name: str) -> float:
    result = _number(value, name, minimum=0.)
    if result == 0.:
        raise ValueError(f"{name} must be positive")
    return result


def canonical_fingerprint(value: Any) -> str:
    """Fingerprint supplied JSON content; this function does no filesystem IO."""
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _rows_by_name(rows: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    result = {}
    for row in rows:
        name = row["parameter"]
        if not isinstance(name, str) or not name or name in result:
            raise ValueError("Parameter names must be unique nonempty strings")
        lo, hi = _number(row["lower"], name), _number(row["upper"], name)
        kind = row["transform"]
        if lo >= hi or kind not in {"log", "discount", "softzero", "linear"}:
            raise ValueError(f"Invalid interval or transform for {name}")
        if kind == "log" and lo <= 0.:
            raise ValueError(f"Log lower bound must be positive for {name}")
        result[name] = {"parameter": name, "lower": lo, "upper": hi, "transform": kind}
    return result


def define_arms(floor_restrictions: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    """Validate supplied frozen restrictions; replace only h_P in share arms.

    The share transforms reproduce run_e1_chain.transform: x=.25*u**2.
    The exponent .86 is a fixed sensitivity, not an externally estimated fact.
    This catalog deliberately does not choose a share utility-unit adjustment.
    """
    floor = _rows_by_name(floor_restrictions)
    if set(floor) != set(FLOOR_NAMES):
        raise ValueError("Expected exactly seven common parameters plus h_P")
    for name, expected in _FLOOR_RESTRICTIONS.items():
        row = floor[name]
        if (row["lower"], row["upper"], row["transform"]) != expected:
            raise ValueError(f"Frozen full restriction changed for {name}")
    shares = {name: floor[name] for name in COMMON_NAMES}
    for name in SHARE_NAMES[len(COMMON_NAMES):]:
        shares[name] = {"parameter": name, "lower": 0., "upper": .25, "transform": "softzero"}
    result = {}
    for arm in ARM_NAMES:
        housing, curvature = arm.split("_")
        names = FLOOR_NAMES if housing == "floor" else SHARE_NAMES
        domain = floor if housing == "floor" else shares
        result[arm] = {
            "arm": arm, "housing": housing,
            "child_reward_exponent": 1. if curvature == "linear" else .86,
            "child_reward_state": "m: current children at home",
            "curvature_status": "linear reference" if curvature == "linear" else "fixed experimental sensitivity",
            "free_parameter_names": list(names), "free_parameter_count": len(names),
            "parameter_restrictions": [copy.deepcopy(domain[name]) for name in names],
            "weighted_target_count": 12, "displayed_target_count": 13,
            "psi_child": "derived separately from completed-fertility normalization 2.1",
            "theta1": "externally fixed; excluded from search",
            "utility_unit_contract": "must be supplied and reviewed separately before launch",
        }
    return result


def unit_to_physical(unit: float, row: Mapping[str, Any]) -> float:
    """Map strict [0,1] coordinates using the frozen transform conventions."""
    u = _number(unit, "unit", minimum=0.)
    if u > 1.:
        raise ValueError("unit must lie in [0,1]")
    checked = _rows_by_name([row])[row["parameter"]]
    lo, hi, kind = checked["lower"], checked["upper"], checked["transform"]
    if u == 0. or u == 1.:
        return lo if u == 0. else hi
    if kind == "log":
        return lo * (hi / lo) ** u
    if kind == "discount":
        return lo + (hi - lo) * (1. - (1. - u) ** 2)
    if kind == "softzero":
        return lo + (hi - lo) * u * u
    return lo + (hi - lo) * u


def physical_to_unit(value: float, row: Mapping[str, Any]) -> float:
    """Strict inverse; never silently clips an out-of-bounds supplied seed."""
    x = _number(value, row["parameter"])
    checked = _rows_by_name([row])[row["parameter"]]
    lo, hi, kind = checked["lower"], checked["upper"], checked["transform"]
    if not lo <= x <= hi:
        raise ValueError(f"{row['parameter']} outside [{lo}, {hi}]")
    if x == lo or x == hi:
        return 0. if x == lo else 1.
    if kind == "log":
        return math.log(x / lo) / math.log(hi / lo)
    q = (x - lo) / (hi - lo)
    if kind == "discount":
        return 1. - math.sqrt(1. - q)
    if kind == "softzero":
        return math.sqrt(q)
    return q


def size_search_budget(*, objective_solve_p90_seconds: float,
                       objective_overhead_seconds: float,
                       workers_per_arm: int, total_limit_seconds: float,
                       repeat_reserve_seconds: float, export_reserve_seconds: float,
                       initial_population: int, de_generations: int | None,
                       objective_timeout_seconds: float) -> dict[str, Any]:
    """Size equal finite arms from explicit timing/concurrency assumptions.

    Four arms run concurrently.  Within each arm the initial bank and every DE
    generation are separate barriers, each requiring ceil(population/workers)
    waves.  Two identical smokes run sequentially within each arm, exercising
    actual objective-loop progression.  Both selected repeats can run concurrently
    when at least two workers are available.  The runtime must use that scheduling
    or add enough reserve.  No observed timing or worker cap is a
    default.  P90 solve-sum plus *assumed* overhead is a planning estimate, never a
    completion guarantee.  A fixed deadline can leave the finite bank unfinished.
    """
    solve_p90 = _positive(objective_solve_p90_seconds, "objective_solve_p90_seconds")
    overhead = _number(objective_overhead_seconds, "objective_overhead_seconds", minimum=0.)
    workers = _count(workers_per_arm, "workers_per_arm", minimum=1)
    total = _positive(total_limit_seconds, "total_limit_seconds")
    reserve = _positive(repeat_reserve_seconds, "repeat_reserve_seconds")
    export = _positive(export_reserve_seconds, "export_reserve_seconds")
    population = _count(initial_population, "initial_population", minimum=1)
    timeout = _positive(objective_timeout_seconds, "objective_timeout_seconds")
    estimate = solve_p90 + overhead
    if timeout < estimate:
        raise ValueError("Objective timeout is below the supplied planning estimate")
    repeat_waves = math.ceil(2 / workers)
    repeat_estimate = repeat_waves * estimate
    if reserve < repeat_estimate:
        raise ValueError("Repeat reserve cannot accommodate both selected repeats at supplied timing")
    smoke_estimate = 2 * estimate
    search_seconds = total - smoke_estimate - reserve - export
    waves_per_generation = math.ceil(population / workers)
    available_waves = max(0, math.floor(search_seconds / estimate))
    max_generations = available_waves // waves_per_generation - 1
    if max_generations < 0:
        raise ValueError("Total limit cannot accommodate the complete initial population and reserves")
    generations = max_generations if de_generations is None else _count(de_generations, "de_generations")
    if generations > max_generations:
        raise ValueError("Requested finite DE generations exceed the supplied time budget")
    if generations and population < 4:
        raise ValueError("DE/rand/1 requires at least four population members")
    population_evaluations = population * (1 + generations)
    # The scientifically verified first smoke is initial slot zero. Reusing it
    # avoids a third unchanged objective while retaining the complete bank.
    search_attempts = population_evaluations - 1
    waves = waves_per_generation * (1 + generations)
    planned = smoke_estimate + waves * estimate + reserve + export
    return {
        "schema": "e5f_four_arm_finite_budget_v1", "status": "preparation_only_no_launch",
        "arm_count": 4, "workers_per_arm": workers, "total_workers": 4 * workers,
        "objective_solve_p90_seconds": solve_p90, "objective_overhead_seconds": overhead,
        "overhead_status": "explicit assumed allowance; excluded from observed solve-sum timing",
        "objective_planning_seconds": estimate, "objective_timeout_seconds": timeout,
        "total_limit_seconds": total, "repeat_reserve_seconds": reserve,
        "export_reserve_seconds": export, "smoke_estimate_seconds": smoke_estimate,
        "initial_population": population, "de_generations": generations,
        "max_de_generations_at_supplied_timing": max_generations,
        "waves_per_generation": waves_per_generation, "search_waves": waves,
        "search_population_evaluations_per_arm": population_evaluations,
        "reused_smoke_seed_cases_per_arm": 1,
        "search_attempts_per_arm": search_attempts, "identical_smokes_per_arm": 2,
        "smoke_schedule": "sequential within each arm; four arms concurrent",
        "selected_repeats_per_arm": 2, "total_attempts_per_arm": search_attempts + 4,
        "total_objective_attempts": 4 * (search_attempts + 4),
        "stationary_calls_per_objective_upper_bound": NORMALIZATION_MAX_STATIONARY_CALLS,
        "stationary_calls_upper_bound": 4 * (search_attempts + 4) * NORMALIZATION_MAX_STATIONARY_CALLS,
        "planned_seconds": planned, "unallocated_seconds": total - planned,
        "objective_timeouts_wave_bound_seconds": (waves + 2 + repeat_waves) * timeout + export,
        "timing_is_completion_guarantee": False,
        "deadline_rule": "stop without retries, replacement slots, or gate relaxation; report unfinished counts",
        "repeat_rule": "compare original selected result separately against each of both repeats",
    }


def validate_target_contract(target_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Preserve and fingerprint complete supplied rows; no new weight choices.

    Values/source agreement with a pinned adopted contract is the caller's
    responsibility.  This structural check does not certify empirical provenance
    or a full-rank identifying Jacobian merely from the count of moments.
    """
    rows = copy.deepcopy(list(target_rows))
    if len(rows) != 13 or {row["id"] for row in rows} != set(TARGET_IDS):
        raise ValueError("Expected all thirteen distinct active target rows")
    for row in rows:
        _number(row["target"], row["id"])
        for field in ("sample", "definition", "model_observation", "uncertainty_status", "mapping_warning"):
            if not row.get(field):
                raise ValueError(f"Missing {field} for {row['id']}")
        source = row.get("source", {})
        for field in ("path", "builder", "record_id", "contract_id"):
            if not source.get(field):
                raise ValueError(f"Missing source {field} for {row['id']}")
        if row["id"] == "initial_normalization":
            if row["target"] != 2.1 or row.get("weight") is not None:
                raise ValueError("psi_child requires the separate unscored normalization 2.1")
        else:
            _positive(row.get("weight"), f"{row['id']} weight")
    return {
        "displayed_rows": 13, "positive_weight_rows": 12,
        "target_rows": rows, "target_rows_sha256": canonical_fingerprint(rows),
        "normalization": {"parameter": "psi_child", "target": 2.1, "weight": None,
                          "status": "externally fixed normalization; parameter derived separately"},
        "free_parameter_counts": {"floor": 8, "shares": 9},
        "identification": "moment count only; no numerical rank or identification claim",
    }


def validate_budget(budget: Mapping[str, Any]) -> None:
    """Reject altered/inconsistent serialized finite counts before bank or DE."""
    keys = ("objective_solve_p90_seconds", "objective_overhead_seconds", "workers_per_arm",
            "total_limit_seconds", "repeat_reserve_seconds", "export_reserve_seconds",
            "initial_population", "de_generations", "objective_timeout_seconds")
    if dict(budget) != size_search_budget(**{key: budget[key] for key in keys}):
        raise ValueError("Budget counts or timing receipts differ from the supplied finite design")


def _rng(master_seed: int, *keys: Any) -> random.Random:
    _count(master_seed, "rng_seed")
    return random.Random(int(canonical_fingerprint([master_seed, *keys]), 16))


def _reflect(value: float) -> float:
    remainder = value % 2.
    return remainder if remainder <= 1. else 2. - remainder


def _domain(arms: Mapping[str, Any], arm: str) -> dict[str, dict[str, Any]]:
    if set(arms) != set(ARM_NAMES):
        raise ValueError("Expected all four arms")
    # Revalidate the common/floor bounds and the catalog before every bank/DE.
    expected = define_arms(arms["floor_linear"]["parameter_restrictions"])
    if dict(arms) != expected:
        raise ValueError("Arm catalog differs from the validated common search contract")
    if arm not in expected:
        raise ValueError(f"Unknown arm {arm!r}")
    return _rows_by_name(expected[arm]["parameter_restrictions"])


def _candidate(arm: str, identifier: str, unit: Mapping[str, float],
               domain: Mapping[str, Any], **metadata: Any) -> dict[str, Any]:
    if set(unit) != set(domain):
        raise ValueError("Candidate coordinate set differs from arm domain")
    return {"id": identifier, "arm": arm, "unit": dict(unit),
            "parameters": {name: unit_to_physical(unit[name], row) for name, row in domain.items()},
            **metadata}


def build_proposal_bank(arms: Mapping[str, Any], *,
                        shared_seeds: Mapping[str, Mapping[str, float]],
                        floor_seeds: Mapping[str, Mapping[str, float]],
                        share_seeds: Mapping[str, Mapping[str, float]],
                        broad_count: int, medium_count: int, local_count: int,
                        medium_unit_scale: float, local_unit_scale: float,
                        rng_seed: int, budget: Mapping[str, Any]) -> dict[str, Any]:
    """Initial seeds plus full-dimensional uniform/medium/local proposals.

    Seed-label sets must agree and all seven common coordinates must be supplied.
    Preference seeds are explicit mappings: {label: {h_P: ...}} and {label:
    {delta_alpha_jump: ..., delta_alpha: ...}}.  Medium/local proposals add an
    independent Gaussian to *every* unit coordinate, reflecting into [0,1]; scale
    means Gaussian standard deviation in that unit interval.  Broad proposals
    independently cover the entire transformed unit cube.  No rejection retries.
    """
    validate_budget(budget)
    if not shared_seeds or any(not isinstance(k, str) or not k for k in shared_seeds):
        raise ValueError("At least one named common seed is required")
    labels = sorted(shared_seeds)
    if set(labels) != set(floor_seeds) or set(labels) != set(share_seeds):
        raise ValueError("All supplied common/preference seed label sets must agree")
    counts = {"broad": _count(broad_count, "broad_count", minimum=1),
              "medium": _count(medium_count, "medium_count", minimum=1),
              "local": _count(local_count, "local_count", minimum=1)}
    medium = _positive(medium_unit_scale, "medium_unit_scale")
    local = _positive(local_unit_scale, "local_unit_scale")
    if not 0. < local < medium <= 1.:
        raise ValueError("Require 0 < local_unit_scale < medium_unit_scale <= 1")
    _count(rng_seed, "rng_seed")
    population = len(labels) + sum(counts.values())
    if population != budget["initial_population"]:
        raise ValueError("Seeds plus proposals must equal the budgeted initial population")
    output = {}
    for arm in ARM_NAMES:
        domain = _domain(arms, arm)
        preference = floor_seeds if arm.startswith("floor_") else share_seeds
        extra_names = {"h_P"} if arm.startswith("floor_") else {"delta_alpha_jump", "delta_alpha"}
        centers, seeds = {}, []
        for index, label in enumerate(labels):
            if set(shared_seeds[label]) != set(COMMON_NAMES) or set(preference[label]) != extra_names:
                raise ValueError(f"Seed {label!r} has missing, unknown, or externally fixed coordinates")
            point = {**shared_seeds[label], **preference[label]}
            centers[label] = {name: physical_to_unit(point[name], row) for name, row in domain.items()}
            proposal = _candidate(arm, f"{arm}_initial_{index:04d}", centers[label], domain,
                                  stage="seed", center_seed=label, generation=0,
                                  selection_generation=0, slot=index)
            # Preserve supplied physical seeds exactly; the inverse/forward path
            # may otherwise introduce last-bit differences on original repeats.
            proposal["parameters"] = {name: _number(point[name], name) for name in domain}
            seeds.append(proposal)
        proposals = list(seeds)
        for stage, count in counts.items():
            for within_stage in range(count):
                index = len(proposals)
                label = labels[within_stage % len(labels)]
                unit = {}
                for name in domain:
                    stream = _rng(rng_seed, "initial", index, name)
                    unit[name] = stream.random() if stage == "broad" else _reflect(
                        centers[label][name] + stream.gauss(0., medium if stage == "medium" else local))
                proposals.append(_candidate(arm, f"{arm}_initial_{index:04d}", unit, domain,
                    stage=stage, center_seed=None if stage == "broad" else label, generation=0,
                    selection_generation=0, slot=index))
        output[arm] = proposals
    return {
        "schema": "e5f_four_arm_transformed_proposal_bank_v1", "rng_seed": rng_seed,
        "initial_population": population, "seed_labels": labels, "counts": counts,
        "scales": {"broad": "independent U[0,1] in every transformed coordinate",
                   "medium": medium, "local": local, "scale_unit": "unit-coordinate Gaussian standard deviation",
                   "boundary_rule": "reflection, with no rejection loop"},
        "random_stream_rule": "SHA256 of JSON [seed, stage, index, parameter]; common names shared across arms",
        "common_parameter_names": list(COMMON_NAMES), "arms": output,
        "budget_sha256": canonical_fingerprint(budget), "arm_catalog_sha256": canonical_fingerprint(arms),
    }


def _completed_scores(population: Sequence[Mapping[str, Any]],
                      scores: Mapping[str, float | None]) -> dict[str, float]:
    identifiers = [row["id"] for row in population]
    if len(set(identifiers)) != len(identifiers) or set(scores) != set(identifiers):
        raise ValueError("Exactly one completed result per unique population member is required")
    # None means a completed rejected objective, not a missing or unrun attempt.
    return {key: math.inf if value is None else _number(value, key, minimum=0.)
            for key, value in scores.items()}


def make_de_generation(arms: Mapping[str, Any], arm: str,
                       population: Sequence[Mapping[str, Any]],
                       completed_scores: Mapping[str, float | None], *,
                       generation: int, budget: Mapping[str, Any], rng_seed: int,
                       mutation_factor: float) -> list[dict[str, Any]]:
    """One finite DE/rand/1 generation, all-coordinate crossover (CR=1).

    Require a complete preceding barrier.  Rejected parents have infinite loss
    but remain population members; their new trials consume ordinary generation
    slots, never retries.  Stop if no parent is valid.  Donor-index streams are
    common across arms; selection depends on each arm's own completed objective.
    """
    validate_budget(budget)
    domain = _domain(arms, arm)
    generation = _count(generation, "generation", minimum=1)
    if generation > budget["de_generations"] or len(population) != budget["initial_population"]:
        raise ValueError("Generation or population exceeds the frozen finite budget")
    if len(population) < 4:
        raise ValueError("DE/rand/1 requires at least four population members")
    factor = _positive(mutation_factor, "mutation_factor")
    if factor > 2.:
        raise ValueError("mutation_factor must not exceed 2")
    losses = _completed_scores(population, completed_scores)
    if all(math.isinf(value) for value in losses.values()):
        raise ValueError("No valid completed parent; stop without replacement search")
    for index, row in enumerate(population):
        if row["arm"] != arm or row["slot"] != index or set(row["unit"]) != set(domain):
            raise ValueError("Population arm, slot, or coordinate contract changed")
        if row.get("generation", 0) >= generation or row.get("selection_generation") != generation - 1:
            raise ValueError("Parents must come from the completed immediately preceding selection barrier")
        if set(row["parameters"]) != set(domain):
            raise ValueError("Physical parameter set differs from the unit coordinates")
        for name, restriction in domain.items():
            predicted = unit_to_physical(row["unit"][name], restriction)
            physical = _number(row["parameters"][name], name)
            if not math.isclose(physical, predicted, rel_tol=2e-14, abs_tol=1e-14):
                raise ValueError("Physical point differs from its transformed coordinate")
    trials = []
    for index, parent in enumerate(population):
        donors = _rng(rng_seed, "de_donors", generation, index).sample(
            [j for j in range(len(population)) if j != index], 3)
        a, b, c = (population[j]["unit"] for j in donors)
        unit = {name: _reflect(a[name] + factor * (b[name] - c[name])) for name in domain}
        trials.append(_candidate(arm, f"{arm}_de{generation:03d}_{index:04d}", unit, domain,
            stage="de_rand_1", generation=generation, selection_generation=generation,
            slot=index, parent_id=parent["id"],
            donor_slots=donors, mutation_factor=factor, crossover_probability=1.))
    return trials


def select_de_generation(population: Sequence[Mapping[str, Any]],
                         parent_scores: Mapping[str, float | None],
                         trials: Sequence[Mapping[str, Any]],
                         trial_scores: Mapping[str, float | None]) -> tuple[list[dict[str, Any]], dict[str, float | None]]:
    """Select only strict improvements after the complete trial barrier.

    Ties and failed trials retain parents.  Returns independent JSON-ready copies
    and their completed scores; no missing receipt is interpreted as rejection.
    """
    old, new = _completed_scores(population, parent_scores), _completed_scores(trials, trial_scores)
    if len(population) != len(trials):
        raise ValueError("A DE generation must contain exactly one trial per parent")
    if not trials or len({row["generation"] for row in trials}) != 1:
        raise ValueError("Require one nonempty, complete generation")
    generation = _count(trials[0]["generation"], "generation", minimum=1)
    survivors, scores = [], {}
    for index, (parent, trial) in enumerate(zip(population, trials)):
        if (trial.get("parent_id") != parent["id"] or trial["arm"] != parent["arm"]
                or parent["slot"] != index or trial["slot"] != index
                or parent["selection_generation"] != generation - 1):
            raise ValueError("DE trial-to-parent pairing changed")
        accept = new[trial["id"]] < old[parent["id"]]
        winner = trial if accept else parent
        survivors.append(copy.deepcopy(winner))
        survivors[-1]["selection_generation"] = generation
        scores[winner["id"]] = trial_scores[winner["id"]] if accept else parent_scores[winner["id"]]
    return survivors, scores
