"""Bounded closed stationary GE at the Phase-A-selected constant renter limit.

This module deliberately does not call the lab's default equilibrium command:
that command clears normalized housing demand. Here a price clears *actual*
birth renewal and population scales demand against the frozen physical supply
curve. Every trial is one fresh fixed-price lifecycle solve. A trial only counts
as an equilibrium after the native calendar, PAYGO, and reporting gates pass.
"""
from __future__ import annotations

import math
import time
import copy
import csv
import signal
from types import SimpleNamespace
from pathlib import Path

import numpy as np

from single_price import solve_fixed_price

RENEWAL_TOL = 1e-6
PAYGO_TOL = 1e-6
NATIVE_L1_TOL = 5e-9


def _require(ok, message):
    if not ok:
        raise RuntimeError(message)


def validate_parameter_estimates(context, rows, actual):
    """Only nine proposed coordinates, dimensions and one derived benefit vary."""
    expected = dict(context["expected_parameters"])
    _require(set(actual) == {r["parameter"] for r in rows}, "Actual parameter row set differs")
    for row in rows:
        key = row["parameter"]
        value, target = float(actual[key]), float(expected[key])
        # Annual beta is a floating-point period/annual round trip in the observer.
        ok = abs(value-target) <= 2e-12 if key == "beta_annual" else value == target
        _require(ok, "Actual proposed parameter differs: " + key)


def _observer_context(context):
    names = ("fp", "prepared", "manifest", "objective", "runtime", "reference")
    missing = [name for name in names if name not in context]
    _require(not missing, "Authenticated frozen observer missing: " + ", ".join(missing))
    return tuple(context[name] for name in names)


def _alarm(_signum, _frame):
    raise TimeoutError("300-second case deadline during native observer")


def _observe_with_deadline(context, live, label, *, final=False):
    end = float(live.get("case_deadline_epoch", time.time() + 300))
    remaining = min(end, float(context["deadline_epoch"])) - time.time()
    _require(remaining > 0, "Case or total deadline before native observer")
    old_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _alarm)
    old_timer = signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        return observe_price(context, live, label, final=final)
    finally:
        signal.setitimer(signal.ITIMER_REAL, *old_timer)
        signal.signal(signal.SIGALRM, old_handler)


def observe_price(context, live, name, *, final=False):
    """Observe the actual native distribution; no further lifecycle solve."""
    fp, prepared, manifest, objective, runtime, reference = _observer_context(context)
    cal = prepared.rt["primitive"].pf.calendar
    P, grid, sd, sol = (live[k] for k in ("P", "b_grid", "sd", "sol"))
    price = np.asarray(live["price"], dtype=float).reshape(-1)
    _require(price.size == 1 and np.isfinite(price[0]) and price[0] > 0, "Invalid scalar price")
    _require(float(P.unsecured_credit_limit) == float(context["selected_d_bar"]), "Credit limit drift")
    _require(float(P.psi_child) == float(context["reference_psi"]), "Child benefit drift")
    _require(np.array_equal(grid, context["b_grid"]), "Wealth grid drift")
    _require(all(float(getattr(P, k)) == float(getattr(context["P"], k)) for k in ("tau_pay", "psi")),
             "Fiscal or sale rule drift")
    _require(np.array_equal(P.H0, context["P"].H0) and np.array_equal(P.xi_supply, context["P"].xi_supply),
             "Physical housing supply curve drift")
    _require(np.array_equal(P.r_bar, context["P"].r_bar) and float(P.pension) == float(context["P"].pension),
             "Supply rent anchor or pension input drift")
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, recon = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(recon["stationary_post_fertility_nesting_l1"], NATIVE_L1_TOL,
                             "Stationary cohort reconstruction")
    runtime.require_abs_gate(recon["stationary_feasibility_projection_mass"], 0.,
                             "Stationary feasibility projection")
    supply = cal.HousingSupplyRule("static-elastic", float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]),
        float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
                             supply_rule=supply, supplied_policy=policy)
    renter_floor = audit_realized_renter_floor(P, grid, policy, ev)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=pre, supply_rule=supply,
                  demographic_seed=reference.get("demographic_seed"))
    out = Path(context["out"]) / "phase_b_ge" / name
    out.mkdir(parents=True, exist_ok=True)
    gates = fp.gates(packet, prepared, out, stationary=True)
    fiscal = gates["fiscal_certificate"]
    paygo = float(fiscal["actual_accounts"]["scaled_pension_budget_residual"])
    _require(abs(paygo) <= PAYGO_TOL, "Actual stationary PAYGO residual fails")
    native = prepared.rt["primitive"].pf.transition
    actual_births = native.calendar_topcode_birth_accounting(
        ev.g_pre, ev.g_post_fertility, float(ev.births), P)["topcode_adjusted_birth_children"]
    entry = float(sol.entry_rate)
    demand = float(np.asarray(ev.demand_by_loc).sum())
    physical_supply = float(np.asarray(ev.supply_by_loc).sum())
    _require(all(map(math.isfinite, (actual_births, entry, demand, physical_supply, paygo))),
             "Nonfinite GE accounting")
    _require(entry > 0 and demand > 0 and physical_supply > 0 and actual_births >= 0,
             "Invalid GE accounting")
    population = physical_supply / demand
    renewal = float(actual_births / (2.1 * entry) - 1)
    result = dict(price=float(price[0]), d_bar=float(P.unsecured_credit_limit),
                  adjusted_births_per_normalized_household=float(actual_births),
                  actual_entry_per_normalized_household=entry,
                  renewal_residual=renewal, population_scale=population,
                  normalized_housing_demand=demand,
                  physical_housing_supply=physical_supply,
                  absolute_housing_demand=population * demand,
                  absolute_housing_residual=population * demand - physical_supply,
                  actual_paygo_residual=paygo, outside_entry=0.0,
                  birth_to_entry_conversion=1 / 2.1,
                  housing_supply_elasticity=float(P.xi_supply[0]),
                  occupied_renter_floor=renter_floor)
    fp.write(out / "closure.json", result)
    if final:
        _require(abs(renewal) <= RENEWAL_TOL, "Birth renewal root fails")
        result["native_population_step"] = native_population_step(
            context, packet, population, actual_births, entry)
        fertility = {k: prepared.rt["observe_initial_fertility"](ev, P, age_projection=k)
                     for k in ("uniform_birth_time", "constant_post_cell")}
        housing = prepared.rt["observe_initial_housing_wealth"](
            ev, P, grid, sd, diagnostic_enabled=True,
            age_projection="uniform_within_age_cell", diagnostic_allow_family_proxies=True,
            include_wealth=True, include_birth_response=True)
        recent = prepared.rt["observe_recent_parent_flow"](
            ev, P, diagnostic_enabled=True, snapshot=prepared.rt["SNAPSHOT"],
            age_projection=prepared.rt["AGE_PROJECTION"], diagnostic_allow_residence_proxy=True,
            input_provenance=dict(case_id=name,
                                  reference_checkpoint_sha256=manifest["checkpoint"]["sha256"]))
        fits = runtime.score_targets(objective, fertility, housing, recent["model_value"],
                                     float(prepared.rt["chain"].extract_moments(sol, P)["tfr"]))
        params = [dict(row) for row in manifest["full_parameter_table"]]
        actual_params = fp.actual_parameters(prepared, P, grid)
        validate_parameter_estimates(context, params, actual_params)
        for row in params:
            row["reference_estimate"], row["estimate"] = row["estimate"], str(float(actual_params[row["parameter"]]))
            key = row["parameter"]
            if key in context["free_coordinates"]:
                row["status"] = "free in exploratory entry calibration"
                lo, hi = float(row["lower"]), float(row["upper"])
                value = float(row["estimate"])
                row["near_bound"] = str(min(value-lo,hi-value) <= .01*(hi-lo))
            elif key == "H0":
                row["status"] = "fixed population scale; not identified by per-household targets"
            elif key == "psi_child":
                row["status"] = "fixed reference benefit; price clears birth renewal"
            elif key == "child_benefit_CRRA_coefficient":
                row["status"] = "derived from fixed benefit and proposed curvature"
        _require(len(fits) == 14 and len(params) == 31, "14 fit/31 parameter rows required")
        fp.table(out / "target_fit.csv", fits)
        fp.table(out / "parameters.csv", params)
        fp.write(out / "observers.json", cal.jsonable(
            dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
        # Native plotting expects supply and demand in the same units. Report
        # S(q)/N alongside normalized demand, but retain absolute S(q) in the
        # solved packet, the fiscal/native gates, and the closure receipt.
        report_packet = dict(packet)
        report_ev = copy.copy(ev)
        report_ev.supply_by_loc = np.asarray(ev.supply_by_loc) / population
        report_packet["evaluation"] = report_ev
        prepared.rt["audit"].standard_diagnostics(report_packet, out, validate_production_young=False)
        plots = sorted(p.name for p in (out / "standard_diagnostics").glob("*.png"))
        _require(plots == sorted(manifest["standard_diagnostic_names"]), "Standard 17 plots required")
        result["standard_plot_count"] = len(plots)
        result["target_fit_rows"] = len(fits)
        result["parameter_rows"] = len(params)
        result["standard_plot_supply_units"] = "physical supply divided by endogenous household population"
        fp.write(out / "closure.json", result)
    return result


def audit_realized_renter_floor(P, grid, policy, ev):
    """Audit the occupied postchoice renter states to which the renter rule applies.

    The tenure index in ``g_current`` and ``bp_pol`` is the realized housing
    choice; a prechoice renter who purchases may legitimately hold mortgage
    debt and is therefore not counted as an unsecured renter violation.
    """
    from small_credit_lab.engine.household import renter_borrowing_floor
    wealth = np.asarray(grid, dtype=float)
    mass = np.asarray(ev.g_current[:, 0], dtype=float)
    saving = np.asarray(policy.bp_pol[:, 0], dtype=float)
    _require(mass.shape == saving.shape and mass.ndim == 6 and mass.shape[0] == wealth.size,
             "Renter state/policy shape differs")
    _require(np.isfinite(mass).all() and (mass >= 0).all(), "Invalid occupied renter mass")
    violation_mass, worst = 0., 0.
    for age in range(mass.shape[2]):
        floor = np.asarray(renter_borrowing_floor(P, wealth, age), dtype=float)
        _require(floor.shape == wealth.shape, "Native renter floor shape differs")
        gap = floor[:, None, None, None, None] - saving[:, :, age, :, :, :]
        occupied = mass[:, :, age, :, :, :] > 0
        _require(np.isfinite(gap[occupied]).all(), "Nonfinite occupied renter saving")
        violation_mass += float(mass[:, :, age, :, :, :][occupied & (gap > 1e-9)].sum())
        if np.any(occupied):
            worst = max(worst, float(np.max(gap[occupied])))
    _require(violation_mass <= 2e-10 and worst <= 1e-9,
             "Occupied renter saving violates constant credit or mortality floor")
    return dict(violation_mass=violation_mass, maximum_shortfall=worst,
                occupied_renter_mass=float(mass.sum()))


def native_population_step(context, packet, population, births, entry):
    """One native 16/20-year step on scaled mass, with no normalization."""
    import copy
    fp, prepared, *_ = _observer_context(context)
    native = prepared.rt["primitive"].pf.transition
    cal = prepared.rt["primitive"].pf.calendar
    ev = copy.copy(packet["evaluation"])
    for field in ("g_pre", "g_post_fertility", "g_current", "g_stay_distribution"):
        value = getattr(ev, field, None)
        if value is not None:
            setattr(ev, field, population * value)
    ev.births *= population
    ev.demand_by_loc = population * ev.demand_by_loc
    queue = native.SplitBirthEntryQueue.constant_prehistory(population * births)
    due, _ = queue.step(population * births)
    _require(len(queue.due_in_16) == 3 and len(queue.due_in_20) == 4,
             "Native 16/20-year entry timing differs")
    _require(abs(due - population * births / 2.1) <= max(1., population) * 2e-10,
             "Native birth-to-entry conversion differs")
    following, _, deaths, mass_res = native.advance_sequential_calendar_distribution(
        ev, np.asarray([due]), packet["parameters"], packet["b_grid"], packet["shared"])
    target = population * packet["stationary_g_pre"]
    raw = float(np.abs(following - target).sum())
    gap = due - population * entry
    corrected = following - target
    corrected[:, :, :, 0] -= gap * cal.entrant_cohort(
        np.asarray([1.]), packet["parameters"], packet["b_grid"])
    adjusted = float(np.abs(corrected).sum())
    _require(adjusted <= max(1., population) * NATIVE_L1_TOL,
             "Native stationary distribution step fails")
    _require(abs(float(mass_res)) <= max(1., population) * 2e-8,
             "Native population mass accounting fails")
    return dict(actual_birth_derived_entry=float(due), required_entry=population * entry,
                entry_gap=float(gap), raw_distribution_l1=raw,
                renewal_adjusted_distribution_l1=adjusted,
                mass_residual=float(mass_res), deaths=float(deaths))


def _bracket(points):
    ordered = sorted(points, key=lambda x: x["price"])
    for left, right in zip(ordered, ordered[1:]):
        if left["renewal_residual"] * right["renewal_residual"] <= 0:
            return left, right
    return None


def run_phase_b(context, phase_a_result, budget):
    """Adaptive diagnostic bracketing, then unchanged safeguarded root/repeat.

    The qref/8 and 8*qref caps are external numerical search limits, not
    restrictions on the economic equilibrium price. Every trial is fresh.
    """
    _observer_context(context)
    d_bar = float(phase_a_result["selected_d_bar"])
    _require(math.isfinite(d_bar) and d_bar >= 0, "Invalid nonnegative constant credit")
    context["selected_d_bar"] = d_bar
    context["reference_psi"] = float(context["P"].psi_child)
    context["deadline_epoch"] = float(budget.deadline_epoch)
    qref = float(np.asarray(context["q_ref"]).reshape(-1)[0])
    lower, upper = qref / 8, 8 * qref
    _require(lower > 0 and math.isfinite(upper), "Invalid reference price bracket")
    _require(int(budget.remaining_lifecycle) >= 2 and time.time() < budget.deadline_epoch,
             "No GE solve and repeat reserve")
    points = []
    phase_b_new = 0
    max_new = int(context.get("phase_b_max_new_lifecycle", 32))
    _require(2 <= max_new <= 32, "Invalid Phase B lifecycle cap")
    max_expansions = 20
    search = dict(method="adaptive_multiplicative_then_safeguarded_secant",
                  bounds_origin_qref=qref, lower=lower, upper=upper,
                  bounds_role="external numerical search caps; not economic restrictions",
                  expansion_factor=1.6, maximum_expansion_attempts=max_expansions,
                  maximum_new_lifecycle=max_new, attempts=[], events=[])
    phase_out = Path(context["out"]) / "phase_b_ge"
    phase_out.mkdir(parents=True, exist_ok=True)

    def record_search():
        context["fp"].write(phase_out / "price_search.json", search)

    def can_search():
        return (int(budget.remaining_lifecycle) >= 2 and phase_b_new < max_new - 1
                and time.time() + budget.stage_deadline_seconds + 400 < budget.deadline_epoch)

    def trial(q, label, *, repeat=False, reason):
        nonlocal phase_b_new
        _require(lower <= q <= upper and time.time() < budget.deadline_epoch,
                 "Price or total deadline exceeded")
        _require(phase_b_new < max_new, "Phase B new-lifecycle cap reached")
        _require(int(budget.remaining_lifecycle) >= (1 if repeat else 2),
                 "Exact repeat reserve would be consumed")
        if not repeat:
            _require(time.time() + budget.stage_deadline_seconds + 400 < budget.deadline_epoch,
                     "No time reserve for selected reporting and exact repeat")
        else:
            _require(time.time() + budget.stage_deadline_seconds < budget.deadline_epoch,
                     "No time for exact repeat")
        attempt = dict(label=label, price=float(q), reason=reason, status="started")
        search["attempts"].append(attempt)
        record_search()
        try:
            live = solve_fixed_price(context, d_bar, q, budget, label,
                                    Path(context["out"]) / "phase_b_ge" / label / "stage")
            phase_b_new += 1
            observed = _observe_with_deadline(context, live, label, final=False)
        except Exception as error:
            attempt.update(status="failed_fatal", error_type=type(error).__name__, error=str(error))
            record_search()
            raise
        attempt.update(status="observed", renewal_residual=float(observed["renewal_residual"]))
        record_search()
        points.append(observed)
        context["fp"].write(phase_out / "latest_completed.json",
                            dict(points=points, remaining_lifecycle=int(budget.remaining_lifecycle),
                                 deadline_epoch=budget.deadline_epoch, price_search=search))
        best = min(points, key=lambda row: abs(row["renewal_residual"]))
        context["fp"].write(phase_out / "best_so_far.json",
                            dict(status="price_trial_not_certified_GE", best=best,
                                 completed_trials=len(points), new_phase_b_lifecycle=phase_b_new))
        return observed, live

    requested_start = context.get("price_start", qref)
    try:
        start = float(requested_start)
    except (TypeError, ValueError):
        start = qref
    if not (math.isfinite(start) and lower <= start <= upper):
        start = qref
    search.update(requested_start=str(requested_start), start=start)
    seed, seed_live = trial(start, "price_start", reason="fresh initial price evaluation")
    if abs(seed["renewal_residual"]) <= RENEWAL_TOL:
        root, root_live = seed, seed_live
    else:
        root = root_live = None
        bracket = _bracket(points)
        first_direction = 1 if seed["renewal_residual"] > 0 else -1
        # A direction is a search ordering, never an acceptance assumption.
        # Defer at most once on an unexpected slope, checking the other side
        # before resuming. All acceptance still uses observed signs and gates.
        directions = [(first_direction, seed, False), (-first_direction, seed, False)]
        expansion = 0
        while root is None and bracket is None and directions and can_search() and expansion < max_expansions:
            direction, previous, resumed = directions.pop(0)
            while root is None and bracket is None and can_search() and expansion < max_expansions:
                proposal = min(upper, previous["price"] * 1.6) if direction > 0 else max(lower, previous["price"] / 1.6)
                if proposal == previous["price"]:
                    search["events"].append(dict(reason="diagnostic_price_cap", direction=direction, price=proposal))
                    record_search()
                    break
                expansion += 1
                item, live = trial(proposal, f"expand_{expansion:02d}",
                                   reason=f"multiplicative expansion direction={direction}; initial residual-directed direction={first_direction}; resumed={resumed}")
                bracket = _bracket(points)
                if abs(item["renewal_residual"]) <= RENEWAL_TOL:
                    root, root_live = item, live
                if root is not None or bracket is not None:
                    break
                slope = (item["renewal_residual"] - previous["renewal_residual"]) / (item["price"] - previous["price"])
                previous = item
                if slope >= 0 and not resumed:
                    search["events"].append(dict(reason="unexpected_nonnegative_slope_try_other_direction", direction=direction, price=item["price"], slope=slope))
                    directions.append((direction, previous, True))
                    record_search()
                    break
        if root is None and bracket is None:
            reason = "budget_or_repeat_reserve" if not can_search() else "finite_expansion_attempt_cap" if expansion >= max_expansions else "both_diagnostic_price_caps"
            search["termination_reason"] = reason
            record_search()
            return dict(status="uncomputed_price_unbracketed", selected_d_bar=d_bar,
                        points=points, price_search=search,
                        remaining_lifecycle=int(budget.remaining_lifecycle))
        iteration = 0
        while root is None and can_search():
            a, b = bracket
            fa, fb = a["renewal_residual"], b["renewal_residual"]
            proposal = (a["price"] * fb - b["price"] * fa) / (fb - fa) if fb != fa else .5 * (a["price"] + b["price"])
            # Keep the new price away from an endpoint where secant stagnates.
            width = b["price"] - a["price"]
            proposal = min(max(proposal, a["price"] + .1 * width), b["price"] - .1 * width)
            iteration += 1
            item, live = trial(proposal, f"root_{iteration:02d}", reason="safeguarded secant/bisection inside observed sign bracket")
            if abs(item["renewal_residual"]) <= RENEWAL_TOL:
                root, root_live = item, live
            else:
                bracket = _bracket(points)
                _require(bracket is not None, "Renewal bracket lost")
    if root is None or time.time() + 400 >= budget.deadline_epoch:
        return dict(status="uncomputed_bounded_budget", selected_d_bar=d_bar,
                    points=points, price_search=search, remaining_lifecycle=int(budget.remaining_lifecycle))
    # Final reporting and the native population step are zero-lifecycle checks.
    certified = _observe_with_deadline(context, root_live, "selected_root", final=True)
    repeat, repeat_live = trial(root["price"], "selected_repeat", repeat=True, reason="fresh exact selected-price repeat")
    repeated = _observe_with_deadline(context, repeat_live, "selected_repeat_final", final=True)
    _require(repeated["renewal_residual"] == certified["renewal_residual"],
             "Exact selected repeat differs in renewal")
    _require(repeated["population_scale"] == certified["population_scale"],
             "Exact selected repeat differs in population")
    _require(set(vars(root_live["sol"])) == set(vars(repeat_live["sol"])),
             "Exact selected repeat solution field set differs")
    for key, value in vars(root_live["sol"]).items():
        twin = getattr(repeat_live["sol"], key, None)
        if isinstance(value, np.ndarray) and value.dtype != object:
            _require(isinstance(twin, np.ndarray) and np.isfinite(value).all()
                     and np.isfinite(twin).all() and np.array_equal(value, twin),
                     "Exact selected repeat differs in solution array: " + key)
    _require(set(vars(root_live["sd"])) == set(vars(repeat_live["sd"])),
             "Exact selected repeat shared field set differs")
    for key, value in vars(root_live["sd"]).items():
        twin = getattr(repeat_live["sd"], key, None)
        if isinstance(value, np.ndarray) and value.dtype != object:
            _require(isinstance(twin, np.ndarray) and np.isfinite(value).all()
                     and np.isfinite(twin).all() and np.array_equal(value, twin),
                     "Exact selected repeat differs in shared array: " + key)
    a = Path(context["out"]) / "phase_b_ge" / "selected_root"
    b = Path(context["out"]) / "phase_b_ge" / "selected_repeat_final"
    for table, expected in (("target_fit.csv", 14), ("parameters.csv", 31)):
        with (a / table).open(newline="") as fa, (b / table).open(newline="") as fb:
            aa, bb = list(csv.DictReader(fa)), list(csv.DictReader(fb))
        _require(len(aa) == len(bb) == expected and aa == bb, "Selected repeat table differs: " + table)
    result = dict(status="passed", selected_d_bar=d_bar, selected_price=root["price"],
                  selected=certified, exact_repeat=repeated, points=points,
                  price_search=search,
                  remaining_lifecycle=int(budget.remaining_lifecycle))
    context["fp"].write(Path(context["out"]) / "phase_b_ge" / "selected.json", result)
    return result


def smoke_phase_b(context, budget):
    """Exact price-loop shape with mock solves and zero lifecycle claims."""
    global solve_fixed_price, observe_price
    _observer_context(context)
    prior_solve, prior_observe = solve_fixed_price, observe_price
    qref = float(context["q_ref"])
    labels = []
    before = int(budget.remaining_lifecycle)
    fake = dict(context)
    fake["price_start"] = qref

    def stage(q):
        return dict(price=np.asarray([q]), sol=SimpleNamespace(a=np.asarray([q])),
                    sd=SimpleNamespace(b=np.asarray([1.])),
                    case_deadline_epoch=time.time() + 300)

    def solve(_context, d_bar, q, _budget, label, _stage_dir):
        _require(float(d_bar) == float(context["selected_d_bar"]), "Mock effective credit drift")
        labels.append(label)
        return stage(q)

    def observe(_context, live, label, *, final=False):
        q = float(live["price"][0])
        if final:
            directory = Path(fake["out"]) / "phase_b_ge" / label
            directory.mkdir(parents=True, exist_ok=True)
            for filename, count, field in (("target_fit.csv", 14, "moment"),
                                           ("parameters.csv", 31, "parameter")):
                with (directory / filename).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=[field, "value"])
                    writer.writeheader()
                    writer.writerows({field: str(i), "value": "1"} for i in range(count))
        return dict(price=q, renewal_residual=(1.02 - q / qref) * .01,
                    population_scale=1.03)

    try:
        solve_fixed_price, observe_price = solve, observe
        result = run_phase_b(fake, dict(selected_d_bar=float(context["selected_d_bar"]), selected_live=stage(qref)), budget)
    finally:
        solve_fixed_price, observe_price = prior_solve, prior_observe
    _require(result["status"] == "passed" and before == budget.remaining_lifecycle,
             "Mock GE loop or zero-lifecycle contract failed")
    _require(labels == ["price_start", "expand_01", "root_01", "root_02", "selected_repeat"],
             "Mock GE loop differed: " + repr(labels))
    return dict(status="passed_mock_zero_lifecycle", selected_d_bar=float(context["selected_d_bar"]), mock_calls=labels,
                selected_price_factor=result["selected_price"] / qref,
                lifecycle_solves=0)
