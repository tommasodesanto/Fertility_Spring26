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
