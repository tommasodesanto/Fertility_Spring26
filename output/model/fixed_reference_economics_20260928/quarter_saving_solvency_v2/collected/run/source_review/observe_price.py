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
    population = 1.0  # Fixed physical population, not housing-clearing scale.
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
        result["native_population_step"] = native_population_step(
            context, packet, population, actual_births, entry)
        result["renewal_status"] = "measured_fixed_price_diagnostic"
        result["housing_market_status"] = "measured_fixed_H0_diagnostic"
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
            if key in context["fixed_coordinates"]:
                row["status"] = "fixed at strict-80 reference estimate for diagnostic"
            elif key in context["free_coordinates"]:
                row["status"] = "free in exploratory entry calibration"
                lo, hi = float(row["lower"]), float(row["upper"])
                value = float(row["estimate"])
                row["near_bound"] = str(min(value-lo,hi-value) <= .01*(hi-lo))
            elif key == "H0":
                row["status"] = "fixed strict-80 physical supply coefficient for diagnostic"
            elif key == "financed_share":
                row["status"] = "experimental fixed policy input"
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
        result["standard_plot_supply_units"] = "physical supply and demand at fixed N0=1"
        fp.write(out / "closure.json", result)
    return result
