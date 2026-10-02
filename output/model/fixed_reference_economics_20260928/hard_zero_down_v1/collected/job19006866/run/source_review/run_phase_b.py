def run_phase_b(context, phase_a_result, budget):
    """One fresh lifecycle solve at the pinned price; no renewal or housing root."""
    _observer_context(context)
    d_bar = float(phase_a_result["selected_d_bar"])
    _require(d_bar >= 0 and math.isfinite(d_bar), "Invalid credit input")
    context["selected_d_bar"] = d_bar
    context["reference_psi"] = float(context["P"].psi_child)
    context["deadline_epoch"] = float(budget.deadline_epoch)
    q = float(context["fixed_price"])
    _require(q > 0 and math.isfinite(q), "Invalid fixed price")
    live = solve_fixed_price(context, d_bar, q, budget, "selected_root",
        Path(context["out"]) / "phase_b_ge" / "selected_root" / "stage")
    observed = _observe_with_deadline(context, live, "selected_root", final=True)
    return dict(status="passed_fixed_price_diagnostic", selected_price=q,
                selected=observed, lifecycle_solves=int(budget.used_lifecycle))
