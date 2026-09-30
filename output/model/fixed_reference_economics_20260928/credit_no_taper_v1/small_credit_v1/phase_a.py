"""Budget-derived positive unsecured credit, followed by one native PE replay.

The bound is necessary at entry, not a certificate of lifetime feasibility.
No entry distribution, fiscal transfer, earnings, preference, or price is changed.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def entrant_budget_bound(context, price_factors=(0.85, 1.0, 1.15)):
    from small_credit_lab.engine import solver
    from small_credit_lab.engine.shared import income_at_state

    P, grid = context["P"], context["b_grid"]
    if not (bool(getattr(P, "exhaustive_saving_control", False))
            and bool(getattr(P, "native_exact_allocation_output", False))
            and str(getattr(P, "interp_method", "linear")) == "linear"
            and float(getattr(P, "rental_wedge_intercept", 0.0)) == 0.0
            and float(getattr(P, "rental_wedge_slope", 0.0)) == 0.0):
        raise RuntimeError("Bound derived only for authenticated exhaustive exact non-wedge renter branch")
    sd = solver.precompute_shared(P, grid)
    conditional = np.asarray(P.fixed_reference_entry_conditional, dtype=float)
    if conditional.shape != (len(grid), len(P.z_grid)):
        raise ValueError("Fixed entrant support has unexpected shape")
    if not np.isfinite(conditional).all() or np.any(conditional < 0):
        raise ValueError("Invalid fixed entrant weights")
    # Occupied entrants begin childless. Read the actual shared arrays: in this
    # fixed eqscale reference childless hb is zero despite P.h_bar_0 being 4.
    childless_c = float(sd.cb_flat[0, 0])
    childless_h = float(sd.hb_flat[0, 0])
    if not np.isfinite(childless_c + childless_h):
        raise ValueError("Invalid shared childless cost")
    rows = []
    for bi, zi in zip(*np.nonzero(conditional > 0)):
        wealth, z = float(grid[bi]), float(P.z_grid[zi])
        income = income_at_state(P, 0, 0, z)
        for factor in price_factors:
            q = float(context["q_ref"]) * float(factor)
            rent = float(P.user_cost_rate) * q
            # The exhaustive exact branch uses positive intratemporal surplus;
            # c_min and 0.01 are output fallbacks, not feasibility thresholds.
            # The optimizer's upper saving candidate is Rv-cb-r*hb-1e-6.
            housing_cost = rent * childless_h
            required = childless_c + housing_cost - float(P.R_gross) * wealth - income
            rows.append(dict(b_index=int(bi), z_index=int(zi), wealth=wealth, z=z,
                             conditional_weight=float(conditional[bi, zi]), price_factor=float(factor),
                             income=float(income), childless_consumption_cost=childless_c,
                             childless_housing_cost=housing_cost,
                             required_d=max(0.0, float(required))))
    if not rows:
        raise ValueError("No occupied entrant cells")
    bound = max(row["required_d"] for row in rows)
    # One-cent operational mesh with the native 1e-6 upper-search buffer.
    mesh = 0.01
    selected = math.ceil((bound + 1e-6 - 1e-14) / mesh) * mesh
    return dict(unrounded_requirement=bound, selected_d_bar=round(selected, 12),
                operational_mesh=mesh, native_upper_search_buffer=1e-6,
                margin_above_requirement=selected-bound,
                robustness_price_factors=[0.85, 1.15],
                bound_price_factors=list(price_factors), childless_cb=childless_c,
                childless_hb=childless_h, occupied_cells=len(rows),
                limiting_cells=[r for r in rows if abs(r["required_d"]-bound)<1e-12],
                scope="necessary entrant budget bound only; c_min/output housing floors do not enter exact branch; native Bellman/KFE is decisive")


def run_phase_a(context, budget):
    from single_price import solve_fixed_price
    out = Path(context["out"])
    bound = entrant_budget_bound(context)
    context["write_json"](out / "entry_budget_bound.json", bound)
    budget.progress("phase_a_bound", bound=bound)
    selected = float(bound["selected_d_bar"])
    # No automatic search or cap widening: at most one selected-cap PE solve.
    stage = solve_fixed_price(context, selected, float(context["q_ref"]), budget,
                              "phase_a_selected_qref", out / "phase_a_selected_qref")
    return dict(status="phase_a_native_bellman_kfe_pass", selected_d_bar=selected,
                selected_stage=str(stage["stage_dir"]), summary=stage["summary"],
                remaining_lifecycle=budget.remaining_lifecycle,
                entry_budget_bound=bound, live_stage=stage, selected_live=stage)
