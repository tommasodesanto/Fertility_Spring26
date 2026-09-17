"""psi_mode: joint -- reproduce the retained corrected_initial calibration's
own (price, psi_child, property-tax rebate) simultaneous root, instead of the
nested price-GE / psi-bisection loop that sandbox/run_ss.py's "root" mode uses.

Traced call sequence (see the assistant's report for the full derivation):
  run.sh -> run_e5f_joint_rebated_initial_scored.py --helper
    run_e5f_rebated_initial_overnight.py --joint run_e5f_joint_rebated_initial_probe.py
  -> candidate()/run_smoke() intercepts the pinned raw driver's subprocess call
     and reroutes it to run_e5f_joint_rebated_initial_probe.raw_mode, which
     monkeypatches driver.solve_balanced_initial_equilibrium ->
     joint_initial_solution (this module's import) and
     driver.calibration.solve_old_steady_state -> a single joint call.
  -> joint_initial_solution (run_e5f_joint_rebated_initial_probe.py:106) roots
     three residuals [housing excess demand / supply, tfr - 2.1,
     rebate residual / max(|revenue|,|outlays|)] over transformed
     [log(price), psi_child, log(transfer)] with a damped Gauss-Newton/line-
     search iteration (solve_three_residual_root, same file, lines 52-104):
     finite-difference steps [0.02, 0.02, 0.05], Levenberg ridge 1e-8, step
     caps [0.30, 0.25, 0.70], backtracking line search (1, 0.5, 0.25),
     bounded evaluations in [4, 20], gates
     housing<=2.5e-5, |tfr-2.1|<=5e-4, rebate<=1e-6.
  -> bind_initial_balanced_pension / certify_initial_pension
     (e5f_stationary_paygo.py, fetched) rebuild the analytic one-market PAYGO
     pension at payroll_tax=0.179 on every fixed-price/psi/transfer evaluation
     and re-certify it on the accepted root.
  -> H0 is NOT re-anchored inside this joint call: the raw driver rebases
     housing supply (xi_supply=0.63) against a *different, unavailable*
     pre-parenthood-utility checkpoint before binding the 9 free structural
     coordinates (H0 among them) from proposal.json -- and that literal
     candidate H0 overwrites whatever the rebase computed. This sandbox
     therefore uses the SAME H0=8.11210048056786 (and the same xi_supply the
     "e5f-floor" profile already sets) that run_ss.py's existing "root" mode
     uses, since that rebase-then-overwrite step is a no-op on H0 itself and
     its unavailable inputs (an "old" pre-utility checkpoint) cannot be
     reproduced locally.
  -> Price is NOT separately anchored either: it is the first joint-root
     coordinate, solved simultaneously with psi_child and the transfer.

Only the fetched copies under
output/model/e5f_final_night_20260913/corrected_initial_template_v6_fetched/
are imported here (never copied into code/); this module and run_ss.py are
the only sandbox files touched. code/model/tools/e5f_social_security.py does
not exist in the live tree (it is a dependency of the fetched
e5f_stationary_paygo.py, itself from the historical corrected_initial_source
commit 70abd4a8): the only local copy found is
tmp/e5f_matched_pf/code/model/tools/e5f_social_security.py, added to
sys.path read-only for import purposes only, exactly as the fetched files
require it.
"""
from __future__ import annotations

import copy
import csv
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

SANDBOX_ROOT = Path(__file__).resolve().parent
MODEL_ROOT = SANDBOX_ROOT.parent
REPO_ROOT = MODEL_ROOT.parents[1]
FETCHED_ROOT = REPO_ROOT / "output/model/e5f_final_night_20260913/corrected_initial_template_v6_fetched"
LEGACY_DEPENDENCY_ROOT = REPO_ROOT / "tmp/e5f_matched_pf/code/model/tools"

PAYROLL_TAX = 0.179  # validate_scientific_contract's pinned PAYGO tax rate
MARGINAL_TOLERANCE = 1e-9
FISCAL_TOLERANCE = 1e-6

# Retained corrected_initial solution (candidate_result.json / target_fit.csv):
# psi and transfer are recorded exactly; the retained equilibrium price is not
# recorded anywhere in output/model/e5f_final_night_20260913/corrected_initial/
# (checked candidate_result.json, summary.json, evaluation_summary.json,
# parameters.csv/parameters_raw.csv -- no price field). Falling back to this
# sandbox's own nested-loop ("root" mode) converged price for the same theta.
RETAINED_PSI_SEED = 0.1489153145785918
RETAINED_TRANSFER_SEED = 0.18660898018995767
FALLBACK_PRICE_SEED = 0.7931310188535463  # output/model/sandbox/baseline/parameters.csv:_solved_price


def _load_fetched() -> tuple[Any, Any]:
    for path in (FETCHED_ROOT, LEGACY_DEPENDENCY_ROOT):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    import e5f_stationary_paygo as paygo  # noqa: E402
    import run_e5f_joint_rebated_initial_probe as joint_probe  # noqa: E402

    return paygo, joint_probe


def _fixed_evaluation(model, bind_pension, price, psi, transfer, *, parameters, b_grid, payroll_tax):
    """Mirror run_e5f_joint_rebated_initial_probe.joint_initial_solution's inner
    `fixed()` closure exactly (same file, lines 116-135), without editing that
    fetched file -- needed here only so an explicit seed can replace its
    internally re-derived transfer_start (see solve_old_steady_state_joint)."""
    P = copy.deepcopy(parameters)
    P.psi_child = float(psi)
    P.property_tax_lump_sum_transfer = float(transfer)
    P, predicted = bind_pension(P, payroll_tax=payroll_tax)
    solution = model.solve_markov_income_at_prices(np.array([price]), P, b_grid, verbose=False, fast_stats=False)
    solution = model.attach_markov_market_accounting(solution, P, b_grid)
    demand, _ = model.markov_market_housing_demand(solution, P, b_grid)
    supply = float(np.asarray(solution.housing_supply).reshape(-1)[0])
    housing = (float(demand[0]) - supply) / max(abs(supply), 1e-12)
    calibration = __import__("importlib").import_module(model.__package__ + ".calibration")
    fertility = float(calibration.extract_moments(solution, P)["tfr"])
    revenue = float(solution.property_tax_revenue)
    outlays = float(solution.property_tax_transfer_outlays)
    rebate = float(solution.property_tax_budget_residual)
    return dict(solution=solution, parameters=P, price=np.array([price]), predicted_pension=predicted,
                completed_fertility=fertility, property_tax_revenue=revenue, rebate_outlays=outlays,
                residual=np.array([housing, fertility - 2.1, rebate / max(abs(revenue), abs(outlays), 1e-12)]))


def solve_old_steady_state_joint(
    chain: Any,
    base_overrides: dict[str, Any],
    *,
    initial_psi: float,
    completed_fertility_target: float,
    completed_fertility_tolerance: float,
    normalize: bool,
    trace_path: Path | None = None,
) -> tuple[Any, Any, Any, float, dict[str, Any]]:
    """Same return contract as run_e5f_transition_calibration.solve_old_steady_state,
    but the psi/price/rebate closure is the fetched joint three-residual root,
    seeded at the retained solution (psi, transfer) and this sandbox's own
    nested-loop converged price (see RETAINED_*_SEED / FALLBACK_PRICE_SEED)."""
    if not normalize:
        raise NotImplementedError("psi_mode: joint only implements the normalized (root) case")
    if abs(float(completed_fertility_target) - 2.1) > 0:
        raise ValueError("Joint adapter is pinned to the 2.1 completed-fertility target")

    paygo, joint_probe = _load_fetched()
    from intergen_eqscale_seq_optimized import solver as model  # noqa: E402

    started = time.perf_counter()
    seed_overrides = dict(base_overrides)
    seed_overrides["psi_child"] = float(initial_psi)
    seed_overrides.setdefault("property_tax_lump_sum_transfer", 0.0)
    _, P0, _ = chain.run_model_cp_dt(dict(seed_overrides), verbose=False)
    if int(P0.I) != 1 or float(P0.tau_H) != 0.04:
        raise ValueError("Joint initial probe requires the saved one-market 1% annual tax economy")
    b_grid = model.make_grid(P0)

    trace_rows: list[list[float]] = []

    def evaluate(transformed):
        price = math.exp(float(transformed[0]))
        psi = float(transformed[1])
        transfer = math.exp(float(transformed[2]))
        result = _fixed_evaluation(model, paygo.bind_initial_balanced_pension, price, psi, transfer,
                                    parameters=P0, b_grid=b_grid, payroll_tax=PAYROLL_TAX)
        residual = result["residual"]
        trace_rows.append([len(trace_rows) + 1, price, psi, transfer,
                            float(residual[0]), float(residual[1]), float(residual[2])])
        if trace_path is not None:
            write_header = not trace_path.exists()
            trace_path.parent.mkdir(parents=True, exist_ok=True)
            with trace_path.open("a", newline="") as handle:
                writer = csv.writer(handle)
                if write_header:
                    writer.writerow(["evaluation", "price", "psi_child", "transfer",
                                      "residual_housing", "residual_fertility_gap", "residual_rebate"])
                writer.writerow(trace_rows[-1])
        return result

    start = np.array([
        math.log(FALLBACK_PRICE_SEED),
        RETAINED_PSI_SEED,
        math.log(RETAINED_TRANSFER_SEED),
    ])
    housing_tolerance = min(float(P0.tol_eq), joint_probe.HOUSING_TOLERANCE)

    try:
        root = joint_probe.solve_three_residual_root(evaluate, start, housing_tolerance=housing_tolerance)
    except TimeoutError:
        seconds = time.perf_counter() - started
        raise TimeoutError(
            f"Joint initial root exhausted its fixed-price evaluation budget after "
            f"{len(trace_rows)} evaluations ({seconds:.1f}s); trace written to {trace_path}"
        ) from None

    selected = root["result"]
    solution, P, price = selected["solution"], selected["parameters"], selected["price"]
    pension = paygo.certify_initial_pension(solution.g, P, marginal_tolerance=MARGINAL_TOLERANCE,
                                            fiscal_tolerance=FISCAL_TOLERANCE)
    seconds = time.perf_counter() - started
    diagnostics = {
        "status": "derived_intercept_joint",
        "psi_child": float(P.psi_child),
        "completed_fertility": float(selected["completed_fertility"]),
        "target": float(completed_fertility_target),
        "absolute_gap": abs(float(selected["completed_fertility"]) - float(completed_fertility_target)),
        "stationary_solves": int(root["evaluations"]) + 1,
        "stationary_solve_seconds": seconds,
        "pension": pension,
        "property_tax_revenue": float(selected["property_tax_revenue"]),
        "rebate_outlays": float(selected["rebate_outlays"]),
        "transfer_period_units": float(P.property_tax_lump_sum_transfer),
        "seed": dict(price=FALLBACK_PRICE_SEED, psi=RETAINED_PSI_SEED, transfer=RETAINED_TRANSFER_SEED),
    }
    return solution, P, price, seconds, diagnostics
