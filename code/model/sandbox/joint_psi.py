"""psi_mode: joint -- reproduce the retained corrected_initial calibration's
own (price, psi_child, property-tax rebate) simultaneous root, instead of the
nested price-GE / psi-bisection loop that sandbox/run_ss.py's "root" mode uses.

Traced call sequence (2026-09-17 re-derivation against
output/model/paper_baseline_sep14/initial_recipe/, the archived recipe that
reproduced the paper's initial steady state exactly on Torch job 17923835):
  run.sh -> run_e5f_joint_rebated_initial_scored.py --helper
    run_e5f_rebated_initial_overnight.py --joint run_e5f_joint_rebated_initial_probe.py
  -> candidate()/run_smoke() intercepts the pinned raw driver's subprocess call
     and reroutes it to run_e5f_joint_rebated_initial_probe.raw_mode, which
     monkeypatches driver.solve_balanced_initial_equilibrium ->
     joint_initial_solution (this module's import) and
     driver.calibration.solve_old_steady_state -> a single joint call.
  -> The raw driver (code/model/tools/run_e5f_initial_revision_probe.py,
     identical on main) builds its parameter object as:
       chain, model = run_e5f_open_population_transition.configure_sequential_model()
         (== audit_closed_reproductive_closure.load_chain(profile="e5f-floor"),
          the SAME chain/profile sandbox/run_ss.py already uses)
       base, _ = e5f_stationary_paygo.rebase_initial_supply(old.parameters,
           asset_prices=old.solution.p_eq, elasticity=0.63)
       base = e5f_parenthood_utility.initialize_parenthood_utility(base)
       candidate = e5f_parenthood_utility.validate_parenthood_candidate(
           structural_candidate, require_complete=True)
       base = e5f_parenthood_utility.bind_parenthood_utility(base, candidate)
       base.psi_child = initial_psi
     `old.parameters` is a historical checkpoint
     (normalized_old.pkl.gz, sha256 dfeb34f0..., under a Torch-scratch-only
     path) that is NOT available locally and cannot be fetched under this
     sandbox's no-cluster-jobs restriction; the already-solved verified
     replay checkpoints that ARE present locally
     (output/model/paper_baseline_sep14/replay_20260917/native_output/raw/
     repetition_0{1,2}/initial_state.pkl.gz) were pickled under Python
     >=3.13/numpy>=2 and cannot be unpickled by this sandbox's pinned
     venv (Python 3.10.10, numpy 1.24.3) either. So `old.parameters` itself
     cannot be reproduced bit-for-bit here.
  -> HOWEVER: a direct attribute-by-attribute diff (2026-09-17, no solve) of
     (a) sandbox/run_ss.py's build_overrides()-assembled parameter object
     against (b) a parameter object built by literally calling
     rebase_initial_supply + initialize_parenthood_utility +
     bind_parenthood_utility on a fresh e5f-floor+e5f_income_entry_overrides()
     base (the best available local proxy for `old.parameters`, chosen
     because both initialize_parenthood_utility's _validate_lifecycle check
     and PARENTHOOD_SEARCH_DOMAIN itself require exactly that base's
     lifecycle/income-grid shape) found ONLY THREE differing attributes out
     of 229, two of which are deterministic (independent of the unavailable
     `old.parameters` values) and were previously missing from
     sandbox/run_ss.py's build_overrides():
       - xi_supply: build_overrides() left the e5f-floor default [1.75];
         rebase_initial_supply ALWAYS resets it to [0.63] regardless of
         `old.parameters`, since it is passed in as a literal argument
         (elasticity=.63 at run_e5f_initial_revision_probe.py:84).
       - c_bar_n: build_overrides() left the e5f-floor default 0.48;
         initialize_parenthood_utility's _FIXED_UTILITY contract ALWAYS
         forces c_bar_n=0.0 (and c_bar_0=0.0, sigma=2.0,
         preference_spec="eqscale", eqscale_form="power",
         child_room_floor=True), regardless of `old.parameters`, because
         the parenthood-only utility is a strict fixed contract
         (e5f_parenthood_utility.py:16-22).
     The third (tenure_choice_kappa: sandbox already forces 0.005 vs a
     fresh base's 0.0 default) is NOT a bug: nothing in the recipe's own
     candidate-binding touches tenure_choice_kappa, so it is only ever set
     by `old.parameters`' own history, and sandbox's existing hard-coded
     0.005 override (run_ss.py, "retained: externally fixed") was already
     independently established from the retained calibration's own record.
  -> This module now calls the *actual* rebase_initial_supply /
     initialize_parenthood_utility / validate_parenthood_candidate /
     bind_parenthood_utility functions (imported, not reimplemented) to
     assemble the joint-root's base parameter object, which repairs both
     xi_supply and c_bar_n. This closes the identifiable part of the gap;
     it does NOT certify a bit-for-bit match, since `old.parameters`' own
     historical values for every OTHER attribute (income grid draw, entry
     wealth calibration, survival schedule, etc., to the extent any of
     those differ from a freshly-built e5f-floor+e5f_income_entry_overrides()
     base) remain unverifiable without the unavailable checkpoint.
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
  -> Price is not separately anchored: it is the first joint-root
     coordinate, solved simultaneously with psi_child and the transfer, from
     the same FALLBACK_PRICE_SEED documented below (the true `old.solution
     .p_eq` seed is unavailable for the same checkpoint reason as above).

Only the fetched copies under
output/model/e5f_final_night_20260913/corrected_initial_template_v6_fetched/
are imported here (never copied into code/); this module and run_ss.py are
the only sandbox files touched. code/model/tools/e5f_social_security.py does
not exist in the live tree (it is a dependency of the fetched
e5f_stationary_paygo.py, itself from the historical corrected_initial_source
commit 70abd4a8): the only local copy found is
tmp/e5f_matched_pf/code/model/tools/e5f_social_security.py, added to
sys.path read-only for import purposes only, exactly as the fetched files
require it. rebase_initial_supply is imported from that same fetched
e5f_stationary_paygo.py (byte-identical to code/model/tools/e5f_stationary_
paygo.py, diffed 2026-09-17); initialize_parenthood_utility /
validate_parenthood_candidate / bind_parenthood_utility are imported from
the live code/model/tools/e5f_parenthood_utility.py, since that file is not
part of the archived initial_recipe/ at all (it lives on main).
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


# Literal override-dict keys produced by run_ss.py's build_overrides() for the
# 9 free structural coordinates (via display_theta_to_overrides), plus the
# hbar_child_rooms=0.0 zero restriction it also sets. These must be stripped
# from the "base" overrides before assembly, since the recipe's own
# bind_parenthood_utility (not a raw override-dict merge) is what binds them.
_CANDIDATE_OVERRIDE_KEYS = frozenset({
    "H0", "chi", "first_birth_fixed_cost", "kappa_fert", "kappa_fert_continuation",
    "theta0", "theta1", "hbar_first_child_jump", "hbar_child_rooms", "beta",
})


def _load_fetched() -> tuple[Any, Any]:
    for path in (FETCHED_ROOT, LEGACY_DEPENDENCY_ROOT):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    import e5f_stationary_paygo as paygo  # noqa: E402
    import run_e5f_joint_rebated_initial_probe as joint_probe  # noqa: E402

    return paygo, joint_probe


def _assemble_unsolved_parameters(model: Any, overrides: dict[str, Any]) -> Any:
    """Build a fully-populated parameter object without a Bellman/GE solve.

    Mirrors intergen_eqscale_seq_optimized/solver.py:run_model_cp_dt's own
    pre-solve steps exactly (setup_parameters -> apply_overrides -> the same
    handful of derived-field finalizations), stopping before make_grid/solve.
    This is the "base" that the archived recipe's old.parameters would have
    supplied, standing in for that unavailable historical checkpoint (see
    this module's docstring)."""
    P = model.setup_parameters()
    P = model.apply_overrides(P, overrides)
    if not hasattr(P, "beta") or P.beta is None:
        P.beta = 1 / (1 + P.rho) if hasattr(P, "rho") else 0.96
    P.rho = 1 / P.beta - 1
    P.rho_hat = P.rho
    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.R_gross = 1 + P.q
    P.phi = np.asarray(P.phi, dtype=float).reshape(-1)
    if P.phi.size == 1:
        P.phi = P.phi.item() * np.ones(P.n_parity)
    if getattr(P, "entry_init_override", None) is not None:
        e0 = np.maximum(np.asarray(P.entry_init_override, dtype=float).reshape(-1), 0)
        P.entry_shares = e0 / e0.sum()
        P.entry_by_loc = (1 / P.J) * P.entry_shares
    if not hasattr(P, "use_stochastic_aging"):
        P.use_stochastic_aging = False
    P = model.finalize_location_choice_spec(P)
    return P


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
    candidate_theta: dict[str, float],
    trace_path: Path | None = None,
) -> tuple[Any, Any, Any, float, dict[str, Any]]:
    """Same return contract as run_e5f_transition_calibration.solve_old_steady_state,
    but the psi/price/rebate closure is the fetched joint three-residual root,
    seeded at the retained solution (psi, transfer) and this sandbox's own
    nested-loop converged price (see RETAINED_*_SEED / FALLBACK_PRICE_SEED).

    `base_overrides` is run_ss.py's fully-merged override dict (chain
    defaults + e5f_income_entry_overrides + any spec deltas + the flattened
    9-coordinate candidate); `candidate_theta` is the same 9 coordinates in
    their DISPLAY names (beta_annual, h_P, ...), as loaded from the retained
    candidate_result.json. The candidate keys are stripped back out of
    base_overrides and rebound canonically via bind_parenthood_utility, per
    this module's docstring."""
    if not normalize:
        raise NotImplementedError("psi_mode: joint only implements the normalized (root) case")
    if abs(float(completed_fertility_target) - 2.1) > 0:
        raise ValueError("Joint adapter is pinned to the 2.1 completed-fertility target")

    paygo, joint_probe = _load_fetched()
    import e5f_parenthood_utility as parenthood  # noqa: E402  (code/model/tools/, live on main)
    from intergen_eqscale_seq_optimized import solver as model  # noqa: E402

    started = time.perf_counter()
    base_overrides_stripped = {
        key: value for key, value in base_overrides.items() if key not in _CANDIDATE_OVERRIDE_KEYS
    }
    base_overrides_stripped.setdefault("property_tax_lump_sum_transfer", 0.0)
    # Guard: base_overrides is run_ss.py's fully-merged dict, which already
    # flattened candidate_theta into these literal keys via
    # display_theta_to_overrides. Since they are about to be dropped in
    # favor of binding candidate_theta canonically below, fail loudly if a
    # spec silently changed one of them instead of only the mechanism
    # switches specs are documented to touch -- a stripped-and-ignored spec
    # override on a structural coordinate would silently change the science.
    implied = dict(candidate_theta)
    if "beta_annual" in implied:
        implied["beta"] = implied.pop("beta_annual") ** 4.0
    if "h_P" in implied:
        implied["hbar_first_child_jump"] = implied.pop("h_P") - float(implied.get("hbar_child_rooms", 0.0))
    implied.setdefault("hbar_child_rooms", 0.0)
    for key in _CANDIDATE_OVERRIDE_KEYS & set(base_overrides):
        if key not in implied or not math.isclose(float(base_overrides[key]), float(implied[key]), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(
                f"psi_mode joint cannot honor a spec override on structural coordinate {key!r} "
                "(base_overrides diverges from candidate_theta); this coordinate is bound "
                "canonically via bind_parenthood_utility(candidate_theta), not from base_overrides."
            )
    base = _assemble_unsolved_parameters(model, base_overrides_stripped)
    # elasticity=0.63 is the recipe's own literal argument
    # (run_e5f_initial_revision_probe.py:84); the anchor price only feeds the
    # H0 recomputation that bind_parenthood_utility immediately overwrites
    # with the candidate's literal H0, so any positive price is equivalent.
    anchor_price = np.asarray(base.r_bar, dtype=float) / float(base.user_cost_rate)
    base, _rebase_info = paygo.rebase_initial_supply(base, asset_prices=anchor_price, elasticity=0.63)
    base = parenthood.initialize_parenthood_utility(base)
    checked_candidate = parenthood.validate_parenthood_candidate(candidate_theta, require_complete=True)
    P0 = parenthood.bind_parenthood_utility(base, checked_candidate)
    P0.psi_child = float(initial_psi)
    P0.property_tax_lump_sum_transfer = float(base_overrides_stripped["property_tax_lump_sum_transfer"])
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
