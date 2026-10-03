"""Play-button driver for a reproducible fixed-price household-model run.

Edit the dictionaries and FIXED_PRICE below, then run this file from an IDE.
Importing it is inert. It uses the authenticated soft-timing ModelPlayground,
solves once at the stated price, and never searches for a market-clearing price.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import hashlib
import os
from pathlib import Path
import signal
import sys
import time
import traceback

np = None

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "code/model/tools"

# EDITABLE selected internal coordinates: these are bound by ModelPlayground
# after the direct native P inputs below, so they take precedence on overlap.
INTERNAL_PARAMETERS = {
    "beta_annual": 0.9671931106058198,  # Annual discount factor (four-year beta is compounded).
    "chi": 1.0972389108161984,  # Owner housing-service premium.
    "first_birth_fixed_cost": 0.35161589615733957,  # Fixed utility cost of the first birth.
    "kappa_fert": 0.12428888389792507,  # First-birth choice shock/logit scale.
    "kappa_fert_continuation": 0.36310164268347417,  # Later-birth attempt choice shock/logit scale.
    "theta0": 0.10630491855239614,  # Bequest utility scale.
    "h_P": 2.504044687713671,  # Physical room floor added at the first child.
    "child_benefit_curvature": 0.10071173914767594,  # Curvature of child benefit by children at home.
    "tenure_choice_kappa": 0.014279685267457497,  # Tenure-choice logit scale.
    "psi_child": 0.17840979194160872,  # Child benefit scale.
}

# EDITABLE external inputs, assigned directly to same-named native P fields.
# Income already includes the pension profile; the script checks that its
# retirement segment is constant before synchronizing pension fields.
EXTERNAL_INPUTS = {
    "sigma": 2.0,  # Consumption curvature.
    "alpha_cons": 0.733,  # Consumption share in the within-period utility aggregator.
    "theta1": 0.008193084126995582,  # Linear fertility utility term.
    "theta_n": 0.0,  # Number-of-children utility curvature parameter.
    "R_gross": 1.08243216,  # Gross annual asset return.
    "delta": 0.05545379079326218,  # Annual housing depreciation rate.
    "tau_H": 0.042393443095490375,  # Annual property-tax rate.
    "psi": 0.06,  # Selling-cost fraction.
    "phi": [0.8, 0.8, 0.8, 0.8],  # Uniform financed share by owner housing rung.
    "unsecured_credit_limit": 0.0,  # Maximum unsecured borrowing (zero here).
    "c_min": 0.04,  # Consumption floor.
    "owner_size_cost": 0.0,  # Owner housing-size adjustment cost.
    "owner_size_cost_power": 2.0,  # Power in the owner size-cost function.
    "owner_size_cost_ref": 6.0,  # Reference owner housing size.
    "retirement_income_z_scale": 0.0,  # Income-type scaling in retirement.
    "fecundity_omega1": 0.02,  # First fecundity profile coefficient.
    "fecundity_omega2": 0.134,  # Second fecundity profile coefficient.
    "property_tax_lump_sum_transfer": 0.0,  # Lump-sum property-tax transfer.
    "H0": [6.757074077757929],  # Fixed housing supply by location.
    "eta_supply": [1.75],  # Housing supply elasticity by location.
    "xi_supply": [0.63],  # Housing supply scale by location.
    "r_bar": [0.16],  # Reference unit rent by location.
    "income": [[  # Disposable income by location and lifecycle period.
        2.650830656801071, 2.650830656801071, 3.4664708588937074,
        3.4664708588937074, 4.078201010463186, 4.078201010463186,
        4.078201010463186, 4.017027995306238, 4.017027995306238,
        4.017027995306238, 3.8131179447830785, 3.8131179447830785,
        0.917784047463731, 0.917784047463731, 0.917784047463731,
        0.917784047463731, 0.917784047463731,
    ]],
    "survival_probs": [1.0] * 12 + [  # One-period survival probabilities, length J - 1.
        0.9391263063710125, 0.9184976343249724,
        0.8849521927812863, 0.8300468061015381,
    ],
}

# Advanced direct assignments to native P fields. Empty by default; there is
# deliberately no hidden name mapping. Unknown fields and shape changes fail.
NATIVE_OVERRIDES = {}

# A single stationary household solve at this fixed price; no price root or fit.
FIXED_PRICE = 0.7266387868818555
RUN_BUDGET_SECONDS = 600
RUNS_ROOT = ROOT / "tmp/model_runs"
DERIVED_NATIVE_FIELDS = {"q", "user_cost_rate", "pension", "pension_by_loc"}
PARAMETER_BOUND_NATIVE_FIELDS = set(INTERNAL_PARAMETERS) | {
    "beta", "rho", "rho_hat", "eps_fert", "hbar_first_child_jump",
    "hbar_child_rooms", "child_room_floor",
}
CACHED_GRID_PROCESS_FIELDS = {
    "earnings_transaction_grid", "b_min", "b_max", "b_grid_power",
    "income_shock_persistence",
}


def _as_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _as_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_as_jsonable(v) for v in value]
    if isinstance(value, (str, bool, int, float)) or value is None:
        return value
    return str(value)


def _assign_native_inputs(P, values, label):
    for name, value in values.items():
        if not hasattr(P, name):
            raise KeyError(f"{label} names unknown native P field {name!r}")
        current = getattr(P, name)
        if isinstance(current, np.ndarray):
            candidate = np.asarray(value, dtype=current.dtype)
            if candidate.shape != current.shape:
                raise ValueError(
                    f"{label}.{name} has shape {candidate.shape}; expected {current.shape}"
                )
            setattr(P, name, candidate.copy())
        elif isinstance(current, (int, float, np.generic)) and not isinstance(current, bool):
            candidate = float(value)
            if not np.isfinite(candidate):
                raise ValueError(f"{label}.{name} must be finite")
            setattr(P, name, candidate)
        else:
            raise TypeError(f"{label}.{name} is not a supported scalar or ndarray field")


def _validate_native_overrides(overrides=None):
    """Reject direct edits that cannot update the authenticated cached grids/process."""
    overrides = NATIVE_OVERRIDES if overrides is None else overrides
    blocked = {
        name for name in overrides
        if name in CACHED_GRID_PROCESS_FIELDS or name.startswith(("b_core_", "b_frac_"))
    }
    if blocked:
        raise ValueError(
            "These native overrides require rebuilding grids/process consistently; "
            "this runner uses saved grids: " + ", ".join(sorted(blocked))
        )


def _validate_phi(phi_values):
    phi = np.asarray(phi_values, dtype=float)
    if phi.shape != (4,) or not np.isfinite(phi).all() or np.any(phi < 0.0) or np.any(phi > 1.0):
        raise ValueError("phi must contain four finite financed shares in [0, 1]")
    if not np.allclose(phi, phi[0], rtol=0.0, atol=1e-12):
        raise ValueError(
            "This runner requires one uniform financed share across all four owner rungs; "
            "nonuniform phi is unsupported because native credit-floor and household "
            "down-payment rules use different phi representations."
        )
    return phi


def _prepare_inputs(model):
    P = model.P
    _validate_native_overrides()
    _assign_native_inputs(P, EXTERNAL_INPUTS, "EXTERNAL_INPUTS")
    _assign_native_inputs(P, NATIVE_OVERRIDES, "NATIVE_OVERRIDES")

    _validate_phi(P.phi)

    income = np.asarray(P.income, dtype=float)
    if income.ndim != 2 or income.shape[1] < 1:
        raise ValueError("income must be a location-by-age array")
    pension_profile = income[:, int(P.J_R):int(P.J)]
    if pension_profile.shape[1] < 1 or income.shape[1] != int(P.J):
        raise ValueError("income age dimension must equal J and include the retirement segment J_R:J")
    if not np.allclose(pension_profile, pension_profile[:, :1], rtol=0.0, atol=1e-12):
        raise ValueError(
            "The supplied retirement income varies by age; this workflow only supports "
            "a constant retirement segment for pension synchronization."
        )
    pension_by_loc = pension_profile[:, 0].copy()
    if np.asarray(P.pension_by_loc).shape != pension_by_loc.shape:
        raise ValueError("pension_by_loc shape does not match income locations")
    P.pension_by_loc = pension_by_loc
    P.pension = float(pension_by_loc[0])

    # Derived native primitives: gross return determines q, and user cost adds
    # depreciation and property tax exactly as in the authenticated baseline.
    P.q = float(P.R_gross) - 1.0
    P.user_cost_rate = P.q + float(P.delta) + float(P.tau_H)
    model.params.clear()
    model.params.update({key: float(value) for key, value in INTERNAL_PARAMETERS.items()})

    baseline = model._authenticated_base_P
    changed = []
    for name in list(EXTERNAL_INPUTS) + list(NATIVE_OVERRIDES) + [
        "pension", "pension_by_loc", "q", "user_cost_rate",
    ]:
        before, after = getattr(baseline, name), getattr(P, name)
        equal = (np.array_equal(np.asarray(before), np.asarray(after), equal_nan=True)
                 if isinstance(before, np.ndarray) or isinstance(after, np.ndarray)
                 else before == after)
        if not equal:
            changed.append(name)
    return changed


def _validate_structure(P, b_grid):
    """Reject grid or core-state dimensions inconsistent with native P."""
    if int(P.Nb) != len(b_grid):
        raise ValueError(f"Nb={P.Nb} but b_grid has {len(b_grid)} nodes")
    if int(P.J) != np.asarray(P.income).shape[1] or int(P.J_R) >= int(P.J):
        raise ValueError("J/J_R do not match the income age dimension")
    if np.asarray(P.survival_probs).shape != (int(P.J) - 1,):
        raise ValueError("survival_probs must have J - 1 entries")
    if len(P.H_own) != int(P.n_house):
        raise ValueError("H_own length does not match n_house")
    if np.asarray(P.Pi_z).shape != (int(P.Nz), int(P.Nz)):
        raise ValueError("Pi_z shape does not match Nz")
    if np.asarray(P.Pi_child).shape != (
        int(P.n_parity), int(P.n_child_states), int(P.n_child_states)
    ):
        raise ValueError("Pi_child shape does not match child-state dimensions")


def _validate_solution_structure(solution, P, b_grid):
    """Check core native policy/distribution tensors against the stated grids."""
    state_shape = (
        int(P.Nb), int(P.n_house) + 1, int(P.n_sub), int(P.J), int(P.Nz),
        int(P.n_parity), int(P.n_child_states),
    )
    for name in ("V", "c_pol", "hR_pol", "bp_pol", "g", "g_stay_distribution",
                 "g_beginning_distribution"):
        value = np.asarray(getattr(solution, name))
        if value.shape != state_shape:
            raise ValueError(f"solution.{name} has shape {value.shape}; expected {state_shape}")
    if np.asarray(solution.b_grid).shape != np.asarray(b_grid).shape or not np.array_equal(
        solution.b_grid, b_grid
    ):
        raise ValueError("solution wealth grid does not match the requested native grid")


def _solver_diagnostics(solution):
    diagnostics = {"meaning": "Finite/value/probability/distribution checks only; no convergence certificate."}
    finite_fields = {}
    probability_fields = {}
    distribution_fields = {}
    for name, value in vars(solution).items():
        if not isinstance(value, np.ndarray) or value.dtype.kind not in "fiu":
            continue
        array = np.asarray(value)
        if array.dtype.kind in "fc":
            finite_fields[name] = {
                "finite_count": int(np.isfinite(array).sum()),
                "element_count": int(array.size),
                "finite_min": float(np.nanmin(array)) if np.isfinite(array).any() else None,
                "finite_max": float(np.nanmax(array)) if np.isfinite(array).any() else None,
            }
        lower_name = name.lower()
        if ("prob" in lower_name or "hazard" in lower_name
                or lower_name.startswith("frac_") or "_rate_by_" in lower_name
                or lower_name == "pop_share"):
            probability_fields[name] = {
                "min": float(np.nanmin(array)), "max": float(np.nanmax(array)),
                "within_unit_interval": bool(
                    np.isfinite(array).all() and np.all(array >= -1e-10)
                    and np.all(array <= 1.0 + 1e-10)
                ),
            }
        if name == "g" or (name.startswith("g") and "distribution" in lower_name):
            distribution_fields[name] = {
                "sum": float(np.nansum(array)),
                "minimum": float(np.nanmin(array)),
                "nonnegative": bool(np.isfinite(array).all() and np.all(array >= -1e-12)),
            }
    diagnostics["finite_arrays"] = finite_fields
    diagnostics["probability_like_arrays"] = probability_fields
    diagnostics["distribution_arrays"] = distribution_fields
    diagnostics["all_reported_finite_arrays_finite"] = all(
        item["finite_count"] == item["element_count"] for item in finite_fields.values()
    )
    diagnostics["all_probability_like_arrays_in_unit_interval"] = all(
        item["within_unit_interval"] for item in probability_fields.values()
    )
    for name in ("adult_entry_stationary_relative_gap",):
        value = getattr(solution, name, None)
        if value is not None:
            diagnostics[name] = _as_jsonable(value)
            diagnostics[name + "_meaning"] = "Renewal/stationary-entry diagnostic; not a convergence certificate."
    return diagnostics


def _print_solver_status(solution, diagnostics, solve_seconds):
    print(f"Household solve returned in {solve_seconds:.1f}s.")
    print(
        "Array checks: finite="
        f"{diagnostics['all_reported_finite_arrays_finite']}; probability/hazard bounds="
        f"{diagnostics['all_probability_like_arrays_in_unit_interval']}"
    )
    for name, values in diagnostics["distribution_arrays"].items():
        print(f"Distribution {name}: sum={values['sum']:.10g}, nonnegative={values['nonnegative']}")
    gap = getattr(solution, "adult_entry_stationary_relative_gap", None)
    if gap is not None:
        print(f"Renewal diagnostic adult_entry_stationary_relative_gap={float(np.asarray(gap).squeeze()):.8g}")
    demand, supply = getattr(solution, "housing_demand", None), getattr(solution, "housing_supply", None)
    if demand is not None and supply is not None:
        residual = np.asarray(demand, dtype=float) - np.asarray(supply, dtype=float)
        print(f"Housing demand minus supply at fixed price={np.array2string(residual, precision=8)}")
    print("These are diagnostics at fixed prices; no market root or convergence certificate is reported.")


class _BudgetExpired(TimeoutError):
    pass


def _on_alarm(_signum, _frame):
    raise _BudgetExpired(f"Run exceeded its {RUN_BUDGET_SECONDS}s budget")


def main():
    # Set thread limits before importing the native runtime.
    global np
    for key in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import numpy as _numpy
    np = _numpy
    import model_playground as playground_module
    from model_playground import ModelPlayground
    from model_run_io import reserve_run_directory, save_run

    if set(INTERNAL_PARAMETERS) != set(playground_module.PARAMETER_ORDER):
        raise ValueError(
            "INTERNAL_PARAMETERS must contain exactly the current PARAMETER_ORDER; "
            f"expected {playground_module.PARAMETER_ORDER}"
        )
    conflicts = (set(NATIVE_OVERRIDES) & DERIVED_NATIVE_FIELDS)
    if conflicts:
        raise ValueError(
            "NATIVE_OVERRIDES cannot set derived fields "
            f"{sorted(conflicts)}; edit R_gross, delta, tau_H, or income instead."
        )
    conflicts = set(NATIVE_OVERRIDES) & PARAMETER_BOUND_NATIVE_FIELDS
    if conflicts:
        raise ValueError(
            "NATIVE_OVERRIDES cannot shadow selected internal-parameter bindings "
            f"{sorted(conflicts)}; edit INTERNAL_PARAMETERS instead."
        )
    _validate_native_overrides()

    run_directory = reserve_run_directory(RUNS_ROOT)
    started_utc = datetime.now(timezone.utc).isoformat()
    start = time.monotonic()
    old_handler = signal.signal(signal.SIGALRM, _on_alarm) if hasattr(signal, "SIGALRM") else None
    if hasattr(signal, "alarm"):
        signal.alarm(RUN_BUDGET_SECONDS)
    try:
        model = ModelPlayground()
        _validate_structure(model.P, model.b_grid)
        changed_fields = _prepare_inputs(model)
        print("Fixed-price household run (not a market equilibrium or calibration).")
        print(f"Fixed price: {FIXED_PRICE:.12g}; threads: 1; total budget: {RUN_BUDGET_SECONDS}s.")
        print("Selected internal parameters:")
        for name in playground_module.PARAMETER_ORDER:
            print(f"  {name} = {INTERNAL_PARAMETERS[name]:.12g}")
        print("Changed native fields vs authenticated P baseline: " +
              (", ".join(changed_fields) if changed_fields else "none"))
        _validate_structure(model.P, model.b_grid)
        print("Solving backward household policies and forward distribution...")
        solve_started = time.monotonic()
        # The literal dictionary remains authoritative; no calibration/root loop.
        result = model.solve(price=FIXED_PRICE)
        solve_seconds = time.monotonic() - solve_started
        _validate_solution_structure(result.solution, result.P, model.b_grid)
        if not hasattr(result.solution, "timing"):
            result.solution.timing = "transaction_inside"
        result.label = "fixed-price stationary household solve; no market-price root/calibration"
        solver_diagnostics = _solver_diagnostics(result.solution)
        _print_solver_status(result.solution, solver_diagnostics, solve_seconds)

        from small_credit_lab.engine import diagnostics as native_diagnostics
        plots_directory = run_directory / "native_diagnostics"
        native_diagnostics.write_diagnostics(result.solution, result.P, plots_directory)
        plot_count = len(list(plots_directory.glob("*.png")))
        if plot_count != 17:
            raise RuntimeError(f"Expected 17 native diagnostic PNGs; found {plot_count}")

        elapsed_before_save = time.monotonic() - start
        reference_record = json.loads(playground_module.SELECTION.read_text(encoding="utf-8"))
        selected_source = ROOT / reference_record["source"]
        internal_changes = {
            key: {"authenticated": float(model._reference_params[key]),
                  "current": float(INTERNAL_PARAMETERS[key])}
            for key in playground_module.PARAMETER_ORDER
            if float(model._reference_params[key]) != float(INTERNAL_PARAMETERS[key])
        }
        metadata = {
            "workflow": "fixed-price household solve",
            "scope": "Authenticated original-timing soft reference; one stationary household solve at a fixed price.",
            "market_price_root_solved": False,
            "calibration_run": False,
            "convergence_certificate": None,
            "solver_diagnostics": solver_diagnostics,
            "started_utc": started_utc,
            "solve_seconds": solve_seconds,
            "elapsed_seconds_before_serialization": elapsed_before_save,
            "run_budget_seconds": RUN_BUDGET_SECONDS,
            "thread_limit": 1,
            "fixed_price": FIXED_PRICE,
            "internal_parameters": INTERNAL_PARAMETERS,
            "external_inputs": EXTERNAL_INPUTS,
            "native_overrides": NATIVE_OVERRIDES,
            "changed_native_fields_vs_authenticated_base": changed_fields,
            "changed_fields_status": "Experimental relative to the authenticated native P baseline.",
            "authenticated_internal_reference": {
                "selection_file": str(playground_module.SELECTION.relative_to(ROOT)),
                "source_file": reference_record["source"],
                "source_sha256": hashlib.sha256(selected_source.read_bytes()).hexdigest(),
                "source_record_sha256": reference_record["source_sha256"],
                "parameters": model._reference_params,
            },
            "internal_parameter_changes_vs_authenticated_reference": internal_changes,
            "runner_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "derived_fields": {
                "q": "R_gross - 1",
                "user_cost_rate": "q + delta + tau_H",
                "pension": "constant income profile over J_R:J; synchronized after validation",
            },
            "native_diagnostic_png_count": plot_count,
            "native_diagnostics_directory": "native_diagnostics",
        }
        def final_metadata():
            return {"elapsed_seconds_through_roundtrip_validation": time.monotonic() - start}

        save_run(result, run_metadata=metadata, runs_root=RUNS_ROOT,
                 run_directory=run_directory, update_latest=True,
                 metadata_finalizer=final_metadata)
        total_elapsed = time.monotonic() - start
        print(f"Saved validated fixed-price run: {run_directory}")
        print(f"Native diagnostics: {plots_directory} ({plot_count} PNGs)")
        print(f"Total elapsed including save and round-trip validation: {total_elapsed:.1f}s.")
        print("No market-price root or calibration was run; fixed-price residuals are diagnostics only.")
        return run_directory
    except Exception:
        failure = {
            "status": "failed",
            "started_utc": started_utc,
            "failed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - start,
            "run_budget_seconds": RUN_BUDGET_SECONDS,
            "error": traceback.format_exc(),
            "latest_pointer_updated": False,
        }
        (run_directory / "failure.json").write_text(
            json.dumps(failure, indent=2) + "\n", encoding="utf-8"
        )
        raise
    finally:
        if hasattr(signal, "alarm"):
            signal.alarm(0)
        if old_handler is not None:
            signal.signal(signal.SIGALRM, old_handler)


if __name__ == "__main__":
    main()
