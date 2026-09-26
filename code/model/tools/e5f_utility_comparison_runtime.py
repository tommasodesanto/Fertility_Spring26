"""Opt-in four-arm runtime on the immutable September 25 nightpair source.

No source default or historical share utility is changed. The share composite
normalization is an experimental economic restriction requiring author review.
Imports/validation, including precompute_shared, must run on Torch.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np

ARMS = {
    "floor_linear": ("floor", 1.0),
    "floor_concave": ("floor", 0.86),
    "shares_linear": ("shares", 1.0),
    "shares_concave": ("shares", 0.86),
}
PENSION_RATIO = 0.2294460118659327
ALPHA0 = 0.733
COMMON = ("H0", "beta_annual", "chi", "first_birth_fixed_cost",
          "kappa_fert", "kappa_fert_continuation", "theta0")
SHARE_COORDINATES = ("delta_alpha_jump", "delta_alpha")


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def import_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def reference_composite_factor(alpha, reference_rent):
    """A(alpha)=K(alpha0,r*)/K(alpha,r*); K=a^a((1-a)/r*)^(1-a).

    At common prices (1,r*), max_{c+r*s=E} A*c^a*s^(1-a)=K(alpha0,r*)*E.
    Hence the scale alone determines compensated material expenditure at r*.
    This is substantive across child states; conditional demands alone do not
    identify it. It leaves childless and constant-alpha utility exactly intact.
    """
    alpha = np.asarray(alpha, dtype=float)
    rent = float(reference_rent)
    if not math.isfinite(rent) or rent <= 0 or not np.isfinite(alpha).all():
        raise ValueError("finite positive reference rent and finite shares required")
    if np.any(alpha <= 0) or np.any(alpha >= 1):
        raise ValueError("expenditure shares must be interior")
    log_k = alpha * np.log(alpha) + (1 - alpha) * np.log((1 - alpha) / rent)
    log_k0 = ALPHA0 * math.log(ALPHA0) + (1 - ALPHA0) * math.log((1 - ALPHA0) / rent)
    # Exact 1.0 avoids a rounding change to the inherited childless utility.
    return np.where(alpha == ALPHA0, 1.0, np.exp(log_k0 - log_k))


def comparison_shared(original_shared, P, grid):
    """Build native arrays then change only the declared preference objects."""
    shared = original_shared(P, grid)
    arm = getattr(P, "utility_comparison_arm", None)
    if arm not in ARMS:
        raise ValueError("explicit four-arm utility specification required")
    housing, exponent = ARMS[arm]
    if (P.child_state_mode != "independent_count" or int(P.n_parity) != 4
            or int(P.n_child_states) != 4 or P.alpha_cons != ALPHA0
            or P.sigma != 2.0 or P.eqscale_form != "power"):
        raise ValueError("four-arm runtime requires the unchanged native lifecycle and material curvature")
    if housing == "shares":
        if P.child_room_floor or P.hbar_first_child_jump != 0 or P.hbar_child_rooms != 0:
            raise ValueError("share arm cannot retain a housing floor")
        factors = reference_composite_factor(shared.alpha_flat, P.utility_reference_rent)
        # Native CRRA is escale * Q^(1-sigma)/(1-sigma).
        shared.escale_flat = shared.escale_flat * factors ** (1.0 - P.sigma)
    if exponent != 1.0:
        benefit = np.zeros_like(shared.psi_v)
        for n in range(int(P.n_parity)):
            for m in range(1, min(n + 1, int(P.n_child_states))):
                benefit[n, m] = P.psi_child * float(m) ** exponent
        shared.psi_v = benefit
        shared.psi_flat = benefit.reshape(1, shared.nc, order="F")
        # The native compression also stores the benefit in type_psi.
        triples = np.column_stack([shared.c_bar.reshape(-1, order="F"),
                                   shared.h_bar.reshape(-1, order="F"),
                                   benefit.reshape(-1, order="F")])
        unique, mapping = np.unique(triples, axis=0, return_inverse=True)
        shared.n_types, shared.type_map = len(unique), mapping
        shared.type_cb, shared.type_hb, shared.type_psi = unique.T
    return shared


def pension_tax_from_demographics(P):
    """Derive tau=rho*N_retired/N_working for equal per-retiree pensions."""
    from e5f_stationary_paygo import stationary_age_income_mass
    if float(getattr(P, "retirement_income_z_scale", 0.0)) != 0.0:
        raise ValueError("adopted pension rule requires equal retiree benefits")
    mass = stationary_age_income_mass(P)
    workers = float(mass[:, :int(P.J_R), :].sum())
    retirees = float(mass[:, int(P.J_R):, :].sum())
    if workers <= 0 or retirees <= 0:
        raise ValueError("positive worker and retiree populations required")
    rate = PENSION_RATIO * retirees / workers
    if not 0 < rate < 1:
        raise ValueError("derived payroll tax is outside (0,1)")
    return rate, dict(pension_to_mean_gross_worker_earnings=PENSION_RATIO,
                      worker_mass=workers, retiree_mass=retirees,
                      retiree_worker_ratio=retirees / workers, payroll_tax=rate,
                      transition_rule="hold baseline tax fixed; equal pension balances actual PAYGO")


def verify_pension_ratio(fiscal):
    """Certify the adopted benefit ratio on solved households, not seed defaults."""
    accounts = fiscal["actual_accounts"]
    mean_earnings = accounts["payroll_tax_base_period"] / accounts["worker_household_mass"]
    ratio = accounts["pension_period_units"] / mean_earnings
    if not math.isfinite(ratio) or abs(ratio - PENSION_RATIO) > 1e-9:
        raise RuntimeError("solved pension/gross-worker-earnings ratio differs from adopted target")
    return ratio


def bind_point(selected, point, arm, old, reference_rent):
    """Preserve the audited binding; replace only the explicitly varied block."""
    from e5f_parenthood_utility import bind_parenthood_utility, validate_parenthood_utility
    housing, exponent = ARMS[arm]
    names = set(COMMON) | ({"h_P"} if housing == "floor" else set(SHARE_COORDINATES))
    if set(point) != names or not all(math.isfinite(float(x)) for x in point.values()):
        raise ValueError("candidate coordinates differ from the declared arm")
    if not 0.94 <= point["beta_annual"] <= 0.99:
        raise ValueError("annual beta must retain the adopted 0.99 cap")
    base_point = {key: point[key] for key in COMMON}
    base_point["h_P"] = point.get("h_P", selected["parameters"].hbar_first_child_jump)
    P = bind_parenthood_utility(selected["parameters"], base_point)
    P.theta1 = old.THETA1
    P.delta = 1.0 - (1.0 - old.ANNUAL_DEP) ** int(P.period_years)
    P.tau_H = old.ANNUAL_PROPERTY_TAX * int(P.period_years)
    P.user_cost_rate = P.q + P.delta + P.tau_H
    validate_parenthood_utility(P)
    if housing == "shares":
        P.child_room_floor = False
        P.hbar_first_child_jump = P.hbar_child_rooms = 0.0
        for name in SHARE_COORDINATES:
            if not 0 <= float(point[name]) <= 0.25:
                raise ValueError("share coordinate outside inherited experimental bounds")
            setattr(P, name, float(point[name]))
    P.adult_entry_clock = "split_birth_vintage"
    P.utility_comparison_arm = arm
    P.utility_child_benefit_exponent = exponent
    P.utility_reference_rent = float(reference_rent)
    # Tax is derived from demographics once for this baseline and then pinned.
    return P


def configure_runtime(pair, old, lock, arm, objective_path, reference_rent):
    """Reuse nightpair's accounting/gates and change the named scientific delta.

    Call in a fresh process per arm. The frozen pair/ancestor/source remain
    untouched; runtime functions reference their existing globals.
    """
    pair.configure(old, "oasi_087510", lock)
    tax, plan, selected, _ = pair.prepare(old, lock)
    old.OBJECTIVE = Path(objective_path)
    old.OBJECTIVE_SHA = file_hash(objective_path)
    objective = json.loads(Path(objective_path).read_text())
    old.FREE = COMMON + (("h_P",) if ARMS[arm][0] == "floor" else SHARE_COORDINATES)
    old.TAX, fiscal_rule = pension_tax_from_demographics(selected["parameters"])
    old.apply_point = lambda saved, point: bind_point(saved, point, arm, old, reference_rent)
    old_setup = old.setup_runtime
    def setup(*args):
        runtime = old_setup(*args)
        model = runtime["model"]
        original = model.precompute_shared
        if getattr(original, "_four_arm_installed", False):
            raise RuntimeError("utility runtime must be installed only once per process")
        def shared(P, grid):
            return comparison_shared(original, P, grid)
        shared._four_arm_installed = True
        shared._four_arm_native_precompute = original
        model.precompute_shared = shared
        return runtime
    old.setup_runtime = setup
    original_parameters = tax.actual_parameters
    def actual(P):
        return dict(original_parameters(P), delta_alpha_jump=P.delta_alpha_jump, delta_alpha=P.delta_alpha)
    tax.actual_parameters = actual
    original_evaluate = old.evaluate_point
    def evaluate(**kwargs):
        receipt = original_evaluate(**kwargs)
        if not 1 <= int(receipt["objective_stationary_solves"]) <= 23:
            raise RuntimeError("native normalizer exceeded the frozen solve-count bound")
        receipt.update(free_count=len(old.FREE), utility_comparison_arm=arm,
                       benefit_exponent=ARMS[arm][1], utility_reference_rent=reference_rent,
                       utility_normalization="experimental compensated expenditure at fixed reference rent",
                       fiscal_rule=fiscal_rule,
                       pension_ratio_actual=verify_pension_ratio(receipt["fiscal"]),
                       estate_observer_warning="all positive estates; SCF target is child-directed")
        path = kwargs["output"] / "parameters.csv"
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            if row["parameter"] == "payroll_tax":
                row["status"] = "derived from adopted pension ratio and baseline demographics"
        rows.extend(dict(parameter=name, estimate=value, lower="", upper="", near_bound="", status=status)
                    for name, value, status in (
                        ("child_benefit_exponent", ARMS[arm][1], "fixed sensitivity; not estimated"),
                        ("utility_reference_rent", reference_rent, "fixed experimental common normalization"),
                        ("pension_to_gross_worker_earnings", PENSION_RATIO, "adopted CPS2007 target")))
        old.table(path, rows)
        old.write(kwargs["output"] / "receipt.json", receipt)
        return receipt
    old.evaluate_point = evaluate
    return tax, plan, selected, objective, fiscal_rule
