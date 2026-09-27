"""Estate-funded entry ledger; no policy, utility, or entry-law mutation.

Net positive estates after the inherited selling cost provisionally finance
positive entrant assets. The residual is a nonutility sink. This accounting
restriction does not settle financial counterparties or physical housing.
Numerical execution and synthetic tests belong on Torch.
"""
from __future__ import annotations

import math

# Capability must accompany an authenticated source contract, never replace it.
SUPPORTS_NATIVE_DUE_STAYER_CREDIT = True


class EstateFundingShortfall(RuntimeError):
    """Only this exception denotes a candidate rejected for estate funding."""

    def __init__(self, ledger):
        self.audit = ledger
        super().__init__(
            "Net estate funding shortfall: "
            f"B={ledger['available_estates_period']:.12g}, "
            f"E={ledger['entry']['positive_financial_endowment']:.12g}"
        )


def policy_mass_branches(evaluation, P):
    """Return disjoint (mass, saving, consumption) branches at current tenure.

    Under DUE the same post-tenure node can contain buyers and inherited
    owners with different saving policies; averaging policies is invalid for
    signed estates and feasibility tests. The split must be owned by this date.
    """
    import numpy as np
    g = np.asarray(evaluation.g_current, dtype=float)
    policy = evaluation.policy
    base = (g, np.asarray(policy.bp_pol), getattr(policy, "c_pol", None))
    if not bool(getattr(P, "native_due_stayer_credit", False)):
        return (base,)
    stay = getattr(evaluation, "g_stay_distribution", None)
    bp = getattr(policy, "bp_pol_stay", None)
    c = getattr(policy, "c_pol_stay", None)
    if stay is None or bp is None or c is None:
        raise ValueError("DUE accounting requires date-owned stayer mass/saving/consumption")
    stay, bp, c = (np.asarray(x, dtype=float) for x in (stay, bp, c))
    if (g.ndim != 7 or any(x.shape != g.shape for x in (stay, bp, c))
            or any(not np.isfinite(x).all() for x in (g, stay, bp, c))
            or np.any(g < 0) or np.any(stay < 0) or np.any(stay[:, 0] != 0)):
        raise ValueError("DUE accounting requires aligned finite owner-only stayer mass")
    other = g - stay
    if np.any(other < -1e-12):
        raise ValueError("DUE stayer mass exceeds current mass")
    # Only floating-point subtraction dust, never substantive mass, is removed.
    other = np.maximum(other, 0.0)
    return ((other, base[1], base[2]), (stay, bp, c))


def branch_estate_accounts(evaluation, P, death, housing_values, selling_cost):
    """Integrate each origin before applying positive/negative estate splits."""
    from audit_e5f_estate_resource_account import signed_accounts
    branches = policy_mass_branches(evaluation, P)
    result = signed_accounts(branches[0][0], branches[0][1], death, housing_values, selling_cost)
    for mass, saving, _ in branches[1:]:
        extra = signed_accounts(mass, saving, death, housing_values, selling_cost)
        for key in result["totals"]:
            result["totals"][key] += extra["totals"][key]
        for row, other in zip(result["by_age"], extra["by_age"]):
            for key in row:
                if key != "age_index":
                    row[key] += other[key]
    return result


def audit(evaluation, P, b_grid, *, next_entrant_cohort=None):
    """Return a ledger, or raise EstateFundingShortfall with its rejected ledger.

    Arrays have axes wealth, tenure, location, age, income, children ever born,
    and children at home. Estates use g_current and post-saving b'; entrants
    use the age-zero renter mass in g_pre. Both flows remain in period units.
    Without next_entrant_cohort, retain the stationary comparison unchanged.
    With it, date-t death estates fund the supplied actual date-t+1 entry cohort
    (axes wealth, tenure, location, income, children ever born, children at home).
    This ledger does not settle financial counterparties or physical housing.
    """
    import numpy as np
    from audit_e5f_estate_resource_account import signed_accounts

    g = np.asarray(evaluation.g_current, dtype=float)
    pre = np.asarray(evaluation.g_pre, dtype=float)
    bp = np.asarray(evaluation.policy.bp_pol, dtype=float)
    bg = np.asarray(b_grid, dtype=float)
    houses = np.asarray(P.H_own, dtype=float)
    prices = np.asarray(evaluation.policy.price, dtype=float).reshape(-1)
    if (g.ndim != 7 or pre.shape != g.shape or bp.shape != g.shape
            or any(size == 0 for size in g.shape)
            or int(P.I) != 1 or g.shape[2] != 1
            or g.shape[3] != int(P.J) or int(P.J) < 1
            or houses.ndim != 1 or g.shape[1] != 1 + houses.size
            or bg.shape != (g.shape[0],) or prices.shape != (1,)):
        raise ValueError("Estate audit requires aligned one-market seven-axis arrays")
    if (not np.isfinite(g).all() or not np.isfinite(pre).all()
            or np.any(g < 0) or np.any(pre < 0)
            or not np.isfinite(bp).all() or not np.isfinite(bg).all()
            or not np.isfinite(houses).all() or np.any(houses <= 0)
            or not np.isfinite(prices).all() or np.any(prices <= 0)):
        raise ValueError("Estate audit inputs must be finite with valid masses/products")
    years = float(P.period_years)
    selling_cost = float(P.psi)
    if not math.isfinite(years) or years <= 0 or not 0 <= selling_cost < 1:
        raise ValueError("Invalid period length or inherited housing selling cost")
    if not bool(getattr(P, "use_postdecision_current_distribution", True)):
        raise ValueError("Estates require post-transaction current mass")
    if (str(getattr(P, "estate_receiver", "none")) != "none"
            or bool(getattr(P, "bequest_net_of_selling_cost", False))
            or any(float(getattr(P, name, 0.0)) != 0.0 for name in
                   ("estate_lump_sum_transfer", "estate_probe_transfer", "estate_tax_rate"))):
        raise ValueError("Entry funding audit requires unchanged donor utility and no adult estate transfers")
    if bool(getattr(P, "use_age_survival", False)):
        survival = np.asarray(P.survival_probs, dtype=float)
        if (survival.shape != (int(P.J) - 1,)
                or not np.isfinite(survival).all()
                or np.any(survival < 0) or np.any(survival > 1)):
            raise ValueError("Invalid age survival schedule")
    else:
        survival = np.ones(int(P.J) - 1)
    death = np.r_[1.0 - survival, 1.0]
    pre_age = pre.sum(axis=(0, 1, 2, 4, 5, 6))
    current_age = g.sum(axis=(0, 1, 2, 4, 5, 6))
    mass_scale = max(1.0, float(pre.sum()), float(g.sum()))
    if pre.sum() <= 0 or np.max(np.abs(pre_age - current_age)) > 1e-10 * mass_scale:
        raise ValueError("Pre-fertility and current distributions must preserve positive age mass")

    estates = branch_estate_accounts(evaluation, P, death, np.r_[0.0, prices[0] * houses], selling_cost)
    dated = next_entrant_cohort is not None
    entrant = (np.asarray(next_entrant_cohort, dtype=float)
               if dated else pre[:, :, :, 0])
    if dated and (entrant.shape != pre[:, :, :, 0].shape
                  or not np.isfinite(entrant).all() or np.any(entrant < 0)):
        raise ValueError("Next entrant cohort must be aligned, finite and nonnegative")
    owner_entry_mass = float(entrant[:, 1:].sum())
    if owner_entry_mass > 1e-10 * mass_scale:
        raise ValueError("Inherited entrants must have renter tenure before transactions")
    # Do not offset positive entrant funding by negative entrant positions.
    entrant_by_asset = entrant.sum(axis=(1, 2, 3, 4, 5))
    entry_mass = float(entrant_by_asset.sum())
    if not dated and entry_mass <= 0:
        raise ValueError("Positive entrant mass is required")
    positive = float(np.dot(entrant_by_asset, np.maximum(bg, 0.0)))
    negative = float(np.dot(entrant_by_asset, np.maximum(-bg, 0.0)))
    entry = dict(mass=entry_mass, positive_financial_endowment=positive,
                 negative_financial_position=negative,
                 signed_financial_endowment=positive - negative,
                 owner_prechoice_mass=owner_entry_mass)
    totals = estates["totals"]
    gross_residual = totals["gross_positive"] - positive
    net_residual = totals["net_positive"] - positive
    tolerance = 1e-10 * max(1.0, totals["net_positive"], positive)
    funded = net_residual >= -tolerance
    ledger = dict(
        audit_id=("estate_funded_dated_entry_provisional_net_v1" if dated
                  else "estate_funded_entry_provisional_net_v1"),
        status="funded" if funded else "funding_shortfall",
        period_years=years, units="model financial units per model period",
        estate=estates, entry=entry,
        available_estates_period=totals["net_positive"],
        available_estates_valuation="provisional_positive_net_estates_after_inherited_selling_cost",
        inherited_selling_cost=selling_cost,
        gross_residual_after_entry_period=gross_residual,
        net_residual_after_entry_period=net_residual,
        residual_sink_period=max(net_residual, 0.0) if funded else None,
        funding_shortfall_period=max(-net_residual, 0.0),
        funding_gate_tolerance=tolerance,
        annual_gross_positive_estates=totals["gross_positive"] / years,
        annual_net_positive_estates=totals["net_positive"] / years,
        annual_positive_entrant_funding=positive / years,
        entrant_minus_death_mass=entry_mass - totals["death_mass"],
        timing=("date-t death flow finances supplied actual date-t+1 entrant cohort" if dated
                else "stationary death flow finances the next entrant cohort; no transition implementation"),
        policy_changes=False, adult_transfers=False, entry_distribution_changes=False,
        donor_bequest_valuation_changes=False,
        certifies_counterparty_or_physical_housing_settlement=False,
        caveats=["Availability valuation is provisional; donor utility remains unchanged.",
                 "Negative estates remain explicit liabilities with an unresolved creditor rule.",
                 "Negative entrant positions do not reduce positive entrant funding.",
                 "This is not a child-directed bequest-flow measurement correction.",
                 "An admissible ledger alone does not certify stationary equilibrium."],
    )
    if not funded:
        raise EstateFundingShortfall(ledger)
    return ledger


def run_self_tests():
    """Small fixtures on Torch or under explicit author-authorized local testing."""
    import os
    import sys
    if not os.environ.get('ALLOW_LOCAL_RUNTIME_TESTS') and (sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID")):
        raise RuntimeError("Run estate audit synthetic tests only in a Torch allocation")
    from types import SimpleNamespace
    import numpy as np

    g = np.zeros((2, 2, 1, 2, 1, 1, 1))
    bp = np.zeros_like(g)
    g[0, 0, 0, 0, 0, 0, 0], g[1, 0, 0, 0, 0, 0, 0] = .4, .6
    g[0, 0, 0, 1, 0, 0, 0] = .25
    g[0, 1, 0, 1, 0, 0, 0] = .5
    g[1, 1, 0, 1, 0, 0, 0] = .25
    bp[0, 0, 0, 1, 0, 0, 0] = -2.
    bp[0, 1, 0, 1, 0, 0, 0] = -9.
    bp[1, 1, 0, 1, 0, 0, 0] = 5.
    P = SimpleNamespace(I=1, J=2, H_own=np.array([10.]), period_years=4.,
                        psi=.2, use_age_survival=True, survival_probs=np.array([1.]))
    evaluation = SimpleNamespace(g_current=g, g_pre=g.copy(),
                                 policy=SimpleNamespace(bp_pol=bp, price=np.array([1.])))
    original = (g.copy(), bp.copy(), evaluation.g_pre.copy())
    result = audit(evaluation, P, np.array([-2., 5.]))
    expected = dict(gross_positive=4.25, gross_negative=.5,
                    net_positive=3.25, net_negative=1., housing_sale_cost=1.5,
                    death_mass=1., negative_estate_death_mass=.75)
    errors = [abs(result["estate"]["totals"][k] - v) for k, v in expected.items()]
    errors += [abs(result["entry"]["positive_financial_endowment"] - 3.),
               abs(result["entry"]["negative_financial_position"] - .8),
               abs(result["net_residual_after_entry_period"] - .25),
               abs(result["annual_net_positive_estates"] - .8125)]
    assert max(errors) < 1e-12
    try:
        audit(evaluation, P, np.array([-2., 6.]))
    except EstateFundingShortfall as exc:
        assert abs(exc.audit["net_residual_after_entry_period"] + .35) < 1e-12
        assert exc.audit["residual_sink_period"] is None
    else:
        raise AssertionError("Funding shortfall was not rejected")
    evaluation.g_pre[0, 0, 0, 0, 0, 0, 0] -= .1
    evaluation.g_pre[0, 1, 0, 0, 0, 0, 0] += .1
    try:
        audit(evaluation, P, np.array([-2., 5.]))
    except ValueError:
        pass
    else:
        raise AssertionError("Owner entrants were not rejected as an invalid input")
    evaluation.g_pre = original[2].copy()
    P.estate_lump_sum_transfer = .1
    try:
        audit(evaluation, P, np.array([-2., 5.]))
    except ValueError:
        pass
    else:
        raise AssertionError("Adult estate transfer was not rejected as an invalid input")
    del P.estate_lump_sum_transfer
    assert np.array_equal(g, original[0]) and np.array_equal(bp, original[1])
    assert np.array_equal(evaluation.g_pre, original[2])
    return dict(status="passed", fixtures=4, maximum_absolute_error=max(errors),
                covers=["signed estates and sale costs", "positive entrant funding without debt offset",
                        "dedicated funding rejection", "invalid owner entry", "invalid adult transfer",
                        "period annualization", "no input mutation"])
