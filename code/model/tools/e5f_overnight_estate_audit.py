"""Estate-funded entry ledger; no policy, utility, or entry-law mutation.

Net positive estates after the inherited selling cost provisionally finance
positive entrant assets. The residual is a nonutility sink. This accounting
restriction does not settle financial counterparties or physical housing.
Numerical execution and synthetic tests belong on Torch.
"""
from __future__ import annotations

import math


class EstateFundingShortfall(RuntimeError):
    """Only this exception denotes a candidate rejected for estate funding."""

    def __init__(self, ledger):
        self.audit = ledger
        super().__init__(
            "Net estate funding shortfall: "
            f"B={ledger['available_estates_period']:.12g}, "
            f"E={ledger['entry']['positive_financial_endowment']:.12g}"
        )


def audit(evaluation, P, b_grid):
    """Return a ledger, or raise EstateFundingShortfall with its rejected ledger.

    Arrays have axes wealth, tenure, location, age, income, children ever born,
    and children at home. Estates use g_current and post-saving b'; entrants
    use the age-zero renter mass in g_pre. Both flows remain in period units.
    The stationary comparison pairs death funding with the next entry cohort;
    it does not implement a dated transition ledger or intraperiod settlement.
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

    estates = signed_accounts(g, bp, death, np.r_[0.0, prices[0] * houses], selling_cost)
    entrant = pre[:, :, :, 0]
    owner_entry_mass = float(entrant[:, 1:].sum())
    if owner_entry_mass > 1e-10 * mass_scale:
        raise ValueError("Inherited entrants must have renter tenure before transactions")
    # Do not offset positive entrant funding by negative entrant positions.
    entrant_by_asset = entrant.sum(axis=(1, 2, 3, 4, 5))
    entry_mass = float(entrant_by_asset.sum())
    if entry_mass <= 0:
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
        audit_id="estate_funded_entry_provisional_net_v1",
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
        timing="stationary death flow finances the next entrant cohort; no transition implementation",
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
    """Small hand-calculated fixtures; executable only within a Torch allocation."""
    import os
    import sys
    if sys.platform != "linux" or not os.environ.get("SLURM_JOB_ID"):
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
