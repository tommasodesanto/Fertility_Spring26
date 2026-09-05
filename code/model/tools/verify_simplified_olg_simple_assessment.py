#!/usr/bin/env python3
"""Construct and check small analytical examples for the simplicity assessment.

No calibration or equilibrium search. Equilibria are constructed from chosen
owner bundles, then checked against the original dated problems, tenure rule,
demography, housing clearing and fiscal budget. Independent optimizations use
all seven dated household choices, without the reduced value functions.

Run from any directory:
    python3 code/model/tools/verify_simplified_olg_simple_assessment.py --plot
Outputs stay in the existing simplified_olg_amendments evidence folder.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, brentq, minimize
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def utility(z, p):
    x, h, n, saving, c2, h2, estate = z
    space = h - p["kappa"] * n
    if min(x, space, n, c2, h2, estate) <= 0:
        return -1e12
    return (math.log(x) + p["alpha"] * math.log(space)
            + p["theta"] * math.log(n)
            + p["beta"] * (math.log(c2) + p["gamma"] * math.log(h2)
                           + p["omega"] * math.log(estate)))


def constructed_equilibrium(tau=0.05, rental_cap_binds=True):
    p = dict(q=.5, phi=.8, alpha=.4, beta=.4, gamma=.3, omega=.4,
             chi=.15, kappa=.5, sigma=.12, nu=2., Hbar=2., P=1., b=.2,
             owner_cap=2., tau=tau)
    p["rental_cap"] = .25 if rental_cap_binds else 1.5
    q, u = p["q"], .5 * (1 + tau)
    rho = 1 + p["beta"] * (1 + p["gamma"] + p["omega"])
    rho_r = 1 + p["beta"] * (1 + p["omega"])
    x, h, n = 1., 1., .75 if rental_cap_binds else .49
    p["theta"] = n * (p["chi"] / x + p["alpha"] * p["kappa"] / (h-p["kappa"]*n))
    w = rho*x + p["chi"]*n + u*h
    if rental_cap_binds:
        hr = p["rental_cap"]
        available = w - (1+q)*u*hr

        def foc(nr):
            xr = (available-p["chi"]*nr) / rho_r
            return (p["theta"]/nr - p["chi"]/xr
                    - p["alpha"]*p["kappa"]/(hr-p["kappa"]*nr))

        nr = brentq(foc, 1e-12, hr/p["kappa"]-1e-12, xtol=1e-14)
        xr = (available-p["chi"]*nr) / rho_r
        hor = hr
    else:
        xr = w / (rho+p["alpha"]+p["theta"])
        nr = p["theta"]*xr / (p["chi"]+p["kappa"]*u)
        hr = p["alpha"]*xr/u+p["kappa"]*nr
        hor = p["gamma"]*p["beta"]*xr/(q*u)
    ho = p["gamma"]*p["beta"]*x/(q*u)
    pi = (1/p["nu"]-nr)/(n-nr)
    assert 0 < pi < 1
    housing_per_pair = pi*(h+ho)+(1-pi)*(hr+hor)
    cohort = p["Hbar"]/housing_per_pair
    T = q*tau*p["P"]*housing_per_pair/2
    p["y"] = w-p["b"]-(1+q)*T
    c2, c2r = p["beta"]*x/q, p["beta"]*xr/q
    e, er = p["omega"]*c2/q, p["omega"]*c2r/q
    saving = p["y"]+p["b"]+T-x-p["chi"]*n-((1-p["phi"])+q*tau)*p["P"]*h
    saving_r = p["y"]+p["b"]+T-xr-p["chi"]*nr-u*hr
    owner = [x, h, n, saving, c2, ho, e]
    renter = [xr, hr, nr, saving_r, c2r, hor, er]
    wo, wr = utility(owner, p), utility(renter, p)
    p["xi_bar"] = p["sigma"]*math.log(pi/(1-pi))+wr-wo
    result = dict(parameters=p, user_cost=u, transfer=T, cohort_mass=cohort,
                  owner_share=pi, owner=owner, renter=renter,
                  rental_cap_binds=rental_cap_binds,
                  owner_value=wo, renter_value=wr)
    return result


def dated_constraints(case, owner):
    """Original linear budgets and physical/financial restrictions, in seven choices."""
    p = case["parameters"]
    q, u, P, T = p["q"], case["user_cost"], p["P"], case["transfer"]
    rows = np.zeros((2, 7))
    rows[0, [0, 2, 3]] = [1, p["chi"], 1]
    rows[1, [3, 4, 5, 6]] = [-1/q, 1, u, q]
    if owner:
        rows[0, 1] = ((1-p["phi"])+q*p["tau"])*P
        rows[1, 1] = p["phi"]*P/q-P
    else:
        rows[0, 1] = u
    rhs = np.array([p["y"]+p["b"]+T, T])
    ineq = [np.array([0., 1., -p["kappa"], 0., 0., 0., 0.])]
    lower, upper = [1e-10], [np.inf]
    if owner:
        ineq += [np.array([0., (1-p["phi"])*P, 0., 0., 0., 0., 0.]),
                 np.array([0., 1., 0., 0., 0., -1., 0.]),
                 np.array([0., 0., 0., 0., 0., -P, 1.])]
        lower += [-np.inf, 0., 0.]
        upper += [p["b"], np.inf, np.inf]
    cap = p["owner_cap"] if owner else p["rental_cap"]
    bounds = Bounds([1e-10, 1e-10, 1e-10, 0., 1e-10, 1e-10, 1e-10],
                    [np.inf, cap, np.inf, np.inf, np.inf,
                     np.inf if owner else cap, np.inf])
    return rows, rhs, np.asarray(ineq), np.array(lower), np.array(upper), bounds


def check_equilibrium(case):
    p, u, T, Y = (case[k] for k in ("parameters", "user_cost", "transfer", "cohort_mass"))
    zo, zr = np.array(case["owner"]), np.array(case["renter"])
    pi = case["owner_share"]
    residuals = []
    optimizer_checks = []
    for name, z, is_owner in (("owner", zo, True), ("renter", zr, False)):
        M, rhs, G, low, high, bounds = dated_constraints(case, is_owner)
        residuals.extend((M@z-rhs).tolist())
        assert np.all(G@z >= low-1e-12) and np.all(G@z <= high+1e-12)
        assert np.all(z >= bounds.lb-1e-12) and np.all(z <= bounds.ub+1e-12)
        residuals.append(p["theta"]/z[2]-p["chi"]/z[0]
                         -p["alpha"]*p["kappa"]/(z[1]-p["kappa"]*z[2]))
        residuals.append(1/z[0]-p["beta"]/(p["q"]*z[4]))

        # Independent optimum: original current and old budgets are imposed together.
        initial = z*np.array([.97, 1.01, .98, 1.02, 1.01, .96, 1.02])
        initial = np.minimum(np.maximum(initial, bounds.lb*2), bounds.ub)
        fit = minimize(lambda zz: -utility(zz, p), initial, method="SLSQP",
                       bounds=bounds,
                       constraints=[LinearConstraint(M, rhs, rhs),
                                    LinearConstraint(G, low, high)],
                       options={"ftol": 1e-12, "maxiter": 600})
        error = float(np.max(np.abs(fit.x-z)))
        budget_error = float(np.max(np.abs(M@fit.x-rhs)))
        assert fit.success, fit.message
        assert error < 2e-5 and budget_error < 1e-9, (name, error, budget_error)
        assert abs(utility(fit.x, p)-utility(z, p)) < 1e-9
        optimizer_checks.append(dict(tenure=name, maximum_choice_error=error,
                                     maximum_budget_error=budget_error))
    residuals.extend([
        pi*zo[2]+(1-pi)*zr[2]-1/p["nu"],
        Y*(pi*(zo[1]+zo[5])+(1-pi)*(zr[1]+zr[5]))-p["Hbar"],
        2*Y*T-p["q"]*p["tau"]*p["P"]*p["Hbar"],
        expit((case["owner_value"]+p["xi_bar"]-case["renter_value"])/p["sigma"])-pi,
        u-((1+p["q"]*p["tau"])*p["P"]-p["q"]*p["P"])
    ])
    mv_y = p["alpha"]*zo[0]/(zo[1]-p["kappa"]*zo[2])
    mv_o = p["gamma"]*zo[4]/zo[5]
    assert mv_y > u and abs(mv_o-u) < 1e-12
    assert zo[3] > 0 and zo[1] < p["owner_cap"]
    assert zo[5] < zo[1] and zo[6] > p["P"]*zo[5]
    if not case["rental_cap_binds"]:
        assert max(zr[1], zr[5]) < p["rental_cap"]
    assert max(abs(v) for v in residuals) < 1e-11
    return dict(maximum_original_equilibrium_residual=max(abs(v) for v in residuals),
                young_housing_value=mv_y, old_housing_value=mv_o,
                optimizer_checks=optimizer_checks)


def check_direct_reallocation(case):
    p, u, T = case["parameters"], case["user_cost"], case["transfer"]
    x, h, n, saving, c2, h2, estate = case["owner"]
    q, P, phi = p["q"], p["P"], p["phi"]
    space = h-p["kappa"]*n
    gap = p["alpha"]*x/space-p["gamma"]*c2/h2
    rows = []
    for eps in (1e-6, 1e-4, 1e-3):
        for ratio in (0., .5, 1.5):
            K = ratio*gap
            D = c2*math.expm1(-p["gamma"]*math.log1p(-eps/h2))
            C = K*eps
            L = D-u*eps
            sn = saving+phi*P*eps-q*P*eps
            young_resources = p["y"]+p["b"]+T
            old_assets = saving/q-phi*P*h/q
            old_assets_new = sn/q-phi*P*(h+eps)/q
            residuals = [
                x+p["chi"]*n-D-C+sn+((1-phi)+q*p["tau"])*P*(h+eps)
                -(young_resources-L-C),
                old_assets_new+P*(h+eps)-old_assets-P*h,
                c2+D+q*estate+u*(h2-eps)-(old_assets+P*h+T+L),
                (x-D-C)+(c2+D)+C-x-c2,
                (h+eps)+(h2-eps)-h-h2,
                math.log1p(D/c2)+p["gamma"]*math.log1p(-eps/h2),
                (sn-saving)+q*P*eps-phi*P*eps
            ]
            gain = math.log1p(-(D+C)/x)+p["alpha"]*math.log1p(eps/space)
            assert max(abs(v) for v in residuals) < 1e-12
            assert (gain > 0) == (ratio < 1)
            assert x-D-C > 0 and sn > 0 and h+eps < p["owner_cap"]
            assert estate >= P*(h2-eps) and h2 < h+eps
            rows.append(dict(epsilon=eps, marginal_cost=K, young_gain=gain,
                             old_gain=residuals[5], maximum_accounting_error=max(abs(v) for v in residuals)))
    return rows


def check_renter_replication(case):
    """Remove only the rental cap from this fixed-price deviation, not from GE."""
    p, u = case["parameters"], case["user_cost"]
    z = np.array(case["owner"], dtype=float)
    z[3] += p["q"]*p["P"]*z[1]-p["phi"]*p["P"]*z[1]
    assert z[3] >= 0
    M, rhs, *_ = dated_constraints(case, False)
    err = float(np.max(abs(M@z-rhs)))
    assert err < 1e-12
    assert abs(utility(z, p)-case["owner_value"]) < 1e-12
    eps = 1e-5
    better = z.copy()
    better[0] -= u*eps
    better[1] += eps
    assert float(np.max(abs(M@better-rhs))) < 1e-12
    gain = utility(better, p)-utility(z, p)
    assert gain > 0
    return dict(replication_budget_error=err, replicated_renter_saving=float(z[3]),
                utility_gain_from_more_rental_space=gain,
                scope="Fixed-price deviation with adequate rental space and no ownership taste; not a cap-removal policy equilibrium.")


def check_fertility():
    cases = [dict(rho=2., alpha=1., chi=1., kappa=.1, theta=1.1,
                  x=1., s=1., n=1., u=.5, p=.5),
             dict(rho=2., alpha=1., chi=1., kappa=.1, theta=1.1,
                  x=1., s=1., n=1., u=.5, p=0.)]
    output = []
    for z in cases:
        rho, alpha, chi, kappa, theta, x, s, n, u, payment = (z[k] for k in
            ("rho", "alpha", "chi", "kappa", "theta", "x", "s", "n", "u", "p"))
        h = s+kappa*n
        w = rho*x+chi*n+u*h
        def fertility(dh):
            hn = h+dh
            wn = w+(u-payment)*dh
            def foc(nn):
                xx = (wn-u*hn-chi*nn)/rho
                return theta/nn-chi/xx-alpha*kappa/(hn-kappa*nn)
            upper = min(hn/kappa, (wn-u*hn)/chi)-1e-10
            return brentq(foc, 1e-10, upper, xtol=1e-14)
        eps = 1e-5
        numerical = (fertility(eps)-fertility(-eps))/(2*eps)
        numerator = alpha*kappa/s**2-chi*payment/(rho*x**2)
        denominator = theta/n**2+chi**2/(rho*x**2)+alpha*kappa**2/s**2
        analytical = numerator/denominator
        assert abs(numerical-analytical) < 1e-8
        output.append(dict(payment=payment, predicted_derivative=analytical,
                           numerical_derivative=numerical))
    return output


def check_original_negative_fertility_example():
    """Embed the negative derivative in the original dated owner problem."""
    p = dict(q=.5, phi=.8, alpha=1., beta=.4, gamma=.5, omega=1.,
             chi=1., kappa=.1, theta=1.1, P=1., b=.22, y=3.33,
             owner_cap=2., rental_cap=.25, tau=0.)
    case = dict(parameters=p, user_cost=.5, transfer=0.)
    z = np.array([1., 1.1, 1., 1.33, .8, .8, 1.6])
    M, rhs, G, low, high, bounds = dated_constraints(case, True)
    assert np.max(abs(M@z-rhs)) < 1e-12
    assert np.all(G@z >= low-1e-12) and np.all(G@z <= high+1e-12)
    initial = z*np.array([.98, .99, 1.01, 1.02, .98, 1.02, .99])
    fit = minimize(lambda zz: -utility(zz, p), initial, method="SLSQP",
                   bounds=bounds,
                   constraints=[LinearConstraint(M, rhs, rhs),
                                LinearConstraint(G, low, high)],
                   options={"ftol": 1e-12, "maxiter": 600})
    error = float(np.max(abs(fit.x-z)))
    assert fit.success and error < 2e-5, (fit.message, error)
    assert abs(utility(fit.x, p)-utility(z, p)) < 1e-9
    assert z[3] > 0 and z[5] < z[1] and z[6] > p["P"]*z[5]
    assert p["alpha"]*z[0]/(z[1]-p["kappa"]*z[2]) > case["user_cost"]
    return dict(parameters=p, owner=z.tolist(), maximum_choice_error=error,
                scope="An optimal constrained owner at fixed prices, not a constructed policy equilibrium.")


def demographic_path():
    """An exact accounting illustration, with prescribed fertility, not a GE solution."""
    t = np.arange(31)
    stock_ratio, decay = 1.2, .6
    log_fertility_ratio = (1-decay)*math.log(stock_ratio)*decay**t
    young = np.exp(math.log(stock_ratio)*(1-decay**t))
    old = np.r_[1., young[:-1]]
    total_households = (young+old)/2
    err = float(np.max(np.abs(young[1:]/young[:-1]-np.exp(log_fertility_ratio[:-1]))))
    assert err < 1e-12
    return t, log_fertility_ratio, total_households, stock_ratio, err


def make_figures(case):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                         "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    p = case["parameters"]
    x, h, n, _, c2, h2, _ = case["owner"]
    s = h-p["kappa"]*n
    def values(eps):
        compensation = c2*np.expm1(-p["gamma"]*np.log1p(-eps/h2))
        return (p["alpha"]*(x-compensation)/(s+eps),
                p["gamma"]*(c2+compensation)/(h2-eps))
    crossing = brentq(lambda e: values(e)[0]-values(e)[1], 0., h2*.4)
    eps = np.linspace(0., crossing*1.9, 301)
    vy, vo = values(eps)
    fig, ax = plt.subplots(figsize=(6.5, 3.0), layout="constrained")
    ax.plot(eps, vy, color="#245578", lw=2, label="Young owner's value")
    ax.plot(eps, vo, color="#a04b37", lw=2, label="Old owner's compensation cost")
    ax.fill_between(eps, vy, vo, where=eps<=crossing, color="#c6dce8", alpha=.7)
    ax.axvline(crossing, color=".5", ls="--", lw=.9)
    ax.scatter([0., 0.], [vy[0], vo[0]], s=25, c=["#245578", "#a04b37"], zorder=5)
    ax.set_xlim(-crossing*.025, eps[-1])
    ax.set_xlabel(r"Housing moved from the old owner to the young owner, $\epsilon$")
    ax.set_ylabel("Consumption per housing unit")
    ax.set_xticks([0., crossing], ["Competitive\nallocation", "Values equal\nfor this pair"])
    ax.legend(frameon=False, fontsize=9, loc="lower center", bbox_to_anchor=(.5, 1.01))
    for ext in ("pdf", "png"):
        fig.savefig(OUT/f"simple_allocation_figure.{ext}", dpi=220)
    plt.close(fig)

    t, log_fertility_ratio, households, scale, _ = demographic_path()
    use = t <= 13
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.35), layout="constrained")
    axs[0].plot(t[use], np.exp(log_fertility_ratio[use]), color="#245578", lw=2)
    axs[0].axhline(1., color=".5", ls="--", lw=1)
    axs[0].set_title("A temporary fertility increase", fontsize=10)
    axs[0].set_ylabel("Fertility relative to baseline")
    axs[0].set_yticks([1.], ["Replacement"])
    axs[1].plot(t[use], households[use], color="#245578", lw=2)
    axs[1].axhline(scale, color=".5", ls="--", lw=1)
    axs[1].set_title("A permanently larger population", fontsize=10)
    axs[1].set_ylabel("Young and old households")
    axs[1].set_yticks([1., scale], [r"$N^*_0$", r"$N^*_1$"])
    for ax in axs:
        ax.set_xlabel("Generations after the change")
        ax.set_xticks([0, 4, 8, 12])
        ax.set_xlim(0, 13)
    for ext in ("pdf", "png"):
        fig.savefig(OUT/f"simple_population_figure.{ext}", dpi=220)
    plt.close(fig)
    return dict(compensated_pair_optimum=crossing,
                population_figure="Conditional demographic identity with prescribed fertility; no equilibrium transition was solved.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cases = [constructed_equilibrium(tau) for tau in (0., .05, .27)]
    cases.append(constructed_equilibrium(0., rental_cap_binds=False))
    report = {"scope": "Analytical equilibrium witnesses, original-problem checks, and conditional figures; no calibration or policy simulation.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "specification_source": "latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex",
              "specification_sha256": hashlib.sha256((ROOT / "latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex").read_bytes()).hexdigest(),
              "equilibria": []}
    for case in cases:
        checks = check_equilibrium(case)
        checks["direct_reallocations"] = check_direct_reallocation(case)
        checks["renter_replication"] = check_renter_replication(case)
        report["equilibria"].append(dict(**case, checks=checks))
    report["fertility_checks"] = check_fertility()
    report["original_negative_fertility_example"] = check_original_negative_fertility_example()
    report["demographic_identity_error"] = demographic_path()[-1]
    # Stationary housing supply scaling: every per-household choice and price stays fixed.
    base = cases[1]
    scale = 1.2
    p = base["parameters"]
    scaled_Y, scaled_H = scale*base["cohort_mass"], scale*p["Hbar"]
    pi = base["owner_share"]
    z, zr = base["owner"], base["renter"]
    per_pair = pi*(z[1]+z[5])+(1-pi)*(zr[1]+zr[5])
    scaled_residual = max(abs(scaled_Y*per_pair-scaled_H),
                          abs(2*scaled_Y*base["transfer"]-p["q"]*p["tau"]*p["P"]*scaled_H))
    assert scaled_residual < 1e-12
    report["stationary_stock_scaling"] = dict(scale=scale, maximum_equilibrium_residual=scaled_residual,
        scope="Corresponding stationary allocations only; no claim of transition existence or convergence.")
    if args.plot:
        report["figures"] = make_figures(base)
    (OUT/"simple_assessment_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"equilibrium_witnesses": len(cases),
                      "independent_original_optimizations": 2*len(cases)+1,
                      "finite_reallocations": 9*len(cases),
                      "maximum_equilibrium_residual": max(c["checks"]["maximum_original_equilibrium_residual"] for c in report["equilibria"]),
                      "maximum_optimizer_choice_error": max(report["original_negative_fertility_example"]["maximum_choice_error"],
                          max(v["maximum_choice_error"] for c in report["equilibria"] for v in c["checks"]["optimizer_checks"])),
                      "output": str(OUT/"simple_assessment_checks.json")}, indent=2))


if __name__ == "__main__":
    main()
