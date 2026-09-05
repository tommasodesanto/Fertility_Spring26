#!/usr/bin/env python3
"""Check the analytical local OLG transition against the original dated model.

This is a bounded theory check, not calibration or a general transition solver.
The analytical cubic is derived symbolically. Finite paths check its derivative
and original budgets; their terminal closure does not prove infinite convergence.
Run --smoke first, then the eight declared cases. Use --stationary-only for the
positive-child-cost population condition, and --figure for the analytical plot.
Use --heterogeneity-only for uniform type constraints and exact aggregation.
Use --general-only for the general root identities and the negative initial response.
All outputs stay in the existing simple-theory evidence folder. The numerical
paths never supply the plotted curves; these come from the analytical formula.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize, root
import sympy as sp

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def parameters(**overrides):
    p = dict(q=.8, phi=.8, b=.2, y=.95, alpha=.4, beta=.4, gamma=.3,
             omega=2., kappa=.5, theta=2/15, nu=2., chi=0., tau=0.,
             rental_cap=.25, owner_cap=2., sigma=1., taste_weight=0.,
             Hbar=1213/928)
    p.update(overrides)
    return p


def household_utility(z, p):
    x, h, n, saving, c2, h2, estate = z
    return (np.log(x)+p["alpha"]*np.log(h-p["kappa"]*n)
            +p["theta"]*np.log(n)
            +p["beta"]*(np.log(c2)+p["gamma"]*np.log(h2)
                          +p["omega"]*np.log(estate)))


def capped_fertility(h, disposable, rho, p):
    """Stable small root of the original fertility FOC; smooth at chi=0."""
    chi, kappa, alpha, theta = (p[k] for k in ("chi", "kappa", "alpha", "theta"))
    middle = kappa*disposable*(theta+alpha)+chi*h*(theta+rho)
    constant = theta*disposable*h
    quadratic = chi*kappa*(theta+rho+alpha)
    return 2*constant/(middle+np.sqrt(middle*middle-4*quadratic*constant))


def young_choices(prices, transfers, p):
    P0, P1, P2 = prices
    T0, T1 = transfers
    q, tau = p["q"], p["tau"]
    u0 = (1+q*tau)*P0-q*P1
    u1 = (1+q*tau)*P1-q*P2
    w = p["y"]+p["b"]+T0+q*T1
    output = {}
    for tenure in ("owner", "renter"):
        own = tenure == "owner"
        h = p["b"]/((1-p["phi"])*P0) if own else p["rental_cap"]
        rho = 1+p["beta"]*(1+p["omega"]+(p["gamma"] if own else 0))
        available = w-u0*h-(0 if own else q*u1*p["rental_cap"])
        n = capped_fertility(h, available, rho, p)
        x = (available-p["chi"]*n)/rho
        c2 = p["beta"]*x/q
        h2 = p["gamma"]*c2/u1 if own else p["rental_cap"]
        estate = p["omega"]*c2/q
        current_cost = ((1-p["phi"])+q*tau)*P0*h if own else u0*h
        saving = p["y"]+p["b"]+T0-x-p["chi"]*n-current_cost
        z = np.array([x, h, n, saving, c2, h2, estate])
        assets = saving/q-(p["phi"]*P0*h/q if own else 0)
        output[tenure] = dict(z=z, assets=assets, utility=household_utility(z, p))
    pi = 1/(1+p["taste_weight"]*np.exp(
        (output["renter"]["utility"]-output["owner"]["utility"])/p["sigma"]))
    output.update(pi=pi, u0=u0, u1=u1,
                  fertility=pi*output["owner"]["z"][2]+(1-pi)*output["renter"]["z"][2],
                  housing=pi*output["owner"]["z"][1]+(1-pi)*output["renter"]["z"][1])
    return output


def complex_jacobian(fun, point):
    size = np.asarray(fun(point)).size
    matrix = np.empty((size, len(point)))
    for j in range(len(point)):
        perturbed = np.asarray(point, dtype=complex)
        perturbed[j] += 1j*1e-25
        matrix[:, j] = np.imag(fun(perturbed))/1e-25
    return matrix


def steady_state(p):
    def equations(v):
        P, T = v
        hh = young_choices([P, P, P], [T, T], p)
        old_h = hh["pi"]*hh["owner"]["z"][5]+(1-hh["pi"])*hh["renter"]["z"][5]
        return np.array([p["nu"]*hh["fertility"]-1,
                         T-.5*p["q"]*p["tau"]*P*(hh["housing"]+old_h)])
    guess = np.array([p["b"]/(1-p["phi"]), .5*p["q"]*p["tau"]*p["Hbar"]])
    fit = root(equations, guess, jac=lambda v: complex_jacobian(equations, v), tol=1e-11)
    assert np.max(abs(equations(fit.x))) < 1e-10, fit.message
    P, T = fit.x
    hh = young_choices([P, P, P], [T, T], p)
    old_h = hh["pi"]*hh["owner"]["z"][5]+(1-hh["pi"])*hh["renter"]["z"][5]
    return dict(price=float(P), transfer=float(T), cohort=float(p["Hbar"]/(hh["housing"]+old_h)),
                households=hh, parameters=p)


def old_choices(assets, inherited, price, next_price, transfer, p, owner):
    q = p["q"]
    u = (1+q*p["tau"])*price-q*next_price
    if owner:
        c2 = (assets+price*inherited+transfer)/(1+p["gamma"]+p["omega"])
        h2 = p["gamma"]*c2/u
    else:
        h2 = p["rental_cap"]
        c2 = (assets+transfer-u*h2)/(1+p["omega"])
    return np.array([c2, h2, p["omega"]*c2/q])


def transition(initial, final, horizon):
    p = final["parameters"]
    def equations(vector, return_rows=False):
        price = np.r_[vector[:horizon+1], final["price"], final["price"]]
        transfer = np.r_[vector[horizon+1:], final["transfer"]]
        Y = initial["cohort"]
        O = Y
        past = initial["households"]
        market, fiscal, rows = [], [], []
        for t in range(horizon+1):
            hh = young_choices(price[t:t+3], transfer[t:t+2], p)
            old = {tenure: old_choices(past[tenure]["assets"], past[tenure]["z"][1],
                                      price[t], price[t+1], transfer[t], p, tenure=="owner")
                   for tenure in ("owner", "renter")}
            old_h = past["pi"]*old["owner"][1]+(1-past["pi"])*old["renter"][1]
            market.append(Y*hh["housing"]+O*old_h-p["Hbar"])
            fiscal.append(transfer[t]-p["q"]*p["tau"]*price[t]*p["Hbar"]/(Y+O))
            if return_rows:
                rows.append(dict(t=t, price=price[t], transfer=transfer[t], young=Y, old=O,
                                 hh=hh, old_choices=old, prices=price[t:t+3], transfers=transfer[t:t+2],
                                 past=past, market_residual=market[-1], fiscal_residual=fiscal[-1]))
            O, Y = Y, p["nu"]*hh["fertility"]*Y
            past = hh
        return rows if return_rows else np.r_[market, fiscal]
    time_index = np.arange(horizon+1)
    price_guess = final["price"]+(initial["price"]-final["price"])*.25**time_index
    guess = np.r_[price_guess, np.full(horizon+1, final["transfer"])]
    fit = root(equations, guess, jac=lambda v: complex_jacobian(equations, v), tol=1e-10)
    residual = float(np.max(abs(equations(fit.x))))
    assert residual < 1e-9, (fit.message, residual)
    return dict(rows=equations(fit.x, True), maximum_residual=residual,
                function_evaluations=int(fit.nfev), horizon=horizon)


def original_constraints(row, p, owner):
    P0, P1, P2 = row["prices"]
    T0, T1 = row["transfers"]
    q = p["q"]
    u0, u1 = row["hh"]["u0"], row["hh"]["u1"]
    M = np.zeros((2, 7))
    M[0, [0, 2, 3]] = [1, p["chi"], 1]
    M[1, [3, 4, 5, 6]] = [-1/q, 1, u1, q]
    M[0, 1] = ((1-p["phi"])+q*p["tau"])*P0 if owner else u0
    if owner:
        M[1, 1] = p["phi"]*P0/q-P1
    rhs = [p["y"]+p["b"]+T0, T1]
    G = [[0., 1., -p["kappa"], 0., 0., 0., 0.]]
    lo, hi = [1e-10], [np.inf]
    if owner:
        G += [[0., (1-p["phi"])*P0, 0., 0., 0., 0., 0.],
              [0., 1., 0., 0., 0., -1., 0.],
              [0., 0., 0., 0., 0., -P2, 1.]]
        lo += [-np.inf, 0., 0.]
        hi += [p["b"], np.inf, np.inf]
    cap = p["owner_cap"] if owner else p["rental_cap"]
    bounds = Bounds([1e-10]*3+[0.]+[1e-10]*3,
                    [np.inf, cap, np.inf, np.inf, np.inf,
                     np.inf if owner else cap, np.inf])
    return M, np.array(rhs), np.array(G), np.array(lo), np.array(hi), bounds


def check_path(path, p, optimize_dates=()):
    max_budget, max_foc = 0., 0.
    minimum_slack = dict(saving=np.inf, retention=np.inf, estate=np.inf, owner_cap=np.inf,
                         owner_housing_wedge=np.inf, renter_young_wedge=np.inf, renter_old_wedge=np.inf,
                         actual_old_renter_wedge=np.inf)
    optimizations = []
    for row in path["rows"]:
        hh, u = row["hh"], row["hh"]["u0"]
        assert u>0 and hh["u1"]>0 and 0<=hh["pi"]<=1
        for tenure in ("owner", "renter"):
            owner = tenure == "owner"
            z = hh[tenure]["z"]
            M, rhs, G, lo, hi, bounds = original_constraints(row, p, owner)
            budget = float(np.max(abs(M@z-rhs)))
            max_budget = max(max_budget, budget)
            assert np.all(G@z>=lo-1e-9) and np.all(G@z<=hi+1e-9)
            assert np.all(z>=bounds.lb-1e-9) and np.all(z<=bounds.ub+1e-9)
            foc = p["theta"]/z[2]-p["chi"]/z[0]-p["alpha"]*p["kappa"]/(z[1]-p["kappa"]*z[2])
            max_foc = max(max_foc, abs(float(foc)))
            minimum_slack["saving"] = min(minimum_slack["saving"], z[3])
            wedge = p["alpha"]*z[0]/(z[1]-p["kappa"]*z[2])-u
            if owner:
                minimum_slack["retention"] = min(minimum_slack["retention"], z[1]-z[5])
                minimum_slack["estate"] = min(minimum_slack["estate"], z[6]-row["prices"][2]*z[5])
                minimum_slack["owner_cap"] = min(minimum_slack["owner_cap"], p["owner_cap"]-z[1])
                minimum_slack["owner_housing_wedge"] = min(minimum_slack["owner_housing_wedge"], wedge)
            else:
                minimum_slack["renter_young_wedge"] = min(minimum_slack["renter_young_wedge"], wedge)
                minimum_slack["renter_old_wedge"] = min(minimum_slack["renter_old_wedge"], p["gamma"]*z[4]/z[5]-hh["u1"])
            if row["t"] in optimize_dates:
                # Scaling only changes numerical units, which matters for tiny rental caps.
                scaled_bounds=Bounds(bounds.lb/z,bounds.ub/z)
                initial=np.minimum(np.array([.98,.99,1.01,1.02,.98,1.02,.99]),scaled_bounds.ub)
                fit = minimize(lambda zz: -household_utility(zz*z, p), initial, method="SLSQP",
                               bounds=scaled_bounds, constraints=[LinearConstraint(M*z, rhs, rhs),
                                                                 LinearConstraint(G*z, lo, hi)],
                               options={"ftol": 1e-12, "maxiter": 600})
                optimum=fit.x*z
                discrepancy = float(np.max(abs(optimum-z)))
                assert fit.success and discrepancy<3e-5, (tenure, row["t"], fit.message, discrepancy)
                assert abs(household_utility(optimum,p)-household_utility(z,p))<1e-9
                optimizations.append(dict(date=row["t"], tenure=tenure, maximum_choice_error=discrepancy))
        # Actual old owners at date zero retain pre-reform assets and title.
        for tenure in ("owner", "renter"):
            old = row["old_choices"][tenure]
            past = row["past"][tenure]
            assert np.all(old>0)
            resources = past["assets"]+row["transfer"]
            if tenure=="owner":
                resources += row["price"]*past["z"][1]
                assert 0<old[1]<past["z"][1]
                assert old[2]>row["prices"][1]*old[1]
            else:
                assert old[1]<=p["rental_cap"]
                minimum_slack["actual_old_renter_wedge"]=min(
                    minimum_slack["actual_old_renter_wedge"],p["gamma"]*old[0]/old[1]-u)
            error = old[0]+u*old[1]+p["q"]*old[2]-resources
            max_budget = max(max_budget, abs(float(error)))
    assert max_budget<1e-10 and max_foc<1e-9
    assert all(v>0 for v in minimum_slack.values()), minimum_slack
    return dict(maximum_original_budget_error=max_budget, maximum_fertility_foc_error=max_foc,
                minimum_branch_margins={k:float(v) for k,v in minimum_slack.items()},
                original_optimizations=optimizations)


def analytical_check():
    z, ym, y0, y1, y2, d, H = sp.symbols("z ym y0 y1 y2 d H", positive=True)
    q, D, w = sp.Rational(4,5), sp.Rational(3,58), sp.Rational(23,20)
    H0 = sp.Rational(1213,928)
    market = y1+D*((w/d-1)/q*ym+y0**2/y1)/(y0/y1-q*y1/y2)-H
    base = {ym:1,y0:1,y1:1,y2:1,d:1,H:H0}
    polynomial = sum(sp.diff(market,v).subs(base)*z**i for i,v in enumerate((ym,y0,y1,y2)))
    wanted = 1140*z**3-3253*z**2+945*z-45
    ratio = sp.simplify(polynomial/wanted)
    assert not ratio.free_symbols and ratio!=0
    intervals = [(sp.Rational(59,1000),sp.Rational(60,1000)),
                 (sp.Rational(261,1000),sp.Rational(262,1000)),
                 (sp.Rational(253,100),sp.Rational(254,100))]
    signs = []
    for left,right in intervals:
        values = [wanted.subs(z,left),wanted.subs(z,right)]
        assert values[0]*values[1]<0
        signs.append([str(left),str(right),str(values[0]),str(values[1])])
    roots = np.sort(np.array([float(v) for v in sp.nroots(wanted)]))
    initial = y1+sp.Rational(1,11)*(-sp.Rational(301,928)+d/y1)/(d/y1-q*d*y1/y2)-H
    initbase={y1:1,y2:1,d:1,H:H0}
    boundary = np.array([float(sp.diff(initial,v).subs(initbase)) for v in (y1,y2,d,H)])
    f1,f2,fd,_ = boundary
    stable = roots[:2]
    matrix=np.array([[1.,1.],[f1*stable[0]+f2*stable[0]**2,f1*stable[1]+f2*stable[1]**2]])
    gradient = 345/1213
    coefficients=np.linalg.solve(matrix,[-gradient,-(f1+f2)*gradient-fd])
    t=np.arange(41)
    young=gradient+coefficients[0]*stable[0]**t+coefficients[1]*stable[1]**t
    increment_floor=-coefficients[0]*(1-stable[0])-coefficients[1]*(1-stable[1])
    assert increment_floor>0 and coefficients[0]>0 and coefficients[1]<0
    # Certify the small positive initial response with rational interval bounds.
    F1,F2,Fd=sp.Rational(33783,10208),-sp.Rational(285,232),sp.Rational(1505,10208)
    J=sp.Rational(345,1213)
    s_lower=intervals[0][0]+intervals[1][0]
    s_upper=intervals[0][1]+intervals[1][1]
    product_lower=(1-intervals[0][1])*(1-intervals[1][1])
    assert F1+F2*s_upper>0
    numerator_lower=-Fd-F2*J*product_lower
    certified_m_lower=numerator_lower/(F1+F2*s_lower)
    assert certified_m_lower>0
    return dict(polynomial=str(wanted), rational_sign_brackets=signs, roots=roots.tolist(),
                initial_boundary_derivatives=boundary.tolist(), boundary_determinant=float(np.linalg.det(matrix)),
                stationary_young_derivative=gradient, coefficients=coefficients.tolist(),
                uniform_first_order_increment_bound=float(increment_floor),
                rational_initial_response_lower_bound=str(certified_m_lower), young_derivative=young.tolist())


def analytical_figure(analytical):
    """Plot equation T6 itself, independently of every finite-path calculation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t=np.arange(9)
    young=np.asarray(analytical["young_derivative"])
    fertility=np.diff(young[:10])/2
    household_population=young[:9]+np.r_[0.,young[:8]]
    terminal=2*analytical["stationary_young_derivative"]
    plt.rcParams.update({"font.family":"serif","font.size":11,
                         "axes.spines.top":False,"axes.spines.right":False,
                         "pdf.fonttype":42,"ps.fonttype":42})
    fig,axes=plt.subplots(1,2,figsize=(9.6,2.85),layout="constrained")
    for ax in axes:
        ax.set_xlabel("Generations after reform")
        ax.set_xlim(-.15,8.2)
        ax.set_xticks([0,2,4,6,8])
        ax.axhline(0,color=".55",lw=.8)
        ax.grid(axis="y",alpha=.17)
    axes[0].plot(t,fertility,color="#255f85",marker="o",ms=4,lw=1.8)
    axes[0].set_title("Fertility")
    axes[0].set_ylabel(r"First-order response, $\partial n_t/\partial d$")
    axes[0].set_ylim(-.005,.105)
    axes[0].set_yticks([0,.04,.08])
    axes[0].annotate("Returns to replacement",xy=(6,fertility[6]),xytext=(3,.05),
                     fontsize=9,arrowprops={"arrowstyle":"-","color":".45"})
    axes[1].plot(t,household_population,color="#996238",marker="o",ms=4,lw=1.8)
    axes[1].axhline(terminal,color="#996238",ls="--",lw=1)
    axes[1].set_title("Young and old households")
    axes[1].set_ylabel(r"First-order response, $\partial N_{\rm hh,t}/\partial d$")
    axes[1].set_ylim(-.03,.68)
    axes[1].set_yticks([0,.2,.4,.6])
    axes[1].text(4,.595,"Larger stationary population",fontsize=9,ha="center")
    paths=[OUT/"local_transition_analytical_figure.pdf",OUT/"local_transition_analytical_figure.png"]
    fig.savefig(paths[0],bbox_inches="tight")
    fig.savefig(paths[1],dpi=180,bbox_inches="tight")
    plt.close(fig)
    return [str(path) for path in paths]


def serializable(value):
    if isinstance(value, dict):
        return {k:serializable(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def stationary_population_checks():
    wchi, ell, d, rho, g, alpha, theta, nu, kappa, chi = sp.symbols(
        "wchi ell d rho g alpha theta nu kappa chi", positive=True)
    X = (wchi-ell*d)/rho
    h = kappa/nu*(nu*(alpha+theta)*X-chi)/(nu*theta*X-chi)
    total_h = h*(1+g*X/d)
    formula = (g*wchi/(d*(d+g*X))-alpha*nu*chi*ell/
               ((nu*theta*X-chi)*(nu*(alpha+theta)*X-chi)))/rho
    identity = sp.factor(-sp.diff(sp.log(total_h),d)-formula)
    assert identity==0
    results=[]
    for name,p in [("zero_child_goods",parameters()),
                   ("positive_child_goods",parameters(chi=.015)),
                   ("small_old_housing_weight",parameters(chi=.015,gamma=.001,rental_cap=.001))]:
        steady=steady_state(p)
        horizon_zero=transition(steady,steady,0)
        checks=check_path(horizon_zero,p,(0,))
        ell0=1-p["q"]
        rho0=1+p["beta"]*(1+p["gamma"]+p["omega"])
        wchi0=p["y"]+p["b"]-p["chi"]/p["nu"]
        d0=p["b"]/(1-p["phi"])
        g0=p["beta"]*p["gamma"]/(p["q"]*ell0)
        subs={wchi:wchi0,ell:ell0,d:d0,rho:rho0,g:g0,alpha:p["alpha"],
              theta:p["theta"],nu:p["nu"],kappa:p["kappa"],chi:p["chi"]}
        predicted=float(formula.subs(subs))
        for expression,actual in [(X,steady["households"]["owner"]["z"][0]),
                                  (h,steady["households"]["owner"]["z"][1]),
                                  (g*X*h/d,steady["households"]["owner"]["z"][5])]:
            assert abs(float(expression.subs(subs))-actual)<1e-9
        step=1e-5
        minus=steady_state(dict(p,phi=1-p["b"]/(d0-step)))
        plus=steady_state(dict(p,phi=1-p["b"]/(d0+step)))
        actual=(np.log(plus["cohort"])-np.log(minus["cohort"]))/2/step
        assert abs(predicted-actual)<1e-8,(name,predicted,actual)
        assert (predicted<0)==(name=="small_old_housing_weight")
        results.append(dict(name=name,parameters=p,log_population_derivative=predicted,
                            independent_stationary_derivative=float(actual),
                            price=steady["price"],young_housing=float(steady["households"]["owner"]["z"][1]),
                            old_housing=float(steady["households"]["owner"]["z"][5]),checks=checks))
    return dict(symbolic_identity=str(identity),cases=results)


def heterogeneity_checks():
    """Check a compact entrant rectangle without assuming away tenure sorting."""
    y,b,P0,P1,P2,phi=sp.symbols("y b P0 P1 P2 phi",positive=True)
    q,rho=sp.Rational(4,5),sp.Rational(58,25)
    d=b/(1-phi)
    h=d/P0
    x=(y+b-d+q*d*P1/P0)/rho
    c2=sp.Rational(2,5)*x/q
    h2=sp.Rational(3,10)*c2/(P1-q*P2)
    saving=y+b-x-(1-phi)*P0*h
    objects=dict(housing=h,fertility=h/2,adult_consumption=x,saving=saving,
                 financial_assets=saving/q-phi*P0*h/q,old_housing=h2,
                 old_consumption=c2,estate=2*c2/q)
    for expression in objects.values():
        assert sp.hessian(expression,(y,b))==sp.zeros(2)
    corners=[dict(y=income,b=wealth) for income in (.94,.96) for wealth in (.195,.205)]
    settings=[("limiting_stationary",parameters(),[1.,1.,1.],[0.,0.]),
              ("limiting_dated",parameters(phi=1-.2/1.005),[1.003,1.004,1.0047],[0.,0.]),
              ("mixed_dated",parameters(phi=1-.2/1.001,chi=.005,tau=.005,taste_weight=.02),
               [1.0,1.0005,1.0007],[.0026,.0026])]
    results=[]
    for name,p,prices,transfers in settings:
        type_rows=[]
        owner_vectors=[]
        actual_old_vectors=[]
        average_reference=young_choices(prices,transfers,p)
        for corner in corners:
            pp=dict(p,**corner)
            initial_p=dict(pp,phi=.8)
            # Predetermined assets and title use pre-reform conditional choices.
            past=young_choices([1.,1.,1.],[transfers[0]]*2,initial_p)
            hh=young_choices(prices,transfers,pp)
            old={tenure:old_choices(past[tenure]["assets"],past[tenure]["z"][1],
                                    prices[0],prices[1],transfers[0],pp,tenure=="owner")
                 for tenure in ("owner","renter")}
            row=dict(t=0,hh=hh,prices=np.asarray(prices),transfers=np.asarray(transfers),
                     past=past,old_choices=old,price=prices[0],transfer=transfers[0])
            checks=check_path(dict(rows=[row]),pp,(0,))
            owner_vectors.append(np.r_[hh["owner"]["z"],hh["owner"]["assets"]])
            actual_old_vectors.append(old["owner"])
            type_rows.append(dict(parameters=pp,owner_choices=hh["owner"]["z"],
                                  owner_share=hh["pi"],checks=checks))
        aggregation_error=None
        if p["chi"]==p["tau"]==p["taste_weight"]==0:
            reference=np.r_[average_reference["owner"]["z"],average_reference["owner"]["assets"]]
            averaging=np.mean(np.asarray(owner_vectors),axis=0)
            refpast=young_choices([1.,1.,1.],[0.,0.],dict(p,phi=.8))
            refold=old_choices(refpast["owner"]["assets"],refpast["owner"]["z"][1],
                               prices[0],prices[1],transfers[0],p,True)
            aggregation_error=float(max(np.max(abs(reference-averaging)),
                                        np.max(abs(refold-np.mean(actual_old_vectors,axis=0)))))
            assert aggregation_error<1e-12
        results.append(dict(name=name,prices=prices,transfers=transfers,
                            exact_aggregation_error=aggregation_error,types=type_rows))
    return dict(scope="Uniform conditional branches and exact limiting aggregation. Dated prices in these checks are diagnostic inputs, not equilibrium paths.",
                mean_income=.95,mean_wealth=.2,support=dict(income=[.94,.96],wealth=[.195,.205]),
                symbolic_hessians="All eight owner quantity/asset maps have identically zero Hessian in income and wealth.",
                cases=results)


def general_transition_checks():
    """Symbolic root/boundary identities, a bounded coefficient grid, and one economy."""
    q,D,C,L,z,ym,y0,y1,y2,r,delta=sp.symbols(
        "q D C L z ym y0 y1 y2 r delta",positive=True)
    ell=1-q
    polynomial=C*q*z**3-(ell-D+C*(1+q))*z**2+(C-2*D)*z+D-ell*C
    market=y1+D*((r-1)/q*ym+y0**2/y1)/(y0/y1-q*y1/y2)
    base={ym:1,y0:1,y1:1,y2:1,r:ell+q*ell*C/D}
    derived=sum(sp.diff(market,v).subs(base)*z**i for i,v in enumerate((ym,y0,y1,y2)))
    assert sp.factor(-ell*derived-polynomial)==0
    R=(ell*C-D)/(q*C)
    assert sp.factor(polynomial.subs(z,1)+ell*(1+C))==0
    assert sp.factor(sp.diff(polynomial,z).subs(z,1)+ell*(2+C))==0
    assert sp.factor(polynomial.subs(z,R)+ell*(ell*C-D)*(q*C*C+ell*C-D)/(C*C*q*q))==0
    assert sp.factor(polynomial.subs(z,-1)-(4*D-ell-C*(3+q)))==0
    v=sp.symbols("v")
    transformed=(1-v)**3*polynomial.subs(z,(1+v)/(1-v))
    expected_transform=(ell+(3+q)*C-4*D)*v**3+(ell+4*D+(7*q-3)*C)*v**2+ell*(C-1)*v-ell*(C+1)
    assert sp.factor(transformed-expected_transform)==0
    initial=y1+L*((C*ell/L-1)/(1+delta)+1/y1)/(1/y1-q*y1/y2)
    initial_base={y1:1,y2:1,delta:0}
    expected=[1+(C*(1+q)-L)/ell,-C*q/ell,L/ell-C]
    for v,target in zip((y1,y2,delta),expected):
        assert sp.factor(sp.diff(initial,v).subs(initial_base)-target)==0
    grid=[]
    for qq in (.05,.2,.4,.5,.8,.95):
        ll=1-qq
        for lf in (.2,.5,.8):
            LL=lf*ll
            for df in (.2,.5,.8):
                DD=df*LL
                for cf in (.05,.1,.25,.5,.9):
                    CC=cf
                    rr=ll+qq*ll*CC/DD
                    coef=[CC*qq,-(ll-DD+CC*(1+qq)),CC-2*DD,DD-ll*CC]
                    roots=np.roots(coef)
                    stable=roots[np.abs(roots)<1]
                    unstable=roots[np.abs(roots)>1]
                    stability_gap=ll+(3+qq)*CC-4*DD
                    assert (len(stable)==2 and len(unstable)==1)==(stability_gap>0)
                    if stability_gap<0:
                        grid.append(dict(q=qq,D=DD,L=LL,C=CC,r=rr,
                                         stability_gap=stability_gap,stable_root_count=len(stable)))
                        continue
                    Lambda=float(unstable[0].real)
                    assert abs(unstable[0].imag)<1e-12 and Lambda>1
                    S=float(np.sum(stable).real)
                    f1,f2,fd=1+(CC*(1+qq)-LL)/ll,-CC*qq/ll,LL/ll-CC
                    denominator=f1+f2*S
                    assert denominator>1-LL/ll>0
                    J=DD*rr/(qq*ll*(1+CC))
                    m=(CC-LL/ll+J*(1+CC)/(Lambda-1))/denominator
                    if CC>=LL/ll:
                        assert m>0
                        exact_test=True
                    else:
                        Z=1+DD*rr/(qq*(LL-ll*CC))
                        test=float(np.polyval(coef,Z))
                        assert (test>0)==(m>0)
                        exact_test=test>0
                    sufficient_gap=ll*CC*(ll+(1+qq)*CC)-LL*(ll-DD+CC)
                    if qq>=.5 and rr>1 and sufficient_gap>0:
                        assert m>0
                    grid.append(dict(q=qq,D=DD,L=LL,C=CC,r=rr,
                                     stability_gap=stability_gap,stable_root_count=len(stable),
                                     maximum_stable_modulus=float(max(abs(stable))),
                                     boundary_factor=denominator,initial_increment=m,
                                     primitive_exact_test=exact_test,sufficient_gap=sufficient_gap))
    pp=parameters(q=.4,y=1.,alpha=1.2,theta=.4,beta=1/12,omega=1.7,
                  rental_cap=.025,Hbar=1.05)
    ss=steady_state(pp)
    assert abs(ss["price"]-1)<1e-10 and abs(ss["cohort"]-1)<1e-10
    assert np.max(abs(ss["households"]["owner"]["z"]-
                      np.array([.48,1.,.5,.52,.1,.05,.425])))<1e-10
    counter_roots=np.roots([.05*.4,-(.6-.02+.05*1.4),.05-2*.02,.02-.6*.05])
    counter_stable=counter_roots[np.abs(counter_roots)<1]
    S=float(np.sum(counter_stable).real)
    f1,f2,fd=1+(.05*1.4-.1)/.6,-.05*.4/.6,.1/.6-.05
    J=2/21
    m=float(((-fd-f2*J*np.prod(1-counter_stable))/(f1+f2*S)).real)
    assert m<0
    cases=[]
    for sign in (-1,1):
        changed=dict(pp,phi=1-.2/(1+sign*1e-4))
        final=steady_state(changed)
        path=transition(ss,final,24)
        checks=check_path(path,changed,(0,2) if sign==1 else ())
        cases.append(dict(sign=sign,parameters=changed,first_fertility=float(path["rows"][0]["hh"]["fertility"]),
                          final_cohort=final["cohort"],maximum_residual=path["maximum_residual"],checks=checks))
    actual=(cases[1]["first_fertility"]-cases[0]["first_fertility"])/2e-4
    assert abs(actual-m/2)<1e-8
    assert cases[1]["first_fertility"]<.5 and cases[1]["final_cohort"]>1
    return dict(scope="Symbolic general recurrence and boundary identities, 270 declared coefficient checks including r below one, and an original-equation counterexample. The grid checks arithmetic, not the proof or an economically calibrated parameter region.",
                symbolic_identities="All six polynomial/transform and three initial-boundary identities are exact.",
                coefficient_cases=grid,
                counterexample=dict(parameters=pp,roots=serializable(counter_roots.real),
                                    root_imaginary_parts=serializable(counter_roots.imag),
                                    initial_population_increment_derivative=m,stationary_population_derivative=J,
                                    predicted_initial_fertility_derivative=m/2,
                                    original_equation_initial_fertility_derivative=actual,cases=cases))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke",action="store_true")
    parser.add_argument("--stationary-only",action="store_true")
    parser.add_argument("--figure",action="store_true")
    parser.add_argument("--heterogeneity-only",action="store_true")
    parser.add_argument("--general-only",action="store_true")
    args=parser.parse_args()
    start=time.monotonic()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.general_only:
        report=general_transition_checks()
        report["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        target=OUT/"local_transition_general_checks.json"
        target.write_text(json.dumps(serializable(report),indent=2)+"\n")
        print(json.dumps(dict(output=str(target),coefficient_cases=len(report["coefficient_cases"]),
                             initial_fertility_derivative=report["counterexample"]["predicted_initial_fertility_derivative"],
                             stationary_population_derivative=report["counterexample"]["stationary_population_derivative"]),indent=2))
        return
    if args.figure:
        print(json.dumps(dict(figures=analytical_figure(analytical_check())),indent=2))
        return
    if args.heterogeneity_only:
        report=heterogeneity_checks()
        report["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        target=OUT/"local_transition_heterogeneity_checks.json"
        target.write_text(json.dumps(serializable(report),indent=2)+"\n")
        print(json.dumps(dict(output=str(target),cases=[dict(name=c["name"],aggregation_error=c["exact_aggregation_error"]) for c in report["cases"]]),indent=2))
        return
    if args.stationary_only:
        report=stationary_population_checks()
        report["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        target=OUT/"local_transition_stationary_checks.json"
        target.write_text(json.dumps(serializable(report),indent=2)+"\n")
        print(json.dumps(dict(output=str(target),cases=[dict(name=c["name"],derivative=c["log_population_derivative"]) for c in report["cases"]]),indent=2))
        return
    analytical=analytical_check()
    base=parameters()
    cases=[("zero",base,base,24),
           ("credit_plus_tiny",base,parameters(phi=1-.2/1.0001),24)]
    if not args.smoke:
        mixed=parameters(chi=.005,tau=.005,taste_weight=.02)
        positive_cost=parameters(chi=.015)
        cases += [("credit_minus_tiny",base,parameters(phi=1-.2/.9999),24),
                  ("credit_finite_24",base,parameters(phi=1-.2/1.005),24),
                  ("credit_finite_40",base,parameters(phi=1-.2/1.005),40),
                  ("stock_finite",base,parameters(Hbar=base["Hbar"]*1.005),24),
                  ("positive_cost_credit",positive_cost,dict(positive_cost,phi=1-.2/1.001),24),
                  ("mixed_credit",mixed,dict(mixed,phi=1-.2/1.001),24)]
    results=[]
    for name,p0,p1,horizon in cases:
        initial,final=steady_state(p0),steady_state(p1)
        path=transition(initial,final,horizon)
        dates=(0,2,10) if name=="credit_finite_24" else (0,2) if name=="mixed_credit" else ()
        check=check_path(path,p1,dates)
        rows=[dict(date=r["t"],price=float(r["price"]),young=float(r["young"]),old=float(r["old"]),
                   fertility=float(r["hh"]["fertility"]),owner_share=float(r["hh"]["pi"])) for r in path["rows"]]
        result=dict(name=name,parameters_before=p0,parameters_after=p1,horizon=horizon,
                    initial_steady={k:float(initial[k]) for k in ("price","transfer","cohort")},
                    final_steady={k:float(final[k]) for k in ("price","transfer","cohort")},
                    maximum_equilibrium_residual=path["maximum_residual"],checks=check,rows=rows)
        results.append(result)
        print(json.dumps(dict(case=name,residual=path["maximum_residual"],
                              first_fertility=rows[0]["fertility"],final_cohort=final["cohort"])),flush=True)
    report=dict(scope="Analytical cubic and original-equation checks. Finite terminally closed paths are supporting evidence, not the convergence proof.",
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                original_specification_sha256=hashlib.sha256((ROOT/"latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex").read_bytes()).hexdigest(),
                analytical=analytical,cases=results,elapsed_seconds=time.monotonic()-start)
    if not args.smoke:
        byname={r["name"]:r for r in results}
        plus=np.array([r["young"] for r in byname["credit_plus_tiny"]["rows"]])
        minus=np.array([r["young"] for r in byname["credit_minus_tiny"]["rows"]])
        derivative=(plus-minus)/.0002
        error=float(np.max(abs(derivative-np.array(analytical["young_derivative"][:25]))))
        assert error<2e-6,error
        short=np.array([[r["price"],r["young"]] for r in byname["credit_finite_24"]["rows"]])
        long=np.array([[r["price"],r["young"]] for r in byname["credit_finite_40"]["rows"][:25]])
        horizon_error=float(np.max(abs(short-long)))
        assert horizon_error<1e-9,horizon_error
        report.update(central_derivative_error=error,horizon_comparison_error=horizon_error)
    target=OUT/("local_transition_smoke.json" if args.smoke else "local_transition_checks.json")
    target.write_text(json.dumps(serializable(report),indent=2)+"\n")
    print(json.dumps(dict(output=str(target),elapsed_seconds=report["elapsed_seconds"]),indent=2))


if __name__=="__main__":
    main()
