#!/usr/bin/env python3
"""Verify the supporting transition extensions without changing the main note.

These are analytical theory checks, not calibration or a finite-horizon
simulation. The receipt separates exact identities from floating-point checks
against the original household equations.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"
from verify_simplified_olg_local_transition import (
    complex_jacobian,
    parameters,
    young_choices,
)


def symbolic_checks():
    q, beta, gamma, omega, alpha, theta, kappa, w, d, h, a = sp.symbols(
        "q beta gamma omega alpha theta kappa w d h a", positive=True
    )
    ell = 1 - q
    rho_o = 1 + beta * (1 + gamma + omega)
    rho_r = 1 + beta * (1 + omega)
    price = d / h
    x_o = (w - ell * d) / rho_o
    x_r = (w - (1 + q) * ell * price * a) / rho_r

    def value(x, house, old_house):
        n = theta * house / (kappa * (alpha + theta))
        return (
            sp.log(x) + alpha * sp.log(house - kappa * n)
            + theta * sp.log(n)
            + beta * (sp.log(beta * x / q) + gamma * sp.log(old_house)
                      + omega * sp.log(beta * omega * x / q**2))
        )

    delta = value(x_o, h, beta * gamma * x_o / (q * ell * price)) - value(x_r, a, a)
    v = (1 + q) * ell * price * a / x_r
    m = ell * d / x_o
    A = alpha + theta + beta * gamma - v
    F = m + beta * gamma - v
    identities = {
        "owner_value_minus_renter_value_h": sp.diff(delta, h) - A / h,
        "owner_value_minus_renter_value_d": sp.diff(delta, d) + F / d,
        "owner_value_minus_renter_value_theta": sp.expand_log(
            sp.diff(delta, theta) - sp.log(h / a), force=True
        ),
    }
    C = beta * gamma * x_o / (q * ell * d)
    identities["old_housing_ratio_derivative"] = d * sp.diff(C, d) + C + beta * gamma / (q * rho_o)
    x, c, k, C0, L = sp.symbols("x c k C L", positive=True)
    inherited_assets_over_d = ell * C0 / L - 1
    initial = x + L * (inherited_assets_over_d * x + k) / (1 - q * x**2 / c) - k * (1 + C0)
    base = {x: 1, c: 1, k: 1}
    identities["actual_old_initial_x"] = sp.diff(initial, x).subs(base) - (1 + (C0 * (1 + q) - L) / ell)
    identities["actual_old_initial_c"] = sp.diff(initial, c).subs(base) + C0 * q / ell
    identities["actual_old_preference_shift"] = sp.diff(initial, k).subs(base) - (L / ell - 1 - C0)
    # Reconstruct the finite recurrence and surprises from the dated budgets.
    K, stock, lag, prev, nxt, yy, assets, titles, d0 = sp.symbols(
        "K stock lag prev nxt yy assets titles d0", positive=True
    )
    Ppast, Pnow, Pfuture = K*d*lag/prev, K*d*prev/x, K*d*x/nxt
    xp = (w-d+q*d*Pnow/Ppast)/rho_o
    clearing = x/K + lag*beta*gamma*xp/(q*(Pnow-q*Pfuture)) - stock
    residual = (K*stock-x)*(prev-q*x*x/nxt) - beta*gamma/rho_o*((w/d-1)/q*lag*x+prev**2)
    identities["finite_interior_housing_clearing"] = clearing*K*(prev-q*x*x/nxt)+residual
    Pnow, Pfuture = K*d*yy/x, K*d*x/nxt
    clearing0 = x/K+L*(assets+Pnow*titles)/(Pnow-q*Pfuture)-stock
    residual0 = (K*stock-x)*(yy-q*x*x/nxt)-L*(assets*x/d+K*yy*titles)
    identities["finite_actual_old_housing_clearing"] = clearing0*K*(yy-q*x*x/nxt)+residual0
    prior_x = (w-d0+q*d0*yy**2/(lag*nxt))/rho_o
    inherited = lag*(w-d0-prior_x)/q
    identities["baseline_generated_assets"] = inherited-((rho_o-1)*(w-d0)*lag/(q*rho_o)-d0*yy**2/(rho_o*nxt))
    identities["finite_interior_credit_forcing"] = sp.diff(residual,d)-beta*gamma/rho_o*w*lag*x/(q*d**2)
    identities["finite_initial_credit_forcing"] = sp.diff(residual0,d)-L*assets*x/d**2
    for name, expression in identities.items():
        assert sp.simplify(expression) == 0, name
    return {name: "exact zero after differentiation of the original budgets/utilities" for name in identities}


def original_margins(p, price=1.0):
    hh = young_choices([price] * 3, [0.0, 0.0], p)
    q, u = p["q"], (1 - p["q"]) * price
    residuals, margins = {}, {}
    for tenure in ("owner", "renter"):
        own = tenure == "owner"
        x, h, n, saving, c2, h2, estate = hh[tenure]["z"]
        cost = (1 - p["phi"]) * price * h if own else u * h
        residuals[tenure + "_young_budget"] = x + p["chi"] * n + saving + cost - p["y"] - p["b"]
        residuals[tenure + "_old_budget"] = c2 + q * estate + u * h2 - hh[tenure]["assets"] - (price * h if own else 0)
        residuals[tenure + "_fertility"] = p["theta"] / n - p["chi"] / x - p["alpha"] * p["kappa"] / (h - p["kappa"] * n)
        margins[tenure + "_positive_bundle"] = min(x, h - p["kappa"] * n, n, saving, c2, h2, estate)
        margins[tenure + "_young_housing_gap"] = p["alpha"] * x / (h - p["kappa"] * n) - u
        if own:
            residuals["down_payment"] = (1 - p["phi"]) * price * h - p["b"]
            margins.update(owner_physical=p["owner_cap"] - h,
                           owner_retention=h - h2, owner_estate=estate - price * h2)
        else:
            margins["old_renter_housing_gap"] = p["gamma"] * c2 / h2 - u
    assert max(abs(z) for z in residuals.values()) < 1e-10, residuals
    assert min(margins.values()) > 0, margins
    return hh, {"maximum_original_equation_residual": max(abs(z) for z in residuals.values()),
                "strict_margins": margins}


def limiting_witnesses():
    cases = {
        "convergent_with_w_less_than_d": parameters(y=.6, alpha=.9, theta=.3, rental_cap=.05, Hbar=277/232),
        "feasible_but_wrong_stable_dimension": parameters(q=.5, y=.35, beta=.5, gamma=1., omega=2.,
                                                         alpha=24., theta=8., rental_cap=.01, Hbar=31/30),
    }
    result = {}
    for name, p in cases.items():
        hh, checks = original_margins(p)
        q, ell = p["q"], 1 - p["q"]
        rho = 1 + p["beta"] * (1 + p["gamma"] + p["omega"])
        D = p["beta"] * p["gamma"] / rho
        r = (p["y"] + p["b"]) / (p["b"] / (1 - p["phi"]))
        C = D * (r - ell) / (q * ell)
        gap = ell + (3 + q) * C - 4 * D
        exact_gap = sp.Rational(847, 1160) if name.startswith("convergent") else sp.Rational(-1, 20)
        rq, rb, rg, ro, ry, rwealth, rphi = [sp.Rational(str(p[key])) for key in ("q", "beta", "gamma", "omega", "y", "b", "phi")]
        rrho = 1 + rb*(1+rg+ro)
        rD = rb*rg/rrho
        rr = (ry+rwealth)/(rwealth/(1-rphi))
        rC = rD*(rr-(1-rq))/(rq*(1-rq))
        assert sp.simplify((1-rq)+(3+rq)*rC-4*rD-exact_gap) == 0
        assert abs(gap - float(exact_gap)) < 1e-13
        roots = np.roots([C*q, -(ell-D+C*(1+q)), C-2*D, D-ell*C])
        assert sum(abs(z) < 1 for z in roots) == (2 if gap > 0 else 1)
        checks.update(parameters=p, r=r, C=C, exact_convergence_gap=str(exact_gap),
                      stable_root_count=int(sum(abs(z) < 1 for z in roots)))
        result[name] = checks
    return result


def mixed_stationary_checks():
    result = []
    for scale in (.25, 1., 4.):
        p = parameters(q=.5, beta=.4, gamma=.3, omega=.4, alpha=.2, theta=7/15,
                       kappa=7/8, y=99/50, Hbar=99/100, rental_cap=.25, sigma=scale)
        base = young_choices([1.]*3, [0., 0.], p)
        p["taste_weight"] = np.exp((base["owner"]["utility"]-base["renter"]["utility"])/scale)
        hh, checks = original_margins(p)
        assert abs(hh["pi"]-.5) < 1e-13 and abs(hh["fertility"]-.5) < 1e-13

        def original_outputs(v):
            h, d, theta = v
            varied = dict(p, phi=1-p["b"]/d, theta=theta)
            H = young_choices([d/h]*3, [0., 0.], varied)
            old_h = H["pi"]*H["owner"]["z"][5] + (1-H["pi"])*p["rental_cap"]
            return np.array([H["fertility"], H["housing"], old_h, H["pi"], d/h])

        J = complex_jacobian(original_outputs, [1., 1., p["theta"]])
        derivatives = {}
        for shock, col, conversion in (("phi", 1, 5.), ("theta", 2, 1.)):
            h_derivative = -J[0, col] / J[0, 0]
            full = (J[:, col] + J[:, 0]*h_derivative)*conversion
            base_output = original_outputs([1., 1., p["theta"]])
            S = base_output[1] + base_output[2]
            assert abs(S-p["Hbar"]) < 1e-12
            N_derivative = -2*p["Hbar"]*(full[1]+full[2])/S**2
            assert abs(full[0]) < 1e-12 and N_derivative > 0 and full[4] > 0
            if shock == "phi":
                assert 0 < h_derivative < 1 and full[3] < 0 and abs(full[1]) < 1e-12 and full[2] < 0
            else:
                assert h_derivative < 0 and full[1] < 0 and full[2] < 0
            derivatives[shock] = dict(N=float(N_derivative), hOwner=float(h_derivative*conversion),
                                     hYoung=float(full[1]), hOld=float(full[2]), pi=float(full[3]), P=float(full[4]))
        checks.update(taste_scale=scale, taste_location=float(-scale*np.log(p["taste_weight"])),
                      owner_share=float(hh["pi"]), derivatives=derivatives)
        result.append(checks)
    return result


def finite_transition_certificate():
    """Exact infinite-sequence bounds; no finite terminal condition."""
    from fractions import Fraction as F

    q, beta, gamma, omega = F(4,5), F(2,5), F(3,10), F(2)
    alpha, nu, kappa, b, y = F(2,5), F(2), F(1,2), F(1,5), F(19,20)
    ell = 1-q
    A = 1+gamma+omega
    rho = 1+beta*A
    D, L, w, Hbar = beta*gamma/rho, gamma/A, y+b, F(1213,928)
    k = F(999,1000)                         # K after the preference shock; Kpre=1.
    theta = alpha*k/(nu/kappa-k)
    S = k*Hbar
    dlo, dhi = F(1), F(101,100)
    elo, ehi = (w/dhi-1)/q, (w/dlo-1)/q
    mB, MB = k-D*ehi/(ell+D)*(1-k), F(1)
    mP, MP = F(997,1000), F(251,250)
    a_pre, h_pre = -F(301,928), F(1)
    report = {'theta1':str(theta), 'phi1_max':str(1-b/dhi),
              'baseline_box':list(map(float,(mB,MB))),
              'policy_box':list(map(float,(mP,MP)))}


    def full_box(name, dmin, dmax, m, M):
        """Uniform in d, all dates, and every neighbor in the population box."""
        emin, emax = (w/dmax-1)/q, (w/dmin-1)/q
        assert emin>0 and 0<m<M<S
        rmin, rmax = (m/M)**2, (M/m)**2
        # E_b>0 gives exact vertex tests for E(m)>=0 and E(M)<=0.
        slacks = {
          'user_cost': m-q*M*M/m,
          'E_b': S-M-2*D*M,
          'self_lower': ell*S-(ell+D)*m-D*emax*M,
          'self_upper': (ell+D)*M+D*emin*m-ell*S,
        }
        assert all(v>=0 for v in slacks.values())
        assert slacks['user_cost']>0 and slacks['E_b']>0
        # Bounds for G=-E_x and the three absolute row coefficients.
        Gmin = m-q*M*M/m+2*q*m*(S-M)/M+D*emin*m
        Gmax = M-q*m*m/M+2*q*M*(S-m)/m+D*emax*M
        ac = (D*emin*m/Gmax, D*emax*M/Gmin)
        bc = ((S-M-2*D*M)/Gmax, (S-m-2*D*m)/Gmin)
        cc = (q*(S-M)*m*m/(M*M*Gmax), q*(S-m)*M*M/(m*m*Gmin))
        row = ac[1]+bc[1]+cc[1]
        assert Gmin>0 and 0<row<1
        # Original young and generated-old owner conditions, checked over d.
        for d in (dmin,dmax):
            xmin=(w-d+q*d*rmin)/rho
            xmax=(w-d+q*d*rmax)/rho
            vals={
              'adult_consumption':xmin,
              'saving':y-xmax,                       # a'=y-x, since closing cost=b.
              'purchase':(alpha+theta)*xmin-d*(1-q*rmin),
              'physical_owner_cap':F(2)-M/(k*m),
              'old_retention':q*d*rmin*(1-q*rmax)-beta*gamma*xmax,
              'old_estate':omega-q*(omega+gamma)*rmax,
            }
            assert all(v>0 for v in vals.values()), (name,d,vals)
            # Each d-dependent inequality above is affine in d, so endpoints suffice.
            for key,val in vals.items():
                slacks[key]=min(slacks.get(key,val),val)
        # Conditional renter caps also stay strict; not needed for the all-owner map.
        rentcap=F(1,4)
        rhoR=1+beta*(1+omega)
        Pmax=k*dmax*M/m
        U=Pmax*(1-q*rmin)
        cash=w-rentcap*(1+q)*U
        renter_margins=(cash,(alpha+theta)*cash-rhoR*rentcap*U,
                        beta*gamma*cash-rhoR*q*rentcap*U)
        assert all(v>0 for v in renter_margins)
        report[name]={'row_bound':float(row),
                     'slacks':{key:float(val) for key,val in slacks.items()},
                     'renter_cap_margins':list(map(float,renter_margins))}
        return Gmin,Gmax,ac,bc,cc

    base_coeffs=full_box('baseline',dlo,dlo,mB,MB)
    pol_coeffs=full_box('policy',dlo,dhi,mP,MP)


    def E0(x,c,Y,assets,d,titles):
        # assets and titles are TOTAL original old claims, not new-policy choices.
        return (S-x)*(Y-q*x*x/c)-L*((assets/d)*x+k*Y*titles)

    # The first preference shock retains the actual Section2 stationary old.
    base_lower=E0(mB,mB,F(1),a_pre,dlo,h_pre)
    base_upper=E0(MB,MB,F(1),a_pre,dlo,h_pre)
    assert base_lower>0 and base_upper<0
    assert base_upper==(k-1)*(ell*(1+D*(w-ell)/(q*ell))-L)

    # At ANY later baseline date, write U=Y_{t-1}, Y=Y_t, V=Y_{t+1}^{baseline}.
    # Actual claims: Atot=cA U - Y^2/(rho V); Htot=Y/k.
    # V is the OLD baseline forecast, even when policy unexpectedly changes prices.
    cA=(rho-1)*(w-dlo)/(q*rho)
    Atot_min=cA*mB-dlo*MB*MB/(rho*mB)
    Atot_max=cA*MB-dlo*mB*mB/(rho*MB)
    assert Atot_min<Atot_max<0
    # E0 increases in Y once Htot=Y/k is substituted, and decreases in Atot.
    assert S-MP-2*L*MB>0
    policy_lower=E0(mP,mP,mB,Atot_max,dhi,mB/k)
    policy_upper=E0(MP,MP,MB,Atot_min,dlo,MB/k)
    assert policy_lower>0 and policy_upper<0
    report['initial_boundary']={
      'baseline_lower':float(base_lower),'baseline_upper':float(base_upper),
      'policy_lower':float(policy_lower),'policy_upper':float(policy_upper),
      'actual_total_assets_range':list(map(float,(Atot_min,Atot_max)))}


    def boundary_coeffs(name,m,M,Ymin,Ymax,amin,amax,dmin,dmax):
        # Here amin,amax bound TOTAL assets for the first-update derivative.
        assert amin<=amax<0
        Gmin=Ymin-q*M*M/m+2*q*m*(S-M)/M+L*amin/dmin
        Gmax=Ymax-q*m*m/M+2*q*M*(S-m)/m+L*amax/dmax
        cc=(q*(S-M)*m*m/(M*M*Gmax), q*(S-m)*M*M/(m*m*Gmin))
        ff=(L*amin*M/(dmin*dmin*Gmin), L*amax*m/(dmax*dmax*Gmax))
        assert Gmin>0 and 0<cc[0]<=cc[1]<1 and ff[0]<=ff[1]<0
        report[name+'_boundary_row']=float(cc[1])
        return cc,ff

    boundary_coeffs('baseline',mB,MB,F(1),F(1),a_pre,a_pre,dlo,dlo)
    c0,f0=boundary_coeffs('policy',mP,MP,mB,MB,Atot_min,Atot_max,dlo,dhi)

    # Original initial-old choices: c2=R/A, h2=L R/u, e=omega R/(q A).
    # Verify R>0, h2<H, e>Pnext*h2 for every actual inherited household.
    rbmin,rbmax=(mB/MB)**2,(MB/mB)**2
    individual_amin=(w-dlo-(w-dlo+q*dlo*rbmax)/rho)/q
    individual_amax=(w-dlo-(w-dlo+q*dlo*rbmin)/rho)/q
    individual_hmin,individual_hmax=mB/(k*MB),MB/(k*mB)
    for name,m,M,Ymin,Ymax,dmin,dmax,amin,amax,hmin,hmax in [
      ('baseline',mB,MB,F(1),F(1),dlo,dlo,a_pre,a_pre,h_pre,h_pre),
      ('policy',mP,MP,mB,MB,dlo,dhi,individual_amin,individual_amax,
       individual_hmin,individual_hmax)]:
        # Independent bounds intentionally allow incompatible corners: conservative.
        Pmin=k*dmin*Ymin/M; Pmax=k*dmax*Ymax/m
        Pnextmax=k*dmax*M/m
        umin=Pmin-q*Pnextmax
        resource_min=amin+Pmin*hmin
        resource_max=amax+Pmax*hmax
        retention=hmin-L*resource_max/umin
        estate=omega*umin-q*gamma*Pnextmax
        assert umin>0 and resource_min>0 and retention>0 and estate>0
        report[name+'_actual_old']={
          'resource_min':float(resource_min), 'retention_slack':float(retention),
          'estate_slack':float(estate)}

    # Differentiate at FIXED actual inherited Y,A,H. z0=0 and, for i>=2,
    # zi=-ai*z(i-2)+bi*z(i-1)+ci*z(i+1)+fi; z1=c0*z2+f0.
    # Uniform coefficient intervals come directly from E and E0, for all dates/d.
    Gmin,Gmax,ac,bc,cc=pol_coeffs
    fc=(D*w*mP*mP/(q*dhi*dhi*Gmax),D*w*MP*MP/(q*dlo*dlo*Gmin))
    row=max(ac[1]+bc[1]+cc[1],c0[1])
    B=max(fc[1],-f0[0])/(1-row)
    assert row<1 and B>0
    # The inverse (I-DT)^-1 exists on bounded sequences, so every zi is in [-B,B].
    # Repeated interval substitution preserves enclosure of the INFINITE solution.
    # The unretained tail remains [-B,B]; no terminal steady-state value is imposed.
    DYADIC=2**80

    def down(x):return F((x.numerator*DYADIC)//x.denominator,DYADIC)
    def up(x):return -down(-x)
    def add(a,b):return down(a[0]+b[0]),up(a[1]+b[1])
    def mul(a,b):
        products=[x*y for x in a for y in b]
        return down(min(products)),up(max(products))
    def neg(a):return -a[1],-a[0]

    N=100
    z=[(F(0),F(0))]+[(-B,B)]*(N+1)
    for _ in range(100):
        new=[(F(0),F(0)),add(f0,mul(c0,z[2]))]
        for i in range(2,N+1):
            new.append(add(add(fc,neg(mul(ac,z[i-2]))),
                           add(mul(bc,z[i-1]),mul(cc,z[i+1]))))
        new.append((-B,B))
        z=new
    assert z[1][0]>F(17,1000) and z[1][1]<F(51,1000)
    report['policy_impact_derivative']={
      'coefficient_row_bound':float(row), 'derivative_norm_bound':float(B),
      'dY1_dd_interval':list(map(float,z[1])),
      'claimed_exact_enclosure':['17/1000','51/1000'],
      'iterations':100,'retained_dates':100,'outward_dyadic_bits':80}

    C0=D*(w/dlo-ell)/(q*ell)
    C1=D*(w/dhi-ell)/(q*ell)
    assert S/(1+C0)==k and S/(1+C1)>k
    report['stationary_young_population']={
      'pre':1.0,'baseline':float(k),'largest_certified_credit':float(S/(1+C1))}
    report["policy_impact_derivative"]["exact_dY1_dd_interval"] = list(map(str, z[1]))
    report["exact_population_boxes"] = {"baseline": list(map(str, (mB, MB))), "policy": list(map(str, (mP, MP)))}
    return report


def main():
    report = {
        "scope": "Supporting theory only. No calibration, finite-horizon simulation, new planner power, or main-note revision.",
        "symbolic_original_equation_checks": symbolic_checks(),
        "finite_transition_certificate": finite_transition_certificate(),
        "limiting_branch_witnesses": limiting_witnesses(),
        "mixed_stationary_original_equation_checks": mixed_stationary_checks(),
    }
    report["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    evidence = [
        "latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex",
        "output/pdf/simplified_olg_amendment_proposal.pdf",
        "output/model/simplified_olg_amendments/theory_slides_misallocation.pdf",
        "output/model/simplified_olg_amendments/combined_transition_figure.pdf",
        "output/model/simplified_olg_amendments/transition_extensions.md",
        "output/model/simplified_olg_amendments/transition_extension_reviews.json",
        "code/model/tools/verify_simplified_olg_local_transition.py",
    ]
    report["evidence_sha256"] = {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in evidence
    }
    target = OUT / "transition_extension_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verification": "passed", "receipt": str(target)}, indent=2))


if __name__ == "__main__":
    main()
