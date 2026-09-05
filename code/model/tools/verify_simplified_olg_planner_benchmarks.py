#!/usr/bin/env python3
"""Independent small algebra checks for the two planner benchmarks.

Run from the repository root with Python 3. No equilibrium or calibration solves.
The direct-allocation check reconstructs conditional household baselines from
original budgets and verifies finite compensated allocations. These constructed
households are not a full equilibrium example. Results are generated under
output/model/simplified_olg_amendments/.
"""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def baseline(tau):
    q, phi, P, P_next = .8, .7, 4., 4.
    y, b, beta, alpha, gamma, omega = 4.8, 2.4, .8, 1., .5, 3.
    chi, kappa, theta = 1., 1., 2.
    T = T_next = 2*q*tau*P
    u = (1+q*tau)*P-q*P_next
    h = b/((1-phi)*P)
    w = y+b+T+q*T_next
    factor = 1+beta*(1+gamma+omega)
    lo, hi = 1e-12, min(h/kappa, (w-u*h)/chi)-1e-12
    for _ in range(100):
        n = (lo+hi)/2
        x = (w-u*h-chi*n)/factor
        foc = theta/n-chi/x-alpha*kappa/(h-kappa*n)
        if foc > 0:
            lo = n
        else:
            hi = n
    n = (lo+hi)/2
    x = (w-u*h-chi*n)/factor
    c = x+chi*n
    saving = y+b+T-c-((1-phi)+q*tau)*P*h
    a_next = saving/q-phi*P*h/q
    old_resources = a_next+P_next*h+T_next
    c_next = old_resources/(1+gamma+omega)
    h_next = gamma*c_next/u
    e_next = omega*c_next/q
    # Interior old seller from the same preference and price environment.
    h_old, H_old = 1., 2.
    c_old = u*h_old/gamma
    e_old = omega*c_old/q
    a_old = c_old+q*e_old+u*h_old-P*H_old-T
    assert saving > 0 and 0 < h_next < h
    assert e_next > P_next*h_next and e_old > P_next*h_old
    assert abs(c_next-beta*x/q) < 2e-14
    assert alpha*x/(h-kappa*n) > u
    return locals()


def direct_case(tau, cost_ratio, eps):
    z = baseline(tau)
    q, phi, P, P_next = (z[k] for k in ('q','phi','P','P_next'))
    x, s = z['x'], z['h']-z['kappa']*z['n']
    mv_y = z['alpha']*x/s
    mv_o = z['gamma']*z['c_old']/z['h_old']
    K = cost_ratio*(mv_y-mv_o)
    D = z['c_old']*math.expm1(-z['gamma']*math.log1p(-eps/z['h_old']))
    C = K*eps
    # L is the net cash compensation above sale/user-cost receipts.
    L = D-z['u']*eps
    c_new = z['c']-D-C
    h_new = z['h']+eps
    saving_new = z['saving']+phi*P*eps-q*P_next*eps
    a_next_new = saving_new/q-phi*P*h_new/q
    current_young = c_new+saving_new+((1-phi)+q*tau)*P*h_new-(z['y']+z['b']+z['T']-L-C)
    future_young = a_next_new+P_next*h_new+z['T_next']-z['old_resources']
    current_old = z['c_old']+D+q*z['e_old']+z['u']*(z['h_old']-eps)-(z['a_old']+P*z['H_old']+z['T']+L)
    goods = c_new+z['c_old']+D+C-z['c']-z['c_old']
    housing = h_new+z['h_old']-eps-z['h']-z['h_old']
    du_old = math.log1p(D/z['c_old'])+z['gamma']*math.log1p(-eps/z['h_old'])
    du_young = math.log1p(-(D+C)/x)+z['alpha']*math.log1p(eps/s)
    derivative = (mv_y-mv_o-K)/x
    residual = max(abs(v) for v in (current_young, future_young, current_old, goods, housing, du_old))
    assert residual < 3e-14
    assert c_new > z['chi']*z['n'] and saving_new > 0
    assert z['h_next'] < h_new and z['e_next'] > P_next*z['h_next']
    assert z['e_old'] > P_next*(z['h_old']-eps)
    assert (1-phi)*P*h_new > z['b']  # Requires relaxed private purchase finance.
    assert (du_young > 0) == (cost_ratio < 1)
    if eps == 1e-7:
        assert abs(du_young/eps-derivative) < 2e-6
    return dict(tau_p=tau, cost_ratio=cost_ratio, epsilon=eps,
                young_utility_change=du_young, old_utility_change=du_old,
                predicted_derivative=derivative, observed_slope=du_young/eps,
                maximum_accounting_residual=residual)



def conditional_cash_check():
    """Exact four-group exercise; fixed tenure, fertility and future prices.

    This omits future market clearing and intermediary-owner welfare. It is
    NOT a constrained Pareto theorem for the complete economy.
    """
    parameters = json.loads((OUT/'transition_verification.json').read_text())['parameters']
    with (OUT/'transition_phi80.csv').open() as f:
        raw = next(csv.DictReader(f))
    r = {k:float(v) for k,v in raw.items() if k not in ('renter_old_branch','owner_old_branch')}
    q, tau, phi = (parameters[k] for k in ('q','tau_p','phi'))
    beta, alpha, gamma, omega = (parameters[k] for k in ('beta','alpha','gamma','omega_b'))
    y, b, chi, kappa = (parameters[k] for k in ('income','liquid_wealth','chi','kappa'))
    A = 1+gamma+omega
    factor = 1+beta*A
    P0, Pnext, T0 = r['price'], r['next_price'], r['transfer']
    # The saved row is stationary to numerical tolerance.
    Tnext = T0
    Hbar, hR = parameters['housing_stock'], parameters['renter_cap']
    Y, O = r['young_mass'], r['old_mass']
    masses = dict(YO=Y*r['young_owner_share'], YR=Y*(1-r['young_owner_share']),
                  OO=O*r['previous_owner_share'], OR=O*(1-r['previous_owner_share']))
    n, h0, ho0 = r['owner_fertility'], r['owner_housing'], r['old_owner_housing']
    u0 = r['owner_user_cost']
    x0 = (y+b+T0+q*Tnext-u0*h0-chi*n)/factor
    co0 = u0*ho0/gamma
    Wo0 = A*co0
    H_old = h0  # stationary inherited owner title
    ao = Wo0-P0*H_old-T0
    ar = r['renter_saving']/q
    Wr0 = ar+T0-u0*hR
    a = 1+q*tau
    gap = alpha*x0/(h0-kappa*n)-u0
    B = gap/((1-phi)*P0)
    kYO = (a*h0-T0/P0+B*(1-phi)*h0)/(1+B)
    kr = a*hR-T0/P0
    kOO = a*ho0-H_old-T0/P0
    C = masses['YO']*kYO+(masses['YR']+masses['OR'])*kr+masses['OO']*kOO
    dYO = (kYO-(1-phi)*h0)/((1-phi)*P0)
    dOO = -(1-gamma/A)*ho0*a/u0
    D = masses['YO']*dYO+masses['OO']*dOO
    etaOO = gamma/(A*u0)
    assert C > 0 and 0 < D/C < etaOO
    rows=[]
    for step in (1e-4, 1e-5, 1e-6):
        P = P0-step
        T = q*tau*P*Hbar/(Y+O)
        u = a*P-q*Pnext
        def owner_value(h):
            ell = (1-phi)*P*h-b
            x = (y+b+T+q*Tnext+ell-u*h-chi*n)/factor
            return factor*math.log(x/x0)+alpha*math.log((h-kappa*n)/(h0-kappa*n))
        lo, hi = .99*h0, 1.01*h0
        assert owner_value(lo) < 0 < owner_value(hi)
        for _ in range(90):
            mid=(lo+hi)/2
            if owner_value(mid)>0:
                hi=mid
            else:
                lo=mid
        h=(lo+hi)/2
        ellYO=(1-phi)*P*h-b
        ellYR=(u-u0)*hR-(T-T0)
        ho=(Hbar-masses['YO']*h-(masses['YR']+masses['OR'])*hR)/masses['OO']
        Wo=A*u*ho/gamma
        ellOO=Wo-ao-P*H_old-T
        ellOR=-(masses['YO']*ellYO+masses['YR']*ellYR+masses['OO']*ellOO)/masses['OR']
        Wr=ar+T+ellOR-u*hR
        duOO=A*math.log(Wo/Wo0)-gamma*math.log(u/u0)
        duOR=(1+omega)*math.log(Wr/Wr0)
        # Current household optimization and positivity conditions.
        x=(y+b+T+q*Tnext+ellYO-u*h-chi*n)/factor
        saving=y+b+T+ellYO-(x+chi*n)-((1-phi)+q*tau)*P*h
        young_old_resources=saving/q-phi*P*h/q+Pnext*h+Tnext
        young_old_c=young_old_resources/A
        young_old_h=gamma*young_old_c/u0
        young_old_e=omega*young_old_c/q
        assert saving>0 and 0 < young_old_h < h < parameters['owner_cap']
        assert young_old_e > Pnext*young_old_h
        assert 0 < ho < H_old and omega*(Wo/A)/q > Pnext*ho
        assert gamma/(hR) > u*(1+omega)/Wr  # old rental cap binds
        assert duOO > 0 and duOR > 0
        residual=max(abs(owner_value(h)),
                     abs(masses['YO']*ellYO+masses['YR']*ellYR+masses['OO']*ellOO+masses['OR']*ellOR),
                     abs(masses['YO']*h+masses['OO']*ho+(masses['YR']+masses['OR'])*hR-Hbar))
        assert residual < 2e-13
        assert h < h0 and ho > ho0
        rows.append(dict(price_change=-step, grants=dict(YO=ellYO,YR=ellYR,OO=ellOO,OR=ellOR),
                         young_owner_housing=h,old_owner_housing=ho,
                         old_owner_utility_gain=duOO,old_renter_utility_gain=duOR,
                         maximum_residual=residual))
    return dict(scope='Four household groups only; fixed fertility/tenure/future prices; intermediary ownership omitted',
                compensation_sum=C,compensated_housing_slope=D,old_housing_cash_slope=etaOO,cases=rows)



def conditional_obstruction_check():
    """A conditional current-date economy, not a stationary OLG equilibrium."""
    q,tau,phi,P,Pnext = .5,.05,.8,1.,1.
    alpha,beta,gamma,omega,chi,kappa = .4,.4,.3,.4,.15,.5
    y,b,theta,T,Tnext = 3.73375,.2,41/240,.0175,.0175
    A=1+gamma+omega
    u=(1+q*tau)*P-q*Pnext
    n,h,x,saving=.5,1.,2.,1.65125
    assert abs(x+chi*n+saving+((1-phi)+q*tau)*P*h-(y+b+T)) < 1e-14
    W=saving/q-phi*P*h/q+Pnext*h+Tnext
    c2=W/A; h2=gamma*c2/u; e2=omega*c2/q
    assert abs(c2-1.6)<1e-14 and h2<h and e2>Pnext*h2
    assert abs(theta/n-chi/x-alpha*kappa/(h-kappa*n)) < 1e-14
    hR=.5
    factor=1+beta*(1+omega)
    wR=y+b+T+q*Tnext-q*u*hR
    lo,hi=1e-10,hR/kappa-1e-10
    for _ in range(90):
        nr=(lo+hi)/2; xr=(wR-u*hR-chi*nr)/factor
        f=theta/nr-chi/xr-alpha*kappa/(hR-kappa*nr)
        if f>0: lo=nr
        else: hi=nr
    nr=(lo+hi)/2; xr=(wR-u*hR-chi*nr)/factor
    ar=y+b+T-xr-chi*nr-u*hR
    cr2=(ar/q+Tnext-u*hR)/(1+omega)
    assert ar>0 and gamma*cr2/u>hR
    aO,HO,aR=1.3625,1.,1.9825
    cO=(aO+P*HO+T)/A;hO=gamma*cO/u;eO=omega*cO/q
    cR=(aR+T-u*hR)/(1+omega);eR=omega*cR/q
    assert abs(cO-1.4)<1e-14 and abs(hO-.8)<1e-14
    assert hO<HO and eO>Pnext*hO and gamma*cR/u>hR
    assert abs(.5*(h+hR+hO+hR)-1.4)<1e-14
    a=1+q*tau;delta=1-phi;t=T/P
    gap=alpha*x/(h-kappa*n)-u;B=gap/(delta*P)
    kYO=(a*h-t+B*delta*h)/(1+B)
    kr=a*hR-t;kOO=a*hO-HO-t
    C=.5*(kYO+2*kr+kOO)
    D=.5*((kYO-delta*h)/(delta*P)-(1-gamma/A)*hO*a/u)
    assert gap>0 and C>0 and D<0
    return dict(scope='Current-date conditional economy; inherited old states not asserted stationary; transfer account covers four household groups',
                q=q,tau_p=tau,phi=phi,price=P,young_housing_value_gap=gap,
                compensation_sum=C,compensated_housing_slope=D,
                young_renter_fertility=nr,young_renter_adult_consumption=xr,
                young_renter_saving=ar,young_owner_future=dict(c=c2,h=h2,e=e2),
                old_owner=dict(c=cO,h=hO,e=eO),old_renter=dict(c=cR,h=hR,e=eR))


def main():
    # Smoke-test the exact case function before the 27-case check.
    direct_case(.02, .5, 1e-7)
    cases = [direct_case(tau, cost, eps)
             for tau in (0., .02, .05)
             for cost in (0., .5, 1.5)
             for eps in (1e-3, 1e-5, 1e-7)]
    source = ROOT/'latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex'
    result = dict(scope='Constructed conditional households; no equilibrium or calibration solve',
                  source=str(source.relative_to(ROOT)),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  input_sha256={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest()
                                for name in ('transition_verification.json','transition_phi80.csv')},
                  direct_cases=cases,
                  conditional_cash=conditional_cash_check(),
                  conditional_obstruction=conditional_obstruction_check(),
                  max_accounting_residual=max(c['maximum_accounting_residual'] for c in cases),
                  positive_gain_cases=sum(c['young_utility_change'] > 0 for c in cases),
                  cost_obstruction_cases=sum(c['young_utility_change'] < 0 for c in cases))
    OUT.mkdir(parents=True, exist_ok=True)
    path=OUT/'planner_benchmark_checks.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='direct_cases'}, indent=2))


if __name__ == '__main__':
    main()
