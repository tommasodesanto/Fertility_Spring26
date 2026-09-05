#!/usr/bin/env python3
"""Independent small algebra checks for the two planner benchmarks.

Run from the repository root with Python 3. No numerical equilibrium or calibration solves.
The direct-allocation check reconstructs conditional household baselines from
original budgets and verifies finite compensated allocations. These constructed
households are not a full equilibrium example. Results are generated under
output/model/simplified_olg_amendments/.
The committed-transfer check constructs a complete analytical equilibrium path.
Use --original-optimizers for independent SciPy checks and --plot for its
supplemental diagnostic (Matplotlib); the previous theory figures are unchanged.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

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


def committed_cash_path_check():
    """Verify a constructed infinite clearing path, with a contraction tail.

    This is an illustrative mathematical check, not a calibration. Fertility
    is fixed household by household. Transfers are committed at both ages;
    a passive initial title owner trades the outside bond and receives the
    government's datewise residual. Its entire present-value account is kept.
    """
    par = json.loads((OUT/'transition_verification.json').read_text())['parameters']
    with (OUT/'transition_phi80.csv').open() as f:
        row = next(csv.DictReader(f))
    q,tau,phi,alpha,beta,gamma,omega = (par[k] for k in
        ('q','tau_p','phi','alpha','beta','gamma','omega_b'))
    y,b,chi,kappa,Hbar,hR = (par[k] for k in
        ('income','liquid_wealth','chi','kappa','housing_stock','renter_cap'))
    N=float(row['young_mass']); share=float(row['young_owner_share'])
    m,mr=N*share,N*(1-share)
    n,nr=float(row['owner_fertility']),float(row['renter_fertility'])
    A=1+gamma+omega; rho=1+beta*A; B=1+beta*(1+omega)
    delta=1-phi; ap=1+q*tau; d=ap-q
    def at_price(P):
        u=d*P; T=q*tau*P*Hbar/(2*N); h=b/(delta*P)
        x=(y+b+(1+q)*T-u*h-chi*n)/rho
        ho=gamma*beta*x/(q*u)
        return u,T,h,x,ho
    # Remove the saved transition row's roundoff in the stationary identities.
    lo,hi=.99*float(row['price']),1.01*float(row['price'])
    for _ in range(100):
        P=(lo+hi)/2; u,T,h,x,ho=at_price(P)
        if m*(h+ho)+2*mr*hR>Hbar: lo=P
        else: hi=P
    P=(lo+hi)/2; u,T,h,x,ho=at_price(P)
    s=h-kappa*n; co=beta*x/q
    xr=(y+b+(1+q)*T-(1+q)*u*hR-chi*nr)/B
    cr=beta*xr/q
    ay=y+b+T-x-chi*n-(delta+q*tau)*P*h
    ar=y+b+T-xr-chi*nr-u*hR
    ao=ay/q-phi*P*h/q
    R=Hbar-m*h
    mv=alpha*x/s; lam=alpha*ho/(rho*s); eta=(1+omega)/A
    assert 0<lam<1 and mv>u and R>0
    slope=m*ho*eta*rho*(mv-u)/(rho-q*alpha*ho/s)
    def branch_value(owner,n_fixed,pt,pn,ut,un,tt,tn,ell_y,ell_o):
        cap=(b+ell_y)/(delta*pt) if owner else hR
        factor=rho if owner else B
        xx=(y+b+tt+q*tn+ell_y+q*ell_o-ut*cap-chi*n_fixed-
            (0 if owner else q*un*hR))/factor
        cc=beta*xx/q
        hh=gamma*cc/un if owner else hR
        ee=omega*cc/q
        aa=y+b+tt+ell_y-xx-chi*n_fixed-((delta+q*tau)*pt if owner else ut)*cap
        assert xx>0 and cap>kappa*n_fixed and aa>0
        assert alpha*xx/(cap-kappa*n_fixed)>ut
        if owner:
            assert hh<cap and ee>pn*hh
        else:
            assert gamma*cc/hR>un
        val=(math.log(xx)+alpha*math.log(cap-kappa*n_fixed)+par['theta_initial']*math.log(n_fixed)+
             beta*(math.log(cc)+gamma*math.log(hh)+omega*math.log(ee)))
        return val,[xx,cap,aa,cc,hh,ee]
    v_o,_=branch_value(True,n,P,P,u,u,T,T,0,0)
    v_r,_=branch_value(False,nr,P,P,u,u,T,T,0,0)
    taste_cutoff=v_r-v_o
    original_optimizer_rows=[]
    def optimize_original(owner,n_fixed,pt,pn,ut,un,tt,tn,ell_y,ell_o,reference):
        import numpy as np
        from scipy.optimize import minimize
        pc=(delta+q*tau)*pt if owner else ut
        matrix=np.array([[1.,pc,1.,0.,0.,0.],
                         [0.,phi*pt/q-pn if owner else 0.,-1/q,1.,un,q]])
        rhs=np.array([y+b+tt+ell_y-chi*n_fixed,tn+ell_o])
        def obj(z):
            xx,hh,aa,cc,hh2,ee=z
            return -(math.log(xx)+alpha*math.log(hh-kappa*n_fixed)+
                     beta*(math.log(cc)+gamma*math.log(hh2)+omega*math.log(ee)))
        def jac(z):
            xx,hh,aa,cc,hh2,ee=z
            return -np.array([1/xx,alpha/(hh-kappa*n_fixed),0,beta/cc,beta*gamma/hh2,beta*omega/ee])
        constraints=[{'type':'eq','fun':lambda z:matrix@z-rhs,'jac':lambda z:matrix}]
        if owner:
            imat=np.array([[0.,1.,0.,0.,-1.,0.],[0.,0.,0.,0.,-pn,1.]])
            constraints.append({'type':'ineq','fun':lambda z:imat@z,'jac':lambda z:imat})
        cap=(b+ell_y)/(delta*pt) if owner else hR
        answer=minimize(obj,np.array(reference)*[1.005,.999,1.01,1.005,.999,1.005],jac=jac,
            method='SLSQP',constraints=constraints,
            bounds=[(1e-9,None),(kappa*n_fixed+1e-8,cap),(0,None),(1e-9,None),
                    (1e-9,None if owner else hR),(1e-9,None)],
            options={'ftol':1e-13,'maxiter':150})
        assert answer.success,answer.message
        error=float(np.max(np.abs(answer.x-reference)))
        assert error<2e-6
        original_optimizer_rows.append(dict(owner=owner,fixed_fertility=n_fixed,
            maximum_choice_discrepancy=error,budget_residual=float(np.max(np.abs(matrix@answer.x-rhs)))))
    cases=[]
    for eps,young_gain in ((1e-3,0.),(1e-4,0.),(1e-5,0.),(1e-4,2e-6)):
        horizon=80
        prices=[P+u*eps/ap]+[P]*(horizon+2)
        rents=[ap*prices[t]-q*prices[t+1] for t in range(horizon+2)]
        rebates=[q*tau*p*Hbar/(2*N) for p in prices]
        old_c0=co*(1+eps)**(gamma/A)
        old_h0=ho*(1+eps)**(-eta)
        hs=[h+ho-old_h0]
        xs=[x*(s/(hs[0]-kappa*n))**(alpha/rho)*math.exp(young_gain/rho)]
        assert xs[0]<x
        for t in range(1,horizon+1):
            hs.append(h+ho*(1-xs[t-1]/x))
            xs.append(x*(s/(hs[t]-kappa*n))**(alpha/rho))
        ly=[delta*prices[t]*hs[t]-b for t in range(horizon+1)]
        lo_cash=[A*old_c0-ao-prices[0]*h-rebates[0]]
        for t in range(horizon):
            lo_cash.append((rho*xs[t]+rents[t]*hs[t]+chi*n-y-b-
                            rebates[t]-q*rebates[t+1]-ly[t])/q)
        lr=[(rents[0]-u)*hR-(rebates[0]-T)]+[0.]*horizon
        lor=lr[:]
        residuals=[]; margins=[]; tenure_margins=[]; path=[]; pv=R*(prices[0]-P)
        for t in range(horizon):
            pt,ut,tt=prices[t],rents[t],rebates[t]
            ht,xt=hs[t],xs[t]
            hot=old_h0 if t==0 else ho*xs[t-1]/x
            cot=old_c0 if t==0 else beta*xs[t-1]/q
            eot=omega*cot/q
            initial_a=ao if t==0 else path[-1]['owner_saving']/q-phi*prices[t-1]*hs[t-1]/q
            inherited_h=h if t==0 else hs[t-1]
            saving=y+b+tt+ly[t]-xt-chi*n-(delta+q*tau)*pt*ht
            renter_saving=y+b+tt+lr[t]-xr-chi*nr-ut*hR
            passive=-m*(ly[t]+lo_cash[t])-mr*(lr[t]+lor[t])
            pv+=q**t*passive
            residuals.extend((m*(ht+hot)+2*mr*hR-Hbar,
                m*(ly[t]+lo_cash[t])+mr*(lr[t]+lor[t])+passive,
                cot+q*eot+ut*hot-initial_a-pt*inherited_h-tt-lo_cash[t],
                beta*xt/q+q*omega*beta*xt/q**2+rents[t+1]*gamma*beta*xt/(q*rents[t+1])-
                (saving/q-phi*pt*ht/q+prices[t+1]*ht+rebates[t+1]+lo_cash[t+1]),
                cr+q*omega*cr/q+ut*hR-
                ((ar/q if t==0 else path[-1]['renter_saving']/q)+tt+lor[t]),
                B*math.log(xr/xr),
                rho*math.log(xt/x)+alpha*math.log((ht-kappa*n)/s)-(young_gain if t==0 else 0),
                delta*pt*ht-b-ly[t]))
            margins.extend((saving,renter_saving,hot,inherited_h-hot,
                eot-prices[t+1]*hot,par['owner_cap']-ht,
                alpha*xt/(ht-kappa*n)-ut,gamma*cr/hR-ut))
            # The grants depend on household identity, not its chosen tenure.
            # A deviation keeps this household's own fixed fertility and cash.
            for identity_owner,n_fixed,ell_y,ell_old in (
                (True,n,ly[t],lo_cash[t+1]),(False,nr,lr[t],lor[t+1])):
                values=[]
                for choice_owner in (False,True):
                    args=(choice_owner,n_fixed,pt,prices[t+1],ut,rents[t+1],tt,rebates[t+1],ell_y,ell_old)
                    vv,ref=branch_value(*args)
                    values.append(vv)
                    if '--original-optimizers' in sys.argv and eps==1e-3 and t in (0,1,3):
                        optimize_original(*args,ref)
                tenure_margins.append(values[1]+taste_cutoff-values[0] if identity_owner else
                                      values[0]-values[1]-taste_cutoff)
            path.append(dict(date=t,price=pt,user_cost=ut,owner_housing=ht,
                old_owner_housing=hot,owner_adult_consumption=xt,
                old_owner_consumption=cot,old_owner_estate=eot,
                housing_residual=m*(ht+hot)+2*mr*hR-Hbar,
                transfer_residual=m*(ly[t]+lo_cash[t])+mr*(lr[t]+lor[t])+passive,
                estate_slack=eot-prices[t+1]*hot,retention_slack=inherited_h-hot,
                owner_saving=saving,renter_saving=renter_saving,
                young_owner_grant=ly[t],old_owner_grant=lo_cash[t],
                young_renter_grant=lr[t],old_renter_grant=lor[t],
                passive_grant=passive))
        exact_resource=-m*((1+omega)*(old_c0-co)+B*sum(q**t*(xs[t]-x) for t in range(horizon)))
        residuals.extend((pv-exact_resource,
            A*math.log(old_c0/co)-gamma*math.log1p(eps)))
        maxres=max(abs(v) for v in residuals)
        assert maxres<3e-12 and min(margins)>0 and min(tenure_margins)>0 and pv>0
        # The analytic contraction proves the infinite tail; truncation is only
        # for displayed receipts, and its bound is recorded separately.
        tail_x_bound=abs(xs[0]-x)*lam**horizon/(1-q*lam)*q**horizon
        case_slope=slope-m*B*x*(young_gain/eps)/(rho*(1-q*lam))
        assert abs(pv/eps-case_slope)<.004*abs(case_slope)
        if eps==1e-3:
            with (OUT/'committed_cash_path.csv').open('w') as f:
                writer=csv.DictWriter(f,fieldnames=path[0],lineterminator='\n');writer.writeheader();writer.writerows(path)
            if '--plot' in sys.argv:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                fig,axes=plt.subplots(3,2,figsize=(11,10),constrained_layout=True)
                short=path[:9]; dates=[z['date'] for z in short]
                def panel(ax,title,series,ylabel):
                    for label,values in series:
                        ax.plot(dates,values,marker='o',markersize=3,label=label)
                    ax.set(title=title,xlabel='Date',ylabel=ylabel)
                    ax.grid(alpha=.2);ax.legend(fontsize=8)
                panel(axes[0,0],'Prices',[(k,[100*(z[col]/base-1) for z in short])
                    for k,col,base in [('House price','price',P),('Service cost','user_cost',u)]],'% from baseline')
                panel(axes[0,1],'Housing by age',[(k,[z[col]-base for z in short])
                    for k,col,base in [('Young owners','owner_housing',h),('Old owners','old_owner_housing',ho)]],'Change per household')
                panel(axes[1,0],'Consumption and estates',[(k,[100*(z[col]/base-1) for z in short])
                    for k,col,base in [('Young adult consumption','owner_adult_consumption',x),
                        ('Old consumption','old_owner_consumption',co),('Old estate','old_owner_estate',omega*co/q)]],'% from baseline')
                panel(axes[1,1],'Cash transfers',[(k,[z[col] for z in short])
                    for k,col in [('Young owners','young_owner_grant'),('Old owners','old_owner_grant'),
                        ('Young renters','young_renter_grant'),('Passive owner (total)','passive_grant')]],'Goods')
                panel(axes[2,0],'Market and fiscal residuals',[(k,[z[col] for z in short])
                    for k,col in [('Housing','housing_residual'),('New transfer budget','transfer_residual')]],'Residual')
                panel(axes[2,1],'Private constraints',[(k,[z[col] for z in short])
                    for k,col in [('Young owner saving','owner_saving'),('Renter saving','renter_saving'),
                        ('Old estate slack','estate_slack'),('Old retention slack','retention_slack')]],'Level or slack')
                fig.suptitle('Committed transfers: supplemental verification of the analytical path\nFixed individual fertility; illustrative, not calibrated',fontsize=13)
                fig.savefig(OUT/'committed_cash_diagnostics.png',dpi=160)
                plt.close(fig)
        cases.append(dict(relative_initial_user_cost_change=eps,
            initial_young_owner_lifetime_utility_gain=young_gain,
            passive_present_value_gain=pv,resource_identity_gain=exact_resource,
            predicted_gain_derivative=case_slope,observed_gain_derivative=pv/eps,
            maximum_original_equation_residual=maxres,minimum_branch_slack=min(margins),
            minimum_taste_threshold_tenure_margin=min(tenure_margins),
            discounted_consumption_tail_bound=tail_x_bound,
            impact_young_housing_gain=hs[0]-h,impact_old_housing_change=old_h0-ho))
    # Independent linear recurrence when all transfers cease after date zero.
    ax=(T-q*P*h)/(rho*x); bx=q*(T+P*h)/(rho*x)
    qa=q/d; qb=bx-ap/d-h/ho; qc=ax
    disc=qb*qb-4*qa*qc
    roots=sorted(((-qb-math.sqrt(disc))/(2*qa),(-qb+math.sqrt(disc))/(2*qa)))
    renter_loading=-((P*hR*(ap-q*roots[0])-T)*(1+q*roots[0]))
    return dict(scope='Committed two-age cash transfers; fixed individual fertility; stationary four-group regime; passive initial title owner with outside bond access',
        baseline=dict(price=P,user_cost=u,owner_housing=h,old_owner_housing=ho,
            owner_adult_consumption=x,young_housing_value=mv,passive_initial_titles=R,
            contraction_factor=lam,owner_mass=m,renter_mass=mr),
        cases=cases,original_optimizer_checks=original_optimizer_rows,
        one_time_tail=dict(roots=roots,renter_utility_numerator_loading=renter_loading,
            qualification='Linear local tail only; not a general no-inefficiency theorem'))


def main():
    # Smoke-test the exact case function before the 27-case check.
    direct_case(.02, .5, 1e-7)
    cases = [direct_case(tau, cost, eps)
             for tau in (0., .02, .05)
             for cost in (0., .5, 1.5)
             for eps in (1e-3, 1e-5, 1e-7)]
    source = ROOT/'latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex'
    result = dict(scope='Conditional checks plus a complete analytical committed-transfer equilibrium path; no calibration or numerical equilibrium solve',
                  source=str(source.relative_to(ROOT)),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  constrained_theorem_sha256=hashlib.sha256((ROOT/'latex/JMP_DS_suggestions/simplified_olg_constrained_efficiency.tex').read_bytes()).hexdigest(),
                  input_sha256={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest()
                                for name in ('transition_verification.json','transition_phi80.csv')},
                  direct_cases=cases,
                  conditional_cash=conditional_cash_check(),
                  conditional_obstruction=conditional_obstruction_check(),
                  committed_cash_path=committed_cash_path_check(),
                  max_accounting_residual=max(c['maximum_accounting_residual'] for c in cases),
                  positive_gain_cases=sum(c['young_utility_change'] > 0 for c in cases),
                  cost_obstruction_cases=sum(c['young_utility_change'] < 0 for c in cases))
    OUT.mkdir(parents=True, exist_ok=True)
    path=OUT/'planner_benchmark_checks.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='direct_cases'}, indent=2))


if __name__ == '__main__':
    main()
