#!/usr/bin/env python3
"""Small, uncalibrated checks for the proposed OLG analytical amendments.

Run from repository root: python3 code/model/tools/verify_simplified_olg_amendments.py
Default: conditional analytical checks. Add --transitions for illustrative
equilibrium checks, or --saved-transitions to inspect existing path files without
solving again. Writes only output/model/simplified_olg_amendments. No calibration.
Scalar fertility choices are solved independently by a bracketed root method.
"""
from pathlib import Path
import csv
import hashlib
from datetime import datetime, timezone
import json
import math
import time
import sys
import signal
from dataclasses import replace, asdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, minimize, least_squares

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def household(rho, alpha, theta, chi, kappa, wealth, price, housing):
    net = wealth - price * housing
    assert net > 0
    upper = min(housing / kappa, net / chi)

    def foc(n):
        x = (net - chi * n) / rho
        s = housing - kappa * n
        return theta / n - chi / x - alpha * kappa / s

    n = brentq(foc, upper * 1e-12, upper * (1 - 1e-12), xtol=1e-14)
    x, s = (net - chi * n) / rho, housing - kappa * n
    delta = theta / n**2 + chi**2 / (rho * x**2) + alpha*kappa**2/s**2
    derivative = (alpha*kappa/s**2 - chi*price/(rho*x**2)) / delta
    return n, x, s, derivative


def primitive_cases(smoke=False):
    rows = []
    # 3*3*3*9*25=6,075 candidates; only binding caps are solved.
    for rho in ((2.0,) if smoke else (1.4, 2.0, 4.6)):
        for alpha in ((1.0,) if smoke else (0.4, 1.0, 2.0)):
            for theta in ((1.1,) if smoke else (0.3, 1.1, 3.0)):
                total = rho + alpha + theta
                threshold = alpha / (rho + alpha)
                maxshare = (alpha + theta) / total
                for share in np.linspace(0.02, maxshare * 0.98, 9):
                    for v in np.geomspace(0.02, 30, 25):
                        free_share = (alpha + theta*v/(1+v)) / total
                        if share >= free_share - 1e-8:
                            continue
                        h = share/v  # wealth=chi=kappa=1, u=v
                        n, x, s, deriv = household(rho, alpha, theta, 1, 1, 1, v, h)
                        step = 1e-5*h
                        plus = household(rho, alpha, theta, 1, 1, 1, v, h+step)[0]
                        minus = household(rho, alpha, theta, 1, 1, 1, v, h-step)[0]
                        fd = (plus-minus)/(2*step)
                        assert abs(fd-deriv) < 3e-6*max(1, abs(deriv)), (fd, deriv)
                        assert alpha*x/s > v
                        sufficient = share >= threshold
                        if sufficient:
                            assert v > alpha/rho
                            assert deriv > 0
                        rows.append(dict(rho=rho, alpha=alpha, theta=theta,
                                         share=share, v=v, housing=h, fertility=n,
                                         derivative=deriv, finite_difference=fd,
                                         primitive_sufficient=sufficient))
    assert any(r['primitive_sufficient'] for r in rows)
    assert any(r['derivative'] < 0 for r in rows)
    return rows


def primitive_full_branch_example():
    q, phi, P, y, b = .8, .7, 4., 4.8, 2.4
    beta, gamma, omega = .8, .5, 3.
    rho, alpha, theta, chi, kappa = 4.6, 1., 2., 1., 1.
    u, w, cap = (1-q)*P, y+b, b/((1-phi)*P)
    n, x, space, deriv = household(rho, alpha, theta, chi, kappa, w, u, cap)
    saving = w-(x+chi*n)-(1-phi)*P*cap
    old_wealth = (saving-phi*P*cap)/q+P*cap
    c2 = old_wealth/(1+gamma+omega)
    h2, estate = gamma*c2/u, omega*c2/q
    share, threshold = u*cap/w, alpha/(rho+alpha)
    assert share > threshold and alpha*x/space > u
    assert saving > 0 and h2 < cap and estate > P*h2
    assert np.allclose([n, x, space, c2, deriv], [1, 1, 1, 1, 19/74], atol=1e-12)
    return dict(q=q, phi=phi, P=P, property_tax=0, y=y, b=b, beta=beta,
                gamma=gamma, omega_B=omega, alpha=alpha, theta=theta,
                chi=chi, kappa=kappa, rho=rho, housing=cap, fertility=n,
                saving=saving, old_c=c2, old_h=h2, old_estate=estate,
                cap_share=share, primitive_threshold=threshold, dn_dh=deriv)


def common_cap_check():
    rho, alpha, theta, chi, kappa = 2.2, 1., 1.1, 1., .1
    q, future_rent, h, wealth_before_rent, u = .8, .5, 1.1, 5., .4

    def solve(cap):
        w = wealth_before_rent-q*q*future_rent*cap
        return household(rho, alpha, theta, chi, kappa, w, u, cap)

    n, x, s, partial = solve(h)
    delta = theta/n**2+chi**2/(rho*x**2)+alpha*kappa**2/s**2
    total = partial-q*q*future_rent*chi/(rho*x*x*delta)
    fd = (solve(h+1e-5)[0]-solve(h-1e-5)[0])/2e-5
    assert abs(fd-total) < 1e-8
    return dict(partial=partial, corrected=total, finite_difference=fd)


def welfare_cases():
    # Construct a conditional owner optimum satisfying the original constraints.
    q, P, future_P, tax, phi = .8, 2., 2., .05, .7
    beta, gamma, omega, alpha, chi, kappa, theta = .8, .5, 3., 1., 1., .1, 1.1
    x, s, n, transfer = 1., 1., 1., .1
    A = 1+gamma+omega
    rho = 1+beta*A
    u = P+q*tax*P-q*future_P
    h, c = s+kappa*n, x+chi*n
    b = (1-phi)*P*h
    wealth = rho*x+chi*n+u*h
    y = wealth-b-(1+q)*transfer
    saving = y+b+transfer-c-((1-phi)+q*tax)*P*h
    old_financial = (saving-phi*P*h)/q
    old_resources = old_financial+P*h+transfer
    c2 = old_resources/A
    h2, e2 = gamma*c2/u, omega*c2/q
    assert saving > 0 and h2 < h and e2 > future_P*h2
    assert abs(c2-beta*x/q) < 1e-12
    assert abs(theta/n-chi/x-alpha*kappa/s) < 1e-12
    assert alpha*x/s > u
    seller_h, seller_H = .8, 1.2
    seller_c = u*seller_h/gamma
    seller_e = omega*seller_c/q
    seller_a = A*seller_c-P*seller_H-transfer
    assert seller_h < seller_H and seller_e > future_P*seller_h
    cost_rate, premium = .1, .15
    rows = []
    for eps in np.geomspace(1e-7, .02, 36):
        seller_c_new = seller_c*(seller_h/(seller_h-eps))**gamma
        exact = seller_c_new-seller_c-u*eps
        for instrument, L in (("exact", exact), ("uniform_premium", premium*eps), ("market_only", 0.)):
            cost = cost_rate*eps
            loan = (P+q*tax*P)*eps+L+cost
            debt_due, sale = loan/q, future_P*eps
            c2_new = c2-(debt_due-sale)
            cj_new = seller_c+u*eps+L
            old_budget_resid = cj_new+q*seller_e+u*(seller_h-eps)-(seller_a+P*seller_H+transfer+L)
            young_gain = alpha*math.log1p(eps/s)+beta*math.log(c2_new/c2)
            seller_gain = math.log(cj_new/seller_c)+gamma*math.log1p(-eps/seller_h)
            assert abs(old_budget_resid) < 1e-12
            assert abs(q*(debt_due-sale)-(u*eps+L+cost)) < 1e-12
            assert c2_new > 0 and seller_e > future_P*(seller_h-eps)
            assert e2 > future_P*h2 and h2 < h+eps
            assert young_gain > 0
            if instrument == "exact":
                assert abs(seller_gain) < 1e-12
            elif instrument == "uniform_premium":
                assert seller_gain > 0
            else:
                assert seller_gain < 1e-14  # second-order loss, near machine zero at tiny eps
            rows.append(dict(epsilon=eps, instrument=instrument, compensation=L,
                             real_cost=cost, young_gain=young_gain, seller_gain=seller_gain,
                             debt_due=debt_due, added_resale=sale))
    baseline = dict(q=q, P=P, future_P=future_P, property_tax=tax, phi=phi,
                    beta=beta, gamma=gamma, omega_B=omega, alpha=alpha, chi=chi,
                    kappa=kappa, theta=theta, y=y, b=b, common_transfer=transfer,
                    c=c, h=h, n=n, saving=saving, old_financial=old_financial,
                    old_c=c2, old_h=h2, old_estate=e2, seller_c=seller_c,
                    seller_h=seller_h, seller_title=seller_H, seller_estate=seller_e,
                    seller_financial=seller_a, usercost=u, m_Y=alpha*x/s)
    return rows, baseline


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def original_household_check(m, tenure, prices, transfers, choice):
    """Maximize original lifetime utility with the two original dated budgets."""
    p = m.P
    P0, P1, P2 = prices
    T0, T1 = transfers
    owner = tenure == 'owner'
    u0 = P0+p.q*p.tau_p*P0-p.q*P1
    u1 = P1+p.q*p.tau_p*P1-p.q*P2
    purchase_cost = ((1-p.phi)+p.q*p.tau_p)*P0 if owner else u0
    cap = min(p.owner_cap, p.liquid_wealth/((1-p.phi)*P0)) if owner else p.renter_cap
    # Variables: usable current goods, children, current housing, saving,
    # old consumption, old housing, estate. Both budget equations are linear.
    future = choice['future_old']
    reference = np.array([choice['usable_consumption'], choice['fertility'],
                          choice['housing'], choice['saving'], future['consumption'],
                          future['housing'], future['estate']])
    young_coeff = np.array([1.,p.chi,purchase_cost,1.,0.,0.,0.])
    old_coeff = np.array([0.,0.,p.gross_return*p.phi*P0-P1 if owner else 0.,
                         -p.gross_return,1.,u1,p.q])
    matrix = np.vstack([young_coeff,old_coeff])
    rhs = np.array([p.income+p.liquid_wealth+T0,T1])

    def obj(z):
        x,n,h,saving,c2,h2,e = z
        space = h-p.kappa*n
        if min(x,n,space,c2,h2,e) <= 0:
            return 1e8
        return -(math.log(x)+p.alpha*math.log(space)+p.theta_final*math.log(n)
                 +p.beta*(math.log(c2)+p.gamma*math.log(h2)+p.omega_b*math.log(e)))

    def grad(z):
        x,n,h,saving,c2,h2,e = z
        space = h-p.kappa*n
        return -np.array([1/x,p.theta_final/n-p.alpha*p.kappa/space,
                          p.alpha/space,0,p.beta/c2,p.beta*p.gamma/h2,p.beta*p.omega_b/e])

    inequalities = [np.array([0.,-p.kappa,1.,0.,0.,0.,0.])]
    if owner:
        inequalities += [np.array([0.,0.,1.,0.,0.,-1.,0.]),
                         np.array([0.,0.,0.,0.,0.,-P2,1.])]
    imat = np.vstack(inequalities)
    bounds = [(1e-10,None),(1e-10,None),(1e-10,cap),(0,None),
              (1e-10,None),(1e-10,None if owner else p.renter_cap),(1e-10,None)]
    answer = minimize(obj,reference*np.array([1.02,.98,.99,1.02,.98,.98,1.02]),
                      jac=grad,bounds=bounds,method='SLSQP',
                      constraints=[{'type':'eq','fun':lambda z:matrix@z-rhs,'jac':lambda z:matrix},
                                   {'type':'ineq','fun':lambda z:imat@z,'jac':lambda z:imat}],
                      options={'ftol':1e-12,'maxiter':300})
    error = float(np.max(np.abs(answer.x-reference)))
    assert answer.success,answer.message
    assert error < 2e-5,(tenure,error)
    assert np.max(np.abs(matrix@answer.x-rhs)) < 1e-9
    return dict(tenure=tenure,max_choice_gap=error,utility_gap=float(obj(answer.x)-obj(reference)),
                max_budget_gap=float(np.max(np.abs(matrix@answer.x-rhs))),iterations=int(answer.nit))


def transition_checks():
    """Policy paths with no gains tax; preserve all historical builder outputs."""
    import build_simplified_olg_mixed_tenure_theory as m
    import build_simplified_olg_saddle_path_diagnostic as local
    saved = m.P
    base = replace(m.Parameters(),tau_g=0.,theta_final=.35)
    started = time.monotonic()
    def stop_for_time(signum, frame):
        raise TimeoutError('Ten-minute transition-check budget reached.')
    previous_alarm_handler = signal.signal(signal.SIGALRM, stop_for_time)
    signal.alarm(600)
    summaries, originals, rows_all, systems = [], [], {}, []
    try:
        m.P = base
        initial = m.solve_steady_state(base.theta_initial,.34)
        for phi in (.80,.81,.85,.90):
            if time.monotonic()-started > 540:
                raise RuntimeError('Transition verification reached its nine-minute stop.')
            case_start = time.monotonic()
            m.P = replace(base,phi=phi)
            final = m.solve_steady_state(base.theta_initial,initial['price'])
            tr = m.transition(initial,final,28)
            rows = tr['rows']; rows_all[phi] = rows
            # Independently reconstruct clearing and both cohort laws.
            max_housing = max(abs(r['young_mass']*r['average_young_housing']+
                                  r['old_mass']*r['average_old_housing']-base.housing_stock) for r in rows)
            max_fiscal = max(abs((r['young_mass']+r['old_mass'])*r['transfer']-
                                 base.q*base.tau_p*r['price']*base.housing_stock) for r in rows)
            for first,second in zip(rows[:-1],rows[1:]):
                assert abs(second['old_mass']-first['young_mass']) < 1e-12
                assert abs(second['young_mass']-base.nu*first['average_fertility']*first['young_mass']) < 1e-12
            assert max_housing < 3e-8 and max_fiscal < 1e-9
            # Compare the analytical choices to a separate constrained optimizer.
            prices = [r['price'] for r in rows]+[final['price']]*2
            transfers = [r['transfer'] for r in rows]+[final['transfer']]
            for date in (0,2,28):
                choices = m.young_mixture(base.theta_initial,*prices[date:date+3],*transfers[date:date+2])
                for tenure in ('renter','owner'):
                    check=original_household_check(m,tenure,prices[date:date+3],transfers[date:date+2],choices[tenure])
                    originals.append(dict(phi=phi,date=date,**check))
            j0,j1=local.pencil(final,False,False)
            jr,kr,meta=local.eliminate_static_transfer(j0,j1)
            stability=local.analyze_pencil(j0,j1)
            # Second derivative method: compare centered differences at half step.
            old_step=local.STEP;local.STEP=old_step/2
            ja,jb=local.pencil(final,False,False);local.STEP=old_step
            derivative_error=max(float(np.max(np.abs(ja-j0))),float(np.max(np.abs(jb-j1))))
            systems.append(dict(phi=phi,**meta,**stability,derivative_step_gap=derivative_error))
            assert stability['stable_count']==7 and stability['finite_unstable_count']==2
            assert meta['reduced_next_jacobian_rank']==9
            assert stability['stable_projection_min_singular_value']>1e-4
            terminal = tr['terminal_state']
            terminal_target = {'young_mass':final['cohort_mass'],'old_mass':final['cohort_mass'],
                               'previous_owner_share':final['young']['owner_share'],
                               'renter_financial_wealth':final['young']['renter']['next_financial_wealth'],
                               'owner_financial_wealth':final['young']['owner']['next_financial_wealth'],
                               'owner_housing':final['young']['owner']['housing']}
            terminal_gap=max(abs(terminal[k]-v)/max(1,abs(v)) for k,v in terminal_target.items())
            summary=dict(phi=phi,price_initial=initial['price'],price_final=final['price'],
                         initial_population=2*initial['cohort_mass'],final_population=2*final['cohort_mass'],
                         population_change_pct=100*(final['cohort_mass']/initial['cohort_mass']-1),
                         impact_fertility=rows[0]['average_fertility'],
                         fertility_min=min(r['average_fertility'] for r in rows),
                         fertility_max=max(r['average_fertility'] for r in rows),
                         max_housing_gap=max_housing,max_fiscal_gap=max_fiscal,
                         max_solver_residual=tr['max_residual'],terminal_relative_gap=terminal_gap,
                         minimum_owner_retention_slack=min(r['old_owner_retention_slack'] for r in rows),
                         minimum_owner_estate_slack=min(r['old_owner_estate_slack'] for r in rows),
                         seconds=time.monotonic()-case_start)
            summaries.append(summary)
            write_csv(OUT/f'transition_phi{round(phi*100)}.csv',rows)
            (OUT/'transition_progress.json').write_text(json.dumps(summaries,indent=2)+'\n')
            print(json.dumps(summary),flush=True)
            if phi==.85:
                selected_final, selected_tr = final,tr
        m.P=replace(base,phi=.85)
        long=m.transition(initial,selected_final,40)
        horizon_gap=max(abs(long['rows'][i][key]-selected_tr['rows'][i][key])
                        for i in range(29) for key in ('price','transfer','average_fertility','young_mass','old_mass'))
        assert horizon_gap < 2e-6,horizon_gap
        write_csv(OUT/'transition_phi85_h40.csv',long['rows'])

        # A different path solver, with an independently written residual loop.
        # This loop writes housing and property-tax balance independently.
        def second_residual(log_paths):
            horizon=28
            prices=np.r_[np.exp(log_paths[:horizon+1]),[selected_final['price']]*2]
            transfers=np.r_[np.exp(log_paths[horizon+1:]),selected_final['transfer']]
            Y=O=initial['cohort_mass'];pi=initial['young']['owner_share']
            renter,owner=initial['young']['renter'],initial['young']['owner']
            market,fiscal=[],[]
            for t in range(horizon+1):
                old_r=m.old_renter_allocation(renter['next_financial_wealth'],prices[t],prices[t+1],transfers[t])
                old_o=m.old_owner_allocation(owner['next_financial_wealth'],owner['housing'],prices[t],prices[t],prices[t+1],transfers[t])
                young=m.young_mixture(base.theta_initial,*prices[t:t+3],*transfers[t:t+2])
                used=Y*young['average_housing']+O*((1-pi)*old_r['housing']+pi*old_o['housing'])
                market.append(used/base.housing_stock-1)
                fiscal.append(((Y+O)*transfers[t]-base.q*base.tau_p*prices[t]*base.housing_stock)/
                              (base.q*base.tau_p*selected_final['price']*base.housing_stock))
                Y,O=base.nu*young['average_fertility']*Y,Y
                pi,renter,owner=young['owner_share'],young['renter'],young['owner']
            return np.array(market+fiscal)

        reference=selected_tr['log_paths']
        assert np.max(np.abs(second_residual(reference)))<1e-8
        alternate=least_squares(second_residual,reference+.001*np.sin(np.arange(len(reference))),
                                xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=25)
        second_gap=float(np.max(np.abs(alternate.x-reference)))
        assert alternate.success and second_gap<2e-6,(alternate.message,second_gap)
        second_max=float(np.max(np.abs(second_residual(alternate.x))))
        assert second_max<1e-8
        write_csv(OUT/'original_household_optimizer_checks.csv',originals)
        (OUT/'transition_stability_checks.json').write_text(json.dumps(systems,indent=2)+'\n')
        result=dict(status='passed',scope='Illustrative credit reforms; zero gains tax, property tax retained',
                    parameters=asdict(base),policy_summaries=summaries,horizon28_vs40_gap=horizon_gap,
                    alternate_solver_log_path_gap=second_gap,alternate_solver_residual=second_max,
                    original_optimizer_max_choice_gap=max(r['max_choice_gap'] for r in originals),
                    original_optimizer_checks=len(originals),elapsed_seconds=time.monotonic()-started)
        (OUT/'transition_verification.json').write_text(json.dumps(result,indent=2)+'\n')
        fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
        for phi in (.80,.81,.85,.90):
            rows=rows_all[phi];t=[r['date'] for r in rows];label=f'{phi:.0%} financed'
            axes[0,0].plot(t,[r['price'] for r in rows],label=label)
            axes[0,1].plot(t,[r['average_fertility'] for r in rows])
            axes[0,2].plot(t,[100*(r['young_mass']+r['old_mass'])/(2*initial['cohort_mass']) for r in rows])
        selected=rows_all[.85];t=[r['date'] for r in selected]
        axes[1,0].plot(t,[r['average_young_housing'] for r in selected],label='Young')
        axes[1,0].plot(t,[r['average_old_housing'] for r in selected],label='Old')
        axes[1,1].plot(t,[r['young_owner_share'] for r in selected],label='Young')
        axes[1,1].plot(t,[r['previous_owner_share'] for r in selected],label='Old')
        axes[1,2].semilogy(t,[max(abs(r['housing_residual']),1e-16) for r in selected],label='Housing')
        axes[1,2].semilogy(t,[max(abs(r['government_residual']),1e-16) for r in selected],label='Government')
        titles=['House price','Fertility (model child units)','Adult households (initial = 100)',
                'Housing per household: 85% policy','Owner share: 85% policy','Equation errors: 85% policy']
        for ax,title in zip(axes.flat,titles):ax.set_title(title,fontsize=10);ax.set_xlabel('Generations after reform')
        for ax in (axes[0,0],axes[1,0],axes[1,1],axes[1,2]):ax.legend(fontsize=8)
        fig.suptitle('Credit access in the illustrative OLG economy; property tax remains 5%',fontsize=12)
        fig.savefig(OUT/'credit_policy_transitions.pdf');fig.savefig(OUT/'credit_policy_transitions.png',dpi=160)
        plt.close(fig)
        print('Transition checks passed; '+json.dumps({k:result[k] for k in ('elapsed_seconds','horizon28_vs40_gap','alternate_solver_log_path_gap','original_optimizer_max_choice_gap')}),flush=True)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_alarm_handler)
        m.P=saved


def accounting_checks():
    """Check child rescaling and stationary supply using independent arithmetic."""
    rho, alpha, theta, chi, kappa, wealth, u, h = 2., 1., 1.1, 1., .1, 3.55, .5, 1.1
    n, x, s, _ = household(rho, alpha, theta, chi, kappa, wealth, u, h)
    errors=[]
    for scale in (.3, 2.1, 4.2, 10.):
        nL, xL, sL, _ = household(rho, alpha, theta, chi/scale, kappa/scale, wealth, u, h)
        errors.extend([abs(nL-scale*n), abs(xL-x), abs(sL-s),
                       abs((2/scale)*nL-2*n)])
        utility_gap = rho*math.log(xL/x)+alpha*math.log(sL/s)+theta*math.log(nL/n)
        assert abs(utility_gap-theta*math.log(scale))<1e-12
    assert max(errors)<1e-12
    # Independent terminal household, person, fiscal, and supply identities.
    records=[]
    for stock,price,hY,hO in ((2.,.4,.8,.5),(3.,.5,.9,.4)):
        cohort=stock/(hY+hO)
        households=2*cohort
        persons=2*cohort+2*cohort+4.2*.5*cohort
        transfer=.5*.05*price*stock/households
        assert abs(persons/households-(2+2+4.2/2)/2)<1e-12
        assert abs(transfer-.5*.05*price*(hY+hO)/2)<1e-12
        records.append(dict(stock=stock,price=price,households=households,mean_housing=(hY+hO)/2))
    lhs=math.log(records[1]['households']/records[0]['households'])
    rhs=math.log(records[1]['stock']/records[0]['stock'])-math.log(records[1]['mean_housing']/records[0]['mean_housing'])
    assert abs(lhs-rhs)<1e-12
    return dict(status='passed',child_scales=[.3,2.1,4.2,10.],
                maximum_rescaling_error=max(errors),supply_identity_error=abs(lhs-rhs))


def saved_transition_checks():
    """Use saved equilibrium paths for the report and population accounting."""
    paths={}
    for phi in (80,81,85,90):
        with (OUT/f'transition_phi{phi}.csv').open() as f:
            paths[phi]=[{k:(v if k.endswith('_branch') else float(v)) for k,v in r.items()}
                        for r in csv.DictReader(f)]
    base=paths[80][0]
    initial=base['young_mass']+base['old_mass']
    decompositions=[];max_product_gap=0.;max_clearing_gap=0.
    for phi,rows in paths.items():
        cumulative=1.
        assert abs(rows[0]['young_mass']-base['young_mass'])<1e-12
        assert abs(rows[0]['old_mass']-base['old_mass'])<1e-12
        for t,r in enumerate(rows):
            max_product_gap=max(max_product_gap,abs(r['young_mass']/paths[80][t]['young_mass']-cumulative))
            max_clearing_gap=max(max_clearing_gap,abs(r['young_mass']*r['average_young_housing']+
                                      r['old_mass']*r['average_old_housing']-2.))
            cumulative*=r['average_fertility']/paths[80][t]['average_fertility']
        r=rows[0]; share0=base['young_owner_share'];share=r['young_owner_share']
        within=(1-share0)*(r['renter_fertility']-base['renter_fertility'])+share0*(r['owner_fertility']-base['owner_fertility'])
        composition=(share-share0)*(r['owner_fertility']-r['renter_fertility'])
        actual=r['average_fertility']-base['average_fertility']
        assert abs(within+composition-actual)<1e-12
        decompositions.append(dict(phi=phi/100,impact_fertility_change=actual,
                                   within_tenure_at_initial_shares=within,tenure_composition=composition))
    assert max_product_gap<1e-8 and max_clearing_gap<1e-8
    assert paths[90][0]['average_fertility']<base['average_fertility']
    result=dict(status='passed',maximum_cohort_product_error=max_product_gap,
                maximum_housing_identity_error=max_clearing_gap,
                impact_decomposition=decompositions)
    (OUT/'population_and_composition_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    # A compact paper figure; retain the full diagnostic packet unchanged.
    fig,axes=plt.subplots(1,2,figsize=(8.,3.1),constrained_layout=True)
    for phi,rows in paths.items():
        t=[r['date'] for r in rows if r['date']<=12]
        short=rows[:len(t)]
        axes[0].plot(t,[r['average_fertility'] for r in short],label=f'{phi}% financed')
        axes[1].plot(t,[100*(r['young_mass']+r['old_mass'])/initial for r in short])
    axes[0].set_title('Fertility (model child units)',fontsize=11)
    axes[1].set_title('Adult households (initial = 100)',fontsize=11)
    for ax in axes:
        ax.set_xlabel('Generations after reform')
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',alpha=.15)
    axes[0].legend(fontsize=10,frameon=False)
    fig.savefig(OUT/'credit_policy_summary.pdf')
    fig.savefig(OUT/'credit_policy_summary.png',dpi=180)
    plt.close(fig)
    print(json.dumps(result,indent=2))


def main():
    started = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    # Smoke-test the exact candidate loop function before its fixed-size grid.
    n, x, s, negative = household(2, 1, 1.1, 1, .1, 3.55, .5, 1.1)
    assert np.allclose([n, x, s, negative], [1, 1, 1, -15/161], atol=1e-12)
    primitive_cases(smoke=True)
    fertility = primitive_cases()
    welfare, baseline = welfare_cases()
    common = common_cap_check()
    primitive_example = primitive_full_branch_example()
    accounting = accounting_checks()
    write_csv(OUT/'fertility_checks.csv', fertility)
    write_csv(OUT/'welfare_checks.csv', welfare)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
    selected = [r for r in fertility if r['rho']==2. and r['alpha']==1. and r['theta']==1.1]
    ax = axes[0]
    pos = [r for r in selected if r['derivative'] > 0]
    neg = [r for r in selected if r['derivative'] <= 0]
    for group, label, color in ((pos, 'Fertility increases', '#167c80'), (neg, 'Fertility decreases', '#bd593e')):
        ax.scatter([r['v'] for r in group], [r['share'] for r in group], s=15, color=color, label=label)
    ax.axhline(1/3, color='black', linestyle='--', linewidth=1, label='Primitive sufficient threshold')
    ax.set(xscale='log', xlabel=r'Space cost relative to goods cost, $\kappa u/\chi$',
           ylabel='Housing spending at cap / lifetime resources', title='Binding-cap household choices')
    ax.legend(fontsize=7, loc='lower right')
    ax = axes[1]
    for instrument, label in (('exact', 'Exact compensation'), ('uniform_premium', 'Uniform seller premium')):
        group = [r for r in welfare if r['instrument']==instrument]
        ax.plot([r['epsilon'] for r in group], [r['young_gain'] for r in group], label=label)
    ax.set(xlabel='Housing transferred', ylabel='Young lifetime utility gain', title='Feasible compensated reallocations')
    ax.legend(fontsize=8)
    fig.suptitle('Analytical checks: constructed examples, not equilibrium or calibration', fontsize=11)
    fig.savefig(OUT/'analytical_checks.pdf')
    fig.savefig(OUT/'analytical_checks.png', dpi=170)
    plt.close(fig)
    summary = dict(status='passed', verified_at_utc=datetime.now(timezone.utc).isoformat(),
                   verification_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   proposal_sha256=hashlib.sha256((ROOT/'latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex').read_bytes()).hexdigest(), scope='Conditional analytical checks; no GE or calibrated result',
                   binding_cases=len(fertility), primitive_sufficient_cases=sum(int(r['primitive_sufficient']) for r in fertility),
                   negative_response_cases=sum(int(r['derivative']<0) for r in fertility),
                   max_finite_difference_error=max(abs(r['derivative']-r['finite_difference']) for r in fertility),
                   welfare_cases=len(welfare), common_cap_check=common, accounting_checks=accounting,
                   welfare_baseline=baseline, primitive_full_branch_example=primitive_example, elapsed_seconds=time.monotonic()-started)
    (OUT/'verification_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('welfare_baseline', 'primitive_full_branch_example')}, indent=2))


if __name__ == '__main__':
    main()
    if '--transitions' in sys.argv:
        transition_checks()
    if '--transitions' in sys.argv or '--saved-transitions' in sys.argv:
        saved_transition_checks()
