#!/usr/bin/env python3
"""Small, uncalibrated checks for the proposed OLG analytical amendments.

Run from repository root: python3 code/model/tools/verify_simplified_olg_amendments.py
Writes only output/model/simplified_olg_amendments. No equilibrium or calibration
is solved. Scalar fertility choices are solved independently by bisection.
"""
from pathlib import Path
import csv
import hashlib
from datetime import datetime, timezone
import json
import math
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

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
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


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
                   welfare_cases=len(welfare), common_cap_check=common,
                   welfare_baseline=baseline, primitive_full_branch_example=primitive_example, elapsed_seconds=time.monotonic()-started)
    (OUT/'verification_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('welfare_baseline', 'primitive_full_branch_example')}, indent=2))


if __name__ == '__main__':
    main()
