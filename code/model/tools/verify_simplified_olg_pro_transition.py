#!/usr/bin/env python3
"""Replay Pro's exact certificate and compare it with the original OLG budgets.

The nine household evaluations are arithmetic/optimality checks, not an
equilibrium simulation. No quantitative model or manuscript is modified.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import itertools
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"
sys.path.insert(0, str(ROOT / "code/model/tools"))
from verify_simplified_olg_mixed_transition import anchor
from verify_simplified_olg_local_transition import (
    young_choices, old_choices, check_path, complex_jacobian,
)
import simplified_olg_pro_certificate as pro
SRC = Path(pro.__file__)

def mid(z):
    z = z.v if isinstance(z, pro.V) else z
    return float((z.lo + z.hi) / 2)

def run_original_checks():
    p0, _, initial = anchor(4)
    primitive_names = {
        "q": "q", "beta": "beta", "alpha": "alpha", "omega": "omega",
        "gamma": "gamma", "chi": "chi", "kap": "kappa",
        "theta0": "theta", "phi0": "phi", "b": "b", "y": "y",
        "stock": "Hbar", "tax": "tau", "a": "rental_cap",
        "nu": "nu", "sigma": "sigma",
    }
    assert all(float(getattr(pro, name)) == p0[key]
               for name, key in primitive_names.items())
    assert p0["owner_cap"] == 2
    r = Q(1, 2000)
    shock = Q(1, 20000)
    cases = [(Q(0), Q(0), 0)]
    cases += list(itertools.product((Q(0), shock), (Q(0), shock), (-1, 1)))
    rows = []
    maximum_policy_error = maximum_tenure_error = 0.
    for case, (dt, df, sign) in enumerate(cases):
        pp = [1 + sign*r, 1 - sign*r, 1 + sign*r]
        yy = {-1: 1-sign*r, 0: 1+sign*r, 1: 1-sign*r}
        ts = [pro.taxcoef*pp[0]/(yy[0]+yy[-1]),
              pro.taxcoef*pp[1]/(yy[1]+yy[0])]
        p = dict(p0, theta=float(pro.theta0-dt), phi=float(pro.phi0+df))
        own = young_choices(list(map(float, pp)), list(map(float, ts)), p)
        exact = pro.young_iv(0, lambda k:pro.V(pp[k]), lambda k:pro.V(yy[k]),
                             pro.V(pro.theta0-dt), pro.V(pro.phi0+df))
        u, un = mid(exact["u"]), mid(exact["un"])
        for tenure, key in (("owner", "o"), ("renter", "r")):
            x, n = mid(exact["x"+key]), mid(exact["n"+key])
            h = mid(exact["h"]) if key=="o" else float(pro.a)
            cost = (((1-p["phi"])+p["q"]*p["tau"])*float(pp[0])*h
                    if key=="o" else u*h)
            saving = p["y"]+p["b"]+float(ts[0])-x-p["chi"]*n-cost
            c2 = p["beta"]*x/p["q"]
            h2 = p["gamma"]*c2/un if key=="o" else float(pro.a)
            estate = p["omega"]*c2/p["q"]
            z = np.array([x, h, n, saving, c2, h2, estate])
            maximum_policy_error = max(maximum_policy_error, float(max(abs(z-own[tenure]["z"]))))
        maximum_tenure_error = max(maximum_tenure_error, abs(mid(exact["pi"])-own["pi"]))
        old = {m:old_choices(initial[m]["assets"], initial[m]["z"][1],
                              float(pp[0]), float(pp[1]), float(ts[0]), p, m=="owner")
               for m in ("owner", "renter")}
        row = dict(t=case, prices=list(map(float,pp)), transfers=list(map(float,ts)),
                   price=float(pp[0]), transfer=float(ts[0]), hh=own,
                   old_choices=old, past=initial)
        result = check_path({"rows":[row]}, p, optimize_dates=(0,8))
        rows.append(dict(case=case, theta_decline=str(dt), credit_change=str(df),
                         price_population_corner=sign, original_checks=result))
    assert maximum_policy_error < 5e-14 and maximum_tenure_error < 5e-14

    # Independent derivative calculation from the existing full household utilities.
    # These are stationary residuals, not a simulated transition.
    def original_stationary(v):
        P, Y, theta, phi = v
        p = dict(p0, theta=theta, phi=phi)
        T = p["q"]*p["tau"]*P*p["Hbar"]/(2*Y)
        h = young_choices([P,P,P], [T,T], p)
        old = h["pi"]*h["owner"]["z"][5]+(1-h["pi"])*h["renter"]["z"][5]
        return np.array([Y*(h["housing"]+old)-p["Hbar"],
                         Y-p["nu"]*Y*h["fertility"]])

    D = complex_jacobian(original_stationary,
                         [1.,1.,float(pro.theta0),float(pro.phi0)])
    L = [[sum(c for (offset, component), c in pro.stencil[k].items() if component==l)
          for l in range(2)] for k in range(2)]
    A = [[sum(c for s,c in pro.AT[k].items() if s%2==l) for l in range(2)]
         for k in range(2)]
    defect = [[Q(k==l)-sum(A[k][j]*L[j][l] for j in range(2))
               for l in range(2)] for k in range(2)]
    defect_norm = max(sum(abs(x) for x in row) for row in defect)
    assert defect_norm < Q(2,10**9)
    reference = pro.residual_iv(3,Q(0),Q(0),Q(0))
    par = np.array([[mid(reference[k].d[name]) for name in ("theta","phi")]
                    for k in range(2)])
    stationary_derivative_error = max(float(np.max(abs(np.array(L,dtype=float)-D[:,:2]))),
                                     float(np.max(abs(par-D[:,2:]))))
    assert stationary_derivative_error < 1e-12
    response = -np.linalg.solve(D[:,:2],D[:,2:])
    assert response[1,0]>2 and response[1,1]>.8

    report = dict(
        scope="Nine independent household evaluations and stationary derivatives; no equilibrium path simulation.",
        original_anchor_primitives_match=True,
        primitive_name_mapping=primitive_names,
        maximum_policy_difference=maximum_policy_error,
        maximum_tenure_probability_difference=maximum_tenure_error,
        household_cases=rows,
        stationary_original_jacobian=D.tolist(),
        maximum_stationary_derivative_difference=stationary_derivative_error,
        stationary_response_rows=["price","young_cohort"],
        stationary_response_columns=["theta","phi"],
        stationary_response=response.tolist(),
        stationary_total_population_response=(2*response[1]).tolist(),
        exact_summed_tail_A=[[str(x) for x in row] for row in A],
        exact_summed_tail_L=[[str(x) for x in row] for row in L],
        summed_tail_inverse_defect=str(defect_norm),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        pro_verifier_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
        original_helper_sha256=hashlib.sha256((ROOT/"code/model/tools/verify_simplified_olg_local_transition.py").read_bytes()).hexdigest(),
        original_anchor_helper_sha256=hashlib.sha256((ROOT/"code/model/tools/verify_simplified_olg_mixed_transition.py").read_bytes()).hexdigest(),
    )
    return report

if __name__ == "__main__":
    certificate = pro.run_certificate()
    expected = json.loads((OUT / "oracle_transition_math_expected.json").read_text())
    comparison = dict(certificate)
    comparison["verifier_sha256"] = expected["verifier_sha256"]
    assert comparison == expected, "Packaged certificate differs beyond its source hash."
    report = run_original_checks()
    report["exact_certificate"] = certificate
    report["matches_original_certificate_except_packaged_source_hash"] = True
    target = OUT / "oracle_transition_math_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verification": "passed", "maximum_each_shock": certificate["maximum_each_shock"], "receipt": str(target)}, indent=2))
