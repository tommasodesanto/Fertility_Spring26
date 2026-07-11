#!/usr/bin/env python3
"""Decisive numeric checks for the bequest-utility audit (spec audit B).

Checks, against intergen_housing_fertility.solver.bequest_utility_vec:
  1. Normalized CRRA form matches the paper display Eq. (quant_bequest)
     term by term at production sigma=2, theta1=0.01.
  2. B(0,n)=0 exactly under normalize_bequest_utility=True, for all n.
  3. Weak monotonicity in W including the max(W,0) clip region
     (negative gross estates map to B(0,n)).
  4. Family-size scaling uses the encoded index n (marginal ratio =
     (1+theta_n*n)/(1+theta_n*0)).
  5. Old (unnormalized) form level gap: old = new - theta0*(1+theta_n n)/theta1
     at sigma=2 (level shift 100*scale at theta1=0.01).
  6. get_completed_fertility is invariant under deterministic child aging
     (bequest scaling unaffected by apply_child_aging at terminal).
"""
from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility.solver import (
    bequest_utility_vec,
    get_completed_fertility,
)

P = SimpleNamespace(
    theta0=0.221498,
    theta_n=1.0711,
    theta1=0.01,
    sigma=2.0,
    normalize_bequest_utility=True,
    estate_tax_rate=0.0,
    estate_tax_exemption=0.0,
    n_child_stages=4,
)

W = np.array([-2.0, -0.5, 0.0, 0.01, 0.1, 0.5, 1.0, 5.0, 20.0])

ok = True
for n in (0, 1, 2):
    got = bequest_utility_vec(W, n, P)
    scale = P.theta0 * (1.0 + P.theta_n * n)
    Wc = np.maximum(W, 0.0)
    want = scale * ((P.theta1 + Wc) ** (1.0 - P.sigma) - P.theta1 ** (1.0 - P.sigma)) / (1.0 - P.sigma)
    if not np.allclose(got, want, rtol=0, atol=1e-13):
        ok = False
        print(f"FORMULA MISMATCH n={n}: {got - want}")
    if abs(got[2]) > 1e-14 or abs(got[0]) > 1e-14:
        ok = False
        print(f"B(0,{n}) != 0: {got[2]}, B(neg,{n})={got[0]}")
    if np.any(np.diff(got) < -1e-14):
        ok = False
        print(f"NON-MONOTONE n={n}: {got}")
print("1-3. normalized formula / B(0,n)=0 / monotone incl clip:", "PASS" if ok else "FAIL")

# 4. marginal scaling ratio
eps = 1e-6
m0 = np.diff(bequest_utility_vec(np.array([0.4, 0.4 + eps]), 0, P))[0]
m2 = np.diff(bequest_utility_vec(np.array([0.4, 0.4 + eps]), 2, P))[0]
ratio = m2 / m0
print(f"4. marginal ratio n=2 vs n=0: {ratio:.6f} (expected {1 + 2 * P.theta_n:.6f})",
      "PASS" if abs(ratio - (1 + 2 * P.theta_n)) < 1e-4 else "FAIL")

# 5. old-form level gap
P_old = SimpleNamespace(**{**P.__dict__, "normalize_bequest_utility": False})
for n in (0, 2):
    new = bequest_utility_vec(np.array([0.5]), n, P)[0]
    old = bequest_utility_vec(np.array([0.5]), n, P_old)[0]
    scale = P.theta0 * (1.0 + P.theta_n * n)
    shift = scale / P.theta1  # theta1^{1-sigma}/(sigma-1) at sigma=2
    print(f"5. n={n}: new={new:.4f} old={old:.4f} gap={new - old:.4f} expected_shift={shift:.4f}",
          "PASS" if abs((new - old) - shift) < 1e-10 else "FAIL")

# 6. completed fertility invariant under deterministic aging map
K = P.n_child_stages
def age(nn, cs):
    if cs == 0:
        return 0
    if cs >= K + 1:
        return cs
    if cs < K:
        return cs + 1
    return 0 if nn == 0 else (K + 1 if nn == 1 else K + 2)

bad = []
for nn in range(3):
    for cs in range(K + 3):
        if get_completed_fertility(nn, cs, P) != get_completed_fertility(nn, age(nn, cs), P):
            # (nn=0, cs in 1..K) is unreachable; flag only reachable states
            if not (nn == 0 and 1 <= cs <= K):
                bad.append((nn, cs))
print("6. nk invariant under child aging (reachable states):", "PASS" if not bad else f"FAIL {bad}")
