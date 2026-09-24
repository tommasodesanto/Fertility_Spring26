"""Exact finite bookkeeping examples. Does not import or solve the model."""
from fractions import Fraction as F
from itertools import product
from math import comb
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
f, mu, delta = F(21, 10), F(2, 9), F('0.602359422009')
c = 1 / f
checks = []

def record(name, value):
    checks.append({'name': name, 'exact': str(value), 'decimal': float(value)})

def ledger(name, H, C, B, X, M, R, rho=F(1), outside=F(0), JH=F(0), JC=F(0)):
    departures = M + R
    E = c * rho * departures + outside
    Cn = C + B - departures + JC
    Hn = H - X + E + JH
    assert Cn >= 0 and Hn >= 0
    assert f * Hn + Cn - (f * H + C) == B - f*X - (1-rho)*departures + f*outside + f*JH + JC
    assert 2 * Hn + Cn - (2 * H + C) == B - 2*X - (1-rho)*departures + (2*c-1)*rho*departures + 2*outside + 2*JH + JC
    for k,v in [('H_next', Hn), ('C_next', Cn), ('entry', E)]:
        record(name + '.' + k, v)
    return Hn, Cn

ledger('no_death_births', F(2), F(0), F(2), F(0), 2*mu, F(0))
ledger('ordinary_realized', F(1), F(1), F(0), F(0), F(1), F(0))
ledger('death_release', F(1), F(2), F(0), F(1), F(0), F(2))
ledger('half_survive', F(1), F(2), F(0), F(1,2), mu, F(1))
ledger('retention_outside_bridge', F(3), F(4), F(1), F(1,2), F(1), F(1), F(4,5), F(1,10), F(1,20), F(1,7))
record('ten_ordinary_two_death_entrants', c*12)
assert c*10 + c*2 == c*12
assert 2*c-1 == -F(1,21)
assert c*f == 1
assert F(1,2)*f == F(21,20)
assert F(20,21)*F(1,2)*f == 1

# Enumerate every subset of departing represented children. The last slot is
# the tagged third child, when tag_present=True. This is an independent
# enumeration check of compressed first-moment propagation z'=K*m'/m*z.
enumerated_cases = 0
for m in range(1,4):
    for tag_present in (False, True):
        for survival in (F(0), F(1,3), F(1)):
            g_to = [F(0)]*(m+1)
            z_to = [F(0)]*(m+1)
            M = F(0)
            for stays in product((0,1), repeat=m):
                p = survival
                for stays_i in stays:
                    p *= (1-mu) if stays_i else mu
                mn = sum(stays)
                zn = delta if tag_present and stays[-1] else F(0)
                g_to[mn] += p
                z_to[mn] += p*zn
                M += p*((m-mn) + (delta if tag_present and not stays[-1] else 0))
            z = delta if tag_present else F(0)
            for mn in range(m+1):
                K = F(comb(m,mn))*(1-mu)**mn*mu**(m-mn)
                assert g_to[mn] == survival*K
                assert z_to[mn] == survival*K*F(mn,m)*z
            R = (1-survival)*(m+z)
            weighted_next = sum(F(mn)*g_to[mn]+z_to[mn] for mn in range(m+1))
            assert m+z == M+R+weighted_next
            assert M == survival*mu*(m+z)
            assert z_to[0] == 0
            enumerated_cases += 1

# Mixing cells does not require a hidden binary decision state: the first
# moment is linear. Two different tagged-status mixtures are transported by
# the same state-dependent choice probability and still obey the formula.
for m in range(1,4):
    for prob_tag in (F(0), F(1,4), F(2,3), F(1)):
        for choice in (F(1,7), F(3,5), F(1)):
            for mn in range(m+1):
                K = F(comb(m,mn))*(1-mu)**mn*mu**(m-mn)
                compressed = choice*K*F(mn,m)*delta*prob_tag
                detailed = choice*prob_tag*K*F(mn,m)*delta + choice*(1-prob_tag)*K*0
                assert compressed == detailed

# Third birth creates only its incremental weight. Earlier departed children
# cannot be retroactively reweighted. Example: two prior births, one remains.
C_before = F(1)
C_after_birth = F(2)+delta
assert C_after_birth-C_before == 1+delta
record('third_birth_weight', 1+delta)
record('third_birth_two_slots_M', mu*(2+delta))
record('third_birth_two_slots_C_next', (1-mu)*(2+delta))
record('third_birth_all_three_M', mu*(3+delta))
record('third_birth_all_three_C_next', (1-mu)*(3+delta))
# Conditional two slots -> one: each slot equally likely, expected tag half.
assert F(1)+delta/2 == (F(1)+F(1)+delta)/2
record('third_birth_two_to_one_conditional_C', 1+delta/2)

# Parent death does not pay an estate a second time through a new household.
estate, recipients, entry_wealth_mean = F(150), F(3), F(7)
transfer = estate / recipients
assert transfer*recipients == estate
record('estate_pool_transfer', transfer)
record('entry_asset_boundary_inflow', c*2*entry_wealth_mean)

# With certain terminal release, lifetime departure probabilities add to one;
# the death and ordinary routes partition each date exactly.
survivals = [F(1), F(1), F(4,5), F(1,2), F(0)]
pending, normal, death = F(1), F(0), F(0)
for s in survivals:
    normal += pending*s*mu
    death += pending*(1-s)
    pending *= s*(1-mu)
assert normal+death == 1 and pending == 0
record('finite_cohort_ordinary_share', normal)
record('finite_cohort_death_share', death)

out = {
    'status': 'passed',
    'scope': 'Exact rational toy bookkeeping and exhaustive subset enumeration; no model imports, model solves, checkpoint loading or equilibrium tests.',
    'enumerated_tagged_transition_cases': enumerated_cases,
    'mixing_cases': 36,
    'checks': checks,
}
(ROOT/'hand_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'], 'enumerated_cases':enumerated_cases, 'scalar_results':len(checks)}))
