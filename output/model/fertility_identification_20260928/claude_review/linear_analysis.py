#!/usr/bin/env python3
"""Read-only linear analysis of the September 28 current-point Jacobian.

Standard library only (no numpy, no model imports). Inputs are the saved
Jacobian CSV, the scaled-SVD JSON, the overnight block0506 primary rescore
table and the NCHS first-birth counts file. Writes linear_analysis_outputs.md
next to this script. Run from anywhere:

    python3 output/model/fertility_identification_20260928/claude_review/linear_analysis.py
"""
import csv, json, math, os, collections, datetime

ROOT = '/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26'
ID = os.path.join(ROOT, 'output/model/fertility_identification_20260928')
JAC = os.path.join(ID, 'lead_review/jacobian.csv')
SVD = os.path.join(ID, 'lead_review/jacobian_scaled_svd.json')
FIT = os.path.join(ROOT, 'output/model/overnight_calibration_20260928/cluster/final_review/target_fit_primary_rescore.csv')
NCHS = os.path.join(ROOT, 'code/data/nchs_natality_timing/first_birth_counts_year_age.csv')
OUT = os.path.join(ID, 'claude_review/linear_analysis_outputs.md')

rows = list(csv.DictReader(open(JAC)))
svd = json.load(open(SVD))
params = svd['parameters']; moms = svd['moments']; scales = svd['column_scales']
ALLM = ['initial_normalization','cps_childlessness','cps_exactly_one','nchs_mean_age','nchs_share30',
        'wealth_earnings','bequest_wealth','old_dispersion','mean_rooms','ownership_30_55',
        'first_birth_rooms','family_rooms','recent_parent_ownership','early_fertility','psi_child']
fit = {r['moment']: r for r in csv.DictReader(open(FIT))}
w0 = {m: float(fit[m]['weight']) for m in moms}
gap = {m: float(fit[m]['gap']) for m in moms}
tgt = {m: float(fit[m]['target']) for m in moms}
BOUNDS = {'H0':(0.2,80),'beta_annual':(0.94,0.99),'chi':(0.1,5),'first_birth_fixed_cost':(0,8),
          'kappa_fert':(0.02,50),'kappa_fert_continuation':(0.02,50),'theta0':(0,8),
          'delta_alpha_jump':(0,0.25),'child_benefit_curvature':(0,0.8),'tenure_choice_kappa':(0.001,0.1)}
SHORT = {'initial_normalization':'norm2.1','cps_childlessness':'childless','cps_exactly_one':'one',
         'nchs_mean_age':'age1','nchs_share30':'share30','wealth_earnings':'W/E','bequest_wealth':'beq',
         'old_dispersion':'oldp90','mean_rooms':'rooms','ownership_30_55':'own','first_birth_rooms':'fbrooms',
         'family_rooms':'famrooms','recent_parent_ownership':'recent','early_fertility':'early','psi_child':'psi'}
PSHORT = {'H0':'H0','beta_annual':'beta','chi':'chi','first_birth_fixed_cost':'fbcost','kappa_fert':'k_first',
          'kappa_fert_continuation':'k_cont','theta0':'theta0','delta_alpha_jump':'dalpha',
          'child_benefit_curvature':'curv','tenure_choice_kappa':'k_ten'}

def J(field):
    d = {(r['parameter'], r['moment']): float(r[field]) for r in rows}
    return {m: [d[(p, m)] * scales[p] for p in params] for m in ALLM}   # per log-unit of parameter
def Jraw(field):
    d = {(r['parameter'], r['moment']): float(r[field]) for r in rows}
    return {m: [d[(p, m)] for p in params] for m in ALLM}

def solve(A, b):
    n = len(A); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c])); M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1): M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def norm(v): return math.sqrt(dot(v, v))
def val(p, x): return scales[p] * math.exp(max(min(x, 50), -50))

L = []
def w(s=''): L.append(s)
w('# Linear analysis outputs (generated %s EDT)' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
w(); w('Inputs: `%s`, `%s`, `%s`, `%s`.' % (JAC, SVD, FIT, NCHS))
w('All derivatives are the saved central differences at the overnight block0506 anchor. "Per log-unit" means the saved derivative multiplied by the anchor parameter value (the contract SVD convention), so entries read as the moment change for a 100 log-percent parameter change; divide by 100 for a 1 percent change.')

# 1. pivot tables
for field, label in [('full_derivative', 'full step'), ('half_derivative', 'half step')]:
    Jl = J(field)
    w(); w('## Jacobian per log-unit parameter (%s)' % label); w()
    w('| parameter | ' + ' | '.join(SHORT[m] for m in ALLM) + ' |'); w('|' + '---|' * (len(ALLM) + 1))
    for i, p in enumerate(params):
        w('| %s | ' % PSHORT[p] + ' | '.join('%+.3g' % Jl[m][i] for m in ALLM) + ' |')
Jf, Jh = J('full_derivative'), J('half_derivative')
w(); w('## Step sensitivity: |full - half| / max(|full|, |half|)'); w()
w('| parameter | ' + ' | '.join(SHORT[m] for m in ALLM) + ' |'); w('|' + '---|' * (len(ALLM) + 1))
for i, p in enumerate(params):
    w('| %s | ' % PSHORT[p] + ' | '.join('%.2f' % (abs(Jf[m][i] - Jh[m][i]) / max(abs(Jf[m][i]), abs(Jh[m][i]), 1e-12)) for m in ALLM) + ' |')
w(); w('Entries above 0.2 mark derivatives that are not stable across the two step sizes (normalization row, old-dispersion row, several bequest entries).')

# 2. implied weights
w(); w('## Primary weights, implied scales and anchor loss contributions'); w()
w('| moment | target | model (block0506) | gap | weight w | 1/sqrt(w) | w*gap^2 |'); w('|---|---:|---:|---:|---:|---:|---:|')
for m in moms:
    w('| %s | %.6g | %.6g | %+.5g | %.6g | %.5g | %.4g |' % (m, tgt[m], tgt[m] + gap[m], gap[m], w0[m], 1 / math.sqrt(w0[m]), w0[m] * gap[m] ** 2))
w('| total | | | | | | %.4f |' % sum(w0[m] * gap[m] ** 2 for m in moms))
w('| total excluding early fertility | | | | | | %.4f |' % (sum(w0[m] * gap[m] ** 2 for m in moms) - w0['early_fertility'] * gap['early_fertility'] ** 2))

# 3. constrained directions
def constrained(Jl, cons):
    g = Jl['early_fertility']; A = [Jl[m] for m in cons]
    AAt = [[dot(a, b) for b in A] for a in A]; lam = solve(AAt, [dot(a, g) for a in A])
    d = [g[j] - sum(lam[i] * A[i][j] for i in range(len(A))) for j in range(len(g))]
    return g, d
w(); w('## Constrained early-fertility directions (which parameter moves raise early fertility while holding named moments fixed, to first order)'); w()
for field, Jl in [('full step', Jf), ('half step', Jh)]:
    for name, cons in [('hold mean first-birth age', ['nchs_mean_age']),
                       ('hold mean age and childlessness', ['nchs_mean_age', 'cps_childlessness']),
                       ('hold mean age, childlessness, exactly-one', ['nchs_mean_age', 'cps_childlessness', 'cps_exactly_one'])]:
        g, d = constrained(Jl, cons)
        gain = dot(g, d) / norm(d); step = 0.10 / dot(g, d); ds = [x * step for x in d]
        w('### %s, %s' % (field, name)); w()
        w('Early-fertility gain per unit log-parameter step along the best feasible direction: **%+.4f** (unconstrained gradient norm %.4f). Log-parameter step for +0.10 early fertility has norm %.2f:' % (gain, norm(g), norm(ds)))
        w(); w('| parameter | log step | multiplier |'); w('|---|---:|---:|')
        for p, x in zip(params, ds): w('| %s | %+.3f | %.3f |' % (p, x, math.exp(x)))
        w(); w('Implied first-order change in every reported moment: ' + ', '.join('%s %+.3g' % (SHORT[m], dot(Jl[m], ds)) for m in ALLM))
        w()

# 4. Gauss-Newton with ridge
w('## Linearized Gauss-Newton steps from block0506 (ridge on log steps as a trust region)'); w()
w('| lane | ridge | lane loss at anchor | linear-predicted lane loss | common-primary rescore | other-moment primary loss | max abs log step | predicted early | predicted mean age | predicted W/E | predicted ownership |'); w('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
gn_points = {}
for mult, lane in [(1, 'primary'), (10, 'early10'), (100, 'early100')]:
    for ridge in [1000.0, 100.0, 10.0]:
        wl = dict(w0); wl['early_fertility'] *= mult; n = len(params)
        A = [[sum(wl[m] * Jf[m][i] * Jf[m][j] for m in moms) + (ridge if i == j else 0) for j in range(n)] for i in range(n)]
        b = [-sum(wl[m] * Jf[m][i] * gap[m] for m in moms) for i in range(n)]
        d = solve(A, b); ng = {m: gap[m] + dot(Jf[m], d) for m in moms}
        L0 = sum(wl[m] * gap[m] ** 2 for m in moms); L1 = sum(wl[m] * ng[m] ** 2 for m in moms); L1p = sum(w0[m] * ng[m] ** 2 for m in moms)
        w('| %s | %g | %.2f | %.2f | %.2f | %.2f | %.2f | %.3f | %.3f | %.2f | %.3f |' % (lane, ridge, L0, L1, L1p, L1p - w0['early_fertility'] * ng['early_fertility'] ** 2, max(abs(x) for x in d), tgt['early_fertility'] + ng['early_fertility'], tgt['nchs_mean_age'] + ng['nchs_mean_age'], tgt['wealth_earnings'] + ng['wealth_earnings'], tgt['ownership_30_55'] + ng['ownership_30_55']))
        gn_points[(lane, ridge)] = {p: val(p, x) for p, x in zip(params, d)}
w(); w('Predicted parameter points (ridge 10):'); w()
w('| parameter | anchor | primary GN | early10 GN | early100 GN | bounds |'); w('|---|---:|---:|---:|---:|---|')
for p in params:
    w('| %s | %.5g | %.5g | %.5g | %.5g | [%g, %g] |' % (p, scales[p], gn_points[('primary', 10.0)][p], gn_points[('early10', 10.0)][p], gn_points[('early100', 10.0)][p], BOUNDS[p][0], BOUNDS[p][1]))

# 5. linear frontier
w(); w('## Linear trade-off frontier: minimum other-moment primary loss for a given early-fertility gain'); w()
w('Minimizes the primary loss of the nine other scored moments (plus ridge 10 on log steps) subject to a first-order early-fertility gain delta. Anchor other-moment loss 12.07.'); w()
w('| early gain | early level | min other-moment loss | mean age | childless | exactly-one | W/E | ownership | rooms | k_first | k_cont | fbcost | beta |'); w('|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
others = [m for m in moms if m != 'early_fertility']; n = len(params)
for delta in [0.03, 0.05, 0.067, 0.10, 0.132, 0.15, 0.20, 0.274]:
    A = [[sum(w0[m] * Jf[m][i] * Jf[m][j] for m in others) + (10.0 if i == j else 0) for j in range(n)] + [Jf['early_fertility'][i]] for i in range(n)]
    A.append(Jf['early_fertility'] + [0.0])
    b = [-sum(w0[m] * Jf[m][i] * gap[m] for m in others) for i in range(n)] + [delta]
    d = solve(A, b)[:n]; ng = {m: gap[m] + dot(Jf[m], d) for m in moms}
    Lo = sum(w0[m] * ng[m] ** 2 for m in others); pt = {p: val(p, x) for p, x in zip(params, d)}
    w('| %+.3f | %.3f | %.1f | %.2f | %.3f | %.3f | %.2f | %.3f | %.2f | %.3f | %.3f | %.3f | %.4f |' % (delta, tgt['early_fertility'] + ng['early_fertility'], Lo, tgt['nchs_mean_age'] + ng['nchs_mean_age'], tgt['cps_childlessness'] + ng['cps_childlessness'], tgt['cps_exactly_one'] + ng['cps_exactly_one'], tgt['wealth_earnings'] + ng['wealth_earnings'], tgt['ownership_30_55'] + ng['ownership_30_55'], tgt['mean_rooms'] + ng['mean_rooms'], pt['kappa_fert'], pt['kappa_fert_continuation'], pt['first_birth_fixed_cost'], pt['beta_annual']))
w(); w('Comparison points from actual solves: early_frontier half first-birth scale (+0.067, other-moment loss 491.5), quarter scale (+0.132, other-moment loss 1398.4); overnight identity winner 0526 (+0.121, other-moment loss 550.5 = 552.79 - 100*0.1528^2).')

# 6. NCHS cell shares and bound
nrows = list(csv.DictReader(open(NCHS)))
byage = collections.Counter(); N = 0.0
for r in nrows:
    if 2003 <= int(r['year']) <= 2006:
        byage[int(r['age'])] += float(r['n_first_births']); N += float(r['n_first_births'])
def share(lo, hi): return sum(v for a, v in byage.items() if lo <= a <= hi) / N
def mid(a):
    if a <= 21: return 20
    if a >= 42: return 44
    return 24 + 4 * ((a - 22) // 4)
cells = [('<=21 (model cell 18-21, includes ages 12-17)', share(0, 21)), ('22-25', share(22, 25)), ('26-29', share(26, 29)), ('30-33', share(30, 33)), ('34-37', share(34, 37)), ('38-41', share(38, 41)), ('42+', share(42, 99))]
w(); w('## NCHS 2003-2006 pooled first births: cell shares and mean-age conventions'); w()
w('| model cell | share of first births |'); w('|---|---:|')
for c, s in cells: w('| %s | %.4f |' % (c, s))
w('| memo: age < 18 | %.4f |' % share(0, 17)); w('| memo: age <= 19 | %.4f |' % share(0, 19)); w('| memo: age >= 30 | %.4f |' % share(30, 99))
w(); w('Mean first-birth age, raw single-year: %.3f. Mean of cell midpoints (contract rule, the actual target): %.3f. Raw mean among mothers 18+: %.3f.' % (sum(a * v for a, v in byage.items()) / N, sum(mid(a) * v for a, v in byage.items()) / N, sum(a * v for a, v in byage.items() if a >= 18) / sum(v for a, v in byage.items() if a >= 18)))
M = 1 - tgt['cps_childlessness']; f1 = share(0, 21); f2 = share(22, 25)
w(); w('## Time-aggregation bound on model early fertility under data-consistent first-birth timing'); w()
w('Model children ever born at [25,26) = P(birth in cell 18-21) + 0.875 * P(birth in cell 22-25), with at most one birth per four-year cell. With eventual-mother share M = 1 - %.4f = %.4f and NCHS first-birth cell shares f1 = %.4f, f2 = %.4f, and p2 the probability that a cell-1 mother has a second birth in cell 2:' % (tgt['cps_childlessness'], M, f1, f2))
w(); w('    early = M*f1 + 0.875*(M*f2 + p2*M*f1) = %.4f + %.4f*p2' % (M * f1 + 0.875 * M * f2, 0.875 * M * f1)); w()
w('| p2 | model-consistent early fertility |'); w('|---:|---:|')
for p2 in [0.0, 0.4, 0.5, 0.6, 0.8, 1.0]: w('| %.1f | %.3f |' % (p2, M * f1 + 0.875 * M * f2 + 0.875 * M * f1 * p2))
w(); w('Target 0.8095; block0506 model 0.5354 (implied p2 about %.2f under data-consistent timing). The p2 = 1 ceiling is %.3f.' % ((0.5354 - (M * f1 + 0.875 * M * f2)) / (0.875 * M * f1), M * f1 + 0.875 * M * f2 + 0.875 * M * f1))

# 7. near-bound rule
w(); w('## Near-bound flag arithmetic'); w()
for p in ['kappa_fert', 'kappa_fert_continuation']:
    lo, hi = BOUNDS[p]; est = scales[p]
    w('- %s = %.4f, bounds [%g, %g]: raw-range rule threshold 0.01*(hi-lo) = %.3f; distance to lower bound %.3f -> flagged %s; in log units the distance to the lower bound is %.2f of a %.2f log range (%.0f percent).' % (p, est, lo, hi, 0.01 * (hi - lo), est - lo, est - lo <= 0.01 * (hi - lo), math.log(est / lo), math.log(hi / lo), 100 * math.log(est / lo) / math.log(hi / lo)))

open(OUT, 'w').write('\n'.join(L) + '\n')
print('wrote', OUT, len(L), 'lines')
