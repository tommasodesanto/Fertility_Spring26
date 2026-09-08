#!/usr/bin/env python3
"""Check exact identities in the separate utilitarian theory note; no model solve."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'output/model/simplified_olg_amendments/utilitarian_checks.json'
checks = []


def zero(name, expression):
    residual = s.factor(expression)
    assert residual == 0, (name, residual)
    checks.append({'name': name, 'exact_residual': str(residual)})


q, P, Pnext, tau, phi, eps = s.symbols('q P Pnext tau phi epsilon', positive=True)
dx, dn, chi, kappa = s.symbols('Delta_x Delta_n chi kappa', positive=True)
alpha, gamma, omega, beta, m, mu, mf = s.symbols(
    'alpha gamma omega beta m mu m_F', positive=True)
price = (1 + q*tau)*P
p, L = price-q*Pnext, price-phi*P
K = 1+gamma+omega
g = gamma/(K*p)

zero('direct young current budget', -q*Pnext*eps + price*eps-p*eps)
zero('direct old current budget', -price*eps+q*Pnext*eps+p*eps)
zero('direct estate preserved', q*Pnext*eps/q-Pnext*eps)
zero('direct young old resources preserved', -Pnext*eps+Pnext*eps)
zero('direct aggregate bond purchases cancel', -q*Pnext*eps+q*Pnext*eps)

transfer = s.symbols('old_cash_transfer', real=True)
old_dc, old_dh, old_de = transfer/K, gamma*transfer/(K*p), omega*transfer/(K*q)
old_dae = q*(old_de-Pnext*old_dh)
zero('old-only redistribution original cash budget', old_dc+old_dae+price*old_dh-transfer)
zero('old-only redistribution estate equation', old_dae/q+Pnext*old_dh-old_de)
for name, quantity in [('consumption',old_dc),('housing',old_dh),('estates',old_de),('financial saving',old_dae)]:
    zero('balanced old-only redistribution preserves aggregate '+name, quantity+quantity.subs(transfer,-transfer))


for label, dh, dc in [('committed fertility', eps, dx),
                      ('restored fertility', eps+kappa*dn, dx+chi*dn)]:
    da = -phi*P*dh/q
    G, J = dc+L*dh, (p-L)*dh/q
    donor = dh/g
    R = G+q*J-donor
    zero(label+' original young budget', dc+q*da+price*dh-G)
    zero(label+' original mortgage covenant', q*da+phi*P*dh)
    zero(label+' original old resources', da+Pnext*dh+J)
    zero(label+' government saving funds future grant', donor+R-G-q*J)
    zero(label+' current housing clearing', dh-g*donor)
    if label == 'restored fertility':
        zero('restored fertility residual includes children',
             R-(dx+(p-1/g)*eps+(chi+kappa*(p-1/g))*dn))

A, p0, L0 = s.symbols('A p L', positive=True)
g0 = gamma/(K*p0)
Lambda = beta*m/q+mu
marg_h = beta*p0*m/q+mu*L0
zero('direct marginal housing gap',
     marg_h-p0*m-((beta/q-1)*p0*m+mu*L0))
zero('transfer welfare decomposition',
     Lambda*A+marg_h-m/g0-mf*(A+p0-1/g0)
     -((Lambda-m)*A+(marg_h-p0*m)+(m-mf)*(A+p0-1/g0)))
R0 = dx+(p0-1/g0)*eps
D0 = eps/g0
zero('old consumption cuts fund present-value resources',
     D0/K+R0/(1+omega)-dx/(1+omega))
zero('consumption plus discounted estates conserved',
     dx-D0/K-R0/(1+omega)-omega*(D0/K+R0/(1+omega)))

x, space, theta, n, Gamma, Delta = s.symbols(
    'x space theta n Gamma Delta', positive=True)
curvA = L0**2/x**2 + alpha/space**2+Gamma*Delta**2
curvB = -L0*chi/x**2 + alpha*kappa/space**2
curvD = chi**2/x**2 + alpha*kappa**2/space**2+theta/n**2
zero('gift housing numerator is strictly positive',
     L0*curvD+chi*curvB
     -(L0*theta/n**2+alpha*kappa*(chi+L0*kappa)/space**2))
zero('gift fertility numerator is strictly positive',
     L0*curvB+chi*curvA
     -(alpha*(chi+L0*kappa)/space**2+chi*Gamma*Delta**2))
f = theta/n-chi/(s.Symbol('c', positive=True)-chi*n)-alpha*kappa/(s.Symbol('h', positive=True)-kappa*n)
c, h = s.symbols('c h', positive=True)
zero('fertility first-order condition strictly decreases in n',
     s.diff(f,n)+theta/n**2+chi**2/(c-chi*n)**2+alpha*kappa**2/(h-kappa*n)**2)

special = {gamma:alpha*(1+omega)}
zero('simple transfer case equal housing responses',
     (g0-alpha/(p0*(1+alpha))).subs(special))
zero('free fertility changes simple-case residual sign',
     (p0/alpha+(p0-1/g0)+(chi+kappa*(p0-1/g0))*dn
      -(chi-kappa*p0/alpha)*dn).subs(special))

# Simplifications identified in the completed Pro review and rederived locally.
Hcap = s.symbols('Hcap', positive=True)
zcap = K*p0*Hcap/gamma
zero('old capped and uncapped cash values meet at cap threshold',
     (1+omega)/(zcap-p0*Hcap)-gamma/(p0*Hcap))
zero('old uncapped cash value at cap threshold', K/zcap-gamma/(p0*Hcap))
ah_free = (alpha/p0+theta*kappa/(chi+p0*kappa))/(1+alpha+theta)
zero('endogenous fertility changes the simple balanced housing response',
     ah_free-alpha/(p0*(1+alpha))
     -theta*(kappa*p0-alpha*chi)/(p0*(chi+kappa*p0)*(1+alpha)*(1+alpha+theta)))
x_new, s_new = s.symbols('x_new s_new', positive=True)
zero('original grant pair profiled fertility derivative',
     chi/x+alpha*kappa/space-chi/x_new-alpha*kappa/s_new
     -(chi*(1/x-1/x_new)+alpha*kappa*(1/space-1/s_new)))

# Independent verification of the full-choice phi=q appendix.
r, t, G, taxQ = s.symbols('r t G Q', positive=True)
E, age_a = 1+alpha+theta, alpha+theta
f = t+t*(1-t)/(theta-age_a*t)
zero('capped fertility inverse solves original condition', theta/t-1/(f-t)-alpha/(1-t))
fp = E/age_a+alpha*theta/(age_a*(theta-age_a*t)**2)
zero('capped fertility inverse first derivative', s.diff(f,t)-fp)
zero('capped fertility inverse second derivative', s.diff(f,t,2)-2*alpha*theta/(theta-age_a*t)**3)
pa_h = (alpha+theta/(r+1))/E
an_chi = theta*r/(E*(r+1))
t_star = theta/(alpha*(r+1)+theta)
q_star = an_chi*fp.subs(t,t_star)
zero('cap threshold tax ratio', q_star-(1-pa_h*(1-1/(alpha*r))))
ah = pa_h/p0
d = ah/g0
zero('cap threshold surplus is strictly positive', q_star+d-1-ah*(1/g0-p0+p0/(alpha*r)))
sqrt_arg = alpha*r/(E*(alpha*r+age_a))
fp_at_upper = E/age_a+alpha/(age_a*theta*sqrt_arg)
zero('explicit upper cap endpoint has unit tax ratio', an_chi*fp_at_upper-1)
rebateS = taxQ+d*G-G
zero('full-choice current government budget', G+rebateS-taxQ-d*G)
free_C = (1-p0*ah)*G-taxQ-d*G/K+rebateS/(1+omega)
free_qE = -omega*d*G/K+omega*rebateS/(1+omega)
zero('full-choice consolidated goods and estate account', free_C+free_qE)
nu, nbar = s.symbols('nu nbar', positive=True)
scale = nu*nbar
zero('common child-cost scaling delivers replacement', nbar/scale-1/nu)
zero('common child-cost scaling preserves child goods', chi*scale*n/scale-chi*n)
zero('common child-cost scaling preserves child housing', kappa*scale*n/scale-kappa*n)

note = ROOT/'latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex'
old = ROOT/'latex/JMP_DS_suggestions/simplified_olg_conventional_finance.tex'
text, prior = note.read_text(), old.read_text()
new_model = text.split(r'\section{Environment}',1)[1].split(r'\section{Equilibrium and the welfare comparison}',1)[0]
old_model = prior.split(r'\section{Environment}',1)[1].split(r'\section{Equilibrium and financing}',1)[0]
new_model = new_model.replace(r'\Needspace{10\baselineskip}'+chr(10), '')
assert new_model == old_model, 'Author-established environment or household notation changed'
checks.append({'name':'environment and all four household problems preserved verbatim', 'pass':True})

receipt = {
    'generated_utc':datetime.now(timezone.utc).isoformat(),
    'method':'Exact symbolic identities plus source-preservation comparison; no numerical equilibrium or parameter-neighborhood proof.',
    'passed':True, 'check_count':len(checks), 'checks':checks,
    'sympy_version':s.__version__,
    'source_sha256':sha256(note.read_bytes()).hexdigest(),
    'driver_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
    'analytical_review_files':[
        'utilitarian_direct_review.md','utilitarian_transfers_review.md',
        'utilitarian_fertility_path_review.md','utilitarian_transfers_hostile_review.md',
        'utilitarian_fertility_hostile_review.md',
        'utilitarian_assembled_transfer_review.md','utilitarian_assembled_fertility_review.md',
        'utilitarian_old_redistribution_review.md','utilitarian_free_choice_transfer_review.md',
        'utilitarian_free_choice_hostile_review.md'],
    'scope_limits':[
        'Identities support the analytical proofs; they are not a substitute for existence or inequality proofs.',
        'Main young-targeted transfer construction assumes committed fertility and tenure.',
        'The full-choice appendix requires phi=q and taste-informed individual targeting; aggregate fertility remains fixed.',
        'No all-date policy fertility sign, endogenous-population welfare ranking, or convergence theorem is certified.',
        'Pro review was retrieved and its capped-funder and exact-grant fertility simplifications independently checked; it predates the full-choice appendix.']}
OUT.write_text(json.dumps(receipt,indent=2)+'\n')
print(f'{len(checks)} exact checks passed; receipt: {OUT}')
