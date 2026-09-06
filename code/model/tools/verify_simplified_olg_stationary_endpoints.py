"""Exact stationary endpoint signs for the proposed two-stage OLG exercise.

Writes a verification receipt under output/model/simplified_olg_amendments.
Uses original household helpers without changing the model or its parameters.
"""
import json
import sys
from pathlib import Path
import numpy as np
import sympy as s
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'output/model/simplified_olg_amendments'
sys.path.insert(0, str(ROOT / 'code/model/tools'))
from verify_simplified_olg_mixed_transition import anchor
from verify_simplified_olg_local_transition import steady_state
R=s.Rational
q=R(1,2); beta=alpha=omega=R(2,5); gamma=R(3,10); chi=R(3,20); kappa=R(1,2)
theta=R(141,400); tau=R(467,9250); phi=R(4,5)
sigma,L=s.symbols('sigma L', positive=True)
dP,dT,dtheta,dphi=s.symbols('dP dT dtheta dphi')
u=R(9717,18500); pi=R(11,21); a=R(1,4); h2=R(1480,3239); S=R(68104,68019)
dh=5*dphi-dP; du=u*dP; dw=(1+q)*dT
ans={}
for own,x,h,n,rho in [(True,s.S.One,s.S.One,R(3,4),R(42,25)), (False,R(99,74),a,R(9,40),R(39,25))]:
    dhm=dh if own else 0
    A=theta/n**2+alpha*kappa**2/(h-kappa*n)**2
    M=s.Matrix([[rho,chi],[-chi/x**2,A]])
    ans[own]=M.inv()*s.Matrix([dw-h*du-u*dhm-(0 if own else q*a*du), alpha*kappa/(h-kappa*n)**2*dhm+dtheta/n])
xO,nO=ans[True]; xR,nR=ans[False]
dWO=dw-du+(alpha/R(5,8)-u)*dh-beta*gamma*dP
dWR=(dw-(1+q)*a*du)/R(99,74)
dpi=pi*(1-pi)/sigma*(dWO-dWR+L*dtheta)
dhy=pi*dh+(1-a)*dpi
dhold=pi*h2*(xO-dP)+(h2-a)*dpi
dS=dhy+dhold
dn=pi*nO+(1-pi)*nR+R(21,40)*dpi
M,rhs=s.linear_eq_to_matrix([dn,dT-q*tau/2*(S*dP+dS)], [dP,dT])
sol=M.inv()*rhs
sub={dP:s.cancel(sol[0]),dT:s.cancel(sol[1])}
outs={'P':sol[0],'T':sol[1], 'Y':-dS.subs(sub)/S,
      'hYoung':dhy.subs(sub), 'hOld':dhold.subs(sub),
      'hOwner':dh.subs(sub), 'hOldOwner':(h2*(xO-dP)).subs(sub), 'pi':dpi.subs(sub)}
expressions={}
for shock,ex in [('theta',{dtheta:1,dphi:0}),('phi',{dtheta:0,dphi:1})]:
    expressions[shock]={k:s.factor(v.subs(ex)) for k,v in outs.items()}
    signs={'P':1,'Y':1,'hYoung':-1 if shock=='theta' else 1,'hOld':-1,
           'hOwner':-1 if shock=='theta' else 1,'hOldOwner':-1}
    for key,wanted in signs.items():
        num,den=s.fraction(expressions[shock][key])
        assert s.degree(num,L)<=1 and s.degree(den,L)==0
        for endpoint in [1,2]:
            assert all(s.sign(c)==wanted for c in s.Poly(num.subs(L,endpoint),sigma).all_coeffs())
            assert all(c>0 for c in s.Poly(den.subs(L,endpoint),sigma).all_coeffs())
num,den=s.fraction(s.factor(M.det()))
assert all(c<0 for c in s.Poly(num,sigma).all_coeffs())
assert all(c>=0 for c in s.Poly(den,sigma).all_coeffs()) and s.Poly(den,sigma).LC()>0
# Since 2<e<3, 1<log(10/3)<2. The affine-in-L coefficient checks
# therefore certify all listed signs for every sigma>0 exactly.
def summary(ss):
    h=ss['households'];p=ss['parameters'];pi=h['pi'];P=ss['price'];T=ss['transfer'];u=(1-p['q']+p['q']*p['tau'])*P
    residuals=[]; margins=[]
    for own,m in [(True,'owner'),(False,'renter')]:
        x,hs,n,saving,c2,hs2,e=h[m]['z']; assets=h[m]['assets']
        residuals.extend([p['theta']/n-p['chi']/x-p['alpha']*p['kappa']/(hs-p['kappa']*n),
                          c2+p['q']*e+u*hs2-assets-(P*hs if own else 0)-T])
        margins.extend([x,hs-p['kappa']*n,n,saving,c2,hs2,e,p['alpha']*x/(hs-p['kappa']*n)-u])
        if own:
            residuals.append((1-p['phi'])*P*hs-p['b'])
            margins.extend([p['owner_cap']-hs,hs-hs2,e-P*hs2])
        else:
            margins.append(p['gamma']*c2/hs2-u)
    assert min(margins)>0 and max(abs(z) for z in residuals)<1e-10
    return {'P':P,'T':T,'Y':ss['cohort'], 'hYoung':h['housing'],
            'hOld':pi*h['owner']['z'][5]+(1-pi)*h['renter']['z'][5],
            'hOwner':h['owner']['z'][1], 'hOldOwner':h['owner']['z'][5], 'pi':pi}, min(margins), max(abs(z) for z in residuals)
checks=[]
for scale in [0.25,1,4,20]:
    p,_,_=anchor(scale)
    for shock in ['theta','phi']:
        eps=1e-6
        low,ml,rl=summary(steady_state(dict(p,**{shock:p[shock]-eps})))
        high,mh,rh=summary(steady_state(dict(p,**{shock:p[shock]+eps})))
        expected={k:float(v.subs({sigma:scale,L:s.log(R(10,3))})) for k,v in expressions[shock].items()}
        observed={k:(high[k]-low[k])/(2*eps) for k in high}
        error=max(abs(expected[k]-observed[k]) for k in high)
        assert error<1e-6
        checks.append({'sigma':scale,'shock':shock,'max_derivative_error':error,
                       'minimum_original_branch_margin':min(ml,mh),'maximum_original_equation_residual':max(rl,rh),
                       'derivatives':expected})
report={'scope':'Stationary only. Local derivative signs for the exact mixed-tenure family at every finite sigma>0. Smoothness gives sufficiently small finite changes with original strict branches. No uniform reform radius or transition proof.',
        'log_enclosure':'1 < log(10/3) < 2, since 2 < e < 3',
        'determinant':str(s.factor(M.det())),
        'expressions':{a:{b:str(c) for b,c in vv.items()} for a,vv in expressions.items()},
        'original_stationary_checks':checks}
(OUT / 'stationary_endpoint_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'exact_sign_assertions':'passed for all sigma>0','stationary_cases':len(checks)*2,
                  'max_derivative_error':max(c['max_derivative_error'] for c in checks),
                  'minimum_original_branch_margin':min(c['minimum_original_branch_margin'] for c in checks),
                  'maximum_original_equation_residual':max(c['maximum_original_equation_residual'] for c in checks),
                  'sigma4':[c for c in checks if c['sigma']==4],
                  'proof_script':str(Path(__file__).resolve()),
                  'proof_receipt':str(OUT / 'stationary_endpoint_checks.json')},indent=2))
