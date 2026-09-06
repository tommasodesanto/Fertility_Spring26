#!/usr/bin/env python3
"""Read-only theta-shock certificate. No output files or repository mutations."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/model/tools'))
import verify_simplified_olg_mixed_transition as m
import numpy as np
import sympy as s
r=s.Rational
sigma=r(4); base=m.exact_linearization(sigma)
J,B=base['J'],base['B'];q=r(1,2);beta=alpha=omega=r(2,5);gamma=r(3,10);chi=r(3,20);kappa=r(1,2)
theta=r(141,400);tau=r(467,9250);H=r(68104,68019);pi=r(11,21);a=r(1,4);u=r(9717,18500)
xO=hO=s.S.One;nO=r(3,4);xR=r(99,74);nR=r(9,40)
rhoO=1+beta*(1+gamma+omega);rhoR=1+beta*(1+omega)
L=s.symbols('L',real=True)
def partials(x,h,n,rho,resource,dtheta=0):
    A=theta/n**2+alpha*kappa**2/(h-kappa*n)**2
    return s.Matrix([[rho,chi],[-chi/x**2,A]]).inv()*s.Matrix([resource,dtheta/n])
xOt,nOt=partials(xO,hO,nO,rhoO,0,1)
xRt,nRt=partials(xR,a,nR,rhoR,0,1)
# Envelope theorem: d(WO-WR)/dtheta = log(nO/nR) at fixed prices/transfers.
pit=pi*(1-pi)/sigma*L
Ft=s.Matrix([(hO-a)*pit,-2*(pi*nOt+(1-pi)*nRt+(nO-nR)*pit)])
Gt=s.Matrix([0,0,0,0,beta*gamma/q*(xO*pit+pi*xOt),-a*pit])
Gv=s.zeros(6,2); Fv=s.zeros(2,2)
for k in range(2):
    dun=s.S.One if k==0 else s.S.Zero
    dYn=s.S.One if k==1 else s.S.Zero
    dw=-q*q*tau*H/4*dYn
    dxO,dnO=partials(xO,hO,nO,rhoO,dw)
    dxR,dnR=partials(xR,a,nR,rhoR,dw-q*a*dun)
    dpi=pi*(1-pi)/sigma*(dw/xO-beta*gamma/u*dun-(dw-q*a*dun)/xR)
    Gv[:,k]=s.Matrix([0,dun,dYn,0,beta*gamma/q*(xO*dpi+pi*dxO),-a*dpi])
    Fv[:,k]=s.Matrix([(hO-a)*dpi,dYn-2*(pi*dnO+(1-pi)*dnR+(nO-nR)*dpi)])
assert Fv==base['implicit']
D=(Gt-Gv*Fv.inv()*Ft).applyfunc(s.cancel)
stat=((s.eye(6)-J).inv()*D).applyfunc(s.cancel)
assert all(s.simplify(e)==0 for e in ((s.eye(6)-J)*stat-D))
# log(10/3) = 2 atanh(7/13): rational partial sum and positive tail bound.
x=r(7,13); terms=70
loglo=2*sum(x**(2*k+1)/r(2*k+1) for k in range(terms))
loghi=loglo+2*x**(2*terms+1)/(r(2*terms+1)*(1-x*x))
logiv=m.I(loglo,loghi)
def as_interval(expr):
    poly=s.Poly(expr,L); assert poly.degree()<=1
    return m.I(poly.coeff_monomial(1))+m.I(poly.coeff_monomial(L))*logiv
Di=[as_interval(e) for e in D];sti=[as_interval(e) for e in stat]
p=base['poly'];z=p.gens[0];five=s.Poly(p.as_expr()/z,z)
assert five.count_roots(-s.oo,s.oo)==3
isolated=five.intervals(eps=r(1,10**35))
assert len(isolated)==3 and all(mult==1 for _,mult in isolated)
neg,small,pos=[m.I(*bounds) for bounds,mult in isolated]
assert neg.hi<-1 and 0<small.lo<small.hi<1 and pos.lo>1
coef=five.monic().all_coeffs()
prod=-m.I(coef[-1])/(neg*small*pos);total=-m.I(coef[1])-neg-small-pos
assert 0<prod.lo<prod.hi<1 and (total**2-4*prod).hi<0
assert r(41,100)**2<prod.lo and prod.hi<r(43,100)**2 and small.hi<r(1,100)
coef=p.monic().all_coeffs();row=s.zeros(1,6);row[0,0]=1;res=s.zeros(1,6)
for k in range(6):
    scalar=sum(coef[j]*z**(k-j) for j in range(k+1))
    res+=scalar*row*J**(5-k)
assert all(s.Poly(e,z).is_zero for e in res*(J-z*s.eye(6))+p.as_expr()*row)
def left(root):
    vals=[m.horner(s.Poly(e,z).all_coeffs(),root) for e in res]
    normal=vals[max(range(6),key=lambda k:abs(vals[k].midpoint()))]
    return [e/normal for e in vals]
lm,lp=left(neg),left(pos)
initial,det=m.solve([list(B.row(i)) for i in range(4)]+[lm,lp],[0]*4+[m.dot(lm,sti),m.dot(lp,sti)])
first=[x+y for x,y in zip(m.mv(J.tolist(),initial),Di)]
fertility=first[2]/2
rel=[x-y for x,y in zip(initial,sti)]
w1=m.mv(J.tolist(),rel);w2=m.mv(J.tolist(),w1)
signal=w2[2]-small*w1[2]
quantities=dict(initial_fertility=fertility,stationary_total_households=sti[2]+sti[3],stationary_price=sti[0],boundary_det=det,complex_mode_signal=signal)
bounds=dict(initial_fertility=(r(84713,100000),r(84714,100000)),stationary_total_households=(r(479728,100000),r(479729,100000)),stationary_price=(r(310315,100000),r(310316,100000)),boundary_det=(r(15282,100000),r(15283,100000)),complex_mode_signal=(r(14668,100000),r(14669,100000)))
for key,val in quantities.items():
    lo,hi=bounds[key]
    assert lo<val.lo and val.hi<hi
print('CERTIFIED rational outward bounds:',{key:[str(lo),str(hi)] for key,(lo,hi) in bounds.items()})
print('Exact positive branch margins:',m.exact_anchor_checks(base)['exact_branch_margins'])
pn,Z,hh=m.anchor(4);numeric=m.linearization(pn,Z,hh);v=Z[[1,2]]
Ftn=m.complex_jacobian(lambda th:m.values(Z,v,dict(pn,theta=th[0]))[0],np.array([pn['theta']])).ravel()
Gtn=m.complex_jacobian(lambda th:m.values(Z,v,dict(pn,theta=th[0]))[1],np.array([pn['theta']])).ravel()
Fte=np.array(Ft.subs(L,s.log(r(10,3))),float).ravel();Gte=np.array(Gt.subs(L,s.log(r(10,3))),float).ravel()
err=max(np.max(abs(Ftn-Fte)),np.max(abs(Gtn-Gte)))
assert err<1e-12
print('Original-helper forcing partial max error:',err)
initial_ss=m.steady_state(pn);results=[]
for delta in [-1e-5,1e-5]:
    pp=dict(pn,theta=pn['theta']+delta)
    final=m.steady_state(pp);path=m.transition(initial_ss,final,24)
    checked=m.check_path(path,pp,optimize_dates=(0,1))
    r0=path['rows'][0]
    assert r0['past'] is initial_ss['households']
    results.append(dict(delta=delta,n0=float(r0['hh']['fertility']),Nstar=float(2*final['cohort']),Pstar=final['price']))
    print('Original-equation check:',dict(delta=delta,max_equilibrium_error=path['maximum_residual'],max_budget_error=checked['maximum_original_budget_error'],max_household_optimization_error=max(z['maximum_choice_error'] for z in checked['original_optimizations'])))
central=dict(initial_fertility=(results[1]['n0']-results[0]['n0'])/2e-5,stationary_total_households=(results[1]['Nstar']-results[0]['Nstar'])/2e-5,stationary_price=(results[1]['Pstar']-results[0]['Pstar'])/2e-5)
for key,value in central.items(): assert abs(value-float(quantities[key].midpoint()))<3e-9
print('Finite central derivatives:',central)
print('Finite small-shock example:',results[0])
print('PASS. Exact root, boundary, shock sign and original-equation checks; no files written.')


# Independent symbolic check of the dated household comparisons.
sa, skappa, schi, stheta, sk, sx, sn, sh, su = s.symbols(
    'alpha kappa chi vartheta k x n h u', positive=True)
sdh, sdu, sdw, sdtheta, sdx, sdn = s.symbols('dh du dw dtheta dx dn')
sspace = sh-skappa*sn
sfert = stheta/sn-schi/sx-sa*skappa/sspace
sbudget = (1+sk)*sdx+schi*sdn+su*sdh+sh*sdu-sdw
sdfert = s.diff(sfert,sx)*sdx+s.diff(sfert,sn)*sdn+s.diff(sfert,sh)*sdh+s.diff(sfert,stheta)*sdtheta
ssol = s.solve([sbudget,sdfert],[sdx,sdn],dict=True)[0][sdn]
sDelta = stheta/sn**2+schi**2/((1+sk)*sx**2)+sa*skappa**2/sspace**2
sclaimed = (sdtheta/sn+(sa*skappa/sspace**2-schi*su/((1+sk)*sx**2))*sdh+schi*(sdw-sh*sdu)/((1+sk)*sx**2))/sDelta
assert s.simplify(ssol-sclaimed)==0
sc2,sh2,sgamma,seps = s.symbols('c2 h2 gamma eps',positive=True)
scomp = sc2*((sh2/(sh2-seps))**sgamma-1)
sgain = s.log(1-scomp/sx)+sa*s.log(1+seps/sspace)
assert s.simplify(s.diff(sgain,seps).subs(seps,0)-(sa*sx/sspace-sgamma*sc2/sh2)/sx)==0
print('PASS. Independent dated fertility differential and compensated-allocation derivative.')
