#!/usr/bin/env python3
"""Exact certificate for an interior mixed-tenure OLG transition.

Run --smoke, then the default bounded verification (four short finite paths).
The proof uses rational matrices, exact root isolation, and rational intervals.
Original dated household equations provide a separate derivative/budget check;
finite terminal closures are arithmetic checks, never convergence proofs.
--certificate-only omits finite paths. --figure writes a supplemental analytical
response using the stable linear solution, never the finite numerical paths.
--family-only certifies stability and stationary signs for all positive taste
dispersion, and the boundary/initial-fertility signs on the full interval [1,4].
The old all-owner verification driver is imported without modification.
"""

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as s
from verify_simplified_olg_local_transition import (
    parameters, young_choices, complex_jacobian, steady_state, transition,
    check_path,
)

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"


def anchor(sigma=1):
    p=parameters(q=.5,beta=.4,alpha=.4,gamma=.3,omega=.4,chi=.15,
                 theta=141/400,tau=467/9250,y=3521349281/1677802000,
                 Hbar=68104/68019,sigma=sigma)
    T=3975571/314587875
    hh=young_choices([1.]*3,[T,T],p)
    pi=11/21
    p['taste_weight']=(1-pi)/pi*np.exp((hh['owner']['utility']-hh['renter']['utility'])/sigma)
    hh=young_choices([1.]*3,[T,T],p)
    assert abs(hh['pi']-pi)<1e-12 and abs(hh['fertility']-.5)<1e-12
    Z=np.array([1.,9717/18500,1.,1.,.24*pi,.25*(1-pi)])
    return p,Z,hh


def values(Z,v,p):
    P,u,Y,O,M,R=Z
    un,Yn=v
    q,tau=p['q'],p['tau']
    Pn=((1+q*tau)*P-u)/q
    Pnn=((1+q*tau)*Pn-un)/q
    T=q*tau*P*p['Hbar']/(Y+O)
    Tn=q*tau*Pn*p['Hbar']/(Yn+Y)
    hh=young_choices([P,Pn,Pnn],[T,Tn],p)
    F=np.array([Y*hh['housing']+M/u+R-p['Hbar'],Yn-p['nu']*Y*hh['fertility']])
    G=np.array([Pn,un,Yn,Y,p['beta']*p['gamma']/q*Y*hh['pi']*hh['owner']['z'][0],
                p['rental_cap']*Y*(1-hh['pi'])])
    return F,G


def linearization(p,Z,hh):
    v=Z[[1,2]]
    F,G=values(Z,v,p)
    assert np.max(abs(F))<1e-12 and np.max(abs(G-Z))<1e-12
    Fz=complex_jacobian(lambda z:values(z,v,p)[0],Z)
    Fv=complex_jacobian(lambda vv:values(Z,vv,p)[0],v)
    Gz=complex_jacobian(lambda z:values(z,v,p)[1],Z)
    Gv=complex_jacobian(lambda vv:values(Z,vv,p)[1],v)
    J=Gz-Gv@np.linalg.solve(Fv,Fz)
    phi=np.array([p['phi']])
    Fphi=complex_jacobian(lambda f:values(Z,v,dict(p,phi=f[0]))[0],phi).ravel()
    Gphi=complex_jacobian(lambda f:values(Z,v,dict(p,phi=f[0]))[1],phi).ravel()
    Dphi=Gphi-Gv@np.linalg.solve(Fv,Fphi)
    derivative=np.linalg.solve(np.eye(6)-J,Dphi)
    pi=hh['pi']; L=p['gamma']/(1+p['gamma']+p['omega'])
    Hold=pi*hh['owner']['z'][1]
    B=np.zeros((4,6)); B[0,2]=1;B[1,3]=1;B[2,5]=1;B[3,4]=1
    B[3,0]=-L*(Hold+pi*p['q']*p['tau']*p['Hbar']/2)
    B[3,2]=B[3,3]=L*pi*p['q']*p['tau']*p['Hbar']/4
    roots,E=np.linalg.eig(J)
    stable=abs(roots)<1
    result=dict(J=J,Fv=Fv,roots=roots,B=B,derivative=derivative,Dphi=Dphi)
    if sum(stable)==4:
        Es=E[:,stable]
        a=np.linalg.solve(B@Es,-B@derivative)
        dz0=derivative+Es@a
        dz1=derivative+J@Es@a
        result.update(boundary=np.linalg.det(B@Es),first_fertility=dz1[2]/2,initial=dz0)
    return result


def exact_linearization(scale=s.S.One):
    R=s.Rational
    q=R(1,2); beta=alpha=omega=R(2,5); gamma=R(3,10); chi=R(3,20);kappa=R(1,2)
    theta=R(141,400);tau=R(467,9250);sigma=scale;phi=R(4,5)
    P=Y=O=Yn=s.S.One;u=R(9717,18500);H=R(68104,68019);pi=R(11,21);a=R(1,4)
    xO=hO=s.S.One;nO=R(3,4);xR=R(99,74);nR=R(9,40)
    M=beta*gamma/q*pi*xO;Rold=a*(1-pi)
    rho=1+beta*(1+gamma+omega);rhoR=1+beta*(1+omega)
    e=[s.eye(9)[:,i] for i in range(9)]
    dP,du,dY,dO,dM,dR,dun,dYn,dphi=e
    dPn=((1+q*tau)*dP-du)/q
    dT=q*tau*H/(Y+O)*dP-q*tau*P*H/(Y+O)**2*(dY+dO)
    dTn=q*tau*H/(Yn+Y)*dPn-q*tau*P*H/(Yn+Y)**2*(dYn+dY)
    dw=dT+q*dTn
    dhO=hO/(1-phi)*dphi-hO/P*dP


    def household_dx_dn(x,h,n,rho,dh,resource):
        A=theta/n**2+alpha*kappa**2/(h-kappa*n)**2
        matrix=s.Matrix([[rho,chi],[-chi/x**2,A]])
        rhs=s.Matrix.vstack(resource.T,(alpha*kappa/(h-kappa*n)**2*dh).T)
        ans=matrix.inv()*rhs
        return ans[0,:].T,ans[1,:].T


    dxO,dnO=household_dx_dn(xO,hO,nO,rho,dhO,dw-hO*du-u*dhO)
    dxR,dnR=household_dx_dn(xR,a,nR,rhoR,s.zeros(9,1),dw-a*du-q*a*dun)
    dWO=dw/xO-hO/xO*du+(alpha/(hO-kappa*nO)-u/xO)*dhO-beta*gamma/u*dun
    dWR=(dw-a*du-q*a*dun)/xR
    dpi=pi*(1-pi)/sigma*(dWO-dWR)
    dhbar=pi*dhO+(hO-a)*dpi
    dnbar=pi*dnO+(1-pi)*dnR+(nO-nR)*dpi
    dMnext=beta*gamma/q*(pi*xO*dY+Y*xO*dpi+Y*pi*dxO)
    dRnext=a*((1-pi)*dY-Y*dpi)
    dF1=(pi*hO+(1-pi)*a)*dY+Y*dhbar+dM/u-M/u**2*du+dR
    dF2=dYn-2*(R(1,2)*dY+Y*dnbar)
    F=s.Matrix.vstack(dF1.T,dF2.T)
    G=s.Matrix.vstack(dPn.T,dun.T,dYn.T,dY.T,dMnext.T,dRnext.T)
    J=G[:,:6]-G[:,6:8]*F[:,6:8].inv()*F[:,:6]
    Q=G[:,8]-G[:,6:8]*F[:,6:8].inv()*F[:,8]
    J=J.applyfunc(s.cancel); Q=Q.applyfunc(s.cancel)
    stationary=((s.eye(6)-J).inv()*Q).applyfunc(s.cancel)
    L=gamma/(1+gamma+omega)
    B=s.zeros(4,6); B[0,2]=1;B[1,3]=1;B[2,5]=1;B[3,4]=1
    B[3,0]=-L*(pi*hO+pi*q*tau*H/2)
    B[3,2]=B[3,3]=L*pi*q*tau*H/4


    z=s.symbols("z")
    return dict(J=J,Q=Q,B=B,stationary=stationary,poly=J.charpoly(z),implicit=F[:,6:8],dWO=dWO,dWR=dWR)

BITS=180
DEN=2**BITS


def frac(x):
    if isinstance(x,F):return x
    if isinstance(x,s.Rational):return F(int(x.p),int(x.q))
    return F(x)


class I:
    def __init__(self,lo,hi=None):
        lo,hi=frac(lo),frac(lo if hi is None else hi)
        self.lo=F((lo*DEN).__floor__(),DEN)
        self.hi=F((hi*DEN).__ceil__(),DEN)
        assert self.lo<=self.hi
    def __add__(self,x):
        x=iv(x);return I(self.lo+x.lo,self.hi+x.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,x):return self+-iv(x)
    def __rsub__(self,x):return iv(x)+-self
    def __mul__(self,x):
        x=iv(x);values=[a*b for a in (self.lo,self.hi) for b in (x.lo,x.hi)]
        return I(min(values),max(values))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=iv(x);assert x.lo>0 or x.hi<0,(float(x.lo),float(x.hi))
        return self*I(1/x.hi,1/x.lo)
    def __rtruediv__(self,x):return iv(x)/self
    def __pow__(self,n):
        assert n>=0
        answer=I(1)
        for _ in range(n):answer=answer*self
        return answer
    def pair(self):return [str(self.lo),str(self.hi)]
    def approx(self):return [float(self.lo),float(self.hi)]
    def midpoint(self):return (self.lo+self.hi)/2


def iv(x):return x if isinstance(x,I) else I(x)
def dot(a,b):return sum((iv(x)*iv(y) for x,y in zip(a,b)),I(0))
def mv(A,x):return [dot(row,x) for row in A]
def solve(A,b):
    A=[[iv(x) for x in row] for row in A];b=[iv(x) for x in b]
    n=len(b);det=I(1)
    for j in range(n):
        k=max(range(j,n),key=lambda i:abs(A[i][j].midpoint()))
        if k!=j:A[j],A[k]=A[k],A[j];b[j],b[k]=b[k],b[j];det=-det
        pivot=A[j][j];assert pivot.lo>0 or pivot.hi<0
        det=det*pivot
        A[j]=[x/pivot for x in A[j]];b[j]=b[j]/pivot
        for i in range(n):
            if i==j:continue
            multiple=A[i][j]
            A[i]=[x-multiple*y for x,y in zip(A[i],A[j])]
            b[i]=b[i]-multiple*b[j]
    return b,det


def horner(coeff,x):
    total=I(0)
    for c in coeff:total=total*x+c
    return total


def certificate(data):
    J,Q,B,stat,p=(data[k] for k in ['J','Q','B','stationary','poly'])
    z=p.gens[0]
    five=s.Poly(p.as_expr()/z,z)
    assert five.count_roots(-s.oo,s.oo)==3
    isolated=five.intervals(eps=s.Rational(1,10**35))
    assert len(isolated)==3 and all(m==1 for bounds,m in isolated)
    roots=[I(*bounds) for bounds,m in isolated]
    neg,small,pos=roots
    assert neg.hi<-1 and 0<small.lo<small.hi<1 and pos.lo>1
    coefficients=five.monic().all_coeffs()
    prod=-I(coefficients[-1])/(neg*small*pos)
    total=-I(coefficients[1])-neg-small-pos
    assert 0<prod.lo<prod.hi<1
    assert -(1+prod.lo)<total.lo<total.hi<1+prod.lo
    assert (total**2-4*prod).hi<0
    assert F(42,100)**2<prod.lo and prod.hi<F(44,100)**2
    assert small.hi<F(1,100)

    # A row of the polynomial resolvent is a left eigenvector at each root.
    coefficients=p.monic().all_coeffs()
    row=s.zeros(1,6);row[0,0]=1
    L=s.zeros(1,6)
    for k in range(6):
        scalar=sum(coefficients[j]*z**(k-j) for j in range(k+1))
        L+=scalar*row*J**(5-k)
    identity=L*(J-z*s.eye(6))+p.as_expr()*row
    assert all(s.Poly(v,z).is_zero for v in identity)
    def left(root):
        values=[horner(s.Poly(v,z).all_coeffs(),root) for v in L]
        at=max(range(6),key=lambda k:abs(values[k].midpoint()))
        normal=values[at]
        return [value/normal for value in values]
    Lminus,Lplus=left(neg),left(pos)
    matrix=[list(B.row(i)) for i in range(4)]+[Lminus,Lplus]
    rhs=[0]*4+[dot(Lminus,list(stat)),dot(Lplus,list(stat))]
    initial,det=solve(matrix,rhs)
    assert det.lo>0 or det.hi<0
    first=[a+b for a,b in zip(mv(J.tolist(),initial),list(Q))]
    fertility=first[2]/2
    assert fertility.lo>0 and stat[2]>0

    # A nonzero component in the dominant complex stable pair.
    relative=[x-y for x,y in zip(initial,list(stat))]
    w1=mv(J.tolist(),relative)
    w2=mv(J.tolist(),w1)
    complex_signal=w2[2]-small*w1[2]
    assert complex_signal.lo>0 or complex_signal.hi<0
    receipt=dict(arithmetic='Exact rational endpoints with outward rounding to dyadic multiples of 2^-180 after every operation.',
                 real_root_intervals=[r.pair() for r in roots],
                 real_root_approx=[r.approx() for r in roots],
                 remaining_quadratic_sum=total.pair(),remaining_quadratic_product=prod.pair(),
                 remaining_quadratic_approx={'sum':total.approx(),'product':prod.approx()},
                 stable_complex_modulus_bound=['21/50','11/25'],
                 boundary_determinant=det.pair(),boundary_determinant_approx=det.approx(),
                 initial_state_derivative=[r.pair() for r in initial],
                 initial_fertility_derivative=fertility.pair(),initial_fertility_approx=fertility.approx(),
                 stationary_cohort_derivative=str(stat[2]),
                 complex_mode_signal=complex_signal.pair(),complex_mode_signal_approx=complex_signal.approx())
    return receipt, initial


def exact_anchor_checks(data):
    """Rational original resource/FOC/inequality checks, both tenures."""
    r=s.Rational
    q=r(1,2); alpha=beta=omega=r(2,5); gamma=r(3,10)
    chi=r(3,20); kappa=r(1,2); theta=r(141,400); phi=r(4,5)
    tax=r(467,9250); b=r(1,5); y=r(3521349281,1677802000)
    stock=r(68104,68019); a=r(1,4); prob=r(11,21)
    u=r(9717,18500); rebate=r(3975571,314587875)
    assert u==1+q*tax-q and rebate==q*tax*stock/2
    result={}; young_h=old_h=mean_n=s.S.Zero
    for own,x,h,n in [(True,s.S.One,s.S.One,r(3,4)),
                      (False,r(99,74),a,r(9,40))]:
        rho=1+beta*(1+omega+(gamma if own else 0))
        assert rho*x+chi*n==y+b+(1+q)*rebate-u*h-(0 if own else q*u*a)
        assert theta/n==chi/x+alpha*kappa/(h-kappa*n)
        c2=beta*x/q; h2=gamma*c2/u if own else a; estate=omega*c2/q
        payment=((1-phi)+q*tax)*h if own else u*h
        saving=y+b+rebate-x-chi*n-payment
        assets=saving/q-(phi*h/q if own else 0)
        assert c2+u*h2+q*estate==assets+(h if own else 0)+rebate
        margins=dict(adult_goods=x,adult_space=h-kappa*n,saving=saving,
                     young_housing_wedge=alpha*x/(h-kappa*n)-u)
        if own:
            assert (1-phi)*h==b
            margins.update(owner_cap=2-h,retention=h-h2,estate=estate-h2,
                           old_housing_above_rental_cap=h2-a)
            assert gamma*c2/h2==u
        else:
            margins.update(old_renter_wedge=gamma*c2/a-u)
        assert all(value>0 for value in margins.values())
        weight=prob if own else 1-prob
        young_h+=weight*h;old_h+=weight*h2;mean_n+=weight*n
        result['owner' if own else 'renter']={key:str(value) for key,value in margins.items()}
    assert young_h+old_h==stock and mean_n==r(1,2)
    assert data['implicit'].det()!=0
    return dict(exact_branch_margins=result,implicit_determinant=str(data['implicit'].det()))


def path_checks(p, smoke=False):
    initial=steady_state(p)
    deltas=[1e-4] if smoke else [-1e-5,1e-5]
    horizons=[12] if smoke else [24,40]
    cases=[]; saved={}
    for delta in deltas:
        final=steady_state(dict(p,phi=p['phi']+delta))
        for horizon in horizons:
            path=transition(initial,final,horizon)
            checked=check_path(path,final['parameters'],optimize_dates=(0,1,horizon))
            cases.append(dict(delta=delta,horizon=horizon,
                              maximum_equilibrium_residual=path['maximum_residual'],
                              initial_fertility=float(path['rows'][0]['hh']['fertility']),
                              final_cohort=final['cohort'],**checked))
            saved[delta,horizon]=(final,path)
    comparisons={}
    if not smoke:
        smaller=saved[1e-5,24][1]['rows'];larger=saved[1e-5,40][1]['rows']
        comparisons['horizon_maximum_price_difference']=float(max(abs(a['price']-b['price']) for a,b in zip(smaller,larger)))
        dY=(saved[1e-5,40][0]['cohort']-saved[-1e-5,40][0]['cohort'])/(2e-5)
        dn=(saved[1e-5,40][1]['rows'][0]['hh']['fertility']-saved[-1e-5,40][1]['rows'][0]['hh']['fertility'])/(2e-5)
        comparisons.update(central_initial_fertility_derivative=float(dn),central_final_cohort_derivative=float(dY))
        assert comparisons['horizon_maximum_price_difference']<1e-10
    return dict(cases=cases,comparisons=comparisons)


def analytical_figure(numeric, out):
    """Evaluate the stable first-order solution; do not iterate unstable noise."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.linalg import schur
    J=numeric['J']; stat=numeric['derivative']
    tri,basis,nstable=schur(J,output='real',sort=lambda re,im:np.hypot(re,im)<1)
    assert nstable==4
    Es=basis[:,:4]; stable=tri[:4,:4]
    coefficients=np.linalg.solve(numeric['B']@Es,-numeric['B']@stat)
    states=[]
    for t in range(21):
        states.append(stat+Es@np.linalg.matrix_power(stable,t)@coefficients)
    states=np.array(states)
    fertility=(states[1:,2]-states[:-1,2])/2
    population=states[:-1,2]+states[:-1,3]
    dates=np.arange(20)
    plt.rcParams.update({'font.family':'serif','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(9.2,3.1),layout='constrained')
    axes[0].axhline(0,color='.45',lw=.7)
    axes[0].plot(dates,fertility,color='#235c8e',marker='o',markersize=3)
    axes[0].set(ylabel=r'$\partial\bar n_t/\partial\phi$',xlabel='Young cohort',title='Fertility relative to replacement',xlim=(0,12))
    axes[1].axhline(2*stat[2],color='.45',lw=.8,ls='--',label='New steady state')
    axes[1].plot(dates,population,color='#235c8e',marker='o',markersize=3)
    axes[1].set(ylabel=r'$\partial(Y_t+O_t)/\partial\phi$',xlabel='Date',title='Change in household population',xlim=(0,12))
    axes[1].legend(frameon=False,fontsize=9,loc='lower right')
    paths=[]
    for ext in ('pdf','png'):
        path=out/f'mixed_transition_analytical_figure.{ext}'
        fig.savefig(path,dpi=190,bbox_inches='tight');paths.append(str(path.relative_to(ROOT)))
    plt.close(fig)
    return dict(paths=paths,source='Stable analytical derivative, with Schur projection to avoid accumulating numerical unstable components.',
                first_order_fertility=fertility.tolist(),first_order_population=population.tolist())


def resolvent(data):
    p,J=data['poly'],data['J'];z=p.gens[0];co=p.monic().all_coeffs()
    row=s.zeros(1,6);row[0,0]=1
    ans=s.zeros(1,6)
    for k in range(6):ans+=sum(co[j]*z**(k-j) for j in range(k+1))*row*J**(5-k)
    return ans.applyfunc(s.expand)

def family_checks():
    start=time.monotonic();sig=s.symbols('sigma',positive=True)
    fam=exact_linearization(sig);J,Q,B,stat,p=(fam[k] for k in ['J','Q','B','stationary','poly']);z=p.gens[0]
    five=s.Poly(p.as_expr()/z,z)
    assert all(s.diff(v,sig,2)==0 for v in J) and all(s.diff(v,sig,2)==0 for v in Q)
    # Full-family polynomial and Routh certificate.
    v=s.symbols('v');C=s.Poly(s.cancel((1-v)**5*five.as_expr().subs(z,(1+v)/(1-v))),v)
    co=C.all_coeffs();routh=[co[0::2],co[1::2]]
    for k in range(4):
        prev,cur=routh[-2:];row=[s.factor((cur[0]*prev[j+1]-prev[0]*cur[j+1])/cur[0]) for j in range(2)]+[s.S.Zero];routh.append(row)
    expected=[-1,-1,-1,1,1,-1]
    for row,wanted in zip(routh,expected):
        num,den=s.fraction(s.cancel(row[0]));n=s.Poly(num,sig);d=s.Poly(den,sig)
        assert all(s.sign(c)==wanted for c in n.all_coeffs()) and all(c>0 for c in d.all_coeffs())
    # Negative discriminant for all sig>0, so the three-real-root count cannot change.
    disc=s.Poly(s.discriminant(five.as_expr(),z),sig)
    assert all(c<0 for c in disc.all_coeffs())
    assert s.Poly(five.as_expr().subs(sig,1),z).count_roots(-s.oo,s.oo)==3
    # The rank-one parameter change gives a left eigenvector polynomial
    # independent of sigma, reducing repeated-parameter interval expansion.
    A=J.subs(sig,0); slope=J.diff(sig); row=slope[1,:]
    column=slope[:,0]/row[0]
    assert slope==column*row
    p0=s.Poly(p.as_expr().subs(sig,0),z);co=p0.monic().all_coeffs()
    L=s.zeros(1,6)
    for k in range(6):L+=sum(co[j]*z**(k-j) for j in range(k+1))*row*A**(5-k)
    L=L.applyfunc(s.expand)
    identity=L*(J-z*s.eye(6))+p.as_expr()*row
    assert all(s.Poly(s.expand(v),z,sig).is_zero for v in identity)
    Lcoeff=[[s.Poly(c,sig).all_coeffs() for c in s.Poly(v,z).all_coeffs()] for v in L]



    def compile_expression(expr):
        num,den=s.fraction(s.cancel(expr))
        return ([frac(c) for c in s.Poly(num,sig).all_coeffs()],
                [frac(c) for c in s.Poly(den,sig).all_coeffs()])
    def evaluate(compiled,x):
        numerator,denominator=compiled
        return horner(numerator,x)/horner(denominator,x)
    def exact_horner(coefficients,x):
        total=F(0)
        for coefficient in coefficients:total=total*x+coefficient
        return total
    jcompiled=[[compile_expression(x) for x in row] for row in J.tolist()]
    qcompiled=[compile_expression(x) for x in Q]
    scompiled=[compile_expression(x) for x in stat]
    pcompiled=[compile_expression(x) for x in five.monic().all_coeffs()]


    endpoint_cache={}
    def endpoint(x):
        if x not in endpoint_cache:
            poly=s.Poly(five.as_expr().subs(sig,s.Rational(x.numerator,x.denominator)),z)
            endpoint_cache[x]=([bounds for bounds,m in poly.intervals(eps=s.Rational(1,10**25))],[frac(c) for c in poly.all_coeffs()])
        return endpoint_cache[x]


    def cell(lo,hi):
        scale=I(lo,hi)
        endpoints=[endpoint(lo)[0],endpoint(hi)[0]]
        roots=[]
        for j,direction in enumerate([1,-1,1]):
            lower=min(e[j][0] for e in endpoints);upper=max(e[j][1] for e in endpoints)
            for edge in [lo,hi]:
                assert direction*exact_horner(endpoint(edge)[1],frac(lower))<0
                assert direction*exact_horner(endpoint(edge)[1],frac(upper))>0
            roots.append(I(lower,upper))
        neg,small,pos=roots
        assert neg.hi<-1 and pos.lo>1 and 0<small.lo<small.hi<F(1,100)
        coefficients=five.monic().all_coeffs()
        prod=-evaluate(pcompiled[-1],scale)/(neg*small*pos)
        total=-evaluate(pcompiled[1],scale)-neg-small-pos
        assert F(40,100)**2<prod.lo and prod.hi<F(45,100)**2
        assert (total**2-4*prod).hi<0
        def left(root):
            values=[horner([horner(c,scale) for c in term],root) for term in Lcoeff]
            at=max(range(6),key=lambda k:abs(values[k].midpoint()))
            return [v/values[at] for v in values]
        lm,lp=left(neg),left(pos)
        si=[evaluate(x,scale) for x in scompiled]
        ji=[[evaluate(x,scale) for x in row] for row in jcompiled]
        qi=[evaluate(x,scale) for x in qcompiled]
        init,det=solve([list(B.row(i)) for i in range(4)]+[lm,lp],[0]*4+[dot(lm,si),dot(lp,si)])
        first=[a+b for a,b in zip(mv(ji,init),qi)]
        fertility=first[2]/2
        assert fertility.lo>0 and si[2].lo>0
        return dict(scale=[str(lo),str(hi)],fertility=fertility.pair(),
                    stationary=si[2].pair(),boundary=det.pair())



    pending=[(F(1),F(4))];passed=[];failures=0
    while pending:
        if time.monotonic()-start>600 or failures>12000:
            raise RuntimeError('Family proof exceeded its 600-second or 12000-subdivision bound; no certificate issued.')
        lo,hi=pending.pop()
        try:passed.append(cell(lo,hi))
        except AssertionError:
            failures+=1;mid=(lo+hi)/2;pending.extend([(mid,hi),(lo,mid)])
        if (len(passed)+failures)%100==0:print('progress',len(passed),len(pending),failures,time.monotonic()-start,'last',float(lo),float(hi),flush=True)
    else:
        pass

    passed.sort(key=lambda cell:frac(cell['scale'][0]))
    assert frac(passed[0]['scale'][0])==1 and frac(passed[-1]['scale'][1])==4
    assert all(frac(a['scale'][1])==frac(b['scale'][0]) for a,b in zip(passed,passed[1:]))
    direction=s.Matrix([*stat,stat[1],stat[2],1])
    stationary_expressions=dict(cohort=s.factor(stat[2]),
                                owner_value=s.factor((fam['dWO'].T*direction)[0]),
                                renter_value=s.factor((fam['dWR'].T*direction)[0]))
    for name,expr in stationary_expressions.items():
        numerator,denominator=s.fraction(expr)
        wanted=1 if name=='cohort' else -1
        assert all(s.sign(c)==wanted for c in s.Poly(numerator,sig).all_coeffs())
        assert all(c>0 for c in s.Poly(denominator,sig).all_coeffs())
    original_comparisons=[]
    for scale in [s.Rational(1,4),1,2,4]:
        param,state,hh=anchor(float(scale)); numeric=linearization(param,state,hh)
        matrix_error=float(np.max(abs(np.array(J.subs(sig,scale),float)-numeric['J'])))
        assert matrix_error<1e-10
        low=steady_state(dict(param,phi=param['phi']-1e-5))
        high=steady_state(dict(param,phi=param['phi']+1e-5))
        errors={}
        for name,expr in stationary_expressions.items():
            if name=='cohort': observed=(high['cohort']-low['cohort'])/(2e-5)
            else:
                tenure='owner' if name=='owner_value' else 'renter'
                observed=(high['households'][tenure]['utility']-low['households'][tenure]['utility'])/(2e-5)
            errors[name]=abs(float(expr.subs(sig,scale))-observed)
            assert errors[name]<1e-6
        original_comparisons.append(dict(scale=str(scale),jacobian_error=matrix_error,stationary_derivative_errors=errors))
    return dict(scope='Four stable/two unstable roots and three real nonzero roots for every positive sigma. Positive stationary cohort and negative values in both tenures for every positive sigma. Transverse initial boundary and positive initial fertility certified for every sigma in [1,4]. The oscillation theorem is certified near sigma=1, not across this whole interval.',
                arithmetic='Rational interval proof; real-root brackets hold throughout each sigma interval because characteristic coefficients are affine in sigma. Rank-one left eigenvector polynomial is independent of sigma and its identity is checked exactly.',
                interval=['1','4'],cells=passed,subdivisions=failures,elapsed_seconds=time.monotonic()-start,
                characteristic_polynomial=str(p.as_expr()),
                routh_first_column=[str(row[0]) for row in routh],routh_signs=expected,
                discriminant_coefficients=[str(c) for c in disc.all_coeffs()],
                stationary_derivatives={key:str(value) for key,value in stationary_expressions.items()},
                original_equation_comparisons=original_comparisons,
                global_bounds={key:[str(min(frac(c[key][0]) for c in passed)),str(max(frac(c[key][1]) for c in passed))] for key in ['fertility','stationary','boundary']})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--certificate-only',action='store_true')
    parser.add_argument('--figure',action='store_true')
    parser.add_argument('--family-only',action='store_true')
    args=parser.parse_args();start=time.monotonic();OUT.mkdir(parents=True,exist_ok=True)
    if args.family_only:
        family=family_checks()
        family['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        family['original_helper_sha256']=hashlib.sha256(Path(__file__).with_name('verify_simplified_olg_local_transition.py').read_bytes()).hexdigest()
        family['created_utc']=datetime.now(timezone.utc).isoformat()
        output=OUT/'mixed_transition_family_certificate.json'
        output.write_text(json.dumps(family,indent=2)+'\n')
        print(json.dumps(dict(receipt=str(output.relative_to(ROOT)),cells=len(family['cells']),elapsed_seconds=family['elapsed_seconds']),indent=2))
        return
    data=exact_linearization(); certified,initial=certificate(data)
    original=exact_anchor_checks(data)
    p,Z,hh=anchor();numeric=linearization(p,Z,hh)
    jac_error=float(np.max(abs(np.array(data['J'],float)-numeric['J'])))
    policy_error=float(np.max(abs(np.array(data['Q'],float).ravel()-numeric['Dphi'])))
    assert jac_error<1e-11 and policy_error<1e-11
    assert all(abs(complex(x).imag)<1e-10 for x in numeric['initial'])
    derivative_error=float(max(abs(float(a.midpoint())-complex(b).real) for a,b in zip(initial,numeric['initial'])))
    assert derivative_error<1e-10
    receipt=dict(created_utc=datetime.now(timezone.utc).isoformat(),
                 proof_scope='Local nonlinear convergence and strict initial-fertility/final-population signs at an interior mixed equilibrium. Finite paths do not prove convergence.',
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 original_helper_sha256=hashlib.sha256(Path(__file__).with_name('verify_simplified_olg_local_transition.py').read_bytes()).hexdigest(),
                 parameters=p,owner_probability=float(hh['pi']),taste_location=float(-p['sigma']*np.log(p['taste_weight'])),
                 certificate=certified,original_equations=original,
                 cross_checks=dict(jacobian_maximum_difference=jac_error,policy_derivative_maximum_difference=policy_error,initial_derivative_maximum_difference=derivative_error),
                 exact_matrices={key:[[str(v) for v in row] for row in data[key].tolist()] for key in ('J','Q','B','stationary','implicit')},
                 characteristic_polynomial=str(data['poly'].as_expr()))
    if not args.certificate_only:
        receipt['finite_original_paths']=path_checks(p,args.smoke)
        comparison=receipt['finite_original_paths']['comparisons']
        if comparison:
            comparison['initial_fertility_derivative_error']=abs(comparison['central_initial_fertility_derivative']-float(iv(certified['initial_fertility_derivative'][0]).midpoint()))
            comparison['final_cohort_derivative_error']=abs(comparison['central_final_cohort_derivative']-float(data['stationary'][2]))
            assert comparison['initial_fertility_derivative_error']<1e-6
            assert comparison['final_cohort_derivative_error']<1e-6
    if args.figure:receipt['figure']=analytical_figure(numeric,OUT)
    receipt['elapsed_seconds']=time.monotonic()-start
    suffix='smoke' if args.smoke else ('certificate' if args.certificate_only else 'checks')
    output=OUT/f'mixed_transition_{suffix}.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(receipt=str(output.relative_to(ROOT)),elapsed_seconds=receipt['elapsed_seconds'],
                         owner_probability=receipt['owner_probability'],cross_checks=receipt['cross_checks'],
                         initial_fertility=certified['initial_fertility_approx'],final_cohort=certified['stationary_cohort_derivative']),indent=2))


if __name__=='__main__':
    main()
