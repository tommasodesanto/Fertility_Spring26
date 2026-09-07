"""Bounded read-only checks for the positive-child-cost stationary result.
Run from the project root: python3 code/model/tools/verify_simplified_olg_positive_costs.py
Imports the original helper; never runs its main or writes repository outputs.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import math
import numpy as np
import verify_simplified_olg_local_transition as original


def analytical(p, xO, nO, xR, nR, h, pi):
    q, alpha, beta, gamma, omega, theta, chi, kappa, a = (p[k] for k in ('q','alpha','beta','gamma','omega','theta','chi','kappa','rental_cap'))
    rhoO, rhoR = 1+beta*(1+gamma+omega), 1+beta*(1+omega)
    d=p['b']/(1-p['phi']); P=d/h; u=(1-q)*P
    C=beta*gamma*xO/(q*u*h); m=u*h/xO; v=(1+q)*u*a/xR
    A=alpha*h/(h-kappa*nO)+beta*gamma-v; F=m+beta*gamma-v
    gO=chi*nO/xO; gR=chi*nR/xR
    tO=kappa*nO/(h-kappa*nO); tR=kappa*nR/(a-kappa*nR)
    KO=gO+alpha*tO*(1+tO)+gO*gO/rhoO
    KR=gR+alpha*tR*(1+tR)+gR*gR/rhoR
    eta=alpha*tO*(1+tO)/KO; lamO=gO*m/(rhoO*KO); lamR=gR*v/(rhoR*KR)
    r=nR/nO; delta=a/h; epsR=(1-pi)/pi*r*lamR
    j=gO/rhoO; Gamma=1-j*eta; D0=1+m/rhoO-j*lamO
    Lambda=(1-delta)/(1-r); E0=(lamO+epsR)/(eta+epsR); E1=F/A
    young_upper=Lambda*(lamO+epsR)+max(Q(0),1-Lambda*(eta+epsR))*E1
    old_lower=C*(D0-max(Q(0),Gamma)*E1)
    Tcondition=1+C*Gamma-(1+C-2*delta)/(1-r)*(eta+epsR)
    U=pi*nO*(eta+epsR); V=pi*nO*(lamO+epsR); W=pi*nO/KO+(1-pi)*nR/KR
    k=pi*(1-pi)/p['sigma']; gap=nO-nR
    E=(V+gap*k*F)/(U+gap*k*A); pidot=k*(A*E-F)
    S=pi*h*(1+C)+(1-pi)*2*a
    Sdot=pi*h*(E+C*(Gamma*E-D0))+(h*(1+C)-2*a)*pidot
    Nphi=-2*p['Hbar']/S**2*Sdot/(1-p['phi'])
    L=math.log(float(nO/nR)); htheta=-(float(W)+float(gap*k)*L)/float(U+gap*k*A)
    pitheta=float(k)*(float(A)*htheta+L)
    Stheta=float(pi*h)*((1+float(C*Gamma))*htheta-float(C*j/KO))+float(h*(1+C)-2*a)*pitheta
    Ntheta=-float(2*p['Hbar']/S**2)*Stheta
    return dict(P=P,u=u,C=C,A=A,F=F,eta=eta,lambdaO=lamO,lambdaR=lamR,Gamma=Gamma,D0=D0,E0=E0,E1=E1,young_upper=young_upper,old_lower=old_lower,credit_margin=old_lower-young_upper,theta_condition=Tcondition,E=E,pidot=pidot,Sdot=Sdot,Nphi=Nphi,Ntheta=Ntheta,Pphi=P*(1-E)/(1-p['phi']),Ptheta=-float(P)*htheta,
                nO_dot=nO*(eta*E-lamO),nR_dot=nR*lamR*(E-1),mean_young_h_phi=(pi*h*E+(h-a)*pidot)/(1-p['phi']),mean_old_h_phi=(pi*C*h*(Gamma*E-D0)+(C*h-a)*pidot)/(1-p['phi']),hO_phi=h*E/(1-p['phi']),h2O_phi=C*h*(Gamma*E-D0)/(1-p['phi']),pi_phi=pidot/(1-p['phi']))


def original_check(p, xO, nO, xR, nR, h, pi, optimize):
    exact=analytical(p,xO,nO,xR,nR,h,pi)
    pp=original.parameters(**{k:float(v) for k,v in p.items()})
    P=float(exact['P']); hh=original.young_choices([P]*3,[0,0],pp)
    assert max(abs(hh['owner']['z'][0]-float(xO)),abs(hh['owner']['z'][2]-float(nO)),abs(hh['renter']['z'][0]-float(xR)),abs(hh['renter']['z'][2]-float(nR)))<1e-12
    diff=hh['owner']['utility']-hh['renter']['utility']
    xi=float(p['sigma'])*math.log(float(pi/(1-pi)))-diff
    pp['taste_weight']=math.exp(-xi/float(p['sigma']))
    ss=original.steady_state(pp);hh=ss['households'];P,T=ss['price'],ss['transfer']
    row=dict(t=0,prices=[P]*3,transfers=[T]*2,hh=hh,old_choices={tenure:original.old_choices(hh[tenure]['assets'],hh[tenure]['z'][1],P,P,T,pp,tenure=='owner') for tenure in ('owner','renter')},past=hh,transfer=T,price=P)
    verification=original.check_path(dict(rows=[row]),pp,optimize_dates=(0,) if optimize else ())
    numeric={}
    for parameter in ('phi','theta'):
        low=original.steady_state(dict(pp,**{parameter:pp[parameter]-1e-6}))
        high=original.steady_state(dict(pp,**{parameter:pp[parameter]+1e-6}))
        numeric[parameter]=dict(N=(2*high['cohort']-2*low['cohort'])/2e-6,P=(high['price']-low['price'])/2e-6)
        assert abs(numeric[parameter]['N']-float(exact['Nphi' if parameter=='phi' else 'Ntheta']))<2e-7
        assert abs(numeric[parameter]['P']-float(exact['Pphi' if parameter=='phi' else 'Ptheta']))<2e-7
    return dict(exact={k:str(v) for k,v in exact.items()},approximate={k:float(v) for k,v in exact.items()},pi=float(pi),xi_location=xi,allocations={tenure:hh[tenure]['z'].tolist() for tenure in ('owner','renter')},original_check=verification,finite_differences=numeric)



def run():
    positive=dict(q=Q(1,2),phi=Q(4,5),b=Q(9717,46250),y=Q(779829,370000),alpha=Q(2,5),beta=Q(2,5),gamma=Q(3,10),omega=Q(2,5),theta=Q(141,400),chi=Q(3,20),kappa=Q(1,2),nu=Q(2),tau=Q(0),rental_cap=Q(1,4),owner_cap=Q(2),sigma=Q(1),Hbar=Q(68104,68019))
    negative=dict(q=Q(4,5),phi=Q(4,5),b=Q(1,5),y=Q(199,62),alpha=Q(2,5),beta=Q(1,10),gamma=Q(1,5),omega=Q(139,155),theta=Q(12,5),chi=Q(2),kappa=Q(1,2),nu=Q(2),tau=Q(0),rental_cap=Q(1,10),owner_cap=Q(2),sigma=Q(1),Hbar=Q(377,664))
    a=analytical(positive,Q(1),Q(3,4),Q(99,74),Q(9,40),Q(1),Q(11,21))
    assert a['E0']<a['E1']<1 and a['credit_margin']>0 and a['theta_condition']>0 and a['C']>positive['rental_cap']
    b=analytical(negative,Q(1),Q(1),Q(51,20),Q(17,100),Q(1),Q(33,83))
    assert b['Nphi']<0 and b['C']>negative['rental_cap']
    assert b['nO_dot']<0 and b['nR_dot']<0
    assert b['Nphi']==-Q(4083925586718912,2831234115298355)
    return dict(positive=original_check(positive,Q(1),Q(3,4),Q(99,74),Q(9,40),Q(1),Q(11,21),True),negative=original_check(negative,Q(1),Q(1),Q(51,20),Q(17,100),Q(1),Q(33,83),True))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
