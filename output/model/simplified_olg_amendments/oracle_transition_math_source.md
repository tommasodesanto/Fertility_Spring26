# Pro verifier as received

Original source SHA-256: 69d3fe3fbc2202d10cecfe9b43f03332b80c332494b8786b53ba38da9208d36f.
The following source is preserved unchanged, including its final newline.

~~~python
#!/usr/bin/env python3
"""Exact half-line certificate for the finite mixed-tenure OLG transition.

Run with Python 3.10 or newer, without -O. Uses only the standard library.
The companion JSON contains fixed dyadic coefficients of an approximate
inverse. Their provenance is immaterial to acceptance: both inverse defects
are checked on the entire half-line by exact rational arithmetic.
There is no simulated equilibrium path or terminal-state condition.
"""
from fractions import Fraction as Q
from math import isqrt, factorial
from pathlib import Path
import json
import hashlib
import sys

q=Q(1,2); beta=alpha=omega=Q(2,5); gamma=Q(3,10)
chi=Q(3,20); kap=Q(1,2); theta0=Q(141,400); phi0=Q(4,5)
b=Q(1,5); y=Q(3521349281,1677802000); stock=Q(68104,68019)
tax=Q(467,9250); taxcoef=q*tax*stock; ap=1+q*tax
rhoO=1+beta*(1+gamma+omega); rhoR=1+beta*(1+omega)
a=Q(1,4); nu=Q(2); sigma=Q(4); pish=Q(11,21)
u0=ap-q; T0=taxcoef/2
xO0=Q(1); xR0=Q(99,74); nO0=Q(3,4); nR0=Q(9,40)
sO0=1-kap*nO0; sR0=a-kap*nR0
LL=gamma/(1+gamma+omega)
aO0=beta*(1+gamma+omega)*xO0/q-1-T0
ar0=beta*(1+omega)*xR0/q+u0*a-T0


class Jet:
 def __init__(self,v,d=None): self.v=Q(v);self.d={} if d is None else {k:Q(z) for k,z in d.items() if z}
 def __add__(self,o):
  o=j(o);d=self.d.copy()
  for k,v in o.d.items():d[k]=d.get(k,Q(0))+v
  return Jet(self.v+o.v,d)
 __radd__=__add__
 def __neg__(self):return Jet(-self.v,{k:-v for k,v in self.d.items()})
 def __sub__(self,o):return self+-j(o)
 def __rsub__(self,o):return j(o)+-self
 def __mul__(self,o):
  o=j(o);d={k:v*o.v for k,v in self.d.items()}
  for k,v in o.d.items():d[k]=d.get(k,Q(0))+v*self.v
  return Jet(self.v*o.v,d)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=j(o);return self*Jet(1/o.v,{k:-v/o.v**2 for k,v in o.d.items()})
 def __rtruediv__(self,o):return j(o)/self
 def __pow__(self,n):return Jet(self.v**n,{k:n*self.v**(n-1)*v for k,v in self.d.items()})

def j(o):return o if isinstance(o,Jet) else Jet(o)

def var(t,k):return Jet(1,{(t,k):1})

def lin_young(t,P,Y):
 u=ap*P(t)-q*P(t+1); un=ap*P(t+1)-q*P(t+2)
 T=taxcoef*P(t)/(Y(t)+Y(t-1));Tn=taxcoef*P(t+1)/(Y(t+1)+Y(t))
 h=b/(1-phi0)/P(t);w=y+b+T+q*Tn
 def fn(h,rho,x,n,B):
  s=h.v-kap*n;An=theta0/n**2+alpha*kap**2/s**2
  z=(B-B.v)*chi/(rho*x*x)+(h-h.v)*alpha*kap/s**2
  nj=Jet(n)+z/(An+chi**2/(rho*x*x))
  xj=(B-chi*nj)/rho
  assert xj.v==x
  return xj,nj
 xo,no=fn(h,rhoO,xO0,nO0,w-u*h)
 xr,nr=fn(j(a),rhoR,xR0,nR0,w-a*u-q*a*un)
 dw=w-w.v;du=u-u.v;dun=un-un.v;dh=h-h.v
 dd=(dw-h.v*du)/xO0+(alpha/sO0-u.v/xO0)*dh-beta*gamma/un.v*dun-(dw-a*du-q*a*dun)/xR0
 pi=Jet(pish)+pish*(1-pish)/sigma*dd
 return dict(u=u,un=un,T=T,h=h,xo=xo,no=no,xr=xr,nr=nr,pi=pi)

def lin_res(t,boundary=False):
 P=lambda k:var(k,0)
 Y=lambda k:var(k-1,1) if (not boundary or k>=1) else j(1)
 h=lin_young(t,P,Y)
 young=Y(t)*(h['pi']*h['h']+(1-h['pi'])*a)
 if boundary and t==0:
  old=pish*LL*(aO0+P(0)+h['T'])/h['u']+(1-pish)*a
 else:
  prev=lin_young(t-1,P,Y)
  old=Y(t-1)*(beta*gamma/q*prev['pi']*prev['xo']/h['u']+a*(1-prev['pi']))
 H=young+old-stock
 D=Y(t+1)-nu*Y(t)*(h['pi']*h['no']+(1-h['pi'])*h['nr'])
 assert H.v==0 and D.v==0,(H.v,D.v)
 return H,D

stencil=[r.d for r in lin_res(0)]
bound=[r.d for r in lin_res(0,True)]

def lentry(t,k,s,l):
 row=bound[k] if t==0 else stencil[k]
 return row.get((s if t==0 else s-t,l),Q(0)) if s>=0 else Q(0)

BITS=90
SCALE=1<<BITS

def floorq(x):return Q((x.numerator*SCALE)//x.denominator,SCALE)

def ceilq(x):return -floorq(-x)

class IV:
 def __init__(self,lo,hi=None):
  if isinstance(lo,IV): self.lo,self.hi=lo.lo,lo.hi;return
  self.lo=Q(lo);self.hi=Q(lo if hi is None else hi)
  assert self.lo<=self.hi
 def __add__(self,o):
  if hasattr(o,"v"):return NotImplemented
  o=iv(o);return IV(floorq(self.lo+o.lo),ceilq(self.hi+o.hi))
 __radd__=__add__
 def __neg__(self):return IV(-self.hi,-self.lo)
 def __sub__(self,o):
  if hasattr(o,"v"):return NotImplemented
  return self+-iv(o)
 def __rsub__(self,o):return iv(o)+-self
 def __mul__(self,o):
  if hasattr(o,"v"):return NotImplemented
  o=iv(o); vals=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
  return IV(floorq(min(vals)),ceilq(max(vals)))
 __rmul__=__mul__
 def __truediv__(self,o):
  if hasattr(o,"v"):return NotImplemented
  o=iv(o);assert o.lo*o.hi>0,('division',o)
  return self*IV(floorq(1/o.hi),ceilq(1/o.lo))
 def __rtruediv__(self,o):return iv(o)/self
 def __pow__(self,n):
  assert isinstance(n,int)
  if n<0:return iv(1)/(self**(-n))
  if n==0:return iv(1)
  vals=[self.lo**n,self.hi**n]
  if n%2==0 and self.lo<=0<=self.hi: vals.append(Q(0))
  return IV(floorq(min(vals)),ceilq(max(vals)))
 def absmax(self):return max(abs(self.lo),abs(self.hi))
 def __repr__(self):return f'[{float(self.lo):.8g},{float(self.hi):.8g}]'
 def sqrt(self):
  assert self.lo>0
  l=isqrt((self.lo.numerator*SCALE*SCALE)//self.lo.denominator)
  h=isqrt((self.hi.numerator*SCALE*SCALE)//self.hi.denominator)+1
  return IV(Q(l,SCALE),Q(h,SCALE))
 def log(self):
  def endpoint(x):
   assert x>0
   z=(x-1)/(x+1);zz=z*z;power=z;ss=Q(0);N=0
   while True:
    ss+=power/(2*N+1); N+=1; power*=zz
    rem=2*abs(power)/((2*N+1)*(1-zz))
    if rem<Q(1,1<<(BITS+10)):break
    assert N<1000
   return IV(floorq(2*ss-rem),ceilq(2*ss+rem))
  return IV(endpoint(self.lo).lo,endpoint(self.hi).hi)
 def exp(self):
  def endpoint(x):
   assert abs(x)<1,('exp range',x)
   ss=Q(1);power=Q(1)
   for n in range(1,34):
    power*=x/n;ss+=power
   # Taylor remainder <= e |x|^34/34!, and e<3.
   rem=3*abs(x)**34/factorial(34)
   return IV(floorq(ss-rem),ceilq(ss+rem))
  return IV(endpoint(self.lo).lo,endpoint(self.hi).hi)

def iv(o):return o if isinstance(o,IV) else IV(o)

def symrad(c,r):return IV(Q(c)-Q(r),Q(c)+Q(r))

class V:
 def __init__(self,v,d=None):self.v=iv(v);self.d={} if d is None else {k:iv(z) for k,z in d.items()}
 def __add__(self,o):
  o=v(o);d=self.d.copy()
  for k,z in o.d.items():d[k]=d.get(k,iv(0))+z
  return V(self.v+o.v,d)
 __radd__=__add__
 def __neg__(self):return V(-self.v,{k:-z for k,z in self.d.items()})
 def __sub__(self,o):return self+-v(o)
 def __rsub__(self,o):return v(o)+-self
 def __mul__(self,o):
  o=v(o);d={k:z*o.v for k,z in self.d.items()}
  for k,z in o.d.items():d[k]=d.get(k,iv(0))+z*self.v
  return V(self.v*o.v,d)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=v(o);return self*V(1/o.v,{k:-z/o.v**2 for k,z in o.d.items()})
 def __rtruediv__(self,o):return v(o)/self
 def __pow__(self,n):return V(self.v**n,{k:n*self.v**(n-1)*z for k,z in self.d.items()})

def v(o):return o if isinstance(o,V) else V(o)

LOG_NR=iv(nO0/nR0).log()

def young_iv(t,P,Y,theta,phi):
 u=ap*P(t)-q*P(t+1);un=ap*P(t+1)-q*P(t+2)
 T=taxcoef*P(t)/(Y(t)+Y(t-1));Tn=taxcoef*P(t+1)/(Y(t+1)+Y(t))
 h=b/(1-phi)/P(t);w=y+b+T+q*Tn
 def fn(h,rho,B):
  hh=h.v;BB=B.v;th=theta.v
  A=chi*kap*(th+rho+alpha)
  C=th*BB*hh
  E=kap*BB*(th+alpha)+chi*hh*(th+rho)
  n=2*C/(E+(E**2-4*A*C).sqrt())
  x=(BB-chi*n)/rho;s=hh-kap*n
  assert x.lo>0 and s.lo>0 and n.lo>0
  den=th/n**2+alpha*kap**2/s**2+chi**2/(rho*x**2)
  keys=B.d.keys()|h.d.keys()|theta.d.keys()
  dn={k:(theta.d.get(k,iv(0))/n+alpha*kap/s**2*h.d.get(k,iv(0))+chi/(rho*x**2)*B.d.get(k,iv(0)))/den for k in keys}
  nn=V(n,dn);xx=(B-chi*nn)/rho
  return xx,nn
 xo,no=fn(h,rhoO,w-u*h)
 xr,nr=fn(v(a),rhoR,w-a*u-q*a*un)
 so=h.v-kap*no.v;sr=iv(a)-kap*nr.v
 dv=(rhoO*(xo.v/xO0).log()-rhoR*(xr.v/xR0).log()+alpha*((so/sO0).log()-(sr/sR0).log())+theta.v*((no.v/nO0).log()-(nr.v/nR0).log())+(theta.v-theta0)*LOG_NR-beta*gamma*(un.v/u0).log())
 ex=(dv/sigma).exp();pi=11*ex/(10+11*ex)
 keys=w.d.keys()|u.d.keys()|un.d.keys()|h.d.keys()|theta.d.keys()
 dd={k:(w.d.get(k,iv(0))-h.v*u.d.get(k,iv(0)))/xo.v+(alpha/so-u.v/xo.v)*h.d.get(k,iv(0))-beta*gamma/un.v*un.d.get(k,iv(0))-(w.d.get(k,iv(0))-a*u.d.get(k,iv(0))-q*a*un.d.get(k,iv(0)))/xr.v+theta.d.get(k,iv(0))*(LOG_NR+(no.v/nO0).log()-(nr.v/nR0).log()) for k in keys}
 pp=V(pi,{k:pi*(1-pi)/sigma*z for k,z in dd.items()})
 return dict(u=u,un=un,T=T,Tn=Tn,h=h,xo=xo,no=no,xr=xr,nr=nr,pi=pp,so=so,sr=sr,w=w)

BASEI=[iv(z) for z in [1,1,pish,aO0,1,ar0]]

def residual_iv(t,r,dt,df,I=None):
 if I is None:I=BASEI
 P=lambda k: V(symrad(1,r),{(k,0):1})
 Y=lambda k: V(symrad(1,r),{(k-1,1):1}) if k>=1 else V(I[0] if k==0 else I[1])
 theta=V(IV(theta0-dt,theta0),{'theta':1});phi=V(IV(phi0,phi0+df),{'phi':1})
 h=young_iv(t,P,Y,theta,phi)
 young=Y(t)*(h['pi']*h['h']+(1-h['pi'])*a)
 if t==0:
  O,oldpi,oldasset,oldtitle=I[1:5]
  old=O*oldpi*LL*(oldasset+P(0)*oldtitle+h['T'])/h['u']+O*(1-oldpi)*a
 else:
  prev=young_iv(t-1,P,Y,theta,phi)
  old=Y(t-1)*(beta*gamma/q*prev['pi']*prev['xo']/h['u']+a*(1-prev['pi']))
 H=young+old-stock
 D=Y(t+1)-nu*Y(t)*(h['pi']*h['no']+(1-h['pi'])*h['nr'])
 return H,D

def get_rows(r,dt,df,I=None):
 out=[]
 for t in range(4):out.append(residual_iv(t,r,dt,df,I))
 return out

def generator_bounds(r,dt):
 # All baseline cohorts, including date zero; independent prices/populations in box.
 P=lambda k:V(symrad(1,r))
 Y=lambda k:V(symrad(1,r))
 h=young_iv(0,P,Y,V(IV(theta0-dt,theta0)),V(phi0))
 # Post-mortgage NET assets; actual forecasts under phi0 are retained.
 assetO=(y+b+h['T']-h['xo']-chi*h['no']-(1+q*tax)*P(0)*h['h'])/q
 assetR=(y+b+h['T']-h['xr']-chi*h['nr']-h['u']*a)/q
 return [symrad(1,r),symrad(1,r),h['pi'].v,assetO.v,h['h'].v,assetR.v]

K=24
ABITS=46
AS=1<<ABITS
PATH=Path(__file__).with_name("finite_mixed_preconditioner.json")

def load_A():
 data=json.load(open(PATH));assert data['K']==K and data['dyadic_bits']==ABITS
 cv=lambda row:{int(s):Q(c,AS) for s,c in row.items()}
 return list(map(cv,data['boundary'])),list(map(cv,data['tail']))

AB,AT=load_A()
assert len(AB)==2*K and len(AT)==2
for i,row in enumerate(AB):
    t=i//2
    assert all(isinstance(s,int) and 0 <= s//2 <= t+K for s in row)
for row in AT:
    assert all(-K <= s//2 <= K for s in row)

def Arow(i):
 t,k=divmod(i,2)
 return AB[i] if t<K else {2*t+s:c for s,c in AT[k].items() if 2*t+s>=0}

def Lrow(i):
 t,k=divmod(i,2)
 row=bound[k] if t==0 else stencil[k]
 return {2*(j if t==0 else j+t)+l:c for (j,l),c in row.items() if (j if t==0 else j+t)>=0}

def exact_operator_checks():
 # Outside this finite set, every product row is a translation of the last
 # interior rows: all operators have finite bandwidth and a constant tail.
 nr=2*(2*K+8)
 def product_error(R,S):
  errs=[]
  for i in range(nr):
   out={i:Q(1)}
   for j,c in R(i).items():
    for k,d in S(j).items():out[k]=out.get(k,Q(0))-c*d
   errs.append(sum(abs(c) for c in out.values()))
  return max(errs),errs
 al,als=product_error(Arow,Lrow);la,las=product_error(Lrow,Arow)
 an=max(sum(abs(z) for z in Arow(i).values()) for i in range(2*K+2))
 return al,la,an,als

def derivative_bounds(rows,al,als):
 # Bounds ||I-A DF|| by rows. Triangle inequality is intentionally conservative.
 diff=[]
 for t,rs in enumerate(rows):
  diff.append([])
  for k,z in enumerate(rs):
   keys={key for key in z.d if isinstance(key,tuple)}|{divmod(j,2) for j in Lrow(2*t+k)}
   e=sum((z.d.get(key,iv(0))-lentry(t,k,*key)).absmax() for key in keys)
   diff[-1].append(e)
 norm=[];forcings={'theta':[],'phi':[]}
 for i in range(2*(K+5)):
  err=als[i] if i<len(als) else al
  for j,c in Arow(i).items():
   t,k=divmod(j,2);err+=abs(c)*diff[min(t,3)][k]
  norm.append(err)
  for par in forcings:
   z=iv(0)
   for j,c in Arow(i).items():
    t,k=divmod(j,2);z-=c*rows[min(t,3)][k].d[par]
   forcings[par].append(z)
 return max(norm),norm,forcings,diff

def household_certificate(R,DT,DF,I):
 P=lambda k:V(symrad(1,R))
 Y=lambda k:V(symrad(1,R))
 h=young_iv(0,P,Y,V(IV(theta0-DT,theta0)),V(IV(phi0,phi0+DF)))
 xo,no,xr,nr=[h[s].v for s in ('xo','no','xr','nr')]
 ho=h['h'].v;u=h['u'].v;un=h['un'].v;T=h['T'].v
 p=P(0).v;p1=P(1).v;p2=P(2).v;ph=IV(phi0,phi0+DF)
 savO=y+b+T-xo-chi*no-((1-ph)+q*tax)*p*ho
 savR=y+b+T-xr-chi*nr-u*a
 h2=beta*gamma*xo/(q*un)
 assert ho.lo>a and h2.lo>a
 c2O=beta*xo/q;c2R=beta*xr/q
 eo=beta*omega*xo/q**2;er=beta*omega*xr/q**2
 out={'price':p,'user_cost':u,'adult_goods_owner':xo,'adult_goods_renter':xr,
      'fertility_owner':no,'fertility_renter':nr,'adult_space_owner':h['so'],'adult_space_renter':h['sr'],
      'saving_owner':savO,'saving_renter':savR,'owner_cap':2-ho,
      'owner_purchase_gap':alpha*xo/h['so']-u,'young_renter_cap_gap':alpha*xr/h['sr']-u,
      'generated_old_retention':ho-h2,'generated_old_estate_gap':eo-p2*h2,
      'generated_old_renter_cap_gap':gamma*c2R/a-un,
      'old_goods_owner':c2O,'old_goods_renter':c2R,'estate_owner':eo,'estate_renter':er}
 TI=taxcoef*p/(I[0]+I[1]);res=I[3]+p*I[4]+TI
 cI=res/(1+gamma+omega);hI=LL*res/u;eI=omega*cI/q
 cRI=(I[5]+TI-u*a)/(1+omega)
 out.update(initial_owner_resources=res,initial_owner_retention=I[4]-hI,initial_owner_estate_gap=eI-p1*hI,
            initial_renter_goods=cRI,initial_renter_cap_gap=gamma*cRI/a-u)
 assert all(z.lo>0 for z in out.values())
 pi=h['pi'].v
 assert Q(523,1000)<pi.lo<pi.hi<Q(525,1000),pi
 assert all(z.lo<=c<=z.hi for z,c in zip(I,[1,1,pish,aO0,1,ar0]))
 return out,pi


def interval_record(z):
    return {"lower":str(z.lo),"upper":str(z.hi)}


def run_certificate():
    if not __debug__:
        raise RuntimeError("Run this verifier without -O; its assertions are proof checks.")
    radius=Q(1,2000)
    baseline_radius=Q(1,5000)
    policy_radius=Q(1,4000)
    shock=Q(1,20000)
    inherited=[IV("0.9998","1.0002"),IV("0.9998","1.0002"),
               IV("0.52349","0.52413"),IV("0.34569","0.34904"),
               IV("0.99980","1.00021"),IV("1.61670","1.61740")]
    generated=generator_bounds(baseline_radius,shock)
    assert all(a.lo<=z.lo<=z.hi<=a.hi for a,z in zip(inherited,generated))

    # Exact reference identities and independent consistency of two derivative
    # implementations: rational implicit differentiation and interval formulas.
    assert theta0/nO0 == chi/xO0+alpha*kap/sO0
    assert theta0/nR0 == chi/xR0+alpha*kap/sR0
    for t in range(4):
        for k,row in enumerate(residual_iv(t,Q(0),Q(0),Q(0),BASEI)):
            assert row.v.lo<=0<=row.v.hi
            for key,z in row.d.items():
                if isinstance(key,tuple):
                    assert z.lo<=lentry(t,k,*key)<=z.hi

    al,la,an,als=exact_operator_checks()
    assert al<Q(2,10**9) and la<Q(3,10**9) and an<Q(401,100)
    rows=get_rows(radius,shock,shock,inherited)
    contraction,row_contractions,forcing,raw_derivative_errors=derivative_bounds(rows,al,als)
    kappa=Q(1,20)
    first_row=Q(2,125)
    C={"theta":Q(10,3),"phi":Q(451,100)}
    assert contraction<kappa and row_contractions[1]<first_row
    for parameter in C:
        assert max(z.absmax() for z in forcing[parameter])<C[parameter]
    impact_forcing={"theta":Q(839,500),"phi":Q(29,250)}
    tail_forcing={"theta":Q(237,100),"phi":Q(109,100)}
    for parameter in C:
        assert forcing[parameter][1].lo>impact_forcing[parameter]
        assert forcing[parameter][-1].lo>tail_forcing[parameter]
    impact={p:impact_forcing[p]-first_row*C[p]/(1-kappa) for p in C}
    endpoint={p:tail_forcing[p]-kappa*C[p]/(1-kappa) for p in C}
    assert impact["theta"]>Q(8,5) and impact["phi"]>Q(1,25)
    assert endpoint["theta"]>2 and endpoint["phi"]>Q(4,5)

    # Uniform residual estimates imply self-mapping on both infinite balls.
    assert C["theta"]*shock+kappa*baseline_radius<baseline_radius
    assert C["phi"]*shock+kappa*policy_radius<policy_radius
    assert baseline_radius+policy_radius<radius
    margins,ownership=household_certificate(radius,shock,shock,inherited)
    assert margins["owner_purchase_gap"].lo>Q(111,1000)
    assert margins["young_renter_cap_gap"].lo>Q(336,100)
    assert margins["saving_owner"].lo>Q(972,1000)
    assert margins["saving_renter"].lo>Q(808,1000)
    for name in ("generated_old_retention","initial_owner_retention"):
        assert margins[name].lo>Q(541,1000)
    for name in ("generated_old_estate_gap","initial_owner_estate_gap"):
        assert margins[name].lo>Q(180,1000)
    for name in ("generated_old_renter_cap_gap","initial_renter_cap_gap"):
        assert margins[name].lo>Q(757,1000)

    result={
        "result":"all exact rational checks passed",
        "scope":"infinite half-line equilibrium operator; not a finite-horizon simulation",
        "primitives":{name:str(globals()[name]) for name in
                      ("q","beta","alpha","omega","gamma","chi","kap","theta0","phi0", "b","y","stock","tax","a","nu","sigma")},
        "taste":"location = 4 log(11/10) - W_O(reference) + W_R(reference); fixed across shocks",
        "maximum_each_shock":str(shock),
        "sequence_radii":{"common":str(radius),"baseline":str(baseline_radius),"policy_about_baseline_tail":str(policy_radius)},
        "inherited_order":["Y0","O0","old_owner_share","old_owner_net_assets","old_owner_purchased_title","old_renter_assets"],
        "inherited_box":[interval_record(z) for z in inherited],
        "generated_baseline_inherited_enclosure":[interval_record(z) for z in generated],
        "exact_operator_bounds":{"I_minus_A_L":str(al),"I_minus_L_A":str(la),"A_norm":str(an),
                                 "nonlinear_contraction":str(contraction),"first_population_row":str(row_contractions[1])},
        "advertised_operator_bounds":{"I_minus_A_L":"2/1000000000","I_minus_L_A":"3/1000000000", "A_norm":"401/100", "nonlinear_contraction":"1/20", "first_population_row":"2/125"},
        "preconditioned_parameter_forcing_bounds":{p:str(C[p]) for p in C},
        "first_population_forcing_intervals":{p:interval_record(forcing[p][1]) for p in C},
        "tail_population_forcing_intervals":{p:interval_record(forcing[p][-1]) for p in C},
        "certified_first_population_derivative_lower":{p:str(impact[p]) for p in C},
        "certified_stationary_young_population_derivative_lower":{p:str(endpoint[p]) for p in C},
        "ownership_probability":interval_record(ownership),
        "household_margins":{name:interval_record(z) for name,z in margins.items()},
        "tail_bound":"sup_{s>=t} ||X_s-X_star||_infty <= (1/1000)*(1/20)**floor(t/27)",
        "arithmetic":"rational endpoints; outward dyadic rounding at 90 bits; rational Taylor remainders for log/exp; integer square-root bounds",
        "operator_representation":"24 exceptional initial block rows; thereafter a fixed +/-24 block stencil; both inverse defects verified exactly through 56 block rows, after which rows are translations",
        "preconditioner_sha256":hashlib.sha256(PATH.read_bytes()).hexdigest(),
        "verifier_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    return result


if __name__=="__main__":
    result=run_certificate()
    target=Path(__file__).with_name("finite_mixed_transition_certificate.json")
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"verification":result["result"],"maximum_each_shock":result["maximum_each_shock"],
                      "receipt":str(target)},indent=2))
~~~
