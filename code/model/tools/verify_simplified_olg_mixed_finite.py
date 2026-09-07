#!/usr/bin/env python3
"""Verify the explicit infinite mixed-tenure transition neighborhood.

Run from the project root: python3 code/model/tools/verify_simplified_olg_mixed_finite.py
The certificate uses rational intervals, not a finite-horizon simulation.
"""
import sys, math, json
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sympy as s
import verify_simplified_olg_mixed_transition as old
BITS=110; DEN=2**BITS

def frac(x):
    if isinstance(x,F):return x
    if isinstance(x,s.Rational):return F(int(x.p),int(x.q))
    return F(x)
class I:
    def __init__(self,lo,hi=None):
        lo=frac(lo);hi=lo if hi is None else frac(hi)
        self.lo=F((lo*DEN).__floor__(),DEN);self.hi=F((hi*DEN).__ceil__(),DEN)
        assert self.lo<=self.hi
    def __add__(self,o):o=iv(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+-iv(o)
    def __rsub__(self,o):return iv(o)+-self
    def __mul__(self,o):
        o=iv(o);p=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)];return I(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=iv(o);assert o.lo>0 or o.hi<0,('zero division',o.pair())
        return self*I(1/o.hi,1/o.lo)
    def __rtruediv__(self,o):return iv(o)/self
    def __pow__(self,n):
        assert n>=0;out=I(1)
        for _ in range(n):out=out*self
        return out
    def sqrt(self):
        assert self.lo>0
        a=math.isqrt((self.lo.numerator*DEN*DEN)//self.lo.denominator)
        b=math.isqrt((self.hi.numerator*DEN*DEN)//self.hi.denominator)+1
        return I(F(a,DEN),F(b,DEN))
    def log(self):
        assert self.lo>0
        z=(self-1)/(self+1);b=z.mag();assert b<F(3,5)
        n=100 if b>F(1,8) else 24
        out=I(0);power=z
        for j in range(n):out+=power/F(2*j+1);power=power*z*z
        tail=2*b**(2*n+1)/(F(2*n+1)*(1-b*b))
        return 2*out+I(-tail,tail)
    def exp(self):
        b=self.mag();assert b<F(1,4),('exp too wide',float(b))
        out=I(1);term=I(1)
        for j in range(1,25):term=term*self/j;out+=term
        tail=2*b**25/math.factorial(25)
        return out+I(-tail,tail)
    def mag(self):return max(abs(self.lo),abs(self.hi))
    def pair(self):return [float(self.lo),float(self.hi)]

def iv(x):return x if isinstance(x,I) else I(x)
NVAR=10
class A:
    def __init__(self,x,der=None):self.x=iv(x);self.d=[I(0)]*NVAR if der is None else der
    @staticmethod
    def variable(x,j):
        d=[I(0)]*NVAR;d[j]=I(1);return A(x,d)
    def __add__(self,o):o=ad(o);return A(self.x+o.x,[a+b for a,b in zip(self.d,o.d)])
    __radd__=__add__
    def __neg__(self):return A(-self.x,[-a for a in self.d])
    def __sub__(self,o):return self+-ad(o)
    def __rsub__(self,o):return ad(o)+-self
    def __mul__(self,o):o=ad(o);return A(self.x*o.x,[a*o.x+self.x*b for a,b in zip(self.d,o.d)])
    __rmul__=__mul__
    def __truediv__(self,o):
        o=ad(o);return A(self.x/o.x,[(a*o.x-self.x*b)/(o.x*o.x) for a,b in zip(self.d,o.d)])
    def __rtruediv__(self,o):return ad(o)/self
    def __pow__(self,n):
        out=A(1)
        for _ in range(n):out=out*self
        return out
    def sqrt(self):v=self.x.sqrt();return A(v,[a/(2*v) for a in self.d])
    def log(self):return A(self.x.log(),[a/self.x for a in self.d])
    def exp(self):v=self.x.exp();return A(v,[a*v for a in self.d])

def ad(x):return x if isinstance(x,A) else A(x)
def matmul(a,b):return [[sum((iv(x)*iv(y) for x,y in zip(row,col)),I(0)) for col in zip(*b)] for row in a]
def matadd(a,b):return [[x+y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]
def negmat(a):return [[-x for x in row] for row in a]
def norm(a):return max(sum(iv(x).mag() for x in row) for row in a)
def exactmat(x):return [[frac(v) for v in row] for row in x.tolist()]
def eye(n):return [[F(int(i==j)) for j in range(n)] for i in range(n)]
def inv2(a):
    det=a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/det,-a[0][1]/det],[-a[1][0]/det,a[0][0]/det]]
def colnormbounds(a,radii):return [sum(iv(v).mag()*r for v,r in zip(row,radii)) for row in a]

q=F(1,2);beta=alpha=omega=F(2,5);gamma=F(3,10);chi=F(3,20);kap=F(1,2)
phi0=F(4,5);theta0=F(141,400);tau=F(467,9250);a=F(1,4);b=F(1,5)
y=F(3521349281,1677802000);Hbar=F(68104,68019);sigma=F(4)
rhoO=1+beta*(1+gamma+omega);rhoR=1+beta*(1+omega);L=gamma/(1+gamma+omega)
pi0=F(11,21);u0=F(9717,18500);T0=F(3975571,314587875)
Z0=[F(1),u0,F(1),F(1),beta*gamma/q*pi0,a*(1-pi0)]
logratio=I(F(3,10)).log()

def evaluate(zrad,vrad,prad):
    vals=[I(c-r,c+r) for c,r in zip(Z0,zrad)]+[I(u0-vrad[0],u0+vrad[0]),I(1-vrad[1],1+vrad[1]),I(phi0-prad[0],phi0+prad[0]),I(theta0-prad[1],theta0+prad[1])]
    P,u,Y,O,M,R,un,Yn,phi,theta=[A.variable(v,j) for j,v in enumerate(vals)]
    Pn=((1+q*tau)*P-u)/q;Pnn=((1+q*tau)*Pn-un)/q
    T=q*tau*P*Hbar/(Y+O);Tn=q*tau*Pn*Hbar/(Yn+Y);w=y+b+T+q*Tn
    hh=[]
    for own in (True,False):
        h=b/((1-phi)*P) if own else A(a)
        rho=rhoO if own else rhoR
        avail=w-u*h-(0 if own else q*un*a)
        mid=kap*avail*(theta+alpha)+chi*h*(theta+rho)
        const=theta*avail*h;quad=chi*kap*(theta+rho+alpha)
        n=2*const/(mid+(mid*mid-4*quad*const).sqrt())
        x=(avail-chi*n)/rho;supp=h-kap*n
        c2=beta*x/q;h2=gamma*c2/un if own else A(a);estate=omega*c2/q
        payment=((1-phi)+q*tau)*P*h if own else u*h
        saving=y+b+T-x-chi*n-payment
        assets=saving/q-(phi*P*h/q if own else 0)
        hh.append(dict(x=x,h=h,n=n,s=supp,c2=c2,h2=h2,estate=estate,saving=saving,assets=assets))
    ho,hr=hh
    # Exact anchor centering fixes BOTH logistic taste parameters, avoiding a
    # rounded irrational taste intercept. These are original W^R-W^O changes.
    dw=rhoR*(hr['x']/F(99,74)).log()-rhoO*ho['x'].log()
    dw+=alpha*((hr['s']/F(11,80)).log()-(ho['s']/F(5,8)).log())
    dw+=theta*((hr['n']/F(9,40)).log()-(ho['n']/F(3,4)).log())
    dw+=(theta-theta0)*logratio+beta*gamma*(un/u0).log()
    pi=1/(1+F(10,11)*(dw/sigma).exp())
    housing=pi*ho['h']+(1-pi)*a;fert=pi*ho['n']+(1-pi)*hr['n']
    FF=[Y*housing+M/u+R-Hbar,Yn-2*Y*fert]
    GG=[Pn,un,Yn,Y,beta*gamma/q*Y*pi*ho['x'],a*Y*(1-pi)]
    return FF,GG,dict(owner=ho,renter=hr,pi=pi,P=P,u=u,Pn=Pn,Pnn=Pnn,T=T,Tn=Tn,Y=Y,Yn=Yn,un=un,phi=phi,theta=theta)

# Numerical vectors only choose rational coordinates. All asserted bounds below
# use exact rational J, Q, Q^-1 and outward rational interval arithmetic.
data=old.exact_linearization(s.Rational(4));J=data['J'];Fv0=data['implicit'];C0=exactmat(Fv0.inv())
# Frozen rational coordinates. Their numerical origin is irrelevant to every
# proof assertion: all block and boundary calculations below are exact.
Q=s.Matrix([[s.Rational(v) for v in row] for row in [
 ['3357591253/1000000000000','0','-202851141137/500000000000','-338025666639/1000000000000','-301293427/125000000000','95323191951/1000000000000'],
 ['1718491519/500000000000','0','-219496015737/500000000000','-238604624561/1000000000000','-1/10','-1/10'],
 ['31955691/10000000000','0','-7786209013/40000000000','-184999390761/500000000000','712867613/1000000000000','-5129959667/500000000000'],
 ['1','1/1000000000000','-1','0','-440449/50000000000','-2473096997/1000000000000'],
 ['-7397734087/100000000000','-525243243243/1000000000000','-2817390331/31250000000','-18101161249/500000000000','47415113/500000000000','606517071/200000000000'],
 ['71168931627/500000000000','1','-134659311767/1000000000000','1100557953/125000000000','-9788597/50000000000','-773123759/1000000000000']]])
Qinv=Q.inv();JJ=Qinv*J*Q;DU=JJ[4:6,4:6];DUinv=DU.inv();QM=exactmat(Q);QI=exactmat(Qinv);DUI=exactmat(DUinv)
B0=s.zeros(4,6);B0[0,2]=1;B0[1,3]=1;B0[2,5]=1;B0[3,4]=1
B0[3,0]=-s.Rational(L)*(s.Rational(pi0)+s.Rational(pi0*q*tau*Hbar/2))
BS0=B0*Q[:,:4];Cbound=-(BS0.inv()*B0*Q[:,4:6])
print('chart exact boundary norm',float(norm(exactmat(Cbound))),flush=True)


def certificate(rad,eps,vrad):
    zr=[sum(abs(x) for x in row)*rad for row in QM]
    FF,GG,hh=evaluate(zr,[vrad,vrad],[eps,eps])
    DF=[f.d for f in FF];DG=[g.d for g in GG]
    Fz=[r[:6] for r in DF];Fv=[r[6:8] for r in DF];Fp=[r[8:] for r in DF]
    Gz=[r[:6] for r in DG];Gv=[r[6:8] for r in DG];Gp=[r[8:] for r in DG]
    inner=matadd(eye(2),negmat(matmul(C0,Fv)));innernorm=norm(inner)
    innerforcing=max(colnormbounds(matmul(matmul(C0,Fz),QM),[rad]*6))+max(colnormbounds(matmul(C0,Fp),[eps,eps]))
    print('inner',float(innernorm),'forcing/radius',float(innerforcing/vrad),flush=True)
    assert innernorm<1 and innerforcing+innernorm*vrad<vrad
    IV=inv2(Fv)
    gZ=matadd(Gz,negmat(matmul(matmul(Gv,IV),Fz)))
    gp=matadd(Gp,negmat(matmul(matmul(Gv,IV),Fp)))
    chart=matmul(matmul(QI,gZ),QM)
    upper=chart[:4]
    low=matmul(DUI,chart[4:])
    for i in range(2):
        for j in range(6):low[i][j]=I(int(j==4+i))-low[i][j]
    outer=max(norm(upper),norm(low)+norm(DUI),norm(exactmat(Cbound)))
    chartp=matmul(QI,gp)
    force=max(max(colnormbounds(chartp[:4],[eps,eps])),max(colnormbounds(matmul(DUI,chartp[4:]),[eps,eps])))
    print('outer',float(outer),'force/radius',float(force/rad),flush=True)
    assert outer<1 and force+outer*rad<rad
    return dict(rad=rad,eps=eps,zrad=zr,vrad=vrad,outer=outer,inner=innernorm,forcing=force,hh=hh,chart=chart,upper=upper,low=low,gp=gp,chartp=chartp)



def private_checks(hh):
    own,rent=hh['owner'],hh['renter']
    margins={}
    for name,h in [('owner',own),('renter',rent)]:
        for key in ('x','s','n','saving','c2','estate'):
            margins[name+'_'+key]=h[key].x.lo
        margins[name+'_young_cap']= (alpha*h['x']/h['s']-hh['u']).x.lo
    margins['owner_physical_cap']=(2-own['h']).x.lo
    margins['owner_retention']=(own['h']-own['h2']).x.lo
    margins['owner_estate_slack']=(own['estate']-hh['Pnn']*own['h2']).x.lo
    margins['owner_h2_above_renter']=(own['h2']-a).x.lo
    margins['renter_old_cap']=(gamma*rent['c2']/a-hh['un']).x.lo
    margins['owner_h_above_renter']=(own['h']-a).x.lo
    for key in ('P','Pn','Pnn','u','un'):margins[key]=hh[key].x.lo
    assert all(v>0 for v in margins.values()),{k:float(v) for k,v in margins.items() if v<=0}
    assert hh['pi'].x.lo>F(52,100) and hh['pi'].x.hi<F(53,100)
    return {k:float(v) for k,v in margins.items()}


def generated_old(base):
    zr=base['zrad']
    # These are actual adjacent equilibrium states, both in the chart ball;
    # use their tighter radii, not the auxiliary forward-solve rectangle.
    FF,GG,hh=evaluate(zr,[zr[1],zr[2]],[0,base['eps']])
    U=hh['Y'].x;Y=hh['Yn'].x;pi=hh['pi'].x
    assets=hh['owner']['assets'].x;title=hh['owner']['h'].x
    return dict(Y=Y,O=U,B=U*pi,A=U*pi*assets,H=U*pi*title,
                aO=assets,hO=title,aR=hh['renter']['assets'].x)


def initial_old():
    saveO=y+b+T0-F(1)-chi*F(3,4)-((1-phi0)+q*tau)
    saveR=y+b+T0-F(99,74)-chi*F(9,40)-u0*a
    aO=(saveO-phi0)/q;aR=saveR/q
    assert L*(pi0*aO+pi0+T0*pi0)==Z0[4]
    return dict(Y=I(1),O=I(1),B=I(pi0),A=I(pi0*aO),H=I(pi0),
                aO=I(aO),hO=I(1),aR=I(aR))


def actual_boundary(states,pol):
    Y,O,B,AA,H=[states[k] for k in ('Y','O','B','A','H')]
    coef=-L*(H+q*tau*Hbar*B/(Y+O))
    BM=[[I(x) for x in row] for row in exactmat(B0)];BM[3][0]=coef
    BS=matmul(BM,[row[:4] for row in QM]);BU=matmul(BM,[row[4:] for row in QM])
    BSinv=exactmat(BS0.inv())
    dev=matadd(eye(4),negmat(matmul(BSinv,BS)))
    err=norm(dev);assert err<1
    # Neumann inverse bound, preserving the exact boundary row cancellations.
    coupling=norm(matmul(BSinv,BU))/(1-err)
    rhs=[Y-1,O-1,a*(O-B)-Z0[5],L*(AA+H+q*tau*Hbar*B/(Y+O))-Z0[4]]
    forcing=max(colnormbounds(matmul(BSinv,[[v] for v in rhs]),[F(1)]))/(1-err)
    assert coupling<1 and forcing+coupling*pol['rad']<pol['rad']
    zr=pol['zrad'];P=I(1-zr[0],1+zr[0]);u=I(u0-zr[1],u0+zr[1])
    Pn=((1+q*tau)*P-u)/q;T=q*tau*P*Hbar/(Y+O)
    co=(states['aO']+P*states['hO']+T)/(1+gamma+omega)
    ho=gamma*co/u;eo=omega*co/q
    cr=(states['aR']+T-u*a)/(1+omega)
    gaps=[co,states['hO']-ho,eo-Pn*ho,cr,gamma*cr/a-u]
    assert all(v.lo>0 for v in gaps),[v.pair() for v in gaps]
    return dict(coupling=coupling,forcing=forcing,err=err,
                old_margins=[float(v.lo) for v in gaps],BM=BM,BS=BS,BU=BU,BSinv=BSinv)



def derivative_certificate(regime,boundary,columns):
    zr=regime['zrad']
    FF,GG,hh=evaluate(zr,[zr[1],zr[2]],[regime['eps'],regime['eps']])
    DF=[f.d for f in FF];DG=[g.d for g in GG]
    Fz=[r[:6] for r in DF];Fv=[r[6:8] for r in DF];Fp=[r[8:] for r in DF]
    Gz=[r[:6] for r in DG];Gv=[r[6:8] for r in DG];Gp=[r[8:] for r in DG]
    IV=inv2(Fv)
    gz=matadd(Gz,negmat(matmul(matmul(Gv,IV),Fz)))
    gp=matadd(Gp,negmat(matmul(matmul(Gv,IV),Fp)))
    chart=matmul(matmul(QI,gz),QM);chartp=matmul(QI,gp)
    upper=chart[:4];low=matmul(DUI,chart[4:])
    for i in range(2):
        for j in range(6):low[i][j]=I(int(j==4+i))-low[i][j]
    fp=chartp[:4]+negmat(matmul(DUI,chartp[4:]))
    # Validated inverse for the varying actual-old boundary; E has norm <1.
    inv0=boundary['BSinv'];EE=matadd(eye(4),negmat(matmul(inv0,boundary['BS'])))
    rr=norm(EE);power=eye(4);approx=[[I(0) for _ in range(4)] for _ in range(4)]
    for j in range(4):approx=matadd(approx,matmul(power,inv0));power=matmul(power,EE)
    tail=rr**4/(1-rr)*norm(inv0)
    approx=[[v+I(-tail,tail) for v in row] for row in approx]
    cb=negmat(matmul(approx,boundary['BU']))
    kk=max(norm(upper),norm(low)+norm(DUI),norm(cb));assert kk<1
    print('ACTUAL derivative contraction',float(kk),flush=True)
    outputs={}
    for j in columns:
        name=['phi','theta'][j]
        forcing=[row[j] for row in fp];bd=max(v.mag() for v in forcing)/(1-kk)
        unknown=[I(-bd,bd) for _ in range(6)];nn=32
        seq=[list(unknown) for _ in range(nn+2)]
        def dot(row,vec):return sum((iv(x)*v for x,v in zip(row,vec)),I(0))
        for it in range(64):
            new=[]
            for t in range(nn+1):
                ss=[dot(row,seq[0][4:]) for row in cb] if t==0 else [dot(row,seq[t-1])+forcing[i] for i,row in enumerate(upper)]
                uu=[dot(row,seq[t])+dot(DUI[i],seq[t+1][4:])+forcing[4+i] for i,row in enumerate(low)]
                new.append(ss+uu)
            new.append(list(unknown));seq=new
        impact=dot(QM[2],seq[1]) # inherited Y0 is held fixed.
        # The constant-state operator has the same interior derivative rows.
        zz=list(unknown)
        for _ in range(100):
            zz=[dot(row,zz)+forcing[i] for i,row in enumerate(upper)]+[
                dot(row,zz)+dot(DUI[i],zz[4:])+forcing[4+i] for i,row in enumerate(low)]
        terminal=dot([QM[2][k]+QM[3][k] for k in range(6)],zz)
        assert impact.lo>0 and terminal.lo>0, (name,impact.pair(),terminal.pair())
        if name=='theta':
            assert impact.lo>F(1694,1000) and terminal.lo>F(4797,1000)
        else:
            assert impact.lo>F(143,1000) and terminal.lo>F(222,100)
        outputs[name]={'impact_dY1':impact.pair(),'terminal_dN':terminal.pair(),
          'impact_exact':[str(impact.lo),str(impact.hi)],'terminal_exact':[str(terminal.lo),str(terminal.hi)],
          'norm_bound':float(bd),'norm_bound_exact':str(bd)}
        print(name,outputs[name],flush=True)
    return outputs


def exact_anchor_check():
    """Exact original stationary identities; interval containment alone is not equality."""
    xO, nO, hO = F(1), F(3,4), F(1)
    xR, nR = F(99,74), F(9,40)
    assert u0 == 1+q*tau-q
    assert 2*T0 == q*tau*Hbar
    for own,x,n,h,rho in ((True,xO,nO,hO,rhoO),(False,xR,nR,a,rhoR)):
        assert theta0/n == chi/x+alpha*kap/(h-kap*n)
        available = y+b+(1+q)*T0-u0*h-(0 if own else q*u0*a)
        assert rho*x+chi*n == available
        payment = ((1-phi0)+q*tau)*h if own else u0*h
        saving = y+b+T0-x-chi*n-payment
        assets = saving/q-(phi0*h/q if own else 0)
        c2, estate = beta*x/q, beta*omega*x/q**2
        h2 = beta*gamma*x/(q*u0) if own else a
        assert c2+q*estate+u0*h2 == assets+(h if own else 0)+T0
        assert saving>0 and alpha*x/(h-kap*n)>u0
        if own:
            assert (1-phi0)*h == b and 2>h>h2 and estate>h2
        else:
            assert gamma*c2/a>u0
    h2O=beta*gamma*xO/(q*u0)
    assert pi0*nO+(1-pi0)*nR == F(1,2)
    assert pi0*(hO+h2O)+(1-pi0)*2*a == Hbar
    assert Z0[4] == pi0*u0*h2O and Z0[5] == (1-pi0)*a
    return "Original stationary budgets, FOCs, constraints, housing, rebates and replacement hold exactly."


def point_check():
    FF,GG,hh=evaluate([F(0)]*6,[F(0)]*2,[F(0)]*2)
    assert all(f.x.lo<=0<=f.x.hi for f in FF)
    assert all(g.x.lo<=z<=g.x.hi for g,z in zip(GG,Z0))
    fv=[r.d[6:8] for r in FF];fz=[r.d[:6] for r in FF];fp=[r.d[8:] for r in FF]
    gv=[r.d[6:8] for r in GG];gz=[r.d[:6] for r in GG];gp=[r.d[8:] for r in GG]
    jac=matadd(gz,negmat(matmul(matmul(gv,inv2(fv)),fz)))
    forcing=matadd(gp,negmat(matmul(matmul(gv,inv2(fv)),fp)))
    for ii in range(2):
        for jj in range(2):
            exact=frac(Fv0[ii,jj]);assert fv[ii][jj].lo<=exact<=fv[ii][jj].hi
    for ii in range(6):
        for jj in range(6):
            exact=frac(J[ii,jj]);assert jac[ii][jj].lo<=exact<=jac[ii][jj].hi
            assert jac[ii][jj].hi-jac[ii][jj].lo<F(1,10**23)
        exact=frac(data['Q'][ii]);assert forcing[ii][0].lo<=exact<=forcing[ii][0].hi
    assert hh['pi'].x.lo<=pi0<=hh['pi'].x.hi
    print('PASS original exact F_v, J, phi forcing, probabilities and anchor equations',flush=True)


def run():
    anchor_result=exact_anchor_check()
    point_check()
    # Fixed, reproducible radii; no numerical search or equilibrium simulation.
    base=certificate(F(1,10**8),F(1,10**11),F(1,10**6))
    pol=certificate(F(1,10**5),F(1,10**8),F(1,10**3))
    marginsB=private_checks(base['hh']);marginsP=private_checks(pol['hh'])
    initial=actual_boundary(initial_old(),base)
    inherited=actual_boundary(generated_old(base),pol)
    print('ACTUAL OLD uniform later boundary',float(inherited['coupling']),float(inherited['forcing']/pol['rad']),inherited['old_margins'],flush=True)
    report={'exact_anchor_check':anchor_result,'baseline':{k:float(base[k]) for k in ('rad','eps','outer','inner','forcing')},
            'policy':{k:float(pol[k]) for k in ('rad','eps','outer','inner','forcing')},
            'policy_inherited':{k:float(inherited[k]) for k in ('coupling','forcing','err')},
            'initial_old_margins':initial['old_margins'],'later_old_margins':inherited['old_margins'],
            'baseline_private_margins':marginsB,'policy_private_margins':marginsP,
            'owner_probability':pol['hh']['pi'].x.pair(),
            'Q_exact':[[str(v) for v in row] for row in QM],
            'baseline_bounds_exact':{k:str(base[k]) for k in ('rad','eps','outer','inner','forcing')},
            'policy_bounds_exact':{k:str(pol[k]) for k in ('rad','eps','outer','inner','forcing')}}
    report['baseline_derivatives']=derivative_certificate(base,initial,[1])
    report['policy_derivatives']=derivative_certificate(pol,inherited,[0])
    print(json.dumps(report,indent=2),flush=True)
    return report

if __name__=='__main__':run()
