"""Experimental joint choice with independent centered logistic H and F shocks.

Q[..., tenure, attempt]; scales need not be ordered. Integration uses the
bounded residual above the exact one-shock envelope, with adaptive quadrature.
"""
import numpy as np
from numba import njit

_X8, _W8 = np.polynomial.legendre.leggauss(8)
_X16, _W16 = np.polynomial.legendre.leggauss(16)

@njit(cache=True)
def logistic(x):
    z = np.exp(-abs(x))
    return 1/(1+z) if x >= 0 else z/(1+z)

@njit(cache=True)
def softplus(x):
    return max(x, 0.) + np.log1p(np.exp(-abs(x)))

@njit(cache=True)
def _node(q, sh, sf, t):
    h = sh*t
    d0 = int(q[1,0]+h > q[0,0])
    d1 = int(q[1,1]+h > q[0,1])
    m0, m1 = q[d0,0]+d0*h, q[d1,1]+d1*h
    out = np.zeros(5)
    if not np.isfinite(m0):
        pf = 1.
    elif not np.isfinite(m1):
        pf = 0.
    else:
        g = (m1-m0)/sf
        pf = logistic(g)
        out[0] = sf*np.log1p(np.exp(-abs(g)))
    out[1+2*d0] = 1-pf
    out[2+2*d1] = pf
    return out * logistic(t)*logistic(-t)

@njit(cache=True)
def _rule(q, sh, sf, lo, hi, nodes, weights):
    out = np.zeros(5)
    for k in range(len(nodes)):
        out += weights[k]*_node(q, sh, sf, (lo+hi)/2+(hi-lo)*nodes[k]/2)
    return out*(hi-lo)/2

@njit(cache=True)
def _integrate(q, sh, sf, lo, hi, tol, depth):
    # Explicit stack avoids cached recursive JIT calls (unsafe on some runtimes).
    left=np.empty(64);right=np.empty(64);budget=np.empty(64);level=np.empty(64,dtype=np.int64)
    left[0]=lo;right[0]=hi;budget[0]=tol;level[0]=depth
    size=1;out=np.zeros(5);evaluations=0
    while size:
        size-=1;lo=left[size];hi=right[size];tol=budget[size];depth=level[size]
        a=_rule(q,sh,sf,lo,hi,_X8,_W8)
        b=_rule(q,sh,sf,lo,hi,_X16,_W16)
        evaluations+=1
        if np.max(np.abs(a-b))<=tol:
            out+=b
        else:
            if depth==0 or evaluations>10000:
                raise RuntimeError('Two-shock quadrature did not converge')
            mid=(lo+hi)/2
            left[size]=lo;right[size]=mid;budget[size]=tol/2;level[size]=depth-1;size+=1
            left[size]=mid;right[size]=hi;budget[size]=tol/2;level[size]=depth-1;size+=1
    return out

@njit(cache=True)
def _one(raw, sh, sf, tol):
    center = np.max(raw)
    if not np.isfinite(center):
        return -np.inf, np.zeros((2,2))
    q = raw-center
    # Exact reductions avoid quadrature in states without a choice margin.
    rows = np.array([np.any(np.isfinite(q[0])),np.any(np.isfinite(q[1]))])
    cols = np.array([np.any(np.isfinite(q[:,0])),np.any(np.isfinite(q[:,1]))])
    p = np.zeros((2,2))
    if not (rows[0] and rows[1]):
        d = int(rows[1])
        if not (cols[0] and cols[1]):
            a = int(cols[1]); p[d,a]=1.
            return center+q[d,a],p
        pf=logistic((q[d,1]-q[d,0])/sf)
        p[d,0]=1-pf;p[d,1]=pf
        return center+max(q[d,0],q[d,1])+sf*np.log1p(np.exp(-abs(q[d,1]-q[d,0])/sf)),p
    if not (cols[0] and cols[1]):
        a=int(cols[1]);ph=logistic((q[1,a]-q[0,a])/sh)
        p[0,a]=1-ph;p[1,a]=ph
        return center+max(q[0,a],q[1,a])+sh*np.log1p(np.exp(-abs(q[1,a]-q[0,a])/sh)),p
    c0,c1=np.max(q[0]),np.max(q[1])
    base=max(c0,c1)+sh*np.log1p(np.exp(-abs(c1-c0)/sh))
    # Omitted residual is bounded by 2*sf*log(2)*logistic(-T).
    T=max(36.,np.log(max(1.,sf)/tol)+5.)
    cuts=np.empty(80); cuts[0]=-T;cuts[1]=T;cuts[2]=0.;n=3
    for a in range(2):
        if np.isfinite(q[0,a]) and np.isfinite(q[1,a]):
            t=(q[0,a]-q[1,a])/sh
            if -T<t<T: cuts[n]=t;n+=1
    cuts[:n].sort()
    original=cuts[:n].copy()
    for k in range(len(original)-1):
        lo,hi=original[k],original[k+1];h=sh*(lo+hi)/2
        d0=int(q[1,0]+h>q[0,0]);d1=int(q[1,1]+h>q[0,1])
        if d0!=d1:
            t=(q[d0,0]-q[d1,1])/(sh*(d1-d0))
            # Resolve thin F boundary layers even if their center is outside.
            for width in (-36.,-16.,-4.,-1.,0.,1.,4.,16.,36.):
                point=t+width*sf/sh
                if lo<point<hi:cuts[n]=point;n+=1
    cuts[:n].sort()
    out=np.zeros(5)
    for k in range(n-1):
        out+=_integrate(q,sh,sf,cuts[k],cuts[k+1],tol/(n-1),20)
    p=out[1:].reshape((2,2))
    # Tail probabilities are the derivative of the exact base envelope.
    threshold=(c0-c1)/sh
    for lo,hi in ((-np.inf,-T),(T,np.inf)):
        own=logistic(-max(lo,threshold))-logistic(-hi) if hi>max(lo,threshold) else 0.
        rent=logistic(min(hi,threshold))-logistic(lo) if min(hi,threshold)>lo else 0.
        p[0,int(q[0,1]>q[0,0])]+=rent
        p[1,int(q[1,1]>q[1,0])]+=own
    return center+base+out[0],p

@njit(cache=True)
def _batch(q, sh, sf, tol):
    v=np.empty(len(q));p=np.empty_like(q)
    for i in range(len(q)):
        v[i],p[i]=_one(q[i],sh,sf[i],tol)
    return v,p

def choose(plans, housing_scale, fertility_scale, tolerance=2e-11):
    q=np.asarray(plans,dtype=float)
    if q.shape[-2:]!=(2,2) or np.any(np.isnan(q)) or np.any(np.isposinf(q)):
        raise ValueError('Expected finite or negative-infinite four-plan values')
    sf=np.broadcast_to(fertility_scale,q.shape[:-2]).astype(float)
    if not np.isfinite(housing_scale) or housing_scale<=0 or np.any(~np.isfinite(sf)) or np.any(sf<=0):
        raise ValueError('Both independent shock scales must be positive')
    if not np.isfinite(tolerance) or tolerance<=0:
        raise ValueError('Positive quadrature tolerance required')
    v,p=_batch(np.ascontiguousarray(q.reshape(-1,2,2)),float(housing_scale),sf.ravel(),float(tolerance))
    occupied=np.any(np.isfinite(q.reshape(-1,4)),axis=1)
    if np.any(~np.isfinite(p)) or np.any(~np.isfinite(v[occupied])):
        raise RuntimeError('Nonfinite two-shock result')
    if np.any(abs(p.sum(axis=(1,2))[occupied]-1)>10*tolerance) or np.any(p<0) or np.any(p>1+10*tolerance):
        raise RuntimeError('Two-shock probability accounting failed')
    return v.reshape(q.shape[:-2]),p.reshape(q.shape)
