"""Focused independent CES mathematics and objective/policy integration checks."""
from pathlib import Path
import sys,unittest
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from refactor_lab.engine import kernels as k

class CESChecks(unittest.TestCase):
    def test_renter_foc_cap_accounting(self):
        a=.733;eta=.487;rho=(eta-1)/eta
        for ec,eh in ((1.,1.),(1.5,1.8)):
            for rent in (.02,.2,1.,2.):
                for S in (.01,.5,10.):
                    u,c,h=k.ces_renter_allocation(S,rent,8.,a,-1.,eta,ec,eh)
                    self.assertAlmostEqual(c+rent*h,S,places=11)
                    muc=k.ces_marginal_c(c,h,a,-1.,eta,ec,eh)
                    Q=a*(c/ec)**rho+(1-a)*(h/eh)**rho
                    muh=(1-a)/eh*(h/eh)**(rho-1)*Q**((-1-rho)/rho)
                    if h<8.-1e-10: self.assertAlmostEqual(muh/muc,rent,places=10)
                    else: self.assertGreaterEqual(muh/muc+1e-10,rent)
                    price=(a**eta*ec**(1-eta)+(1-a)**eta*(rent*eh)**(1-eta))**(1/(1-eta))
                    if h<8.-1e-10:self.assertAlmostEqual(u,-price/S,places=10)
    def test_cd_limit_lambda_zero(self):
        a=.733;c=.8;h=3.;e=1.6
        cd=(c**a*h**(1-a))**-1/-1
        self.assertAlmostEqual(k.ces_flow(c,h,a,-1.,1.,1.,1.),cd,places=12)
        self.assertAlmostEqual(k.ces_flow(c,h,a,-1.,.487,e,e),e*k.ces_flow(c,h,a,-1.,.487,1.,1.),places=12)
    def test_owner_marginal(self):
        c=2.;h=6.;dc=1e-5
        derivative=(k.ces_flow(c+dc,h,.733,-1.,.487,1.3,1.6)-k.ces_flow(c-dc,h,.733,-1.,.487,1.3,1.6))/(2*dc)
        self.assertAlmostEqual(derivative,k.ces_marginal_c(c,h,.733,-1.,.487,1.3,1.6),places=8)
    def test_exhaustive_against_independent_interval_optimizer(self):
        bg=np.array([0.,.05,.2,.6,1.2,2.]);V=np.array([-.9,-.8,-.65,-.4,-.3,-.2])
        for owner in (False,True):
            for resources in (.4,2.,10.):
                oc=.1 if owner else 0.;hi=resources-oc-1e-6;lo=0.
                bp,value=k.exhaustive_saving_scalar(lo,hi,resources,V,bg,.3,0.,0.,.12,3.,.733,-1.,.97,1.,oc,1.,owner,.487,1.4,1.8,5.)
                def f(x):
                    if owner:u=k.ces_flow(resources-oc-x,5.,.733,-1.,.487,1.4,1.8)
                    else:u=k.ces_renter_allocation(resources-x,.3,3.,.733,-1.,.487,1.4,1.8)[0]
                    return u+.12+.97*np.interp(x,bg,V)
                points=sorted(set([lo,hi]+[float(x) for x in bg if lo<x<hi]))
                best=max(f(x) for x in points)
                for left,right in zip(points[:-1],points[1:]):
                    opt=minimize_scalar(lambda x:-f(x),bounds=(left,right),method='bounded',options={'xatol':1e-12})
                    best=max(best,-opt.fun)
                self.assertAlmostEqual(value,best,places=8)
                self.assertAlmostEqual(value,f(bp),places=11)
    def test_full_block_returns_objective_allocations(self):
        bg=np.array([0.,.1,.5]);V=np.zeros((3,1));zero=np.zeros(1);one=np.ones(1);ec=np.array([1.3]);eh=np.array([1.6]);R=np.array([.3,1.,2.]);psi=np.array([.15]);zV=np.zeros_like(V)
        v,bp,c,h=k.full_renter_block_kernel(R,R,V,zV,0,bg,zero,zero,psi,zero,np.array([.733]),one,.3,8.,1e-3,0.,0.,.733,-1.,.97,0.,0.,.381966,.618034,1e-3,1,zero,0,0,0.,0.,6.,True,None,0.,.487,ec,eh)
        np.testing.assert_allclose(c+.3*h+bp,R[:,None],rtol=0,atol=1e-12)
        for i in range(3):self.assertAlmostEqual(v[i,0],k.ces_flow(c[i,0],h[i,0],.733,-1.,.487,ec[0],eh[0])+.15,places=11)
        v,bp,c=k.full_owner_block_kernel(R,R,V,zV,0,bg,zero,zero,psi,zero,np.array([.733]),one,zero,.1,4.,1.,1.2,1e-3,.733,-1.,.97,0.,0.,.381966,.618034,1e-3,0,1,zero,0,0,0,0.,True,False,-np.inf,.487,ec,eh)
        np.testing.assert_allclose(c+bp+.1,R[:,None],rtol=0,atol=1e-12)
        for i in range(3):self.assertAlmostEqual(v[i,0],k.ces_flow(c[i,0],4.8,.733,-1.,.487,ec[0],eh[0])+.15,places=11)

if __name__=='__main__':unittest.main()
