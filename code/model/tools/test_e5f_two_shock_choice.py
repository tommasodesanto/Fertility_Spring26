"""Independent adaptive integration and choice-envelope checks."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):os.environ.setdefault(key,"1")
import sys, unittest, itertools, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'code/model'))
from intergen_eqscale_seq_optimized.two_shock_choice import choose


def reference(q,sh,sf):
    # Direct conditional softplus, integrated by QUADPACK, independent of
    # the implementation's bounded-residual/envelope decomposition.
    center=np.max(q)
    if not np.isfinite(center):return -np.inf,np.zeros((2,2))
    q=q-center
    points=[-45.,0.,45.]
    for a in range(2):
        if np.all(np.isfinite(q[:,a])):
            points.append(float(np.clip((q[0,a]-q[1,a])/sh,-45,45)))
    for d0,d1 in itertools.product(range(2),repeat=2):
        if d0!=d1 and np.isfinite(q[d0,0]) and np.isfinite(q[d1,1]):
            root=(q[d0,0]-q[d1,1])/(sh*(d1-d0))
            points.extend(np.clip(root+np.array([-40,-10,-2,0,2,10,40])*sf/sh,-45,45))
    points=sorted(set(points))
    def node(t,k):
        u=q+np.array([0.,sh*t])[:,None]
        ds=np.argmax(u,axis=0);m=np.max(u,axis=0)
        if not np.isfinite(m[0]):pf=1.;v=m[1]
        elif not np.isfinite(m[1]):pf=0.;v=m[0]
        else:pf=expit((m[1]-m[0])/sf);v=max(m)+sf*np.log1p(np.exp(-abs(m[1]-m[0])/sf))
        p=np.zeros((2,2));p[ds[0],0]=1-pf;p[ds[1],1]=pf
        return (v if k==0 else p.ravel()[k-1])*expit(t)*expit(-t)
    out=[sum(quad(node,lo,hi,args=(k,),epsabs=2e-12,epsrel=2e-12,limit=150)[0] for lo,hi in zip(points[:-1],points[1:])) for k in range(5)]
    return center+out[0],np.array(out[1:]).reshape(2,2)

class TwoShockTests(unittest.TestCase):
    def test_reference_all_feasibility_masks_and_scale_orders(self):
        rng=np.random.default_rng(231)
        for sh,sf in ((.005,2.17),(2.,.005),(.7,.4)):
            for mask in range(16):
                q=rng.normal(size=(2,2));q.ravel()[[not(mask&(1<<i)) for i in range(4)]]=-np.inf
                v,p=choose(q,sh,sf);rv,rp=reference(q,sh,sf)
                np.testing.assert_allclose(v,rv,atol=2e-9,rtol=0)
                np.testing.assert_allclose(p,rp,atol=2e-9,rtol=0)
    def test_independent_limit(self):
        q=np.array([[0.,.3],[-.2,.1]])
        v,p=choose(q,.4,.8)
        ph,pf=expit(-.2/.4),expit(.3/.8)
        np.testing.assert_allclose(p,np.outer([1-ph,ph],[1-pf,pf]),atol=2e-11,rtol=0)
        self.assertAlmostEqual(float(v),.4*np.logaddexp(0,-.2/.4)+.8*np.logaddexp(0,.3/.8),places=10)
    def test_value_gradient_symmetry_translation_refinement(self):
        rng=np.random.default_rng(436)
        for sh,sf in ((.005,2.),(2.,.005),(.7,.4)):
            for _ in range(5):
                q=rng.normal(size=(2,2));v,p=choose(q,sh,sf)
                vt,pt=choose(q.T,sf,sh)
                np.testing.assert_allclose(v,vt,atol=2e-10,rtol=0)
                np.testing.assert_allclose(p,pt.T,atol=2e-10,rtol=0)
                vv,pp=choose(q-1e5,sh,sf)
                np.testing.assert_allclose(vv+1e5,v,atol=5e-10,rtol=0)
                np.testing.assert_allclose(pp,p,atol=5e-10,rtol=0)
                vr,pr=choose(q,sh,sf,tolerance=2e-13)
                np.testing.assert_allclose(v,vr,atol=1e-10,rtol=0)
                np.testing.assert_allclose(p,pr,atol=1e-10,rtol=0)
                for d,a in itertools.product(range(2),repeat=2):
                    up=q.copy();down=q.copy();up[d,a]+=1e-5;down[d,a]-=1e-5
                    derivative=(choose(up,sh,sf)[0]-choose(down,sh,sf)[0])/2e-5
                    self.assertAlmostEqual(float(derivative),p[d,a],delta=2e-7)
    def test_bellman_uses_parity_scales_and_exact_joint_menu(self):
        from types import SimpleNamespace
        from intergen_eqscale_seq_optimized import joint_nested as joint
        P=SimpleNamespace(E_loc=np.array([0.]),mu_stay=0.,tenure_choice_kappa=.13,
            two_shock_choice=True,kappa_fert=.8,
            kappa_fert_continuation=None,A_f_start=1,A_f_end=2,
            n_parity=3,first_birth_fixed_cost=.2)
        q=np.array([[[[[.1,.2],[.3,.4],[.1,.7]],[[.3,.9],[.8,.5],[.3,.6]],[[.6,.4],[.9,.7],[.5,.9]]]]])
        # shape: wealth=1, tenure state=1, parity=3, children=3, chosen tenure=2
        replies=iter(((q[...,0],np.zeros(q.shape[:-1],dtype=int)),(q[...,1],np.ones(q.shape[:-1],dtype=int))))
        def kernel(*args):return next(replies)
        value,prob,_,wait=joint.bellman_block(np.zeros((1,2)),(np.array([False]),),P,0,np.array([.6]),kernel)
        plans=np.full(q.shape+(2,),-np.inf);plans[...,0]=q
        for n in range(2):
            for c in range(n+1):plans[...,n,c,:,1]=.4*q[...,n,c,:]+.6*(q[...,n+1,c+1,:]-(.2 if n==0 else 0.))
        expected,pr=choose(plans,.13,.8)
        np.testing.assert_allclose(value,expected,atol=1e-12)
        np.testing.assert_allclose(prob,pr,atol=1e-12)
        np.testing.assert_allclose(wait,joint.logsum_prob(q,.13)[1],atol=1e-12)

    def test_thin_layer_and_missing_diagonal(self):
        for ratio in (1e-4,1e-6):
            q=np.array([[0.,-np.inf],[-np.inf,.123]])
            v,p=choose(q,1.,ratio);rv,rp=reference(q,1.,ratio)
            np.testing.assert_allclose(v,rv,atol=2e-10,rtol=0)
            np.testing.assert_allclose(p,rp,atol=2e-10,rtol=0)

if __name__=='__main__':unittest.main()
