"""Independent smooth optimizers verify the exhaustive analytic saving oracle."""
import unittest
import numpy as np
from scipy.optimize import minimize_scalar

from run_e5f_independent_numerical_audit import segment_oracle
from intergen_eqscale_seq_optimized.kernels import eval_owner_scalar, eval_renter_scalar
import run_e5f_global_saving_quantification as quant


class SegmentOracleTests(unittest.TestCase):
    def test_global_blocks_dominate_local_and_balance_interior_synthetic_budgets(self):
        bg=np.array([-3.,-2.,-1.,0.,1.,2.,4.])
        R=1.08*bg+4.
        continuation=np.cumsum(np.random.default_rng(19).normal(.1,.3,(len(bg),2)),axis=0)
        previous=np.zeros_like(continuation)
        cb=np.array([.2,.3]);hb=np.array([.5,.8]);psi=np.zeros(2);gb=np.array([0.,.1])
        al=np.array([.7,.72]);es=np.array([1.,1.15])
        base=(R,R,continuation,previous,0,bg,cb,hb,psi,gb,al,es)
        tail=(.7,-1.,.9,.5,.4,(3-np.sqrt(5))/2,(np.sqrt(5)-1)/2,1e-3)
        rental_args=base+(.3,5.,.04,.2,.5)+tail
        old=quant.LOCAL_RENTER(*rental_args);new=quant.global_renter(*rental_args)
        self.assertTrue(np.all(new[0]>=old[0]-1e-10))
        np.testing.assert_allclose(new[2]+.3*new[3]+new[1],R[:,None]+np.clip(gb[None,:]-R[:,None],0,gb[None,:]),atol=1e-12,rtol=0)
        owner_args=base+(np.array([-2.,-2.]),.4,2.,1.,1.05,.04)+tail+(1,)
        old=quant.LOCAL_OWNER(*owner_args);new=quant.global_owner(*owner_args)
        self.assertTrue(np.all(new[0]>=old[0]-1e-10))
        np.testing.assert_allclose(new[2]+.4+new[1],R[:,None]+np.clip(gb[None,:]-R[:,None],0,gb[None,:]),atol=1e-12,rtol=0)

    def test_multipeaked_linear_continuation_and_rental_cap(self):
        rng=np.random.default_rng(1905)
        for owner in (False,True):
            for case in range(6):
                bg=np.array([-3.,-1.,0.,.25,1.,2.,4.,8.])
                continuation=np.cumsum(rng.normal(.2,.7,len(bg)))
                lo,hi=-1.5,7.5
                R,ri,hb,cb,pc,hmax,al,oms,beta,es,oc,Ko=9.,.6,.4,.1,.2,3.,.73,-1.,.96,1.18,.7,.82
                bp,value=segment_oracle(lo,hi,R,continuation,bg,ri,hb,cb,pc,hmax,al,oms,beta,es,oc,Ko,owner)
                if owner:
                    f=lambda x: float(eval_owner_scalar(x,R,continuation,bg,oc,cb,pc,Ko,al,oms,beta,es))
                else:
                    Kr=(al**al*((1-al)/ri)**(1-al))**oms
                    f=lambda x: float(eval_renter_scalar(x,R,continuation,bg,cb+ri*hb,pc,ri*(hmax-hb)/(1-al),cb,ri,hmax,hmax-hb,Kr,al,oms,beta,es))
                points=[lo,hi]+[x for x in bg if lo<x<hi]
                if not owner:
                    kink=R-cb-ri*hb-ri*(hmax-hb)/(1-al)
                    if lo<kink<hi: points.append(kink)
                points=sorted(points)
                reference=max(f(x) for x in points)
                for a,b in zip(points[:-1],points[1:]):
                    fit=minimize_scalar(lambda x:-f(x),bounds=(a,b),method="bounded",options={"xatol":1e-12})
                    reference=max(reference,-fit.fun)
                with self.subTest(owner=owner,case=case):
                    self.assertLessEqual(abs(value-reference),2e-9)
                    self.assertAlmostEqual(value,f(bp),places=12)


if __name__=="__main__": unittest.main()
