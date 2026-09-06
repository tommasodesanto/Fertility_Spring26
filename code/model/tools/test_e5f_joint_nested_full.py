"""Independent reference checks for the experimental joint operator."""
import sys
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
from intergen_eqscale_seq_optimized import joint_nested as joint
from e5f_joint_nested_choice import choose, scatter_joint_block

class JointFullTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(2345)
        self.P=SimpleNamespace(I=1,J=2,n_parity=4,n_child_states=4,A_f_start=1,A_f_end=1,
            fecundity_omega1=.4,fecundity_omega2=0.,age_start=20.,da=28.,fecundity_terminal_age=45.)
    def fixture(self):
        shape=(5,3,1,2,2,4,4)
        plans=self.rng.normal(size=shape+(2,2))
        plans[:,:,:,1,...,1]=-np.inf
        plans[...,3,:,:,1]=-np.inf
        _,prob=choose(plans,.7,.4)
        prod=np.zeros(shape+(2,),dtype=int);prod[...,1]=self.rng.integers(1,3,size=shape)
        _,wait=joint.logsum_prob(plans[...,0],.7)
        obj=SimpleNamespace(probabilities=prob,products=prod,wait_probabilities=wait)
        g=self.rng.random(shape)
        for n in range(4):g[...,n,n+1:]=0
        return g/g.sum(),obj
    def test_nested_vs_independent_reference(self):
        plans=self.rng.normal(size=(17,2,2));plans[1,1,:]=-np.inf
        for lam in (1,.5,.001):
            inner,pa=joint.logsum_prob(plans,.8*lam)
            value,pt=joint.logsum_prob(inner,.8)
            expected,prob=choose(plans,.8,lam)
            np.testing.assert_array_equal(value,expected)
            np.testing.assert_array_equal(pa*pt[...,None],prob)
    def test_factor_is_exact_direct_scatter(self):
        g,obj=self.fixture(); post,kernel,born,attempt,risk=joint.factor_distribution(g,obj,self.P)
        nb,nt=g.shape[:2]
        idx=np.broadcast_to(np.arange(nb-1).tolist()+[nb-2],(nt,nt,4,4,nb)).copy()
        wt=np.zeros_like(idx,dtype=float);wt[...,-1]=1
        for j in range(2):
            for z in range(2):
                current=np.zeros((nb,nt,4,4));before=np.zeros_like(current)
                for n in range(4):
                    for c in range(n+1):
                        dest=(n+1,c+1) if n<3 and j==0 else None
                        c1=obj.products[:,:,0,j,z,n+1,c+1] if dest else obj.products[:,:,0,j,z,n,c]
                        direct,pst,_=scatter_joint_block(g[:,:,0,j,z,n,c],obj.probabilities[:,:,0,j,z,n,c],
                            obj.products[:,:,0,j,z,n,c],c1,.6 if j==0 else 0,dest,(n,c),idx,wt)
                        current+=direct;before+=pst
                compressed=np.sum(post[:,:,0,j,z,...,None]*kernel[:,:,0,j,z],axis=1).transpose(0,3,1,2)
                np.testing.assert_allclose(current,compressed,rtol=0,atol=2e-17)
                np.testing.assert_allclose(before,post[:,:,0,j,z],rtol=0,atol=2e-17)
        self.assertAlmostEqual(post.sum(),1,places=14)
        self.assertAlmostEqual(born.sum(),.6*attempt.sum(),places=14)
    def test_selected_origin_tenure_preserved(self):
        g,obj=self.fixture()
        t,kt,*_=joint.factor_distribution(g,obj,self.P,'first_birth_treated')
        c,kc,*_=joint.factor_distribution(g,obj,self.P,'first_birth_control')
        self.assertAlmostEqual(t.sum(),c.sum(),places=14)
        # Within-tenure outcome products may change, but renter/owner split cannot.
        mt=t[...,None]*kt;mc=c[...,None]*kc
        np.testing.assert_allclose(mt[...,0].sum(),mc[...,0].sum(),rtol=0,atol=1e-16)
        self.assertEqual(c[...,1:,:].sum(),0)
    def test_action_marginals_have_exact_support(self):
        probability=np.array([[[.7,.1],[.20000000000000015,0.]]])
        # A tiny sum-roundoff in wait is removed without changing birth mass.
        marginal=joint.action_marginals(probability)
        self.assertTrue(np.all((marginal>=0)&(marginal<=1)))
        np.testing.assert_array_equal(marginal.sum(axis=-1),[1.])
        self.assertEqual(marginal[0,1],.1)
        with self.assertRaises(RuntimeError):joint.action_marginals(probability*4)

    def test_wait_has_no_births(self):
        g,obj=self.fixture();post,k,b,a,r=joint.factor_distribution(g,obj,self.P,'wait')
        np.testing.assert_allclose(post,g,rtol=0,atol=1e-17)
        self.assertEqual(b.sum(),0)

if __name__=='__main__':unittest.main()
