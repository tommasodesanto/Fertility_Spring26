"""Small accounting fixtures; no equilibrium solves. Run with NUMBA_DISABLE_JIT=1."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intergen_eqscale_seq_optimized import solver as s
from intergen_eqscale_seq_optimized import kernels as k
from intergen_eqscale_seq_optimized.utils import make_grid


class NativeCurrentHousehold(unittest.TestCase):
    def setUp(self):
        self.bg = np.array([-2., -1., 0., 1., 2.])
        self.P = NS(I=1, z_grid=np.array([.5, 1.5]), Nb=5, b_max=2.,
                    fixed_reference_entry_grid=self.bg.copy(),
                    earnings_transaction_grid=self.bg.copy(),
                    fixed_reference_entry_conditional=np.array([[.2, 0], [.8, 0], [0, .3], [0, .7], [0, 0.]]))

    def test_explicit_grid_and_joint_entry(self):
        s.configure_current_household_contract(self.P)
        grid = make_grid(self.P)
        np.testing.assert_array_equal(grid, self.bg)
        self.assertFalse(np.shares_memory(grid, self.P.earnings_transaction_grid))
        for iz, z in enumerate(self.P.z_grid):
            ix, w = s.entry_wealth_grid_weights(grid, self.P, z_value=z)
            restored = np.zeros(5); restored[ix] = w
            np.testing.assert_array_equal(restored, self.P.fixed_reference_entry_conditional[:, iz])
        with self.assertRaises(ValueError):
            s.entry_wealth_grid_weights(grid, self.P, z_value=1.)
        with self.assertRaises(ValueError):
            s.entry_wealth_grid_weights(grid, self.P, j=1, z_value=.5)
        with self.assertRaises(AssertionError):
            s.entry_wealth_grid_weights(grid + .01, self.P, z_value=.5)

    def test_entry_rejects_invalid_mass(self):
        s.configure_current_household_contract(self.P)
        self.P.fixed_reference_entry_conditional[0, 0] = -.1
        with self.assertRaises(ValueError):
            s.entry_wealth_grid_weights(self.bg, self.P, z_value=.5)

    def test_forward_is_exact_transaction_without_collateral_transfer(self):
        P = NS(native_purchase_income=True)
        hc = np.array([[0., 1.]])
        he = .94 * hc
        phi = np.ones((1, 2, 1, 1)) * .8
        bd = np.zeros((1, 1, 2, 2), bool)
        grant = np.zeros_like(phi)
        ix, wt = s.build_forward_tenure_transition_maps(P, self.bg, hc, he, phi, bd, grant)
        for old in range(2):
            for new in range(2):
                x = self.bg if old == new else self.bg + he[0, old] - hc[0, new]
                mapped = (1-wt[0,old,new,0,0])*self.bg[ix[0,old,new,0,0]] + wt[0,old,new,0,0]*self.bg[ix[0,old,new,0,0]+1]
                np.testing.assert_allclose(mapped, np.clip(x, -2., 2.), atol=1e-15)
        # At b=0 purchase wealth is -1, not collateral-clipped to -.8.
        self.assertEqual(float((1-wt[0,0,1,0,0,2])*self.bg[ix[0,0,1,0,0,2]] + wt[0,0,1,0,0,2]*self.bg[ix[0,0,1,0,0,2]+1]), -1.)

    def test_support_rejects_sale_above_grid(self):
        values = np.zeros((5, 2, 1, 1, 1)); values[:, 0] = 10.
        hc = np.array([[0., 1.]]); he = .94*hc
        floor = np.full((1,2,1,1), -.8); dp = np.full_like(floor,.2)
        bd = np.zeros((1,1,2,2),bool); grant=np.zeros_like(floor)
        _, legacy = k.tenure_choice_kernel(values,self.bg,he,hc,dp,floor,bd,grant,values)
        _, exact = k.tenure_choice_kernel(values,self.bg,he,hc,dp,floor,bd,grant,values,False,True)
        self.assertEqual(legacy[-1,1,0,0,0],0)
        self.assertEqual(exact[-1,1,0,0,0],1)
        _, _, probs = k.tenure_logit_kernel(values,self.bg,he,hc,dp,floor,bd,grant,.1,values,True)
        self.assertEqual(probs[-1,1,0,0,0,0],0.)

    def test_purchase_income_changes_screen_without_adding_wealth(self):
        values=np.zeros((5,2,1,1,1)); values[:,1]=10.
        hc=np.array([[0.,1.]]); he=.94*hc
        floor=np.full((1,2,1,1),-.8); dp=np.full_like(floor,.2)
        bd=np.zeros((1,1,2,2),bool); grant=np.zeros_like(floor)
        _, legacy=k.tenure_choice_kernel(values,self.bg,he,hc,dp,floor,bd,grant,values)
        income_over_R=.3
        _, exact=k.tenure_choice_kernel(values,self.bg,he,hc,dp-income_over_R,np.maximum(floor-income_over_R,self.bg[0]),bd,grant,values,False,True)
        self.assertEqual(legacy[2,0,0,0,0],0)
        self.assertEqual(exact[2,0,0,0,0],1)

    def test_exact_owner_allocation_preserves_value_saving(self):
        bg=np.array([0., .01, .1]); v=np.zeros((3,1)); z=np.zeros(1); one=np.ones(1)
        args=(np.full(3,.001),np.full(3,.001),v,v,0,bg,z,z,z,z,.7*one,one,z,0.,1.,1.,1.,.1,.7,-1.,.9,0.,0.,.381966,.618034,1e-8,0,1)
        legacy=k.full_owner_block_kernel(*args)
        exact=k.full_owner_block_kernel(*args,exact_allocation_output=True)
        np.testing.assert_array_equal(legacy[0],exact[0]);np.testing.assert_array_equal(legacy[1],exact[1])
        np.testing.assert_allclose(exact[2]+exact[1],.001,atol=1e-15)
        self.assertTrue(np.all(legacy[2]>=.1))

    def test_exact_renter_allocation_with_deep_continuation(self):
        bg=np.array([0., .01, .1]); v=np.full((3,1),-2e9); z=np.zeros(1); one=np.ones(1)
        args=(np.full(3,.001),np.full(3,.001),v,np.zeros_like(v),0,bg,
              z,z,z,z,.7*one,one,1.,2.,.1,0.,0.,.7,-1.,.9,0.,0.,.381966,.618034,1e-8,1)
        legacy=k.full_renter_block_kernel(*args)
        exact=k.full_renter_block_kernel(*args,exact_allocation_output=True)
        np.testing.assert_array_equal(legacy[0],exact[0]); np.testing.assert_array_equal(legacy[1],exact[1])
        np.testing.assert_allclose(exact[1]+exact[2]+exact[3],.001,atol=1e-15)
        self.assertTrue(np.any(legacy[2]+legacy[3] > .001))

    def test_owner_zero_rollover_enforces_collateral_floor(self):
        bg=np.array([-2.,-1.,0.,1.]); v=np.zeros((4,1)); z=np.zeros(1); one=np.ones(1)
        args=(np.full(4,1.),np.full(4,1.),v,v,0,bg,z,z,z,z,.7*one,one,-.8*one,
              0.,1.,1.,1.,.1,.7,-1.,.9)
        remainder=(.381966,.618034,1e-8,0,1)
        legacy=k.full_owner_block_kernel(*args,1.,0.,*remainder)
        exact=k.full_owner_block_kernel(*args,0.,0.,*remainder)
        self.assertLess(legacy[1][0,0],-.8)
        self.assertTrue(np.all(exact[1]>=-.8))


if __name__ == '__main__':
    unittest.main()
