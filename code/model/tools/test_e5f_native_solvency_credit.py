"""Synthetic support/operator tests only; no household or equilibrium solves."""
import inspect
import sys
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intergen_eqscale_seq_optimized import solver, kernels
from unittest.mock import patch


class NativeSolvencyTests(unittest.TestCase):
    def test_default_off_is_noop(self):
        self.assertFalse(solver.validate_native_solvency_mode(SimpleNamespace()))
        self.assertFalse(solver.validate_native_solvency_mode(SimpleNamespace(native_solvency_credit=False)))
        with self.assertRaises(ValueError):
            solver.validate_native_solvency_mode(SimpleNamespace(native_solvency_credit=True))

    def test_supported_mode_requires_all_native_flags(self):
        p = SimpleNamespace(native_solvency_credit=True,native_purchase_income=True,
            native_fixed_reference_entry=True,native_explicit_transaction_grid=True,
            native_exact_allocation_output=True,exhaustive_saving_control=True)
        with patch.object(solver,'parent_age_maturation_active',return_value=False), \
             patch.object(solver,'mortgage_stay_floor_active',return_value=False), \
             patch.object(solver,'bequest_utility_net_active',return_value=False), \
             patch.object(solver,'estate_receiver_active',return_value=False):
            self.assertTrue(solver.validate_native_solvency_mode(p))
            p.native_exact_allocation_output=False
            with self.assertRaises(ValueError):
                solver.validate_native_solvency_mode(p)

    def test_feasible_boundary_and_cutoff_margin(self):
        cutoff=solver.DEAD_VALUE_CUTOFF
        grid=np.array([-2.,-1.,0.,1.])
        v=np.array([[cutoff, -1e10],[cutoff+1, -1e10],[0.,-1e10],[1.,-1e10]])
        floors, dead=solver.native_solvency_support_floor(v,grid)
        np.testing.assert_array_equal(floors,[-1.,1.])
        np.testing.assert_array_equal(dead,[False,True])
        v[2,0]=-1e10
        with self.assertRaises(ValueError):
            solver.native_solvency_support_floor(v,grid)

    def test_current_price_net_death_estate(self):
        # Gross house price 10, sale cost .06: debt -9.4 is exactly solvent.
        bg=np.array([-9.4000001,-(1.-.06)*10.,-9.,0.])
        np.testing.assert_array_equal(solver.native_solvency_death_mask(bg,10.,.06),
                                     [True,False,False,False])
        np.testing.assert_array_equal(solver.native_solvency_death_mask(bg,0.,.06),
                                     [True,True,True,False])

    def test_dated_continuation_changes_support(self):
        grid=np.array([-1.,0.,1.]); death=np.zeros((3,1))
        stationary=np.ones((2,3,1))
        dated=stationary.copy();dated[1,0,0]=-1e10
        old=solver.native_solvency_continuation(stationary,[.999999, .000001],death,1.)
        new=solver.native_solvency_continuation(dated,[.999999, .000001],death,1.)
        self.assertEqual(solver.native_solvency_support_floor(old,grid)[0][0],-1.)
        self.assertEqual(solver.native_solvency_support_floor(new,grid)[0][0],0.)
        no_reach=solver.native_solvency_continuation(dated,[1.,0.],death,1.)
        np.testing.assert_array_equal(no_reach,old)

    def test_possible_death_and_certain_death(self):
        death=np.array([[-1e10],[2.]])
        nextv=np.array([[[1.],[-1e10]]])
        mixed=solver.native_solvency_continuation(nextv,[1.],death,.99999)
        np.testing.assert_array_equal(mixed,np.full((2,1),-1e10))
        certain=solver.native_solvency_continuation(nextv,[1.],death,0.)
        np.testing.assert_array_equal(certain,death)

    def test_renter_argument_added_at_end(self):
        fn=getattr(kernels.full_renter_block_kernel,'py_func',kernels.full_renter_block_kernel)
        names=list(inspect.signature(fn).parameters)
        self.assertEqual(names[25:33],['exhaustive_saving','yadj_v','pen_on',
            'wedge_on','w0','w1','hk','exact_allocation_output'])
        self.assertEqual(names[33],'natural_floor_v')
        self.assertIsNone(inspect.signature(fn).parameters['natural_floor_v'].default)

if __name__=='__main__':
    unittest.main()
