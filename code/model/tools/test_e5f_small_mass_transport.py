"""Mass conservation and scale invariance for rare matched birth cohorts."""
from pathlib import Path
from types import SimpleNamespace as NS
import sys
from unittest import TestCase,main
from unittest.mock import Mock
import numpy as np

MODEL=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(MODEL),str(MODEL/'tools')]
from intergen_eqscale_seq_optimized import solver,utils
import run_e5f_transition_calibration as calibration


def fixture(compiled=False):
    b=np.array([0.,1.]);shape=(2,2,1,2,1,2,2)
    P=NS(n_house=1,I=1,n_parity=2,n_child_states=2,n_child_stages=0,use_numba_scatter=compiled)
    shared=NS(nc=4)
    loc=np.ones((2,2,1,1,2,1,2,2));ten=np.zeros(shape+(2,))
    ten[...,0]=1-1e-7;ten[...,1]=1e-7
    bp=np.broadcast_to(b.reshape(2,1,1,1,1,1,1),shape).copy()
    lidx=np.zeros((1,2,2),dtype=np.int64);lwt=np.broadcast_to(b,lidx.shape).copy()
    tidx=np.zeros((1,2,2,2,2,2),dtype=np.int64);twt=np.broadcast_to(b,tidx.shape).copy()
    g=np.zeros((2,2,1,1,2,2));g[0,0,0,0,1,1]=1
    choices=np.zeros(shape,dtype=np.int64)
    args=(0,loc,choices,ten,bp,P,b,shared,lidx,lwt,tidx,twt,False,None,np.ones((1,1)))
    current=(loc,choices,ten,lidx,lwt,tidx,twt)
    expected=g.copy();expected[0,0,0,0,1,1]=1-1e-7;expected[0,1,0,0,1,1]=1e-7
    return g,args,current,expected


class SmallMassTransportTests(TestCase):
    def test_default_retains_original_behavior(self):
        g,args,_,_=fixture()
        ordinary=solver.advance_cohort_one_period_markov_income(3.55e-9*g,*args)
        explicit=solver.advance_cohort_one_period_markov_income(3.55e-9*g,*args,mass_pruning_tolerance=1e-15)
        np.testing.assert_array_equal(ordinary,explicit)
        self.assertEqual(float(ordinary[:,1].sum()),0.)
        self.assertGreater(abs(float(ordinary.sum())/3.55e-9-1),5e-9)

    def test_precise_transition_matches_closed_form_at_all_scales(self):
        for compiled in ([False,True] if solver.NUMBA_AVAILABLE else [False]):
            g,args,_,expected=fixture(compiled)
            for mass in [1.,3.55e-9,1e-13,1e-20]:
                result=solver.advance_cohort_one_period_markov_income(mass*g,*args,mass_pruning_tolerance=0.)
                np.testing.assert_allclose(result/mass,expected,rtol=0,atol=3e-16)
                self.assertGreater(float(result[:,1].sum()),0.)
                self.assertLess(abs(float(result.sum())/mass-1),3e-16)

    def test_precise_current_choice_matches_closed_form_at_all_scales(self):
        for compiled in ([False,True] if solver.NUMBA_AVAILABLE else [False]):
            g,_,args,expected=fixture(compiled)
            for mass in [1.,3.55e-9,1e-13,1e-20]:
                cross=np.zeros((2,2,1,2,1,2,2));cross[:,:,:,0,:,:,:]=mass*g
                result=solver.realize_current_cross_section(cross,*args,use_compiled_scatter=compiled,mass_pruning_tolerance=0.)
                np.testing.assert_allclose(result[:,:,:,0,:,:,:]/mass,expected,rtol=0,atol=3e-16)
                self.assertEqual(float(result[:,:,:,1,:,:,:].sum()),0.)

    def test_numpy_scatter_keeps_tiny_positive_columns(self):
        mass=np.array([[.2,1e-21],[.8,3e-21]])
        idx=np.zeros_like(mass,dtype=np.int64);weights=np.array([[0.,0.],[1.,1.]])
        np.testing.assert_array_equal(utils.scatter_redistribute_cols(idx,weights,mass,2,mass_pruning_tolerance=0.),mass)
        np.testing.assert_array_equal(utils.scatter_redistribute_cols_sameidx(idx[:,0],weights[:,0],mass,2,mass_pruning_tolerance=0.),mass)

    def test_passing_transport_never_retries_and_is_exact(self):
        original=np.array([.2,.8])*(1-2e-10);retry=Mock(side_effect=AssertionError('unneeded retry'))
        expected,old_gate=calibration.normalize_distribution_mass_roundoff(original,expected_mass=1.,stage='test')
        actual,gate=calibration.normalize_branch_transport_mass(original,expected_mass=1.,stage='test',retry_without_pruning=retry)
        np.testing.assert_array_equal(actual,expected);self.assertEqual(gate,old_gate);retry.assert_not_called()

    def test_failed_transport_retry_preserves_original_relative_gate(self):
        g,args,_,_=fixture();mass=3.55e-9
        first=solver.advance_cohort_one_period_markov_income(mass*g,*args)
        retry=lambda:solver.advance_cohort_one_period_markov_income(mass*g,*args,mass_pruning_tolerance=0.)
        result,gate=calibration.normalize_branch_transport_mass(first,expected_mass=mass,stage='test',retry_without_pruning=retry)
        self.assertLess(abs(float(result.sum())/mass-1),3e-16)
        self.assertEqual(gate['relative_tolerance'],5e-9)
        self.assertGreater(gate['initial_relative_gap'],gate['relative_tolerance'])
        with self.assertRaisesRegex(RuntimeError,'mass gate failed'):
            calibration.normalize_branch_transport_mass(first,expected_mass=mass,stage='test',retry_without_pruning=lambda:.9*mass*g)

    def test_invalid_original_mass_never_retries(self):
        for values,expected in [(np.array([np.nan]),1.),(np.array([-1.]),1.),(np.array([1.]),np.nan),(np.array([1.]),-1.)]:
            retry=Mock(side_effect=AssertionError('invalid-input retry'))
            with self.assertRaises(RuntimeError):
                calibration.normalize_branch_transport_mass(values,expected_mass=expected,stage='test',retry_without_pruning=retry)
            retry.assert_not_called()

    def test_invalid_retried_distribution_is_rejected(self):
        for retried in [np.array([-.1,1.1]),np.array([np.nan,0.]),np.ones((1,1))]:
            with self.assertRaisesRegex(RuntimeError,'invalid distribution'):
                calibration.normalize_branch_transport_mass(np.array([.9,0.]),expected_mass=1.,stage='test',retry_without_pruning=lambda:retried)

    def test_default_call_still_rejects_mass_failure(self):
        with self.assertRaisesRegex(RuntimeError,'mass gate failed'):
            calibration.normalize_branch_transport_mass(np.array([.9]),expected_mass=1.,stage='test')


if __name__=='__main__':main()
