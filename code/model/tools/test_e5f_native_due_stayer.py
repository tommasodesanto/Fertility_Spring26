"""Small DUE stayer fixtures; no household life-cycle or equilibrium solve."""
import sys
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intergen_eqscale_seq_optimized import solver as s, kernels as k

class DueStayer(unittest.TestCase):
    def test_floor_and_interest_units(self):
        b=np.array([-80.,-72.,-60.,10.])
        np.testing.assert_array_equal(s.native_due_owner_floor(b,-72.),[-80.,-72.,-72.,-72.])
        # Principal -80 cannot become -84 just because its interest is due.
        R=1.05; income=10.; saving=-80.
        self.assertEqual(R*b[0]+income-saving,6.)
        self.assertGreater(saving, R*b[0])

    def test_death_solvency_is_separate(self):
        P=NS(J=3,use_age_survival=True,survival_probs=[1.,.9],psi=.06)
        self.assertEqual(s.native_due_death_floor(P,0,70.,1.),-np.inf)
        for j in (1,2):
            floor=s.native_due_death_floor(P,j,70.,1.)
            self.assertAlmostEqual(float(s.native_due_owner_floor(-80.,-56.,death_floor=floor)),-65.8)

    def test_native_default_and_explicit_flag(self):
        P=NS(J=3,native_purchase_income=True,owner_ltv_multipliers=np.ones(4))
        np.testing.assert_array_equal(s.owner_borrowing_floor(P,[-2.,0.],-.8,0),[-.8,-.8])
        with self.assertRaises(ValueError): s.owner_borrowing_floor(P,-2.,-.8,0,stay_on=True)
        P.native_due_stayer_credit=True
        np.testing.assert_array_equal(s.owner_borrowing_floor(P,[-2.,0.],-.8,0,stay_on=True),[-2.,-.8])
        np.testing.assert_array_equal(s.owner_borrowing_floor(P,[-2.,0.],-.8,0),[-.8,-.8])

    def test_owner_kernel_distinct_buyer_and_stayer(self):
        bg=np.array([-2.,-1.,0.,1.]); v=np.zeros((4,1)); z=np.zeros(1); one=np.ones(1)
        args=(np.full(4,1.),np.full(4,1.),v,v,0,bg,z,z,z,z,.7*one,one,-.8*one,
              0.,1.,1.,1.,.1,.7,-1.,.9,0.,0.,.381966,.618034,1e-8,0,1)
        default=k.full_owner_block_kernel(*args,exact_allocation_output=True)
        off=k.full_owner_block_kernel(*args,exact_allocation_output=True,due_stayer=False)
        for x,y in zip(default,off): np.testing.assert_array_equal(x,y)
        due=k.full_owner_block_kernel(*args,exact_allocation_output=True,due_stayer=True)
        death=k.full_owner_block_kernel(*args,exact_allocation_output=True,due_stayer=True,due_death_floor=-.9)
        self.assertGreaterEqual(default[1][0,0],-.8)
        self.assertLess(due[1][0,0],-.8)
        self.assertGreaterEqual(death[1][0,0],-.9)
        for out in (default,due,death): np.testing.assert_allclose(out[1]+out[2],1.,atol=1e-14)
        np.testing.assert_array_equal(default[1][2:],due[1][2:])

    def test_actual_savings_fallback_rejects_due(self):
        # Exercise the stage itself before it can reach the incomplete fallback.
        with self.assertRaisesRegex(NotImplementedError, "explicit death solvency"):
            s._savings_stage(None, NS(native_due_stayer_credit=True), np.array([-1.,0.]),
                             None, NS(use_full_kernel=False), None, 0, 1., 0., 0., None,
                             stay_floor=True)

    def test_packed_solution_retains_stayer_payload(self):
        shape=(1,2,1,1,1,1,1)
        a=np.zeros(shape); stay=np.ones(shape)*.2; bs=a-.8; cs=a+.3
        P=NS(I=1,n_house=1,H_own=[1.],user_cost_rate=.1,H0=1.,r_bar=.1,xi_supply=0.,
             _bp_pol_stay=bs,_c_pol_stay=cs)
        st=NS(g_stay_distribution=stay)
        with patch.object(s,"income_transition_values",return_value=(np.ones(1),np.ones(1),np.ones((1,1)))), patch.object(s,"housing_demand_normalizer",return_value=1.):
            sol=s.pack_solution_markov_income(a,a,a,a,a,None,a,a,a,a,st,np.ones(1),np.ones(1),P)
        self.assertIs(sol.g_stay_distribution,stay)
        self.assertIs(sol.bp_pol_stay,bs)
        self.assertIs(sol.c_pol_stay,cs)

    def test_origin_specific_estate_moment(self):
        shape=(1,2,1,1,1,1,1)
        mass=np.zeros(shape); mass[:,1]=1.
        bp=np.zeros(shape); bp[:,1]=-.5
        stay=np.zeros(shape); stay[:,1]=.25
        bs=bp.copy(); bs[:,1]=-.8
        P=NS(J=1,I=1,J_R=1,z_grid=[1.],period_years=4.,tau_pay=0.,
             income=np.ones((1,1))*4.,H_own=[1.],native_due_stayer_credit=True,
             _g_stay_distribution=stay,_bp_pol_stay=bs)
        st=NS()
        s.add_aggregate_wealth_bequest_flow_moments(st,mass,mass,bp,P,np.array([0.]),np.array([1.]))
        self.assertAlmostEqual(st.annual_bequest_flow,(.75*.5+.25*.2)/4.)

    def test_origin_mass_keeps_stayers_only(self):
        g=np.ones((2,3,1,1,1,1,1))
        lp=np.ones((2,3,1,1,1,1,1,1))
        tc=np.ones_like(g,dtype=np.int16)
        tp=np.zeros(g.shape+(3,)); tp[...,0]=.25; tp[...,1]=.5; tp[...,2]=.25
        out=s.realize_stayer_cross_section(g,lp,tc,tp)
        self.assertTrue(np.all(out[:,0]==0.))
        self.assertTrue(np.all(out[:,1]==.5))
        self.assertTrue(np.all(out[:,2]==.25))
        # A location departure is not grandfathered ownership.
        lp[:,1]=.2
        self.assertTrue(np.all(s.realize_stayer_cross_section(g,lp,tc,tp)[:,1]==.1))
        deterministic=s.realize_stayer_cross_section(g,lp,tc,None)
        self.assertTrue(np.all(deterministic[:,2]==0.))

if __name__=='__main__': unittest.main()
