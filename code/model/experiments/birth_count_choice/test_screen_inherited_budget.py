import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location("screen", Path(__file__).with_name("screen_inherited_budget.py"))
screen = importlib.util.module_from_spec(spec); spec.loader.exec_module(screen)

def packet(*, b=0., old=0, last=False, hb=1., phi=.8, strict=False):
    P = NS(b_grid=np.array([-2., 0.]), n_house=1, I=1, J=2, Nz=1, H_own=np.array([2.]), income=np.array([[1.,1.]]), z_grid=np.array([1.]), J_R=2, R_gross=1.1, psi=.1, delta=.05, tau_H=.05, pension=1., retirement_income_z_scale=0., native_purchase_income=True, native_due_stayer_credit=True, unsecured_credit_limit=0., child_room_floor=strict, hbar_child_rooms=1., hbar_first_child_jump=0., owner_h_bar_scale=1., hR_max=2., owner_ltv_multipliers=np.ones(1), owner_size_cost=0., child_earnings_penalty=False, rental_wedge=False, use_age_survival=False)
    S = NS(phi_choice=np.full((1,2,2,1), phi), cb_flat=np.array([.2,.2]), hb_flat=np.array([hb,hb]), gb_flat=np.zeros(2))
    g=np.zeros((2,2,1,2,1,2,1)); g[1 if b == 0 else 0,old,0,1 if last else 0,0,1,0]=1.
    return {"parameters":P,"b_grid":P.b_grid,"shared":S,"stationary_g_pre":g}

class TestStayWitness(unittest.TestCase):
    def test_owner_same_tenure_has_no_transaction_charge(self):
        out=screen.screen(packet(b=0.,old=1),q=1.,r=1.,pension=1.,scale=1.)
        self.assertEqual(out["witnessed_stay_mass"],1.)
    def test_due_floor_is_minimum_of_current_debt_and_collateral(self):
        p=packet(b=-2.,old=1,phi=.8); p["parameters"].income[0,0]=.5
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.)
        self.assertAlmostEqual(out["top_unresolved_margins"][0]["margin"], 1.1*-2+.5-(-2.)-(.2+.2))
    def test_last_age_death_floor_binds_independently(self):
        p=packet(b=-2.,old=1,last=True); p["parameters"].income[0,1]=.1
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.)
        expected=1.1*(-2)+.1-(-1.8)-(.2+.2)
        self.assertAlmostEqual(out["top_unresolved_margins"][0]["margin"],expected)
    def test_zero_credit_renter_threshold_and_limit(self):
        p=packet(b=-2.,old=0); p["parameters"].income[0,0]=2.2 # exact margin is zero
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.); self.assertEqual(out["unresolved_mass"],1.)
        p=packet(b=0.,old=0,hb=2.); out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.); self.assertEqual(out["unresolved_mass"],1.)
    def test_f_order_children_and_no_mass_mutation(self):
        p=packet(b=0.,old=0); p["shared"].cb_flat=np.array([.1,9.]); p["shared"].hb_flat=np.array([.1,9.]); before=p["stationary_g_pre"].copy()
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.); self.assertEqual(out["unresolved_mass"],1.); self.assertTrue(np.array_equal(before,p["stationary_g_pre"]))
    def test_changed_pension_rejected(self):
        with self.assertRaises(ValueError): screen.screen(packet(),q=1.,r=1.,pension=1.1,scale=1.)
    def test_owner_exit_only_witness(self):
        p=packet(b=-2.,old=1,hb=.1); p["parameters"].income[0,0]=1.05; p["parameters"].delta=.25; p["shared"].cb_flat[:]=.5
        R, q, psi, H, b, income, cbar, r, hbar = 1.1, 1., .1, 2., -2., 1.05, .5, 1., .1
        stay_margin=R*b+income-max(p["parameters"].b_grid[0],min(b,-p["shared"].phi_choice[0,1,1,0]*q*H))-(cbar+(p["parameters"].delta+p["parameters"].tau_H)*q*H)
        x=b+(1-psi)*q*H/R
        exit_margin=R*x+income-max(p["parameters"].b_grid[0],0)-cbar-r*hbar
        self.assertAlmostEqual(stay_margin,-.25)
        self.assertAlmostEqual(x,-.36363636363636354)
        self.assertAlmostEqual(exit_margin,.05)
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.)
        self.assertEqual(out["witnessed_stay_mass"],0.)
        self.assertEqual(out["witnessed_exit_only_mass"],1.)
    def test_exit_grid_domain_and_strict_owner_room(self):
        p=packet(b=0.,old=1,hb=2.,strict=True)
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.)
        self.assertEqual(out["witnessed_stay_mass"],0.)
        p["parameters"].H_own[0]=3.
        out=screen.screen(p,q=1.,r=1.,pension=1.,scale=1.)
        self.assertEqual(out["witnessed_stay_mass"],1.)

if __name__ == "__main__": unittest.main()
