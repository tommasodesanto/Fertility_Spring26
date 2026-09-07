"""Global saving checks against independent dense objectives and audited oracle."""
import ast
import sys
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'code/model'))
from intergen_eqscale_seq_optimized import kernels as k


def plain(fn):
    return getattr(fn, 'py_func', fn)


def oracle_reference():
    source = ast.parse((ROOT/'code/model/tools/run_e5f_independent_numerical_audit.py').read_text())
    fn = next(n for n in source.body if isinstance(n, ast.FunctionDef) and n.name == 'segment_oracle')
    fn.decorator_list = []
    namespace = dict(np=np, eval_owner_scalar=plain(k.eval_owner_scalar), eval_renter_scalar=plain(k.eval_renter_scalar))
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<audited_oracle>', 'exec'), namespace)
    return namespace['segment_oracle']


class ExhaustiveSavingTests(unittest.TestCase):
    def test_strict_interpolation_rejects_positive_dead_endpoint_weight(self):
        bg=np.array([0.,1.,2.]);v=np.array([-1e10,-10.,-5.])
        self.assertGreater(k._interp_with_clip(bg,v,.99),-1e9)
        self.assertEqual(k._interp_with_clip(bg,v,.99,True),-1e10)
        for x in (0.,1.,2.,3.,-1.):
            self.assertEqual(k._interp_with_clip(bg,v,x,True),k._interp_with_clip(bg,v,x))
        live=np.array([-20.,-10.,-5.])
        for x in np.linspace(-1.,3.,51):
            self.assertEqual(k._interp_with_clip(bg,live,x,True),k._interp_with_clip(bg,live,x))

    def test_strict_tenure_kernel_rejects_unsupported_sale_buy_and_switch(self):
        bg=np.array([0.,1.,2.]);v=np.full((3,3,1,1,1),-1e10)
        heq=np.array([[0.,.99,2.]]);hc=np.array([[0.,1.01,2.]])
        dp=(.2*hc)[:,:,None,None];bm=(-.8*hc)[:,:,None,None]
        birth=np.zeros((1,1,3,3),dtype=np.bool_);grants=np.zeros((1,3,1,1))
        args=(bg,heq,hc,dp,bm,birth,grants)
        for destination,origin_b,origin_tenure in ((0,0,1),(1,2,0),(1,0,2)):
            values=v.copy();values[:,destination,0,0,0]=[-1e10,-10.,-5.]
            legacy,_=k.tenure_choice_kernel(values,*args)
            strict,_=k.tenure_choice_kernel(values,*args,True)
            idx=(origin_b,origin_tenure,0,0,0)
            self.assertGreater(legacy[idx],-1e9)
            self.assertEqual(strict[idx],-1e10)

    def setUp(self):
        self.bg = np.array([-3., -2., 0., .3, 1., 3., 6.])
        self.v = np.array([0., 1., 0., 4., 2., 7., 5.])

    def arguments(self, owner, oms=-.7):
        return (-2.8, 7.79, 8., self.v, self.bg, .7, .3, .2, .04,
                2., .64, oms, .96, 1.13, .8, 1.2, owner)

    def independent_values(self, args, grid):
        lo, hi, r, v, bg, rent, hb, cb, pc, hmax, al, oms, beta, es, cost, Ko, owner = args
        if owner:
            c = r-cost-cb-grid
            feasible = c > 1e-10
            utility = es*Ko*np.maximum(c, 1e-10)**(al*oms)/oms+pc
        else:
            surplus = r-cb-rent*hb-grid
            feasible = surplus > 1e-10
            h = np.minimum((1-al)*surplus/rent, hmax-hb)
            c = r-cb-rent*(hb+h)-grid
            utility = es*(np.maximum(c,1e-10)**al*np.maximum(h,1e-10)**(1-al))**oms/oms+pc
        return np.where(feasible, utility+beta*np.interp(grid,bg,v), -1e10)

    def test_global_maximum_dominates_dense_grid(self):
        rng = np.random.default_rng(721)
        for owner in (False, True):
            for oms in (-1.2, -.3, .4):
                for _ in range(5):
                    args = list(self.arguments(owner, oms)); args[3] = rng.normal(2.,3.,self.bg.size)
                    bp, val = plain(k.exhaustive_saving_scalar)(*args)
                    grid = np.unique(np.r_[np.linspace(args[0],args[1],6001), self.bg[(self.bg>args[0])&(self.bg<args[1])]])
                    self.assertGreaterEqual(val+1e-10, self.independent_values(args,grid).max())
                    self.assertAlmostEqual(val,float(self.independent_values(args,np.array([bp]))[0]),places=10)
                    self.assertTrue(args[0] <= bp <= args[1])

    def test_exact_audited_candidate_order_and_values(self):
        reference = oracle_reference()
        for owner in (False, True):
            for oms in (-.7,.3):
                args = self.arguments(owner,oms)
                self.assertEqual(plain(k.exhaustive_saving_scalar)(*args),reference(*args))

    def test_clipped_continuation_and_collapsed_interval(self):
        for owner in (False,True):
            args = list(self.arguments(owner));args[0]=6.2;args[1]=6.6
            bp, _ = plain(k.exhaustive_saving_scalar)(*args)
            self.assertEqual(bp,6.2)
            args[1]=args[0]
            self.assertEqual(plain(k.exhaustive_saving_scalar)(*args)[0],6.2)

    def test_borrowing_limit_is_an_explicit_candidate(self):
        for owner in (False,True):
            args = list(self.arguments(owner));args[3]=np.zeros_like(self.v)
            self.assertEqual(plain(k.exhaustive_saving_scalar)(*args)[0],args[0])

    def test_renter_cap_boundary_candidate(self):
        args=list(self.arguments(False))
        rent,hb,cb=args[5:8];hmax,al,oms,beta,es=args[9:14]
        surplus=rent*(hmax-hb)/(1-al);point=args[2]-cb-rent*hb-surplus
        Kr=(al**al*((1-al)/rent)**(1-al))**oms
        slope=es*Kr*surplus**(oms-1)/beta
        args[3]=slope*self.bg
        bp,val=plain(k.exhaustive_saving_scalar)(*args)
        self.assertAlmostEqual(bp,point,places=11)

    def test_rejects_uncovered_objective(self):
        args=list(self.arguments(False));args[11]=0.
        with self.assertRaises(ValueError):plain(k.exhaustive_saving_scalar)(*args)

    def test_renter_reports_budget_and_objective_below_legacy_floors(self):
        # Independent intratemporal checks, including the housing cap,
        # nonzero subsistence, a means-tested transfer and state-specific alpha.
        for resources, rent, hmax, cb, hb, al, grant in (
                (.03, 1., 2., 0., 0., .7, 0.),
                (.10, 10., 2., .01, .002, .8, 0.),
                (.20, 1., .025, .01, .02, .6, .1)):
            bg=np.array([0., 1.]); ones=np.ones(1); zeros=np.zeros((2,1))
            args=(np.full(2,resources),np.zeros(2),zeros,zeros,0,bg,
                  ones*cb,ones*hb,ones*.04,ones*grant,ones*al,ones*1.13,
                  rent,hmax,.04,0.,0.,.64,-.7,.96,0.,0.,.381966,.618034,1e-3,1)
            value,saving,c,h=k.full_renter_block_kernel(*args)
            available=resources+grant
            surplus=available-cb-rent*hb-saving
            expected_h=hb+np.minimum((1-al)*surplus/rent,hmax-hb)
            expected_c=available-rent*expected_h-saving
            np.testing.assert_allclose(h,expected_h,rtol=0,atol=1e-16)
            np.testing.assert_allclose(c,expected_c,rtol=0,atol=1e-16)
            utility=1.13*((c-cb)**al*(h-hb)**(1-al))**(-.7)/(-.7)+.04
            np.testing.assert_allclose(value,utility,rtol=2e-15,atol=1e-13)
        args=list(self.arguments(False));args[1]=args[0]-1
        with self.assertRaises(ValueError):plain(k.exhaustive_saving_scalar)(*args)


if __name__ == '__main__':unittest.main()
