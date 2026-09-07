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
        args=list(self.arguments(False));args[1]=args[0]-1
        with self.assertRaises(ValueError):plain(k.exhaustive_saving_scalar)(*args)


if __name__ == '__main__':unittest.main()
