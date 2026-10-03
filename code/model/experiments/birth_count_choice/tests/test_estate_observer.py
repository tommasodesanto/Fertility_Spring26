"""No-solve frozen-observer estate checks, including retained real arrays."""
from pathlib import Path
import importlib
import importlib.util
import subprocess
import sys
import unittest
from types import SimpleNamespace
import numpy as np

EXPERIMENT = Path(__file__).resolve().parents[1]
ROOT = EXPERIMENT.parents[3]
sys.path[:0] = [str(EXPERIMENT.parent), str(ROOT/'code/model')]
adapter = importlib.import_module('birth_count_choice.model.estate_observer_adapter')
distribution = importlib.import_module('birth_count_choice.model.engine.distribution')


class EstateObserverTests(unittest.TestCase):
    def test_toy_changes_only_estate_statistics_and_default_exact(self):
        legacy = importlib.import_module('production.engine.distribution')
        mass = np.ones((1,2,1,1,1,1,1))
        bp = np.array([7.,-40.]).reshape(mass.shape)
        P = SimpleNamespace(J=1,J_R=1,I=1,z_grid=np.ones(1),period_years=4.,
            income=np.ones((1,1)),H_own=np.array([100.]),psi=.06,
            bequest_net_of_selling_cost=True,estate_flow_net_of_selling_cost=True)
        args = (mass,mass,bp,P,np.array([40.]),np.ones(1))
        old, new = SimpleNamespace(),SimpleNamespace()
        legacy.add_aggregate_wealth_bequest_flow_moments(old,*args)
        adapter.aggregate_moments(distribution,legacy,new,*args)
        self.assertEqual(new.annual_bequest_flow,(7.+54.)/4.)
        for k in vars(old).keys()-{'annual_bequest_flow','annual_bequest_flow_to_aggregate_wealth'}:
            np.testing.assert_array_equal(getattr(new,k),getattr(old,k))
        P.bequest_net_of_selling_cost=P.estate_flow_net_of_selling_cost=False
        default=SimpleNamespace()
        adapter.aggregate_moments(distribution,legacy,default,*args)
        for k in vars(old):np.testing.assert_array_equal(getattr(default,k),getattr(old,k))

    def test_authenticated_helper_real_cached_net_flow_and_other_rows(self):
        # Legacy model imports must not contaminate reporting authentication in
        # the parent unittest process, which requires a fresh native import.
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--cached-check'],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def _check_authenticated_helper_real_cached_net_flow_and_other_rows(self):
        from birth_count_choice.model.storage import load_case
        case=ROOT/'output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T211930831130Z_6809d185'
        cached,_=load_case(case)
        P,sol=cached.P,cached.solution
        P._g_stay_distribution=sol.g_stay_distribution
        P._bp_pol_stay=sol.bp_pol_stay
        ev=SimpleNamespace(g_post_fertility=sol.g_beginning_distribution,g_current=sol.g,
            policy=SimpleNamespace(bp_pol=sol.bp_pol,price=sol.p_eq))
        path=ROOT/'tmp/e5f_overnight_local_20260927/portable/calibration_code_integration_20260927_v2/source/code/model/tools/e5f_initial_housing_observer.py'
        spec=importlib.util.spec_from_file_location('estate_test_actual_frozen_housing',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        rows={name:module._row(name,'test') for name in module.MOMENT_NAMES}
        ages=module.uniform_age_cell_overlap(P,76.,85.)
        module._wealth_diagnostics(rows,ev,P,cached.b_grid,ages)
        old=rows
        receipt=adapter.adapt_observer(module.observe_initial_housing_wealth,distribution)
        new={name:module._row(name,'test') for name in module.MOMENT_NAMES}
        module._wealth_diagnostics(new,ev,P,cached.b_grid,ages)
        key='annual_bequest_flow_to_aggregate_wealth'
        self.assertAlmostEqual(new[key]['model_value'],sol.annual_bequest_flow_to_aggregate_wealth,places=13)
        self.assertAlmostEqual(new[key]['model_value'],.006361342219724652,places=13)
        self.assertGreater(old[key]['model_value'],new[key]['model_value'])
        for name in old.keys()-{key}:self.assertEqual(new[name],old[name])
        self.assertEqual(new[key]['denominator'],old[key]['denominator'])
        self.assertEqual(receipt['original_file_sha256'],adapter.SOURCE_PIN)
        self.assertTrue(receipt['age_geometry_and_state_weights_unchanged'])
        print('cached estate-A observer: gross=',old[key]['model_value'],'net=',new[key]['model_value'])


if __name__=='__main__':
    if sys.argv[1:] == ['--cached-check']:
        EstateObserverTests()._check_authenticated_helper_real_cached_net_flow_and_other_rows()
    else:
        unittest.main()
