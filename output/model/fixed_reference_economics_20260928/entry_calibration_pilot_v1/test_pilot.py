"""Pure input and mocked-controller checks: zero native model evaluations."""
import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
import inputs,runner
class PilotChecks(unittest.TestCase):
    def test_three_laws(self):
        P,g=inputs.proposal();means=[]
        for arm in inputs.ARMS:
            Q,r=inputs.entry(P,g,arm);self.assertEqual(r['clipped_draw_mass'],0)
            np.testing.assert_allclose(Q.fixed_reference_entry_conditional.sum(0),1,atol=2e-15)
            if arm!='empirical_credit':self.assertEqual(r['negative_wealth_share'],0)
            if arm=='zero_wealth':self.assertEqual(r['mean_wealth'],0)
            else:means.append(r['mean_wealth'])
        self.assertAlmostEqual(*means,places=12)
    def test_all_nine_bindings(self):
        P,g=inputs.proposal();seed,bounds,_=inputs.seed_and_bounds()
        for key in inputs.PARAMETERS:
            point=dict(seed);point[key]+=runner.step_sizes(seed,bounds)[key];Q=inputs.bind(P,point,bounds)
            if key=='beta_annual':self.assertEqual(Q.beta,point[key]**P.period_years);self.assertEqual(Q.rho,1/Q.beta-1)
            else:self.assertEqual(getattr(Q,key),point[key])
            self.assertEqual(Q.eps_fert,Q.kappa_fert)
            np.testing.assert_array_equal(Q.H0,P.H0);self.assertEqual(Q.psi_child,P.psi_child)
    def test_target_contract(self):
        self.assertEqual(len(runner.PLAN['target_contract']),14)
        self.assertEqual(sum(x['role']=='scored' for x in runner.PLAN['target_contract']),10)
    def test_mock_exact_loop(self):
        seed,bounds,_=inputs.seed_and_bounds()
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);result=runner.search(out,seed,bounds,runner.mocked_evaluator(out,seed,bounds),runner.time.time()+3600,mock=True)
            self.assertEqual(result['completed_full_ge'],15);self.assertEqual(result['selected_repeats'],2);self.assertEqual(result['lifecycle_solves'],0)
            self.assertEqual(result['identification']['rank'],9)
    def test_failed_probe_never_zero_derivative(self):
        seed,bounds,_=inputs.seed_and_bounds()
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);good=runner.mocked_evaluator(out,seed,bounds)
            def evaluator(label,p):
                return dict(status='inadmissible_numerical',lifecycle_solves=0) if 'probe' in label else good(label,p)
            result=runner.search(out,seed,bounds,evaluator,runner.time.time()+3600,mock=True)
            self.assertEqual(result['identification']['status'],'incomplete_jacobian')
            self.assertEqual(result['completed_full_ge'],5)
    def test_baseline_failure_blocks(self):
        seed,bounds,_=inputs.seed_and_bounds()
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(RuntimeError,'Baseline full GE failed'):
                runner.search(Path(temp),seed,bounds,lambda *_:dict(status='inadmissible_numerical'),runner.time.time()+3600,mock=True)
    def test_nonnegative_credit_phaseB(self):
        P,_=inputs.proposal()
        for arm in inputs.ARMS:
            with tempfile.TemporaryDirectory() as temp:
                result=runner.phase_b_mock_checks(Path(temp),P,arm,runner.time.time()+3600)
                self.assertEqual(result['lifecycle_solves'],0)
                self.assertEqual(result['rejected_actual_parameter_drifts'],31)
    def test_rank_deficient_skips_GN(self):
        seed,bounds,_=inputs.seed_and_bounds()
        with tempfile.TemporaryDirectory() as temp:
            def evaluator(label,p):return dict(status='passed',residual=[1.]*10,lifecycle_solves=0)
            result=runner.search(Path(temp),seed,bounds,evaluator,runner.time.time()+3600,mock=True)
            self.assertEqual(result['identification']['status'],'underidentified_local_jacobian')
            self.assertEqual(result['completed_full_ge'],13)
            cases=json.loads((Path(temp)/'cases.json').read_text())
            self.assertFalse(any(x['kind']=='damped_Gauss_Newton' for x in cases))
if __name__=='__main__':unittest.main()
