"""Controller, distribution and binding tests; no native lifecycle solves."""
import copy,json,tempfile,time,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import inputs,runner
LANE='nonnegative_mean_120x9'
class RoundChecks(unittest.TestCase):
    def setUp(self):self.seed,self.bounds,_=inputs.seed_and_bounds(LANE)
    def test_lane_seeds_and_provenance(self):
        self.assertEqual(inputs.LANES[LANE]['seed'],inputs.LANES['nonnegative_mean_160x15']['seed'])
        for lane,cfg in inputs.LANES.items():
            inputs.check_point(cfg['seed'],self.bounds)
            self.assertEqual(cfg['seed_provenance']['full_ge_repeats'],2)
            self.assertEqual(len(cfg['seed_provenance']['file_sha256']),4)
    def test_actual_grids_and_entry_laws(self):
        means=[]
        for lane,cfg in inputs.LANES.items():
            P,g=inputs.proposal(lane);self.assertEqual([P.Nb,P.Nz],cfg['dimensions'])
            self.assertFalse(any(k.startswith('_') for k in vars(P)))
            Q,r=inputs.entry(P,g,cfg['arm']);self.assertEqual(r['clipped_draw_mass'],0)
            np.testing.assert_allclose(Q.fixed_reference_entry_conditional.sum(0),1,rtol=0,atol=2e-15)
            if cfg['arm']=='nonnegative_mean':self.assertEqual(r['negative_wealth_share'],0)
            means.append(r['mean_wealth'])
        self.assertAlmostEqual(means[0],means[1],places=12)
        # Income-discretization changes can alter annual entrant income across grids.
        self.assertTrue(all(np.isfinite(means)))
    def test_bind_all_nine_fixed_objects(self):
        P,_=inputs.proposal(LANE)
        for key in inputs.PARAMETERS:
            point=dict(self.seed);point[key]+=runner.step_sizes(self.seed,self.bounds)[key]
            Q=inputs.bind(P,point,self.bounds)
            self.assertEqual(Q.eps_fert,Q.kappa_fert)
            self.assertEqual(Q.beta,point['beta_annual']**P.period_years)
            self.assertEqual(Q.rho,1/Q.beta-1)
            np.testing.assert_array_equal(Q.H0,P.H0);self.assertEqual(Q.psi_child,P.psi_child)
            if key!='beta_annual':self.assertEqual(getattr(Q,key),point[key])
    def test_target_count_and_dimensions(self):
        self.assertEqual(len(runner.PLAN['target_contract']),14)
        self.assertEqual(sum(r['role']=='scored' for r in runner.PLAN['target_contract']),10)
        for dims in [(120,9),(160,15)]:
            x=runner.expected_parameters(self.seed,dims)
            self.assertEqual([x['wealth_grid_nodes'],x['income_states']],list(dims))
    def test_exact_loop_each_lane_and_center_updates(self):
        for lane,cfg in inputs.LANES.items():
            seed,bounds,_=inputs.seed_and_bounds(lane)
            with tempfile.TemporaryDirectory() as temp:
                out=Path(temp);result=runner.search(out,seed,bounds,runner.mocked_evaluator(out,seed,bounds,lane),time.time()+14400,mock=True,lane=lane)
                n=min(runner.PLAN['maximum_rounds'][cfg['size']],int(np.ceil((runner.PLAN['maximum_full_ge'][cfg['size']]-4)/12)))
                self.assertEqual(result['rounds_completed'],n)
                self.assertEqual(result['completed_full_ge'],runner.PLAN['maximum_full_ge'][cfg['size']])
                self.assertEqual(result['selected_repeats'],2)
                self.assertEqual(result['lifecycle_solves'],0)
                rounds=json.loads((out/'rounds.json').read_text())
                self.assertNotEqual(rounds[0]['center_parameters'],rounds[1]['center_parameters'])
                self.assertTrue(all(r['rank']==9 for r in rounds[:-1]));self.assertEqual(rounds[-1]['status'],'incomplete_jacobian')
    def test_failed_probe_not_zero_derivative(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);good=runner.mocked_evaluator(out,self.seed,self.bounds,LANE)
            def ev(label,p,end):return dict(status='inadmissible_numerical',lifecycle_solves=0) if 'probe' in label else good(label,p,end)
            r=runner.search(out,self.seed,self.bounds,ev,time.time()+14400,mock=True)
            self.assertEqual(r['search_stop_reason'],'incomplete_jacobian');self.assertEqual(r['completed_full_ge'],5)
            self.assertNotIn('rank',r['identification']);self.assertEqual(r['selected_repeats'],2)
    def test_rank_deficient_skips_GN(self):
        with tempfile.TemporaryDirectory() as temp:
            r=runner.search(Path(temp),self.seed,self.bounds,lambda *_:dict(status='passed',residual=[1.]*10,lifecycle_solves=0),time.time()+14400,mock=True)
            self.assertEqual(r['search_stop_reason'],'rank_deficient');self.assertEqual(r['completed_full_ge'],13)
    def test_failed_baseline_and_repeat_halt(self):
        for where in ['baseline','baseline_repeat']:
            with tempfile.TemporaryDirectory() as temp:
                out=Path(temp);good=runner.mocked_evaluator(out,self.seed,self.bounds,LANE)
                def ev(label,p,end):return dict(status='inadmissible_numerical') if label=='000_baseline' and where=='baseline' or label=='001_baseline_repeat' and where=='baseline_repeat' else good(label,p,end)
                with self.assertRaisesRegex(RuntimeError,'Baseline'):runner.search(out,self.seed,self.bounds,ev,time.time()+14400,mock=True)
    def test_selected_repeat_mismatch_halts(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);good=runner.mocked_evaluator(out,self.seed,self.bounds,LANE)
            def ev(label,p,end):
                r=good(label,p,end)
                if 'selected_repeat' in label:r['residual'][0]+=1
                return r
            with self.assertRaisesRegex(RuntimeError,'Selected full GE repeat failed'):runner.search(out,self.seed,self.bounds,ev,time.time()+14400,mock=True)
    def test_accounting_failure_propagates(self):
        with tempfile.TemporaryDirectory() as temp:
            def ev(*_):raise RuntimeError('Accounting drift')
            with self.assertRaisesRegex(RuntimeError,'Accounting drift'):runner.search(Path(temp),self.seed,self.bounds,ev,time.time()+14400,mock=True)
    def test_budget_reserve_and_guard_classifier(self):
        for lane in inputs.LANES:
            self.assertGreaterEqual(runner.reserve_seconds(lane,800),710+2*1.3*800)
        now=time.time();budget=SimpleNamespace(deadline_epoch=now+650,stage_deadline_seconds=300)
        self.assertTrue(runner.is_search_budget_exit(RuntimeError('No time reserve for selected reporting and exact repeat'),budget,now+14400))
        self.assertFalse(runner.is_search_budget_exit(RuntimeError('Accounting drift'),budget,now+14400))
        self.assertFalse(runner.is_search_budget_exit(RuntimeError('No time reserve for selected reporting and exact repeat'),budget,budget.deadline_epoch))
    def test_reserved_search_exit_still_repeats_best(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);good=runner.mocked_evaluator(out,self.seed,self.bounds,LANE);deadline=time.time()+14400;ends=[]
            def ev(label,p,end):
                ends.append((label,end))
                if 'probe' in label:return dict(status='budget_exhausted',lifecycle_solves=0)
                return good(label,p,end)
            with patch.object(runner,'compare_repeated',return_value={'status':'mock_comparison'}):r=runner.search(out,self.seed,self.bounds,ev,deadline,mock=False)
            self.assertEqual(r['selected_repeats'],2)
            self.assertTrue(all(end<deadline for label,end in ends if 'selected_repeat' not in label))
            self.assertTrue(all(end==deadline for label,end in ends if 'selected_repeat' in label))
    def test_native_phaseB_mock_binding_31_each_grid(self):
        for lane in inputs.LANES:
            P,_=inputs.proposal(lane)
            with tempfile.TemporaryDirectory() as temp:
                r=runner.phase_b_mock_checks(Path(temp),P,lane,time.time()+300)
                self.assertEqual(r['rejected_actual_parameter_drifts'],31);self.assertEqual(r['lifecycle_solves'],0)
if __name__=='__main__':unittest.main()
