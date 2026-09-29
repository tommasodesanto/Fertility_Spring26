"""Torch-only synthetic checks of the exact restricted controller loop."""
from __future__ import annotations
import argparse,csv,math,sys,tempfile,time,unittest
from pathlib import Path
import numpy as np
import search

def fixture(root,behavior='ok'):
    names=list(search.PARAMETERS)
    point=dict(H0=6.,beta_annual=.97,chi=1.,kappa_fert=.2,kappa_fert_continuation=.3,
        theta0=.12,delta_alpha_jump=.13,child_benefit_curvature=.06,tenure_choice_kappa=.012)
    bounds={name:[0.,10.] for name in names}
    bounds.update(beta_annual=[.94,.99],delta_alpha_jump=[0.,.25],child_benefit_curvature=[0.,.8],tenure_choice_kappa=[.001,.1])
    scale=search.scales(point,names)
    target={name:point[name]+.1*scale[i] for i,name in enumerate(names)}
    fit=root/'target_fit.csv'
    with fit.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['moment','weight','gap']);writer.writeheader()
        for i in range(10):writer.writerow(dict(moment=f'm{i}',weight=1,gap=-.1 if i<9 else 0.))
    old=dict(center=point,target=target,scale=scale.tolist(),behavior=behavior)
    search.write(root/'fixture.json',old)
    lane=dict(anchor_point=point,initial_point=point,anchor_fixed_cost=.35270914196085973,initial_psi=.12,
        anchor_case=str(root),anchor_loss=.09,source_fingerprint='source',target_fingerprint='targets')
    return dict(schema='zero_first_birth_cost_v1',synthetic=True,parameters=names,scored_moments=[f'm{i}' for i in range(10)],
        bounds=bounds,fixed_parameters={'first_birth_fixed_cost':0.0},lanes={'zero_cost':lane},
        worker_command=[sys.executable,str(Path(__file__).resolve()),'--fake-worker','--fixture',str(root/'fixture.json')],
        hard_end_epoch=time.time()+120,budget=dict(total_seconds=100,final_reserve_seconds=5,case_seconds=2,
            maximum_stationary_solves=8,max_evaluations=16))

def fake_worker(args):
    req=search.read(args.request);fixture=search.read(args.fixture);out=args.output;out.mkdir()
    identity={name:req[name] for name in ('candidate_id','lane','config_sha256','source_fingerprint','target_fingerprint')}
    role=req['role'];behavior=fixture['behavior']
    if behavior=='zero_censor' and role=='zero_cost_center':
        search.write(out/'FAILURE.json',dict(identity,status='censored',authenticated=True,model_evaluations=8));return 2
    if behavior=='bad_probe' and role=='jac0_2':
        search.write(out/'FAILURE.json',dict(identity,status='inadmissible',authenticated=True));return 2
    if behavior=='anchor_mismatch' and role=='anchor_replay':offset=.2
    else:offset=0.
    names=list(search.PARAMETERS);point=req['point'];scale=fixture['scale']
    residuals=[(point[name]-fixture['target'][name])/scale[i]+offset for i,name in enumerate(names)]+[0.]
    case=out/'case';case.mkdir();artifact=case/'synthetic.json';search.write(artifact,dict(synthetic=True))
    checkpoint=case/'initial_state.pkl.gz';checkpoint.write_bytes(b'synthetic\n');digest=search.sha(checkpoint)
    search.write(case/'receipt.json',dict(identity,case_checkpoint_sha256=digest))
    search.write(case/'scientific_identity.json',dict(identity,checkpoint_sha256=digest))
    success=dict(identity,status='passed',loss=float(np.dot(residuals,residuals)),residuals=residuals,
        point=point,fixed_parameters=req['fixed_parameters'],psi=.12,case_path=str(case.resolve()),
        model_evaluations=1,elapsed_seconds=0.,checkpoint_sha256=digest,
        artifacts=[dict(path=str(artifact.resolve()),sha256=search.sha(artifact))])
    search.write(out/'SUCCESS.json',success);return 0

class ExactLoop(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def run_case(self,behavior):
        c=fixture(self.root,behavior)
        stream=search.Stream(c,'zero_cost',self.root/'run','synthetic')
        return stream,stream.run()
    def test_complete_sixteen_case_loop(self):
        stream,final=self.run_case('ok')
        self.assertEqual(final['attempts'],16)
        self.assertEqual(final['repeat_count'],2)
        self.assertTrue(final['numerical_repeat_screens_passed'])
        self.assertEqual([r['role'] for r in stream.records[:2]],['anchor_replay','zero_cost_center'])
        self.assertEqual(stream.records[0]['result']['fixed_parameters']['first_birth_fixed_cost'],.35270914196085973)
        self.assertTrue(all(r['result']['fixed_parameters']['first_birth_fixed_cost']==0 for r in stream.records[1:] if r['status']=='success'))
        jac=search.read(stream.output/'jacobian_0.json')
        self.assertEqual(jac['valid_probe_count'],9);self.assertEqual(jac['diagnostics']['numerical_rank'],9)
        np.testing.assert_allclose(jac['matrix'][:9],np.eye(9),atol=1e-10)
        self.assertEqual(final['selected']['point'].keys(),stream.spec['initial_point'].keys())
    def test_zero_center_censor_stops_before_refit(self):
        stream,final=self.run_case('zero_censor')
        self.assertEqual(final['attempts'],2);self.assertEqual(final['stop_reason'],'zero_cost_center_unavailable')
        self.assertFalse((stream.output/'jacobian_0.json').exists())
    def test_anchor_mismatch_stops_before_cost_change(self):
        stream,final=self.run_case('anchor_mismatch')
        self.assertEqual(final['attempts'],1);self.assertEqual(final['stop_reason'],'anchor_replay_mismatch')
    def test_failed_probe_never_imputed(self):
        stream,final=self.run_case('bad_probe')
        jac=search.read(stream.output/'jacobian_0.json')
        self.assertFalse(jac['complete']);self.assertNotIn('matrix',jac)
        self.assertFalse(any(r['role'].startswith('gn') for r in stream.records))
    def test_free_count_and_fixed_cost(self):
        c=fixture(self.root)
        self.assertEqual(len(c['parameters']),9);self.assertNotIn('first_birth_fixed_cost',c['parameters'])
        with self.assertRaises(RuntimeError):search.check_point(dict(c['lanes']['zero_cost']['initial_point'],first_birth_fixed_cost=0),c['parameters'],c['bounds'])

if __name__=='__main__':
    if '--fake-worker' in sys.argv:
        p=argparse.ArgumentParser();p.add_argument('--fake-worker',action='store_true');p.add_argument('--fixture',type=Path)
        p.add_argument('--request',type=Path);p.add_argument('--output',type=Path)
        sys.exit(fake_worker(p.parse_args()))
    import os
    search.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ExactLoop)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():sys.exit(1)
    here=Path(__file__).resolve().parent
    files={name:search.sha(here/name) for name in ('worker.py','search.py','prepare.py','tests_search.py','run.sh')}
    search.write(here/'TESTS.json',dict(status='passed',test_count=result.testsRun,synthetic=True,
        model_imports=0,model_solves=0,source_hashes=files))
