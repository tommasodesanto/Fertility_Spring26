"""Torch-only synthetic tests for the exact finite subprocess search loop."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np
import search
import worker


def configuration(root,behavior='linear'):
    names=list(search.PARAMETERS)
    center=dict(H0=6.,beta_annual=.96,chi=1.,first_birth_fixed_cost=.6,kappa_fert=.2,
        kappa_fert_continuation=.3,theta0=.12,delta_alpha_jump=.13,child_benefit_curvature=.06,tenure_choice_kappa=.012)
    bounds={k:[0.,10.] for k in names};bounds.update(beta_annual=[.94,.99],delta_alpha_jump=[0.,.25],child_benefit_curvature=[0.,.8],tenure_choice_kappa=[.001,.1])
    scale=search.scales(center,names)
    target={k:center[k]+float(scale[i])*.1 for i,k in enumerate(names)}
    spec=dict(initial_point=center,initial_psi=.13,source_fingerprint='synthetic_one',target_fingerprint='ten_weighted_targets',seed=2609281)
    return dict(schema='two_stream_search_v1',synthetic=True,parameters=names,scored_moments=[f'm{i}' for i in range(10)],bounds=bounds,
        lanes={'one_birth':spec,'two_birth':dict(spec,source_fingerprint='synthetic_two',seed=2609282)},
        worker_command=[sys.executable,str(Path(__file__).resolve()),'--fake-worker','--behavior',behavior,'--fixture',str(root/'fixture.json')],
        hard_end_epoch=time.time()+120,budget=dict(total_seconds=90,final_reserve_seconds=4,case_seconds=2,maximum_stationary_solves=8,max_evaluations=36),
        fixture=dict(center=center,scale=scale.tolist(),target=target))


def fake_worker(args):
    request=search.read(args.request);out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    fixture=search.read(args.fixture);names=list(search.PARAMETERS)
    identity={k:request[k] for k in ('candidate_id','lane','config_sha256','source_fingerprint','target_fingerprint')}
    behavior=args.behavior
    if behavior=='timeout' or (behavior=='two_timeouts' and request['role']!='initial_replay' and not request['role'].startswith('final_repeat')):
        process=subprocess.Popen([sys.executable,'-c','import time;time.sleep(120)'])
        search.write(out/'grandchild.json',dict(pid=process.pid));time.sleep(120)
    if behavior=='fatal_timeout':
        search.write(out/'FAILURE.json',dict(identity,status='fatal',authenticated=True,model_evaluations=8))
        time.sleep(120)
    if behavior=='mismatched_fatal_timeout':
        wrong_identity=dict(identity,candidate_id='wrong_candidate')
        search.write(out/'FAILURE.json',dict(wrong_identity,status='fatal',authenticated=True,model_evaluations=8))
        time.sleep(120)
    if behavior in ('fatal','inadmissible','budget_censored'):
        status='censored' if behavior=='budget_censored' else behavior
        search.write(out/'FAILURE.json',dict(identity,status=status,authenticated=True,model_evaluations=8));return 2
    if behavior=='incomplete_jac' and request['role']=='jac0_3':
        search.write(out/'FAILURE.json',dict(identity,status='inadmissible',authenticated=True));return 2
    point=request['point'];scale=np.array(fixture['scale'])
    target=fixture['target'];r=np.array([(point[k]-target[k])/scale[i] for i,k in enumerate(names)])
    # Each lane has its own normalization derivative; using the other lane's
    # derivative would predict the wrong starting benefit in the GN proposal.
    multiplier=1. if request['lane']=='one_birth' else -2.
    displacement=np.array([(point[k]-fixture['center'][k])/scale[i] for i,k in enumerate(names)])
    psi=.13+multiplier*.001*float(displacement.sum())
    case=out/'case';case.mkdir();artifact=case/'synthetic.json';search.write(artifact,dict(synthetic=True))
    checkpoint=case/'initial_state.pkl.gz';checkpoint.write_bytes(b'synthetic-checkpoint-v1\n')
    checkpoint_sha256=search.sha(checkpoint)
    search.write(case/'receipt.json',dict(identity,case_checkpoint_sha256=checkpoint_sha256))
    search.write(case/'scientific_identity.json',dict(identity,checkpoint_sha256=checkpoint_sha256))
    receipt=dict(identity,status='passed',loss=float(r@r),residuals=r.tolist(),point=point,psi=psi,case_path=str(case.resolve()),synthetic=True,
         model_evaluations=1,elapsed_seconds=0.,checkpoint_sha256=checkpoint_sha256,
         artifacts=[dict(path=str(artifact.resolve()),sha256=search.sha(artifact))])
    if behavior=='wrong_fingerprint':receipt['source_fingerprint']='wrong'
    search.write(out/'SUCCESS.json',receipt);return 0


class SearchTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def make(self,behavior='linear',lane='one_birth'):
        c=configuration(self.root,behavior);search.write(self.root/'fixture.json',c.pop('fixture'))
        return search.Stream(c,lane,self.root/('out_'+lane),'synthetic_config'),c

    def test_ten_coordinate_bounds_and_actual_steps(self):
        stream,c=self.make();center=dict(stream.spec['initial_point'])
        for i,k in enumerate(c['parameters']):center[k]=c['bounds'][k][i%2]
        scale,design=search.probe_design(center,c['parameters'],c['bounds'])
        self.assertEqual(len(design),10)
        for i,row in enumerate(design):
            search.check_point(row['point'],c['parameters'],c['bounds'])
            self.assertNotEqual(row['scaled_step'],0.)
            self.assertAlmostEqual(row['scaled_step'],(row['point'][row['parameter']]-center[row['parameter']])/scale[i])
        with self.assertRaises(RuntimeError):search.check_point(dict(center,extra=1),c['parameters'],c['bounds'])

    def test_exact_linear_loop_improves_and_keeps_final_slots(self):
        stream,c=self.make();result=stream.run()
        self.assertLess(result['selected']['loss'],.001)
        self.assertLessEqual(result['attempts'],36)
        self.assertEqual(result['repeat_count'],2)
        self.assertTrue(result['numerical_repeat_screens_passed'])
        self.assertEqual(stream.records[-1]['role'],'final_repeat1')
        self.assertEqual(stream.records[-2]['role'],'final_repeat0')
        first=search.read(stream.output/'jacobian_0.json')
        self.assertTrue(first['complete']);self.assertEqual(first['valid_probe_count'],10)
        np.testing.assert_allclose(first['matrix'],np.eye(10),atol=1e-10)
        self.assertEqual(first['diagnostics']['numerical_rank'],10)
        self.assertLessEqual(max(abs(x) for p in first['proposals'] for x in p['scaled_step']),1.+1e-12)
        self.assertEqual(first['center']['point'],stream.spec['initial_point'])

    def test_fresh_warm_derivative_per_lane(self):
        centers=[];designs=[];probes=[]
        for lane,mult in [('one_birth',1.),('two_birth',-2.)]:
            root=self.root/lane;root.mkdir();c=configuration(root);fixture=c.pop('fixture');search.write(root/'fixture.json',fixture)
            stream=search.Stream(c,lane,root/'out','synthetic_config')
            center=stream.evaluate(stream.spec['initial_point'],stream.spec['initial_psi'],'initial_replay')
            scale,design=search.probe_design(center['point'],stream.names,stream.bounds)
            cases=[stream.evaluate(p['point'],center['psi'],f'probe{i}') for i,p in enumerate(design)]
            J,dpsi=search.jacobian(center,cases,design)
            np.testing.assert_allclose(dpsi,np.full(10,.001*mult),atol=1e-11)
            proposals,info=search.gn_proposals(center,stream.names,stream.bounds,scale,J,dpsi)
            for p in proposals:
                expected=center['psi']+float(dpsi@p['scaled_step'])
                self.assertAlmostEqual(p['initial_psi'],expected,places=12)
            centers.append(dpsi)
        self.assertFalse(np.array_equal(*centers))

    def test_incomplete_probe_does_not_make_jacobian(self):
        stream,c=self.make('incomplete_jac');result=stream.run()
        first=search.read(stream.output/'jacobian_0.json')
        self.assertFalse(first['complete']);self.assertNotIn('matrix',first)
        self.assertFalse(any(r['role'].startswith('gn') for r in stream.records))
        self.assertTrue(any(r['role'].startswith('explore') for r in stream.records))
        self.assertEqual(result['repeat_count'],2)

    def test_fatal_stops_without_retry(self):
        stream,c=self.make('fatal');result=stream.run()
        self.assertEqual(result['attempts'],1)
        self.assertEqual(result['stop_reason'],'fatal_integrity_or_scientific_failure')

    def test_fingerprint_mismatch_stops(self):
        stream,c=self.make('wrong_fingerprint');result=stream.run()
        self.assertEqual(result['attempts'],1)
        self.assertEqual(stream.records[0]['status'],'fatal')

    def test_owned_timeout_kills_descendants(self):
        stream,c=self.make('timeout')
        result=stream.run();self.assertEqual(stream.records[0]['status'],'censored')
        child=search.read(stream.output/stream.records[0]['candidate_id']/'grandchild.json')['pid']
        stat=Path('/proc')/str(child)/'stat'
        # SIGKILL reaches the owned process group synchronously, but /proc may
        # briefly expose the child before the kernel reaps it.
        deadline=time.monotonic()+2
        while stat.exists() and stat.read_text().split()[2]!='Z' and time.monotonic()<deadline:
            time.sleep(.05)
        self.assertTrue(not stat.exists() or stat.read_text().split()[2]=='Z')

    def test_authenticated_fatal_written_before_timeout_stops_stream(self):
        stream,c=self.make('fatal_timeout');result=stream.run()
        self.assertEqual(result['attempts'],1)
        self.assertEqual(stream.records[0]['status'],'fatal')
        self.assertEqual(result['stop_reason'],'fatal_integrity_or_scientific_failure')

    def test_mismatched_failure_written_before_timeout_stops_stream(self):
        stream,c=self.make('mismatched_fatal_timeout');result=stream.run()
        self.assertEqual(result['attempts'],1)
        self.assertEqual(stream.records[0]['status'],'fatal')
        self.assertEqual(result['stop_reason'],'fatal_integrity_or_scientific_failure')

    def test_target_accounting_allows_one_ulp_but_rejects_economic_gap(self):
        expected=.1+.2
        self.assertTrue(worker.same_accounting_value(math.nextafter(expected,math.inf),expected))
        self.assertFalse(worker.same_accounting_value(expected+1e-6,expected))

    def test_two_consecutive_censors_stop_chain(self):
        stream,c=self.make('two_timeouts')
        result=stream.run()
        self.assertEqual(result['attempts'],5)
        self.assertEqual(result['stop_reason'],'two_consecutive_censored_or_inadmissible')
        self.assertEqual(sum(r['status']=='censored' for r in stream.records),2)
        self.assertEqual(result['repeat_count'],2)

    def test_authenticated_solve_cap_is_censor(self):
        stream,c=self.make('budget_censored');stream.run()
        self.assertEqual(stream.records[0]['status'],'censored')

    def test_actual_time_guard_and_final_reserve(self):
        stream,c=self.make();stream.cutoff=time.time()+.1
        self.assertFalse(stream.can_start())
        self.assertTrue(stream.can_start(final=True))
        stream.records=[{}]*34
        self.assertFalse(stream.can_start())
        self.assertTrue(stream.can_start(final=True))
        stream.records=[{}]*36;self.assertFalse(stream.can_start(final=True))
        stream.records=[];stream.end=time.time()+.1;self.assertFalse(stream.can_start(final=True))

    def test_rank_deficiency_is_diagnostic(self):
        stream,c=self.make();center=dict(point=stream.spec['initial_point'],psi=.13,residuals=[1.]*10)
        proposals,diag=search.gn_proposals(center,stream.names,stream.bounds,search.scales(center['point'],stream.names),np.zeros((10,10)),np.zeros(10))
        self.assertEqual(diag['numerical_rank'],0);self.assertEqual(proposals,[])

    def test_partial_rank_jacobians_remain_provisional(self):
        stream,c=self.make();center=dict(point=stream.spec['initial_point'],psi=.13,residuals=[1.]*10)
        for rank in range(1,10):
            J=np.zeros((10,10));J[np.arange(rank),np.arange(rank)]=1.
            _,diag=search.gn_proposals(center,stream.names,stream.bounds,search.scales(center['point'],stream.names),J,np.zeros(10))
            self.assertEqual(diag['numerical_rank'],rank)
            self.assertIn('not statistical',diag['interpretation'])

    def test_source_and_expiring_launch_approval(self):
        source=Path(search.__file__).resolve()
        config=dict(pins=[dict(path=str(source),sha256=search.sha(source))],synthetic=False,
                    lanes={'one_birth':dict(source_fingerprint='synthetic_one'),
                           'two_birth':dict(source_fingerprint='synthetic_two')})
        synthetic=self.root/'synthetic.json';search.write(synthetic,dict(status='passed',config_sha256='config',synthetic=True))
        integrations={}
        for lane in ('one_birth','two_birth'):
            path=self.root/f'integration_{lane}.json';search.write(path,dict(status='passed',config_sha256='config',lane=lane,real_model_evaluations=2,source_fingerprint='synthetic_'+lane.split('_')[0]))
            integrations[lane]=dict(path=str(path),sha256=search.sha(path))
        approval_path=self.root/'approval.json'
        approval=dict(schema='two_stream_launch_approval_v1',status='approved_two_stream_overnight',config_sha256='config',allowed_lanes=['one_birth'],
           not_before_epoch=time.time()-1,expires_epoch=time.time()+20,source_fingerprints={'one_birth':'synthetic_one','two_birth':'synthetic_two'},
           synthetic_receipt=dict(path=str(synthetic),sha256=search.sha(synthetic)),integration_receipts=integrations)
        search.write(approval_path,approval)
        with patch.dict(os.environ,dict(TWO_STREAM_LAUNCH_APPROVAL=str(approval_path),EXPECTED_TWO_STREAM_LAUNCH_APPROVAL_SHA256=search.sha(approval_path))):
            search.verify_launch(config,'one_birth','config')
            with self.assertRaises(RuntimeError):search.verify_launch(config,'two_birth','config')
            approval['expires_epoch']=time.time()-1;search.write(approval_path,approval)
            os.environ['EXPECTED_TWO_STREAM_LAUNCH_APPROVAL_SHA256']=search.sha(approval_path)
            with self.assertRaises(RuntimeError):search.verify_launch(config,'one_birth','config')
        config['pins'][0]['sha256']='bad'
        with self.assertRaises(RuntimeError):search.verify_launch(config,'one_birth','config')


if __name__=='__main__':
    if '--fake-worker' in sys.argv:
        parser=argparse.ArgumentParser();parser.add_argument('--fake-worker',action='store_true');parser.add_argument('--behavior');parser.add_argument('--fixture',type=Path);parser.add_argument('--request',type=Path);parser.add_argument('--output',type=Path)
        sys.exit(fake_worker(parser.parse_args()))
    search.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Synthetic tests belong on Torch')
    unittest.main()
