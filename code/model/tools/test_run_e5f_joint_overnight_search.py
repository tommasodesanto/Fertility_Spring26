"""Pure orchestration checks; no model or cluster work is launched."""
import copy
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest import mock

import run_e5f_joint_overnight_search as search


class OvernightTests(unittest.TestCase):
    def test_de_deterministic_bounded_all_dimensions_free(self):
        pop=[[i/64+j/1000 for j in range(11)] for i in range(32)]
        a=search.de_trial(pop,list(range(32)),random.Random(17))
        self.assertEqual(a,search.de_trial(pop,list(range(32)),random.Random(17)))
        self.assertTrue(all(0<=x<=1 for u in a for x in u))
        self.assertTrue(all(any(u[j]!=p[j] for u,p in zip(a,pop)) for j in range(11)))

    def test_bound_repair_uses_parent_not_clipped_endpoint(self):
        self.assertEqual(search.repair([.2,.8],[-2,3]),[.1,.9])

    def test_payload_removes_stale_loss_and_preserves_preference_change(self):
        domain=search.adapter.SEARCH_DOMAIN
        seed={'old_psi_child':.3,'best_candidate':{'theta':{'beta':.99,'psi_child':0},
              'new_psi_child':-.02,'transition_loss':1.},'loss':1.}
        unit=[.5]*11
        p=search.unit_to_payload(seed,domain,unit,'test')
        self.assertNotIn('loss',p);self.assertNotIn('transition_loss',p['best_candidate'])
        for d,u in zip(domain,unit):
            v=search.transform(u,d)
            if d['name']=='beta_annual':actual=p['best_candidate']['theta']['beta']**.25
            elif d['name']=='psi_child_change_2023':actual=p['best_candidate']['new_psi_child']-p['old_psi_child']
            else:actual=p['best_candidate']['theta'][d['name']]
            self.assertAlmostEqual(actual,v,places=13)
        with self.assertRaises(ValueError):search.unit_to_payload(seed,domain,[2]*11,'bad')

    def test_full_360_budget_reserves_reproductions(self):
        self.assertEqual(2+32+8*32+2*(22+12)+2,360)
        self.assertTrue(search.has_budget(346,12))
        self.assertFalse(search.has_budget(347,12))
        self.assertTrue(search.has_budget(358,2,repeats=0))

    def test_accept_only_actual_strict_improvements_in_matching_slots(self):
        p=[[.1]*11,[.2]*11];loss=[2.,3.]
        new=[({'id':2,'loss':4.},[.8]*11),({'id':1,'loss':1.},[.9]*11)]
        result,values=search.select_population(copy.deepcopy(p),loss[:],new)
        self.assertEqual(result,[[.9]*11,[.2]*11]);self.assertEqual(values,[1.,3.])
        with self.assertRaises(RuntimeError):search.select_population(p,loss,new[:1])

    def test_polish_compares_original_anchor_and_ranks_actual_gains(self):
        center=[.5]*11
        a=center[:];a[0]+=.0025
        b=center[:];b[1]-=.0025
        pairs=[({'label':'coordinate_0_plus','loss':8.},a),
               ({'label':'coordinate_1_minus','loss':7.},b)]
        d,order=search.polish_directions(center,10.,pairs)
        self.assertGreater(d[0],0);self.assertLess(d[1],0);self.assertEqual(order[:2],[1,0])
        vectors,_=search.joint_polish_vectors(center,d,order,{tuple(center)})
        self.assertTrue(any(u[0]>.5 and u[1]<.5 for u in vectors))
        self.assertEqual(len(vectors),len(set(map(tuple,vectors))))

    def test_exact_repeat_copies_original_generator_bytes(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);prior=root/'prior';prior.mkdir()
            center=prior/'center.json';center.write_text('{"original": 1.0000000000000002}\n')
            candidate=prior/'task_001';(candidate/'cases/task_001').mkdir(parents=True)
            for name,body in [('summary.json','{"best_candidate":{"candidate":"task_001"}}'),
                              ('target_fit_long.csv','x\n1\n'),('parameter_table.csv','x\n1\n')]:
                (candidate/name).write_text(body)
            (candidate/'cases/task_001/transition_path.csv').write_text('x\n1\n')
            original={'id':1,'center':'center.json','center_sha256':search.adapter.digest(center),
                      'panel_task_id':1,'panel_size':1,'panel_seed':7,'panel_design':'mixed','radius':.005,'output':'task_001'}
            plan=prior/'plan.json';plan.write_text('{}')
            run=search.Search.__new__(search.Search);run.root=root/'search';run.root.mkdir()
            run.c={'base_plan':{},'controller_sha256':'test'};run.finish_epoch=2000;run.search_epoch=1000
            row={'id':1,'plan':str(plan),'plan_sha256':'test','summary':str(candidate/'summary.json')}
            with mock.patch.object(search.adapter,'load_plan',return_value={'cases':[original]}):
                generated,_=run.new_plan('repeats',repeat=(row,[.5]*11))
            data=json.loads(generated.read_text())
            self.assertEqual(data['repeat_origin']['center_sha256'],search.adapter.digest(center))
            for case in data['cases']:
                self.assertEqual((generated.parent/case['center']).read_bytes(),center.read_bytes())
                self.assertEqual(case['panel_seed'],7)

    def test_failed_child_stops_before_any_queued_case_starts(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);plan=root/'plan.json';plan.write_text('{}')
            run=search.Search.__new__(search.Search)
            run.c={'max_workers':1};run.phase='';run.completed=0
            run.write_state=mock.Mock();run.kill_owned=mock.Mock();run.run_child=mock.Mock(side_effect=RuntimeError('failed solve'))
            cases=[{'id':i} for i in (1,2,3)]
            with mock.patch.object(search.adapter,'load_plan',return_value={'stage':'test','cases':cases}),\
                 mock.patch.object(search.adapter,'write_json'):
                run.root=root
                with self.assertRaisesRegex(RuntimeError,'failed solve'):run.run_batch(plan,'test')
            self.assertEqual(run.run_child.call_count,1)
            run.kill_owned.assert_called_once()


    def test_complete_search_flow_reaches_original_generator_repeats(self):
        self.exercise_complete_search(False)

    def test_recovery_probes_every_parameter_each_round_and_repeats_winner(self):
        self.exercise_complete_search(True)

    def exercise_complete_search(self, recovery):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw)
            run=search.Search.__new__(search.Search)
            run.c={'seed':123,'search_domain':search.adapter.SEARCH_DOMAIN,'output_root':str(root),'max_histories':360}
            if recovery:
                run.c.update(schema='e5f_joint_local_recovery_v1', max_histories=210, polish_rounds=6,
                             polish_radius=.00125, minimum_polish_radius=.0003125)
            run.root=root;run.completed=4 if recovery else 2;run.ledger=[];run.seed={'panel_design':{'unit_vector':[.5]*11}}
            run.best=({'loss':100.,'id':1,'summary':str(root/'seed/summary.json')},[.5]*11)
            run.require_smoke=mock.Mock();run.stage_fits=lambda n,repeats=False: search.has_budget(run.completed,n,maximum=run.c['max_histories'],repeats=0 if repeats else 2)
            phases=[]
            def evaluated(stage,vectors,labels):
                phases.append((stage,len(vectors)));answer=[]
                if recovery and stage.endswith('_coordinates'):
                    self.assertEqual(len(vectors),22)
                    for j in range(11):
                        for i in (2*j,2*j+1):
                            self.assertNotEqual(vectors[i][j],run.best[1][j])
                            self.assertEqual([v for k,v in enumerate(vectors[i]) if k!=j],
                                             [v for k,v in enumerate(run.best[1]) if k!=j])
                for i,(unit,label) in enumerate(zip(vectors,labels),1):
                    run.completed+=1
                    row={'id':i,'label':label,'loss':100.-run.completed,'summary':str(root/'selected/summary.json')}
                    pair=(row,list(unit));run.ledger.append(pair);answer.append(pair);run.best=pair
                return answer
            run.evaluate=evaluated
            selected_repeat=[]
            def make_repeat(stage,repeat):
                selected_repeat.append(copy.deepcopy(repeat));return root/'repeat_plan.json','hash'
            run.new_plan=make_repeat
            def repeats(*_):
                phases.append(('final_repeats',2));run.completed+=2;rows=[]
                for i in (1,2):
                    dest=root/f'repeat_{i}';dest.mkdir()
                    (dest/'case_receipt.json').write_text(json.dumps({'reference':{'exact_twelve_row_fit':True}}))
                    rows.append(({'summary':str(dest/'summary.json'),'id':i},[.5]*11))
                return rows
            run.run_batch=repeats;run.finish_report=mock.Mock()
            with mock.patch.object(search,'verify_graph_match',return_value=17):run.search()
            self.assertEqual(phases[0],('polish_1_coordinates',22) if recovery else ('initial_population',32))
            self.assertEqual(sum(name.startswith('generation_') for name,_ in phases),0 if recovery else 8)
            self.assertEqual(sum(name.endswith('_coordinates') for name,_ in phases),6 if recovery else 2)
            self.assertEqual(phases[-1],('final_repeats',2));self.assertLessEqual(run.completed,210 if recovery else 360)
            self.assertEqual(len(selected_repeat),1)
            run.finish_report.assert_called_once_with('complete_verified')

if __name__=='__main__':unittest.main()
