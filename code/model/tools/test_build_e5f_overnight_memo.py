"""Synthetic reporting tests; none of these values is a model result."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import build_e5f_overnight_memo as memo


def save(path, value):
    path.write_text(json.dumps(value))


def csvfile(path, values):
    with path.open('w', newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)


def fixture(root, completed=True):
    objective={'target_rows':[
        {'restriction_id':name,'target':2.1 if name=='initial_normalization' else 1.,
         'actual_weight':None if name=='initial_normalization' else 100. if name=='mean_rooms' else 1.}
        for name in memo.MOMENTS],
        'parameter_restrictions':[{'parameter':name,'lower':0.,'upper':1.}
                                  for name in memo.PARAMETERS if name!='psi_child']}
    save(root/'objective.json',objective)
    index={'primary_objective':'objective.json','cases':[],
           'status_note':'SYNTHETIC QA FIXTURE - not an empirical or model result.',
           'source_manifest_sha256':'synthetic-shared-source',
           'expected_workers':{'cluster':24,'local':10}}
    if completed:
        for name,weighting,changed,gap in [('primary_case','primary','mean_rooms',.1),('identity_case','identity','cps_childlessness',.5)]:
            case=root/name;case.mkdir()
            fits=[]
            for target in objective['target_rows']:
                moment=target['restriction_id'];g=gap if moment==changed else 0.
                weight=target['actual_weight']
                if weighting=='identity' and weight is not None:weight=1.
                fits.append(dict(moment=moment,target=target['target'],model=target['target']+g,gap=g,
                                 weight='' if weight is None else weight,
                                 loss_contribution='' if weight is None else weight*g*g))
            csvfile(case/'target_fit.csv',fits)
            parameters=[dict(parameter=p,estimate=.5,lower=0.,upper=1.,near_bound=False,status='free')
                        for p in memo.PARAMETERS if p!='psi_child']
            parameters.append(dict(parameter='psi_child',estimate=.2,lower='',upper='',near_bound='',status='normalized'))
            csvfile(case/'parameters.csv',parameters)
            save(case/'receipt.json',dict(status=memo.SUCCESS,loss=sum(r['loss_contribution'] or 0 for r in fits),
                normalization={'psi_child':.2},objective_stationary_solves=2,objective_stationary_solve_seconds=240,
                source_manifest_sha256='synthetic-shared-source',case_checkpoint_sha256='synthetic-checkpoint'))
            save(case/'stationary_solves.json',[{'status':'completed'},{'status':'completed'}])
            index['cases'].append({'path':name,'weighting':weighting,'status':'success'})
    save(root/'index.json',index)
    return root/'index.json'


class MemoTests(unittest.TestCase):
    def test_early_weight_experiment_checks_every_weight(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=fixture(root);index=memo.read(path)
            entry=index['cases'][0];entry['weighting']='early_fertility_3000'
            case=root/'primary_case';fits=list(memo.rows(case/'target_fit.csv','moment').values())
            for row in fits:
                if row['moment']=='early_fertility':row['weight']='3000.0'
            csvfile(case/'target_fit.csv',fits);save(path,index)
            targets,restrictions=memo.target_spec(memo.read(root/'objective.json'))
            checked=memo.validate_case(dict(entry,path=str(case)),targets,restrictions)
            self.assertEqual(checked['weighting'],'early_fertility_3000')
            changed=next(r for r in fits if r['moment']=='mean_rooms')
            changed['weight']='1';changed['loss_contribution']=float(changed['gap'])**2
            csvfile(case/'target_fit.csv',fits)
            receipt=memo.read(case/'receipt.json');receipt['loss']=sum(float(r['loss_contribution'] or 0) for r in fits)
            save(case/'receipt.json',receipt)
            with self.assertRaisesRegex(ValueError,'unexpected weights'):
                memo.validate_case(dict(entry,path=str(case)),targets,restrictions)

    def test_pending_keeps_every_target_and_parameter_without_values(self):
        with tempfile.TemporaryDirectory() as folder:
            result=memo.collect(fixture(Path(folder),False))
            self.assertEqual(result['status'],'no_completed_primary_point')
            self.assertEqual(len(result['target_fit']),14)
            self.assertEqual(len(result['parameters']),10)
            self.assertTrue(all(r['model'] is None for r in result['target_fit']))

    def test_common_weights_and_separate_primary_selection(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);result=memo.collect(fixture(root))
            self.assertTrue(result['selected']['path'].endswith('primary_case'))
            self.assertAlmostEqual(result['selected']['primary_loss'],1.)
            self.assertAlmostEqual(result['groups']['identity']['best_under_primary_weights'],.25)
            self.assertEqual(result['counts']['success'],2)
            self.assertEqual(len(result['parameters']),10)

    def test_mismatched_target_is_excluded_not_rescored(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);index=fixture(root)
            fits=list(memo.rows(root/'identity_case/target_fit.csv','moment').values())
            fits[1]['target']='999';csvfile(root/'identity_case/target_fit.csv',fits)
            result=memo.collect(index)
            self.assertEqual(result['counts']['collector_rejected'],1)
            self.assertEqual(result['groups']['identity']['completed'],0)

    def test_source_mismatch_cannot_win_common_comparison(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);index=fixture(root)
            receipt=memo.read(root/'identity_case/receipt.json');receipt['source_manifest_sha256']='other'
            save(root/'identity_case/receipt.json',receipt)
            result=memo.collect(index)
            self.assertEqual(result['counts']['comparison_excluded'],1)
            self.assertEqual(result['groups']['identity']['completed'],0)

    def test_fixed_economic_change_cannot_be_treated_as_weight_experiment(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);index=fixture(root)
            receipt=memo.read(root/'identity_case/receipt.json')
            receipt['retained_owner_grid']=[1.,2.,3.]
            save(root/'identity_case/receipt.json',receipt)
            result=memo.collect(index)
            self.assertEqual(result['counts']['comparison_excluded'],1)
            self.assertEqual(result['groups']['identity']['completed'],0)

    def test_incomplete_solve_ledger_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);index=fixture(root)
            save(root/'primary_case/stationary_solves.json',[{'status':'failed'}])
            result=memo.collect(index)
            self.assertIsNone(result['selected'])
            self.assertEqual(result['counts']['collector_rejected'],1)

    def test_controller_records_and_duplicate_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=fixture(root);index=memo.read(path)
            controller=root/'controller';controller.mkdir()
            save(controller/'checkpoint.json',{'records':[{'case':'a','case_path':str(root/'primary_case'),
                                                         'status':'success','loss':1.}]})
            index['controllers']=[{'path':'controller','weighting':'primary'}];save(path,index)
            result=memo.collect(path)
            self.assertEqual(result['counts']['success'],2)
            self.assertEqual(len(result['controller_states']),1)

    def test_complete_supporting_outputs(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);index=fixture(root)
            memo.build(index,root/'report',pdf=False)
            self.assertEqual(len(memo.rows(root/'report/target_fit.csv','moment')),14)
            self.assertEqual(len(memo.rows(root/'report/parameters.csv','parameter')),10)
            self.assertIn('not statistical z-scores',(root/'report/fit_overview.svg').read_text())


if __name__=='__main__':unittest.main()
