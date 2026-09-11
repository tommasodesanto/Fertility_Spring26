import math
import tempfile
import json
from pathlib import Path
import unittest
import numpy as np
import run_search as s

class SearchTests(unittest.TestCase):
    def test_only_known_equilibrium_failure_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'raw').mkdir();(p/'preflight.json').write_text('{}')
            f={'error_type':'RuntimeError','phase':'stationary_equilibrium','error':'Initial housing equilibrium failed its unchanged strict gate'}
            (p/'raw/failure.json').write_text(json.dumps(f))
            self.assertEqual(s.failure_status(p),'rejected_equilibrium')
            f['error']='Fingerprint mismatch';(p/'raw/failure.json').write_text(json.dumps(f))
            self.assertEqual(s.failure_status(p),'failed')
    def test_missing_preflight_is_always_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(s.failure_status(tmp),'failed')
    def test_parameter_transform_roundtrip(self):
        for n,v in [('beta_annual',.995),('H0',7.5),('theta0',.44)]:
            self.assertAlmostEqual(s.inverse(n,s.transform(n,v)),v)
    def test_bounds_handle_reversed_beta(self):
        names=['beta_annual','theta0'];r={'beta_annual':{'lower':.94,'upper':.9995},'theta0':{'lower':0,'upper':8}}
        x=np.array([s.transform('beta_annual',.995),math.log(.4)])
        lo,hi=s.bounds(names,r,x,10)
        self.assertAlmostEqual(s.inverse(names[0],x[0]+lo[0]),.9995)
        self.assertAlmostEqual(s.inverse(names[0],x[0]+hi[0]),.94)
        self.assertTrue(np.all(lo<0));self.assertTrue(np.all(hi>0))
    def test_residual_keeps_all_scored_rows_only(self):
        score={'target_fit':[{'scored':False},{'scored':True,'gap':3,'actual_weight':4},{'scored':True,'gap':-1,'actual_weight':9}]}
        np.testing.assert_array_equal(s.residual(score),[6,-3])
    def test_exact_batch_loop_writes_each_completion(self):
        completed=[]
        result=s.batch([1,2,3],lambda i:dict(case_id=i,status='verified',loss=i*i),completed.append,2)
        self.assertEqual(sorted(x['case_id'] for x in completed),[1,2,3])
        self.assertEqual(len(result),3)
    def test_failed_case_is_preserved(self):
        completed=[]
        result=s.batch([1,2],lambda i:dict(case_id=i,status='failed' if i==2 else 'verified'),completed.append,2)
        self.assertEqual(sum(r['status']=='failed' for r in result),1)
        self.assertEqual(len(completed),2)
    def test_worker_exception_propagates(self):
        def worker(i):raise ValueError('Fingerprint mismatch')
        with self.assertRaisesRegex(ValueError,'Fingerprint'):
            s.batch([1],worker,lambda r:None,1)

if __name__=='__main__':unittest.main()
