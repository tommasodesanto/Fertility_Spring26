import copy,unittest
from run_e5f_local_continuation import validate_plan,request
class TestPlan(unittest.TestCase):
    def setUp(self):
        self.obj={'parameter_restrictions':[{'parameter':'x','lower':0.,'upper':1.}]}
        self.plan=dict(workers=2,total_seconds=1740,search_seconds=1080,case_seconds=480,point={'x':.5},cases=[dict(id='trial',point={'x':.6})],authorization='Author daytime local authorization',supersedes_overnight_deadline=True)
    def test_valid(self):validate_plan(self.plan,self.obj)
    def test_bad_budgets(self):
        for key,value in [('workers',3),('total_seconds',1800),('search_seconds',1200),('case_seconds',481),('supersedes_overnight_deadline',False)]:
            with self.subTest(key=key),self.assertRaises(AssertionError):
                p=copy.deepcopy(self.plan);p[key]=value;validate_plan(p,self.obj)
    def test_invalid_candidates(self):
        for cases in [[],[dict(id='../escape',point={'x':.2})],[dict(id='repeat_1',point={'x':.2})],[dict(id='trial',point={'x':1.1})],[dict(id='trial',point={'x':float('nan')})],[dict(id='trial',point={'x':.1})]*2]:
            with self.subTest(cases=cases),self.assertRaises(AssertionError):
                p=copy.deepcopy(self.plan);p['cases']=cases;validate_plan(p,self.obj)
    def test_authenticated_request(self):
        class Driver:
            canon=staticmethod(lambda p:'pointpin')
            candidate_id=staticmethod(lambda c,p:'candidatepin')
            norm_inputs=staticmethod(lambda c:{'fixed':True})
        c=dict(files={'driver':{'sha256':'sourcepin'}},objective={'sha256':'targetpin'})
        r=request(Driver,c,'contractpin',dict(id='trial',point={'x':.2}),'search',100.,False)
        self.assertEqual(r['context']['source_sha256'],'sourcepin');self.assertEqual(r['deadline_epoch'],100.)
        self.assertEqual(r['scientific_candidate_id'],'candidatepin');self.assertFalse(r['graphs'])
if __name__=='__main__':unittest.main()
