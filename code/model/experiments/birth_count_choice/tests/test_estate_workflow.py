"""Zero-solve workflow/contract propagation and full rescoring tests."""
import copy,json,subprocess,sys,tempfile,time,unittest
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(BASE))
from model import calibration
from model.estate_contract import apply_experiment_flags,experiment_flags,contract,rescore_rows,make_evaluator
from model.inputs import load_inputs

class EstateWorkflowTests(unittest.TestCase):
    def test_whitelist_and_caps(self):
        for cap in (1,2,3):
            P,_=load_inputs();apply_experiment_flags(P,experiment_flags(cap))
            self.assertEqual(P.birth_count_choice_cap,cap)
            self.assertTrue(P.bequest_net_of_selling_cost and P.estate_flow_net_of_selling_cost)
        for cap in (0,4,True):
            with self.assertRaises(ValueError):experiment_flags(cap)
        with self.assertRaises(ValueError):apply_experiment_flags(None,{'something_else':True})
        flags=experiment_flags(1);flags['estate_flow_net_of_selling_cost']=False
        with self.assertRaises(ValueError):apply_experiment_flags(None,flags)

    def test_rescore_only_changes_wealth_target(self):
        anchor=json.loads(calibration.ANCHOR.read_text())
        old=anchor['target_fit'];new,residual=rescore_rows(old)
        self.assertEqual(len(new),14);self.assertEqual(residual.shape,(10,))
        for a,b in zip(old,new):
            if a['moment']!='wealth_earnings':self.assertEqual(a,b)
            else:
                self.assertEqual(b['target'],'4.45838713455674')
                self.assertEqual(a['model'],b['model']);self.assertEqual(a['weight'],b['weight'])
        self.assertEqual(contract()[3]['beta_annual'],(.93,.99))

    def test_adapter_carries_flags_and_supplied_arrays(self):
        anchor=json.loads(calibration.ANCHOR.read_text())
        source=calibration.ANCHOR.parent/'native_postcheck/selected_postcheck/phase_b_ge/selected_root'
        P,grid=load_inputs();P.adapter_probe=np.array([123.]);calls=[]
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);report=out/'report';report.mkdir()
            for name in ('target_fit.csv','parameters.csv'):(report/name).write_bytes((source/name).read_bytes())
            closure=json.loads((source/'closure.json').read_text())
            def mock(Q,b,**kwargs):
                self.assertEqual(Q.birth_count_choice_cap,1)
                self.assertTrue(Q.bequest_net_of_selling_cost and Q.estate_flow_net_of_selling_cost)
                self.assertIsNot(Q.adapter_probe,P.adapter_probe)
                np.testing.assert_array_equal(b,grid)
                self.assertEqual(kwargs['closure'],'population_one');calls.append(kwargs)
                return dict(report_directory=str(report),closure=closure,price=anchor['selected']['price'],lifecycle_solves=2)
            end=time.time()+120
            tp,wp=contract()[1:3]
            evaluate=make_evaluator(out,'none',P,grid,end,birth_cap=1,solver=mock,target_fingerprint=tp,weight_fingerprint=wp)
            receipt=evaluate('mock',anchor['selected']['parameters'],end)
            self.assertEqual(len(calls),1);self.assertEqual(len(receipt['target_fit']),14)
            self.assertEqual(len(receipt['parameter_table']),31)
            beta=next(r for r in receipt['parameter_table'] if r['parameter']=='beta_annual')
            self.assertEqual(beta['lower'],'.93')
            self.assertEqual(receipt['target_fingerprint'],tp)
            self.assertTrue((report/'target_fit_new_contract.csv').is_file())
            self.assertEqual((report/'target_fit.csv').read_bytes(),(source/'target_fit.csv').read_bytes())
            json.dumps(receipt,allow_nan=False)
        self.assertFalse(hasattr(P,'birth_count_choice_cap'))

    def test_cli_is_inert(self):
        result=subprocess.run([sys.executable,str(BASE/'run_estate.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0);self.assertIn('--birth-cap',result.stdout)

if __name__=='__main__':unittest.main()
