import unittest, csv, json
from pathlib import Path
import numpy as np
from inspect_e5f_saved_households import weighted_quantile,sample_index
class SavedHouseholdsTest(unittest.TestCase):
    def test_zero_atom_is_not_interpolated(self):
        self.assertEqual(weighted_quantile(np.array([-1,0,3]),np.array([.2,.7,.1]),[.1,.5,.95]),[-1,0,3])
    def test_quantile_is_order_independent(self):
        self.assertEqual(weighted_quantile(np.array([3,-1,0]),np.array([.1,.2,.7]),[.1,.5,.95]),[-1,0,3])
    def test_sampler_preserves_joint_support(self):
        p=np.zeros((2,3));p[1,2]=1
        self.assertEqual(sample_index(p,np.random.default_rng(1)),(1,2))
    def test_sampler_rejects_mass_loss(self):
        with self.assertRaises(AssertionError):sample_index(np.array([.1,.1]),np.random.default_rng(1))
    def test_saved_native_factorization_and_paths(self):
        root=Path(__file__).resolve().parents[3]/'output/model/daytime_calibration_20260927/households'
        if not (root/'simulation_verification.json').exists():self.skipTest('Run authenticated saved-state audit first')
        receipt=json.loads((root/'simulation_verification.json').read_text())
        self.assertEqual(receipt['agents'],3)
        self.assertEqual(receipt['native_solves'],0)
        self.assertLess(max(receipt['transition_factorization_l1']),1e-11)
        with (root/'three_household_lives.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),receipt['row_count'])
        self.assertEqual({r['agent'] for r in rows},{'1','2','3'})
        for agent in ['1','2','3']:
            rr=[r for r in rows if r['agent']==agent]
            self.assertTrue(np.all(np.diff([float(r['age']) for r in rr])==4))
            self.assertTrue(np.all(np.diff([int(r['children_capped_three']) for r in rr])>=0))
            self.assertTrue(all(float(r['consumption'])>0 for r in rr))
if __name__=='__main__':unittest.main()
