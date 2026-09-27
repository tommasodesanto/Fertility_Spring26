import hashlib,json,tempfile,unittest
from pathlib import Path
import numpy as np
import e5f_recent_parent_flow_observer as observer
import e5f_current_transition_runtime as runtime
class DeadTailTests(unittest.TestCase):
    def fixture(self,mass=1e-14,value=-1e10):
        post=np.array([1.,mass]).reshape(2,1,1,1,1,1,1)
        loc=np.array([1.,0.]).reshape(2,1,1,1,1,1,1,1)
        values=np.array([1.,value]).reshape(post.shape)
        return loc,post,values
    def test_small_dead_mass_retained_exactly(self):
        a,b,v=self.fixture();original=[x.copy() for x in (a,b,v)]
        error,e=observer._audit_location_lottery(a,b,v,allow_dead_tail=True)
        self.assertEqual(error,1.);self.assertEqual(e['offending_mass'],1e-14);self.assertTrue(e['accepted_dead_tail'])
        for x,y in zip((a,b,v),original):np.testing.assert_array_equal(x,y)
    def test_total_mass_gate(self):
        with self.assertRaises(ValueError):observer._audit_location_lottery(*self.fixture(1.01e-12),allow_dead_tail=True)
    def test_small_live_invalid_rejected(self):
        with self.assertRaises(ValueError):observer._audit_location_lottery(*self.fixture(value=0.),allow_dead_tail=True)
    def test_default_still_rejects_dead_tail(self):
        with self.assertRaises(ValueError):observer._audit_location_lottery(*self.fixture())
    def test_default_valid_ignores_unneeded_value(self):
        a,b,v=self.fixture();a[:]=1
        self.assertEqual(observer._audit_location_lottery(a,b,None)[0],0.)
    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):observer._audit_location_lottery(*self.fixture(value=np.nan),allow_dead_tail=True)
    def test_frozen_measurement_reconstruction(self):
        main=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
        c=json.loads((main/'tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/production_contract.json').read_text())
        frozen=Path(c['source_root'])/'code/model/tools/e5f_recent_parent_flow_observer.py'
        runtime.verify_current_recent_observer(Path(observer.__file__),frozen)
        with tempfile.TemporaryDirectory() as d:
            changed=Path(d)/'observer.py';changed.write_text(Path(observer.__file__).read_text().replace('PRUNING_TOLERANCE = 1e-15','PRUNING_TOLERANCE = 1e-14'))
            with self.assertRaises(RuntimeError):runtime.verify_current_recent_observer(changed,frozen)
if __name__=='__main__':unittest.main()
