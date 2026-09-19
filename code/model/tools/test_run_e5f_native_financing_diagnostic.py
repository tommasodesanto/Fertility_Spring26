import sys, unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_native_financing_diagnostic as d

def policy(value=1.): return SimpleNamespace(**{n:np.array([value]) for n in d.POLICY_NAMES})
class TestNativeFinancing(unittest.TestCase):
 def test_missing_arrays_blocks_execution(self):
  with self.assertRaisesRegex(ValueError,'missing mandatory'): d.policy_arrays(SimpleNamespace(V=np.ones(1)))
 def test_exact_baseline_mismatch_blocks_treatment(self):
  with self.assertRaises(AssertionError): d.compare_exact(policy(),policy(1+2e-10))
 def test_mass_mapping_gate_uses_same_pre_mass(self):
  ev=SimpleNamespace(g_pre=np.array([1.]),g_post_fertility=np.array([1.]),g_current=np.array([.9]),births=np.array([0.]),policy=policy())
  with self.assertRaisesRegex(ValueError,'mass gate'): d.gates(ev)
 def test_probability_gate(self):
  p=policy();p.fert_probs=np.array([1.01]);ev=SimpleNamespace(g_pre=np.array([1.]),g_post_fertility=np.array([1.]),g_current=np.array([1.]),births=np.array([0.]),policy=p)
  with self.assertRaisesRegex(ValueError,'probability'): d.gates(ev)
 def test_contract_hash_fail_closed(self):
  with self.assertRaisesRegex(ValueError,'checkpoint contract'): d.validate_contract(Path(__file__),Path(__file__).parent,Path(__file__).parent)
 def test_cases_have_requested_contrasts(self):
  self.assertEqual(d.CASES['baseline'],(.8,0.));self.assertEqual(d.CASES['mortgage_only'],(1.,0.));self.assertEqual(d.CASES['unsecured_only'],(.8,5.));self.assertEqual(d.CASES['both'],(1.,5.))
if __name__=='__main__': unittest.main()
