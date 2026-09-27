"""No model imports/solves: authentication and economic-contract guards."""
import tempfile
import unittest
import copy
import numpy as np
from pathlib import Path
from types import SimpleNamespace as NS
import e5f_current_transition_runtime as r

class RuntimeContract(unittest.TestCase):
    def test_economic_contract_retains_fixed_psi_and_tax(self):
        result=r.economic_contract(NS(psi_child=.14,tau_pay=.08,property_tax_lump_sum_transfer=0.))
        self.assertFalse(result['preferences']['renormalize_fertility'])
        self.assertEqual(result['paygo']['payroll_tax'],.08)
        self.assertEqual(result['natural_credit']['classification'],'experimental_not_enabled')
        self.assertEqual(result['population_and_geography']['classification'],'outstanding_for_transition')

    def test_rebates_are_not_inherited_silently(self):
        with self.assertRaises(ValueError):
            r.economic_contract(NS(psi_child=.14,tau_pay=.08,property_tax_lump_sum_transfer=.1))

    def test_wrong_source_rejected(self):
        with self.assertRaises(RuntimeError):
            r.require_current_model(NS(__file__='/tmp/frozen/solver.py'))

    def test_hash_and_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'evidence.json'
            r.write(path,{'status':'no solve'})
            self.assertEqual(r.read(path),{'status':'no solve'})
            self.assertEqual(len(r.sha(path)),64)

    def test_only_reviewed_reporting_extension_is_allowed(self):
        frozen=r.PORTABLE/'tools_v4/e5f_calibration_runtime.py'
        current=r.ROOT/'code/model/tools/e5f_calibration_runtime.py'
        r.verify_current_reporter(current,frozen)
        with tempfile.TemporaryDirectory() as directory:
            bad=Path(directory)/'bad.py';bad.write_text(current.read_text()+'\n# unreviewed extra change\n')
            with self.assertRaises(RuntimeError): r.verify_current_reporter(bad,frozen)

    def test_array_census_retains_inactive_differences(self):
        old=dict(b_grid=np.array([0.,1.]), shared=NS(v=np.ones(2)),
                 solution=NS(V=np.array([2.,3.])), evaluation=NS(g_current=np.array([1.,0.])),
                 stationary_g_pre=np.array([1.,0.]))
        new=copy.deepcopy(old);new['solution'].V[1]=99.
        report=r.compare_arrays(old,new)
        row=report['arrays']['solution.V']
        self.assertEqual(row['max_abs'],96.)
        self.assertEqual(row['mass_weighted_abs'],0.)
        self.assertFalse(row['exact'])
        self.assertEqual(report['status'],'review_required')

if __name__=='__main__': unittest.main()
