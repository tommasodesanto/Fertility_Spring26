import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest
import numpy as np
import run_e5f_current_transition_smoke as s

class SmokeTests(unittest.TestCase):
    def test_prices_change_and_constant_control(self):
        p=s.trial_prices(1.,.8,'natural')
        self.assertTrue(.8<p[1]<p[0]<1.)
        np.testing.assert_equal(s.trial_prices(1.,1.,'reference'),[1.,1.])
        with self.assertRaises(ValueError):s.trial_prices(1.,1.,'natural')

    def test_exact_comparison_native_queues(self):
        def obj(x=0.):
            queue=NS(due_in_16=(.5,)*3,due_in_20=(.5,)*4)
            return NS(prices=np.array([1.,2.]),rents=np.ones(2),values=[np.array([3.])],
                      terminal_state=NS(g_pre=np.array([1.+x]),scheduled_entries=queue,scheduled_raw_entries=queue),
                      rows=[dict(moment=2.,kind='diagnostic')])
        self.assertTrue(all(v==0 for v in s.exact_difference(obj(),obj()).values()))
        self.assertGreater(s.exact_difference(obj(),obj(.1))['terminal_g'],0)

    def test_accounting_gates_explicit(self):
        a=dict(budget={'budget_excess_mass':0.},purchase={'maximum_occupied_transaction_wealth_error':0.,'transaction_outside_grid_mass':0.},
               estate={'status':'funded','audit_id':'estate_funded_dated_entry_provisional_net_v1','estate':{'totals':{'net_negative':0.}}})
        self.assertTrue(all(s.audit_gates(a).values()))
        a['estate']['audit_id']='stationary'
        self.assertFalse(s.audit_gates(a)['estate_next_cohort'])
        a['purchase']['transaction_outside_grid_mass']=1e-6
        self.assertFalse(s.audit_gates(a)['transaction_outside_grid_mass'])

    def test_approval_rejects_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory);cp=d/'checkpoint';cp.write_bytes(b'checkpoint')
            source=d/'source.py';source.write_text('x=1')
            receipt=d/'receipt.json';s.write(receipt,{'status':'native_reference_replay_pass'})
            approval=d/'approval.json';a=dict(approved=True,current_source_files={str(source):s.sha(source)},
                reference_checkpoint_sha256=s.sha(cp),baseline_receipt_path=str(receipt),baseline_receipt_sha256=s.sha(receipt))
            s.write(approval,a);self.assertTrue(s.validate_approval(approval,cp)['approved'])
            source.write_text('x=2')
            with self.assertRaises(ValueError):s.validate_approval(approval,cp)

    def test_terminal_approval_and_immutable_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory);case=d/'case';case.mkdir()
            checkpoint=case/'initial_state.pkl.gz';checkpoint.write_bytes(b'native')
            receipt=case/'receipt.json';s.write(receipt,{'status':'native_terminal_replay_verified'})
            core=d/'intergen_eqscale_seq_optimized'/'solver.py';core.parent.mkdir();core.write_text('x=1')
            pins={str(core):s.sha(core)}
            baseline=d/'baseline.json';s.write(baseline,{'current_source_files':pins})
            approval=d/'terminal.json';a=dict(approved=True,status='native_terminal_replay_verified',
                current_source_files=pins,terminal_checkpoint_sha256=s.sha(checkpoint),
                terminal_receipt_path=str(receipt),terminal_receipt_sha256=s.sha(receipt))
            s.write(approval,a);self.assertTrue(s.validate_terminal_approval(approval,case,baseline)['approved'])
            s.archive_sources(d/'archive',pins)
            manifest=s.read(d/'archive/source_snapshot_manifest.json')
            self.assertEqual(s.sha(manifest[str(core)]['snapshot']),s.sha(core))
            a['status']='unverified';s.write(approval,a)
            with self.assertRaises(ValueError):s.validate_terminal_approval(approval,case,baseline)

    def test_nonfinite_replay_rejected(self):
        q=NS(due_in_16=(0.,)*3,due_in_20=(0.,)*4)
        x=NS(prices=np.ones(2),rents=np.ones(2),values=[np.ones(2)],
             terminal_state=NS(g_pre=np.ones(1),scheduled_entries=q,scheduled_raw_entries=q),rows=[dict(x=float('nan'))])
        with self.assertRaises(ValueError):s.exact_difference(x,x)

if __name__=='__main__':unittest.main()
