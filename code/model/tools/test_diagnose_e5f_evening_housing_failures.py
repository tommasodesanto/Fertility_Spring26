"""Pure Torch process/pin tests. No model import, model solve, or cluster submit.

Run in a small Torch Slurm allocation:
  python -m unittest -v test_diagnose_e5f_evening_housing_failures
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import diagnose_e5f_evening_housing_failures as diagnostic

TORCH_TEST = sys.platform == 'linux' and bool(os.environ.get('SLURM_JOB_ID'))


@unittest.skipUnless(TORCH_TEST, 'Torch Slurm only; no local process tests')
class Pins(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.contract=self.root/'contract.json';self.contract.write_text('{}')
        self.checkpoint=self.root/'checkpoint.json'
        record=dict(case='initial_0041_block',status='inadmissible',design='broad_coverage',lane='block',
                    error=dict(error='Initial housing equilibrium failed its unchanged strict gate',
                               traceback='initial = evaluate(float(initial_psi))',
                               context={'contract_sha256':diagnostic.digest(self.contract)}))
        self.checkpoint.write_text(json.dumps({'records':[record]}))
        def pin(path):return dict(path=str(path),sha256=diagnostic.digest(path))
        self.plan=dict(status='lead_approved',case_id=record['case'],arms=[1.,.5],case_cap_seconds=1.,
                       parent_cap_seconds=30.,global_end_epoch=time.time()+60.,
                       runner=pin(Path(diagnostic.__file__).resolve()),contract=pin(self.contract),
                       checkpoint=pin(self.checkpoint))
        self.path=self.root/'plan.json'
        self.args=SimpleNamespace(plan=self.path,plan_sha256='')
        self.save_plan()

    def save_plan(self):
        self.path.write_text(json.dumps(self.plan));self.args.plan_sha256=diagnostic.digest(self.path)

    def test_valid_small_contract_has_no_model_import(self):
        forbidden_before={k for k in sys.modules if k.startswith('intergen_')}
        p,_,r=diagnostic.authenticate(self.args)
        self.assertEqual(r['case'],'initial_0041_block')
        self.assertEqual(p['arms'],[1.,.5])
        self.assertEqual({k for k in sys.modules if k.startswith('intergen_')},forbidden_before)

    def test_changed_plan_refused(self):
        self.path.write_text(self.path.read_text()+' ')
        with self.assertRaisesRegex(RuntimeError,'Plan pin differs'):diagnostic.authenticate(self.args)

    def test_changed_contract_refused(self):
        self.contract.write_text('{"changed":true}')
        with self.assertRaisesRegex(RuntimeError,'Bad pin: contract'):diagnostic.authenticate(self.args)

    def test_changed_checkpoint_refused(self):
        self.checkpoint.write_text('{"records":[]}')
        with self.assertRaisesRegex(RuntimeError,'Bad pin: checkpoint'):diagnostic.authenticate(self.args)

    def test_unapproved_refused(self):
        self.plan['status']='awaiting_review';self.save_plan()
        with self.assertRaisesRegex(RuntimeError,'awaits explicit lead review'):diagnostic.authenticate(self.args)

    def test_budget_extension_refused(self):
        self.plan['case_cap_seconds']=1801;self.save_plan()
        with self.assertRaisesRegex(RuntimeError,'Budget enlarged'):diagnostic.authenticate(self.args)

    def test_ancestry_contract_identity_refused_even_if_rehashed(self):
        data=json.loads(self.checkpoint.read_text())
        data['records'][0]['error']['context']['contract_sha256']='0'*64
        self.checkpoint.write_text(json.dumps(data))
        self.plan['checkpoint']['sha256']=diagnostic.digest(self.checkpoint);self.save_plan()
        with self.assertRaisesRegex(RuntimeError,'not the original failed-case contract'):diagnostic.authenticate(self.args)


@unittest.skipUnless(TORCH_TEST, 'Torch Slurm only; no local process tests')
class ParentSupervision(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.args=SimpleNamespace(plan=self.root/'unused_plan.json',plan_sha256='synthetic',output=self.root/'output')
        self.plan=dict(case_cap_seconds=1.,parent_cap_seconds=30.,global_end_epoch=time.time()+60.,arms=[1.,.5])

    def test_parent_kills_descendant_after_leader_terminates(self):
        real_popen=subprocess.Popen
        children=[];groups=[]
        grandchild="import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)"
        def launch(cmd,**kwargs):
            self.assertTrue(kwargs['start_new_session'])
            for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
                self.assertEqual(kwargs['env'][name],'1')
            arm=cmd[cmd.index('--arm')+1];pidfile=self.root/f'grandchild_{arm}.pid'
            # Parent retains default SIGTERM behavior; its descendant ignores it.
            synthetic=("import subprocess,sys,time,pathlib; "
                       f"p=subprocess.Popen([sys.executable,'-c',{grandchild!r}]); "
                       f"pathlib.Path({str(pidfile)!r}).write_text(str(p.pid)); time.sleep(60)")
            proc=real_popen([sys.executable,'-c',synthetic],**kwargs)
            groups.append(proc.pid);children.append(pidfile);return proc
        try:
            with patch.object(diagnostic,'authenticate',return_value=(self.plan,{},{})),patch.object(diagnostic.subprocess,'Popen',side_effect=launch):
                diagnostic.parent(self.args)
            result=json.loads((self.args.output/'complete.json').read_text())
            self.assertEqual(len(result['cases']),2)
            self.assertTrue(all(r['status']=='censored_timeout' for r in result['cases']))
            self.assertTrue(all(r['returncode']==-signal.SIGTERM for r in result['cases']))
            self.assertLess(result['elapsed'],10.)
            for file in children:
                self.assertTrue(file.exists(),'Synthetic descendant never started')
                stat=Path('/proc')/file.read_text()/'stat'
                deadline=time.monotonic()+2.
                while stat.exists() and stat.read_text().split()[2]!='Z' and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue(not stat.exists() or stat.read_text().split()[2]=='Z','Descendant survived owned group cleanup')
        finally:
            for pid in groups:
                try:os.killpg(pid,signal.SIGKILL)
                except ProcessLookupError:pass

    def test_reserved_parent_deadline_prevents_dispatch(self):
        self.plan['global_end_epoch']=time.time()+1.
        with patch.object(diagnostic,'authenticate',return_value=(self.plan,{},{})),patch.object(diagnostic.subprocess,'Popen') as launch:
            diagnostic.parent(self.args)
        launch.assert_not_called()
        rows=json.loads((self.args.output/'complete.json').read_text())['cases']
        self.assertEqual([r['status'] for r in rows],['unrun_parent_budget']*2)

    def test_normal_child_finish_is_not_timeout(self):
        real_popen=subprocess.Popen
        def launch(cmd,**kwargs):return real_popen([sys.executable,'-c','raise SystemExit(0)'],**kwargs)
        with patch.object(diagnostic,'authenticate',return_value=(self.plan,{},{})),patch.object(diagnostic.subprocess,'Popen',side_effect=launch):
            diagnostic.parent(self.args)
        rows=json.loads((self.args.output/'complete.json').read_text())['cases']
        self.assertTrue(all(r['status']=='child_finished' and r['returncode']==0 for r in rows))


if __name__=='__main__':unittest.main()
