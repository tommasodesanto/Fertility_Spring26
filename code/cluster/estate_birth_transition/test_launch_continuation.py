"""Bounded launcher argument/pin checks; subprocess creation is always mocked."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock

sys.path.insert(0,str(Path(__file__).resolve().parent))
import launch_continuation as m


class LauncherTests(unittest.TestCase):
    def fixture(self, directory, mode='smoke', wall=21600):
        root=Path(directory)
        driver=root/'driver.py';driver.write_text('# test only\n')
        manifest=root/'manifest.json';manifest.write_text('{}')
        args=['--python',sys.executable,'--driver',str(driver),'--manifest',str(manifest),
            '--mode',mode,'--output',str(root/'launch'),'--wall-seconds',str(wall),'--memory-gib','24']
        return args,root

    def reject(self,args):
        with patch.object(m.subprocess,'Popen') as children,contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as exc:m.main(args)
        self.assertEqual(exc.exception.code,2);children.assert_not_called()

    def test_both_modes_reject_excess_wall_budget_without_children(self):
        for mode in ('smoke','run'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                args,root=self.fixture(d,mode,21601)
                if mode=='run':
                    receipt=root/'smoke.json';receipt.write_text('{}')
                    args+=['--smoke-receipt-pin',json.dumps(dict(path=str(receipt),sha256=m.sha(receipt)))]
                self.reject(args);self.assertFalse((root/'launch').exists())

    def test_zero_wall_budget_without_children(self):
        with tempfile.TemporaryDirectory() as d:
            args,_=self.fixture(d,wall=0);self.reject(args)

    def test_run_requires_smoke_pin_without_children(self):
        with tempfile.TemporaryDirectory() as d:
            args,_=self.fixture(d,mode='run');self.reject(args)

    def test_smoke_rejects_receipt_pin_without_children(self):
        with tempfile.TemporaryDirectory() as d:
            args,_=self.fixture(d);args+=['--smoke-receipt-pin','{}'];self.reject(args)

    def test_smoke_accepts_explicit_six_hour_bound_with_mocked_children(self):
        with tempfile.TemporaryDirectory() as d:
            args,root=self.fixture(d)
            with patch.object(m.subprocess,'Popen',return_value=MagicMock(pid=1234)) as children,contextlib.redirect_stdout(io.StringIO()):
                m.main(args)
            self.assertEqual(children.call_count,2)
            spec=json.loads((root/'launch/invocation.json').read_text())
            self.assertEqual(spec['wall_seconds'],21600);self.assertEqual(spec['memory_gib'],24)
            self.assertIsNone(spec['smoke_receipt_pin'])

    def test_receipt_pin_exact_fields_hash_and_content_remain_required(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'smoke.json';path.write_text('{}')
            pin=dict(path=str(path),sha256=m.sha(path))
            self.assertEqual(m.validate_pin(pin),pin)
            with self.assertRaisesRegex(ValueError,'exactly'):m.validate_pin(dict(pin,extra=True))
            path.write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError,'SHA-256 differs'):m.validate_pin(pin)
            args,_=self.fixture(d,mode='run');args+=['--smoke-receipt-pin',json.dumps(pin)]
            self.reject(args)

    def test_manager_rechecks_pins_before_any_child(self):
        for field in ('driver_sha256','smoke_receipt_pin'):
            with self.subTest(field=field),tempfile.TemporaryDirectory() as d:
                args,root=self.fixture(d)
                receipt=root/'smoke.json';receipt.write_text('{}')
                spec=dict(nonce='nonce',driver=str(root/'driver.py'),manifest=str(root/'manifest.json'),
                    driver_sha256=m.sha(root/'driver.py'),manifest_sha256=m.sha(root/'manifest.json'),
                    smoke_receipt_pin=dict(path=str(receipt),sha256=m.sha(receipt)))
                if field=='driver_sha256':spec[field]='changed'
                else:receipt.write_text('{"changed":true}')
                spec_path=root/'invocation.json';spec_path.write_text(json.dumps(spec))
                with patch.object(m.subprocess,'Popen') as children,patch.object(m.signal,'signal'):
                    m.manager(spec_path,'nonce')
                children.assert_not_called()
                terminal=json.loads((root/'manager_terminal.json').read_text())
                self.assertIsNone(terminal['worker_pid']);self.assertIn('manager_error',terminal['reason'])

    def test_memory_cap_and_self_test_without_children(self):
        with tempfile.TemporaryDirectory() as d:
            args,_=self.fixture(d);args[-1]='25';self.reject(args)
        with patch.object(m.subprocess,'Popen') as children,contextlib.redirect_stdout(io.StringIO()):m.self_test()
        children.assert_not_called()


if __name__=='__main__':unittest.main()
