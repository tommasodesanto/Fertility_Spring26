"""Subprocess lifecycle checks; full-stage scientific smoke runs on Torch."""
import tempfile
import threading
import sys
import unittest
from pathlib import Path
from run_e5f_simple_fertility_search import run_case_process

class ProcessLifecycle(unittest.TestCase):
    def call(self, code, timeout):
        with tempfile.TemporaryDirectory() as tmp:
            active=set()
            log=Path(tmp)/'case.log'
            result=run_case_process([sys.executable,'-c',code],log,timeout,active,threading.Lock())
            self.assertEqual(active,set())
            return result,log.read_text()
    def test_success_writes_case_log(self):
        result,log=self.call("print('receipt')",3)
        self.assertEqual(result,dict(returncode=0,timeout=False))
        self.assertEqual(log.strip(),'receipt')
    def test_failure_is_preserved(self):
        result,_=self.call('raise SystemExit(7)',3)
        self.assertEqual(result,dict(returncode=7,timeout=False))
    def test_timeout_kills_and_reaps_process(self):
        result,_=self.call('import time; time.sleep(60)',.1)
        self.assertEqual(result,dict(returncode=124,timeout=True))

if __name__=='__main__':unittest.main()
