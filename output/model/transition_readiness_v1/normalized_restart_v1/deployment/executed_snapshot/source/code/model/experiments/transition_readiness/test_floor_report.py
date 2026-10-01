"""Saved-packet validation tests; no native runtime or model is constructed."""
import json
from pathlib import Path
import tempfile
import types
import unittest
import floor_report as report

class ReportTests(unittest.TestCase):
    def test_input_hash_mismatch_fails_before_json_or_pickle_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'record.json'; path.write_text('not JSON')
            with self.assertRaisesRegex(RuntimeError, 'Input hash differs'):
                report.load_record(path, '0'*64)

    def test_packet_hash_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'record.json'; packet = Path(tmp)/'packet'; packet.write_bytes(b'bad')
            path.write_text(json.dumps(dict(rows=[dict(period=0,calendar_year=2007)],
                diagnostic_packets=[dict(path=str(packet),sha256='0'*64,period=0,dated_rent=1.)],
                market_residual=[0],fiscal_residual=[0],gates=dict(mass=True))))
            with self.assertRaisesRegex(RuntimeError,'Input hash differs'):
                report.load_record(path,report.sha(path))

    def test_original_date_set_enforced(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'record.json'
            path.write_text(json.dumps(dict(rows=[{}]*6,diagnostic_packets=[dict(period=0),dict(period=2),dict(period=5)])))
            with self.assertRaisesRegex(RuntimeError,'first/middle/last'):
                report.load_record(path,report.sha(path))

    def test_extra_png_rejected_and_hashes_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'date_000'/'standard_diagnostics';folder.mkdir(parents=True)
            names=[f'{i:02d}.png' for i in range(17)]
            for name in names:(folder/name).write_bytes(name.encode())
            receipts=dict(sampled_dates=[dict(period=0,plots={name:report.sha(folder/name) for name in names})])
            packets=[dict(period=0)];rows=[dict(calendar_year=2007)]
            self.assertEqual(report.verify_outputs(tmp,receipts,names,packets,rows)[0]['calendar_year'],2007)
            (folder/'extra.png').write_bytes(b'extra')
            with self.assertRaisesRegex(RuntimeError,'17-plot'):
                report.verify_outputs(tmp,receipts,names,packets,rows)
            (folder/'extra.png').unlink();(folder/names[0]).write_bytes(b'changed')
            with self.assertRaisesRegex(RuntimeError,'PNG hash'):
                report.verify_outputs(tmp,receipts,names,packets,rows)

    def test_reporting_solve_guard_restores_actual_callable(self):
        original=lambda: 'not called'
        runtime=types.SimpleNamespace(model=types.SimpleNamespace(solve_bellman_full_markov_income=original))
        with report.forbid_solves(runtime):
            with self.assertRaisesRegex(RuntimeError,'forbidden'):
                runtime.model.solve_bellman_full_markov_income()
        self.assertIs(runtime.model.solve_bellman_full_markov_income,original)

if __name__=='__main__':unittest.main()
