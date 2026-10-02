"""Zero-solve checks against the running broad-array's actual slot-0 contract schema."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
os.environ['PURCHASE_PACKET_ROOT'] = str(PACKET)
spec = importlib.util.spec_from_file_location('purchase_scan_remote', HERE / 'scan_remote.py')
scan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan)
DESIGN = json.loads((PACKET / 'broader_regions_v1/design.json').read_text())
FIXTURE = HERE / 'testdata/region_slot0_search_contract.json'


class RegionReceiptTest(unittest.TestCase):
    def test_running_claimed_and_transient_slot(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            (base / 'source').mkdir()
            shutil.copy2(PACKET / 'broader_regions_v1/design.json', base / 'source/design.json')
            root = base / 'results/slot_0'
            root.mkdir(parents=True)
            shutil.copy2(FIXTURE, root / 'search_region_contract.json')
            original_root = scan.REGION_ROOT
            original_censor = os.environ.get('PURCHASE_REVIEWED_CENSOR')
            scan.REGION_ROOT = base
            os.environ['PURCHASE_REVIEWED_CENSOR'] = str(PACKET / 'restart_controller_v2/controller.py')
            try:
                ready = scan.region_provenance(root, 0, DESIGN)
                self.assertEqual((ready['original_chain_index'], ready['region_slot']), (0, 0))
                (root / 'postcheck').mkdir()
                (root / 'postcheck/completed.json').write_text('{"status":"selected_numerically_verified"}\n')
                with self.assertRaisesRegex(ValueError, 'Unbound broad-region'):
                    scan.region_provenance(root, 0, DESIGN)
                search = json.loads(FIXTURE.read_text())
                post = dict(search, stage='postcheck', optimizer_source_sha256=None)
                (root / 'postcheck_region_contract.json').write_text(json.dumps(post))
                (root / 'search').mkdir()
                shutil.copy2(HERE / 'testdata/region_slot0_native_search_contract.json',
                             root / 'search/search_contract.json')
                shutil.copy2(HERE / 'testdata/region_slot0_native_input_contract.json',
                             root / 'search/input_contract.json')
                (root / 'search/search_completed.json').write_text('{"objective_calls":80}\n')
                pending = scan.region_provenance(root, 0, DESIGN)
                self.assertEqual(pending['region_state'], 'awaiting_terminal')
                (root / 'launcher_terminal.json').write_text(json.dumps(dict(
                    slot=0, start_epoch=search['deadline_epoch']-4500,
                    deadline_epoch=search['deadline_epoch'], absolute_deadline_epoch=1790932500,
                    maximum_objective_calls=80, final_reserve_seconds=900, exit_code=1)))
                with self.assertRaisesRegex(ValueError, 'lacks successful'):
                    scan.region_provenance(root, 0, DESIGN)
            finally:
                scan.REGION_ROOT = original_root
                if original_censor is None:
                    os.environ.pop('PURCHASE_REVIEWED_CENSOR', None)
                else:
                    os.environ['PURCHASE_REVIEWED_CENSOR'] = original_censor


if __name__ == '__main__':
    unittest.main()
