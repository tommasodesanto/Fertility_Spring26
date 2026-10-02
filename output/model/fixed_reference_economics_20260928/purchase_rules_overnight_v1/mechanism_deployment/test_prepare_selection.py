"""Zero-solve routing checks for all collector winner origins."""
from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).with_name('prepare_selection.py')
spec = importlib.util.spec_from_file_location('prepare_selection', MODULE)
selection = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selection)


class PhysicalSourceTest(unittest.TestCase):
    def test_four_origins_preserve_chain_and_source(self):
        with tempfile.TemporaryDirectory() as temp:
            original = selection.ORIGINS.copy()
            try:
                for index, origin in enumerate(original):
                    if origin == 'torch_regions':
                        continue
                    selection.ORIGINS[origin] = str(Path(temp) / origin)
                    chain = 48 + index
                    root = Path(temp) / origin / f'chain_{chain}'
                    selected = dict(origin=origin, chain=chain, remote_root=str(root),
                                    remote_report=str(root / 'postcheck' / selection.REPORT),
                                    remote_arrays=str(root / 'postcheck' / selection.ARRAY))
                    if origin.endswith('restart'):
                        selected['parent_remote_root'] = str(Path(temp) / 'parent' / f'chain_{chain}')
                    self.assertEqual(selection.physical_source(selected), (origin, str(root)))
                    wrong = dict(selected, remote_root=str(Path(temp) / 'torch' / f'chain_{chain}'))
                    if origin != 'torch':
                        with self.assertRaises(RuntimeError):
                            selection.physical_source(wrong)
                    with self.assertRaises(RuntimeError):
                        selection.physical_source(dict(selected, chain=chain+1))
                with self.assertRaises(RuntimeError):
                    selection.physical_source(dict(selected, origin='invented'))
                selection.ORIGINS['torch_regions'] = str(Path(temp) / 'torch_regions')
                slot = 8
                chain = 24
                root = Path(temp) / 'torch_regions' / f'slot_{slot}'
                chosen = dict(origin='torch_regions', source_run='regions', arm='quarter',
                              chain=chain, original_chain_index=chain, region_slot=slot,
                              region_state='ready', region_design_sha256=selection.sha(
                                  selection.PACKET / 'broader_regions_v1/design.json'),
                              region_search_contract_sha256='a'*64,
                              region_postcheck_contract_sha256='b'*64,
                              region_launcher_terminal_sha256='c'*64,
                              remote_root=str(root), remote_report=str(root / 'postcheck' / selection.REPORT),
                              remote_arrays=str(root / 'postcheck' / selection.ARRAY))
                self.assertEqual(selection.physical_source(chosen), ('torch_regions', str(root)))
                for drift in (dict(chain=8), dict(region_slot=9), dict(region_state='awaiting_terminal'),
                              dict(remote_root=str(root.parent / 'chain_24'))):
                    with self.assertRaises(RuntimeError):
                        selection.physical_source(dict(chosen, **drift))
            finally:
                selection.ORIGINS.clear()
                selection.ORIGINS.update(original)


if __name__ == '__main__':
    unittest.main()
