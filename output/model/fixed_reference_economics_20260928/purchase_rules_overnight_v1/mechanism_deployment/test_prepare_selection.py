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
            finally:
                selection.ORIGINS.clear()
                selection.ORIGINS.update(original)


if __name__ == '__main__':
    unittest.main()
