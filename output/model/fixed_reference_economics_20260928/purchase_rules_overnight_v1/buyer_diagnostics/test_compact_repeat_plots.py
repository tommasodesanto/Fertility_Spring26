"""Root-only plot contract: compact selected repeat has no second plot set."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_selected import require_root_plots


class CompactRepeatPlotTest(unittest.TestCase):
    def test_root_17_repeat_0_passes_root_16_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root = base / 'selected_root'
            repeat = base / 'selected_repeat'
            (root / 'standard_diagnostics').mkdir(parents=True)
            repeat.mkdir()
            for index in range(17):
                (root / 'standard_diagnostics' / f'{index:02d}.png').touch()
            self.assertEqual(len(list(repeat.rglob('*.png'))), 0)
            require_root_plots(root)
            (root / 'standard_diagnostics/16.png').unlink()
            with self.assertRaisesRegex(RuntimeError, '17-plot'):
                require_root_plots(root)


if __name__ == '__main__':
    unittest.main()
