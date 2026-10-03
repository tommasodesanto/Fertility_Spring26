#!/usr/bin/env python3
"""Compatibility entry for the archived June diagnostic runner.

Edit historical settings in calibration_archive/model_legacy_20261003/run_intergen_model.py.
Current production runs use code/model/run_model.py.
"""
from pathlib import Path as _Path

_archived_runner = (
    _Path(__file__).resolve().parents[2]
    / "calibration_archive/model_legacy_20261003/run_intergen_model.py"
)
# Execute in this namespace so historical Spyder globals remain visible. Import
# keeps the historical __name__ guard; direct execution still invokes main().
__file__ = str(_archived_runner)
exec(compile(_archived_runner.read_bytes(), __file__, "exec"), globals())
