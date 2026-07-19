"""Stable test entry point that does not depend on pytest collection."""

from __future__ import annotations

import importlib
import inspect
import os
import sys
import unittest
from pathlib import Path


for thread_variable in (
    "NUMBA_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
):
    os.environ.setdefault(thread_variable, "1")


def main() -> int:
    tests_dir = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(str(tests_dir))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        return 1

    direct_count = 0
    for path in sorted(tests_dir.glob("test_*.py")):
        module = importlib.import_module(
            f"intergen_housing_fertility_optimized.tests.{path.stem}"
        )
        for name, function in sorted(vars(module).items()):
            if (
                name.startswith("test_")
                and inspect.isfunction(function)
                and function.__module__ == module.__name__
                and not inspect.signature(function).parameters
            ):
                function()
                direct_count += 1
    print(f"Direct top-level tests: {direct_count} passed")
    print(f"Total tests: {result.testsRun + direct_count} passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
