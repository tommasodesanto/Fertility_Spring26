# NumPy 2 checkpoint module names mapped to NumPy 1.24 equivalents.
# Array contents and numerical functions are unchanged.
import importlib
import sys
import numpy
if int(numpy.__version__.split('.')[0]) < 2:
    for suffix in ('', '.multiarray', '.numeric', '.umath', '._multiarray_umath'):
        sys.modules['numpy._core' + suffix] = importlib.import_module('numpy.core' + suffix)
# Python 3.13 moved concrete Path classes; retain equivalent local classes.
import pathlib
sys.modules.setdefault('pathlib._local', pathlib)
