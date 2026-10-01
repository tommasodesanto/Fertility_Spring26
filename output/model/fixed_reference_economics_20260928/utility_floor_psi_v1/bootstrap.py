"""Reuse the reviewed read-only local compatibility overlay without changing it."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'utility_floor_round2_v1/local_run/bootstrap.py'
source = OLD.read_text()
tail = "runpy.run_path(str(HERE/'runner_local.py'),run_name='__main__')"
assert source.count(tail) == 1
source = source.replace(tail, "sys.path.insert(0, str(HERE)); runpy.run_path(" + repr(str(HERE/'run_psi.py')) + ", run_name='__main__')")
exec(compile(source, str(OLD), 'exec'), {'__file__': str(OLD), '__name__': '__main__'})
