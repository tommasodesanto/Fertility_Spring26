"""Activate the existing reviewed read-only frozen-source overlay before imports."""
from pathlib import Path
BOOTSTRAP=Path(__file__).resolve().parents[3]/'utility_floor_round2_v1/local_run/bootstrap.py'
source=BOOTSTRAP.read_text()
tail="runpy.run_path(str(HERE/'runner_local.py'),run_name='__main__')"
assert source.count(tail)==1
exec(compile(source.replace(tail,'pass'),str(BOOTSTRAP),'exec'),
     {'__file__':str(BOOTSTRAP),'__name__':'_mechanism_readonly_overlay'})
