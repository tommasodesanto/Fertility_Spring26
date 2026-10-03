"""Run one pinned entry point with the reviewed, read-only local source overlay."""
import hashlib
import runpy
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
overlay = ROOT / 'output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/bootstrap.py'
source = overlay.read_text()
prefix = source.split("if '--preflight-context' in sys.argv:", 1)[0]
assert 'DIGESTS=' in prefix and 'Frozen overlay write forbidden' in prefix
scope = {'__file__': str(overlay), '__name__': 'authenticated_local_overlay'}
exec(compile(prefix, str(overlay), 'exec'), scope)
entry = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(entry.parent))
sys.argv = [str(entry)] + sys.argv[2:]
if '--deadline-epoch' in sys.argv:
    index = sys.argv.index('--deadline-epoch') + 1
    if sys.argv[index] == 'START_PLUS_1200':
        sys.argv[index] = str(time.time() + 1200)
runpy.run_path(str(entry), run_name='__main__')
