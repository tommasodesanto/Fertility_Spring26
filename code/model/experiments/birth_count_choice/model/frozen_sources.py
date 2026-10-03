"""Read-only exact-byte recovery for existing missing frozen dependency paths."""
from pathlib import Path
import hashlib
import io
import json
import sys


def install_recovered_sources():
    marker = '_birth_count_exact_frozen_recovery'
    if getattr(sys, marker, False):
        return
    base = Path(__file__).resolve().parents[1] / 'frozen_sources'
    data = json.loads((base / 'mapping.json').read_text())
    mapping = data.get('mapping', data)
    redirects = {}
    for original, entry in mapping.items():
        target = (base / entry['recovered_path']).resolve()
        if not target.is_relative_to(base.resolve()):
            raise RuntimeError('Recovered source escapes experiment directory')
        if hashlib.sha256(target.read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError('Recovered frozen source digest differs: ' + original)
        redirects[Path(original).resolve()] = target
    old_path_open, old_io_open = Path.open, io.open
    def target_for(file, mode):
        if isinstance(file, (str, Path)) and not any(x in mode for x in ('w', 'a', 'x', '+')):
            path = Path(file).resolve()
            if path in redirects and not path.exists():
                return redirects[path]
        return file
    def path_open(self, mode='r', *args, **kwargs):
        return old_path_open(target_for(self, mode), mode, *args, **kwargs)
    def io_open(file, mode='r', *args, **kwargs):
        return old_io_open(target_for(file, mode), mode, *args, **kwargs)
    Path.open, io.open = path_open, io_open
    package_root = base / 'calibration_archive/model_legacy_20261003'
    sys.path.insert(0, str(package_root))
    setattr(sys, marker, True)
    return {'mapping': str(base / 'mapping.json'), 'exact_recovered_files': len(redirects),
            'original_working_tree_preserved': True}
