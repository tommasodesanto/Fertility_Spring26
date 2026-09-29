"""Install the authenticated solver transformation in all loaded solver aliases."""
from __future__ import annotations

import difflib
import hashlib
import importlib
import json
import sys
from pathlib import Path

import patch_solver


def install(out, evaluator, root, graft):
    path = patch_solver.SOLVER
    original = (root/path).read_text()
    changed = patch_solver.patch_text(original)
    generated = out/'effective_sources'/path
    generated.parent.mkdir(parents=True, exist_ok=True)
    generated.write_text(changed)
    generated.with_suffix('.diff').write_text(''.join(difflib.unified_diff(
        original.splitlines(True), changed.splitlines(True),
        fromfile=path, tofile=str(generated))))
    # Evaluation may hold aliases imported before this overlay. Graft into each
    # matching live module so those aliases retain their original identity.
    modules = []
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if not filename or Path(filename).name != 'solver.py':
            continue
        candidate = Path(filename).resolve()
        if candidate == (root/path).resolve():
            modules.append(module)
        elif candidate.is_file() and candidate.read_text() == original:
            # The native runtime can keep a byte-identical frozen private copy.
            modules.append(module)
    if not modules:
        model = evaluator.rt['model']
        modules = [importlib.import_module(model.__package__ + '.solver')]
    changed_names = []
    for module in modules:
        changed_names.extend(graft(module, original, changed, generated))
    if not any(name == 'solve_bellman_full_markov_income' for name in changed_names):
        raise RuntimeError('Sequential Bellman function was not grafted')
    digest = lambda value: hashlib.sha256(value.encode()).hexdigest()
    manifest = dict(path=path, original_source=str(root/path),
                    original_sha256=digest(original), effective_sha256=digest(changed),
                    effective_source=str(generated),
                    loaded_module_names=[module.__name__ for module in modules],
                    grafted_names=changed_names)
    (out/'effective_source_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return manifest
