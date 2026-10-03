"""Build an immutable, source-only Torch package for the matched soft-timing search."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKETS = Path('output/model/fixed_reference_economics_20260928')
DEPLOY = ROOT / PACKETS / 'soft_timing_calibration_20261002_v1/deployment'
NORMAL = PACKETS / 'normalized_calibration_v2'
TIMING = PACKETS / 'purchase_timing_sandbox_v1'
SOFT = PACKETS / 'soft_timing_review_v1'
SANDBOX = Path('code/model/experiments/purchase_timing_sandbox')
PREVIOUS = PACKETS / 'purchase_rules_overnight_v1/previous_soft_checkpoints/previous_soft_best_before_stop.json'
REMOTE = '/scratch/td2248/projects/soft_timing_calibration_20261002_v2'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    DEPLOY.mkdir(parents=True, exist_ok=True)
    pins = json.loads((ROOT / NORMAL / 'source_pins.json').read_text())
    for rel, expected in pins.items():
        path = ROOT / rel
        if not path.is_file() or digest(path) != expected:
            raise SystemExit(f'Original source pin drift: {rel}')
    selection = json.loads((ROOT / SOFT / 'soft_selected.json').read_text())
    source = ROOT / selection['source']
    if digest(source) != selection['source_sha256'] or source.resolve() != (ROOT / PREVIOUS).resolve():
        raise SystemExit('Selected soft checkpoint source drift')
    manifest = json.loads((ROOT / TIMING / 'manifest.json').read_text())
    plan = json.loads((ROOT / NORMAL / 'plan.json').read_text())
    for key in ('target_fingerprint', 'weight_fingerprint'):
        if manifest[key] != plan[key]:
            raise SystemExit(f'Mixed {key}')
    if selection['selected']['weight_fingerprint'] != plan['weight_fingerprint']:
        raise SystemExit('Selected soft weight fingerprint drift')
    selected_contract = [{k: row[k] for k in ('moment', 'target', 'weight', 'role')}
                         for row in selection['selected']['target_fit']]
    if selected_contract != plan['base_target_contract']:
        raise SystemExit('Selected soft target contract drift')
    for pair, expected in manifest['source_pairs'].items():
        original, sandbox = pair.split('|')
        if digest(ROOT / original) != expected['original'] or digest(ROOT / sandbox) != expected['sandbox']:
            raise SystemExit(f'Timing source drift: {pair}')
    files = set(map(Path, pins))
    files.update(p.relative_to(ROOT) for p in (ROOT / NORMAL).iterdir() if p.is_file() and p.suffix in ('.py', '.json'))
    files.update(p.relative_to(ROOT) for p in (ROOT / TIMING).iterdir() if p.is_file() and p.suffix in ('.py', '.json'))
    files.update(p.relative_to(ROOT) for p in (ROOT / SANDBOX / 'source').rglob('*')
                 if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'
                 and p.name != '.DS_Store')
    files.update((SANDBOX / 'calibrate.py', SANDBOX / 'evaluate_selected.py', SANDBOX / 'local_entry.py',
                  SOFT / 'soft_selected.json', PREVIOUS))
    files.add(PACKETS / 'soft_timing_calibration_20261002_v1/driver_plan.json')
    for rel in files:
        if not (ROOT / rel).is_file():
            raise SystemExit(f'Missing stage input: {rel}')
    entrypoints = {p.name: digest(p) for p in sorted(HERE.glob('*.sh'))}
    entrypoints.update({name: digest(HERE / name) for name in ('verify_stage.py', 'verify_smoke_gate.py', 'collect_torch.py')})
    inventory = {'files': {str(rel): digest(ROOT / rel) for rel in sorted(files)},
                 'entrypoints': entrypoints,
                 'target_fingerprint': plan['target_fingerprint'],
                 'weight_fingerprint': plan['weight_fingerprint'],
                 'selected_source_sha256': selection['source_sha256'],
                 'remote_root': REMOTE, 'source_prefix': 'source/',
                 'original_source_pin_count': len(pins), 'no_cache_or_results': True}
    (DEPLOY / 'inventory.json').write_text(json.dumps(inventory, indent=2, sort_keys=True) + '\n')
    archive = DEPLOY / 'stage.tar.gz'
    entries = [(Path('source') / rel, ROOT / rel) for rel in sorted(files)]
    entries += [(Path('inventory.json'), DEPLOY / 'inventory.json')]
    entries += [(Path(p.name), p) for p in sorted(HERE.glob('*.sh'))]
    entries += [(Path('verify_stage.py'), HERE / 'verify_stage.py')]
    entries += [(Path('verify_smoke_gate.py'), HERE / 'verify_smoke_gate.py')]
    entries += [(Path('collect_torch.py'), HERE / 'collect_torch.py')]
    with archive.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode='w') as tar:
                for name, path in entries:
                    blob = path.read_bytes()
                    info = tarfile.TarInfo(str(name))
                    info.size = len(blob)
                    info.mtime = 0
                    info.mode = 0o755 if path.suffix == '.sh' else 0o644
                    tar.addfile(info, io.BytesIO(blob))
    receipt = {'archive': str(archive), 'sha256': digest(archive), 'bytes': archive.stat().st_size,
               'source_files': len(files), 'original_source_pins': len(pins),
               'target_fingerprint': plan['target_fingerprint'],
               'weight_fingerprint': plan['weight_fingerprint'], 'selected_source_sha256': selection['source_sha256']}
    (DEPLOY / 'stage_receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
