"""Build the dated-policy source overlay after the mechanism source is frozen."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = HERE.parents[4]
TRANSITION = ROOT / 'code/model/experiments/transition_readiness'
FIT_PLAN = ROOT / 'output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json'
PACKET_REL = str(PACKET.relative_to(ROOT))


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def files():
    chosen = set(PACKET.glob('*.py')) | set(PACKET.glob('*.json'))
    chosen |= set((PACKET / 'center').glob('*'))
    chosen |= {p for p in (PACKET / 'engines').rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    chosen |= {p for p in (PACKET / 'mechanism').glob('*') if p.suffix in ('.py', '.md')}
    chosen |= {p for p in HERE.glob('*') if p.suffix in ('.py', '.sh', '.json', '.md')
               and p.name not in ('inventory.json', 'stage_receipt.json', 'submission_receipt.json', 'stage_verification.json')}
    chosen |= {p for p in TRANSITION.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    chosen.add(FIT_PLAN)
    chosen |= {ROOT / rel for rel in json.loads((PACKET / 'manifest.json').read_text())['sha256']}
    return {str(p.relative_to(ROOT)): p.read_bytes() for p in sorted(chosen) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=HERE / 'mechanism_stage.tar.gz')
    args = parser.parse_args()
    existing = files()
    pins = json.loads((PACKET / 'manifest.json').read_text())['sha256']
    for rel, expected in pins.items():
        if rel not in existing or sha(existing[rel]) != expected:
            raise SystemExit('Calibration source pin drift: ' + rel)
    for rel, expected in json.loads((PACKET / 'engine_pins.json').read_text()).items():
        key = PACKET_REL + '/' + rel
        if key not in existing or sha(existing[key]) != expected:
            raise SystemExit('Isolated engine pin drift: ' + key)
    for rel in (PACKET_REL + '/results/.mount_point', PACKET_REL + '/collection/readout/.mount_point'):
        existing[rel] = b'mounted at policy runtime\n'
    inventory = {'schema': 'purchase_mechanism_stage_v1', 'files': {rel: sha(blob) for rel, blob in sorted(existing.items())},
                 'calibration_manifest_sha256': sha((PACKET / 'manifest.json').read_bytes()),
                 'mechanism_run_case_sha256': sha((PACKET / 'mechanism/run_case.py').read_bytes()),
                 'selection_mounted_at_runtime': True, 'checkpoint_payloads_included': False}
    (HERE / 'inventory.json').write_text(json.dumps(inventory, indent=2, sort_keys=True) + '\n')
    existing['inventory.json'] = (HERE / 'inventory.json').read_bytes()
    existing['launch_torch.sh'] = (HERE / 'launch_torch.sh').read_bytes()
    existing['submit_torch.sh'] = (HERE / 'submit_torch.sh').read_bytes()
    with args.archive.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w') as tar:
                for rel, blob in sorted(existing.items()):
                    item = tarfile.TarInfo('source/' + rel if rel not in ('inventory.json', 'launch_torch.sh', 'submit_torch.sh') else rel)
                    item.size = len(blob)
                    item.mtime = 0
                    item.mode = 0o755 if rel.endswith('.sh') else 0o644
                    tar.addfile(item, io.BytesIO(blob))
    receipt = {'archive': str(args.archive), 'sha256': sha(args.archive.read_bytes()),
               'bytes': args.archive.stat().st_size, 'source_files': len(inventory['files']),
               'mechanism_run_case_sha256': inventory['mechanism_run_case_sha256']}
    (HERE / 'stage_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
