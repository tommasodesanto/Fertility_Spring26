"""Authenticate the local frozen context and build an evaluator without solving."""
from __future__ import annotations
import argparse
import hashlib
import json
import runpy
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = PACKET.parents[3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chain', type=int, choices=range(58), required=True)
    args = ap.parse_args()
    controller = HERE/'run_local_psi.py' if args.chain >= 48 else PACKET/'run_psi.py'
    module = runpy.run_path(str(controller), run_name='_purchase_rules_module')
    native, inputs = module['native'], module['inputs']
    native.verify_sources()
    for rel, digest in json.loads((PACKET/'source_pins.json').read_text()).items():
        assert inputs.sha(ROOT/rel) == digest, rel
    for rel, digest in json.loads((PACKET/'engine_pins.json').read_text()).items():
        assert inputs.sha(PACKET/rel) == digest, rel
    own = PACKET/'manifest.json'
    if own.is_file():
        for rel, digest in json.loads(own.read_text())['sha256'].items():
            assert inputs.sha(ROOT/rel) == digest, rel
    else:
        raise FileNotFoundError('Overnight source manifest is not staged yet')
    config = module['CONFIG']
    lane = f"floor_s{module['SEEDS'][args.chain]}"
    design = config['nearby_starts'][args.chain]
    seed, bounds, _ = inputs.seed_and_bounds(lane)
    seed = dict(design['parameters'])
    bounds = dict(bounds)
    bounds['psi_child'] = tuple(config['psi_bounds'])
    bounds['h_P'] = (.1, 2.6)
    coordinates = tuple(inputs.parameters(lane)) + ('psi_child',)
    inputs.LANES[lane].update(seed=seed, bounds=bounds, free_coordinates=list(coordinates))
    P, grid = inputs.proposal(lane)
    P, _ = inputs.entry(P, grid, 'nonnegative_mean')
    import numpy as np
    assert np.asarray(P.phi).shape == (4,)
    P.phi = np.full_like(np.asarray(P.phi, dtype=float), .8)
    arm = ('hard' if args.chain < 24 else 'quarter') if args.chain < 48 else ('hard' if args.chain < 53 else 'quarter')
    P.experimental_purchase_saving_fraction = .25 if arm == 'quarter' else 1.0
    with tempfile.TemporaryDirectory(prefix='auth_', dir=HERE) as temp:
        out = Path(temp)
        Q = native.utility_checks(P, grid, lane, out)
        module['normalized_objective'].make_evaluator(
            out, lane, Q, grid, __import__('time').time() + 900.,
            float(design['initial_price']), native_runner=native,
            exploratory=False,
        )
        assert (out/'native_initializer_verification.json').is_file()
        initializer = json.loads((out/'native_initializer_verification.json').read_text())
        assert initializer['status'] == 'passed_zero_solve'
    receipt = dict(status='authenticated_zero_solve_context', arm=arm, chain=args.chain,
                   lifecycle_solves=0, engine_files=len(json.loads((PACKET/'engine_pins.json').read_text())),
                   multistart_files=len(json.loads((PACKET/'source_pins.json').read_text())),
                   packet_manifest_sha256=hashlib.sha256(own.read_bytes()).hexdigest(),
                   frozen_overlay_files=2)
    path = HERE/f'preflight_chain{args.chain}.json'
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
