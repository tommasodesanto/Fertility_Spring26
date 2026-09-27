#!/usr/bin/env python3
"""Regenerate the unchanged 17 diagnostic figures from an authenticated case.

Imports the pinned runtime and loads saved arrays; performs no model solve.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import pickle
import sys
from pathlib import Path


def main(a):
    def sha(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    assert sha(a.contract) == os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']
    c = json.loads(a.contract.read_text())
    assert sha(c['files']['driver']['path']) == c['files']['driver']['sha256']
    spec = importlib.util.spec_from_file_location('calibration_diagnostic_controller',
                                                c['files']['driver']['path'])
    driver = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = driver; spec.loader.exec_module(driver)
    c, objective = driver.verify(a.contract); driver.verify_execution(c)
    receipt = driver.read(a.case / 'receipt.json')
    assert receipt['target_system_sha256'] == c['objective']['sha256']
    checkpoint = a.case / 'initial_state.pkl.gz'
    assert driver.sha(checkpoint) == receipt['case_checkpoint_sha256']
    a.output.mkdir(parents=True, exist_ok=False)
    (a.output / 'runtime').mkdir()
    _, _, _, runtime, _, _, _, _ = driver.setup(c, objective, receipt['point'], a.output / 'runtime')
    with gzip.open(checkpoint, 'rb') as stream:
        packet = pickle.load(stream)
    runtime['audit'].standard_diagnostics(packet, a.output, validate_production_young=False)
    assert len(list((a.output / 'standard_diagnostics').glob('*.png'))) == 17
    driver.write(a.output / 'diagnostic_receipt.json', {
        'status': 'saved_case_diagnostics_regenerated', 'case': str(a.case.resolve()),
        'checkpoint_sha256': receipt['case_checkpoint_sha256'],
        'contract_sha256': driver.sha(a.contract), 'native_solves': 0})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract', type=Path, required=True)
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args())
