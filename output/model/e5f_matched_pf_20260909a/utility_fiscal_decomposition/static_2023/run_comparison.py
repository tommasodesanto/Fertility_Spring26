"""Two repeated stationary diagnostics at the retained 2023 preference.

Uses the frozen initial driver without editing the model. This is a stationary
counterpart, not a reproduction of the inherited 2023 transition distribution
or a new calibration. Full old target values, weights and gaps are retained.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

FINGERPRINT = '3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['old_balanced', 'new_balanced'], required=True)
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    manifest = json.loads((folder / 'manifest.json').read_text())
    for name, pin in manifest['bundle_sha256'].items():
        if sha(folder / name) != pin:
            raise ValueError('Bundle hash mismatch: ' + name)
    contract = json.loads((folder / 'contract.json').read_text())
    original = Path(contract['selected_2023_summary'])
    if sha(original) != contract['selected_2023_summary_sha256']:
        raise ValueError('Original 2023 summary changed')
    selected = json.loads(original.read_text())
    if selected['best_candidate']['new_psi_child'] != contract['initial_psi']:
        raise ValueError('Preference differs from the retained 2023 value')
    from intergen_eqscale_seq_optimized.e5_young_ownership_profile import e5_target_system_for_profile
    targets = e5_target_system_for_profile('baseline')
    if targets.fingerprint != FINGERPRINT or selected['target_fingerprint'] != FINGERPRINT:
        raise ValueError('Old target/weight contract changed')
    targets.require_identified(11 if args.case == 'old_balanced' else 10)
    source = Path(contract['source_root'])
    for name, pin in contract['source_sha256'].items():
        if sha(source / name) != pin:
            raise ValueError('Model source changed: ' + name)
    out = folder / 'results' / args.case
    out.parent.mkdir(exist_ok=True)
    if out.exists():
        raise ValueError('Never overwrite an existing case')
    started = time.monotonic()
    subprocess.run([sys.executable, str(source / 'code/model/tools/run_e5f_initial_revision_probe.py'),
                    '--contract', str(folder / 'contract.json'), '--contract-sha256', sha(folder / 'contract.json'),
                    '--case', args.case, '--output', str(out)], check=True, timeout=780)
    summary = json.loads((out / 'summary.json').read_text())
    if summary['repetitions'] != 2 or summary['normalized'] or summary['stationary_solves'] != 2:
        raise ValueError('Expected two fresh unnormalized stationary repetitions')
    fits = []
    losses = []
    for rep in [1, 2]:
        result = json.loads((out / f'repetition_{rep:02d}/summary.json').read_text())
        moments = result['legacy_stationary_moments']
        loss = targets.loss(moments)
        losses.append(loss)
        for name, target, weight in zip(targets.moment_names, targets.target_values, targets.weights):
            value = float(moments[name])
            fits.append(dict(repetition=rep, moment=name, target=target, model=value,
                             gap=value-target, weight=weight, loss_contribution=weight*(value-target)**2,
                             target_fingerprint=FINGERPRINT, interpretation='fixed-parameter stationary diagnostic'))
    if losses[0] != losses[1]:
        raise ValueError('Repeated complete objective differs')
    write_csv(out / 'target_fit.csv', fits)
    receipt = dict(case=args.case, status='repeated_static_2023_preference_diagnostic',
                   elapsed_seconds=time.monotonic()-started, loss=losses[0], exact_repetitions=2,
                   preference_value=contract['initial_psi'], target_fingerprint=FINGERPRINT,
                   calibrated=False, historical_transition_reproduced=False,
                   scope='Stationary counterpart at retained 2023 preference; same inherited structural parameters; balanced pensions',
                   target_fit_sha256=sha(out/'target_fit.csv'))
    (out / 'comparison_receipt.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
