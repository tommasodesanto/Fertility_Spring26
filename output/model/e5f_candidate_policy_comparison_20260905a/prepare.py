#!/usr/bin/env python3
"""Create a policy input view of the verified fixed-parameter history replay."""
import copy
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
REMOTE = Path('/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a')
FROZEN = ROOT / 'tmp/e5f_policy_cleanup_20260905a'
REPLAY = ROOT / 'output/model/e5f_policy_cleanup_verification_20260905a/history'
REFERENCE = ROOT / 'output/model/e5f_morning_refinement_20260905a/round2/task_009'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


if __name__ == '__main__':
    if (OUT / 'plan.json').exists():
        raise FileExistsError('Refusing to overwrite an existing run plan')
    receipt = read(REPLAY.parent / 'history_receipt.json')
    summary = read(REPLAY / 'summary.json')
    reference = read(REFERENCE / 'summary.json')
    assert receipt['status'] == 'pass'
    assert sha(REPLAY / 'summary.json') == receipt['summary_sha256']
    assert sha(REFERENCE / 'summary.json') == receipt['reference_sha256']
    assert receipt['exact_numeric_history_entries'] == 253
    for name in ['target_fit_long.csv', 'parameter_table.csv']:
        assert rows(REPLAY / name) == rows(REFERENCE / name), name
        shutil.copyfile(REPLAY / name, OUT / name)
    for key in ['theta', 'new_psi_child', 'transition_loss']:
        assert summary['best_candidate'][key] == reference['best_candidate'][key], key
    report = copy.deepcopy(summary)
    report['schema'] = 'verified_fixed_parameter_replay_selection_v1'
    report['status'] = 'complete'
    report['scope'] = 'Policy input view of an already verified replay; not a new calibration collector or production promotion'
    report['best_candidate'].update(old_psi_child=summary['old_psi_child'],
        old_completed_fertility=summary['old_model_completed_fertility'])
    report['selection_provenance'] = dict(reference_summary_sha256=sha(REFERENCE / 'summary.json'),
        replay_summary_sha256=sha(REPLAY / 'summary.json'), verification_receipt_sha256=sha(REPLAY.parent / 'history_receipt.json'),
        original_round2_report_sha256=sha(REFERENCE.parent / 'report/summary.json'), production_promoted=False)
    selected = OUT / 'selected_replay_report.json'
    write(selected, report)
    remote_out = REMOTE / OUT.relative_to(ROOT)
    remote_replay = REMOTE / REPLAY.relative_to(ROOT)
    case = remote_replay / 'cases' / summary['best_candidate']['candidate']
    transition = REPLAY / 'cases' / summary['best_candidate']['candidate'] / 'transition_path.csv'
    renewal_sha = hashlib.sha256(json.dumps(summary['renewal_accounting_contract'], sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    arguments = {
        '--selected-report': str(remote_out / selected.name), '--selected-task-summary': str(remote_replay / 'summary.json'),
        '--selected-case-dir': str(case), '--selected-case-transition': str(case / 'transition_path.csv'),
        '--source': summary['source'], '--expected-report-sha256': sha(selected),
        '--expected-task-summary-sha256': sha(REPLAY / 'summary.json'),
        '--expected-case-transition-sha256': sha(transition), '--expected-source-sha256': summary['source_sha256'],
        '--expected-target-fingerprint': summary['target_fingerprint'],
        '--expected-code-bundle-sha256': summary['code_fingerprints']['bundle_sha256'],
        '--expected-renewal-contract-sha256': renewal_sha,
        '--expected-scientific-contract-sha256': 'collector-style-not-applicable',
        '--expected-selection-sha256': 'collector-style-not-applicable',
    }
    input_hashes = {str(remote_out / selected.name):sha(selected), str(remote_replay / 'summary.json'):sha(REPLAY / 'summary.json'),
        str(case / 'transition_path.csv'):sha(transition), summary['source']:summary['source_sha256']}
    source_hashes = {}
    names = ['run_e5f_post2023_policy_mechanisms.py','run_e5f_post2023_no_policy_continuations.py',
        'collect_e5f_post2023_policy_mechanisms.py','run_e5f_post2023_rebated_property_tax_smoke.py',
        'run_e5f_rebated_tax_shapley_diagnosis.py','run_e5f_independent_numerical_audit.py']
    for name in names:
        relative = 'code/model/tools/' + name
        source_hashes[relative] = sha(FROZEN / relative)
    for relative in ['code/model/tools/run_e5f_candidate_policy_comparison.py','code/cluster/run_e5f_candidate_policy_comparison.sh']:
        source_hashes[relative] = sha(ROOT / relative)
    plan = dict(schema='e5f_candidate_policy_comparison_v1', output_root=str(remote_out),
        candidate_loss=summary['best_candidate']['transition_loss'], code_bundle_sha256=summary['code_fingerprints']['bundle_sha256'],
        target_fingerprint=summary['target_fingerprint'], selected_report_sha256=sha(selected),
        source_hashes=source_hashes, input_hashes=input_hashes,
        common_arguments=[item for pair in arguments.items() for item in pair],
        launch_deadline_epoch=datetime(2026,9,5,22,50,tzinfo=timezone.utc).timestamp(),
        author_window_end_utc='2026-09-06T00:11:00Z', production_promoted=False,
        scope='same calibrated candidate and inherited finite-horizon policy contract; repaired source exactly verified against original history',
        stages={'smoke':{'cases':5,'dates_per_case':2,'minutes_per_case':40},
                'full':{'cases':5,'dates_per_case':11,'minutes_per_case':55},
                'rebate':{'impact_equilibria':3,'component_cells':8,'minutes':55}})
    write(OUT / 'plan.json', plan)
    print('Prepared immutable plan:', sha(OUT / 'plan.json'))
