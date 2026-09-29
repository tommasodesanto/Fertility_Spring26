"""Torch-only pin preparation. Does not import a model or grant main-run approval."""
import hashlib
import json
import os
import sys
from pathlib import Path

assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[3]


def read(p): return json.loads(Path(p).read_text())
def pin(p): return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    destination = BASE/'contract.json'
    if destination.exists():
        raise FileExistsError('Versioned contract already exists; never silently re-pin')
    parent = BASE.parent/'contract_v1/contract.json'
    reference = BASE.parent/'resume_v1/selected_export/primary'
    audit = BASE.parent/'measurement_audit_v1'
    c = read(parent)
    plan = read(audit/'bounded_experiment_plan.json')
    assert pin(parent)['sha256'] == '68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
    result = dict(schema='e5f_numerical_pair_v1', status='prepared_not_approved_for_main',
        reference_label=plan['reference_label'], original_contract=pin(parent),
        plan=pin(audit/'bounded_experiment_plan.json'),
        predictions=pin(audit/'proposed_step_predictions.csv'),
        reference_receipt=pin(reference/'receipt.json'),
        reference_fit=pin(reference/'target_fit.csv'),
        reference_parameters=pin(reference/'parameters.csv'),
        files={name: pin(BASE/name) for name in ('runner.py','adapter.py','prepare.py','run.sh')},
        original_source_manifest=c['source_manifest'],
        target_weight_fingerprint=c['lanes']['primary']['target_weight_fingerprint'],
        objective=c['lanes']['primary']['objective'],
        economic_changes=plan['economic_changes'], specification_changes=[],
        numerical_change='Only the initial psi guess differs across the same ten-coordinate trial point.',
        scored_screen_definition='max(sqrt(actual_weight)*abs(arm0_moment-arm1_moment)) <= 0.01',
        promotion='No automatic promotion; separate final exact repeats required',
        main_model_solves_launched=0)
    destination.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps(dict(contract=pin(destination), model_solves=0)))


if __name__ == '__main__': main()
