#!/usr/bin/env python3
"""Torch-only review-pending lock; never marks a launch approved."""
import hashlib,json
from pathlib import Path
BASE=Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
WORK=BASE/'nightpair_20260925_v1'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()
runtime_files=['run_pair.py','run_task.sh','submit.sh','render_pair.py','prepare_bank.py',
    'build_inputs.py','build_lock.py','tools/e5f_earnings_wealth_contract.py',
    'tools/run_e5f_preference_share_candidate.py']
tax_driver=BASE/'paygo_tax_comparison_20260924/run_paygo_two_rate.py'
lock=dict(status='preflight_pending_review',schema='paired_night_launch_lock_v1',
    run_budget_seconds=21600,search_reserve_seconds=4500,export_reserve_seconds=900,
    per_objective_cap_seconds=3100,workers_per_arm=20,points_per_worker_max=18,
    rates={'greaney_179':.179,'oasi_087510':.08751017424959717},
    first_birth_target_exact=1.465,
    objective_sha256=sha(WORK/'inputs/objective.json'),
    proposed_contract_sha256=sha(WORK/'inputs/proposed_common_contract.json'),
    proposal_bank_sha256=sha(WORK/'inputs/proposal_bank.json'),
    source_manifest_sha256=sha(WORK/'inputs/source_manifest.json'),
    ancestor_sha256=sha(WORK/'ancestor_commute.py'),
    tax_driver_sha256=sha(tax_driver),
    runtime_file_sha256={rel:sha(WORK/rel) for rel in runtime_files},
    repeat_model_absolute_tolerance=0.0,repeat_loss_absolute_tolerance=0.0,
    submission='blocked_until_lead_review_and_explicit_launch_go',
    note='No model solve has run under this new objective; the old checkpoint is a common seed only.')
path=WORK/'inputs/launch_lock.json'
if path.exists():raise RuntimeError('refuse overwrite launch lock')
path.write_text(json.dumps(lock,sort_keys=True,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(path=str(path),sha256=sha(path),status=lock['status']),sort_keys=True))
