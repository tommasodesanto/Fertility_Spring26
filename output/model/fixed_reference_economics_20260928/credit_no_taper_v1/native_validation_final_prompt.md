# Codex worker task — five concrete remaining harness defects

Full project startup, bounded latest context reads only; no giant status dump.
10-minute worker_fast task. Own ONLY native_validation_v3 under this packet.
Preserve v1/v2. Copy four small prepared files from v2 as successor; fix the below.
No other agents, no git, no Slurm submission, no heavy local work or cleanup.

1. The ACTUAL native runtime has require_current_model in
code/model/tools/e5f_current_transition_runtime.py:68–76. It rejects the overlay
parameters path because it is outside the native package. Prepare an EXPLICIT
isolated identity-guard extension in this harness, installed before parent
authenticate, never by altering module.__file__ or changing original source.
It must require the original native solver path/hash, reject ANY other mixed
module path, and allow precisely intergen_eqscale_seq_optimized.parameters
at plan.overlay_parameters_path with exact OVERLAY_SHA. Source/contract pins
remain checked by existing authenticate/runtime.setup/verify_sources. No
numerical/economic/scientific gate may be changed. Keep original guard when
overlay absent; disclose this sole authenticated-module exception in receipts.
Lead approves that identity extension specification; this is not permission
to waive any gate or economic requirement. The actual prepared.rt.model must
use the registered overlay's setup_parameters and unsecured_debt_floor, verify
that after authentication. Test guard rejects another external package module.

2. self-test imports e5f_overnight_estate_audit without tools on sys.path.
Add ROOT/code/model/tools before that helper is imported, pin its source as
already done by runtime/plan. Verify the test paths on Torch only.
3. Synthetic fixture currently has owner stayer mass .4 but zero total owner
mass. Set g owner entries to include stay. This is a fixture consistency bug,
not permission to change policy_mass_branches or real incidence accounting.
4. Synthetic fixture's renter saving is -1 at terminal age, violating the new
zero estate bound. Set terminal renter saving to zero. Add a separate failing
fixture to assert negative terminal saving is rejected.
5. Actual-controller timeout mock injects clock=deadline, which fails the
initial RuntimeError gate rather than TimeoutError expected. Provide controlled
clock progression: initial checks before deadline, then while child still alive
clock reaches case deadline and raises TimeoutError. Test actual termination
path (no sleep > a few seconds). Avoid spin polling; a tiny sleeper is fine.

Review launcher now complete in v2 and preserve entry-time1200 budget, two
360s fresh cases,1CPU24GiB, original160grid, fixedprices/fiscal/psi, no new credit.
For self-test use only imports, controller fixtures and shaped arrays, no
checkpoint reads/lifecycle/KFE/GE. Add option --self-test-only to launcher so
lead can run a 5-minute preflight before production; no lifecycle before tested
preflight receipt hash. Production requires successful preflight receipt pin
in plan AND identical tested source/plan scientific fields. Do not create a
plan circular pin dependency; preflight can be an explicit launcher argument
with SHA independent of immutable scientific plan. Refuse duplicate outputs.
Immutable staging recipe must pin all files and guard extension. No old budget
restart, reference switch, calibration, or positive credit design.

Return concise corrected diff/files/commands, all checks actually run. Do NOT
claim all review items solved just from syntax. Stop on unknown premise.
