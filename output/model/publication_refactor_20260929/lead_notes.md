# Lead checks during implementation

- Exact requested Claude model availability was confirmed from CLI metadata:
  `claude-opus-5-5`, first-party Claude Max authentication. No fallback model.
- Manifest checkpoint identity is resolved: repeat0212 (`b15ba92d...`), distinct
  from the seed/original and other repeat identities. Do not diagnose a mismatch
  from their differing hashes alone.
- Historical full closed-credit GE is 886.527 seconds for six lifecycle calls
  including one repeat; it is not the older September26 multi-psi objective.
- Lead inspected `kernels.py:881-947`: exhaustive saving checks all segment
  endpoints and interior FOC candidates; static renter demand is already
  analytical. Preserve candidate ordering and strict-improvement tie breaking.
- Potential later exact optimization: continuation slopes and FOC powers are
  constant across current-wealth rows within a fixed continuation column.
  `full_renter_block_kernel`/owner scan currently recomputes them in scalar calls.
  Hoisting them requires exact equivalence tests and profiling; this is a
  hypothesis, not authority to change optimization completeness or formulas.
- The one-market pension has an existing analytical age/income recursion in
  `code/model/tools/e5f_stationary_paygo.py`, certified against the actual solved
  distribution. Preserve this shortcut rather than introducing a fiscal loop.

## Borrowing implementation reviewed and compiled checks passed

The other chat's source is
`output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/overlay/`.
Its README and `overrides.json` define `unsecured_credit_limit=None` for legacy
reproduction and an explicit finite scalar for the correction. The planned
scalar is zero. This directory is read-only to this refactor task.

The lead reviewed its complete diff in parameters.py, solver.py and kernels.py.
The scalar replaces the old renter floor (it is not combined with the old
floor); current mortality/terminal solvency raises it to zero when required.
Both tenure kernels and their Python fallback apply the raw sale balance >=0
gate before any clipping. Buyer and owner-stayer kernel/call sites are unchanged.
Native natural-credit mode and unsupported legacy/factored Bellman routes reject
the scalar rather than silently applying another contract. This agrees with the
author-selected rule. The author subsequently authorized small one-core local
tests. The other chat cancelled pending Torch job18843355 and recorded compiled
PASS in `verification_local_v1/receipt.json`: 5.33 seconds, 219,807,744 bytes
peak RSS, no lifecycle solve. The discarded first draft is explicitly unaccepted.

## Lead local component checks

Under the new one-core instruction, `local_smoke_v1/` records 11 passing small
tests and 4,000 byte-identical compiled original/indexed saving choices on the
actual 160-point wealth grid. The compiled pass including import/JIT took
1.65 seconds and about181 MB. A fresh engine import loaded no legacy engine
modules. These checks establish component behavior only; no full-solution
equivalence or GE speedup is yet established.

## Input export

Torch job18844678 uses the independently pinned three-file `source_export`
stage, the read-only September28 reference, one core,24 GB and a10-minute cap.
The physical repeat0212 checkpoint is used because `selected_export/primary`
is an absolute Mac-path symlink that does not resolve before container binding.
The checkpoint SHA remains unchanged. This job only exports inputs and oracle
arrays; it does not solve or recalibrate.

Export completed in13 seconds (about1.8 GB peak RSS); see
`export_receipt.json`. Bundle SHA is
`427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7`,
reference price0.7898695017462086. There are68 top-level exported arrays including
`stationary_g_pre`. Sol's completed `sol_array_coverage.md` establishes that
the relevant fixed-price control instead has113 nested array paths: one grid,
20 shared arrays,67 solution arrays,24 evaluation arrays and stationary_g_pre.
The earlier107 count belongs to another packet scope. Preserve the113-path
control comparison, including lab-owned shared precomputation. The14 saved
underscore parameter fields are overwritten before reads on the complete
Bellman-to-KFE path; verify that claim with fresh parameters in the replay.

## Staged verification

- Component job18845050 stopped in pytest collection: a default path expression
  was evaluated eagerly despite REFACTOR_ROOT. No model code ran. Luna fixed
  only that expression and checked the explicit-root and local fallback cases.
- New immutable stage `lab_src_pass3b` is pinned by
  `pass3b_SOURCE_SHA256SUMS` (manifest8541f1b6...). Component job18845564 has a
  10-minute cap and24 GiB. An attempted reduction to8 GiB while pending was
  rejected by Slurm; the original request remains in force.
- Preliminary fixed-price job18845813 depends on that component job succeeding.
  It retains one CPU/24 GiB and a15-minute cap for two fresh lifecycle replays
  plus oracle checks. Historical per-solve time is about107 seconds. It cannot
  yet certify the20 shared arrays, which pass4 is correcting separately.

## Local execution and same-machine validation

The author explicitly permits individual one-core model runs locally, including
overnight. AGENTS.md/CLAUDE.md and memory now record this. Machine is Apple M5
Pro, 48 GB, arm64. Other chats still own upstream borrowing-runtime validation.
The original model source hashes match the frozen 13-module source inventory.

Preliminary Torch replay18845813 completed: 25 component tests pass; two fresh
lifecycle calls took74.35/64.56 seconds; 67 solution arrays, full14/31 tables
and17plot hashes match the reference. Its20 shared arrays came from the old
precompute (not an independent lab check); final local tests save lab shared.
All99 solver definitions were subsequently split without source/bytecode change.

Existing local Python3.10/NumPy1.24 could not read NumPy2 checkpoint objects;
an isolated NumPy2 Python3.10 runtime then failed on Python3.13 pathlib objects.
No model ran in either failed attempt. The separate Python3.13.15 environment
(local_env_v1/venv313) reads the authenticated checkpoint successfully; its input
bundle and reference-array SHA values exactly match the Torch export. Local
export took1.58 seconds, about1.24GB max RSS. No pickle compatibility shim used.

Baseline and candidate will run in that identical local environment, one at a
time. Required limits remain2 fixed-price solves per candidate stage,18 at-price
calls/900 solve seconds per GE,2700 seconds per paired benchmark, sampled12GiB
RSS. Same-machine old/new equality is mandatory; Mac/Torch differences are
separate diagnostics. No recalibration or closure normalization is authorized.

## Independent replay and full local GE — September 30

Torch18850069 passed all four independent fixed-price certificates: scalar and
indexed engines, two fresh repetitions each. Each includes113exact finite array
paths,14target rows,31parameter rows and17unchanged plot hashes. Lab precompute
supplies the20shared arrays. `fp_pair_gridfix_18850069/COLLECTION.md` records the
collected evidence; three representative standard plots were visually inspected.
The previously known high-wealth policy shapes are retained, not repaired.

`native_local_pair_v1/` records a same-machine native fullGE old/indexed pair:
152.650099166 versus99.546702958 solve seconds (34.7877% less), four household
solves and five forward-distribution passes each, fresh per-engine caches,
one core on AppleM5Pro. Strict90-array comparison and effective public-parameter
comparison passed. Externalwall156.000/103.098s; sampledRSS1.865/1.909GiB. No
source pin was bypassed. This timing excludes historical tables/17plot reporting.
Absolute renewal gap1.70281219775e-6 is reported, not claimed below1e-6: the
lead1e-6 threshold measures only its difference from reference7.91996696621e-7.

Torch matched GE/reporting job18851943 uses the immutable indexed_src_gridfix
engine and frozen original,1CPU/12GiB; maximum18at-price calls/900solve seconds
each;2700s phase/3000s scheduler ceilings; no retries. Separate strict comparator
SHA141deab226f35abb1981cf0473edf05709b67b768ed80eb0fec7d841366700e8
requires identical complete native array key sets. At dispatch, final reporting
is pending. Opus final package pass is permitted to promote only the exact
tested kernel and move verification machinery, with no numerical edits.

Other chat's newest read confirms compiled upstream borrowing smoke passed
and queued its sole unset-credit baseline control18849552. It reconfirms the
two zero-credit entrant cells are infeasible. We do not duplicate that work.
