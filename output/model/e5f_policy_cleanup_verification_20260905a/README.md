# Focused calendar correctness and bloat cleanup

The author requested obvious cleanup, with correctness first. This packet
verifies the policy-reuse repair and removal of repeated boundary evaluations.
The separate code review and its PDF now include a repair addendum; the package
README has a short reading map and clearly labels the historical June runner.
No economic primitive, target, weight, bound, grid, saving optimizer or
acceptance tolerance changes. No calibration or policy result is promoted.

## What passed

- Nine focused regression checks in `code/model/tools/test_calendar_policy_reuse.py`
  and eighteen existing checks in `code/model/tools/test_e5f_transition_accounting.py`.
  The local project interpreter used `NUMBA_DISABLE_JIT=1`; Torch uses the
  Anaconda 2025.06 runtime without disabling compilation.
- Five full compiled Bellman replays in `replay_final/summary.json`: production
  2023, the morning review candidate, and the previously saved baseline, supply
  and dependent-child credit diagnostic impacts. All saved policy arrays,
  continuation probabilities, distributions, births and market residuals agree
  exactly. Reuse after scratch-state corruption and pickle reload also agrees
  exactly. One fresh Bellman solve runs per state; this is not a policy search.
- Each state has the unchanged seventeen-graph diagnostic packet. All 85 graphs
  match the earlier repaired-source pass exactly. The candidate and three
  diagnostic packets also match their previously reviewed PNGs byte for byte;
  the seventeen production-state plots were visually inspected separately.
- The full five-date historical replay passes in `history_receipt.json`: every
  field of all twelve target rows, every parameter/bound row, and 253 saved
  numeric historical entries match the selected morning reference exactly.
  `history/target_fit_long.csv` and `history/parameter_table.csv` are the complete
  tables. Its unchanged terminal diagnostic packet is the exactly reproduced
  `replay_final/review_candidate/standard_diagnostics/`.

The final scientific bundle is
`33167d84113e2bd38d9ee48dcd9ab0403790348610d998d4032fb8c1797ad3e3`.
`manifest_final.json` records all 51 source hashes, helpers and five input hashes;
its SHA-256 is
`b183d9b0392c66fd1bcca903ab0ee962d117247419495116f9fc7f92b33590d1`.
The scientific changes are the calendar policy object, sequential fertility
argument, and calibration readers. The household solver is byte-unchanged.
The current default-off young-ownership profile predates this patch and explains
the other difference from the frozen fifty-file production bundle.

## Bounded execution and retained failures

| Torch job | Purpose | Result |
|---|---|---|
| 17004380 | Regression checks and initial replay preflight | Checks passed; checkpoint location wrong; zero model solves |
| 17004396 | Five compiled saved-state replays before reviewer edge-case guard | All five passed |
| 17004440 | Final nine tests, eighteen accounting checks and five compiled replays | All passed |
| 17004441 | Full-history reference preflight | Reference project name had an extra `a`; zero model solves |
| 17004455 | Same one-case history, corrected reference location | Complete exact reproduction passed |

Each numerical job requests one CPU, 8 GB and a hard twenty-minute limit under
`torch_pr_570_general` / `cpu_short`. Five fixed-state solves were expected to
take a few minutes including compilation and graphs. The historical check
budgeted one candidate: about five old-equilibrium normalization evaluations
plus fifty dated Bellman solves, roughly ten minutes. No parameter perturbation,
calibration extension, numerical retry, or additional candidate is authorized
by these replay scripts. Every saved-state case writes a completed-case receipt;
the historical driver writes old-normalization and dated-period progress.

The independent reviewer found an out-of-range initial-price edge case after
the first successful replay. The final patch rejects that input before solving;
its ninth regression proves no evaluation occurs. `pre_price_guard.patch`
reconstructs the two earlier source files from the final code, and
`manifest_v2.json` identifies that earlier pass. Apply the reverse patch only
in an isolated copy, using `git apply --unidiff-zero pre_price_guard.patch`.
The final replay supersedes it. `manifest.json`, the first submission, and the
initial history contract retain the two preflight location failures rather than
presenting them as successful checks. The corrected references have unchanged
content hashes.

## Reproduce

Use the final repaired source, with imports pinned by `PYTHONPATH`, and the
original hash-pinned checkpoints. The isolated Torch snapshot is:

`/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/`

Local snapshot: `tmp/e5f_policy_cleanup_20260905a/`. It was built from the tracked
`b38e932` model tree plus only this patch's explicitly owned files; unrelated
working-tree changes were excluded. Shared empirical data is a read-only input
through the snapshot's `code/data` symlink. Original production snapshots and
receipts remain untouched. `submit_final.sh` and `submit_history.sh` retain the
exact runtime and commands; use a fresh destination for reruns, since drivers
refuse populated outputs. `history_contract.json` records the full historical
command, reference hashes, count and stop criteria. The ordinary production
adapter remains frozen to its original source fingerprint and was not relaxed.

The replay never silently fills an old incomplete policy. It solves a fresh
policy at the saved price and first requires every policy array to reproduce;
only then is the new policy reused or serialized. It checks the full forward
operator with original inherited inputs where retained, or the saved prepared
input otherwise, and requires exact saved distributions in either case.

To refresh the complete code-review PDF without solving the model:

```sh
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_full_code_correctness_efficiency_review_20260905.md --output output/pdf/e5f_full_code_correctness_efficiency_review.pdf --date '5 September 2026' --heading 'FERTILITY / CODE CORRECTNESS AND EFFICIENCY'
```

## Limits and next cleanup

The boundary failure probe falls from nineteen evaluations at two prices to two
evaluations; this is not an overall runtime benchmark. Each retained policy now
owns about 23.5 MB of continuation probabilities at the production dimensions.
Parameter consistency beyond price remains a caller responsibility; the
independent reviewer traced the maintained preference, policy and fiscal
invalidation paths. The misleading verbose fertility formula, final duplicate
bisection evaluation, broader module split and memory optimization remain open.
The economic and identification issues in the full review are unaffected.

Git retains the source patch, complete historical target/parameter ledgers and
compact replay receipts. The full graphs and repeated diagnostic tables stay in
the local output folder and on Torch; their recorded hashes and replay command
make them reproducible without duplicating them in the source backup.
