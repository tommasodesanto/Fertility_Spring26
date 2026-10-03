# Context refresh — October 2–3, 2026

This archive preserves the complete status and memory replaced during the
author-authorized context refresh. It is historical evidence, read on demand.
The live files are [CALIBRATION_STATUS.md](../../CALIBRATION_STATUS.md) and
[memory/AGENT_MEMORY.md](../../memory/AGENT_MEMORY.md).

## What was preserved

- [CALIBRATION_STATUS.before.md](CALIBRATION_STATUS.before.md): exact first
  working-copy capture, including uncommitted additions, at October 2 23:49:27
  New York. It has 17,973 lines and approximately 148,000 words.
- [AGENT_MEMORY.before.md](AGENT_MEMORY.before.md): complete previous durable
  memory, 3,048 lines and 26,582 words.
- [daily_2026-10-02.before.md](daily_2026-10-02.before.md): the startup daily note.
- [manifest.json](manifest.json): original paths, resolved backing paths, byte
  counts, capture time and SHA-256 checksums. These originals are immutable.
- [cutover_manifest.json](cutover_manifest.json) records the final pre-cutover
  hashes and deployed replacements. [The late-update patch](CALIBRATION_STATUS.late_updates.patch)
  reconstructs the final previous status exactly from the immutable first
  snapshot; that reconstruction was hash-verified. Concurrent chats continued
  working during review. Any changed external memory/daily file is separately
  captured before replacement.
- [Refreshed memory](AGENT_MEMORY.refreshed.md) and
  [the short October 3 daily note](daily_2026-10-03.refreshed.md) have tracked copies here:
  active `memory/` is a symlink to the external NightlyMemory backing directory,
  rather than Git-tracked files. The symlink is preserved.
- The installed nightly-script original is backed up before deployment. Only
  the reviewed prompt rules change; its schedule and execution logic stay intact.

No historical model outputs, failed runs, paper baselines, transcript collections
or manuscript files were moved or deleted. `SESSION_DIARY.md` is unchanged.

## Reconciliation decisions

The replacement status follows the current soft selected-point receipts and the
current author instructions. It distinguishes that working calibration from the
frozen September 14 paper baseline, the September 28 stationary export, the
normalized-v1 transition point, and hard/quarter historical comparisons.

Several older labels were superseded by later evidence:

| Earlier claim | Reconciled state and source |
|---|---|
| Soft timing production pending; eight chains planned | Passed v2 native smokes, array 19086987 submitted, then author-expanded array 19087556: 48 total chains. Structured deployment/submission receipts establish this; old preparation prose remains historical. |
| No grid extension or refinement tests completed | Six fixed-price checks completed in `asset_grid_diagnosis_v1/README.md`. Full grid/target/GE convergence remains unverified. |
| Negative entrant cells block zero renter credit | The current provisional nonnegative mean-preserving mapping has no negative wealth cells. Its empirical and estate-recipient limitations remain open. |
| Split age-16/20 entry queue unimplemented | Current `adult_entry.py` and selected-point closure implement the queue; literal person/genealogy equivalence remains an approximation. |
| First-birth target 0.600 or 0.720 remains current | Working contract uses author-selected 1.465. The reviewed A2h builder uses calendar-year event windows; the model observer still approximates the empirical estimator. |
| Child-directed bequest calculation pending / no builder retained | The SCF calculation is complete and a builder exists. Recipient, spouse and creditor mapping and builder-to-historical-receipt identity remain separate limits. |
| Permanent hard-64 policy outcome unknown | Recovered evidence records failure. All permanent dated paths failed; accepted terminal steady states are not accepted transition paths. |
| Old main-text proposals belong in suggestions directory | Current `AGENTS.md` governs `latex/JMP_DS_mock/`; author manuscript restrictions remain intact. |
| Early four-shock draft disabled, test-only | Later author authorization permitted fitting successive surprises; historical attempts failed before producing an accepted fit. Preserve the requested information timing rather than reverting to an announced path. |

Economic changes were not adopted by this refresh. Existing targets, weights,
parameters, optimizer gates, cluster budgets and manuscript wording were not
changed. The deployment state is reported as of retrieved receipts, not as a new
live scheduler check.

## Evidence and retention review

[recent_chat_source_index.json](recent_chat_source_index.json) records 18 retrieved
chat excerpt files, exact chat IDs/titles, excerpt dates and hashes. Retrieval
covered the recent interest-timing, explorer/grid, fertility scaling, down-payment,
calibration normalization, policy mechanism, entrant mapping, PSID income-control,
refactor, Corina slides, transition and utility discussions. Older excerpts were
used where a recent chat referenced an earlier author decision. The Pro discussion
was treated as advice motivating an experiment, not as author adoption or an
independently verified literature review.

[Target provenance review](target_provenance_review.json) records all 14 current
rows, authoritative empirical builders/records, definitions, sample and estimator,
uncertainty where available, weights, model counterparts and caveats. Historical
source records/activation flags are retained explicitly as provenance; current
weights and roles come from the selected soft CSV. September 27 builder/calendar
corrections and the accepted national ACS source are carried forward. This file
does not mint a new calibration contract or claim exact model-data equivalence.

[Memory retention review](memory_retention_review.md) maps the 85 headings of the
previous memory to retained/merged durable rules, live-status material or archived
chronology. The lead corrected the draft review to retain the explicit manuscript
drafting preference, the adjacent-price numerical gotcha and the historical
successive-surprise contract in its proper location. Every original paragraph
remains available in the immutable snapshot.

## Validation and future upkeep

[verification.json](verification.json) checks original hashes, complete fit/parameter counts,
loss arithmetic, candidate/reference links, mirrored `AGENTS.md`/`CLAUDE.md`,
nightly-script syntax, installed-script identity and cutover preservation.

Current state belongs in the consolidated live status, durable rules/preferences
in active memory, and chronological detail in daily notes or experiment READMEs.
The revised agent and nightly-memory instructions require replacing superseded
state, identifying sources and verification times, and preserving unresolved
decisions. Advisory word budgets are a signal to review structure, not permission
to truncate evidence or silently adopt experimental economics.
