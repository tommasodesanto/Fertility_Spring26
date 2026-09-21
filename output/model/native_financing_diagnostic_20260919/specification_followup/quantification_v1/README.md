# Quantification and specification follow-up — September 20 evening

## Authorization and priority

Tommaso renewed authorization for overnight research, bounded experiments and
Claude Max discussion, with an updated view tomorrow. His priority is to settle
the model choices and calibration before drawing conclusions from policies.
This is a new follow-up after the earlier September 20 diagnostic batch closed.
The earlier decision packet is evidence; it is not a frozen adopted baseline.

Deliverable: one recommended minimal baseline, a table distinguishing economic
choices, externally measured inputs, estimated parameters and numerical
approximations, a reviewed diagnosis of the existing fit, and at most two
well-defined next experiments. Aim for a morning report by September 21 at
09:00 America/New_York, explicitly labeling any blocked or unrun work.

## Active bounded work

1. **Claude Max specification review:** authenticated first-party subscription,
   Fable/max effort; fork of the earlier review, read-only tools, at most 40
   turns and 1,200 seconds. Prompt, launcher, stream and receipt live in
   `claude_review/`. The lead assesses the economic claims; the review does not
   adopt choices. No automatic restart on failure.
2. **Saved calibration trade-offs:** one Luna worker, 15 minutes, owns only
   `code/model/tools/analyze_e5f_saved_income_fit_tradeoffs.py` and `saved_fit/`.
   Reduce the 96 already evaluated proposals to verified complete fit and
   parameter tables, descriptive economic-block losses and an observed
   rooms/ownership trade-off plot. No new solves, policy conclusions, inferred
   Jacobians from adaptive proposals, or claims of an attainable frontier.
3. **Reusable calibration run scope:** one read-only Luna worker, 12 minutes,
   identifies the pinned driver, source, objective, runtime and smallest useful
   numerical design. The lead must write and review a concrete launch contract
   before submission.

## Cluster status and numerical limits

**Update:** the author restored Torch authentication and the normal queue probe
now succeeds. A [controlled calibration panel](sensitivity/DESIGN.md) is in
preparation; inspect a submission receipt before calling it launched.
The following failure is retained as history.

At 23:05 EDT September 20, `code/cluster/torch.sh status` failed SSH authentication
with `Permission denied (gssapi-keyex,gssapi-with-mic,password,keyboard-interactive)`.
The author was asked to refresh the usual login. **No new cluster job has been
submitted.** Local review and saved-output analysis can proceed. Do not work
around credentials or imply that a prepared experiment is running.

If access is restored, numerical work is limited to specification/measurement,
fit or local sensitivity, not a financing or policy sweep. Initial ceiling:
32 full objective evaluations, including smoke controls and selected-point
repetitions; at most eight single-threaded workers. This is an upper envelope,
not an instruction to consume it. A full objective can require several nested
stationary solves; the actual launch manifest must count both and use measured
runtime. No new batch after 02:00 September 21; finish by 07:30 for collection.
No laptop model-search substitution when cluster access is unavailable.

Before any launch: exact source/input/target fingerprints, parameter bounds,
measurement and entry definitions, actual case list or algorithm, solve-count
estimate, per-case and global deadlines, loop smoke and its outputs, progress
at each case or five minutes, latest and best summaries, and stop criteria.
Production must depend on successful smoke without laptop chaining. Stop on
source/target/accounting or unexpected errors; expected rejected candidates
stay explicit. No retry, gate relaxation or expanded scope without a reviewed
changed hypothesis/method. Cancel only never-started blocked dependents.

Retain the existing complete objective for comparable fit diagnostics; any
proposed alternative measurement system needs a separate contract and cannot
be compared by scalar loss as if identical. Near-bound flags must use actual
plan bounds (especially annual beta <= 0.99). Selected candidates need full
13-target/17-parameter tables, source/objective fingerprints, the unchanged
17-plot packet and two exact repetitions. Nothing is adopted automatically.

## Morning report and stop

Write the lead-reviewed recommendation in `morning_view.md` and index its
supporting evidence here. Update canonical status and the parent README after
meaningful completion. Preserve all unrelated dirty work; commit/push only
owned files. The morning follow-up should notify once with the completed
recommendation and any concrete blocker, then delete itself. If the laptop is
asleep, local review and the scheduled follow-up wait until it is available;
only jobs already submitted to Torch continue independently.

The hourly follow-up is `quantification-morning-review`, created for this task.
It should deliver the morning report and then delete itself.

## Reviewed external advice

Claude Max completed its bounded review; see [memo](claude_review/final.md) and
[lead assessment](claude_review/lead_review.md). Additional family ownership
preferences and architectural-impossibility claims were not accepted.

## Submitted execution

Smoke **18153070**, production **18153071**, submitted 23:28 EDT.
[Immutable launch manifest](sensitivity/launch_manifest.json),
[submission](sensitivity/submission.json), and [design](sensitivity/DESIGN.md).
Production depends on successful smoke and requires no laptop chaining.
Three local tests and an independent mocked full loop passed; the native smoke
is the remaining numerical gate. The controlled panel is the sole new model
batch; no additional experiment is preapproved for automatic launch.
Read the [working morning view](morning_view.md) and update it from these jobs.
