# Overnight dated financing experiments: launch handoff

**October 2, 08:23 New York progress.** The 48/64-date extension's four 80%
control cases, both temporary 100% financing 48-date cases, and the
quarter-saving temporary 64-date case passed all numerical gates against
matched controls. In the first period, hard-rule 48-date
first-birth flow fell from 0.04980374667823033 to 0.04953396493829134
(0.5416896477%); its aggregate first-birth hazard fell by 0.08487845
percentage points. Quarter-saving flow fell from 0.04974551614189717 to
0.0495313221352156 (0.4305795241%); its hazard fell by 0.06709008 percentage
points at 48 dates. At 64 dates, quarter-saving first-period flow fell from
0.04974551614189717 to 0.04953228420886753 (0.4286455334%); its hazard
fell by 0.06678874 percentage points. Five other policy cases were running,
with no new failures and checkpoints less than 13 minutes old. The hard-rule
64-date check and all permanent cases remain pending. The quarter 48- and
64-date first-period responses are close; a full-state horizon-overlap check
has not been completed. See
[hard T48](accepted_hard_temporary_h48_0722.json),
[quarter T48](accepted_quarter_temporary_h48_0752.json), and
[quarter T64](accepted_quarter_temporary_h64_0822.json), plus
[progress_0822.json](progress_0822.json). The underlying hard and quarter
fits remain fresh-postchecked at losses 97.0112198 and 51.5560360,
respectively, without optimizer convergence certificates. They are
experimental selected points, not adopted paper baselines.
The date-zero buyer diagnostics for the accepted 48-date hard and quarter
paths passed without model solves; see [dated v5](../buyer_diagnostics/dated_v5/README.md)
and [dated v6](../buyer_diagnostics/dated_v6/README.md).

The following is the original, now historical launch handoff. At preparation,
this was a **prepared, unsubmitted** Torch launcher for the diagnostic plan:
two separately fitted 80% purchase-rule baselines (hard and quarter-saving),
each compared with a 100% financed-share change for one date and permanently.
It runs control, temporary and permanent paths at both 12 and 16 dates: 12
dated cases. The 80% baseline parameters, housing-supply coefficient \(H_0\),
earnings, entry, fiscal objects, child preference and all calibration targets
remain those of each selected rule. The mechanism source in `../mechanism/`
owns the economics and numerical roots; this deployment does not modify it.

Each case has one CPU, 32 GiB, one numerical thread, at most four hours from its actual start, at most 1,024 native policy calls, and an absolute stop at **10:00 New York on October 2, 2026** (epoch `1790949600`). With a 05:30 start, the four-hour cap ends near 09:30. A later start has less time. These are limits, not completion guarantees. Recent single fixed-price household solves took roughly 2–4 minutes before warm reuse; dated mappings may use additional one-date solves, so 1,024 calls is a generous ceiling rather than an expected count. The selected fit and native one-date smoke must pass first. Per-case `run.log`, native mapping progress, latest/best receipts and a terminal or failure receipt remain readable during execution. No automatic retry or contract fallback is configured.

The source overlay uses the same frozen calibration runtime, 120×9 inputs and isolated hard/quarter engines as the running calibration. It adds the current `code/model/experiments/transition_readiness/` source and the pinned dated root controls at `output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json`. The mechanism source is SHA-pinned in `inventory.json` when staged. The chosen postchecks and selected-winner JSONs are mounted read-only from a separate mechanism snapshot. Original Torch, Torch restart, broad-region Torch, local and local restart results stay in their physical locations.

The selected winners come only from `../collection/collect.py`, after a **fresh selected-point postcheck** for each arm. The collector rejects mixed target and weight fingerprints and invalid full 14-target/31-parameter/17-plot reports. `prepare_selection.py` checks the collector's physical source and parent provenance, verifies all postcheck file hashes and the complete target/parameter/plot identity, then copies each winner to immutable `purchase_mechanism_v1/selected_postchecks/chain_N/postcheck`. It publishes the selection manifest only after both copies match. A changed winner requires a new named snapshot version; it cannot replace an existing chain snapshot or selection manifest. A search checkpoint alone is never a policy baseline.

Original review and launch sequence (historical):

1. Confirm `../mechanism/run_case.py`, `selected_runtime.py`, `integration.py` and `dated_phi.py` have the final reviewed hashes in the stage inventory. Check that the permanent terminal and dated price/pension roots use the original numerical gates and fixed fitted \(H_0\).
2. Run the mechanism's zero-solve tests, then stage this source-only packet with `bash mechanism_deployment/stage_torch.sh`. Compare the local and remote archive SHA-256. Execute `preflight_torch.sh` to authenticate both mounted source variants without a model solve.
3. After both 80% calibrations have fresh postchecks, run `collection/collect.py --fetch-reports`; inspect its `selected_hard.json`, `selected_quarter.json`, full target-fit and parameter tables, standard plots, and native arrays. Run `python -m unittest mechanism_deployment/test_prepare_selection.py`, then `python mechanism_deployment/prepare_selection.py --apply` from this packet; inspect `selection_snapshot/manifest.json` and the SHA-verified remote `selected_postchecks/` roots.
4. Submit with `ssh torch 'bash /scratch/td2248/projects/purchase_mechanism_v1/submit_torch.sh'` **after lead review**. It submits two independent one-date native baseline smokes, then the 12-case array with an `afterok` dependency. A failed smoke leaves all production cases pending or cancelled. Confirm scheduler startup and review the first receipt.
5. Collect all 12 terminal cases. Compare first-birth flow and at-risk mass at date 0 and later against each arm's same-horizon 80% control; require 12/16 overlap and report terminal, market, fiscal and population/entry-queue gates. Buyer closing-LTV and matched financial-access diagnostics are a separate saved-policy packet and must not be called an observed mortgage LTV or a birth response.

The exact launch matrix and limits are in `plan.json`. The 10:00 deadline is author requested. Failed or incomplete cases remain explicit; they are not silently replaced by an easier closure.

The approved **source-only v2 staging** completed on October 2. The local and Torch SHA-256 values match for mechanism archive `4a3f0e0a1efe2861c78630e33c280afc963a6d9af56bb5805138b656539fcff5` and buyer archive `aaa7fdad1fc06a5eb24af654a46801e1304f60de17f216c6b652ac766efb51e3`. The [deployment receipt](deployment_v2_receipt.json) records the inventory hashes. Both mechanism source preflights and both buyer financial-map preflights passed with zero model solves; no selection manifest or model job was published. The old archives remain intact. To repeat the two preflights on the staged source, use:

```sh
ssh -o BatchMode=yes torch 'bash /scratch/td2248/projects/purchase_mechanism_v1/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism_deployment/preflight_torch.sh'
ssh -o BatchMode=yes torch 'bash /scratch/td2248/projects/purchase_buyer_diagnostics_v2/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/preflight_torch.sh'
```

The reviewed October 2 source is staged separately at
`/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a` (archive
SHA-256 `051c42a6423fe02396295566f921ca1781693170306c67cfaf4791c046ae2f7a`,
inventory SHA-256 `7cd28f68e96a18b24b3bea9e9b709038576444de20032a9538a82c37574247a3`).
Both zero-solve source preflights passed. The separate, unchanged v1 selection
store contains final best-available experimental checkpoints, hard Torch restart
chain 11 and quarter local chain 54; its manifest SHA-256 is
`19db8dfcc928c4cdea5810f70c4e62a9de65372bd038646d0bc2d49532913879`.
`launch_torch_reviewed.sh` and `submit_torch_reviewed.sh` route reviewed source
and results to the new root while mounting that v1 selection store read-only;
the exact route and source hashes are in `reviewed_routing_receipt.json`.
Final one-date native smoke array 19022272 and the dependent 12-case array
19022273 were submitted under `afterok:19022272`. These are experimental
best-available fits, not accepted paper baselines; case acceptance still depends
on the saved native receipts and numerical gates.

The separate 48/64-date numerical-horizon extension is prepared at
`/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1`. Its
source archive SHA-256 is `94a1962be2ba0142a1bfc71228c4f18ee16e48e968d9f8f26538aae2586f34d2`
and inventory SHA-256 is `a4559bf682de6e15461e6e35df595769e515d6b1c5c236900ff8fc46b9dbc24e`.
Both arm-specific zero-native preflights passed. The extension preserves the
same economics and acceptance gates, uses the v1 selection store read-only,
and writes distinct results. Its 12-case routing, per-case four-hour/1,024-call
limits, absolute 10:00 ET deadline, and external script hashes are recorded in
`extension_routing_receipt.json`. It has not been submitted as of preparation.
