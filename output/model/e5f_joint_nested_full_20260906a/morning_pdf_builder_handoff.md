Rebuilt the evidence-driven report builder: [build_e5f_joint_nested_review.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_review.py).

It now:

- Requires a complete selected candidate and validates all 12 fit rows, 11 free + 3 fixed/derived parameters, arithmetic, numerical gates, and original artifact hashes.
- Accepts optional policy/search-verification/narrative inputs.
- Validates same-anchor policy provenance, four cases, all dates, 17 graphs/date, gates, and recorded effects against paths.
- Distinguishes two-date smoke from complete 44-date paths.
- Requires an explicit label for smoke inputs and never calls a candidate a final calibrated optimum.
- Writes an adjacent verification JSON with `visual_review_pending`.

Example future use:

```bash
python3 code/model/tools/build_e5f_joint_nested_review.py \
  --selected-dir /path/to/selected/task_001 \
  --policy-results /path/to/policy_results \
  --search-verification /path/to/final_verification.json \
  --narrative /path/to/lead_narrative.json \
  --output /path/to/morning_readout.pdf
```

Smoke preview generated and rendered (21 pages); its receipt records source checks and explicitly labels it as a two-date smoke: :codex-file-citation{path="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/tmp/joint_pdf_builder_qa/e5f_joint_nested_smoke_preview_v2.pdf" purpose="output"}

Verified: syntax compilation, PDF rendering, missing-input rejection, unlabeled-smoke rejection, and receipt/hash validation. No scientific code, outputs, commits, or cluster state changed.