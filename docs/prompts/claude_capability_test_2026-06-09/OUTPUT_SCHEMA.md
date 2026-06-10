# Output Schema

Claude should produce one main report and any optional scratch artifacts under:

```text
output/claude_capability_test_2026-06-09/
```

The main report should be:

```text
output/claude_capability_test_2026-06-09/claude_full_pass_report.md
```

Required report sections:

1. `Executive Judgment`
   - five to ten bullet points maximum
   - strongest current asset
   - most serious weakness
   - highest-return next action

2. `Evidence Log`
   - files read
   - commands run
   - tests or smoke checks run
   - checks not run and why

3. `Theory Audit`
   - primitives and timing
   - household problems
   - competitive equilibrium definition
   - constrained planner problem
   - proposition-by-proposition assessment
   - notation/style issues
   - corrected theorem/proposition statements where needed

4. `Math Derivation Check`
   - derive the household FOCs from the stated problem
   - derive the planner FOCs from the stated planner problem
   - compare CE and planner conditions
   - identify precisely which wedge is an externality, which is a private
     constraint, and which requires a policy instrument or transfer assumption
   - state whether the claimed efficiency result is valid as written

5. `Codebase Audit`
   - active model map
   - equation-to-code correspondence for key objects
   - likely bugs or inconsistencies with file and line references
   - missing tests or diagnostics
   - smallest useful smoke tests to add

6. `Calibration Audit`
   - target table: target, code location, live value, model value if cheaply
     available
   - objective and weighting assessment
   - identification and parameter substitutability
   - stale or incomparable outputs
   - next calibration experiment with budget, checkpoints, and stop criteria

7. `Project-Level Diagnosis`
   - organization/source-of-truth issues
   - documentation gaps
   - reproducibility risks
   - file hygiene issues that matter for future work

8. `Ranked Action Plan`
   - P0/P1/P2 issues
   - one-day actions
   - one-week actions
   - items that should not be touched yet

Optional scratch artifacts:

- `theory_clean_room_note.tex`: a clean 2--4 page standalone theory note if
  Claude can improve the compact analytical result without rewriting the whole
  draft.
- `code_findings.csv`: machine-readable list of code findings with path, line,
  severity, and recommendation.
- `verification_log.md`: fuller command output summaries.
- `patch_suggestions.diff`: proposed patch only; do not apply to active files
  unless the user explicitly approves.
