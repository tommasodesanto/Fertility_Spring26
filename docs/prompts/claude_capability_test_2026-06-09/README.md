# Claude Capability Test, 2026-06-09

This microfolder is a self-contained harness for testing a new Claude model on
the fertility project. It is intentionally prompt-and-rubric only: it points
Claude to the live repository sources instead of copying large code or paper
bundles into this folder.

Recommended use:

1. Open Claude in the repository root:
   `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.
2. Paste the contents of `MASTER_PROMPT.md`.
3. Tell Claude to use its strongest available reasoning mode.
4. Keep Claude's generated outputs under
   `output/claude_capability_test_2026-06-09/` unless you explicitly approve
   edits to active paper or model files.

Files:

- `MASTER_PROMPT.md`: copy-paste prompt for Claude.
- `CONTEXT_MANIFEST.md`: required startup files, active artifacts, and
  high-value paths to inspect.
- `OUTPUT_SCHEMA.md`: expected report structure and deliverables.
- `SCORING_RUBRIC.md`: a rubric for judging whether the model gave a genuinely
  useful research-engineering pass.
