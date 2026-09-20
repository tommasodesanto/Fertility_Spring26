# Credit policy retention collection

- Smoke job `18078707`: **FAILED** after 1:24, exit `1:0`, peak RSS `2420980K`.
- Production job `18078708`: **CANCELLED** by dependency after smoke failure; no restart or cancellation was performed here.
- Pre-case verification passed: checkpoint, frozen solver/parameters, copied source files, helper hashes, plan, summary, and wrapper hash gates all reported `ok=true` (`1102` hash records).
- The household solve, cohort outputs, numerical checks, retained policy arrays and 17 plots were produced before the final receipt failed. No completed cross-dose policy comparison was produced. Large arrays remain remote.

The failure occurred during final import-origin receipt assembly: `import origin violation: __mp_main__ loaded from .../run_e5f_credit_policy_retention.py`. The alias pointed to the exact pinned driver. This attempt consumed one household solve; its production stage never ran. The subsequent v2 experiment documents the focused import-alias correction and preserves all scientific gates.
