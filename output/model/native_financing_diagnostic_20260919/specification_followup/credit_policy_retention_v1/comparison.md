# Credit policy retention collection

- Smoke job `18078707`: **FAILED** after 1:24, exit `1:0`, peak RSS `2420980K`.
- Production job `18078708`: **CANCELLED** by dependency after smoke failure; no restart or cancellation was performed here.
- Pre-case verification passed: checkpoint, frozen solver/parameters, copied source files, helper hashes, plan, summary, and wrapper hash gates all reported `ok=true` (`1102` hash records).
- No policy-retention comparison, population gate, numeric gate, or 17-plot manifest was produced. Large arrays were not downloaded.

The failure occurred before any retained policy result: `import origin violation: __mp_main__ loaded from .../run_e5f_credit_policy_retention.py`. The dependency state is therefore recorded as unavailable rather than interpreted.
