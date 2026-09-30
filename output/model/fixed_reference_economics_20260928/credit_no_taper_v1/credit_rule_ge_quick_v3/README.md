# Credit-rule GE quick v1 — prepared, not submitted

This packet is the closed stationary GE add-on for the two companion debt rules:
`ours` removes the renter taper/caps; `author` applies the same rule plus the
strict owner-to-renter negative-raw-proceeds gate. It retains the reference
preferences (including \(\psi\)), entry distributions, income, age/birth
timing, fiscal rule, absolute \(H_s(q)\) curve and `H0`, original 160-node
grid, LTV/grandfathering, actual PAYGO/nonnegative-estate gates, and native
16/20 queue. There is no fertility normalization or physical-stock closure.

At each price it sets \(N=H_s(q)/d(q)\) and uses the actual stationary
`adult_entry_adjusted_birth_children / entry_rate` renewal residual. The root
starts at \(q_0\), then tries +3% if births are excessive and -3% otherwise;
it expands in that direction to 8% and 20%, then uses safeguarded log-secant
inside a bracket. It stops without a fallback. Each case has at most eight new
lifecycle solves, 180 seconds per solve, one CPU/24 GiB, and a hard stop one
minute before `2026-09-30T01:00:25Z`. A root is explicitly **preliminary**:
no exact repeat is permitted under this cap.

`rule_interface.py` is the only bridge to the companion packet. It imports its
`weights` function and its `strict_tenure.py`; it neither recreates a solver
nor introduces economic fields. Each candidate writes `preliminary_numbers`,
`latest_completed`, `best_so_far`, `progress`, and bracket/failure evidence.
The 14-fit and 31-parameter tables are written per successful candidate; a
candidate satisfying renewal/PAYGO/native-queue gates renders the unchanged
17 standard plots only if time remains.

Stage this directory and the unmodified sibling `credit_rule_quick_v1` together
under `/scratch/td2248/projects/fixed_reference_credit_rule_ge_quick_20260929/source_packet/`:
the two resulting children must be named `credit_rule_ge_quick_v1` and
`credit_rule_quick_v1`. The launcher refuses either missing source.
The prepared (not submitted) command is:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/credit_rule_ge_quick_v1/launch_credit_rule_ge_quick.sh
```

Local pure-algebra check: `python -m pytest test_root_math.py`. It does not
import the lifecycle runtime or perform a model solve.
