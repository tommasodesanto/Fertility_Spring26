# Model-matched wealth local calibration, October 3

Experimental arm: adopted soft purchase financing plus post-interest housing
transactions. The only target change is `wealth_earnings`, 6.92658379107299
to 4.45838713455674, retaining weight 7.595098472533724 as a controlled
sensitivity. The target comes from the pooled 2005/07 PSID matched-numerator
calculation in `output/model/wealth_numerator_match_20261002/`; a new standard
error has not been estimated. The bequest target remains unchanged.

`start_plan.json` pins ten distinct historical/nearby parameter starts,
source-checkpoint identity, bounds, and the new complete target/weight hashes.
`native_smoke_chain0/` runs two exact optimizer evaluations and full selected
native verification. The native `target_fit.csv` uses the old empirical target
for traceability. The authoritative experimental table is
`target_fit_new_contract.csv` beside it; `wealth_rescore_receipt.json` checks
the one-row loss replacement. Native reports retain 31 parameter rows and 17
standard plots. Plot files are diagnostics of model states and markets; any
target markers in inherited plots refer to their original source and do not
override the experimental target table.

`overnight_run/` preserves an initial one-chain launch stopped within two
minutes when its schedule was found to risk starving four queued chains. It
is not a completed calibration and will not be retried automatically. The
corrected `overnight_two_wave/` run starts at most five simultaneous one-thread
local processes and requires seven GiB available memory before starting
another. Ten starts run in two waves, each with at most 2 hours 45 minutes and
250 objective calls;
the shared hard stop is 07:00 New York time, with 1800 seconds reserved for
fresh native verification. Queued chains that cannot start within the budget
are reported as such, never silently retried. Search losses remain provisional
until selected native verification passes.
