# CES normalized-share overnight calibration

This is a prepared, non-adopted four-chain bounded Nelder--Mead search for the author-authorized normalized CES-share experiment. It starts from post-interest soft chain 13 and retains the old wealth target (6.92658379107299), timing, earnings, and all other inherited economics. The experiment has 14 target rows, 11 scored moments, and 11 free coordinates: the nine inherited coordinates plus `delta_alpha_jump` and `delta_alpha`, each in `[0,.25]`. The added `family_rooms` target is 0.38509964969278165 with weight 280.52808370152104.

The share rule is `alpha(m)=.733` when childless and `clip(.733-delta_alpha_jump-delta_alpha*m,.05,.95)` for parents. The normalized material denominator applies in all states and `h_P=0`. The existing birth menu, raw utility costs, no-estateA restriction, grids, fixed inputs, timing, and earnings remain unchanged. There is no `r*` correction or `alpha0` numerator, and no added birth shock or cost rescaling. The `family_rooms` weight is the inherited 42-metro bootstrap weight because national uncertainty is unavailable; the model-dependent-child observer remains a proxy.

V5 source inventory SHA-256 is
`39fe1c3316919fb3d8672030afda6fb5fcc29aef0969181d206c9fc0757ca127` (3,008
files). Its three-context preflight, smoke **19132940**, full native gate, and
collector passed. Four-chain production array **19133352** (`0-3%4`) was
submitted October 3 at 22:31:16 EDT. At 22:32:52 EDT all four tasks were
verified `RUNNING` on cs604/cs606/cs633 ([launch health](../../../output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt5/launch_health.json)). All four searches have finished and passed fresh native verification; [full final fits and parameters](../../../output/model/experiments/ces_normalized_shares/overnight_v1/final_results/README.md) are retained. The experimental calibration has not been adopted.

V4 smoke **19132298** passed its two-GE native gate and collector, then Slurm
reported `FAILED 1:0` because its EXIT receipt function lacked an `os` import.
V5 smoke likewise passed its numerical gate and collector, then failed in the same terminal bookkeeping.
V5 leaves model, search and budget source unchanged; a separate reviewed
launcher repairs only EXIT-receipt import, quoting and valid-JSON newline. Four
local and Torch fixtures passed and preserved JSON validity and exit codes. V5 derives
from the v3 parent by hardlink, preserving the six path differences; its three
reference files are in the derived source, so no dependency overlay is needed.
The [experiment packet](../../../output/model/experiments/ces_normalized_shares/overnight_v1/README.md)
links the smoke, launcher, source and submission receipts.

Per chain: 1 CPU, 24 GiB, six hours, at most 500 objective calls and 32
lifecycle solves per GE. A native GE may start only with more than 2,700 seconds
left: a 900-second search-GE window plus the 1,800-second final-verification
reserve. Warm native GEs took 177 seconds; cold GEs took 267–270 seconds. The
two-candidate smoke plus fresh postcheck took about 12.5 minutes; observed speed implies roughly 70–105 GEs
per chain, so wall time is expected to bind before the 500-call cap. An
authenticated numerical candidate that cannot be bracketed receives a
`1e12` penalty; a bounded-budget stop postchecks the best candidate. Other
errors are terminal without retry or fallback. Native verification compares
11 coordinates, full target and parameter tables, all 17 diagnostic plot
hashes, and an exact repeat.

## Fixed-price mortgage-financing diagnostic

The author-requested two-case diagnostic uses the verified best overnight chain, holds its price and all fitted preferences fixed, and changes the uniformly financed mortgage share from 80% to 95%. Run `bash code/cluster/ces_normalized_shares_calibration/credit_diagnostic.sh`. The driver authenticates the frozen v5 stage and the chain-1 JSON, checks the reached utility arrays and that only phi differs, and saves the standard 17 graphs and full arrays for each case. One CPU, 24 GiB and a 30-minute cap. See the [diagnostic packet](../../../output/model/experiments/ces_normalized_shares/credit_diagnostic_v1/README.md) for status, interpretation and limitations.
