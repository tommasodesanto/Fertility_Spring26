# Policy impact check from the verified simultaneous-choice calibration

Question: did the improved calibration also change the fertility response to
family credit and housing supply? This diagnostic uses the selected overnight
candidate (loss 23.791955301663187), without further calibration. Its full
fit and parameter bounds are in the September 8 morning review in the sibling
overnight output folder.

The test holds the inherited 2023 pre-choice population fixed within this
calibration and clears each policy housing market. Baseline replay comes first;
then baseline, 95% LTV for households with dependent children, and a 20% upward
shift in the housing supply schedule. Existing policy definitions are reused.
The tax rate remains 1% annually, no grant or rebate is added, and the dated
supply elasticity remains externally fixed at 0.63. The housing taste scale
remains externally fixed at 0.005. Other estimated parameters and household
choice equations remain unchanged. Wealth grid: 120 nodes.

This is a 2023 impact diagnostic, not a new calibrated target system, welfare
comparison, 2063 forecast or production promotion. No future entry, migration,
retention or fiscal closure is selected by this test. Those must be reconciled
before any full transition. No unsolicited figures or monitoring automation.

Execution: one baseline smoke (fresh fixed-price replay and a cleared baseline),
then three independent cleared policy cases only if smoke passes. This is five
evaluations in total, each market potentially requiring multiple household
solves. Recent full five-date calibrated histories took about 33–44 minutes;
we do not assume policy solve times are identical. Each CLI case has a hard
30-minute cap and Slurm a 35-minute cap, one core and 24 GiB. Smoke and policy
array are separate jobs linked by successful completion. No automatic retries
beyond the unchanged market solver's existing two-pass routine. Failed gates
block comparison; no checks are relaxed. Each case writes a minute heartbeat,
latest completed summary, frozen-best description, quantities, checkpoint and
budget/probability/value diagnostics. The smoke must verify the same execution
and output loop used by policy cases before the array can start.

Report impact percentage changes in births per household, mean rooms and price,
and percentage-point changes in ownership, each relative to the new-model
baseline. Historical policy numbers remain historical until a matched reference
is explicitly checked. Larger fertility effects are a hypothesis to test.

Submitted September 8: smoke **17222216**, dependent three-case policy array
**17222217**, and comparison collector **17222218**. Outcomes are pending;
these are queued tests, not evidence of stronger policy effects.
