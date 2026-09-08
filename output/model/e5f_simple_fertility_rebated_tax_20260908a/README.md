# Priority: rebated property-tax reform under the new benchmark

The author corrected the policy priority to the property-tax experiment.
The maintained comparison is annual property tax 1% with an equal rebate
versus 2% with an equal rebate, using the same inherited 2023 population from
the verified simultaneous fertility-nest calibration. The unrebated 1% fitted
benchmark is a separate reproduction control. Full calibration fits and bounds
remain in the sibling overnight folder's morning review.

The existing coupled root jointly solves the housing price and equal transfer.
Property-tax revenue covers occupied rental and owner housing services at the
asset price; the transfer divides revenue by current household decision units,
not children or owners alone. Four-year model tax rates are .04 and .08.
The dated housing supply schedule and elasticity .63 are unchanged. New nest
flags, housing taste scale .005 and all estimated parameters remain fixed.
This is one-date impact analysis; no future entry, migration or retention law
is selected and no welfare or 2063 claim is made.

First verify the completed unrebated baseline smoke, then run the 1% rebated
case through the entire coupled-root, restored-parameter, fresh fixed-price
replay, audit, checkpoint and receipt loop. Only success permits the 2% reform
through that same loop. Restore the selected transfer explicitly because the
root caches trial evaluations and can leave the mutable parameters at a different
trial. Independently recompute the fiscal ledger after restoration and replay.

Root gates: normalized joint residual at most 1e-4 and absolute fiscal gap at
most 2.5e-5; market residual at most 2e-4. Existing stronger mass (2e-10), budget
violating mass (2e-10), probability, feasibility (1e-6), and occupied value
monotonicity checks remain enforced. Fresh replay must reproduce solved arrays,
distributions and births within 2e-10. No tolerances are relaxed.

Two coupled markets plus one fresh fixed-price replay per market; each coupled
root uses its existing 16-iteration cap and can cache up to roughly 99 trial
evaluations. Recent baseline replay plus market solve took 59 seconds; tax roots
may take several minutes. Each process has a 30-minute hard cap and a 35-minute
Slurm limit. Allocate eight CPUs and 128 GiB per tax job to meet scheduler memory
allocation requirements; numerical computation remains one thread. This memory
reservation accommodates cached policy arrays. Jobs run in dependency order, with no automatic retries after
failure. A comparison collector runs only after both pass. Heartbeats every
30 seconds, latest completed and frozen-best summaries, checkpoints and fiscal
ledgers preserve progress and failures. No figures or monitoring automation.

The first readout reports net changes in births per household, housing services,
asset prices, ownership, young ownership and equal transfers between the two
rebated equilibria. It does not attribute those net changes to separate causal
components; that requires the separate eight-cell tax/price/rebate decomposition.

Submitted September 8: rebated-baseline smoke **17222674**, tax2 reform
**17222675**, collector **17222676**. All were queued at submission; no new
tax effect is established yet. The first single-CPU/high-memory request was
rejected before job creation; the eight-CPU allocation passed scheduler validation.
