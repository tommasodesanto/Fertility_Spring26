# Inherited 2023 property-tax continuations

Author authorized an unexpected permanent increase from 1% to 2% annual
property tax in the inherited 2023 economy. Two independent Torch pipelines
solve matched 1% and 2% continuations; all property-tax receipts are rebated
equally to household heads and PAYGO balances at payroll tax 0.179.
Preferences stay at 0.09221854783921073. Each path has 100 four-year dates,
2023--2419, and a separately verified stationary continuation in 2423.

Both arms inherit the saved 2023 pre-choice distribution behind the current
one-permanent-shock iteration3 figures. The four adjusted/raw birth vintages
are carried from 2007/2011/2015/2019, each divided by 2.1. There is no
immigration, population rescaling or historical age reweighting. This starting
history is unconverged: policy differences remain provisional even if the
subsequent finite paths converge. The announced-four-shock run remains separate.

The supply curve stays fixed in asset-price coordinates, with elasticity 0.63.
The policy changes the annual tax from .01 to .02 (four-year .04 to .08) and
updates the native user-cost identity. No structural parameter, housing-supply
scale, target, uncertainty estimate or numerical gate changes.

Each pipeline: exact saved-2023 reconstruction and native budget checks;
fresh baseline endpoint verification or new 2% endpoint solve; two-date joint
root-loop smoke with its own fixed-price household continuation; then up to
two eight-evaluation rounds of the 100-date price/pension/rebate root.
The second round resumes the previous numerical best/Jacobian only if needed
and time remains. Up to 24 terminal evaluations/one hour, 25-minute smoke,
eight-hour path budget and ten-hour per-job hard limit. At about 32 minutes
per full mapping, allow roughly 4--9 hours per independent arm.

One CPU/32 GiB per job, native numerical threads fixed to one. Each completed
mapping saves rows, fertility, native gates, terminal distances and PNG/PDF.
The first long mapping includes the full standard diagnostic gallery. Per-arm
latest-completed/best-so-far, minute heartbeat, phases and failure receipts
remain available. One arm failing does not stop the other. Terminal approach
and sensitivity to the chosen horizon remain separate acceptance conditions.

Driver: `code/cluster/run_e5f_inherited_2023_tax_long.py`.
Remote batch: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/inherited_2023_tax_20260913a`.
Local submission and pinned manifest are saved alongside this note.

## Expanded to both histories

The author requests policy results for both one-shock and four-shock histories.
The existing four-shock baseline is reused. A queued recovery will freeze its
best completed104-date evaluation after its job ends, recover its own2023
household state and birth queues, and automatically launch only the2% tax arm.
The2% terminal equilibrium is common across histories and is freshly audited
before reuse. Each policy comparison must stay within its own inherited history.

The one-shock2% job17722295 passed terminal verification and its two-mapping
smoke and entered the100-date transition. One-shock1% job17722294 failed strict
saved-value reproduction after recomputing the unchanged user-cost expression.
Retry17722962 in sibling batch `inherited_2023_tax_20260913b` preserves the
exact stored baseline floating-point value. Scientific gates are unchanged.
