# Winner31 fixed-credit-policy distribution accounting

At the winner31 prescribed price \(q_0=0.719168368828958\), changing from the inherited reference pre-fertility distribution to the expanded-credit stationary distribution reduces birth events by **0.008851894096 per household per four-year period** under the *same saved expanded-credit fertility policy*. This is the distribution component of the previously verified three-corner accounting, not a new policy or equilibrium comparison.

Define a group \(g\) by age, income state, inherited housing tenure/product, location, children ever born before the period, and children currently at home. Let \(M_{kg}\) be group mass, \(p_{kgb}\) its conditional distribution across net financial asset node \(b\), and \(f_C(g,b)\) the saved expanded-credit birth expectation at that state. For baseline \(B\) and credit cohort \(C\), the symmetric two-order decomposition on groups present in both distributions is

\[
W_g=\tfrac12(M_{Cg}+M_{Bg})\sum_b f_C(g,b)(p_{Cgb}-p_{Bgb}),\qquad
G_g=\tfrac12(M_{Cg}-M_{Bg})\sum_b f_C(g,b)(p_{Cgb}+p_{Bgb}).
\]

Here \(W_g\) records changes in financial-asset distributions *within* a fixed group, and \(G_g\) records changes in group mass. If a group were absent on one side, its observed birth flow would be reported separately as unmatched support without imputing a missing conditional rate. All 4,599 occupied groups were present in both distributions, so that term is zero. The model has one location state in this packet.

| Outcome | Within-group asset distribution \(W\) | Group mass \(G\) | Unmatched groups | Total distribution term |
|---|---:|---:|---:|---:|
| All birth events | −0.008259773500 | −0.000592120596 | 0 | −0.008851894096 |
| First births | −0.007221411500 | +0.001817705633 | 0 | −0.005403705867 |

The within-group asset-distribution term accounts for 93.3% of the total birth-flow decline in this ordering. For first births, the group-mass term offsets part of a larger within-group decline. This localizes the accounting shift; it does **not** show that debt, mortgage equity extraction, or saving caused the fertility change. Net financial assets include mortgage positions, and the stationary distribution is jointly endogenous to prior fertility, tenure and saving choices.

| Age cell starts at | All-birth asset term | All-birth group-mass term | All-birth total | First-birth asset term | First-birth group-mass term | First-birth total |
|---|---:|---:|---:|---:|---:|---:|
| 18 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |
| 22 | −0.001637 | +0.000155 | −0.001482 | −0.001570 | −0.000221 | −0.001791 |
| 26 | −0.002299 | +0.000174 | −0.002125 | −0.002137 | +0.000202 | −0.001935 |
| 30 | −0.002053 | +0.000049 | −0.002004 | −0.001814 | +0.000639 | −0.001176 |
| 34 | −0.001342 | −0.000185 | −0.001527 | −0.001090 | +0.000671 | −0.000419 |
| 38 | −0.000643 | −0.000431 | −0.001074 | −0.000438 | +0.000358 | −0.000080 |
| 42 | −0.000286 | −0.000355 | −0.000641 | −0.000171 | +0.000169 | −0.000002 |

All table entries are birth events per normalized household per four-year period. Ages 46 and above have zero births under the saved policy. The age-18 distribution component is zero because the entrant distribution is fixed across the two arms.

The remote pass read only the existing winner31 q0 inherited PRE array and credit saved solution; it did not download arrays. It authenticated the same 31-field parameter table and source binding as the prior six-cell run, reconstructed credit PRE with POST nesting L1 below \(1.5\times10^{-16}\), and used the installed native birth operator on origin family cells. The reference PRE hash matched exactly. One CPU and one thread were used; execution took 87.23 seconds, peak memory and solver traps are recorded in `receipt.json`. Household, Bellman, lifecycle, optimizer and GE solve counts are zero. Both PRE masses equal one within floating-point precision; 4,599 groups have positive mass on both sides; the maximum credit-policy state birth-rate difference on overlapping positive-mass states is \(3.33\times10^{-16}\). First-birth and all-birth group and age sums reproduce the saved scalar distribution components to numerical precision (`checks.json`). The original credit solution remains support-limited: unoccupied alternatives and grid convergence are uncertified.

`groups.csv` contains the complete \(g\)-level decomposition; `age.csv` contains the age summaries; `receipt.json` and `checks.json` give the exact identity and reconciliation gates. `reduce.py` is the isolated saved-data reduction, and `run_remote.sh` records its one bounded execution. No new economic intervention was run.

For the repository backup, the two delivered CSVs were converted from CRLF to LF line endings only. `checks.json` hashes the delivered files; the executed reducer and numerical rows are unchanged.

## PRE identity and scope across packets

This reduction pins the **original six-cell response** reference PRE array, SHA-256 `c6589a65ad74e624579c7abe2e4069ab4e07712f299496deefda6ec7a6d5e48a`. Its provenance is [`collected/responses/q0_reference_inherited_states.json`](../../collected/responses/q0_reference_inherited_states.json): the reference-price-1.0 case's reconstructed stationary pre-fertility distribution. The original [`fixed_price_responses.py`](../../fixed_price_responses.py) reconstructs that distribution at line 210 and saves it at lines 256–260. The baseline and credit q0 impact closures also carry the same inherited-distribution hash. `receipt.json` records the **array-content** hash checked in the remote pass, not merely an NPZ-file checksum.

The later mortgage and matched-preference packets use a **separately generated** reference PRE array, SHA-256 `459ab9229e7ce58376bda8c933b583ec059104ede31c57627ce7368c41e5a51b`, recorded in [`purchase_ltv_v1/local_run/retry5/results/q0_reference_inherited_states.json`](../../purchase_ltv_v1/local_run/retry5/results/q0_reference_inherited_states.json) and the [matched-preference note](../matched_preference_credit_v1/README.md). Both metadata files name a reference-price-1.0 reconstructed stationary PRE and report the same shape, but their content hashes differ. Some aggregate outcomes match to numerical roundoff in the retained reports; **elementwise or bitwise equivalence of these two PRE arrays has not been established**. The wealth decomposition here must not be represented as an array-matched decomposition for the mortgage or matched-preference packet.
