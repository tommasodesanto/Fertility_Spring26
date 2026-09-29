# Age 25/26 and teenage birth diagnostic

This saved-data packet tests two proposed explanations for the young number-of-children gap. Torch jobs 18833274 and 18833422 passed. Both used one CPU, at most 4 GB and a 10-minute limit; neither imported or solved the model or read checkpoints. The original and two-birth selected candidates are experimental and do not replace the block0506 reference. No target, weight, parameter, observer, or source model was changed.

The model observer interpolates the [22,26) cell at completed CPS interview age 25 ([25,26)) with post-birth weight 0.875. At age 26 ([26,27)), it interpolates the next [26,30) cell with post-birth weight 0.125. The same ages are matched on the data side. Counts are children ever born capped at three; conditional counts divide by motherhood share.

| Age | Series | Data children/woman | Model children/woman | Data motherhood | Model motherhood | Data children/mother | Model children/mother |
|---:|---|---:|---:|---:|---:|---:|---:|
| 25 | Reference | .810 | .535 | .457 | .450 | 1.770 | 1.190 |
| 25 | Original selected | .810 | .530 | .457 | .448 | 1.770 | 1.185 |
| 25 | Two-birth selected | .810 | .606 | .457 | .438 | 1.770 | 1.385 |
| 26 | Reference | .923 | .615 | .524 | .496 | 1.762 | 1.240 |
| 26 | Original selected | .923 | .609 | .524 | .494 | 1.762 | 1.234 |
| 26 | Two-birth selected | .923 | .689 | .524 | .484 | 1.762 | 1.424 |

Moving the **matched comparison** to age 26 widens the children-per-woman gap: original selected −.279 to −.313; two-birth selected −.203 to −.234. Comparing model age 26 with data age 25 would instead change the object being fitted and is not the active target. The main miss at both matched ages remains children among mothers.

NCHS pooled 2003–2006 period **first births** include 511,141 at maternal ages below 18 (7.731% of 6,611,269), and 1,368,228 below 20 (20.695%). The latter includes ages 18 and 19, which are within the model's age support. These are shares of first births, not the share of women with pre-entry births, all births, or children per woman.

An independent read of the authenticated June 2004/2006 CPS partitions finds children ever born capped at three of .084 at completed age 17 (motherhood .055; 1,973 women), .107 at age 18 (.077; 1,803 women), and .199 at age 19 (.154; 1,781 women). It reproduces the approved age-25 .810 target to numerical precision. Age 17 is [17,18), so .084 is a nearby pre-entry stock comparison, not exactly the stock at the eighteenth birthday. These CPS ages are cross sections of different cohorts. We cannot subtract .084 from the age-25 gap as a cohort contribution or infer the effect of adding an entrant child state without a cohort-matched birth-history measure and a model experiment.

Machine-readable evidence: `age25_26_fit.csv`, `age_tail_summary.json`, and `teen_entry_stock.json`. Sources and SHA-256 checks are recorded in the JSON receipts. `age_tail_diagnostic.py` (one level up) and `teen_entry_stock.py` in this folder reproduce the calculations on Torch.
