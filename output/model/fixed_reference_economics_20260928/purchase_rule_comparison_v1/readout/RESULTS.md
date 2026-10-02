# Four fixed-price purchase-rule diagnostics

All cases hold the strict-80 price, $H_0$, $N_0=1$, and ten calibrated coordinates fixed. This is a partial-equilibrium comparison, not recalibration.

| Case | Purchase eligibility | Saving rule | Gate | Loss |
| --- | --- | --- | --- | ---: |
| hard80 | strict wealth-only; $\phi=0.8$ | standard | passed | 217.20688594292585 |
| hard100 | strict wealth-only; $\phi=1$ | standard | passed | 256.03658589324806 |
| quarter80 | strict wealth-only; $\phi=0.8$ | quarter-saving constraint | passed | 100.26620364054745 |
| quarter100 | strict wealth-only; $\phi=1$ | quarter-saving constraint | passed | 192.40916254556313 |

The ten estimated coordinates are fixed at their strict-80 values in all four cases. `delta_alpha_jump` is externally fixed at zero; any inherited source-table status calling it free is corrected in the display table. The failed first $\phi=1$ attempts reached the negative-estate production gate before target-fit and parameter reports were serialized. Their missing cells are blank, not estimated. Successful isolated solvency retries populate those columns when collected.

Complete 14-moment target fit: [target_fit.csv](target_fit.csv). Complete 31-parameter comparison: [parameters.csv](parameters.csv). Source hashes and gate states: [verification.json](verification.json).
