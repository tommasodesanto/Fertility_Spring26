# Fixed current-child-value shift diagnostic

This is a saved-policy, zero-solve diagnostic, not a recalibration or a new
equilibrium. It raises the current value of successfully having the first
child by `delta` in both the baseline and permanent-both-100 credit arms, while
holding their saved continuation values and every other policy fixed. Since
the saved first-child flow benefit at `m=1` is `psi_child`, rows correspond to
fixed flow levels `psi_child + delta` in `[0.01, 0.10, baseline, 0.25, 0.50]`;
they are not full recalibrations of `psi_child`.

For positive interior attempt probabilities, the prescribed update is
`logit(p_try_new) = logit(p_try) + pi * delta / kappa_fert`. Saved probabilities
of exactly zero or one stay at that endpoint under any finite odds shift; this
diagnostic does not invent probabilities for them. The retained weight is the
saved young, childless renter PRE mass, including states with zero probabilities.

Inputs are the authenticated baseline and permanent-both-100 saved policy
arrays and the common `q0_reference_inherited_states.npz` PRE distribution
used by `ablate_current_two.py`. There are no model solves, no new prices, no
changes to the inherited PRE distribution, and no adoption claim. The delta=0
row is required to reproduce the previously reported credit-minus-baseline
first-birth response of `-0.051680757225407` percentage points.

## Results

The probability columns are expected first births per unit of the unchanged
young, childless renter PRE mass (0.1863761723 total). Positive and negative
mass columns decompose the credit-minus-baseline first-birth mass across PRE
states; their sum is the net mass.

| Current child flow value | Baseline probability | Credit probability | Credit minus baseline (pp) | Positive mass | Negative mass | Net mass |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.11782479 | 0.11767267 | -0.01521187 | 0.000027246 | -0.000055597 | -0.000028351 |
| 0.10 | 0.17298593 | 0.17266466 | -0.03212722 | 0.000030804 | -0.000090681 | -0.000059877 |
| 0.171561928 (baseline) | 0.22237400 | 0.22185719 | -0.05168076 | 0.000029924 | -0.000126244 | -0.000096321 |
| 0.25 | 0.27831321 | 0.27754167 | -0.07715411 | 0.000025691 | -0.000169488 | -0.000143797 |
| 0.50 | 0.43936843 | 0.43794295 | -0.14254866 | 0.000009689 | -0.000275365 | -0.000265677 |

Raising the current child benefit raises first births in both saved-policy
arms. It does not reverse the within-state credit gap, and the aggregate
credit-minus-baseline response remains negative throughout these five rows.
This fixed-current-value exercise does not rule out a full `psi_child`
recalibration changing future choices and continuation values. The underlying
full 14-target, 31-parameter, and 17-plot packets are linked in the
[parent saved-value README](../README.md#retained-artifacts).

Run with `python run_child_value_shift.py`. Results, input hashes and numerical
checks are written to `child_value_shift.csv` and `summary.json`.
