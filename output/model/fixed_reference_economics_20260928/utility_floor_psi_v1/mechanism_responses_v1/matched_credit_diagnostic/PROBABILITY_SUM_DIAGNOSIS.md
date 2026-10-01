# Saved probability-sum diagnosis

The one authorized metadata pass completed in 0.420 seconds with maximum
resident memory 75,460 KiB. It loaded only saved first-birth probabilities and
the baseline pre-fertility weights on Torch. No model calls, array downloads
or Slurm jobs occurred. Full compact evidence is `probability_sum_diagnosis.json`.

The reference-credit array has 22 nonzero action sums whose absolute error
exceeds 1e-12; the largest error is 1.343143221e-7. All 22 states have exactly
zero baseline childless pre-fertility mass. Every age has zero occupied failing
states and zero affected probability mass. The maximum error on occupied
reference states is 3.441691376e-15. The expanded-credit array has no failing
state anywhere; its maximum error is 1.887379142e-15. Both arrays are float64,
finite and within [0,1]; the unused actions are exactly zero, so summing all
four actions does not explain the discrepancy.

The largest reference error occurs at wealth-index12, inherited-tenure-index3,
location0, age-index6 and income-index2. Its wait/try probabilities are
0.9999990463261383 and 8.193595394991998e-7, summing to 0.9999998656856779.
Its baseline childless pre-fertility mass is exactly zero. This identifies a
saved-probability discrepancy outside the weighted diagnostic's occupied
support; it does not establish its computational cause.

Proposed rule for lead review: preserve the 1e-12 tolerance and global
finiteness, bounds and action-identity checks, but require the sum gate on
positive baseline pre-fertility matched-state mass. Count and disclose the
unoccupied violations separately. The recovery mask would still require
valid sums, interior probabilities and finite non-dead inclusive values,
without normalization or clipping. This changes the gate's support to match
the question; it does not certify global support or weaken precision.

The original reducer remains unchanged and its failed attempt is preserved.
No full reduction retry has been made.
