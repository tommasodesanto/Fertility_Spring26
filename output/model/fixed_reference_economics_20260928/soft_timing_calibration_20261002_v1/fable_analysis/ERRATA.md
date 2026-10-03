# Corrections to the immutable Fable launch prompt

The active Fable session was dispatched at 06:06 New York on October 3. Its saved `PROMPT.md` is hashed in `launch_receipt.json` and has **not** been altered. Apply these corrections when reviewing its analysis and, if the agent sees this file during its run, in its deliverable:

1. The early-fertility empirical moment is mean children ever born (capped at three) among women whose **completed interview age is 25**, corresponding to ages $[25,26)$, not the instant of the 25th birthday. The model uses uniform within-period birth-time interpolation in $[22,26)$ with post-birth weight $0.875$. Source: `calibration_archive/context_refresh_20261002/target_provenance_review.json`, early-fertility record, especially the definition and warning.
2. The third experimental arm changes **only the empirical wealth/earnings calibration target** from $6.92658379107299$ to $4.45838713455674$, retaining its weight $7.595098472533724$. The entrant wealth distribution is unchanged. Source: `output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/README.md` and `collection/RESULTS.md`. Call it the **new wealth-target arm**, not an alternative-entry-wealth arm.

Do not interpret a difference between the target systems as an interest-timing or entry-distribution mechanism. These are wording/source-identity corrections, not adopted economic changes.
