# Attempt status: blocked before incidence extraction

The single authorized Torch accounting job, **18832560**, failed after 36 seconds with exit code `1:0` (peak RSS 2,081,212 KiB). No renter-taper CSV, owner-to-renter sale table, entrant summary, or PASS receipt was produced. The checkpoint was authenticated and the staged driver/launcher SHA256 checks passed before runtime setup.

The failure is the gate at `audit_renter_taper.py` line 62: it requires every element of `policy.loc_probs` to equal one. The saved policy array contains at least one non-one element, so execution stopped before any incidence calculation. That global check may include unreachable/infeasible states; the proposed next diagnostic is to inspect location probabilities only on positive `g_post_fertility` states and verify unit stay probability there. This correction is untested and was not applied to the immutable `source_v1`; no second job was submitted.

Artifacts here are the compact [Slurm log](slurm-18832560.out) and [scheduler record](sacct-18832560.txt). The original source and its hash manifest remain preserved in `../source_v1/` and on Torch. No model solves ran and no economic change was made.
