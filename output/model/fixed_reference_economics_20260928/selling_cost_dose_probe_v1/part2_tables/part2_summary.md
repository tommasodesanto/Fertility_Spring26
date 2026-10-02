# Part 2 lock-in diagnostic — childless (n=0, m=0), age cells 0-4

Cases (saved arrays, zero new solves): psi06_phi08, psi06_phi10, psi00_phi08, psi00_phi10.
At-risk mass g_pre: reconstructed exactly as the observer does --
entrant cohort at j=0 from saved entry_by_loc; advance of saved
post-fertility g_beginning_distribution otherwise (saved loc/tenure/
saving policy, rebuilt location/tenure maps).
Attempt prob: at-risk-mass-weighted mean of saved fert_probs parity-0
try slice [...,1] (second axis = beginning tenure state to) -- the
exact engine birth-flow factor. Birth prob = attempt x pi_j; per-cell
pi_j in the prob CSV.
Fecundity vector pi_j (j=0..16): 0.98, 0.965817, 0.941576, 0.900144, 0.82933, 0.708298, 0.501436, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.
Childless mass at (n=0, cs>=1) in saved post-fertility arrays is 0.0
in all 4 cases (checked).
Engine flags: due_stayer=True, entry_censor=True, independent_child_maturation=True, joint_nested_choice=False, parent_age_maturation=False, readiness_gate=False, sequential_births=True, use_age_survival=True, use_stochastic_aging=True
Tenure labels: to=0 renter; to=1..5 owner rungs H_own=(2,4,6,8,10).
Verification per case in part2_verification.csv: post-fertility
reconstruction L1 vs saved arrays; implied aggregate first births
vs observer first_birth_flow; entry mass vs saved entry rate.

## Post-birth tenure transitions

No post-birth-branch tenure policy is saved: the only tenure arrays are pre-tenure policies by beginning tenure state (tenure_choice/tenure_probs indexed by pre-tenure parity/child-state). Keep/upsize/downsize/to-rent shares after the birth branch cannot be computed from saved arrays.
Tenure-related saved keys: tenure_choice, tenure_probs.
Full per-case key inventories are identical: True

No economic interpretation per TASK.
