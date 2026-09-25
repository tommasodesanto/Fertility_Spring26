# Adopted maturation and adult-entry implementation: review packet

This is a source patch and hand-accounting receipt for the September 24 author decision. No local model import, test, solve, compilation, rendering, or bulk hashing was run on the Mac. The new path is explicit opt-in for the next closed, one-market calibration; existing saved results and the September 14 reference were not changed. `maturation.patch` contains the exact tracked-file hunks and new files. `independent_review.json` is the separate source-review receipt.

## Decision and units

Let (C_t) be children born during model date (t), counting the excess represented by the 3+ top bin. The household birth counter records (R_t), the explicit 0/1/2/3 births. If (T_t) households enter 3+ and the mean in that bin is (w_3), then (C_t=R_t+(w_3-3)T_t). Let (P_t=C_t/2.1) be **potential entrant households** from that birth cohort. The conversion is applied once before queueing. The household dependency state (m) still thins independently with probability (2/9) per four-year period; it never creates an extra queue entry.

At date (t), the transition's two waiting queues pop half of cohort (t-3) and half of cohort (t-4). The popped flow is injected at date (t+1):

\[
 B_{t+1}=\tfrac12 P_{t-3}+\tfrac12 P_{t-4}.
\]

Thus a date-0 birth impulse supplies entrant households at dates 4 and 5, 16 and 20 years later. The stationary prehistory initializes all three short-queue and four long-queue slots to half of the old (P). For constant (C_t=C), every popped flow is (P=C/2.1), and pending stock is (3.5P). With queue stock (Q_t), its exact mass identity is (Q_{t+1}=Q_t+P_t-B_{t+1}). In the pre-existing **open** transition, the already specified outside flow (M) and retention \(\rho\) act afterward, (E_{t+1}=M+\rho B_{t+1}). The new timing switch changes neither (M) nor \(\rho\). In the next **closed normalized** calibration, (M=0), \(\rho=1\), and stationary replacement requires (E=B=C/2.1); no new migration or geographic valve is introduced.

The closed stationary head distribution keeps the existing age-survival and total-mass normalization. Since the fertility intercept is normalized to 2.1 births per entrant household and survival is one through reproductive ages, the constant-birth cohort law has zero growth and (B=E). Split timing does not change this constant-birth stationary shape, but the stationary entry **accounting basis changes**: the old measured dependent-departure flow times 0.5 is no longer treated as birth-cohort renewal. The opt-in Markov solve reports `adult_entry_adjusted_birth_children`, `adult_entry_potential_total`, and `adult_entry_stationary_residual`. The final calibration driver must call `require_closed_stationary_renewal(E,C,tolerance)` **after** completing fertility-intercept normalization. The pinned commute objective's existing fertility tolerance is `0.0005` child units; the helper checks `abs(C/E - 2.1) <= 0.0005`, exactly equivalent to `abs(E-B) <= E*0.0005/2.1`. Intermediate fertility candidates need only report their residual.

## Selected receipt arithmetic, not a new solve

The September 24 selected receipt has (R=0.11627493581786193), (T=0.02218833829783723), (w_3=3.602359422009), so (C=0.12964029045028488). It gives (B=C/2.1=0.0617334716429928), versus existing normalized entrant mass (E=0.06173345618094337). The difference (E-B=-1.54620494255\times10^{-8}) is explained by the saved completed-fertility normalization (2.1000005259789205), within the existing tolerance. Old dependency-based mature flow is (0.05515116139576231), a different object. At this (B), each steady-state queue arm pops (0.0308667358214964) household units per period and combined pending stock is (0.216067150750475).

The active source uses age 18 plus four-year cells and `A_f_end=7`, so the last fertility cell is age 42. The selected serialized survival schedule is one through the age-62 transition; first mortality acts in the 66-to-70 transition. Offspring from age-42 births have both scheduled entries by parent age 62. Parent death still removes any remaining modeled dependency burden and triggers no added adult entry. This is a timing check, not a genealogical or literal support account.

## Files, activation, and verification

- `code/model/intergen_eqscale_seq_optimized/adult_entry.py`: pure birth/top-bin conversion, split queue, and final closed-renewal gate.
- `code/model/intergen_eqscale_seq_optimized/parameters.py`: default-off `adult_entry_clock='child_departure'`; opt in with `'split_birth_vintage'`.
- `code/model/intergen_eqscale_seq_optimized/solver.py`: active Markov fast/full KFE birth basis and residual; fails for a nonclosed/non-one-market or non-four-year mode. Household DP and child thinning stay unchanged.
- `code/model/tools/run_e5f_open_population_transition.py`: explicit `--adult-entry-timing split-16-20` on its existing birth-vintage path; default remains `legacy-20`. This driver retains its separate previously defined open-population closure.
- `code/model/intergen_eqscale_seq_optimized/tests/test_adult_entry.py` and `code/model/tools/test_run_e5f_open_population_adult_entry.py`: birth impulse, constant births/prehistory, queue mass, parent-death independence, top bin, final closed renewal, and unchanged old 20-year caller.

Torch should run the two focused test files using its existing source environment, then check the selected native fast/full Markov accounting and the post-normalization final gate. The new overnight driver must explicitly set the stationary clock override; a saved checkpoint with the old law is a seed, not a new adopted equilibrium. The open transition requires its explicit CLI switch to exercise split timing. No new model result is claimed here.
