All decisive claims are now verified in source. Writing the assessment.

**Recommendation.** Keep the birth-vintage queue as the adult-entry clock; do not couple entry to the geometric departure law. Read \(m\) as economic dependency, not co-residence. Change only forward accounting, labels and guards. This is my proposal, not an authorized change. Startup limit: no `memory/daily` file exists, so I used `memory/AGENT_MEMORY.md` and `CALIBRATION_STATUS.md`. Code claims were checked in source; the 42/66 age support and flow magnitudes come from the archived September 12 receipt, not a reopened selected checkpoint.

## What is inconsistent, and what is not

**Verified in source.** Each dependent leaves independently with \(\mu=1/4.5=2/9\) per period, newborn included in the first draw (`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/parameters.py:1229`). Both dated drivers apply survival first and pass survivors only to the child kernel, so dependents of exiting households are deleted with no entry (`code/model/tools/run_e5f_open_population_transition.py:867`). Entry pops the oldest of four slots, appends \(B^w_t/2.1\), and inserts next date, so a birth at \(t\) enters at \(t+20\) years; the queue already carries the top-bin weight (same file, lines 351 and 1219). The stationary solver's own counter uses conversion 0.5 and ordinary departures only (`solver.py:5653`). It is a third object and must remain diagnostic.

**Question 1.** Nobody is lost. Two ledgers each conserve: dependents \(C_{t+1}=C_t+B_t-M_t-R_t\) and the pre-entry register \(N_{t+1}=N_t+B^w_t/2.1-A_t\). What is conflated is identity: leaving \(m\) and entering are different events for the same birth. In the stationary state this creates a gap stock, departed but not entered and using no resources, and an overlap stock, still in \(m\) and already a household. With \(B=0.115\) per household-period from the receipt:

| Object | Child units per household |
|---|---:|
| Dependents in \(m\) | 0.498 |
| Overlap, aged 20 or more | about 0.13 |
| Gap, departed before 20 | about 0.21 |

These are composition errors of the geometric law, not leaks. They are a permissible approximation under the dependency reading and a defect only under a literal co-residence reading. In both readings the defect sits in the \(m\) law the author accepted, not in the queue.

**Question 2.** Parental death adds nothing under the reviewed support. The last birth is at parent age 42; that vintage pops when the parent is 58 and enters at 62; mortality first applies at 66. Every dependent released at death has already entered, so a deduplicated death-entry stream is identically zero and deleting those dependents is correct. The receipt's \(R=0.006\) per period is overlap being resolved, not orphans. The queue is also robust to younger mortality: the vintage enters on schedule whatever happens to the parent; only care during the gap would be unmodeled. Steady state does not distinguish the rules, since \(F=B\) gives \(E=B^w/2.1\) either way with the same household age distribution. The choice is purely about transition timing. The existing arithmetic establishes that coupling sends 22% of a birth cohort into adult households after 4 years and 63% before 18, with cumulative extra entries by year 20 of 2.5 times the queue's after a permanent birth increase. Price, fertility and welfare paths under either rule are unmeasured.

## Which rule, and why

**Question 3.** The entry clock carries fertility into housing demand and into the next generation's fertility, so it should approximate the observed household-formation lag. Real household formation concentrates at ages 18 to 30 and is essentially absent before 16. The queue matches that; a clock with median 12 years and standard deviation 15.9 years does not. Coupling would front-load every population and housing response to a fertility shock by more than a decade and shorten the generation length by construction. For a paper whose population headlines are already small, near 1 to 3%, that is a distortion, not an approximation.

The coordinator oscillated because each turn optimized a different criterion: one-clock bookkeeping, then the timing distribution, then "lifecycle consistency" heard as "same person, same clock". Lifecycle consistency is chronological: born at \(t\), adult near \(t+20\). The queue has it; the geometric \(m\) law is the inconsistent object. Moving entry onto \(m\)'s clock satisfies the letter and violates the substance.

Worked timeline for a birth in 2000 to a parent aged 30:

| Date | Queue, as implemented | Coupled |
|---|---|---|
| 2004 | 22% leave \(m\); no entry | 22% become age-18 households |
| 2012 | 53% have left; no entry | 53% are households |
| 2020 | vintage enters, \(1/2.1\) per birth; 28% still in \(m\) | 72% entered |
| 2036 onward, parent 66+ | remaining dependents deleted at death, no entry | same, already entered |

## Equations, alternatives, next steps

**Question 4, proposed.** Let \(T_t\) be the flow entering the 3+ group, \(\delta=0.602\), \(B^w_t=B_t+\delta T_t\); \(M_t\) ordinary departures from survivors; \(R_t\) death release; \(X_t\) exits.
\[
C_{t+1}=C_t+B_t-M_t-R_t,\qquad
N_{t+1}=N_t+\tfrac{1}{2.1}B^w_t-A_t,\qquad
E_{t+1}=\rho_tA_t+I_t,
\]
with \(A_t\) the vintage born five periods earlier and no \(R_t\) term. Entrants start at age index 0 with \(n=m=0\) and the existing wealth and earnings draws. Estates are generated at death and paid through the existing 45 to 65 pool; entrants receive nothing new. Terminal tests are unchanged. The factor \(1/2.1\) is a replacement normalization: \(Q_t=2.1H_t+C_t\) is a replacement-unit index, and literal residents \(2H_t+C_t\) are not claimed. \(C\) is literal with 3+ counted as three; the register is weighted. Label both.

**Question 5.** Three alternatives, all within the current household state:

1. **Coupling.** Forward-only, but it needs a tagged third-child moment to preserve \(B^w\), reconstruction of that moment in initial states, a new terminal criterion, and reruns of every transition. Steady state unchanged; transitions distorted as above. Roughly 1 to 2 engineering days plus reruns.
2. **Parent-age hazard for \(m\)**, already coded and default off. Ends dependency by parent age 62, shrinking gap and overlap without new states, but it changes the household's child-cost profile, so it alters values and requires recalibration. Not today's step.
3. **Spread queue**, equal thirds at 16, 20 and 24 years. Forward-only and cheap; no evidence favors it now.

Only option 2 touches household information or values; options 1 and 3 change entrant resources over time.

**Question 6.** Cost: hours, no solve. Record the ledgers above in `CALIBRATION_STATUS.md` and the theory note. Extend the read-only measurement script to report \(B,B^w,M,R,X,A,E\) and the gap and overlap stocks per date. Add a guard so neither a death-entry stream nor the 0.5-conversion counter can reach production entry. On the selected parameter object, recheck last fertile age, first mortality age and queue length, which I verified only on the archived receipt. Optionally run a fixed-policy entry-path comparison from the saved birth path to record the size of the rejected change.

Not established: literal person conservation; co-residence realism of \(m\); dated wiring of the estate pool if active; behavior under younger mortality. Evidence that would change the recommendation: a household state making \(m\) chronological, currently ruled out; observed 2007 to 2023 head-age marginals showing the 20-year lag is too long, which would argue for a shorter or spread queue rather than coupling; or a mechanism in which co-resident children's housing needs drive results, which would argue for option 2.
