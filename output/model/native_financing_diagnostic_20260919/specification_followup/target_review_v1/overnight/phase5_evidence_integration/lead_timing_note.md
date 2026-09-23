# Conditional timing algebra; no adopted observer change

The native operator measures the treated/control difference after one full model transition. Its origin branches use identical pre-choice states, assign a first birth in one branch, and hold the other childless. At destination it allows continuation births only in the treated branch. It neither samples a calendar birth date within the four-year cell nor runs the data regression.

For illustration only, add piecewise-constant housing and a calendar birth offset $U\in[0,4)$. Let $D_0$ and $D_1$ be treatment effects at origin and destination; impose zero effect before the origin. If $a=\Pr(U<1)$ and the timing law is independent of effect heterogeneity, the event-$-1$ to $+3$ effect is

\[
aD_0+(1-a)(D_1-D_0)=(1-a)D_1+(2a-1)D_0.
\]

Uniform timing gives $.75D_1-.5D_0$, not $D_1$ in general. But Opus's claim that equality is impossible with a nonzero $D_0$ is false: $U=0$ and $D_0=D_1=1$ yields equality. Equality at special parameter values cannot establish a generally correct observer.

Under the same additional assumptions, the symmetric $-2/+2$ window gives $.5D_1$, and $-4/+4$ gives $D_1$. These identities do not prove that either empirical cohort regression identifies the model intervention, that within-period housing is constant, or that birth timing can be shifted independently of housing/fertility decisions. A longer prebirth window may capture anticipation, partnership and selection absent from the matched intervention. Cohort support, empirical weights, household selection and subsequent births must be addressed separately.

Therefore compute alternative empirical contrasts as diagnostics, but do not silently rescale the model or replace the target. Specify the within-period observation rule and its economic interpretation first. The fixed running comparison remains unchanged.
