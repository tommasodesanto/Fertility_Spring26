# Lead resolution of independent design review

The reviewer read the initial specification before the new source files appeared.
It is a mathematical/specification review, not an independent code certification.

1. Location: the reference has I=1. The experiment rejects multiple locations and includes the original amenity minus staying-cost constant in each conditional tenure value. No omitted spatial alternative remains.
2. Control: the spec now gives both sequential inclusive-value formulas and calls this a fertility-first sequential-logit control. It is not a valid reversed GEV when lambda<1 and is never described as one.
3. Expected maximum: the spec now explicitly labels V as the mean-zero-shock expected maximum. The displayed uncentered law differs by the common Euler-constant times outer scale. Choices are unchanged; no welfare comparison is reported.
4. Feasibility: plan_values branches on conception endpoints and maps the original dead sentinel to negative infinity before mixing outcomes. Unit tests cover zero, one and mixed conception, missing outcomes and first-birth costs.

Lead source review checks capture call order, the existing deterministic housing maximum and transaction maps, current-family to successful-birth mapping, and direct correlated scattering. Runtime acceptance additionally requires two exact reference replays, probability/mass/birth identities and newly occupied household budgets. Production adoption and identification remain open.
