# Simplified OLG theory: hostile audit and resolution

Date: 2026-08-31

This note records the final independent audit of the paper-facing simplified
OLG theory. The audited source snapshot is identified by Git blob hashes:

- section: `de710478cea65b072bc97b893e686a923cbbfba8`
- proof appendix: `338d7ccccbf5aa7a6ca79d964735f7728cb6bf7a`
- numerical driver: `65b71c9559d66da9df1d7f214ee1b3620db8101d`
- claim ledger: `9e20d8822ba78d0d7294e7877d463b247bdf2489`

## Final verdict

The household algebra, old-age active sets, mortgage cancellation, conditional
young policies, mixed-tenure aggregation, steady-state equations, conditional
steady-state existence theorem, intergenerational marginal-value gap,
compensated reallocation result, and fertility derivative pass. Independent
numerical checks reproduce household first-order conditions to below
`3.4e-16`, steady-state residuals below `1.7e-12`, and the reported maximum
scaled transition-truncation residual of `2.07e-10`.

The dynamic result is deliberately split into three claims:

1. The paper defines an exact infinite-horizon transition equilibrium as an
   equilibrium sequence converging between two positive steady states.
2. Directional implicit-function and Poincare--Miranda arguments establish
   roots of terminally closed finite-horizon truncations under their stated
   sufficient conditions.
3. A limit of those roots is an exact infinite-horizon transition only under
   the compactness, continuity, fixed-date convergence, and uniform-tail
   conditions in Proposition `prop:olg_app_transition_limit`.

The computed horizon-28 path is therefore labeled a high-accuracy transition
approximation, not an exact infinite-horizon equilibrium. Its small terminal
state gap and horizon-24 versus horizon-28 stability are supporting diagnostics;
they do not prove the uniform-tail condition.

## Issues found and resolved

- The original finite-horizon language imposed final steady-state prices,
  transfers, and continuation values but did not impose the full generated
  terminal state. Every affected theorem, caption, computation description,
  README entry, and claim-ledger row now says “terminally closed truncation” or
  “transition approximation.”
- The main young problems now state `a' >= 0`, matching their KKT systems.
- The numerical wedge is now called the gross marginal-value gap. The text
  states separately that a compensated improvement requires this gap to exceed
  the real transfer cost `K`.
- The verification output now reports owner-size slack. Its minimum is
  `1.204597`, establishing that the displayed owner access wedge is purely the
  financing wedge in the numerical construction.
- The frozen-branch Jacobian is described only as a finite-root conditioning
  diagnostic, not evidence for the directional theorem at a zero-size shock.
- Distribution-measure punctuation and conflicting uses of `sigma` and `rho`
  in the simplified theory were removed.

## Scope boundary

The integrated paper draft contains a pre-existing inconsistency between two
sets of quantitative counterfactual numbers (`2.1/2.7` percent versus
`3.4/4.9` percent birth effects, with corresponding price changes). This lies
outside the simplified-theory audit and has not been silently resolved. The
live quantitative model and counterfactual contract remain governed by
`CALIBRATION_STATUS.md`.
