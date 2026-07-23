# Wealth and bequest timing repair: advisor recap

## What was wrong

The July 22 calibration measured three wealth-related moments using a hybrid
balance sheet: beginning-of-period liquid wealth \(b_t\) was paired with the
household's newly chosen tenure \(\tau'_t\). This can count a newly purchased
house without subtracting its purchase price, or retain an owner's mortgage
debt after the house has been sold. The reported loss \(0.021954\) is therefore
not a valid calibration result.

For example, a renter with \(b_t=2\) who buys a house worth \(10\) has
post-transaction liquid wealth \(-8\) and total wealth \(2\), not \(12\).
Likewise, an owner with \(b_t=-8\) who sells that house for net proceeds \(9.4\)
enters the renter branch with \(1.4\), not the inherited debt \(-8\).

The household transaction map itself already applies purchases and sales. The
error was in the statistics layer that assembled the calibration moments from
the simulated distribution.

## Correct within-period accounting

The model timeline is:

\[
(b_t,\tau_t)
\longrightarrow \text{family/location/tenure choices and housing transaction}
\longrightarrow \widetilde b_t
\longrightarrow \text{consumption and } b'_t
\longrightarrow \text{survival or death}.
\]

Accordingly:

- Living-household cross-sectional wealth is
  \[
  W_t=b_t+pH_{\tau_t}.
  \]
- The PSID \(p90/p50\) target is also a cross-sectional wealth moment because
  its sample contains living reference persons aged 76--84. It is not a sample
  of decedents.
- If death occurs after the current choices, the transferred estate is
  \[
  E_t=\max\{b'_t+pH_{\tau'_t},0\}.
  \]
  The annual bequest flow weights this object by the age-specific death
  probability and divides the four-year flow by four.

The population code multiplies the next-age transition by survival before
propagating choices. Because survival depends only on age, that multiplication
commutes with the choice transition for survivors. It does not justify using
beginning-of-period wealth for the bequest flow, which must still be evaluated
at post-saving \(b'_t\).

## Verified repaired readout

At the previously reported parameter vector, two strict repaired solves are
bit-identical:

| Moment | Target | Repaired model |
|---|---:|---:|
| Wealth / annual after-tax earnings | 6.9000 | 6.0104 |
| Annual bequest flow / wealth | 0.0088 | 0.00571 |
| Living-old wealth \(p90/p50\), ages 76--84 | 3.4481 | 1.9112 |

The repaired total loss is \(0.349288\). The living-old dispersion and bequest
flow rows contribute \(0.198682\) and \(0.123677\), respectively.

All accounting/unit tests pass in both current-model packages. The repaired
14-by-14 local Jacobian solves all 29 perturbation cases and is numerically full
rank, but its condition number is \(14{,}224\). Identification is therefore
highly joint and weak; numerical full rank is not evidence of a clean
one-parameter/one-moment mapping.

## What can and cannot be claimed

We can defend the corrected measurement contract, the two repeated readouts,
and the conclusion that the old apparent wealth/bequest fit was spurious. We
cannot defend the old parameter estimates or loss. Older M4/M5 results that use
the same wealth-statistics layer require re-audit, although this does not by
itself invalidate their non-wealth moments.

A timing-repaired diagnostic recalibration is running on Torch (array
`14658852`, strict collector `14658853`). It should be interpreted only after
the collector produces two exact strict repeats and the complete target and
parameter tables are reviewed.
