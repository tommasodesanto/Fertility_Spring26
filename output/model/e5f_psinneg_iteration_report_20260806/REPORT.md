# Nonnegative-child-utility iteration: calibration and policy report

## Executive conclusion

This iteration asks whether the repaired sequential-fertility model can replace
the ad hoc child expenditure-share tilts with a more transparent restriction:
each child currently at home requires a minimum quantity of housing services.
The answer is mixed. The child-room floor is economically easier to explain and
fits the parent housing gradient better, but it does not repair the central
fertility-timing or wealth-distribution failures. Under the funded policy
experiment, it produces nearly the same equilibrium responses as the tilt
model: completed fertility rises by 0.15--0.67 percent, while total births rise
by 9.1--10.8 percent because city population expands. The baseline entry level
now matches the approved 16.9 percent outside-origin entrant share exactly.
The policy accounting is correct; the entry elasticity is still externally
fixed and therefore remains outstanding.

Neither model is promoted. The pre-E6 maintained model remains unchanged.

## What changed economically

Parity, (n), remains children ever born. The existing child state, (m), is
the number currently at home. A birth maps ((n,m)) to ((n+1,m+1)). Each
at-home child independently leaves in a four-year period with probability
(2/9), so time at home is 18 years in expectation. No state dimension was
added.

The floor model replaces the two child expenditure-share tilts with

\[
  h_t \geq m_t\,\bar h,
\]

where the estimate is \(\bar h=0.6136\) housing-service units per child at
home. Thus the minimum is 0.614, 1.227, or 1.841 units with one, two, or three
children at home. At the estimate, all discrete owner rungs remain feasible;
the restriction mainly disciplines low-service renter choices. Direct child
flow utility is restricted to \(\psi\geq0\); both the floor and tilt estimates
put it at, or effectively at, zero.

## Calibration verification

Both models use the same 12 targets and canonical weights. The doubled search
contains 16 chains per model. All 32 chains are eligible, and the certified
winners have strict exact repeats with zero loss and moment differences. An
independent solve reproduces each certified moment within (2.2\times10^{-15}).

- Child-room floor: loss 297.1112; market residual (1.05\times10^{-5}).
- Expenditure-share tilts: loss 299.2700; market residual (2.41\times10^{-5}).

The 0.72 percent loss difference is not the substantive reason to prefer the
floor. Its substantive advantage is interpretability; its substantive cost is
excess aggregate housing demand.

## Complete target fit

| Moment | Target | Floor model | Floor gap | Weight | Floor contribution | Tilt model |
|---|---:|---:|---:|---:|---:|---:|
| Completed fertility | 1.918 | 1.96667 | +2.54% | 1,425.74 | 3.377 | 1.97565 |
| Childlessness | 0.188 | 0.103896 | -44.74% | 17,180.7 | 121.529 | 0.107981 |
| Mean age at first birth | 25.3106 | 25.5707 | +1.03% | 16 | 1.083 | 25.7014 |
| Share of first births at 30+ | 0.270062 | 0.327463 | +21.25% | 15,625 | 51.482 | 0.335458 |
| Housing increment, zero to one child | 0.664435 | 0.690190 | +3.88% | 906.056 | 0.601 | 0.687920 |
| Rooms: 3+ minus 1--2 children, ages 30--55 | 0.367700 | 0.359137 | -2.33% | 2,958.51 | 0.217 | 0.341796 |
| Parent-minus-nonparent ownership gap | 0.167662 | 0.161950 | -3.41% | 14,229.6 | 0.464 | 0.172796 |
| Prime-age ownership rate | 0.575472 | 0.582480 | +1.22% | 1,207.85 | 0.059 | 0.598047 |
| Mean occupied rooms, ages 18--85 | 5.77997 | 6.47800 | +12.08% | 11.9732 | 5.834 | 5.94846 |
| Wealth / annual gross labor earnings | 6.8731 | 5.79412 | -15.70% | 6.28767 | 7.320 | 5.78765 |
| Annual bequest flow / wealth | 0.0088 | 0.00871298 | -0.99% | 5,165,289 | 0.039 | 0.00885423 |
| Old-age wealth p90/p50 | 3.44811 | 2.08970 | -39.40% | 56.9598 | 105.106 | 2.09674 |

The floor matches the two family-size housing gradients well. The binding
failures are childlessness, the share of late first births, and upper-tail
wealth dispersion. The supplemental age distribution confirms too many first
births at 18--25 and 38--45 and too few at 26--33.

## Parameters and bounds

### Child-room floor model: 9 estimated parameters

| Parameter | Economic role | Estimate | Bound | Near bound? |
|---|---|---:|---:|---|
| \(\beta\) | Annual discount factor | 0.994634 | [0.94, 0.9995] | No |
| \(\psi\) | Direct flow value per child at home | 0 | [0, 3] | Lower |
| \(\kappa_0\) | First-birth choice-shock scale | 2.98684 | [0.02, 50] | No |
| \(\kappa_+\) | Later-birth choice-shock scale | 1.71332 | [0.02, 50] | No |
| \(\chi\) | Owner housing-service premium | 1.03694 | [0.1, 5] | No |
| \(H_0\) | Housing-supply scale | 9.87303 | [0.2, 80] | No |
| \(\theta_0\) | Bequest-motive strength | 0.159351 | [0, 8] | Lower 2% |
| \(\theta_1\) | Bequest curvature/shift | 0.093730 | [0.02, 16] | Lower 2% |
| \(\bar h\) | Room floor per child at home | 0.613605 | [0.1, 1.8] | No |

The expenditure-share control instead estimates
\(\Delta\alpha=0.016919\in[0,0.25]\) and
\(\Delta\alpha^{jump}=0.020344\in[0,0.25]\), and omits \(\bar h\). Its other
estimates are \(\beta=0.993728\), \(\psi=0.000258\),
\(\kappa_0=3.56626\), \(\kappa_+=1.77449\), \(\chi=1.04204\),
\(H_0=7.43944\), \(\theta_0=0.164510\), and \(\theta_1=0.086022\).
The complete machine-readable parameter table contains every bound and flag.

## Correctly funded policy experiment

The policy baseline is a 1 percent annual property tax rebated lump sum. The
counterfactuals are (i) a 2 percent tax with the revenue rebated, and (ii) the
same 2 percent tax plus a 0.4 purchase grant for homes with at least six rooms.
Every case balances the government budget.

At the baseline, the model entry probability is not imposed arbitrarily. The
approved empirical normalization is that 16.9 percent of entrants originate
outside the city. For each calibrated candidate, the driver computes

\[
 q^*=\frac{1-s^{E,\mathrm{out}}}{B_0/E_0},
 \qquad s^{E,\mathrm{out}}=0.169,
\]

and then recovers the outside value and outside potential-entrant flow so that
stationary population has scale one. This gives $q^*=0.969225$ for the floor
and $q^*=0.969844$ for the tilt control. The two recovered outside objects and
the entry taste scale 2 are held fixed. For each counterfactual, the house price
and transfer clear jointly, while entry determines population scale:

\[
 q=\sum_z\omega_z\Lambda\!\left(2[V_z-\bar V^{out}]\right),
 \qquad
 S E_0=q\left[M+S B_0\right].
\]

Here $E_0$ is required entrant flow per unit population, $B_0$ is mature
city-born entrant flow, and $M$ is the fixed outside potential-entrant flow.
The baseline identity recovers $S=1$ and the 16.9 percent share to numerical
precision in both models.

### General-equilibrium policy effects

All changes are relative to each model's funded 1 percent baseline.

| Model | Policy | TFR | House price | Population | Total births | Entry probability |
|---|---|---:|---:|---:|---:|---:|
| Floor | 2% tax | +0.294% | -15.050% | +9.980% | +10.264% | 0.9821 |
| Floor | 2% tax + grant | +0.674% | -15.196% | +10.104% | +10.771% | 0.9794 |
| Tilts | 2% tax | +0.146% | -15.305% | +8.988% | +9.128% | 0.9825 |
| Tilts | 2% tax + grant | +0.432% | -15.075% | +9.062% | +9.487% | 0.9804 |

For the floor model, TFR levels are 1.98015 in the baseline, 1.98598 under the
2 percent tax, and 1.99351 with the grant. The corresponding house prices are
0.77609, 0.65929, and 0.65815. The balanced period transfers are 0.20913,
0.35455, and 0.30096.

### Fixed-population decomposition

Holding population fixed, the floor model gives TFR changes of +0.395 percent
under the 2 percent tax and +0.788 percent with the grant; house prices fall
18.08 and 18.40 percent. Endogenous entry raises housing demand, attenuating
the price decline and therefore attenuating the fertility response. These
fixed-population rows are decompositions, not counterfactual equilibria.

### Accounting audit

All six full counterfactuals pass. Across the two models, the maximum absolute
fixed-population housing-market residual is $2.23\times10^{-5}$, the maximum
fixed-population budget residual is $2.38\times10^{-5}$, the maximum joint
housing-market residual is $2.35\times10^{-5}$, and the maximum joint budget
residual is $4.89\times10^{-6}$. The accounting identity reproduces the 16.9
percent outside-origin share within $5.4\times10^{-13}$ and baseline scale
one to machine precision; the numerical joint root is within $2.1\times10^{-5}$.

## Interpretation and decision

The room floor succeeds as an economically legible replacement for the two
share tilts, and it strengthens the grant's fertility response slightly. It
does not solve the model's low childlessness, late-fertility, or wealth-tail
problems. The corrected entry normalization cuts the previously reported
19--22 percent total-birth response roughly in half, to 9--11 percent. Most of
that response is still population entry rather than fertility per household.

### Outstanding policy objects

- The baseline entry **level is closed** by the 16.9 percent empirical share.
- The entry taste scale remains fixed at 2; its empirical discipline is
  outstanding. No sensitivity grid is run in this iteration.
- The local-born retention weight remains fixed at 1; its empirical discipline
  is outstanding.
- The distinction between a local reallocation experiment and a national
  fertility experiment remains outstanding.

The older $q=0.5$ entry-adjusted policy rows are withdrawn. Their
fixed-population decomposition rows are unchanged and remain valid.

## Figures and complete files

- `target_relative_gap_comparison.png`: all 12 signed proportional gaps.
- `funded_policy_comparison.png`: the four principal policy responses.
- `target_fit_comparison_full.csv`: full target, model, gap, weight, and loss data.
- `parameter_comparison_full.csv`: every estimate, bound, and near-bound flag.
- `policy_comparison_full.csv`: every level, change, entry, price, fiscal, and market result.
- Floor diagnostics: `../intergen_e5f_child_room_floor_psinneg_extended_20260806/diagnostic_packet/`.
- Tilt diagnostics: `../eqscale_seq_e5_maturation_repair_psinneg_extended_20260806/diagnostic_packet/`.
- Floor policy details: `../intergen_e5f_child_room_floor_psinneg_policy_empirical_entry_20260806/`.
- Tilt policy details: `../eqscale_seq_e5_maturation_repair_psinneg_policy_empirical_entry_20260806/`.
