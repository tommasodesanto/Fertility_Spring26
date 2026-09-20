# Persistent plus iid-transitory earnings candidate

This is a machine-readable diagnostic candidate only. It removes the measured
permanent component, keeps the age profile and pension external, and passes
household gross residual earnings through the flat payroll-tax model without
additional HSV progressivity compression. No solve, recalibration, cluster
run, or adoption decision is implied.

The source is the PSID reference-person earnings packet
`code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/`:
ages 25--60, positive IW and positive RP/spouse combined real gross labor
earnings, 1984--2019; age/year FE residualization; exact calendar-year
autocovariances at 13 lags; 199 person-cluster bootstrap and 10% covariance
correlation shrinkage. The nested no-fixed specification has 13 moments, 3
free parameters, objective `39.5085474236852`, versus `14.3178085110745` for
the fixed-effect specification.

The nested fit is reused directly. The diagnostic rho is recovered from its
lag-2/lag-1 fitted covariance ratio (`0.9703033301891563`), with persistent
variance `0.692826827789733` and transitory variance `0.3444421950282869`.
The optimized 15-state constructor is in
`code/model/tools/build_persistent_transitory_income_candidate.py`; its
transitory input is a period log SD. The coarse four-year diagnostic maps the
annual transitory variance to
`log(1 + (exp(V_trans_annual)-1)/4) = 0.09785297241466638`, preserving the
first two level moments of four independent annual mean-one shocks. This is a
coarse-period assumption, not an empirical estimate or exact lognormal
aggregation of the persistent path.

The full 13-row empirical covariance, bootstrap SE, nested fitted value, gap,
and pair count are in `earnings_comparison.csv`; coefficients and source
hashes are in `candidate.json`.
The constructor's numerical checks cover the stationary mean, row sums,
stationarity, log-income variance, and first two autocovariances. The iid input
is a log standard deviation; the three-node quadrature matches its log variance.
The failed first cohort job stopped on its baseline arm before candidate income
was evaluated. It provides no earnings-model evidence.
