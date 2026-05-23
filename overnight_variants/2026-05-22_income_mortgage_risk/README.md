# Income Risk / Mortgage Account Prototype

This folder is an isolated overnight branch test. It contains a copied
`dt_cp_model/` package and does not modify the live implementation under
`code/model/dt_cp_model`.

The smoke driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_smoke.py --quiet
```

The second-pass scenario driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v2_scenarios.py --quiet --nb 40
```

The third-pass HANK earnings-risk driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v3_hank_z.py --quiet --nb 30 --nz 3
```

The full-equilibrium HANK earnings-risk driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v4_hank_z_ge.py --quiet --nb 30 --nz 7 --max-iter-eq 35
```

The outside-option closure HANK earnings-risk driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v5_hank_z_outside_closure.py --quiet --nb 30 --nz 7 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --baseline-max-iter-eq 35 --max-iter-eq 60 --normalization-passes 2 --scale-target-tol 1e-5 --tol-eq 5e-4
```

The V5 `Nz=5` exploratory equilibrium and figure driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python plot_income_mortgage_risk_v5_hank_z_outside_closure_nz5.py --quiet --nb 30 --nz 5 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --max-iter-eq 60 --tol-eq 5e-4
```

It writes separate `Nz=5` results, diagnostics, report, log, and figures so the
accepted `Nz=7` validation files are not overwritten.

The V5 `Nz=5` directional calibration audit is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_v5_nz5_directional_calibration.py --quiet --nb 30 --nz 5 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --max-iter-eq 35 --tol-eq 5e-4
```

It writes small one-at-a-time and targeted-combination calibration probes
without overwriting the accepted V5 files.

The V5 `Nz=5` parallel global search driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_v5_nz5_global_search.py --run-tag v5_nz5_global_test --results-dir global_search_v5_nz5 --workers 8 --budget-sec 28800 --max-evals 1000000 --seed 20260523 --nb 30 --nz 5 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --max-iter-eq 35 --tol-eq 5e-4 --global-prob 0.75
```

For an overnight laptop run, use the `nohup caffeinate` command in
`REPORT_V5_NZ5_GLOBAL_SEARCH_PLAN.md` so the process survives the terminal and
prevents sleep.

The full-equilibrium figure driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python plot_income_mortgage_risk_v4_hank_z_ge.py --quiet --nb 30 --nz 3 --max-iter-eq 35
```

It writes the equilibrium and sorting figure set to
`figures_v4_hank_z_ge/`.

The borrowing-wedge validation diagnostic is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_hank_z_borrowing_wedge_diagnostics.py --quiet --nb 30 --nz 3 --max-iter-eq 35
```

## Implemented Objects

- Earnings grid \(z\) as a true discrete state. The current HANK pass uses a
  Rouwenhorst approximation with configurable `Nz`, persistence `rho_z`, and
  unconditional log-earnings dispersion `sigma_z`.
- Mortgage-account grid \(\mu\) with five encoded states:
  no mortgage/good credit, no mortgage/bad credit, active low/good, active
  high/good, and active low/bad.
- Diagnostic mortgage objects: amortization rate, good/bad payment rates,
  good/bad financed shares, default penalty, and bad-credit recovery.
- A copied-solver helper that attaches an augmented diagnostic distribution
  over \((z,\mu)\) after one benchmark equilibrium solve.

## Important Limitation

The first overnight implementation did not put \(z\) and \(\mu\) inside the
Bellman recursion. V3 and later correct this for \(z\): the HANK earnings
state is now structural in values, policies, and the forward distribution.
The mortgage-account state \(\mu\) remains diagnostic only and should not be
treated as an accepted structural default/mortgage implementation.

The generated `REPORT.md` states this explicitly and gives the branch verdict.

## Second-Pass V2

`run_income_mortgage_risk_v2_scenarios.py` improves on V1 by making the new
objects affect solved choices in a finite set of partial-equilibrium scenarios.
The script solves 15 economies: 3 earnings states times 5 mortgage-account
states. Earnings scale \(y_a(z)\), bad credit lowers the financed share
\(\phi_g\), and active mortgage accounts add a compact owner-payment wedge
\(\rho_g m\).

This still is not a full structural implementation because there is no single
Bellman recursion over \((b,d,i,a,n,s,z,\mu)\). See `REPORT_V2.md`,
`results_income_mortgage_risk_v2.csv`, and
`diagnostics_income_mortgage_risk_v2.csv`.

## Third-Pass V3 HANK-z

`run_income_mortgage_risk_v3_hank_z.py` implements the standard HANK-style
idiosyncratic earnings process as a true state. The copied solver carries
value functions, policies, fertility probabilities, location probabilities,
tenure choices, and the forward distribution over \((b,d,i,a,n,s,z)\). The
earnings state follows the finite Markov chain reported in
`diagnostics_income_mortgage_risk_v3_hank_z.csv`, and continuation values
average over \(\Pi_z\) after child aging.

This is the correction to the earlier Branch 1 mistake: \(z\) is no longer a
post-solve label or a scenario loop. The output remains a partial-equilibrium
prototype at copied benchmark prices, and \(\mu\) is not yet structural. Adding
\(\mu\) requires tenure-dependent account transitions for purchase,
amortization, sale, and default. See `REPORT_V3_HANK_Z.md`,
`results_income_mortgage_risk_v3_hank_z.csv`, and
`diagnostics_income_mortgage_risk_v3_hank_z.csv`.

## Fourth-Pass V4 HANK-z Full Equilibrium

`run_income_mortgage_risk_v4_hank_z_ge.py` runs the same structural \(z\)-state
household problem inside the copied equilibrium price/entry loop. This
supersedes the V3 fixed-price result as the relevant Branch 1 smoke evidence.
On the coarse `Nb=30`, `Nz=3` run, it accepted at strict tolerance in 15
iterations and runs in the expensive-but-workable range after the 2026-05-22
compiled HANK-\(z\) forward-pass fix. The solver still reports the qualitative
cost category as expensive because this is a new state dimension and the final
full-statistics pass is intentionally retained. See `REPORT_V4_HANK_Z_GE.md`,
`results_income_mortgage_risk_v4_hank_z_ge.csv`,
`diagnostics_income_mortgage_risk_v4_hank_z_ge.csv`, and
`diagnostics_income_mortgage_risk_v4_hank_z_ge_trace.csv`.

## Borrowing-Wedge Diagnostic

`run_hank_z_borrowing_wedge_diagnostics.py` reruns the V4 HANK-\(z\) GE
prototype and asks whether the existing \(\phi\)-based down-payment and
borrowing-floor wedge is active. It does not add a mortgage/default account.

The `Nb=30`, `Nz=3` diagnostic accepted at the same strict GE tolerance as V4.
High-\(z\) renters have a renter-to-owner purchase share of 0.091 versus 0.015
for low-\(z\) renters, and starter down-payment feasibility is 0.590 versus
0.336. Low-\(z\) owners have more mass close to the borrowing floor, 0.089
versus 0.029 for high-\(z\). This validates that the existing collateral and
liquidity wedge is being used. It is not a default or refinancing diagnostic.

See `REPORT_HANK_Z_BORROWING_WEDGE.md`,
`results_hank_z_borrowing_wedge.csv`,
`diagnostics_hank_z_borrowing_wedge.csv`, and
`hank_z_borrowing_wedge.log`.

## Fifth-Pass V5 HANK-z Outside-Option Closure

`run_income_mortgage_risk_v5_hank_z_outside_closure.py` switches the HANK-\(z\)
GE loop to the paper-facing outside-option scale closure. The benchmark now
imposes the \(S=1\) normalization directly inside the copied GE loop. At each
candidate composition it computes \(E_0(p)\) and \(B_0(p)\), calibrates
\(\bar W^E\) to the target \(q^E\) at current entry values, and sets
\[
M=E_0(p)/q^E(p)-B_0(p).
\]
The reported final \(M\) and outside value are the objects to hold fixed in
counterfactuals, where the fixed-\(M\) scale equation is
\[
S E_0(p)=q^E(p)\left[M+S B_0(p)\right].
\]

The accepted `Nb=30`, `Nz=7` run uses \(\rho_z=0.95\), unconditional
\(\sigma_z=0.35\), and a separate entry scale \(\kappa_E=10^6\). A separate
\(\kappa_E\) is required because the outside option compares lifetime-utility
levels; reusing the incumbent location scale \(\kappa_\ell\) made the outside
probability jump to zero or one and destabilized the price loop. The accepted
run no longer uses an outer normalization pass; it converged at strict
tolerance with final scale \(S=1\), \(q^E=0.9\), outside probability `0.1`,
and residual outside-born flow \(M=0.00578025\). It remains yellow because the
un-recalibrated moments are economically poor. See
`REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md`,
`results_income_mortgage_risk_v5_hank_z_outside_closure.csv`,
`diagnostics_income_mortgage_risk_v5_hank_z_outside_closure.csv`,
`diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_closure.csv`, and
`diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_trace.csv`.

## V5 Speed Audit

`run_v5_speed_audit.py` measures the benchmark-normalized V5 closure without
overwriting the accepted V5 result files. The fixed-iteration audit writes
`speed_audit_v5_benchmark_normalized.csv`,
`speed_audit_v5_benchmark_normalized.json`, and
`REPORT_V5_SPEED_AUDIT.md`. A separate same-process warm accepted reference
writes `speed_audit_v5_warm_accepted.csv/json`.

The audit shows that the cold `Nb=30`, `Nz=7` accepted solve took `620.13`
seconds, while the same accepted solve after a same-process warm-up took
`200.82` seconds. The warm solve is still too slow for calibration: the
remaining cost is split across the full HANK-\(z\) Bellman loop, HANK-\(z\)
forward distribution, and the final full-statistics pass.

## V5 Nz=5 Equilibrium Figures

`plot_income_mortgage_risk_v5_hank_z_outside_closure_nz5.py` reruns the V5
benchmark-normalized outside-option closure at `Nb=30`, `Nz=5`, \(\rho_z=0.95\),
unconditional \(\sigma_z=0.35\), and \(\kappa_E=10^6\). It accepted at strict
GE tolerance in 12 iterations with final GE error `3.46432e-4`. The full solve
and plot packet build took `83.10` seconds. The final benchmark normalization
objects are \(S=1\), \(q^E=0.9\), outside probability `0.1`, residual
outside-born flow \(M=0.00580643\), and outside value `-1504137.922586`.

The figure packet is
`figures_v5_hank_z_outside_closure_nz5/HANK_Z_OUTSIDE_CLOSURE_NZ5_FIGURE_PACKET.pdf`.
The run-specific report is `REPORT_V5_HANK_Z_OUTSIDE_CLOSURE_NZ5.md`; moment
values are in `results_income_mortgage_risk_v5_hank_z_outside_closure_nz5.csv`.

## V5 Nz=5 Directional Calibration

`run_v5_nz5_directional_calibration.py` probes whether the bad V5 HANK-\(z\)
moments move under small parameter changes before running a full optimizer.
The best probe so far is `alpha_high_fertility_mid_finance`, which sets
`alpha_cons=0.80`, `kappa_fert=4.0`, and `phi=0.90`. It accepted strict GE in
16 iterations, took `107.73` seconds, and cut the `Nz=5` loss from `311.30` to
`166.72`, improving 16 of 19 target gaps.

This is promising but not calibrated: TFR moves near target, ownership becomes
positive-gradient, room moments improve, and geography improves modestly, but
first-birth timing is still far too late, childlessness is still too high, the
fertility gradient remains wrong-signed, and young liquid wealth worsens. See
`REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md`,
`v5_nz5_directional_calibration.csv`, and
`v5_nz5_directional_calibration_round2.csv`.

## V5 Nz=5 Parallel Global Search

`run_v5_nz5_global_search.py` runs a true global randomized search over 19
parameters, including the original structural calibration parameters,
geography shifters, `alpha_cons`, `phi`, `hR_max`, and `h_own_max`. The
objective is the full 19-moment weighted SMM loss plus GE/scale penalties. It
parallelizes candidates across worker processes and writes every completed
evaluation to `evaluations.jsonl` and `evaluations.csv`, with continuously
updated `best.json`, `best_summary.md`, `latest.json`, `status.json`, and
`heartbeat.txt`.

The intended 8-core laptop command uses 8 workers with one numerical thread per
worker for 8 hours. See `REPORT_V5_NZ5_GLOBAL_SEARCH_PLAN.md`.
