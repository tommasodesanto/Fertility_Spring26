# Experimental CES with separate family scales

This private snapshot changes stationary material preferences and estimates ten coordinates. It is not adopted, and the floor-search seed is a provisional initializer, not a verified CES calibration. Objective and reported policies use the same integrated scalar/indexed/full kernels.

For nonhousing consumption $c$, physical renter rooms $h$ and owner services $\chi h$, write $s=h$ or $\chi h$, respectively. Current children at home are $m$, $e_c=((2+.7m)/2)^{.7}$ and $e_h=e_c(1+\lambda_h 1\{m>0\})$.

$$X=[\alpha(c/e_c)^\rho+(1-\alpha)(s/e_h)^\rho]^{1/\rho},\quad \rho=(\eta-1)/\eta,$$
$$u=X^{1-\sigma}/(1-\sigma)+\psi_{child}m^{1-\kappa}.$$

The inherited direct-child-benefit convention and additive constants remain unchanged. Fixed $\sigma=2$, $\alpha=.733$ and $\eta=.487$; $\alpha$ is a primitive CES weight, so it does not enforce a 26.7% housing expenditure share. No housing or consumption floor, compensation factor, or varying CES weight is used. CES uses inherited model units with no reference-rent compensation. Physical product choice, housing costs, financing, entrant wealth, timing and population/renewal closure remain inherited.

The renter interior satisfies $h/c=[(1-\alpha)\alpha^{-1}(e_c/e_h)^\rho/r]^{\eta}$ and $c+rh=S$. If this allocation exceeds the inherited rental cap, $h=h_{max}$ and $c=S-rh_{max}$. CES expenditure price is $[\alpha^\eta e_c^{1-\eta}+(1-\alpha)^\eta(re_h)^{1-\eta}]^{1/(1-\eta)}$. Owner utility fixes services at $\chi h$.

Saving remains exhaustive over every wealth-grid continuation interval and the renter cap kink. Each interval has linear continuation. Concavity makes $u_c=\beta V'$ the unique interior optimum if bracketed; 44 deterministic bisections find it. All interval endpoints remain evaluated. Existing credit, solvency and timing blocks are retained. The original Cobb-Douglas branch is preserved when `ces_eta=0`. Guarded CES operation requires the native indexed/full kernel and rejects floor/multiplier inputs.

Sources: original reference `code/model/refactor_lab`; initializer `utility_floor_psi_v1/chain_11/results/0026_nm/proposed_parameters.json`; fixed elasticity from Li, Liu, Yang and Yao, *Housing Over Time and Over the Life Cycle*, Table 2. The search packet at `output/model/fixed_reference_economics_20260928/utility_ces_scales_v1` owns native import pins and reporting metadata.

Focused verification: `NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/anaconda3/bin/python -m unittest discover -s code/model/experiments/stationary_ces_scales/tests`. Seven engine tests cover CES first-order conditions, capped KKT, allocation budgets, independent interval optimization, owner marginal utility, CD limit/regression and current-child scale/benefit binding. Controller tests separately exercise zero-lifecycle search and repeat gates. No numerical model solve is part of these tests.

The completed one-hour-budget experiment (Torch job 18921574) finished in
48 minutes 38 seconds. Three verified full equilibria cover the same initializer;
no new optimizer proposal completed an equilibrium. The original-weight loss
is 3638.67943503527. All 14 target rows, 31 native parameter rows and 17 standard
diagnostic plot hashes match across the final baseline verification. This is
an incomplete calibration, not an optimized CES outcome or an adopted model.
The exact launched controller is retained unchanged: its `completed_full_ge`
field counts attempt records, including 16 budget-refused proposals, and its
internal penalty is not a computed model loss. The output packet's
`collected/coverage.json` records 19 attempts, three computed equilibria at one
unique parameter vector, 39 lifecycle solves and no numerical rejections.

The corrected deployed archive has SHA256
`bc61bbba879458598d05d250c1b8aa22f02002693b8264af77e757d5a245d561`.
At source-backup closeout, all 18 engine Python files and `run_ces.py` match
that archive byte for byte. Copied refactor verification receipts describe the
inherited reference; they do not certify CES. The focused CES checks and native
initializer/reporting receipts described above provide the CES verification.
Generated Python/Numba caches and operating-system metadata are not source.
