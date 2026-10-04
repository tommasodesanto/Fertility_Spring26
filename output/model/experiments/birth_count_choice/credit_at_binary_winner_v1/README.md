# Fixed-price credit relaxation at the verified Estate-A binary winner

This is one matched **partial-equilibrium** diagnostic at the verified chain-1 one-birth Estate-A point. Both arms use the same ten coordinates, fixed housing coefficient $H_0=6.1245552591467405$, and price $p=0.7811670615311468$. The only changed economic input is the uniform financed share $\phi:0.8\to1.0$. The Estate-A sale-cost rule, revised interest timing, entry and fiscal rules, target 4.45838713455674, weights, and numerical gates are unchanged. No price root or calibration was run.

**Acceptance:** The $\phi=0.8$ baseline passed the native final observer and all inherited gates. Its 14 target model values match the verified selected packet within $7.11\times10^{-15}$, and its ten scored residuals within $8.40\times10^{-14}$ (existing gate $10^{-10}$). Its new-contract loss is 21.2754133610716, matching the selected value 21.275413361071312. The $\phi=1.0$ lifecycle and stationary distribution solved and wrote 17 standard figures, but the unchanged **negative-estate production gate failed**. It is **unaccepted**; the numbers below are provisional saved-solution diagnostics, not an accepted counterfactual or target fit. No relaxed-arm 14-row target table is claimed.

| Fixed-price quantity | $\phi=0.8$ accepted baseline | $\phi=1.0$ provisional | Change |
|---|---:|---:|---:|
| Expected explicit births per normalized household | 0.115230677006 | 0.112736774921 | -2.164% |
| Ownership, all households | 0.696172645339 | 0.818967491400 | 12.279 percentage points |
| Motherhood at completed interview age 25 | 0.446065877195 | 0.435396631197 | -1.067 percentage points |
| Children ever born at age 25, conditional on motherhood | 1.232071010625 | 1.235715129770 | +0.003644 |
| Children ever born at age 25, all households | 0.549584836121 | 0.538026204621 | -0.011559 |

Age 25 uses the active observer projection: 12.5% of the pre-birth and 87.5% of the post-birth distribution in the age-22–25 cell. The baseline calculation exactly reproduces the accepted `early_fertility` observer value. [The age-25 CDF](age25_children_cdf.png) labels the $\phi=1$ curve provisional; [the full age-cell shares](children_ever_born_by_age.csv) retain pre/post timing.

The relaxed arm has net-negative death estates **1.1913129273e-08** per model period, versus the unchanged $10^{-10}$ gate. A death mass **1.27086870611e-07** of total **0.0617334561807** is affected (share 2.06e-06). The material amount is at age 66 in the two-room owner product: the mean shortfall per affected decedent is 0.093740047384, equal to sale cost $\psi pH=0.06\times0.7811670615311468\times2$. At full financing, debt can reach $pH$, while sale proceeds net of selling cost are $(1-\psi)pH$. This explains the gate failure without implying a large aggregate creditor loss. The [read-only reconstructed ledger](phi_100/estate_from_saved_arrays.json) provides totals, ages, and the branch-correct tenure split; the baseline reconstruction matches its accepted native ledger exactly. The [original reporter failure](phi_100/reporting_failure.json) remains intact.

The source is the verified Estate-A binary search receipt SHA-256 `2e104f260050b802370de3d8bdd0d0792c8a7af9e13be849fadfdd6171106f71`. [Prepared inputs](prepared.json) pin the complete experiment code and both input fingerprints, confirm `phi` is the only difference between arms, and record propagation to `shared.phi_choice` and `shared.phi_state`. [Launch](launch_receipt.json), [terminal](terminal_receipt.json), [baseline fit gate](baseline_fit_gate.json), and [read-only summary](diagnostic_summary.json) are retained. The 17 standard plots for each arm are in [baseline](phi_080/standard_diagnostics/) and [provisional relaxed](phi_100/standard_diagnostics/); their housing-market panels are fixed-price diagnostics and do not show equilibrium clearing. Large saved solution arrays remain local and are not included in Git.

## Complete accepted baseline target fit

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| initial_normalization | 2.1 | 2.100000132416349 | 1.324163489968555e-07 |  |  |
| cps_childlessness | 0.19827875100684264 | 0.19937082965426206 | 0.00109207864741942 | 35532.3042455214 | 0.04237709711010588 |
| cps_exactly_one | 0.21365532522014702 | 0.21983341718038182 | 0.006178091960234805 | 26952.820824310795 | 1.0287573737888578 |
| nchs_mean_age | 25.976263860992496 | 25.973635968522142 | -0.0026278924703539985 | 139.82806784479274 | 0.0009656273046881534 |
| nchs_share30 | 0.2492780130410667 | 0.23830064918452423 | -0.010977363856542466 | 0.0 | 0.0 |
| wealth_earnings | 4.45838713455674 | 5.029110096565928 | 0.5707229620091878 | 7.595098472533724 | 2.4739111666101308 |
| bequest_wealth | 0.007291023472616158 | 0.006460953164947379 | -0.0008300703076687789 | 5165289.256198346 | 3.55897063880858 |
| old_dispersion | 3.51593508651872 | 3.249580850854957 | -0.26635423566376293 | 0.0 | 0.0 |
| mean_rooms | 5.729434240102641 | 5.651369141589238 | -0.07806509851340326 | 128.02070205233477 | 0.7801785911672395 |
| ownership_30_55 | 0.6762604168538028 | 0.6794060981371007 | 0.003145681283297952 | 2339.3623724673616 | 0.02314871759988371 |
| first_birth_rooms | 1.465 | 1.2462912282603833 | -0.21870877173961678 | 137.5652749002964 | 6.580232268624658 |
| family_rooms | 0.38509964969278165 | 0.30252636167737723 | -0.08257328801540442 | 0.0 | 0.0 |
| recent_parent_ownership | 0.12760836356692162 | 0.12865865900004736 | 0.0010502954331257364 | 27055.822957508266 | 0.029845832863430875 |
| early_fertility | 0.8095276384290021 | 0.5495848361206093 | -0.25994280230839284 | 100.0 | 6.75702604719402 |

All 31 accepted baseline parameter records, including estimates and advisory bounds, are in [parameters.csv](phi_080/reporting/phase_b_ge/phi_080/parameters.csv). The ten coordinates retain their Estate-A search bounds; the native table shows the reference advisory beta bounds 0.94–0.99, while the actual Estate-A search bound was 0.93–0.99. The relaxed arm holds all 31 numeric values fixed.

## Reproduce without another solve

`output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/postprocess_credit_at_binary_winner.py` reconstructs the ledger and age-25 distributions from saved arrays. The original two-solve command is `NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/run_credit_at_binary_winner.py --run`, with an external 1,500-second process limit; it refuses duplicate arm directories. `--prepare` authenticates inputs with zero solves.
