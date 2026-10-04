# Net-estate-solvent fixed-price credit diagnostic at the verified Estate-A winner

This v2 experiment holds the verified one-birth Estate-A winner at price $p=0.7811670615311468$, housing coefficient $H_0=6.1245552591467405$, and the same ten estimated coordinates. The matched arms differ only in the uniform financed share $\phi=0.8$ versus $\phi=1.0$. No price root or recalibration was run.

**Experimental economic change from v1:** A household facing possible death must choose post-saving liquid assets $b\prime\geq-(1-\psi)pH$, where $H$ is its chosen owner-room product and $\psi=0.06$ is the Estate-A selling cost; renters have $H=0$ and require $b\prime\geq0$. This is maximized with the existing credit and asset-grid floors for buyers, tenure changers, and stayers in both full and fast savings paths. At mortality ages, nominal $\phi=1$ can therefore finance at most $(1-\psi)=0.94$ of the house value through secured debt. The rule also applies in the terminal period. Earnings, initial wealth/income distributions, timing, transfers, preferences, targets, weights, and numerical gates are unchanged. This is experimental, not an adopted production default.

**Acceptance:** Both arms completed native partial-equilibrium (PE) observation and the unchanged estate, purchase, fiscal/PAYGO, household-budget, stationary-distribution, and feasibility gates. Both audited net-negative death estates equal zero; the [relaxed native gate ledger](phi_100/reporting/phase_b_ge/phi_100/gates.json) and [read-only reconstruction](phi_100/estate_from_saved_arrays.json) agree exactly. This is a fixed-price PE comparison: the $\phi=1$ renewal residual is $-0.0216386158176$, so it is **not a general-equilibrium result**. The final 14-target reporter requires the GE renewal root and was appropriately run only for the 0.8 arm. No 14-row relaxed target fit or relaxed objective loss is claimed.

| Fixed-price quantity | $\phi=0.8$ baseline | $\phi=1.0$ net-estate-solvent PE | Change |
|---|---:|---:|---:|
| Expected explicit births per normalized household | 0.115230677006 | 0.112736504723 | -2.164504% |
| Ownership, all households | 0.696172645339 | 0.818965240629 | 12.279260 percentage points |
| Motherhood at completed interview age 25 | 0.446065877195 | 0.435396244658 | -1.066963 percentage points |
| Children ever born at age 25 conditional on motherhood | 1.232071010625 | 1.235715218499 | +0.003644208 |
| Children ever born at age 25, all households | 0.549584836121 | 0.538025765601 | -0.011559071 |
| Adjusted births per normalized household in native PE observer | 0.129640266154 | 0.126835022243 |  |
| Net-negative death estates per model period | 0 | 0 | 0 in both |

The age-25 measure projects the age-22–25 cell with $0.125$ pre-birth and $0.875$ post-birth weights, exactly matching the active baseline observer. See the [age-25 CDF](age25_children_cdf.png), [full age-cell count shares](children_ever_born_by_age.csv), and [read-only diagnostic summary](diagnostic_summary.json). The v1 unconstrained $\phi=1$ solution failed the estate gate with $1.1913129273\times10^{-8}$ net-negative estates per period; its raw births and ownership were 0.112736774921 and 0.818967491400. The v2 constraint removes that failure while changing these aggregate quantities by only $-2.70\times10^{-7}$ and $-2.25\times10^{-6}$, respectively.

**Baseline replay:** All 14 native target values match the verified selected packet within 7.11e-15; ten scored residuals agree within 8.4e-14 against the existing $10^{-10}$ gate. New-contract baseline loss is 21.2754133610716, matching the selected 21.275413361071312. All 31 baseline parameter records, estimates, advisory bounds, and bound flags are in [parameters.csv](phi_080/reporting/phase_b_ge/phi_080/parameters.csv); all 31 numeric values are held in the relaxed arm. The actual Estate-A search beta bound was 0.93–0.99 even though the native parameter table retains a 0.94–0.99 reference advisory bound.

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

## Source identity and reproduction

The source Estate-A binary search receipt has SHA-256 `2e104f260050b802370de3d8bdd0d0792c8a7af9e13be849fadfdd6171106f71`. [Prepared inputs](prepared.json) pin the source search receipt, selected target and weight fingerprints, winner coordinates, both effective input fingerprints, `shared.phi_choice` and `shared.phi_state` propagation, and hashes of all experiment model and engine Python sources used in the solve. This is a derivative engine correction from pre-change Git commit `a2a3d4c9`; the older `engine_inventory.json` describes ancestry and is not a claim of byte-identical v2 behavior. [Launch](launch_receipt.json) and [terminal](terminal_receipt.json) receipts record the successful fresh run; an initial interrupted attempt was retained locally and is not included in the result. [Baseline native fit](phi_080/reporting/phase_b_ge/phi_080/target_fit_new_contract.csv), [baseline estate ledger](phi_080/estate_from_saved_arrays.json), [relaxed estate ledger](phi_100/estate_from_saved_arrays.json), and 17 [baseline](phi_080/standard_diagnostics/) and 17 [relaxed](phi_100/standard_diagnostics/) standard plots are provided. Price and market-clearing plots should be read as fixed-price diagnostics. Large saved solution arrays remain local and are excluded from Git.

From the repository root, use `NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/postprocess_credit_at_binary_winner.py --output-version v2` to reproduce the ledgers and age-count reductions from saved arrays without a model solve. The two-solve driver is `NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/run_credit_at_binary_winner.py --run --output-version v2` and an external 1,500-second deadline; it refuses duplicate arm directories. Each native solve has a 600-second limit.
