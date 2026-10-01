# National benchmark calibration at normalized household population one

Prepared for the author-confirmed October 1 contract. The national 2007 benchmark
has $N_0=1$. Housing supply is $S(p)=H_0(up/\bar r)^\xi$; $H_0$ is internally
calibrated, with the existing national AHS2007 mean occupied rooms target for ages
18–85, 5.729434240102641. There is no price/rent calibration target or external
anchor for $H_0$.

The implementation profiles the same calibration. For each ten-coordinate
proposal, including free child-benefit level $\psi$, the unchanged native price
root sets $B(p)/(2.1E(p))=1$. Native per-household housing demand $h(p)$ determines
$H_0=h(p)/(up/\bar r)^\xi$. The mean-rooms moment remains scored with its original
weight: this formula clears supply at $N_0=1$ and does not force rooms to its
target. The original $H_0\in[0.2,80]$ restriction is checked at accepted roots,
not intermediate bracketing prices. A violation is an explicit constraint
rejection with the actual lifecycle-call count.

The adapter preserves household policies and distributions and performs no extra
lifecycle solve. Before native gates/reporting it copies $P$ with normalized
$H_0$, reconstructs supply from that coefficient, and updates evaluation supply,
solution supply, aggregate and market residuals, expected parameter tables,
closure receipts, retained selected/repeat arrays and the unchanged 17-plot set.
It checks $P.N_{target}=1$, actual native household mass and solution mass within
the existing $5\times10^{-9}$ tolerance, and a finite positive supply factor.
Renewal root tolerances, fiscal/native gates, model code and SMM scoring are
unchanged. Both fast search and complete native final-report routes use the
isolated adapter, with reviewable literal source diffs written by each process.

Economics relative to the preceding free-$\psi$ floor search are unchanged:
physical first-child housing floor, $D=0$, nonnegative mean-preserving entry,
120×9 grid, owned rooms 2/4/6/8/10, original targets/weights and coordinate
bounds. Fixed survival is one before reproduction ends; actual birth renewal
and completed fertility therefore use the same topcode correction
3.602359422009. No experimental point is adopted by this packet.

The authenticated incumbent is loss 30.99909262396404 at price
0.7188481029285891. `center/` retains its complete 14-target/31-parameter tables
and verification receipt. `plan.json` contains the exact incumbent start plus
23 distinct deterministic modest joint neighbors (seed 20261005); all chains
use original base weights. Each chain gets 7200 seconds from actual launcher
start, at most 100 objective calls, a 900-second selected-native reserve, one
CPU/24 GiB/one thread and the existing maximum 32 lifecycle calls per candidate.
Maximum production work is 2400 exploratory GEs plus 24 selected native
postchecks. Observed median GE times 128.38–166.81 seconds imply about 37–49
cases per chain before the reserve; the time budget can bind first. The nominal
2400-case work is 85.6–111.2 CPU hours but the 24 production clocks cap allocated
CPU time at 48 hours. No automatic retries or extensions.

Prepared verification:

```sh
PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/normalized_calibration_v1/test_configuration.py
```

This four-second check performs zero model solves. `configuration_check.json`
records normalization arithmetic, unit-mass/factor checks, accepted-root bounds,
full/fast source variants, all 24 actual mocked 100-call controller loops,
reserve behavior, constraint accounting and fatal propagation.

The production gate is one Torch native smoke at the exact incumbent, with full
selected-root/repeat gates, complete 14/31/17 artifacts, and all target-fit values
and non-$H_0$ estimates matching the authenticated incumbent within $10^{-10}$:

```sh
python run_psi.py --chain 0 --out NEW_SMOKE_DIRECTORY --deadline-epoch ACTUAL_START_PLUS_7200 --smoke-only
```

Use the deployed full package path and one-thread launcher environment. Native
smoke runtime is not yet measured; the preceding native repeat is expected to
require minutes, with the unchanged 300-second per-lifecycle limit and 32-call
cap. Production searches use `--fast-objective`; native postchecks use
`--verify-only SEARCH_DIRECTORY/search_completed.json`. Deployment owns the
Torch source-only overlay, staging and gated submission. No local native model
solve was run during preparation.
