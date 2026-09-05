# Candidate policy comparison for September 14

Author authorization: September 5, 2026, three-hour independent work window
starting 21:11 UTC, ending 00:11 UTC September 6. The author specifically asks
that apparent inconsistencies be traced through prior chats before diagnosis.

## Fixed question and scope

Does the verified morning candidate change the existing supply, family-credit
and property-tax conclusions? Compare with retained September 4 task_010, using
the same twelve targets, eleven free parameters, policy definitions and finite
2023--2063 closure. No new calibration search or production promotion. The
rebated tax calculation is a separate 2023 comparison of rebated 2% versus
rebated 1%; it is not relative to the unrebated status quo.

Use the already repaired source snapshot at
`/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/`.
Its history replay exactly reproduces the morning candidate's twelve target
rows, parameter table and 253 numeric historical entries. The archived original
source, active working checkout and separate structural-cleanup task stay intact.
The new observer copies each solved state for checkpoints and the unchanged
seventeen-graph packet. It does not change any model or policy equation.

## Recoverable run budget

- Smoke: all five policy cases, 2023 and 2027, parallel at most five. About
  5--10 minutes each based on the completed historical solves; hard 40 minutes.
- Full: the same five cases, eleven dates through 2063, after every smoke gate
  passes. Prior full-case runtimes were 807--974 seconds. Allow 20--30 minutes
  including the observer; hard 55 minutes per case. At most five simultaneous
  tasks. Approximately 5 old equilibria, 20 historical dates and 55 policy dates;
  the price solver may perform several Bellman evaluations at each date.
- Separate tax calculation: three impact equilibria then eight fixed-component
  cells for the exact Shapley decomposition; hard 55 minutes, launched only
  after the smoke. No new long-horizon rebated or perfect-foresight run.
- One run at each stage; numerical failure stops dependent work. A diagnosed
  preflight problem may be corrected explicitly, retaining its failed receipt.
  No automatic parameter search or wider policy menu. Latest launch 22:50 UTC;
  finish collection and the readout within the three-hour window.
- Heartbeats every minute and a completed-date receipt at every policy date;
  latest-completed and fixed-best summaries stay readable. No heartbeat for
  thirty minutes triggers investigation. Every policy date saves a compressed
  checkpoint and the stable seventeen diagnostic graphs. The helper's historical
  `lifecycle_2023.csv` filename is retained inside an explicitly dated folder.

## Contract reconciliation before launch

| Object | Classification | Maintained convention |
|---|---|---|
| Eleven preferences/supply parameters | Estimated | Verified morning candidate, unchanged twelve-row objective |
| Dated supply elasticity | Externally fixed | 0.63 |
| Old-equilibrium supply elasticity | Outstanding economic interpretation | Inherited 1.75; same as both verified calibrations |
| Date-0 supply normalization | Empirically normalized initialization | Clear reweighted 2007 household distribution at inherited price |
| Tenure taste dispersion | Externally fixed | 0.005 |
| Old fertility level | Externally fixed replacement normalization | Derive preference intercept for completed fertility 2.1 |
| Historical household mass and age distribution | Empirically normalized | Census HH-3 totals, national ACS head-age shares through 2023 |
| Birth-to-household conversion | Externally fixed; population interpretation outstanding | Adjusted births / 2.1, twenty-year delay, four waiting slots |
| Post-2023 entry/geography | Externally fixed mechanism closure | Closed national comparison, outside inflow M=0, retention rho=1 |
| Outside-origin share in old accounting | Externally fixed legacy restriction; empirical interpretation outstanding | 0.169 retained; no new estimated migration response |
| Baseline tax and fiscal scope | Externally fixed | Annual 1%; receipts discarded in five-case packet |
| Family credit | Externally fixed policy experiment | Dependent-child households may finance 95%; no purchase grant |
| Rebated-tax comparison | Externally fixed fiscal experiment | Equal contemporaneous rebate, rebated 2% versus rebated 1% |
| Expectations and horizon | Externally fixed mechanism experiment | Temporary-equilibrium sequence through 2063; no welfare or terminal steady-state claim |

This classifies the inherited benchmark transparently; it does not close the
outstanding items or promote it to the separate coherent-person-law branch.
The historical evidence table distinguishes author decisions from agent records.

## Completed results and numerical follow-up

All jobs completed with exit 0:0: smoke array `17009623` (6m28s--7m09s), full
array `17009732` (16m52s--17m48s), tax/Shapley `17009733` (14m08s), and the
six fixed-price saving checks `17009779` (2m10s). Total allocated CPU time is
2.27 hours. Both existing collectors passed. All 55 full dates reproduce their
smoke prefixes exactly, retain a common 2023 state, clear markets, pass mass and
feasibility gates, and never use a grid fallback.

Candidate births per household, relative to its own no-policy path:

| Policy | 2023 | 2063 | Total births, 2063 |
|---|---:|---:|---:|
| Supply schedule +20% | +1.343280% | +2.140588% | +2.567723% |
| Dependent-child LTV 95% | +0.017992% | +0.074818% | +0.082364% |
| Combined | +1.365519% | +2.244344% | +2.682037% |
| Unrebated property tax 1% to 2% | -1.144991% | -1.698748% | -2.019930% |

The rebated 2% versus rebated 1% impact decomposition is -1.376585% tax,
+0.301906% capitalization, +1.574415% rebate, net +0.499736%. The selected
rebate and parameter-state transfers agree exactly in both equilibria. The
possible cached-state concern raised by inspection was not realized and was
not declared a new solver bug. No calibration, target or production promotion.

The full value-monotonicity screen flags 28/55 dates, with maximum exposed
pre-choice mass share 1.773586e-6 (combined, 2031). The eleven tax packets flag
four, maximum 2.581302e-6 (a mixed component cell). Three smoke distributions
were re-solved with the local and previously reviewed exhaustive saving methods:
local arrays/distributions reproduce exactly; global occupied values dominate
and all occupied value drops disappear in those three tests. The largest
fixed-price births difference is 0.000345945%. This is not an equilibrium-path
or Shapley error bound. The largest reporting-budget-excess share is 2.1554e-8.

`artifact_manifest.json` verifies all 1,368 original checkpoint/graph hashes on
Torch and records 2,401 output files. `collection_receipt.json` verifies 2,319
local files; all 82 large checkpoints remain on Torch. There are 1,394 standard
graphs across 76 policy/tax packets and six saving-check packets. Impact and
terminal five-policy packets and all three tax-equilibrium packets were visually
inspected; the report's supplemental comparison does not replace that stable set.

The rebate uses `rebate_plan.json`, a child of the original plan. Its only
addition is a read-only observer for the three equilibria and eight component
cells, plus a selected-transfer/parameter-state consistency assertion. Mixed
component cells are fixed-component diagnostics, not clearing equilibria.
Historical and source reconciliation is in `historical_reconciliation.md`,
including the verified reporting-only policy-driver diff from Git history.

## Reproduction and evidence

`prepare.py` produces the transparent selected-replay report and immutable plan
from existing verified inputs. It does not change the candidate summary or invent
a new calibration search. `plan.json` pins sources, input hashes, selected history,
target fingerprint, policies and deadline. The existing policy driver and
collector retain their gates. `run_e5f_candidate_policy_comparison.py` adds only
observation and exact smoke-prefix checks. `run_e5f_candidate_policy_comparison.sh`
is the bounded Torch launcher. Submission and collection receipts are saved here.

The selected input view uses the repaired replay's `task_001` receipt, which
exactly reproduces morning `round2/task_009`; it is not a third calibrated point.
The immutable plans retain their launch deadlines and completed output paths.
A new numerical reproduction needs a new explicitly scoped output/plan, rather
than overwriting these results or weakening their no-overwrite checks.

Large state checkpoints are retained on Torch. To collect the graphs and all
small evidence files into this named local folder:

```sh
rsync -az --exclude='*.pkl.gz' torch:/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/output/model/e5f_candidate_policy_comparison_20260905a/ output/model/e5f_candidate_policy_comparison_20260905a/
```

The numerical folders and immutable plans are not overwritten by report generation.
The documented graph regeneration command was run: all seventeen baseline-2063
PNGs reproduce exactly (`graph_regeneration_probe.json`).

After collecting the saved outputs, from the local repository root:

```sh
python3 output/model/e5f_candidate_policy_comparison_20260905a/verify_artifacts.py --local
python3 code/model/tools/build_e5f_candidate_policy_comparison_report.py
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_candidate_policy_comparison_review.md --output output/pdf/e5f_candidate_policy_comparison_review.pdf --date '5 September 2026' --heading 'FERTILITY / CALIBRATION AND POLICY' --no-source-index
```

For a single saved solution, the existing function
`run_e5f_independent_numerical_audit.standard_diagnostics(packet, out,
validate_production_young=False)` regenerates all seventeen standard graphs.
On Torch, in the frozen project root, one command is:

```sh
PYTHONPATH=code/model:code/model/tools python3 - <<'PYTHON'
from pathlib import Path
import hashlib, json
import run_e5f_independent_numerical_audit as audit
folder = Path('output/model/e5f_candidate_policy_comparison_20260905a')
relative = 'full_diagnostics/baseline/2063/dated_state.pkl.gz'
checkpoint = folder / relative
assert hashlib.sha256(checkpoint.read_bytes()).hexdigest() == json.loads((folder/'artifact_manifest.json').read_text())['files'][relative]
_, configured = audit.transition.configure_sequential_model()
assert audit.baseline.calibration.code_fingerprint_contract(configured)['bundle_sha256'] == json.loads((folder/'plan.json').read_text())['code_bundle_sha256']
out = Path('tmp/candidate_baseline_2063_regenerated')
out.mkdir(parents=True, exist_ok=False)
audit.standard_diagnostics(audit.load_checkpoint(checkpoint), out, validate_production_young=False)
PYTHON
```

The concise seven-page PDF is `output/pdf/e5f_candidate_policy_comparison_review.pdf`;
its first page provides the decision readout, and pages 5--6 contain complete
fit and parameter tables. Full numerical precision is in the source CSVs.
