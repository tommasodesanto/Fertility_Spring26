# HANDOFF: Independent evaluation of the M4 calibration process and result

Repository root:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

I want an independent, critical assessment of how the M4 standard-bequest
calibration went. This is an audit, not a continuation of the search. Do not
change files, launch jobs, rerun a broad calibration, commit, or push. Small
read-only commands and narrowly targeted numerical checks are allowed. Do not
accept the existing status notes or the prior agent's conclusions without
checking the underlying artifacts.

## 1. Mandatory orientation

Read the project context in the required order:

1. `memory/AGENT_MEMORY.md`
2. the latest `memory/daily/YYYY-MM-DD.md`
3. `CALIBRATION_STATUS.md`
4. `code/model/README.md`
5. `SESSION_DIARY.md` only if needed to resolve chronology

Then read:

- `AGENTS.md`
- `docs/style/econ_writing_style_guide.md` only if assessing paper language
- `docs/prompts/HANDOFF_m4_standard_bequest_20260716.md`, but treat it as
  historical evidence: much of that file predates the user's correction that
  `theta1` must be calibrated internally.
- `/Users/tommasodesanto/.codex/attachments/a3f46331-bc03-4801-bf18-a1dc1eabcf1f/pasted-text.txt`
  for the preceding chat and the user's intended contract.

First run `git status -sb` and preserve the dirty worktree. This is a read-only
review.

## 2. Final intended M4 contract

The final contract—not the earlier fixed-`theta1` proposal—is:

\[
B(W)=\theta_0
\frac{(\theta_1+W)^{1-\sigma}-\theta_1^{1-\sigma}}{1-\sigma},
\qquad W=\max\{b+pH,0\},\qquad \sigma=2.
\]

In code this should use the existing `linear_child_scale` specification with
`theta_n=0`, `normalize_bequest_utility=True`, SSA age survival enabled, and no
owner-LTV taper. `tenure_choice_kappa=0` and `theta_n=0` are external
restrictions. Both `theta0` and `theta1` are estimated; `theta1=0.25` was only
one optimizer start.

There are 14 target moments and 13 free parameters: the 11 clean-frontier
parameters plus `theta0` and `theta1`. The two late-life targets are:

- estate wealth/income median, ages 76--84:
  target `6.50131577436537`, SE `0.2319582116443`, weight
  `18.585767349158665`;
- nonhousing wealth/income median, ages 65--75:
  target `1.90821154211154`, SE `0.109272216032637`, weight
  `83.74916751466371`.

The predeclared fit gates were:

1. two identical strict exact evaluations;
2. estate median within one SE;
3. nonhousing median positive and within two SEs;
4. loss on the 12 established moments no greater than `4.42`;
5. free `theta0` weakly improves on the exact `theta0=0` nested seed;
6. an identification assessment sufficient to justify estimating both bequest
   parameters internally.

Audit whether this contract is actually implemented and whether the last item
was ever defined sharply enough to constitute an ex ante gate.

## 3. Chronology to verify, not assume

The reported chronology is:

1. Six initial production chains finished cleanly but selected an exact winner
   with loss `12.6540621356`, `theta0=0.24062`, and `theta1=0.82286`. This
   failed the established-block gate.
2. The exact nested `theta0=0` reference had loss `307.6847778879`.
3. The first collector failed because it compared `0.25` with
   `0.25000000000000006`. The collector was changed to tolerate only a
   four-ULP transform roundtrip; exact moment-repeat requirements were not
   relaxed.
4. A fixed-`theta1` profile found much better points:

   | Fixed `theta1` | `theta0` | Exact loss |
   |---:|---:|---:|
   | 0.421710877 | 0.247497728 | 4.468698900 |
   | 0.589073251 | 0.259821410 | 4.571466636 |
   | 0.822855927 | -- | 12.654062136 |
   | 1.149418812 | -- | 22.466503693 |
   | 1.605583144 | -- | 80.980478092 |

5. A full 13-dimensional polish starting from the `theta1=0.421710877`
   profile point returned the final exact loss `4.401508692086581`, with
   `theta0=0.247497728306924` and `theta1=0.4217108770366293`. The reported two
   strict repeats are bit-identical.
6. The final fit decomposition is:
   established 12-moment loss `3.6497057301`, estate contribution
   `0.0002041425`, and nonhousing contribution `0.7515988195`.

Determine whether this is a defensible search-and-promotion procedure or
whether the successful point is too dependent on an inadequately explored
profile/local basin. In particular, the fixed profile did not include a cell
below `theta1=0.4217`. Assess whether the final free polish supplies credible
evidence of an interior optimum or merely stayed at its starting value because
the direction is extremely flat.

## 4. Canonical artifacts

Inspect at least:

- `output/model/intergen_standard_bequest_recalibration_20260716/final_report/results.json`
- `output/model/intergen_standard_bequest_recalibration_20260716/final_report/target_fit_full.csv`
- `output/model/intergen_standard_bequest_recalibration_20260716/final_report/parameter_table_full.csv`
- `output/model/intergen_standard_bequest_recalibration_20260716/final_report/acceptance_criteria.csv`
- `output/model/intergen_standard_bequest_recalibration_20260716/final_report/untargeted_diagnostics.csv`
- `output/model/intergen_standard_bequest_recalibration_20260716/theta1_profile_report_final/`
- `output/model/intergen_standard_bequest_recalibration_20260716/identification_final/summary.json`
- `output/model/intergen_standard_bequest_recalibration_20260716/identification_final/weakest_direction.csv`
- `output/model/intergen_standard_bequest_recalibration_20260716/distribution_diagnostic_final/README.md`
- `output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet_final/README.md`
- the two PNGs under
  `output/model/intergen_standard_bequest_recalibration_20260716/diagnostic_packet_final/classic_draft/`

Trace the artifacts back to the active code:

- `code/model/intergen_housing_fertility/calibration.py`
- `code/model/intergen_housing_fertility/parameters.py`
- the bequest implementation in the active solver/kernels
- `code/model/tools/run_intergen_bequest_exit_chain.py`
- `code/model/tools/collect_intergen_internal_bequest_recalibration.py`
- `code/model/tools/collect_intergen_standard_bequest_theta1_profile.py`
- `code/model/tools/audit_intergen_bequest_exit_jacobian.py`
- the relevant cluster submission scripts
- `code/model/intergen_housing_fertility/tests/test_bequest_target_moments.py`

Verify the mathematical bequest object, estate definition, target measurement,
parameter domain, fixed restrictions, objective weights, strict repeat logic,
and provenance of every reported final number.

## 5. Questions the review must answer

### A. Contract and numerical integrity

- Does the code implement the final 14-moment/13-parameter M4 contract exactly?
- Is `theta1` genuinely free everywhere that matters, rather than accidentally
  fixed or inherited from a profile seed?
- Does `W=max{b+pH,0}` in the mathematical description match the exact object
  used in terminal utility and the measured estate moment?
- Are the target values, age windows, denominators, bootstrap weights, and
  signs consistent between empirical construction, calibration code, and
  collector?
- Are the final result, exact repeats, plots, Jacobian, and distribution packet
  all generated from the same parameter vector and target system?
- Was the four-ULP collector fix appropriately narrow, or could it hide a real
  contract mismatch?

Check three possible strictness gaps reported by a separate read-only audit:

1. the production-chain collector may validate arm, target set, and counts but
   not the literal M4 domains, fixed restrictions, evaluators, or six distinct
   starts/seeds;
2. nested-reference validation may omit `active_domain` and
   `fixed_parameters`;
3. a focused domain test may compare against constants imported from the same
   module rather than literal required tuples, allowing coordinated drift to
   pass.

Confirm or reject each concern from the code and explain whether it affects the
reported result or only future reproducibility.

### B. Search quality

- Why did all six initial chains stop at loss `12.65` while the subsequent
  profile immediately found `4.47`?
- Was the original search design inadequate, was there a seeding/geometry
  problem, or is there evidence of a driver/collector discrepancy?
- Is the final `4.4015` point sufficiently validated to freeze, or is one
  narrowly designed low-`theta1` profile/check indispensable?
- Do not recommend a broad new search by reflex. State the smallest decisive
  check, if any, and what outcome would change the conclusion.

### C. Identification

The final Jacobian reportedly has 14 rows and 13 columns. Under SMM weighting,
its relative-threshold rank is `9/13` at `1e-2` and `12/13` at `1e-3`, with
condition number about `3705`. The weakest singular vector loads about `0.996`
on `theta1`; the `theta1` column norm is far smaller than the other parameter
columns. The nearby fixed-profile cell at `theta1=0.5891` raises loss by only
about `0.103` relative to the fixed-profile minimum.

- Can `theta1` honestly be described as internally calibrated?
- Can it be described as separately or precisely identified?
- Is the system formally overidentified but effectively rank-deficient?
- Should the paper report a point estimate, a profile/set, an external
  restriction, or a sensitivity analysis?
- Which existing or additional moment would economically identify `theta1`
  separately from `theta0`, if any?

### D. Economic fit and paper readiness

Do not judge the model only by the targeted scalar loss. Discuss at least:

- aggregate ownership: target `0.575`, model `0.715`;
- ages 25--34 ownership: target `0.341`, model `0.344`;
- old-age ownership: model about `0.964` and a visibly too-steep lifecycle
  ownership path;
- estate median: target `6.501`, model `6.505`;
- nonhousing median at ages 65--75: target `1.908`, model `2.003`;
- ages 76--84 estate p90/p50: PSID `5.291`, model `1.763`;
- ages 76--84 nonhousing p50/p75/p90: PSID
  `2.368/10.336/26.272`, model `0/0/0`;
- top-estate-decile housing share: PSID `0.186`, model `0.996`;
- degenerate model retirement income in the current distribution diagnostic.

Explain whether passing the two targeted wealth medians while failing the
late-life balance-sheet distribution is a tolerable limitation, a reason to
reframe the claims, or evidence that the mechanism remains substantively
wrong. Distinguish what must be fixed before showing an advisor, before a
paper draft, and before quantitative counterfactuals.

### E. Comparability and communication

The older M1 headline loss around `6.86` used a different full target system.
On the 12 shared established moments with M4 weights, M1 has loss about `3.841`
and M4 has `3.650`. Evaluate whether prior descriptions of M4 as returning the
established block to “M1 quality” are accurate. Flag any place where losses
from incompatible objectives, samples, or target sets were compared as though
they were directly comparable. Confirm that no loss was mistakenly described
as a percentage.

### F. Was fixing `kappa_t` a contract failure?

Tommaso has now flagged that he expected `kappa_t` to be calibrated and is
unhappy that the reported M4 instead fixes it at zero. Treat this as a central
audit question, not an innocuous modeling convention.

The repository has an important notation collision that must be resolved
explicitly:

1. several numerical/model documents write the tenure-choice logit scale
   `tenure_choice_kappa` as \(\kappa_t\); positive values smooth the rent/own/
   housing-product choice and zero gives a deterministic hard argmax;
2. a July D12 experiment also used `kappa_t` for a proposed child time-cost
   parameter, which was reportedly rejected at zero and is not part of M4;
3. `kappa_fert` is a different parameter: the fertility-logit scale, estimated
   in M4 at about `2.14745`.

Do not conflate these three objects. For tenure `kappa_t`, answer:

- Was fixing `tenure_choice_kappa=0` explicitly requested or approved by
  Tommaso for M4, or was it silently inherited from the M1 “clean frontier”?
  Reconstruct this from the transcript, status notes, run contract, and code.
- Earlier production specifications searched `tenure_choice_kappa` over a
  positive interval (and some documentation used `0.01`). Why and when was it
  demoted to an external zero restriction? Was that decision based on economic
  evidence, numerical convenience, a prior boundary estimate, or an attempt to
  keep the parameter count below the moment count?
- Did the written proposal's phrase “11 clean-frontier parameters” make this
  restriction sufficiently transparent, or should the agent have presented
  the complete free/fixed parameter table before launch? Assess whether the
  user could reasonably have believed `kappa_t` was still being calibrated.
- Is fixing `kappa_t=0` defensible economically? Distinguish a structural
  random-utility scale from a numerical smoothing device. Check how the active
  Bellman equation and forward distribution use it; if positive, it changes
  realized tenure/product probabilities and is not merely a solver trick.
- Could the zero restriction be materially related to M4's main fit failure:
  aggregate ownership `0.715` versus `0.575` and old-age ownership around
  `0.964`? Use existing `kappa_t` sweeps or cached diagnostics if available;
  do not assert a direction without evidence.
- Does fixing `kappa_t` make the reported M4 objective artificially favorable
  or simply define a different model? Would freeing it plausibly improve the
  ownership path while damaging housing-rung, fertility, or wealth moments?
- If tenure `kappa_t` is restored, count the free parameters and informative
  moments exactly. The current M4 has 14 moments and 13 free parameters but an
  effective weighted Jacobian rank of only 12 at the `1e-3` threshold. Explain
  why merely adding a fourteenth parameter to fourteen moments does not solve
  identification. Name the moment or external restriction that would identify
  tenure smoothing separately from housing preferences, discounting, and the
  bequest block. Consider tenure transition hazards or the age profile of
  ownership, not only the aggregate ownership level.
- Decide whether the current M4 can still be called the requested calibration,
  should be relabeled “M4 conditional on deterministic tenure,” or should be
  treated as incomplete until `kappa_t` is addressed.
- If a corrective computation is necessary, specify the smallest decisive
  design: exact parameter set, added/replaced moment if any, bounds, starting
  points, comparison baseline, runtime class, and pass/fail criterion. Do not
  recommend an open-ended recalibration.

For the separate child time-cost parameter that was also called `kappa_t`,
verify the D12 rejection before accepting it. State whether that experiment has
any bearing on the tenure-choice scale (it should not unless the code or notes
actually link them).

## 6. Deliverable

Lead with a one-paragraph verdict using one of these labels:

- **A: numerically and scientifically ready**;
- **B: numerically accepted, scientifically usable with explicit caveats**;
- **C: provisional; one bounded check is required**;
- **D: not credible; reopen the calibration or model specification**.

Then provide:

1. a verified chronology of what happened;
2. what was done well;
3. mistakes, confusion, or avoidable process failures;
4. a full 14-moment target/model/gap/weight/contribution table;
5. a full 13-parameter estimate/bound/near-bound table;
6. the search-quality verdict;
7. the identification verdict for `theta0` and `theta1` separately;
8. a separate verdict on whether fixing tenure `kappa_t=0` was authorized,
   scientifically defensible, and material to the result;
9. the economic-fit and paper-readiness verdict;
10. a short prioritized action list divided into:
   - required before showing the advisor,
   - required before using M4 in the paper,
   - optional robustness work;
11. the exact two or three sentences Tommaso should use to describe this
    calibration honestly to an advisor.

Every substantive claim should cite an exact local file and, when possible, a
line number, CSV row, JSON key, or code symbol. Separate verified facts from
your interpretation. If canonical files disagree, stop and make that
discrepancy the headline rather than silently choosing one.
