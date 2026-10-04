# Original Claude scratch evidence recovered

**Recovered October 4, 2026.** No model run, cluster action, source change, adoption or manuscript edit occurred. Original files were copied read-only and SHA-256 checked. Recovery is complete for the bounded original scratch source; executed runtime identity remains incomplete for fixed-price experiments.

## Source identity

The suggested Desktop scratch directory exists but is empty. Its CLI session `b8403d74-a3d8-42bd-b700-d8f8b9bffa5c` is an unrelated GeoComply uninstall conversation. The permitted one-tier search of local session indexes/recent project sessions found the relevant session:

- Session: `/Users/tommasodesanto/.claude/projects/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/e6ffcec8-ecc2-4f4e-80e3-a9983d1b3762.jsonl`.
- Original scratch: `/private/tmp/claude-501/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/e6ffcec8-ecc2-4f4e-80e3-a9983d1b3762/scratchpad`.
- Preservation: `preserved/` includes all 28 scratch top-level scripts/results and the complete 56 MB `ge_runs/` tree. `SOURCE_MANIFEST.json` records every original path, SHA-256, bytes and modification time. The filtered transcript retains matching script records and their tool results with original line numbers.

## What is established by recovery

The saved credit tables support a small **fixed-price** aggregate birth effect at the claimed chain-13 price. The baseline paired case raises financed share from .80 to .95: explicit births change −0.250%, while ownership at ages22–25 rises11.419 percentage points. Child room requirements of0/.5/1 and renter caps6/4 give birth effects from−0.250% to+.978%; tightening the cap and adding rooms are additional experimental economic changes. Dropping the two-room owner rung gives much smaller ownership changes than the general11–15point summary. See `CLAIMS.csv` (12rows) and `ARITHMETIC_CHECK.json`; they report source recovery and arithmetic, not numerical certification.

The strongest native GE attempt did **not** pass: cap4, one additional room per child, re-normalized child value .3222148912, with each financed-share arm failed its dated-budget gate at the first trial price. Bad masses were1.632e−7/5.261e−8 and maximum excesses.223/.083. The exact failed reports are preserved in `ge_runs/ge_cap4_hb1_phi80/` and `phi95/`, as are `ge_strongest_case_summary.json` and its log. The original baseline native GE passed, including an exact repeat, at price.7760569760205563 and population1.

`ge_diagnostic_root.py` expressly skips those reporting/budget audits and uses its own secant root. Its default .95 result is price.7723762061 and household scale.9852667009. These results are uncertified diagnostics; small renewal residuals do not erase the retained budget failures. No accepted transition appears in the recovered experiment. The stationary birth-renewal condition uses adjusted births, while the headline script birth measure is the sum of explicit first/second/third flows; these are distinct.

The tax experiment is also its own diagnostic root. Doubled tax without/with rebate yields gross-supply population.998510795/1.244434855. Replacing the supply argument with net-of-tax rent **algebraically at the same solved prices and household policies** yields.874114729/1.089401177, reproducing−12.6%/+8.9%. No separate net-supply household/root solve was recovered or needed by this reaccounting. Rebated rooms fall8.591% in levels, or8.982% in logs. The root stops after20iterations without a final residual assertion, and the saved tax JSON does not retain residuals. The fiscal rebate is iterated six times and remains a diagnostic fiscal closure.

## Executed input and measurement contract

All fixed-price scripts import mutable repository `production.inputs.load_inputs` and `DEFAULT_PRICE`, then call `production.equilibrium.solve_at_price`. The relevant session displays baseline exact comparisons to chain13; the baseline full native report retains 31parameter rows,14target rows,87native arrays, frozen observer identities and source snapshots. The baseline parameter table reports psi_child=.17892072066041628, financed_share=.8, annual real rate.02 and H0=6.40569359569417. Its source snapshot is under `ge_runs/ge_default/runtime_auth/runtime_preparation/native_preparation/source_snapshot/`.

Fixed-price scripts did **not** serialize their complete final P object, grid, all runtime source pins or paired policy/distribution arrays. Their scripted explicit overrides and summary values are recovered exactly, but the complete executed runtime cannot be authenticated from source recovery alone. Do not label them fully reproducible certified cases. The snapshot in the baseline native GE is evidence for that case, not proof that every earlier/later import used identical bytes.

| Family | Explicit experimental changes beyond named chain13 | Measurement/acceptance |
|---|---|---|
| credit_space_2x2 | phi.8/.95; renter cap6/4; optionally owner menu4/6/8/10 replacing2/4/6/8/10 | Fixed price; births are explicit flows; age25 uses.125pre+.875post; ownership uses realized post-tenure `g`; no native report gates. |
| credit_rooms_per_child | phi.8/.95; cap6/4; child-room floor0/.5/1 | Same timing; summary mass/NaN/censor diagnostics only. |
| credit_renormalized_psi | child-room floor1; cap4/6; psi re-normalized.3222148912/.2977576805; phi.8/.95 | Births restored near baseline at fixed price, not calibrated target fit or accepted GE. |
| credit_earnings_penalty | child earnings penalty in working ages plus phi/cap/child rooms | Script recovered; result JSON absent in original scratch; no numerical claim can be recovered from that filename alone. |
| would_be_parents | phi.99; rental cap12; both | This tests99% financing, not exact elimination of down payments. Reference childless-renter weights reconstructed from post-fertility mass. |
| unsecured_credit | renter limit0/.37/.74/1.48 in model wealth units | Fixed price; completed stock measured at j8; no gates. |
| gradient_fixes | transfer floorG0+Gn*m; kappa scaling | Transfer policies unfinanced; fixed-price experimental preference/transfer changes. |
| gradient_fixes_renorm;credit_in_benefit_world | G0=.75,Gn=.3/.6 or kappa x2; psi re-normalized; benefit-worldG0=.75,Gn=.6,psi=.0618127463; phi.95/.99 or unsecured.74 or cap12 | Bundled changes, not only a transfer switch; fixed price and no full calibration. |
| ge_strongest_case | cap4,child-room1,psi.3222148912 with phi.8/.95 | Both modified arms fail native dated-budget gates. |
| ge_diagnostic_root | same modified contract plus default phi.8/.95 | Explicit audit-gate bypass, diagnostic secant root, no production adoption. |
| property_tax_ss | tau_H.0423934430954904→.0847868861909808 per4year period; optional extra-revenue lump-sum | Annual1.059836%→2.119672%; diagnostic secant/fiscal iteration, supplier net accounting is separate arithmetic. |

`fertility_by_income.py` reads the cached CPSJune2024 file at `code/data/cps_fertility/cache/jun24pub.csv`, selects women, children0..5 and positivePWSSWGT, caps children at3, and selects ages24–26/40–44. Income categories areHEFAMINC<=10,11–14,>=15, described as<$40k/$40–100k/>=$100k. Own-family young results usePRFAMREL1/2. Model results use `g_beginning_distribution` post-fertility at j1/j6 and grouped current persistent Markov earnings states (not permanent types) by the midpoint of cumulative state mass, giving36.3%/27.3%/36.3%. These are not matched samples/income concepts; the measurement worker owns the independent check. Original transcript line836 includes NumPy floating-point warnings; line846 asserts an exact-sum cross-check, but no separate cross-check script/output was recovered.

`birth_gap_decomposition.py` reconstructs first-birth value gaps from logit probabilities and converts value differences into wealth units with a forward wealth derivative. It weights occupied childless inherited renters using reconstructed reference pre-fertility mass, filters slopes/probabilities, and reports ages18–21/22–25/26–29 and income states4/5/6. The “.02” and “50times” statements summarize selected state comparisons, not a recovered single population statistic. No paired underlying arrays were saved.

## Confidence and remaining limits

High confidence that original scripts/results and native gate failures are preserved byte-for-byte: all194original files passed SHA checks. High confidence in reported arithmetic from saved summaries, without fresh model solves. Incomplete confidence in complete fixed-price runtime identity and validity because those cases retained no serialized inputs, source pins, policy arrays or native acceptance packet. Original output is evidence to review, not an accepted model contract. The recovery stop condition is satisfied; no broader archive search is warranted.
