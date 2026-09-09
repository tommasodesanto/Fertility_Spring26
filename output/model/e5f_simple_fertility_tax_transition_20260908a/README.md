# Rebated property-tax paths: completed morning packet

**September 9 clarification:** this packet does not complete the author's intended perfect-foresight policy exercise. It uses temporary equilibria and the household birth queue. Earlier verified perfect-foresight/person-demographic work remains preserved but was not integrated with the new fertility-nest specification for these runs. Resolve that handoff before treating this packet as the final quantitative transition.

Both paths completed all11dates2023–2063. Annual1% versus2%, each rebating its own revenue equally per household. Selected simultaneous fertility-nest calibration fixed at loss23.791955; all12targets and11estimated parameters unchanged. No production promotion.

## Read first

- September 9 household/child-count review: `households_children_review/OVERVIEW.md`, with verified policy counts, fresh ACS household-share validation, child-state/entry accounting and unresolved empirical group mappings. This qualifies demographic and family-group interpretation without changing the numerical experiment.
- Morning PDF: `../../pdf/rebated_property_tax_morning_20260909.pdf` (27pages:8readout/calibration pages,19standard diagnostic appendix pages).
- Complete374standard PNGs and receipts: `../../pdf/rebated_tax_complete_graphs_20260909.zip`.
- `morning_packet/READOUT_VERIFIED.md`: all11dated effects.
- `morning_packet/report_manifest.json`: complete report content, target/parameter tables, source/image hashes.
- `morning_packet/verification.json`: independent numerical/result checks.
- `morning_packet/graph_verification.json`: all22dated graph packets and374images verified.
- `results/`: complete collected tables/receipts;52large checkpoints remain onTorch.
- `graphs/v2/full/`: complete final standard graphs. Original `graphs/smoke/` is superseded and must not be used.

## Main result

By2063, births per household rise1.049670%, total births1.251306%, young ownership2.154084percentage points, and young-parent ownership1.745674points. Young rooms fall4.497144%; young-parent rooms fall4.183808%. Thus the fertility gain persists while occupied housing space contracts. Initial purchase prices fall4.801637%, yet ongoing rent/user cost rises18.177669%. The verified impact decomposition is in `../e5f_simple_fertility_tax_channels_20260908a/results/`.

Young uses age nodes26,30,34; exact annual-age ACS alignment is unresolved. Parent groups have dependent children and their means include composition effects. Birth flows differ from the cohort completed-fertility calibration moment.

## Scientific contract and outstanding items

Same verified selected2023inherited population,2019post-advance entryqueue and dated supply elasticity0.63, without reanchoring. Closed national household-unit diagnostic: no outside entries, full retention, adjusted births/2.1, twenty-year lag. No future Census reweighting, person-demographic replacement, perfect-foresight transition or welfare calculation. Prices are perceived permanent at each temporary-equilibrium date. Household formation/headship and the historical-to-endogenous-entry handoff remain unresolved.

Calibration still misses first-birth housing:0.443859versus0.720246rooms. Late-life ownership is almost universal. `income_audit/` independently confirms positive income–housing gradients. `owner_shape_audit/` traces ownership dips to rental/four-room value crossings, with six-room ownership notyet attractive. All inspected origin states havezero mass; grid sensitivity remains untested. The extraction wrote verified CSVs but metadata serialization failed; see its explicit partial receipt. No numerical problem was repaired by relabeling graphs.

Frozen contractSHA256: `ce80b6ad241bec9556d2f3d3cdccb9f89e0ccfdd5b2e80ddcd167ed5156b5b68`.
Scientific bundleSHA256: `4199e948c5f3625c4a2af106623344ddd8f0b032262f26a8d3973223f5bd63c8`.
Targetfingerprint: `3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`.

## Execution and validation

Smoke17250629 completed two exact-loop dates per policy; full17250630 completed eleven dates per policy; collector17250631 completed. Total26coupled roots plus26fresh fixed-price replays. Each branch2CPUs/48GiB with one numericalthread;4hour stage cap and30minutes per date. Checkpoints/queues/latest summaries per date and heartbeat30seconds. Numerical failures would stop the run without retry or relaxed tolerance. Both paths completed within the planned35–45minute compute envelope, with queue time additional.

Independent verification checked26dated packets,324collected artifact hashes and1,639conditions. All11effect rows recompute exactly. Max market residual1.50843e-5 (gate2e-4); fiscal imbalance2.42821e-5 (gate2.5e-5); mass error1.11022e-15. Occupied budget-violation masszero. Both2023impact equilibria and2023/27smoke paths reproduced. Initial states/queues agree; policy household mass cannot differ before2043 and passes that gate.

Graph smoke17251384 and full17251389 completed. Each graph job1CPU/24GiB/20minutes, max4parallel, zero numericalsolves. The obsolete pending17250948 was cancelled after reporting QA. Reporting-only adapter corrects all-owner-product/conception probabilities and pre-choice first-birth risk weights; preserves17standard filenames; verifies native operator identities and unchanged inputs. Original scientific code and numerical manifests stay frozen. Source2803e5fa on `codex/fertility-nest-computation`; PDF builderd122eb52.

## Reproduce and inspect

The isolated worktree is `tmp/e5f_fertility_nest_compute_20260907a` under projectroot. PDF renderer: `code/model/tools/build_e5f_rebated_tax_morning_report.py --manifest <absolute morning_packet/report_manifest.json> --output <new absolute PDF path>`. Use bundled document Python; the manifest pins all content/images. Supplemental source recipe is retained in `morning_packet/build_manifest_and_supplemental.py`; graph verifier/packager in `morning_packet/verify_and_package_graphs.py`. The committed report manifest is authoritative for final captions and formatting changes.

Complete graph regeneration from saved states uses `code/model/tools/run_e5f_simple_fertility_tax_graphs.py` and the recorded `graphs_v2_full.sbatch`, under the frozen Torch snapshot `/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a`. Use a fresh graph output root; never write plots inside numerical run folders. All scientific/target contract hashes and source pins must pass before reads. No new calibration is required to regenerate diagnostics.

## Monitoring completion

The author explicitly requested a goal, overnight monitoring, graphs and sleep prevention. Heartbeat `overnight-rebated-tax-results-and-graphs` is paused because the completed PDF is ready for delivery. CaffeinatePID64940 prevents idle/system sleep for12hours from its recorded start; it expires automatically. Keep power connected and lidopen for the assertion to serve its purpose. No additional computation is needed merely to fill the night.
