# Verified recent-parent diagnostic and prepared follow-up batches

Job **17362130 passed** in 30 scheduler seconds (27.324 driver seconds), using 916,792 KB maximum resident memory. All 13 tests passed with compiled transport forced in the startup-test child. The two full-grid observer files are byte-identical and their canonical result hashes agree. There were zero Bellman or equilibrium solves. Independent read-only verification rehashed all 641 archived source files and the original baseline checkpoint; no mismatch occurred.

The baseline selected-birth ownership diagnostic is **−0.05324358594536549**, or −5.324 percentage points. The unchanged ACS empirical reference is +0.16289550916123285, or +16.290 points. This is an opposite-signed descriptive comparison under the declared synchronized snapshot, age-exposure and residence proxy. The exact empirical model value remains **null**, as do the actual weight and loss contribution. This observation does not resolve annual ACS oldest-child or resident-adult-child alignment, and no policy inference is made.

| Group | Owner numerator | Household denominator | Ownership |
|---|---:|---:|---:|
| Actual current birth from an empty-dependent home | 0.0153983163943 | 0.0308594476346 | 49.898224% |
| All currently empty-dependent homes | 0.0933447301889 | 0.169033618710 | 55.222583% |
| First-birth contribution | 0.00616273141557 | 0.0137950933593 | 44.673358% |
| Continuation-birth contribution | 0.00923558497868 | 0.0170643542752 | 54.122089% |
| Empty never-parent contribution | 0.0341118940326 | 0.0761868930613 | 44.773967% |
| Empty former-parent contribution | 0.0592328361564 | 0.0928467256489 | 63.796365% |

The complete six groups, annual-age vectors, accounting checks and approximation warnings are preserved in `output/baseline_pass_01.json`, its identical second pass, `groups.csv`, and `verification.json`. The largest recorded transport/mass error is approximately 2.23e-13, below the unchanged 2e-10 observer tolerance. Original archived source, input and empirical-definition identities are preserved in `executed_bundle/contract.json`; the collected files have independent remote/local hashes in `transport_manifest.json`.

Verified original output hashes:

- Summary: `4dc8a76f0f27bf96ddbe0218610ac023fd41438fe7f973ba0389c0af5466a4ff`.
- Receipt: `a1299de0b9a645fd5ab6c946672ba3dc2ff3152c3ea24b1285466cd0ca44976a`.
- Each observation file: `7696fe049e1a7603ff8a262c3a70f16386b480108d1a76885af4417faae9046f`.

## Prepared only: remaining checkpoints

`prepared_batches/batch_manifest.json` indexes four bundles. Three bundles (`panel_01`, `panel_02`, `panel_03`) contain six remaining sensitivity checkpoints each, with one observer pass per case. They use the **unchanged successfully executed reader**, the same 641-file observer source manifest from `70abd4a8`, and every original case's 634-file `7e872053` policy-source manifest, checkpoint, parent contract, output contract, collection receipt and numerical gates. Their case identities have not been rewritten.

The fourth bundle, `joint_smoke`, selects repetition 02 of original joint job **17360699**, checkpoint SHA `2b99ac5195f3c22de96350397c196b4d32d4f09fe7058292a3dfebf076064116`. Its external reader has a narrowly adapted parent-validation function for the native **two-repetition joint smoke**. It reuses the exact original panel validator and checks both repetitions, original operator/fiscal/market gates, exact early/price/normalization reproduction, and original checkpoint claims. It retains `joint_smoke_analysis_47` and original source `7e872053`; it does not label the joint smoke as a sensitivity-panel observation. The helper's hash matches the original approved joint run plan. Existing graph hashes are retained in the verified original receipt; graph bytes are not re-read by this observation task.

All four bundles pass local archive/source, light-input and original-gate preflight, Python syntax and shell syntax. Their checkpoints will be rehashed remotely before loading. **No follow-up batch has been submitted or observed.** No result is extrapolated for the other 18 points or the joint candidate.

Each script requests one CPU, 8 GB, five Slurm minutes and automatic partition selection, with a 220-second work budget, 20-second reporting reserve and 240-second hard stop. For each six-case batch, assigning the entire measured two-pass baseline runtime to every single-pass checkpoint gives a conservative planning allowance of 163.95 seconds. This fits the existing limit with about 56 seconds of work-budget headroom; it is an estimate rather than a runtime guarantee. The same compiled startup test is retained for every batch.

Stage each bundle directory as `batches/<bundle>/` beneath `/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8`, reusing the already verified archive. The submission scripts carry their exact new contract hashes. Lead review and submission remain separate steps. The unchanged reader reports `exact_reproduction_passed=false` for one-pass batches because it did not request a second pass; this is **not a failed comparison**. Baseline two-pass reproducibility is independently verified and pinned in every bundle.
