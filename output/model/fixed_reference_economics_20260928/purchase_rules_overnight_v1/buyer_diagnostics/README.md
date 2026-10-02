# Saved-solution buyer diagnostics

`summarize_saved_buyers.py` reads the native selected-point `stage/solution_arrays.npz` on Torch. It does not call a Bellman or equilibrium solver. Supply the owner unit-size ladder and sale-cost fraction from the authenticated run; the script checks array dimensions and tenure probability sums. Example:

```bash
/share/apps/anaconda3/2025.06/bin/python summarize_saved_buyers.py \
  --stage /scratch/td2248/projects/strict_purchase_sandbox_v1/results/onecase/run/strict_purchase/phase_b_ge/selected_repeat/stage/solution_arrays.npz \
  --out /scratch/td2248/projects/strict_purchase_sandbox_v1/buyer_diagnostics.json \
  --owner-rungs 2,4,6,8,10 --sale-cost .06
```

The buyer population is the post-fertility, pre-tenure household distribution multiplied by its destination-tenure probability. Renter-to-owner and owner-to-different-rung flows are reported separately; same-rung owner stayers are excluded. “Renter-to-owner” does not establish lifetime first ownership. The model state ages 26, 30 and 34 implement the existing 25–34 young-owner measurement window on its four-year grid. The reported closing ratio is \(\max(0,Q-A)/Q\), where \(Q=pH\) and \(A=b+(1-\psi)pH_{old}\) for an owner origin. Since the model records net liquid assets rather than a separate gross mortgage and cash account, this is an **implied net funding ratio**, not observed mortgage LTV.

The [strict-80 fixed-coordinate smoke](strict80_fixed_coordinate_smoke.json) ran successfully on the saved native repeat arrays, with no solves. Among young renter-to-owner buyers its implied net ratio has median 0.6545 and 90th percentile 0.7887; the mass-weighted share above 80% is zero. This is an old fixed-coordinate diagnostic, not either overnight recalibrated baseline.

Additional old fixed-coordinate saved-array checks at \(\phi=1\): [hard-100](hard100_fixed_coordinate.json) has zero buyer flow with negative closing cash \(A<0\) at all ages. [Quarter-100](quarter100_fixed_coordinate.json) has 0.0013183 owner-switcher flow mass with \(A<0\), or 1.146% of its owner-switcher buyer mass; renter-to-owner negative-\(A\) flow is zero. At young ages 26, 30 and 34, quarter-100 has 0.0006930 negative-\(A\) owner-switcher flow. This supports a choice-set difference for negative-\(A\) switchers; it does not isolate its contribution to fertility or ownership because the two saved solutions have different distributions.

For a first-birth denominator, provide an exact `first_birth_flow` array in an NPZ with the same dimensions as the saved distribution. The overnight selected stage now retains `distribution.g_pre`; pass that same NPZ as `--pre-birth` and pass the exact `get_fecundity_by_age(P)` vector as `--fecundity-by-age`. The script validates that reconstructed first-birth mass is contained in the post-birth distribution. A synthetic branch-flow and exclusion check passed on Torch. With `--hard-phi .8` or `1`, it reports the share of origin-renter first-birth flow that fails the closing cash screen for every owner rung. This is a **down-payment exclusion** measure. It does not imply that a household that rents wanted to buy, nor does passing the cash screen guarantee feasible consumption or an ending estate.

The quarter-saving buyer restriction depends on next-period saving \(b'-x\), so the saved starting cash and observed tenure probability do not by themselves identify who was unable to buy. A quarter-rule exclusion report requires branch-specific feasible owner choices from the forward-only replay or the Bellman candidate feasibility map. Do not report all new-parent renters as constrained.

`financial_access.py::matched_first_birth_access` is a postcheck hook for that financial-access calculation. It accepts the live native parameter object, shared precomputation, prices, pre-birth distribution and fertility probabilities. It computes the same origin-renter first-birth state weights for \(\phi=.8\) and \(1\), and asks whether any owner product passes the rule-specific purchase screen, transaction-grid support, physical housing floor, saving/death/grid lower bound, and budget with positive consumption surplus. It uses the engine's income, family earnings adjustment, housing-stage context, and fecundity helpers. It reports the stricter \(c_{min}\) margin separately because the owner's numerical utility kernel itself rejects only surplus \(\leq10^{-10}\). It does not evaluate continuation values or choice preferences, so its label is financial access, not desired ownership or caused births. A zero-solve execution with the authenticated older input bundle passed on Torch and found no reverse access loss in the inspected age/income cell; the hook has not yet run on an overnight selected postcheck. Review that execution before numerical reporting.

Once a chain's selected postcheck has passed, `run_selected.py --arm hard|quarter --completed <chain/postcheck/completed.json> --out <new-directory>` authenticates the fit through `mechanism/selected_runtime.py`, reconstructs its exact selected parameters, and writes `buyer_net_closing_ratio.json` and `matched_first_birth_financial_access.json` from its saved repeat stage. It refuses a nonterminal or mismatched selected point and flags any observed first-birth owner choice falling outside the financial map. `run_dated.py --engine-root <matching isolated engine> --packet <date_NNN/diagnostic_packet.pkl.gz> --rule hard|quarter --observed-phi .8|1 --out <new JSON>` applies the same matched-state financial map to a trusted dated policy packet. The dated mechanism driver's `mapping.json` supplies the actual first-birth response; these files supply the access mechanism diagnostic. Neither helper launches a model solve. They must be staged separately from the immutable calibration source archive and smoke-tested on the first genuine selected/dating packet.

## Separate Torch addendum and exact later commands

The original separate addendum root is `/scratch/td2248/projects/purchase_buyer_diagnostics_v1`; preserve it unchanged. The selection-snapshot-compatible addendum stages under `/scratch/td2248/projects/purchase_buyer_diagnostics_v2`. `preflight_torch.sh` checks the SHA-256 of every staged buyer, mechanism, floor and base source, then runs the native financial map with the authenticated input bundle for both purchase rules and **zero model solves**. It does not submit jobs. From the repository root:

```sh
bash output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/stage_torch.sh
ssh -o BatchMode=yes torch 'bash /scratch/td2248/projects/purchase_buyer_diagnostics_v2/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/preflight_torch.sh'
```

After `purchase_mechanism_v1/selection/manifest.json` points to two actually postchecked winners, submit one saved-policy readout per arm. The wrapper mounts the immutable `purchase_mechanism_v1/selected_postchecks/` snapshot, keeps each winner's original physical source and restart-parent provenance, and verifies the chosen JSON, completed receipt, full report and native-array SHA-256 values before entering the native runtime. It then requires both exact ROOT/REPEAT reports, 14 targets, 31 parameters, 17 plots and the retained repeat arrays. It writes only under its new result directory.

```sh
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=hard /scratch/td2248/projects/purchase_buyer_diagnostics_v2/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=quarter /scratch/td2248/projects/purchase_buyer_diagnostics_v2/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
```

After a dated mechanism case has `run/completed.json` with `status=passed`, the same readout wrapper resolves its accepted mapping, checks the selected date packet exists, and reads the model's observed financed share for that date. For example, the date-zero control packet for the hard rule:

```sh
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=dated,BUYER_ARM=hard,BUYER_CASE=case_00_hard_control_h12,BUYER_DATE=date_000 /scratch/td2248/projects/purchase_buyer_diagnostics_v2/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
```

Substitute the actual accepted case name and `date_000`, midpoint or final date. The wrapper refuses missing or failed case receipts. It runs one CPU, 12 GiB, one thread and 20 minutes at most. It never starts a Bellman, equilibrium or transition solve. Existing calibration and mechanism inventories are read-only and remain unchanged.

## Reviewed-runtime v3 addendum

The current source-only addendum is `/scratch/td2248/projects/purchase_buyer_diagnostics_v3`, with archive `buyer_diagnostics_stage_v3.tar.gz`. It authenticates the reviewed mechanism source under `/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a` (including the compact repeat, import and renderer fixes). It reads the eventual immutable winner JSONs and copied postchecks from the separate `/scratch/td2248/projects/purchase_mechanism_v1/selection` and `/selected_postchecks` data store. The prior buyer v1/v2 archives and results remain intact.

From the repository root, stage and run the zero-solve v3 preflight:

```sh
bash output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/stage_torch.sh
ssh -o BatchMode=yes torch 'bash /scratch/td2248/projects/purchase_buyer_diagnostics_v3/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/preflight_torch.sh'
```

Only after the lead publishes final postchecked winners and accepts the policy source, the exact selected readout commands are:

```sh
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=hard /scratch/td2248/projects/purchase_buyer_diagnostics_v3/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=quarter /scratch/td2248/projects/purchase_buyer_diagnostics_v3/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
```

## Corrected v4 saved-policy readout

The v3 selected jobs stopped at a plot-location check before writing diagnostics: the native `selected_root` report has all 17 standard plots, while `selected_repeat` saves a compact closure and state arrays without a second plot set. `run_selected.py` now checks the 17 plots in `selected_root`; `selected_runtime` checks the compact closure against the full root and authenticates the saved arrays through its native replay. Preserve the failed v3 jobs as evidence. The corrected source stages in a separate `/scratch/td2248/projects/purchase_buyer_diagnostics_v4`, retaining the reviewed mechanism source and v1 immutable selection data routing.

```sh
bash output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/stage_torch.sh
ssh -o BatchMode=yes torch 'bash /scratch/td2248/projects/purchase_buyer_diagnostics_v4/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/preflight_torch.sh'
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=hard /scratch/td2248/projects/purchase_buyer_diagnostics_v4/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
ssh -o BatchMode=yes torch 'sbatch --export=ALL,BUYER_MODE=selected,BUYER_ARM=quarter /scratch/td2248/projects/purchase_buyer_diagnostics_v4/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh'
```
