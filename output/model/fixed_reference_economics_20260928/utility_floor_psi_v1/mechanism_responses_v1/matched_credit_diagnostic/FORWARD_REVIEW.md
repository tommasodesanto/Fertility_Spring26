# Saved-policy forward reconstruction for review

No remote forward reconstruction has been executed. The scalar result is
already saved in `AGGREGATE_FLOW_DECOMPOSITION.md` and its JSON.

`forward_flow_decomposition.py` authenticates the original frozen mechanism
runtime and actual 31-parameter binding. It copies the two already validated
credit parameter objects and loads only the saved fields required by native
`policy_from_solution` and deterministic cohort reconstruction. Solver entry
points are replaced by functions that immediately fail. Every dated evaluation
supplies its saved policy. No household, GE or calibration solve is permitted.

First it reconstructs baseline PRE using the native age-forward operator,
requires exact array and hash equality with the retained baseline checkpoint,
requires zero feasibility projection, and reproduces all native birth-order
and total flows. Only then does it reconstruct credit PRE. Both reconstructions
must re-nest into their saved POST distribution at the existing 5e-9 criterion.
Both reproduce saved uniform-clock age40–44 childlessness and topcode-adjusted
birth children divided by entry. Flow comparison tolerance remains 2e-10.

The actual installed birth function is `apply_sequential_fertility`, installed
by `e5f_current_transition_runtime.py:239`. It snapshots the first/subsequent
risk pools before any births, multiplies attempt probabilities by fecundability,
and allows at most one birth per household in a four-year period. Birth order
is measured by cumulative children-count differences across PRE and POST.
To recover contributions by origin children-ever-born and children-at-home
states, the reduction masks each origin family cell and calls this exact
operator. Birth-stage wealth and age do not change, so contributions can then
be aggregated by age and negative/nonnegative inherited wealth. Grouped flows
must sum to each unmasked retained native total.

The three corners are baseline policy/baseline PRE, credit policy/baseline PRE,
and credit policy/credit PRE. Their policy-first identity is exact accounting;
the reverse ordering is not computed. Results retain risk exposure alongside
flows. They use all saved attempt probabilities, including endpoints; the
earlier value-inversion exclusions do not apply. Zero projection is required;
no occupied mass is deliberately removed. The original support-limited native
credit interpretation remains: unoccupied alternatives and grid convergence
are not certified. This does not decompose direct and continuation utility.

`check_forward_flow.py` extracts the exact native birth-function body for a
small fixture, without importing its solver/runtime. It passed fecundability,
no same-period first-to-second chaining, birth-order replay/mismatch rejection,
family/wealth contribution add-up, the three-corner identity and solver traps.
The receipt is `forward_toy_receipt.json`.

Proposed execution is one SSH process using the original launcher's exact
read-only project/source/bundle binds inside the same Apptainer image. Stream
the reviewed script through stdin; bind original results read-only and create
one isolated analysis work directory for runtime authentication/compiler
metadata. Collect compact stdout JSON only. Use one numerical thread and an
external 180-second timeout; stop on failure without retry or extension.

Observed runtime evidence: the successful mechanism launch began at epoch
1790830426; its runtime zero-LC receipt was checked at 1790830458, about32seconds
later. All six native cells and diagnostics completed within205seconds. Saved
metadata reads took under1second and about75MiB. This forward-only pass excludes
household solving and plots;180seconds is plausible, although new forward-kernel
compilation is unmeasured and can cause a bounded failure. Selected policy and
distribution arrays are estimated below250MiB, with native imports/temporary
evaluations bringing working memory roughly below1GiB; this is an estimate,
not a measured forward-run claim. The same16GiB ceiling is ample. No runtime
or memory budget extension is proposed.
