# Handoff: implement the E5F child-room-floor arm on the repaired E5 architecture

Use this entire document as the task. Scope is implementation, tests, and
launcher preparation only. **Do not submit anything to the cluster.** The lead
reviews the diff line-by-line against this specification before any
submission. Hand back when the deliverables at the bottom are ready.

Repository: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`
Package: `code/model/intergen_eqscale_seq_optimized/`
Base: commit `7a402cc` (independent child maturation + funded policy closure).

## Do not touch

- Torch jobs `15398062` (array) and `15398070` (collector) are the running
  control arm. Do not cancel, modify, or write into their output directories.
- The maintained E5b defaults, the twelve targets and weights in
  `e5_profile.py`, the maturation-repair behavior, and the policy driver
  `tools/run_e5_repaired_policy_with_entry.py` are all out of scope tonight.
- Everything you add must be default-off behind an explicit switch, with a
  bitwise-nesting guarantee (below).

## Economic specification

The arm replaces the child expenditure-share tilt with a child room
requirement, on the repaired (n, m) architecture, where n = children ever
born (parity axis) and m = children currently at home (the repaired
child-state axis, `child_state_mode="independent_count"`).

Under `preference_spec="eqscale"` the current code zeroes the Stone-Geary
floors and applies the share tilt (solver.py:2158-2196, specifically
c_bar/h_bar zeroed at 2176-2177 and the alpha tilt at 2178-2180). The E5F arm
changes exactly one block:

1. **Housing floor, keyed off m only:**
   \[ \bar h(m) = \bar h_{child} \cdot m \]
   populated into the existing `h_bar[nn, cs]` slots (the Stone-Geary
   machinery already consumed by the kernels), where the cs axis is the
   at-home count m under the repair. Adults and childless states carry no
   floor: \(\bar h(0) = 0\). No first-birth jump term. `c_bar` stays zero
   everywhere.
2. **Tilt retired:** with the floor active, `delta_alpha = 0` and
   `delta_alpha_jump = 0`, and both leave the search domain.
   `alpha_cons = 0.733` stays external, so \(\alpha(n,m) \equiv \alpha_0\).
3. **Unchanged:** the SSK equivalence scale (`eqscale_form="power"`), the
   child utility flow `psi_child`, both fertility logit scales, the fecundity
   schedule, all finance objects (phi = 0.80, down-payment tests), the
   maturation repair, and the twelve-target system with its weights.
4. **Keying report (do not change, just document):** state precisely which of
   the existing child-linked preference objects (tilt, escale, psi_v) key off
   parity n versus at-home count m in the post-`7a402cc` code, with file:line
   references. If any keying looks inconsistent with "needs follow the
   children at home, joy follows children ever born," flag it in the report;
   do not silently change it.

Utility in the floor arm, schematically (sigma = 2, renter case):
\[
u = e(\cdot)\,\frac{\left[c^{\alpha_0}\,(h-\bar h(m))^{1-\alpha_0}\right]^{1-\sigma}}{1-\sigma} + \psi\,n
\]
with the owner case applying the existing chi-on-residual-services convention
to \(H - \bar h(m)\), exactly as the dormant Stone-Geary path already does.

**Feasibility by construction (assert, don't hope):**

- Domain for the new parameter: `hbar_child_rooms` in **[0.10, 1.80], log
  transform**. Upper bound chosen so the maximum floor (3 at-home children
  under L4 top-coding: 3 x 1.80 = 5.4) stays strictly below the renter cap
  `hR_max = 6.0` — every family can always rent legally sized housing, so no
  state is ever infeasible through the floor. Assert
  `n_child_max * hbar_child_rooms < hR_max` at setup and refuse to run
  otherwise.
- Owner rungs with \(H - \bar h(m) \le 0\) (e.g. the 2-room rung for a family
  with two at-home children) must resolve through the existing
  net-service floor / infeasible-branch handling: deeply negative value, no
  NaN, no positive stationary mass. Add a test.

## Identification

Nine free parameters (ten minus the two tilt coordinates plus
`hbar_child_rooms`) against the unchanged twelve targets — overidentified by
count. The two rooms moments (`housing_increment_0to1`,
`prime30_55_parent_3plus_minus_1to2_mean_rooms`) switch from identifying the
tilt to identifying the floor level and its per-child linearity. State this in
the arm metadata. Do not touch targets or weights.

## Implementation pattern

Follow the existing conventions exactly:

- New module `e5f_floor_profile.py` mirroring `e5_maturation_repair.py`:
  overrides dict (`child_room_floor: True`, `hbar_child_rooms` default 0.0 =
  off), metadata dict, and the E5F domain (the E5 domain with the two tilt
  rows removed and the `hbar_child_rooms` row added).
- Env-gate wiring in `run_e1_chain.py` in the style of the `E3_L4`/`E5`
  gates: `E5F=1` selects the floor arm on top of the repaired-E5 chain
  configuration.
- Preference population in the `precompute_shared` eqscale branch
  (solver.py:2158-2196): when the floor switch is on, fill
  `h_bar[nn, cs] = hbar_child_rooms * m(cs)` instead of zero, keep
  `c_bar = 0`, force the tilt contributions to zero. Touch nothing on the
  default path.
- **Bitwise-nesting guarantee:** with the switch off (or
  `hbar_child_rooms = 0`), solver outputs must be bit-identical to the
  post-repair path. This is the same guarantee pattern as
  `fecundity_omega1 == 0` (parameters.py:742-748) and
  `kappa_fert_continuation = None`.

## Tests (add to `tests/`, all must pass alongside the full existing suite)

a. Nesting: floor switch off and floor = 0 both reproduce the repaired-E5
   solve bit-for-bit (V, policies, stationary distribution, all 12 moments).
b. Floor mechanics: at a state with m = 1, renter utility is computed on
   (h - hbar_child_rooms); marginal utility diverges as h approaches the
   floor from above; h at or below the floor is never chosen with positive
   probability.
c. m-keying: after a child matures (m falls by one at fixed n), the floor
   falls by exactly hbar_child_rooms; a household with n = 2, m = 0 carries
   zero floor.
d. Owner short rung: a rung with H - floor <= 0 gets the infeasible-grade
   value, no NaN anywhere, zero stationary mass on it.
e. Domain guard: setup refuses hbar_child_rooms * m_max >= hR_max.
f. Adults: childless and empty-nest states carry zero floor at every age.

## Smokes and launchers (prepare, do not submit)

1. Local exact-loop smoke in the established pattern: a handful of theta
   draws through the E5F chain path at the strict evaluator (J=17, Nb=120,
   max_iter_eq=40, tol_eq=2.5e-5), each solved twice, losses bit-identical
   across repeats. Include one draw at the domain upper bound to exercise the
   short-rung handling. Report per-solve wall time.
2. Launcher pair mirroring `code/cluster/submit_intergen_e5_maturation_repair.sh`
   and its collector: 8 chains, same walltime and account, distinct run tag
   (suggest `intergen_e5f_child_room_floor_20260805`), plus a 2-chain torch
   smoke script in the house pattern. Write them; do not `sbatch` anything.

## Deliverables (hand back to the lead)

1. Diff map: every changed/added file with line references, one sentence per
   hunk.
2. The keying report from item 4 of the specification.
3. Any place the eqscale branch assumed `h_bar == 0` and how you handled it.
4. Full test-suite output (existing + new), local smoke log with timings.
5. The prepared launcher paths.
6. Nothing submitted to Torch; no edits outside the package, its tests, and
   `code/cluster/` launcher files; no target, weight, default, or
   status-file changes.

The lead will verify the preference-block diff line-by-line against the
specification above before authorizing the cluster smoke and the overnight
array.
