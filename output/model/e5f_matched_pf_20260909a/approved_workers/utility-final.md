# Parenthood utility adapter: completed bounded implementation

Owned source files, in the isolated `codex/balanced-social-security` worktree:

- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/e5f_parenthood_utility.py`
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/test_e5f_parenthood_utility.py`

No core solver, model production checkout, manuscript, other worker file, launch, commit or push was changed by this worker. The recovered CLI-worker draft was inspected after the lead confirmed its process had stopped; this implementation replaces its invalid h_P override and incomplete tests.

## Integration API

1. `initialize_parenthood_utility(Pold)` explicitly performs the one-time migration on a deep copy: `hbar_first_child_jump = old jump + old slope`; `hbar_child_rooms = 0.0` exactly. It binds sigma 2, constant alpha (externally fixed at 0.733), power equivalence scale, and zero consumption floors. It requires the retained four-year, independent-child-count lifecycle and does not rebuild it. All other primitive fields and arrays are preserved.
2. `validate_parenthood_utility(P)` is the launch/reload gate. It rejects any nonzero slope, including 1e-300, changed fixed utility/lifecycle fields, nonfinite inputs and a floor at or above the rental-service cap. **Do not use the explicit legacy migration as a reload validator.**
3. `bind_parenthood_utility(P, candidate=None, *, copy_parameters=True)` accepts an already validated new-utility P and a partial or full physical-coordinate dictionary. It returns a deep copy by default. It translates `h_P` to the existing jump field and `beta_annual` to `beta_annual**4`; beta aliases rho/rho_hat and the fertility-scale alias eps_fert update coherently. It never sends an `h_P` field to the old solver and never calls the general `apply_overrides`, which also rebuilds fiscal income. In-place mode is explicit and validates the whole proposed change before mutation.
4. `validate_parenthood_candidate(candidate, require_complete=True)` requires exactly the nine structural coordinates. Without that flag, partial dictionaries support bounded coordinate probes. Either mode rejects the old slope/jump names, psi_child and all other fixed/unknown fields; numerical bounds are checked.
5. `PARENTHOOD_SEARCH_DOMAIN`, `PARENTHOOD_SEARCH_NAMES` and `parenthood_utility_metadata()` supply the domain and serializable launch restrictions. The eight unchanged bounds/transforms come directly from the existing income-entry domain. The new h_P coordinate has approved numerical bounds [0.1, 2.3] and the previous positive floor's log transform. Existing transform functions and old domains are untouched.
6. `parenthood_utility_overrides(P, candidate)` provides just the translated field mapping when needed; `initial_parenthood_requirement(Pold)` computes only the legacy first-child mapping.

The nine coordinates are beta_annual, kappa_fert, kappa_fert_continuation, chi, H0, theta0, theta1, first_birth_fixed_cost, h_P. psi_child remains outside this structural vector for the separate normalization. Binding deliberately preserves its current value so the frozen-psi comparison is possible. The adapter is an explicit utility contract; target fingerprints, fiscal balancing, supply rebasing and normalizer wiring remain the lead's integration scope.

## Evidence and checks

From the isolated worktree:

```sh
NUMBA_DISABLE_JIT=1 PYTHONPATH=code/model:code/model/tools /usr/bin/python3 -m unittest test_e5f_parenthood_utility -v
```

**11 tests passed in 0.069 seconds** (0.25-second command runtime). Both new files also passed Python AST syntax parsing. Scoped status confirms only these two new utility files are owned here; the workspace diff whitespace check passed.

Tests use the real setup_parameters object, configure_child_state_process, precompute_shared, and real full_renter_block_kernel/full_owner_block_kernel in Python mode. No utility/type/kernel function is mocked. They check:

- all feasible lifetime-parity/current-dependent-child combinations m=0,1,2,3, with psi_child=0 and 0.19;
- exactly equal parent floors while per-state equivalence scales remain distinct;
- actual maturation probability to m=0 and disappearance of floor, child reward and scale increment;
- preservation of old first-child floors and every unrelated field/array, including deliberately externally bound fiscal income and pensions;
- annual beta raised to four, its aliases, candidate H0 shape, no original-object mutation;
- strict rejection of stale candidate keys, even microscopic slope changes, fixed-field corruption, missing full coordinates and bound violations;
- real renter choices satisfy s=(1-alpha)X/r+alpha*h_P when interior, actual CRRA values equal -e/Q+psi*m, rental-cap behavior, exact budgets and the binding zero-saving boundary;
- floor-affordability rejection for renters and strict owner housing-floor rejection, with real feasible owner CRRA values, budgets and retained ownership service premium.

The existing solver's legacy `type_map` collapses parent floor/reward triples at psi=0, but its full state-specific `escale_flat` does not collapse: the tests exhibit this exact state. Source inspection of solve_bellman_full_markov_income confirms that the actual Markov path passes full cb/hb/psi/alpha/escale arrays to the kernels, rather than rebuilding scale from that legacy triple map. No solver change is needed for this utility shape.

## Remaining verification

No compiled lifecycle, stationary equilibrium, normalization, transition or fiscal loop ran. The local Python is 3.9 with NumPy 2.0.2 and no Numba. Its generic old apply_overrides helper contains a Python-3.10-only union isinstance expression, so the fixture sets scalar primitives on the actual setup P and calls the real child-state constructor directly; no runtime or source workaround was applied. The same test module can be run with enabled Numba on the cluster for compiled tiny-block evidence; a full compiled household/initial/fiscal loop and launch/reload integration remain necessary before production use. Existing numerical and empirical gates are unchanged.
