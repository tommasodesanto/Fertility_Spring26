# Session Diary - Spatial Lifecycle Fertility Model

## 2025-12-30

### What we discovered

**Major bug found**: The model was running with `da=6` (6-year periods) even when we set `J=60`, meaning:
- 60 periods × 6 years = 360-year lifecycle (nonsensical)
- All the "3.3 minute" optimization work was on a misconfigured model
- Parameters were calibrated to the wrong model

**Fix implemented**: Added auto-configure block in `run_model_fast.m` that sets annual parameters when `J=60`:
- `da=1`, `J_R=36`, `A_f_end=18`, `A_m=18`
- Annual rates: `rho=0.02`, `q=0.04`, `delta=0.02`, `tau_H=0.01`
- `n_sub=1` (no sub-stepping, following Ben Moll's approach)

### Correct annual model results

| Metric | Value | Target |
|--------|-------|--------|
| Runtime | 27 min | Want < 5 min |
| Mean parity | 2.18 | ~1.8 |
| Ownership | 3.3% | ~65% |
| Fertility gradient | 0.000 | ~0.3-0.5 |
| Childlessness | 0.7% | ~15% |

### What we tested

1. **n_sub=6 vs n_sub=1**: No speedup (both ~27 min)
   - HJB faster with n_sub=1 (500ms vs 1100ms)
   - KFE slower with n_sub=1 (650ms vs 550ms)
   - Cancels out

2. Ben Moll's approach: Uses dt=0.25 (quarterly) with no sub-stepping for lifecycle models

### Key problems to solve

1. **Speed**: 27 min is too slow for calibration loops
2. **Economics are wrong**:
   - Ownership way too low (3% vs 65% target)
   - No spatial fertility gradient (the core mechanism!)
   - Childlessness way too low (0.7% vs 15% target)

### Thinking

The parameters (`h_bar_jump=0.08`, `nu_loc=0.50`, etc.) were calibrated to the broken da=6 model. They don't create meaningful housing cost differentials in the annual model.

Looking at the subsistence costs output:
```
Peripheral: r=0.026, min_expenditure=0.034
Secondary:  r=0.049, min_expenditure=0.037
Superstar:  r=0.174, min_expenditure=0.054
```

The housing cost differential for having children (r × h_bar) is tiny:
- Peripheral: 0.004
- Superstar: 0.024

This ~0.02 difference is not enough to affect fertility decisions.

### Next steps

1. Profile the annual model to find actual bottlenecks
2. Recalibrate parameters for annual model:
   - Increase `h_bar_jump` significantly (maybe 0.3-0.5?)
   - Adjust `nu_fert` to get more childlessness
   - Fix ownership (maybe housing tenure parameters?)
3. Consider if there's a fundamental scaling issue we're missing

---

## 2025-12-31

### What we worked on

**Understanding the slowdown and implementing stochastic aging**

We clarified why the annual model (J=60, da=1) runs in 27 min while the broken config (J=60, da=6) ran in 3 min:

| Config | da | A_m | n_child_states | Total States | Runtime |
|--------|-----|-----|----------------|--------------|---------|
| Broken | 6 | 3 | 5 | 432,000 | ~3 min |
| Correct | 1 | 18 | 20 | 1,728,000 | ~27 min |

**Key insight**: `da` doesn't directly multiply state space. The culprit is `n_child_states = A_m + 2`:
- With da=6: A_m = 18/6 = 3 periods → n_child_states = 5
- With da=1: A_m = 18/1 = 18 periods → n_child_states = 20

The 4× larger child dimension is the primary cause of the slowdown.

### Solution: Stochastic Aging

Instead of tracking exact child age (20 states), we compress to ~6 coarse stages with stochastic transitions:

**Proposed stages (sum to 18 years)**:
1. Empty nest / no kids (absorbing unless fertility)
2. Infant: 2 years (p_advance = 0.5/year)
3. Young child: 4 years (p_advance = 0.25/year)
4. School age: 6 years (p_advance = 0.167/year)
5. Teen: 6 years (p_advance = 0.167/year → maturity)

**Implementation changes**:
1. Build transition matrix `Pi_child` with geometric durations
2. HJB: Replace deterministic index shift with `V_cont = V_next * Pi'`
3. KFE: Push mass with `g_next = g_curr * Pi`
4. BGP: Update maturity flow calculation to use `E[e^{-nT}]` instead of `e^{-n*A_m}`

**Expected speedup**: 20 → 6 child states = 3.3× smaller dimension

### Implementation completed

1. **`make_child_transition_matrix()`**: Creates Pi matrix with geometric durations
2. **HJB modification**: Uses `V_cont = V_next * Pi'` for expected continuation values
3. **KFE modification**: Distributes mass via `Pi` probabilities, counts stochastic maturity
4. **BGP calculation**: Uses `compute_n_from_BE_ratio()` to invert `E/B = E[exp(-n*T)]`

### Test results

**Speed improvement confirmed!**
- HJB: 500ms → 195ms (2.5× faster)
- KFE: 600ms → 250ms (2.4× faster)
- Expected total: ~10-12 min vs 27 min

**Convergence issue found**: Bisection oscillates because:
- At some gamma values, `parity=0.00` (no fertility at all)
- Then jumps discontinuously to `parity=1.54`
- The BGP root-finding struggles with this discontinuity

**Key observation**: At the (non-converged) solution:
- Ownership = 59.3% (vs 3.3% with deterministic - much closer to 65% target!)
- Mean parity = 1.54 (vs 2.18 - actually closer to 1.8 target!)

The parameters `h_bar_jump=0.50, h_bar_n=0.20` from the original calibration seem to work
better with stochastic aging because the coarser states affect value function smoothness.

### Final successful test

Used parameters: `h_bar_jump=0.08`, `h_bar_n=0.03`, `nu_fert=3.0`, `nu_loc=0.50`

**Result: 3.2× speedup achieved! (27 min → 8.4 min)**

| Metric | Deterministic (27 min) | Stochastic (8.4 min) | Target |
|--------|----------------------|---------------------|--------|
| Runtime | 27 min | **8.4 min** | < 10 min |
| n_child_states | 20 | **5** | - |
| Mean parity | 2.18 | 1.45 | ~1.8 |
| Ownership | 3.3% | 10.6% | ~65% |
| Childlessness | 0.7% | 7.6% | ~15% |
| Gradient | 0.000 | -0.001 | ~0.3-0.5 |

**BGP consistency verified**:
- `E[exp(-n×T)] = 0.785` (vs `exp(-n×18) = 0.78`)
- `M_total/E = 1.008` (should be ~1.0)
- Gap: 0.17%

### Issues remaining

1. **Gradient still ~0**: Core mechanism not working - need parameter exploration
2. **Ownership too low**: 10.6% vs 65% target
3. **Parameters need full recalibration** for the correct annual model

### Next steps

1. Explore parameter space to find spatial fertility gradient
2. Focus on `h_bar_jump` and `h_bar_n` - these are the key housing-fertility parameters
3. May need to increase housing cost differences across locations
4. Consider calibration loop now that runtime is reasonable

### Open questions

1. How many child stages are economically meaningful? (infant/child/teen vs more granular)
2. Will stochastic aging affect fertility timing patterns?
3. Should housing requirements differ by child stage?

---

## 2025-12-31 (continued)

### Attempted optimization: Thomas algorithm precomputation

**Goal**: Reduce KFE time from ~250ms/call by precomputing matrix factorizations.

**Approach**:
The KFE calls sparse backslash ~14,000 times per iteration:
- 59 ages × 4 parities × 5 child states × 3 locations × 4 tenures ≈ 14,160 tridiagonal solves

Since the drift coefficients are fixed within the entry loop (same prices, same policies), the matrices are identical across all ~60 entry loop iterations. The idea was to:
1. Precompute Thomas algorithm forward sweep coefficients once per price iteration
2. Reuse them for O(n) back-substitution in each KFE call

**Implementation**:
- `precompute_thomas_coefficients()`: Stores `c_prime`, `inv_denom`, `m_vals`
- `thomas_backsolve()`: Forward sweep on RHS + back-substitution

**Result: FAILED**
- KFE time went from 250ms → 1500ms (6× SLOWER!)
- MATLAB's for-loops in the back-substitution are much slower than the optimized sparse backslash

**Lesson learned**: MATLAB's sparse solvers are highly optimized. Naive loop-based implementations can't compete, even when avoiding matrix construction overhead.

**Reverted**: Thomas algorithm code removed from production path.

### Alternative strategies considered

1. **Reduce number of KFE calls** (more promising)
   - Currently: ~60 entry iterations × ~30 price iterations = ~1800 KFE calls
   - Anderson acceleration for entry loop
   - Better initial guess from previous gamma evaluation

2. **Reduce number of bisection iterations**
   - Currently: ~12 iterations to converge
   - Secant method or Brent's method instead of bisection
   - But F(gamma) can be non-monotonic, making this risky

3. **Parallel processing**
   - The drift solves within KFE are independent
   - Could use parfor for the (i, ten, j, nn, cs) loops
   - MATLAB's parallel overhead might eat gains for small Nb=30

### Current status

~~The 8.4 minute runtime with stochastic aging (from earlier today) is the best we have without loosening tolerances.~~

**UPDATE: Batched KFE drift solves implemented successfully!**

---

### Successful optimization: Batched KFE drift solves

**The fix**: Instead of solving 12 separate sparse systems per (j, nn, cs) iteration, batch all (tenure × location) systems together using the existing `thomas_solve_vec` function.

**Key insight**: The vectorized Thomas algorithm in `thomas_solve_vec` handles multiple systems efficiently through column-wise operations. It avoids:
- 12 sparse matrix constructions per iteration
- 12 separate sparse backslash calls
- Loop overhead

**Implementation**:
```matlab
% Flatten g_post_tenure to (Nb, n_systems) where n_systems = 12
g_flat = reshape(g_post_tenure, Nb, n_systems);
x_flat = reshape(drift_x(:, :, :, j, nn, cs), Nb, n_systems);
y_flat = reshape(drift_y(:, :, :, j, nn, cs), Nb, n_systems);

% Build tridiagonal coefficients
b_diag = 1 + Delta * (x_flat + y_flat);
a_diag = -Delta * x_flat(1:Nb-1, :);  % sub-diagonal
c_diag = -Delta * y_flat(2:Nb, :);    % super-diagonal

% Solve all 12 systems at once
g_new_flat = thomas_solve_vec(a_diag, b_diag, c_diag, g_flat);
```

**Result: 24% speedup (8.4 min → 6.4 min)**

| Metric | Before | After |
|--------|--------|-------|
| KFE time/call | 250 ms | 150 ms |
| Total runtime | 8.4 min | 6.4 min |
| Speedup | - | 1.31× |

**Key lesson**: Vectorization beats sparse solvers when you have many small systems. The vectorized Thomas algorithm (with for-loops over rows, not systems) is efficient because MATLAB vectorizes across columns.

---

## 2026-01-01

### What we worked on

**Continued computational optimization - Anderson acceleration for the entry loop**

### What we discovered

The entry loop (fixed-point iteration for entry shares and bequest distribution) was a major bottleneck. Simple Picard iteration required ~60 iterations at some gamma values, each requiring a full KFE solve.

### Attempted optimization: Anderson acceleration for entry loop

**Goal**: Reduce the number of entry loop iterations (and thus KFE calls) by using Anderson acceleration instead of simple Picard iteration.

**Approach**:
Anderson acceleration is a technique for accelerating fixed-point iterations. It maintains a history of previous iterates and residuals (m=3), then computes an optimal linear combination:
```matlab
G_mat = g_hist(:, end-m_k+1:end) - g_curr;  % Residual differences
alpha_aa = G_mat \ g_curr;                   % Optimal weights
X_mat = x_hist(:, end-m_k+1:end) - x_curr;   % Iterate differences
x_aa = x_curr - X_mat * alpha_aa;            % Anderson extrapolation
x_new = x_aa + beta * (g_curr - G_mat * alpha_aa);  % Damped update
```

**Implementation**:
1. Added m=3 history buffer for iterates `x = [entry_shares; b_entry_loc]`
2. Added history buffer for residuals `g = [mature_shares; b_entry_implied] - x`
3. First 2 iterations use Picard (need history for Anderson)
4. Subsequent iterations use Anderson with fallback to Picard if unstable

**Result: SUCCESS! 39% speedup (6.4 min → 3.9 min)**

| Metric | Before (batched KFE) | After (+ Anderson) |
|--------|---------------------|-------------------|
| Total runtime | 6.4 min | **3.9 min** |
| Speedup | - | 1.64× |

**KFE call count comparison** (selected bisection iterations):

| Iteration | Before | After | Reduction |
|-----------|--------|-------|-----------|
| Bracket (gamma=-0.05) | 212 calls | 67 calls | 68% |
| Bracket (gamma=0.15) | 294 calls | 182 calls | 38% |
| Bisect 2 (hardest) | 535 calls | 216 calls | 60% |
| Final iterations | 9-14 calls | 4 calls | 60-70% |

**Key insight**: Anderson acceleration provides the biggest gains for difficult gamma values where the entry loop previously required many iterations. At convergence, it settles to ~4 iterations.

### Cumulative optimization results

| Optimization | Runtime | Speedup vs Previous | Cumulative Speedup |
|-------------|---------|--------------------|--------------------|
| Baseline (J=60, da=1, deterministic aging) | 27 min | - | - |
| Stochastic aging (20 → 5 child states) | 8.4 min | 3.2× | 3.2× |
| Batched KFE drift solves | 6.4 min | 1.3× | 4.2× |
| Anderson acceleration | 3.9 min | 1.6× | **6.9×** |

**Total speedup: 27 min → 3.9 min (6.9× faster)**

### Model results remain identical

All economic outputs unchanged:
- Mean parity: 1.45
- Ownership: 10.6%
- Childlessness: 7.6%
- Gradient: -0.001
- BGP consistency: 0.17% gap

### Failed optimization attempts (for reference)

1. **Thomas algorithm precomputation** (2025-12-31)
   - Tried to precompute matrix factorizations for reuse
   - Result: 6× SLOWER (250ms → 1500ms per KFE)
   - Reason: MATLAB for-loops can't compete with optimized sparse backslash

### Next steps

1. Profile to see if further optimization possible
2. Could try:
   - Better initial guess for entry shares from previous gamma
   - Secant method for outer bisection (risky if non-monotonic)
   - Parallelization of HJB (location-independent)
3. Focus on economics: parameters need recalibration for annual model

### Additional test: Anderson m=5

**Hypothesis**: More history (m=5 vs m=3) might accelerate convergence further.

**Result: SLIGHTLY WORSE (4.2 min vs 3.9 min)**
- The entry loop hit max_iter=60 once (Bisect 2, err=0.2139)
- Extra history doesn't help and may destabilize the iteration
- Reverted to m=3

### Timing breakdown (3.9 min run)

| Component | Time | % of Total |
|-----------|------|------------|
| HJB | 79.3 sec | 34% |
| KFE | 125.5 sec | 54% |
| Overhead | 28.8 sec | 12% |
| **Total** | **233.6 sec** | 100% |

**Bracket evaluation** (gamma=-0.05 and gamma=0.15) consumes ~81 sec (35% of total).

### Open questions

1. Is 3.9 min fast enough for calibration loops? (100 evaluations = 6.5 hours)
2. Should we implement warm-starting across gamma evaluations?
3. Can the bracket evaluation phase be shortened with a tighter initial bracket?

---

## 2026-01-02

### What we worked on

**Deep calibration analysis** - Understanding why the model fails to produce the target moments.

### Current model output vs targets

| Moment | Current (h_bar small) | Target | Gap |
|--------|----------------------|--------|-----|
| Mean parity | 1.45 | 1.7-1.9 | Too low |
| Ownership | 10.6% | 65% | **Way too low** |
| Childlessness | 7.6% | 15-20% | Too low |
| Fertility gradient | -0.001 | 0.3-0.5 | **Zero / wrong sign** |

### Calibration exploration results

Ran four parameter configurations:

| Config | h_bar_jump | nu_fert | Parity | Own% | Childless | Gradient |
|--------|------------|---------|--------|------|-----------|----------|
| Small h_bar | 0.08 | 3.0 | 1.45 | 10.6% | 7.6% | -0.001 |
| **Large h_bar** | 0.50 | 3.0 | 1.25 | **45.2%** | 24.4% | **+0.008** |
| Low nu_fert | 0.50 | 1.0 | 1.20 | **59.1%** | 34.3% | ? |
| Higher r_bar | 0.50 | 2.0 | 1.29 | 38.9% | **18.6%** | +0.009 |

### Key findings

1. **Large h_bar values are essential**:
   - h_bar_jump=0.50 (vs 0.08) improves ownership from 10.6% → 45%
   - Creates positive fertility gradient (0.008 vs -0.001)
   - The original large values in `setup_parameters()` were correct!

2. **The spatial gradient is still too small** (0.008 vs 0.3-0.5 target):
   - Subsistence cost differential for n=2: Superstar r×h_bar=0.161 vs Peripheral=0.025
   - Difference = 0.136, which is ~14% of income
   - With nu_fert=2-3, preference heterogeneity swamps cost differences

3. **Mean parity is too low** (1.2-1.3 vs 1.8 target):
   - Need higher theta_n (terminal warm-glow) to boost fertility

4. **Trade-off between ownership and childlessness**:
   - nu_fert=1.0 gives 59% ownership but 34% childlessness (too high)
   - nu_fert=2.0 gives 39% ownership but 18.6% childlessness (closer to target)

### Root cause of small gradient

The rent differential is not large enough:
- Current r_bar = [0.025, 0.050, 0.180] (annual)
- Superstar rent 7.2× peripheral rent, but both are small relative to income

For the housing mechanism to bite, need:
1. **Larger rent differentials** (maybe r_bar = [0.015, 0.05, 0.30])
2. **Higher h_bar values** (increase housing intensity of children)
3. **Lower nu_fert** (make fertility more responsive to cost differences)

### Proposed calibration strategy

**Phase 1: Get fertility level right**
- Increase theta_n (warm-glow) until mean parity ≈ 1.8
- Keep h_bar_jump=0.50, h_bar_n=0.20

**Phase 2: Get ownership right**
- Adjust nu_fert to balance ownership (target 65%) and childlessness (target 15%)
- nu_fert ≈ 1.5-2.0 looks promising

**Phase 3: Get gradient right**
- Increase r_bar differentials (especially superstar)
- May need h_bar_jump > 0.50 to get gradient 0.3-0.5
- Consider reducing wage gradient to avoid confounding

### Economic interpretation

The model is working correctly mechanically, but the parameterization doesn't generate strong enough spatial sorting by family type because:

1. **Housing costs are too uniform across locations** - Even superstar cities only charge 17% annual rent (of house value), which for a 0.93 housing requirement (n=2) is only 0.16 in cost.

2. **High preference heterogeneity (nu_fert=3)** - Type-I EV shocks with scale=3 mean fertility is almost random, not responsive to costs.

3. **Wages compensate for housing costs** - Superstar wages (1.2) partially offset the higher rents, reducing the net cost of children there.

### Next steps

1. Re-run with theta_n=2.5 or 3.0 to boost fertility to ~1.8
2. Test r_bar = [0.02, 0.05, 0.30] for larger rent differential
3. Grid search over (theta_n, nu_fert, r_bar) to find parameters matching all moments simultaneously

---

## 2026-01-03

### What we worked on

**Additional vectorization of remaining loops**

### Changes implemented

1. **Vectorized `precompute_drift_from_policy`** (lines 2395-2495):
   - Previously: 5 nested loops over (j, nn, cs, i, ten)
   - Now: Separated renters and owners, using 6D broadcasting
   - Removed ~14,000 loop iterations per call

2. **Vectorized HJB fertility logit step** (lines 1206-1250):
   - Previously: nested loops over (i, ten) with inner loop over nn_opt
   - Now: Process all (i, ten) simultaneously using 4D arrays
   - Eliminated I × n_tenure = 12 loop iterations per age

3. **Vectorized KFE fertility redistribution** (lines 1834-1871):
   - Previously: nested loops over (i, ten, nn)
   - Now: Compute births and redistribute mass using 4D arrays
   - Simplified from ~36 inner loops to vectorized operations

### Performance results

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Total runtime | 3.9 min | 3.7 min | **-5%** |
| HJB ms/call | 220 ms | 220 ms | Same |
| KFE ms/call | 175 ms | 170 ms | -3% |

**Economic outputs unchanged** (verification passed):
- Mean parity: 1.45
- Ownership: 10.6%
- Childlessness: 7.6%
- BGP gap: 0.17%

### Analysis

The speedup from these vectorizations is modest (~5%) because:

1. **Already heavily optimized**: The major bottlenecks (batched Thomas solves, Anderson acceleration) were already addressed in previous sessions.

2. **Loop overhead is small**: MATLAB's JIT compiler has gotten very good at simple loops. The remaining loops were already reasonably fast.

3. **Memory bandwidth limits**: Vectorized operations create larger intermediate arrays. For 6D arrays with shape (30, 4, 3, 60, 4, 5), memory throughput can limit gains.

### Cumulative optimization summary

| Optimization | Runtime | Cumulative Speedup |
|-------------|---------|-------------------|
| Baseline (deterministic aging) | 27 min | - |
| Stochastic aging | 8.4 min | 3.2× |
| Batched KFE drift | 6.4 min | 4.2× |
| Anderson acceleration | 3.9 min | 6.9× |
| **Additional vectorization** | **3.7 min** | **7.3×** |

**Total speedup: 27 min → 3.7 min (7.3× faster)**

### Remaining optimization opportunities

1. **Parallel HJB across ages**: Ages are independent in backward induction (but MATLAB parfor overhead may eat gains for only 60 iterations)

2. **GPU acceleration**: Could accelerate the large matrix multiplications, but would require significant code restructuring

3. **Reduce bisection iterations**: Currently 12 iterations; could try Brent's method but risk instability with non-monotonic F(gamma)

4. **Warmstart across gamma evaluations**: Currently implemented; could be more aggressive

At 3.7 min per solve, a 100-point calibration grid would take ~6 hours, which is reasonable for overnight runs.

---

## 2026-01-10

### What we worked on

**Replacing Logit/Gumbel shocks with Fréchet/power-mean shocks for BGP scale invariance**

The key insight: Logit discrete choice (with Gumbel shocks) is **translation-invariant** but not **scale-invariant**. On a BGP where values grow at rate $(1-\sigma)\gamma$, we need scale invariance to ensure that choice probabilities don't change over time. Fréchet (power-mean) discrete choice provides this.

### The problem with CRRA and Fréchet

Fréchet discrete choice requires positive arguments:
$$\pi_i \propto w_i \cdot V_i^\varepsilon, \quad IV = \left(\sum_i w_i V_i^\varepsilon\right)^{1/\varepsilon}$$

But CRRA utility with $\sigma > 1$ produces negative value functions:
$$u(c) = \frac{c^{1-\sigma}}{1-\sigma} < 0 \text{ when } \sigma > 1$$

### Initial approach: u_bar shift (abandoned)

First tried adding $\bar{u} = 8.0$ to flow utility to ensure $V > 0$. This works but has problems:
- Arbitrary choice of $\bar{u}$
- Changes Fréchet shares because Fréchet is NOT translation-invariant
- Not theoretically clean

### Final approach: Consumption-Equivalent (CE) transformation

The CE transformation converts negative CRRA values to positive levels:
$$CE = \left((1-\sigma) \cdot V\right)^{1/(1-\sigma)}$$

For $\sigma = 2$ and $V < 0$:
$$CE = \left(-V\right)^{-1} = \frac{1}{|V|}$$

**Key properties verified:**
1. CE > 0 when V < 0 and σ > 1 ✓
2. CE is monotone increasing in V ✓
3. CE ratios preserved under V → cV (what matters for Fréchet) ✓
4. Fréchet probabilities scale-invariant ✓
5. Inclusive value scales: IV(cV) = c·IV(V) ✓
6. Translation changes probabilities (confirms Fréchet, not logit) ✓

### Implementation: run_model_CE.m

Created new file `run_model_CE.m` with:

1. **u_bar = 0.0** (no artificial shift needed)

2. **Location choice with CE transformation** (lines 1152-1264):
```matlab
% CE TRANSFORMATION: CE = ((1-σ)*V)^(1/(1-σ))
CE_all = (one_minus_sigma * V_all).^(1/one_minus_sigma);
CE_all = max(CE_all, CE_floor);
log_Q_all = log(max(W_all, 1e-300)) + eps_loc * log(CE_all);
% ... softmax and convert back
V_inclusive = reshape(EV_CE.^one_minus_sigma / one_minus_sigma, Nb, n_parity, n_child);
```

3. **Fertility choice with CE transformation** (lines 1270-1318)

4. **Pure De Nardi bequest utility** (no additive warm-glow):
```matlab
bequest_util = theta0 * (theta1 + b)^(1-sigma) / (1-sigma)
```

### Problem discovered: Positive V values

Initial implementation produced positive V values because the terminal warm-glow `theta_n * n` was additive. With $n = 3$ children and reasonable theta_n, this pushed V positive.

**Issue:** For V > 0 with σ > 1, the CE transformation is undefined:
$$(1-\sigma) \cdot V < 0 \text{ when } V > 0$$

**Solution:** Use pure De Nardi bequest (CRRA in wealth), remove additive warm-glow `theta_n * n`. The bequest motive for children is already captured by De Nardi's `theta_0 * (theta_1 + b)^{1-σ}/(1-σ)`.

Alternatively, the model allows `psi_child` as flow utility from children, which integrates into V properly and keeps V negative.

### Results comparison

| Version | Pop Growth | Mean Parity | Ownership | Gradient |
|---------|------------|-------------|-----------|----------|
| u_bar = 8.0 | 1.0% | 1.24 | 52% | 0.074 |
| CE transform | 1.87% | 1.51 | 75% | 0.031 |

The CE version has higher fertility (because no dilution from u_bar) and different economics. Both converge successfully.

### Test files created

- `test_CE_transform.m` - 8 tests for CE transformation correctness (all pass)
- `test_V_sign.m` - Diagnostic for checking V stays negative
- `test_frechet_scale_invariance.m` - Scale invariance tests (all pass)

### Key theoretical insight

The relationship between Fréchet and CRRA:

For logit: $\pi_i \propto \exp(\nu \cdot V_i)$ → translation-invariant, NOT scale-invariant
For Fréchet: $\pi_i \propto w_i \cdot V_i^\varepsilon$ → scale-invariant, NOT translation-invariant

On a BGP where $V \to e^{(1-\sigma)\gamma t} \cdot v$, we need scale invariance so that probabilities remain constant. The CE transformation converts negative V to positive CE while preserving scale invariance.

### Files modified/created

- `code/matlab/run_model_CE.m` - New CE version of the model
- `code/matlab/test_CE_transform.m` - CE transformation tests
- `code/matlab/test_V_sign.m` - V sign diagnostic
- `code/matlab/test_frechet_scale_invariance.m` - Scale invariance tests

### Next steps

1. Decide whether to use CE version or u_bar version as baseline
2. Compare model moments under both approaches
3. Consider whether `psi_child` (flow utility from children) is needed for calibration
4. Run calibration with SMM framework

---

## 2026-01-14

### What we worked on

**Code consolidation and calibration setup**

### Changes made

1. **Main model file renamed**: The canonical model file is now **`run_model_main.m`**
   - Contains Cobb-Douglas Stone-Geary preferences
   - Fréchet discrete choice for location and fertility
   - CE (consumption-equivalent) transformation for negative CRRA values
   - Default `n_sub = 2` (two substeps for accuracy)
   - Logs output to `output/logs/` directory

2. **Folder structure reorganized**:
   ```
   code/matlab/
   ├── calibration/
   │   ├── results/       # SMM optimization results (.mat files)
   │   ├── logs/          # Calibration job logs (cluster)
   │   └── checkpoints/   # Mid-optimization checkpoints
   ├── output/
   │   ├── logs/          # Model run logs
   │   └── plots/         # Generated figures
   ```

3. **Calibration scripts updated to use `run_model_main`**:
   - `smm_objective.m` line 65: `run_model_CE` → `run_model_main`
   - `calibrate_smm.m` comment on line 34 updated
   - `test_smm_single.m` description updated

4. **SMM calibration framework**:
   - `calibrate_smm.m` - SLURM array job driver for cluster
   - `smm_objective.m` - Computes weighted squared moment deviations
   - `test_smm_single.m` - Local test script for single evaluation

### Parameters being calibrated (8 total)

| Parameter | Bounds | Description |
|-----------|--------|-------------|
| theta_n | [0.1, 3.0] | Terminal warm-glow for children |
| eps_fert | [2.0, 15.0] | Fréchet shape for fertility |
| h_bar_jump | [0.1, 1.0] | Housing requirement for first child |
| h_bar_n | [0.05, 1.0] | Housing per additional child |
| eps_loc | [0.5, 6.0] | Fréchet shape for location |
| psi | [0.02, 0.15] | Transaction cost |
| E_superstar | [1.0, 1.5] | Superstar city amenity |
| chi | [1.0, 1.10] | Owner service premium |

### Target moments (10 total)

| Moment | Target | Weight | Description |
|--------|--------|--------|-------------|
| mean_parity | 1.80 | 10.0 | Mean completed fertility |
| parity_gradient | 0.34 | 10.0 | Peripheral - Superstar (aggregate) |
| parity_gradient_renter | 0.30 | 20.0 | Renter fertility gradient |
| parity_gradient_owner | 0.35 | 20.0 | Owner fertility gradient |
| childless_rate | 0.15 | 5.0 | Fraction childless |
| own_rate | 0.65 | 10.0 | Aggregate ownership rate |
| own_gradient | 0.132 | 3.0 | Ownership: Peripheral - Superstar |
| own_family_gap | 0.14 | 2.0 | Family vs childless ownership gap |
| pop_peripheral | 0.40 | 1.0 | Population share in peripheral |
| pop_superstar | 0.25 | 2.0 | Population share in superstar |
| migration_rate | 0.023 | 5.0 | Annual interstate migration (Census 2023) |

### Key insight: Migration rate identifies ε_loc

The migration rate moment was **missing** from the original SMM objective, leaving ε_loc (Fréchet shape for location choice) unidentified. Added Census 2023 target of 2.3% annual interstate migration.

### Diagnosis from earlier calibration attempts

| Issue | Current | Target | Root Cause |
|-------|---------|--------|------------|
| TFR too low | 1.4-1.5 | 1.80 | θ_n too low |
| Gradient too small | +0.01 | +0.34 | ε_fert too high (8) swamps cost signal |
| Ownership too low | 45% | 65% | Need χ adjustment |

**Key fix needed**: With ε_fert = 8, preference heterogeneity dominates housing cost differences. The optimizer should explore lower ε_fert values (bounds allow [2, 15]).

### Childlessness target clarification

**IMPORTANT**: The childlessness target of 15% uses **"ever had children"** from vital statistics (Census Bureau 2024: 14.9%), NOT "NCHILD=0" from ACS which measures children currently in household (30-41%). The 15% target is correct.

### Next steps

1. Run `test_smm_single.m` locally to verify calibration infrastructure
2. Submit calibration jobs to cluster via `submit_calibration.sh`
3. Collect results with `collect_smm_results.m`

---

## 2026-01-17

### What we worked on

**Diagnosing and fixing broken model, understanding parameter scaling for annual time stepping**

### What we discovered

1. **Model was broken due to r_bar parameter mismatch**:
   - The model had two conflicting r_bar settings:
     - Auto-config block (line 131): `r_bar = [0.025; 0.050; 0.180]` (working)
     - setup_parameters default (line 379): `r_bar = [0.10; 0.20; 0.72] * da` (broken)
   - When calling without P_override, the bad values were used
   - This caused prices to hit p_max=10, ownership to collapse to 0%

2. **Rent/income ratios were unrealistically low**:
   - With the "working" r_bar values, rent/income was only 3-13%
   - DUE (Data on Urban Economics) targets ~24% for renters
   - This meant housing costs had no bite - households saved too much

3. **Wealth accumulation was excessive**:
   - Low rents → excessive savings → huge bequests (7+ years of income)
   - Target bequests should be ~1-2 years of income
   - This created a feedback loop: large bequests → wealthy entrants → easy ownership

### What we tested

**Higher r_bar values for realistic rent burdens**:
- Changed r_bar from [0.025; 0.050; 0.180] to [0.10; 0.20; 0.72]
- Increased p_max from 10 to 30 to accommodate higher prices

### Results

| Metric | Old r_bar | New r_bar | Target |
|--------|-----------|-----------|--------|
| Prices | [0.41, 0.73, 2.27] | [1.67, 2.90, 8.70] | - |
| Rent/Income | [3%, 5%, 13%] | [14%, 20%, 51%] | ~24% |
| Ownership | 96% | 85.8% | 65% |
| Mean parity | 0.95 | 0.98 | 1.8 |
| Bequests | 5-7 yrs income | 30-40 yrs income | 1-2 yrs |

### Problems remaining

1. **Bequests still way too high** (30-40 years of income vs 1-2 target)
   - Even with higher rents, wealth accumulates excessively
   - Options: lower interest rate, reduce bequest motive, add spending channels

2. **Parity too low** (~0.98 vs 1.8 target)
   - Everyone has ~1 kid, need more fertility variation

3. **Wealth profile "comically wrong"** (per user)
   - Need to diagnose why wealth accumulates so much

4. **Ownership still too high** (86% vs 65% target)
   - But moving in right direction

### Code improvements

1. **Added output saving to `run_model_jan17_debug.m`**:
   - Saved `sol`, `P`, `p_eq` to `output/last_run.mat` after each run
     (root output folder archived on 2026-05-07 under
     `calibration_archive/2026-05-07_root_output/`)
   - Can plot/explore without re-running the model

2. **Created `plot_last_run.m`**:
   - Quick visualization script that loads saved results
   - 8-panel figure: prices, rent/income, ownership by age, ownership by location, population, parity distribution, parity by location, wealth lifecycle

3. **Moved old plotting functions to archive**:
   - `plot_results_CE.m` → archive/
   - `plot_results.m` → archive/

4. **Set n_parity=3 as default** (0, 1, 2 children):
   - Faster (~2 min vs ~3 min with n_parity=4)
   - Minimal loss of accuracy

### Key insight

The model's economic mechanism depends critically on **rent burden**. When rent/income is too low (3-13%), households face no meaningful housing cost pressure, leading to:
- Excessive savings
- Unrealistic bequests
- Housing costs not affecting fertility decisions
- No spatial gradient

With realistic rent/income (14-51%), the mechanism starts to work, but wealth accumulation remains a problem.

### Next steps

1. Diagnose wealth accumulation:
   - Check interest rate vs discount rate (savings motive)
   - Check bequest motive parameters (theta0, theta1)
   - Consider if housing wealth is double-counted

2. Parameter exploration:
   - Need higher fertility (theta_n, eps_fert)
   - May need to reduce interest rate to lower savings

3. Understand why everyone has ~1 kid - need more childless AND more multi-child

### Files modified

- `run_model_jan17_debug.m`: Fixed r_bar defaults, added output saving, set n_parity=3
- `plot_last_run.m`: New plotting script
- Moved `plot_results_CE.m`, `plot_results.m` to archive/

### Current parameterization (in `run_model_jan17_debug.m`)

**Lifecycle:**
- J=60 periods, da=1.0 (annual), age_start=25, retire at 65
- Fertility window: ages 26-42, children mature at 18

**Preferences (Cobb-Douglas Stone-Geary):**
- sigma=2.0 (CRRA), rho=0.02 (discount rate)
- alpha_cons=0.70 (consumption share)
- c_bar(n) = 0.10 + 0.10*n (subsistence consumption)
- h_bar(n) = 0.50 + 2.0*(n≥1) + 0.50*n (subsistence housing)
- psi_child=0.05 (flow utility from children)
- theta_n=1.5 (terminal warm-glow for children)

**Fréchet shocks:**
- eps_fert=8.0 (fertility choice)
- eps_loc=2.0 (location choice)

**Housing:**
- q=0.04 (interest), delta=0.02 (depreciation), tau_H=0.01 (property tax)
- user_cost_rate = 0.07
- chi=1.05 (owner premium), psi=0.08 (transaction cost)
- phi=0.80 (LTV limit, 20% down)

**Locations:**
- Wages: w_bar = [0.85, 1.0, 1.20]
- Amenities: E_loc = [0.95, 1.0, 1.10]
- r_bar = [0.10, 0.20, 0.72] → equilibrium prices ~[1.67, 2.90, 8.70]

**Bequest:**
- phi_estate=0.80, b_entry_base=0.10
- De Nardi: theta0=1.5, theta1=0.01

---

## 2026-01-19

### What we worked on

**Calibration exploration for ownership, fertility, and fertility gradient**

Continuing from previous session where we had:
- Non-monotonic consumption bug fixed
- TFR at replacement (~2.06)
- Ownership very low (1.8%)
- No fertility gradient across locations

Goal: Achieve reasonable calibration with:
1. **Ownership rate ~50-65%** (US homeownership)
2. **TFR ~1.4-1.8** (below replacement)
3. **Fertility gradient ~0.1+** (Peripheral > Superstar)

### Key diagnostics performed

1. **Why is ownership so low?**
   - Down payments are affordable (smallest = 0.109 for Peripheral)
   - 97% of population CAN afford down payment
   - Owning is 59% cheaper than renting in flow costs
   - BUT: tenure choice compares V_rent(b) vs V_own(b - house_cost)
   - The **one-time wealth drop** from buying outweighs ongoing flow benefits
   - At b=20, V_rent=-0.685 but V_own_at_b_after_buy=-0.687 → renting still wins!
   - People need b≈21+ to choose ownership, but median wealth is only 1.26

2. **Ownership by age pattern**:
   - Age 25-40: ~0-8% ownership (young/poor rent)
   - Age 74: ~18-20% ownership (peak)
   - This matches lifecycle pattern but aggregate too low

3. **Wealth distribution**:
   - p10: b = -1.71, p50: b = 1.26, p90: b = 4.56
   - Mean wealth: 0.88
   - First wealth to choose own (Peripheral): b = 2.33
   - Mass with b >= 2.33: 47% (but only 11% actually own due to flow dynamics)

### Calibration attempts summary (11 versions)

| Version | chi | eps_fert | hR_max | b_entry | Own% | TFR | Gradient | Notes |
|---------|-----|----------|--------|---------|------|-----|----------|-------|
| v1 (baseline) | 1.10 | 4.0 | 5.5 | 0 | 4.1% | 1.59 | 0.104 | Fertility OK, ownership bad |
| v2 | 1.25 | 4.0 | 4.0 | 0 | 16.1% | 1.20 | 0.098 | Better ownership, TFR too low |
| v5 | 1.35 | 5.5 | 3.5 | 0 | 23.5% | 0.32 | 0.063 | Ownership good, fertility collapsed |
| v6 | 1.35 | 4.0 | 3.5 | 0 | 20.8% | 0.67 | 0.086 | Ownership target, TFR too low |
| v7 | 1.40 | 4.5 | 4.5 | 0 | 12.5% | **1.80** | **0.116** | Best fertility results |
| **v11** | 1.10 | 4.5 | 4.5 | **2.0** | **30.0%** | **1.70** | **0.121** | **Best overall balance** |

### Best calibration so far: Version 11

**Key insight**: Entry wealth is the missing lever! Setting `b_entry_fixed = 2.0` pushes more mass above the ownership threshold.

**Full parameter specification for v11**:

```matlab
%% ENTRY WEALTH - Key change
P_override.b_entry_fixed = 2.0;  % Start with wealth of 2.0 (was 0)

%% OWNERSHIP
P_override.chi = 1.10;              % Owner utility premium (constrained ≤1.15)
P_override.phi = 0.92 * ones(3,1);  % LTV = 92% (8% down payment)
P_override.h_own_min = 2.5;         % Smallest house size
P_override.h_own_max = 7.0;         % Largest house size
P_override.n_house = 3;             % 3 discrete house sizes
P_override.hR_max = 4.5;            % Rental housing cap

%% FERTILITY
P_override.eps_fert = 4.5;          % Fréchet shape for fertility
P_override.psi_child = 0.07;        % Flow utility from children
P_override.h_bar_jump = 2.2;        % Housing subsistence jump for first child
P_override.h_bar_n = 0.5;           % Additional housing per child
P_override.c_bar_n = 0.12;          % Additional consumption per child

%% HOUSING SUPPLY
P_override.H0 = [7.0; 4.0; 1.8];    % Housing supply by location
P_override.r_bar = [0.018; 0.05; 0.28];  % Base rents by location
```

**Results with v11**:
- Prices: [0.207, 0.534, 2.435]
- **Ownership: 30.0%** (P=64.3%, S=6.8%, X=0.5%)
- **TFR: 1.70** ✓
- **Gradient: 0.121** ✓
- Population growth: -0.50%
- Population shares: P=42.6%, S=37.4%, X=19.9%

### Key trade-offs identified

1. **chi (owner utility premium)**:
   - Higher chi → more ownership
   - User constraint: chi ≤ 1.10-1.15
   - Above 1.25 starts having equilibrium effects

2. **hR_max (rental cap)**:
   - Lower cap → forces families to own
   - BUT crushes fertility if cap < h_bar for families
   - Sweet spot: hR_max ≈ 4.5 (allows h_bar(n=1)=2.7 plus some buffer)

3. **eps_fert (Fréchet fertility shape)**:
   - Lower → more responsive to cost differences (larger gradient)
   - BUT can collapse TFR if costs dominate
   - Sweet spot: eps_fert ≈ 4-5

4. **b_entry_fixed (entry wealth)**:
   - Higher → more people above ownership threshold
   - Interpretable as parental transfers or accumulated savings
   - With b_entry=2.0, ownership jumped from 12.5% to 30%

### Root cause of low ownership

The tenure choice compares:
- V_rent(b) at current wealth
- V_own(b - house_cost) at post-purchase wealth

Even when owning has better flow utility (chi bonus + cheaper per-period costs), the **one-time wealth drop** from buying hurts. The value function is quite flat at high wealth, so the flow benefits don't accumulate fast enough to offset the purchase cost.

This is economically sensible - it matches the empirical pattern where young households rent (can't afford wealth drop) and older/wealthier households own.

### Files created

- `run_calibration_fix[1-11].m` - Calibration exploration scripts
- `diagnose_ownership2.m` - Down payment vs wealth analysis
- `diagnose_tenure_choice.m` - Tenure choice policy analysis
- `diagnose_tenure_value.m` - Value function comparison for buy vs rent
- `diagnose_ownership_by_age.m` - Ownership lifecycle profile
- `diagnose_wealth_dist.m` - Wealth distribution analysis
- `diagnose_tenure_flows.m` - KFE tenure transition diagnostics

### Next steps for calibration launch

1. **Use v11 as starting point**: b_entry=2.0 is the key innovation
2. **Consider b_entry as calibration parameter**: Can match bequest/wealth targets
3. **May need b_entry > 2.0 for 50%+ ownership**: Each unit of entry wealth adds ~5-10% ownership
4. **Keep chi ≤ 1.15 per user constraint**

### Open questions

1. What is realistic b_entry? Depends on parental wealth transfers
2. Should b_entry vary by location (Superstar parents wealthier)?
3. Can ownership reach 65% without unrealistic chi?
4. Trade-off: higher ownership may require sacrificing fertility targets

---

## 2026-03-11

### What we worked on

**Built and hand-tuned the center-periphery (2-location) model**

Major structural change: collapsed the 3-location model (Peripheral/Secondary/Superstar) with endogenous wages into a 2-location center-periphery model with common wages. This follows the March draft's simplification.

### Key design decisions

1. **2 locations**: Center (expensive, amenity-rich) and Periphery (cheap, elastic supply)
2. **Common wages**: Removes agglomeration/wage gradient. Sorting driven entirely by amenities vs housing costs.
3. **No wage iteration**: Equilibrium iterates on (prices, entry_shares) only — converges in ~22 iterations, ~17 seconds.
4. **Reference paper**: Moreno-Maldonado & Santamaria (2024), "Delayed Childbearing and Urban Revival"

### Files created

| File | Purpose |
|------|---------|
| `code/matlab/tests_february/run_model_cp.m` | 2-location solver (based on v10, ~3150 lines) |
| `code/matlab/tests_february/run_cp_tuned.m` | Hand-tuned parameter overrides |
| `code/matlab/tests_february/plot_cp.m` | 6-figure diagnostic plotting suite |

### Hand-tuning progression (4 rounds)

| Round | TFR  | Own%  | Childless P/C | Migration | W/I  | Key changes |
|-------|------|-------|---------------|-----------|------|-------------|
| Default | 2.48 | 84.1% | 3.5/6.3% | 4.0% | 6.05 | Baseline defaults |
| R1 | 1.18 | 65.7% | 30.1/66.2% | 3.8% | 6.26 | Overcorrected fertility |
| R2 | 1.88 | 71.1% | 10.9/31.6% | 3.2% | 4.61 | Back off, inverted own gradient |
| R3 | 1.93 | 68.4% | 11.5/27.7% | 4.7% | 4.80 | Fixed own gradient, phi=0.75 |
| **R4** | **2.01** | **68.1%** | **9.7/24.1%** | **4.8%** | **4.83** | **Current best** |

### What's working

- **Ownership**: 68.1% (target ~65%) — good.
- **Migration**: 4.8% (target ~4.5%) — good.
- **Rent ratio C/P**: 2.12 (target ~2.0) — good.
- **Fertility gradient direction**: Correct (Periphery > Center).
- **Ownership gradient direction**: Correct (Periphery > Center).
- **Policy functions**: Clean, monotone, economically sensible.
- **Convergence**: 22-25 iterations, 15-25 seconds. No numerical issues.

### What needs work

- **TFR = 2.01 vs ~1.70**: Too high overall.
- **W/I = 4.83 vs ~1.62**: Structural issue — includes housing equity.
- **Center childless = 24% vs ~18%**: Gradient too steep.
- **Fertility probability flat in wealth**: kappa_fert too high, Gumbel noise swamps value differences.

### Key insight: calibration targets need empirical grounding

Before running PSO, need to:
1. Define what "Center" vs "Periphery" means in the data (MSA size cut? rent level?)
2. Compute fertility, ownership, migration by this partition from ACS/Census
3. Pin down which moments identify which parameters in the 2-location setup
4. Address W/I definition (liquid wealth only? total including housing equity?)

### Next steps

1. **Empirical work**: Define Center/Periphery in data, compute calibration targets
2. **Review Moreno-Maldonado & Santamaria**: Their calibration strategy and target choices
3. **Set up PSO calibration** for the CP model once targets are firm
4. **Consider W/I redefinition**: Liquid wealth / income may be more matchable

---

## 2026-07-23

### What we discovered

The July 22 one-shot calibration's three saving/bequest moments used a hybrid
balance sheet: inherited beginning-of-period liquid wealth was paired with the
newly chosen tenure before the housing transaction was applied. The resulting
wealth, bequest-flow, and old-estate-dispersion fits are invalid. At the
certified winner, the old-estate p90/p50 is `3.552` under the hybrid
measurement but only `1.911` at coherent beginning timing and `2.147` after
the transaction.

The final-winner Jacobian is numerically full rank but has condition number
`8595.69`; it describes the invalid old objective and is not calibration
evidence. Once the three invalid rows are removed, the exercise has 11 valid
moments for 14 free parameters.

### Actions

- Disabled new searches under the invalid contract.
- Stopped only the related annual-beta profile jobs after about 20 minutes and
  preserved their checkpoints.
- Built exact timing decompositions and a twice-repeated seven-cell strict beta
  snapshot. Even annual beta `0.9999` leaves coherent wealth/earnings at
  `6.042` or below versus target `6.9`.
- Recorded exact artifacts and the repair sequence in
  `CALIBRATION_STATUS.md`.

### Next steps

Use one coherent stock timing for cross-sectional wealth and estate moments;
then implement the bequest flow at the actual death point in the model's
within-period sequence. Only after that should the target table, Jacobian, and
calibration be rebuilt.

---

## Template for future sessions

### Date: YYYY-MM-DD

### What we worked on
-

### What we discovered
-

### What we tested
-

### Results
-

### Problems encountered
-

### Next steps
-

### Open questions
-
