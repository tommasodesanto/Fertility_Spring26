# Code Audit Prompt

Use this prompt in a fresh Claude chat. Attach the file `code/matlab/tests_february/run_model_v8_decoupled.m` (2645 lines, single self-contained MATLAB script).

---

## Prompt

I need a comprehensive audit of this MATLAB code. It implements a **structural spatial lifecycle model** with endogenous fertility and housing tenure choice, solved in continuous time (HJB/KFE) with discrete age periods. The model is used for an academic economics paper targeting a top journal.

**The code is a single file containing all functions.** Please read it in full and provide a detailed audit covering the sections below. Be specific — cite line numbers, variable names, and equations. Flag anything that looks wrong, suspicious, or inconsistent.

### 1. Model Structure & State Space
- What are all the state variables? (continuous, discrete)
- What are the dimensions of the value function V and distribution g?
- Are the grid dimensions consistent across HJB, KFE, and aggregation?
- Is the state space correctly indexed everywhere? (Watch for off-by-one errors in tenure, parity, child-state indices.)

### 2. HJB Solver (`solve_hjb_cobb_douglas`)
- **Terminal value**: Is the boundary condition at the last age correct? Does it handle all tenure × parity × child-state combinations?
- **Consumption policy**: Is the Euler equation / FOC correctly implemented? Check the Stone-Geary Cobb-Douglas utility specification — is the first-order condition right?
- **Housing rental choice**: Is the optimal renter housing `hR_pol` derived correctly from the FOC?
- **Tenure choice**: Is the argmax over {rent, own_H1, own_H2, ...} correctly comparing flow utilities + continuation values?
- **Location choice**: Is the Gumbel/Logit formula correct? Is the inclusive value computed correctly?
- **Fertility choice**: Same — is the Logit choice over parity levels correct? Are the continuation values for different parities right?
- **Finite differences**: Are forward/backward differences computed correctly? Is the upwind scheme right (forward drift → forward difference, backward drift → backward difference)?
- **Tridiagonal solve**: Is the implicit scheme set up correctly? Are the A matrix coefficients right?
- **Drift computation**: Is the savings drift `s = y - c - r*h` computed correctly for each tenure type?

### 3. KFE Solver (`solve_kfe`)
- **Transition logic**: At each age boundary, how are discrete choices (location, tenure, fertility) applied? Is mass conserved?
- **Entry distribution**: How are newborns injected? Is the entry wealth and location distribution correct?
- **Movers**: When households change location, do they correctly become renters? Is equity correctly liquidated?
- **Stochastic aging**: If children "mature" stochastically, is the transition matrix between child states correct?
- **Mass conservation**: Does total mass sum to the expected amount at each age? Are there any leaks?

### 4. Equilibrium Loop (`solve_equilibrium`)
- Is the damped iteration on (prices, wages, entry_shares) correct?
- Is the housing market clearing condition right? (Supply curve, demand aggregation)
- Is the wage agglomeration formula correct?
- Is the convergence criterion appropriate?

### 5. Aggregation & Moments
- **Ownership rate**: Correctly computed from the distribution g?
- **TFR / mean parity**: Correctly computed?
- **Population shares**: Correctly computed?
- **Housing demand**: Aggregated correctly from renter policies + owner fixed sizes?
- **Wealth-to-income ratio**: Does it include housing equity or just liquid wealth?
- **Migration rate**: Correctly computed from location choice probabilities?

### 6. Economics Consistency
- Does the model satisfy Walras' law? (Are all budget constraints consistent?)
- Is the borrowing constraint correctly enforced? (Renters b≥0, owners b≥-(1-ψ)pH)
- Is the down payment constraint correctly enforced?
- Are bequests handled correctly at the terminal age?
- Is the pension system correctly specified? (Fixed pension, payroll tax on workers)

### 7. Potential Bugs & Red Flags
- Any variables that are defined but never used?
- Any hardcoded values that should be parameters?
- Any places where array dimensions might not match?
- Any numerical issues (division by zero guards, NaN propagation)?
- Any places where the code does something different from what the comments say?

### 8. Known Issues to Investigate
These are issues we've already observed in the model output. For each, trace the code to identify the root cause:
- **Ownership collapses from ~66% at age 54 to ~1% at age 64.** Almost everyone sells their house at retirement. Is this because movers become renters? Or is there a bug in the tenure choice at retirement?
- **Consumption has a discontinuous jump near wealth = 0.** Is this a grid artifact from the multi-segment wealth grid? Or a genuine kink from the borrowing constraint?
- **Savings drift is non-monotonic and discontinuous.** Related to the consumption issue?

Please be thorough and specific. I'd rather have false positives (flagging something that turns out to be fine) than miss a real bug.
