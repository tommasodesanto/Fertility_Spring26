# Calibration Strategy Review - Spatial Fertility Model

## Context
I'm working on a structural spatial lifecycle model with endogenous fertility and housing tenure choice. The goal is to explain why fertility differs systematically across locations with different housing costs. I need help reviewing and improving my calibration strategy, with the primary goal of achieving **below-replacement fertility** (negative population growth).

---

## Model Overview

### Agents
- Continuous-time OLG model solved in discrete age periods (ages 25-84, J=60 periods)
- Households are couples (2 adults) making joint decisions
- Heterogeneous in: wealth (continuous), number of children (discrete: 0, 1, 2+), location (discrete), housing tenure (rent vs own)

### Geography (3 locations)
1. **Peripheral**: Low wages, low housing costs, family-friendly amenities
2. **Secondary**: Baseline/normalized (wages=1, amenity=1, housing supply=1)
3. **Superstar**: High wages, high housing costs, career amenities, scarce housing

### Key Decisions
1. **Fertility**: Irreversible choice to have children (only during fertile ages 25-44)
2. **Location**: Where to live (with Fréchet taste shocks)
3. **Tenure**: Rent vs own (with transaction costs, down payment constraint)
4. **Consumption/Saving**: Standard lifecycle with borrowing constraint

### Preferences
Stone-Geary Cobb-Douglas utility:
```
u(c, h, n) = α log(c - c̄(n)) + (1-α) log(h - h̄(n)) + ψ_child * n
```
where:
- `α = 0.7` (externally calibrated → 30% housing expenditure share)
- `c̄(n)` = subsistence consumption (increases with children)
- `h̄(n) = h̄_base + h̄_jump * 1{n≥1} + h̄_n * n` = subsistence housing (discrete jump for first child)
- `ψ_child` = flow utility from children

Terminal value includes warm-glow for children: `θ_n * n`

### Fertility Choice
Fréchet discrete choice over {have child, don't have child}:
- Lower `eps_fert` → more children (less randomness in choice)
- Higher `theta_n` → more children (stronger bequest/warm-glow motive)
- Higher `psi_child` → more children (flow utility)
- Higher `h_bar_jump` → fewer children (housing cost of kids)

### Housing Market
- Owners: Pay price `p`, get service flow `χ * h` where `χ > 1` (ownership premium)
- Renters: Pay rent `r = user_cost_rate * p`, capped at `h_R_max`
- Down payment constraint: `b ≥ -φ * p * h` where `φ = 0.80` (20% down payment)
- Housing supply: `H_s(p) = H̄ * p^η` where `η` varies by location (Saiz elasticities)

### General Equilibrium
- Housing market clears in each location
- Wages determined by productivity + agglomeration: `w = A * L^ξ` where `ξ = 0.04`
- Stationary population distribution (births = deaths + migration flows balance)

---

## Current Calibration Results

### Best Result (Loss = 43.58, from overnight cluster run)

**Parameters found:**
| Parameter | Value | Bounds | Description |
|-----------|-------|--------|-------------|
| theta_n | 3.01 | [2, 8] | Terminal warm-glow |
| eps_fert | 3.08 | [1.5, 6] | Fertility Fréchet shape |
| h_bar_jump | 1.69 | [0.5, 3] | Housing jump for 1st child |
| h_bar_n | 0.63 | [0.2, 1] | Housing per additional child |
| psi_child | 0.13 | [0.05, 0.4] | Flow utility from children |
| eps_loc | 1.91 | [1, 4] | Location Fréchet shape |
| psi | 0.09 | [0.02, 0.15] | Transaction cost |
| chi | 1.11 | [1.05, 1.3] | Ownership premium |
| E_peripheral | 1.10 | [0.8, 1.2] | Peripheral amenity |
| E_superstar | 1.18 | [1.0, 1.5] | Superstar amenity |
| A_prod_P | 0.77 | [0.7, 1] | Peripheral productivity |
| A_prod_X | 1.27 | [1, 1.3] | Superstar productivity |
| H_bar_P | 1.89 | [1, 2.5] | Peripheral housing supply |
| rho | 0.031 | [0.01, 0.05] | Discount rate |
| b_entry | 1.28 | [0.3, 1.5] | Entry wealth (hit upper bound!) |

**Moment comparison:**
| Moment | Target | Model | Error |
|--------|--------|-------|-------|
| TFR (mean_parity × 2) | 1.62 | 2.26 | +39% ❌ |
| Fertility gradient (P-S) | 0.26 | 0.10 | -62% ❌ |
| Childless rate | 15% | 17% | OK |
| Ownership rate | 65% | 76% | +16% |
| Ownership gradient (P-S) | 10.4pp | 34.3pp | +230% ❌ |
| Family-childless own gap | +12pp | -4.6pp | Wrong sign! ❌ |
| Pop share Peripheral | 51% | 45% | OK |
| Pop share Superstar | 22% | 21% | OK |
| Wage ratio P/S | 0.84 | 0.90 | OK |
| Wage ratio X/S | 1.09 | 1.13 | OK |
| Rent ratio X/S | 1.43 | 2.87 | +101% ❌ |

---

## Key Problems

### 1. Fertility Too High (TFR = 2.26 vs target 1.62)
The model produces **positive population growth** (+0.25%) when the US has **below-replacement fertility**. The calibration cannot push fertility low enough even with:
- Low `psi_child` (0.13)
- Moderate `theta_n` (3.0)
- Moderate `eps_fert` (3.08)

**Question**: What's creating the strong intrinsic fertility incentive? Is it the terminal value specification? The Fréchet choice structure?

### 2. Fertility Gradient Inverted
Data shows Peripheral > Secondary > Superstar fertility (people have more kids in cheaper places).
Model shows the **opposite**: Superstar has highest fertility!

**Hypothesis**: High wages in Superstar dominate housing cost effects, making children more affordable there in terms of foregone consumption.

### 3. Ownership-Fertility Relationship Wrong
Data: Families own more than childless (positive gap).
Model: Childless own more than families (negative gap).

**Hypothesis**: In the model, having children consumes resources that would otherwise go to down payment. The housing cost of children crowds out ownership.

### 4. Superstar Rents Too High
Model rent ratio 2.87 vs data 1.43. Housing supply in Superstar (`H_bar_X = 0.6`, fixed) may be too constrained.

### 5. Entry Wealth Hit Upper Bound
The optimizer pushed `b_entry` to 1.28 (near max 1.5), suggesting young agents are severely wealth-constrained. This might be compensating for other model deficiencies.

---

## Externally Calibrated Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| α (consumption share) | 0.70 | → 30% housing expenditure |
| φ (LTV ratio) | 0.80 | Standard US mortgage (20% down) |
| q (interest rate) | 0.04 | External |
| ξ (agglomeration) | 0.04 | Combes et al. |
| η_supply | [2.0, 1.5, 0.8] | Saiz (2010) - P, S, X |
| H_bar_X | 0.6 | Fixed (scarce in Superstar) |

---

## Target Moments and Weights

**Fertility targets:**
- TFR = 1.62 (US 2023) → mean_parity = 0.81 (weight 10)
- Fertility gradient P-S = 0.13 (weight 10)
- Fertility gradient by tenure: owners 0.12, renters 0.08 (weight 20 each)
- Childless rate = 15% (weight 5)

**Ownership targets:**
- Aggregate ownership = 65% (weight 10)
- Ownership gradient P-S = 10.4pp (weight 3)
- Family-childless gap = +12pp (weight 2)

**Population:**
- Pop share Peripheral = 51% (weight 1)
- Pop share Superstar = 22% (weight 2)

**GE equilibrium:**
- Wage ratio P/S = 0.84 (weight 5)
- Wage ratio X/S = 1.09 (weight 5)
- Rent ratio X/S = 1.43 (weight 5)

**Intertemporal:**
- Wealth-to-income ratio = 1.62 (weight 5) - for calibrating ρ

---

## Questions for Review

1. **How to achieve below-replacement fertility?**
   - Should we add explicit costs of children (childcare, time costs)?
   - Is the terminal warm-glow specification too strong?
   - Should fertility have a direct utility cost?

2. **How to get the fertility gradient right?**
   - Housing costs should make children more expensive in Superstar
   - But wages are also higher there - is the income effect dominating?
   - Do we need location-specific fertility costs?

3. **Is the model identified?**
   - With 15 parameters and ~14 moments, are we overparameterized?
   - Which parameters are well-identified vs poorly identified?
   - Should we externally calibrate more parameters?

4. **What's the key mechanism?**
   - The paper's thesis is that housing costs affect fertility through the space requirement channel
   - Is `h_bar_jump` (housing requirement for first child) doing this work?
   - What moments actually identify `h_bar_jump` vs `theta_n` vs `psi_child`?

5. **Should we simplify?**
   - Two parities (0, 1+) instead of three (0, 1, 2+)?
   - Exogenous wages instead of GE?
   - Fixed locations (no migration)?

---

## Code Files

- Main model: `run_model_jan19.m`
- Calibration driver: `calibration/calibrate_smm.m`
- Objective function: `calibration/smm_objective.m`
- Results collection: `calibration/collect_results.m`

---

## What I Need Help With

1. Review the calibration strategy - are the targets and weights sensible?
2. Diagnose why fertility is too high and has wrong gradient
3. Suggest model modifications to achieve below-replacement fertility with correct spatial gradient
4. Identify which parameters to externally calibrate vs estimate
5. Think about identification - what variation identifies each parameter?
