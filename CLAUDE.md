# CLAUDE.md - Project Instructions

## Project Overview

This is an academic research project developing a **structural spatial lifecycle model with endogenous fertility and housing tenure choice**. The goal is to explain why fertility differs systematically across locations with different housing costs, how housing markets mediate fertility decisions, and the resulting spatial sorting of households by family type.

**Discipline**: Quantitative Macroeconomics / Spatial Economics / Urban Economics  
**Implementation**: MATLAB  
**Target outlet**: Top economics journal (AER, QJE, Econometrica, REStud)

---

## Research Question

> How do housing costs affect fertility decisions, and what are the aggregate implications for spatial sorting and demographic patterns across cities?

The model integrates three literatures:
1. **Spatial equilibrium** (Rosen-Roback tradition)
2. **Lifecycle fertility** (Becker quality-quantity, Doepke-Kindermann)
3. **Housing tenure choice** (Henderson-Ioannides)

---

## Standards for Assistance

### Academic Rigor

All work must meet publication standards for a top economics journal:

- **Theoretical precision**: Every equation must be correctly specified. Define all variables, parameters, and functional forms explicitly. Use standard notation from the literature.

- **Mathematical accuracy**: Derivations must be error-free. When implementing equations in code, verify they match the theoretical specification exactly.

- **Empirical grounding**: Calibration targets and parameter values must be justified with citations to data sources or published estimates. Flag uncertain or unverified targets.

- **Reproducibility**: Code must be clean, documented, and deterministic.

### When Assisting with This Project

1. **Cite relevant literature** when discussing economic mechanisms or calibration choices (e.g., "Following Doepke & Kindermann (2019)..." or "Standard in the housing literature, see Sommer et al. (2013)...")

2. **Use LaTeX notation** for mathematical expressions in discussions and documentation

3. **Question assumptions** - if a calibration target, parameter bound, or modeling choice seems unjustified or inconsistent with the literature, flag it explicitly rather than proceeding

4. **Verify numerical results** - if model output violates basic economic logic (non-monotonic value functions, negative probabilities, implausible moments), investigate the cause rather than accept the output

5. **Maintain terminological precision**:
   - Distinguish "fertility" (flow/hazard) from "completed fertility" (stock)
   - Distinguish "childlessness" (extensive margin) from "number of children" (intensive margin)
   - Distinguish "housing services" from "housing wealth" from "house prices"

6. **Think about identification** - when calibrating or estimating, consider what variation identifies each parameter and whether targets are informative

---

## Model Structure

### Overview

Continuous-time OLG model solved in discrete age periods. Agents are heterogeneous in:
- Wealth (continuous state, endogenous)
- Number of children (discrete state, endogenous, irreversible)
- Location (discrete choice)
- Housing tenure (discrete choice: rent vs own)

### Geography

Three location types representing the spectrum from rural/peripheral to superstar metros:
- **Peripheral**: Low wages, low housing costs, family-friendly amenities
- **Secondary**: Baseline/normalized
- **Superstar**: High wages, high housing costs, career amenities, scarce housing supply

### Key Economic Mechanisms

1. **Housing-fertility interaction**: Children require space. Subsistence housing requirements increase with family size. In expensive locations, this creates a higher effective cost of children.

2. **Spatial sorting**: Households select into locations based on wages, housing costs, and amenities. Family status affects optimal location choice.

3. **Tenure choice**: Rent vs own tradeoff depends on wealth, borrowing constraints, and expected tenure length.

4. **Lifecycle dynamics**: Fertility is concentrated at young ages; wealth accumulates over the lifecycle; location and tenure may change.

### Preferences

Stone-Geary utility with subsistence requirements for both consumption and housing that depend on family size. Key feature: nonlinear (discrete jump) increase in housing requirement for first child.

### Equilibrium

Stationary spatial equilibrium where:
- Housing markets clear in each location
- Population distribution is stationary given fertility and migration flows
- Prices adjust endogenously

---

## Calibration Philosophy

### Targets
Calibration targets should be:
- Empirically documented with clear data sources
- Relevant to the mechanisms the model is designed to capture
- Informative for identifying the parameters they're matched to

### Parameters
Distinguish between:
- **Externally calibrated**: Set from external data/literature (e.g., interest rates, depreciation)
- **Internally calibrated**: Chosen to match model moments to data targets
- **Normalized**: Set to pin down scale or units

### Validation
After calibration, examine:
- Untargeted moments (out-of-sample fit)
- Comparative statics (do responses go in expected directions?)
- Sensitivity to key parameters

---

## Diagnosing Numerical Issues

If results look wrong, check in order:
1. **Value function monotonicity**: V must be strictly increasing in wealth
2. **Market clearing**: Housing demand must equal supply in each location
3. **Probability bounds**: All choice probabilities in [0,1]
4. **Boundary behavior**: Check solutions at grid edges (poorest, richest, oldest)
5. **Convergence**: Did price iteration and value function iteration converge?

---

## References

Key papers in the literature that inform this project:

**Spatial equilibrium**: Rosen (1979), Roback (1982), Moretti (2011)

**Fertility economics**: Becker (1960), Becker & Lewis (1973), Doepke & Kindermann (2019)

**Housing and lifecycle**: Henderson & Ioannides (1983), Sommer, Sullivan & Verbrugge (2013)

**Superstar cities**: Gyourko, Mayer & Sinai (2013), Hsieh & Moretti (2019)

---

*Last updated: December 2025*
