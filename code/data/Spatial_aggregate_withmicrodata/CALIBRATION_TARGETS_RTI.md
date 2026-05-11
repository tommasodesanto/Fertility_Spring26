# Calibration Targets: RTI-Based Location Classification

**Date**: January 2026
**Data Source**: ACS microdata 2005-2023 (IPUMS), 10% sample
**Classification**: Rent-to-Income (RTI) terciles

---

## Location Classification Method

Rather than using arbitrary population thresholds (e.g., MSA >= 2M for "superstar"), we classify locations by **housing affordability** measured as rent-to-income ratio (RTI = annual rent / household income).

**Rationale**: The model's mechanism is about housing costs affecting fertility. Classifying by RTI directly captures this mechanism.

### RTI Tercile Cutoffs (Population-Weighted)

| Location Type | RTI Range | Interpretation |
|---------------|-----------|----------------|
| **Peripheral** | RTI < 0.224 | Low housing cost burden |
| **Secondary** | 0.224 <= RTI < 0.240 | Medium housing cost burden |
| **Superstar** | RTI >= 0.240 | High housing cost burden |

---

## Summary of Calibration Targets

### Population Shares

| Location | Share | Note |
|----------|-------|------|
| Peripheral | 51% | Low-cost areas (incl. non-metro) |
| Secondary | 27% | Medium-cost metros |
| Superstar | 22% | High-cost metros |

### Total Fertility Rate (TFR)

Computed using FERTYR (births in last 12 months), summing ASFR across 5-year age groups.

| Location | TFR |
|----------|-----|
| Peripheral | 1.95 |
| Secondary | 1.82 |
| Superstar | 1.69 |
| **Gradient (P - S)** | **+0.26** |

### TFR by Tenure

| Tenure | Peripheral | Secondary | Superstar | Gradient |
|--------|-----------|-----------|-----------|----------|
| Owners | 2.00 | 1.90 | 1.71 | **+0.29** |
| Renters | 1.87 | 1.72 | 1.67 | **+0.20** |

**Note**: Owners show a *larger* gradient than renters. This may reflect selection (renters in superstars already negatively selected on fertility).

### Homeownership Rate

| Location | Ownership Rate |
|----------|----------------|
| Peripheral | 68.7% |
| Secondary | 65.5% |
| Superstar | 58.3% |
| **Gradient (P - S)** | **+10.4 pp** |

### Family-Childless Ownership Gap

| Location | Gap (Family - Childless) |
|----------|--------------------------|
| Peripheral | +12.1 pp |
| Secondary | +14.2 pp |
| Superstar | +10.6 pp |
| **Average** | **+12 pp** |

### Childlessness

| Source | Value | Note |
|--------|-------|------|
| Census Bureau vital stats | **15%** | Women 45-50, "ever had children" |
| ACS NCHILD = 0 | 30-35% | Confounded by children leaving home |

**Target the 15% from vital statistics**, not ACS NCHILD.

---

## The TFR vs. NCHILD Paradox

| Measure | Peripheral | Superstar | Gradient | Sign |
|---------|-----------|-----------|----------|------|
| **TFR** (flow) | 1.95 | 1.69 | +0.26 | Positive |
| **NCHILD** (stock) | 1.23 | 1.38 | -0.15 | **Negative** |

**Interpretation**:
- TFR shows actual fertility is *higher* in cheap areas (housing cost mechanism)
- NCHILD shows *more* children in household in expensive areas (selection + children leaving home earlier in cheap areas)

The model should generate **both** patterns.

---

## Wage and Rent Levels

### Median Individual Earnings (Workers 25-60)

| Location | Median Earnings | Relative (vs Secondary) |
|----------|-----------------|-------------------------|
| Peripheral | $32,000 | 0.84 |
| Secondary | $38,020 | 1.00 |
| Superstar | $41,400 | **1.09** |

### Median Monthly Rent

| Location | Median Rent | Relative (vs Secondary) |
|----------|-------------|-------------------------|
| Peripheral | $550 | 0.71 |
| Secondary | $770 | 1.00 |
| Superstar | $1,100 | **1.43** |

**Key insight**: Superstars have only 9% higher wages but 43% higher rent. This is the housing cost burden channel.

---

## Model Parameters (Suggested)

Based on relative values (Secondary = 1.00):

| Parameter | Peripheral | Secondary | Superstar |
|-----------|-----------|-----------|-----------|
| Wage (w_bar) | 0.84 | 1.00 | 1.09 |
| Rent (r_bar) | 0.71 | 1.00 | 1.43 |

---

## Files Generated

All saved to `calibration_targets_output/`:

- `tfr_by_location_rti.csv`
- `tfr_by_tenure_location_rti.csv`
- `ownership_by_location_rti.csv`
- `population_shares_rti.csv`
- `income_by_location_rti.csv`
- `nchild_by_location_rti.csv`
- `family_ownership_gap_rti.csv`

---

## Script

To regenerate these targets, run:
```r
source("compute_rti_targets.R")
```
