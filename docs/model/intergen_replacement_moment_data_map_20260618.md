# One-Market Intergen Replacement-Moment Data Map

Date: 2026-06-18

## Status

This is a candidate data map, not a final target system. All candidate moments
below must be reaudited before they enter the SMM objective: source file,
sample restriction, household unit, weights, formula, room/bedroom convention,
wealth definition, and timing must be checked against the model object. The
current concern is that some hard targets may be empirically real but not the
right counterparts for the one-market intergenerational model.

No target should be dropped or demoted without preserving identification. If a
moment is removed from the hard objective, replace it with a moment for the same
parameter block, fix the affected parameter externally, or state explicitly that
the calibration is underidentified.

## Recommended Data Conventions

Use ACS/MMS household-head extracts as the primary source for housing, tenure,
rooms, and renter costs. They are closest to the current ownership target
construction. Use AHS as a robustness and stock/menu source, especially for
national tenure-size distributions.

Use PSID for lifecycle wealth, completed-fertility links, birth-related housing
increments, and down-payment/access moments. For old-age bequest moments, prefer
nonhousing or liquid wealth as the model-consistent estate object; total net
worth should be a robustness check because the theory treats retained housing as
old consumption that later returns to the stock/passive pool, not as additional
entrant wealth.

Use rooms as the primary housing-service unit for the current model. Bedrooms
are useful robustness checks and may be cleaner for some AHS supply objects, but
they are not the main calibration unit unless the model unit is changed.

## Candidate Replacement Map

| Weak current object | Candidate replacement | Empirical formula | Main source | Parameter block | Status |
|---|---|---|---|---|---|
| `prime_childless_owner_median_rooms` | Prime-age childless owner room-bin shares | \(E[1\{rooms \ge 6\}\mid O,n=0,30\le age\le55]\), plus adjacent bin shares | `code/data/mms_center_periphery/output_family_size_supply/acs_family_size_supply_cells.csv`; likely needs no-location collapse | \(\chi,\bar h_0\), owner housing menu | Build collapsed target and recheck sample |
| `prime_childless_owner_median_rooms` | Owner-renter room gap among childless prime-age households | \(E[h\mid O,n=0]-E[h\mid R,n=0]\) | `code/data/mms_center_periphery/output_family_size_supply/acs_parent_income_sorting_summary.csv` | \(\chi,\bar h_0\), tenure service wedge | Existing columns available; needs one-market collapse |
| `housing_user_cost_share` | Renter rent-to-income, by tenure/parent/age cell | \(12\,rent/hhincome\) for renters, with trimming and head weights | ACS/MMS scripts already compute rent burdens in `analyze_mms_fertility_moves.R` and `analyze_income_fertility_cross_section.R` | \(\alpha,\chi\), housing demand | Renter side feasible; sample must be aligned |
| `housing_user_cost_share` | Owner imputed user-cost-to-income | \(u^O VALUEH/hhincome\), with explicit user-cost convention \(u^O\) | ACS `VALUEH` path in `code/data/mms_center_periphery/validate_acs_home_value_scf.R` | \(\alpha,\chi\), owner user cost | Needs construction; do not use until audited |
| `old_age_own_rate` | Old nonhousing/liquid wealth-to-income | \(E[a/y\mid 65\le age\le75]\) and median analog | PSID wealth scripts under `code/data/psid_followup_mar2026/` | \(\theta_0,\beta\), old saving/bequest strength | Needs new old-age collapse |
| `old_age_parent_childless_gap` | Parent-childless old wealth gap | \(E[a/y\mid children>0,old]-E[a/y\mid children=0,old]\) | PSID completed-fertility and wealth variables | \(\theta_n,\theta_0\), child-linked bequests | Needs new old-age collapse |
| `liquid_wealth_to_income` | Age-profile slope of liquid wealth | \(E[a/y\mid45\le age\le55]-E[a/y\mid25\le age\le35]\) | `code/data/psid_followup_mar2026/output/entry_wealth_v1/` plus wealth scripts | \(\beta,b_0\), saving vs entry wealth | Partly available; needs common definition |
| `tfr` and `childless_rate` only | Fertility composition/parity share | Valid share by completed children, if household-unit convention matches model | PSID fertility-wealth outputs and ACS fertility cross-section | \(\psi_{\mathrm{child}},\kappa_n,\bar c_n\) | Candidate only; must avoid second-birth hazard confusion |
| External finance parameters if made internal | Down-payment shortfall or first-purchase LTV moment | Share short at target LTV/PTI or observed first-purchase down-payment share | `code/data/psid_followup_mar2026/output/dp_*` and `hacamo_table8_replica_v1` | \(\phi\), PTI/access parameters | Use only if finance parameters enter SMM |

## Immediate Construction Plan

1. Build an ACS/MMS one-market housing target CSV that collapses the household
   head sample across location but preserves tenure, child status, age window,
   rooms, room-bin shares, and renter rent-to-income.
2. Build a PSID old-age wealth target CSV for ages 65--75 using completed
   children, with separate liquid/nonhousing wealth and total net worth columns.
3. Build an owner user-cost diagnostic only after fixing the user-cost
   convention and validating ACS `VALUEH` against the intended price object.
4. Assemble a candidate target ledger with target value, standard error or
   sample size, source script, and exact formula for every candidate moment.
5. Rerun the sensitivity/Jacobian audit on any proposed replacement target set
   before launching another production-style search.

The next objective is measurement discipline, not fitting a new objective
quickly. A replacement target system is not valid until the data audit and
identification audit both pass.

## Extracted Values

Candidate replacement values were extracted later the same day and summarized
in `docs/model/intergen_candidate_replacement_targets_20260618.md`. Those
values remain candidate targets only; they do not supersede the audit
requirements above.
