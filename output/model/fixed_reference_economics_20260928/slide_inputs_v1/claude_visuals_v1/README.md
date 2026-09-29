# Frozen-model visual storyboard and prototypes

**2007 stationary reference — block0506, September 28 verified export**.

The author requested Claude's visualization work while discussing economics with the lead. Claude Sonnet produced two supplemental prototypes from saved small CSVs. No model, lifecycle or equilibrium solves ran. All rendering and numeric checks ran on Torch. The 17 standard diagnostic plots remain unchanged; no manuscript or deck was edited.

## Reviewed prototypes

Use **reviewed_v5**. Following the author's readability feedback, the captions use plain language and keep grid bookkeeping here. Torch **18827963** completed 0:0 in nine seconds (1 CPU, 4 GiB, three-minute cap, zero model solves). The lead verified that all 37 numerical values and units match the scientifically reviewed v3. Both actual PNGs were visually inspected; v5 PNGs are byte-identical to the inspected v4. [Verification](reviewed_v5/lead_verification.json) and [manifest](reviewed_v5/actual_output/manifest.json) retain the evidence.

- [Contributors to additional first births](reviewed_v5/actual_output/birth_response_contributors.pdf): four-year fixed-price impact, grouped by inherited age, tenure and net financial wealth. Values are contributions per 1,000 **all initial households**, not responses per household within each group. First births contribute 83.4% of the total birth increase. Source: `../../credit_v1/summary_v1/impact_birth_decomposition.csv`.
- [Fixed-price cohort outcomes and stationary GE](reviewed_v5/actual_output/fixed_price_vs_ge.pdf): completed fertility, ownership and flow-weighted first-birth age. The fixed-price control and credit comparison use the same 262-node grid; the original exported 160-node baseline is not silently substituted. These separate conditional-cohort outcomes and stationary endpoints are not a transition path. Source: `../rendered_output_v2/borrowing_comparison.csv`.

The scripts and staged inputs are in [reviewed_v5](reviewed_v5/); [plotted_data.csv](reviewed_v5/actual_output/plotted_data.csv) records exact displayed units and levels. Frozen remote root: `/scratch/td2248/projects/fixed_reference_claude_visuals_20260929/reviewed_v5/`. Renderer SHA-256: `9ac4904ca6ad088b3984d33bc0fa5dd84c847aa5e2d5fc9c611973f1369bce94`. These are economically checked prototypes for choosing figures. The next priority is standard policy curves against financial wealth: saving with the credit floor, housing choice, and birth probability, with other household circumstances held fixed and occupied wealth ranges marked.

## Four visual questions

1. **Who contributes to the credit response?** The first prototype answers this from occupied-state weights. Next pair contributions with responses per at-risk household to distinguish group size from sensitivity. This does not separately identify purchase-finance and renter-credit channels.
2. **Which changes survive equilibrium price adjustment?** The second prototype shows fertility returning to replacement while tenure and timing still change. Existing supply figures separately show household population and prices. Adjustment speed remains uncomputed.
3. **Where are reference credit limits binding?** An age/tenure view can use `../../constraints_v1/supplemental_constraints_by_age_tenure.csv` and its receipt. The recorded native credit component is the baseline limit before grid/death maxima, not the natural-solvency floor; component binding can overlap other limits. Purchase exclusion is a different object.
4. **More borrowers or greater debt per borrower?** `../../credit_v1/summary_v1/next_saving_debt.csv` distinguishes debt participation from net financial debt. Its `debt_per_branch_household` denominator includes all households in the branch, including zero debt. It is neither debt per indebted household nor new lending. Separating these denominators is necessary before interpreting the response.

Full 14-row fits, 31-row parameters and standard plots remain linked in [credit](../../credit_v1/README.md), [credit GE](../../credit_ge_v1/README.md) and [supply](../../supply_v1/README.md). Preferences including psi and other economic primitives are frozen; finite-grid and estate-counterparty caveats remain.

## Preserved review history

Initial v1 jobs18820670/18820725 and v2 job18820896 are **superseded, not approved for use**. Review corrected an unsupported83.7% annotation, percent versus percentage-point labels, cohort versus impact wording, debt denominators, the comparison's numerical baseline, and CSV quoting/units. Original outputs, prompts and Claude responses remain preserved. The v2 remote sibling directory is recorded in its prior receipt; no existing files were cleaned or deleted. V3 (18821362) remains the numerical comparator. V4 (18827443) simplified captions but introduced three empty CSV unit fields through a renamed display label. V5 repaired that unit mapping in new immutable source and regenerated outputs. Text labels intentionally differ from v3; all numerical values and units match. No output CSV was substituted after rendering.
