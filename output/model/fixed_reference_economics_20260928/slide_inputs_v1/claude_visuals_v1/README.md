# Frozen-model visual storyboard and prototypes

**2007 stationary reference — block0506, September 28 verified export**.

The author requested Claude's visualization work while discussing economics with the lead. Claude Sonnet produced two supplemental prototypes from saved small CSVs. No model, lifecycle or equilibrium solves ran. All rendering and numeric checks ran on Torch. The 17 standard diagnostic plots remain unchanged; no manuscript or deck was edited.

## Reviewed prototypes

Use **reviewed_v3** only. The lead checked the source, all 37 plotted numerical rows against original CSVs, source/input/output hashes, and both actual PNGs. Torch **18821362** completed 0:0 in eight seconds (1 CPU, 4 GiB, three-minute cap). [Verification](reviewed_v3/lead_verification.json) and [manifest](reviewed_v3/actual_output/manifest.json) retain the evidence.

- [Contributors to additional first births](reviewed_v3/actual_output/birth_response_contributors.pdf): four-year fixed-price impact, grouped by inherited age, tenure and net financial wealth. Values are contributions per 1,000 **all initial households**, not responses per household within each group. First births contribute 83.4% of the total birth increase. Source: `../../credit_v1/summary_v1/impact_birth_decomposition.csv`.
- [Fixed-price cohort outcomes and stationary GE](reviewed_v3/actual_output/fixed_price_vs_ge.pdf): completed fertility, ownership and flow-weighted first-birth age. The fixed-price control and credit comparison use the same 262-node grid; the original exported 160-node baseline is not silently substituted. These separate conditional-cohort outcomes and stationary endpoints are not a transition path. Source: `../rendered_output_v2/borrowing_comparison.csv`.

The scripts and staged inputs are in [reviewed_v3](reviewed_v3/); [plotted_data.csv](reviewed_v3/actual_output/plotted_data.csv) records exact displayed units and levels. Frozen remote root: `/scratch/td2248/projects/fixed_reference_claude_visuals_20260929/reviewed_v3/`. Renderer SHA-256: `8bb3e19cf93a7c4c807dfe22618e1b002ac308ba36013cfbe58fba600938b581`. These are economically checked prototypes for choosing figures; presentation wording can be shortened when integrating the chosen figures into the existing deck.

## Four visual questions

1. **Who contributes to the credit response?** The first prototype answers this from occupied-state weights. Next pair contributions with responses per at-risk household to distinguish group size from sensitivity. This does not separately identify purchase-finance and renter-credit channels.
2. **Which changes survive equilibrium price adjustment?** The second prototype shows fertility returning to replacement while tenure and timing still change. Existing supply figures separately show household population and prices. Adjustment speed remains uncomputed.
3. **Where are reference credit limits binding?** An age/tenure view can use `../../constraints_v1/supplemental_constraints_by_age_tenure.csv` and its receipt. The recorded native credit component is the baseline limit before grid/death maxima, not the natural-solvency floor; component binding can overlap other limits. Purchase exclusion is a different object.
4. **More borrowers or greater debt per borrower?** `../../credit_v1/summary_v1/next_saving_debt.csv` distinguishes debt participation from net financial debt. Its `debt_per_branch_household` denominator includes all households in the branch, including zero debt. It is neither debt per indebted household nor new lending. Separating these denominators is necessary before interpreting the response.

Full 14-row fits, 31-row parameters and standard plots remain linked in [credit](../../credit_v1/README.md), [credit GE](../../credit_ge_v1/README.md) and [supply](../../supply_v1/README.md). Preferences including psi and other economic primitives are frozen; finite-grid and estate-counterparty caveats remain.

## Preserved review history

Initial v1 jobs18820670/18820725 and v2 job18820896 are **superseded, not approved for use**. Review corrected an unsupported83.7% annotation, percent versus percentage-point labels, cohort versus impact wording, debt denominators, the comparison's numerical baseline, and CSV quoting/units. Original outputs, prompts and Claude responses remain preserved. The v2 remote sibling directory is recorded in its prior receipt; no existing files were cleaned or deleted. V3 is the final reviewed version.
