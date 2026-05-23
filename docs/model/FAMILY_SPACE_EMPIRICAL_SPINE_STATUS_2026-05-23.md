# Family-Space Empirical Spine Status

This note maps the current empirical evidence to the proposed paper spine:

> Does the shortage of family-capable housing inside high-amenity urban
> locations reduce fertility and reallocate families across tenure and space?

The clean object is family-space access. Buying, renting, and moving are
access technologies. The paper should not start by assuming that the missing
object is rental housing or owner housing.

## Fact 1: Children Raise Housing Demand

Status: available.

Source:
`code/data/psid_followup_mar2026/output/no_location_family_space_packet_20260523/PSID_NO_LOCATION_FAMILY_SPACE_PACKET.md`

Current facts:

- First birth raises rooms by 0.664 rooms by +3 years.
- Second birth raises rooms by 0.704 rooms by +3 years.
- Moved-for-size rises from 0.076 at pre-birth event time -2 to 0.202 in
  post-birth years 0-3.
- Low-slack pre-birth renters expand rooms by 1.451 rooms and move at a high
  rate by +3.

Interpretation:

These facts validate the family-space transition. They do not identify
location-specific moves because this packet deliberately does not use
location-specific PSID.

## Fact 2: Family-Capable Units Are Scarce Or Costly In Central Locations

Status: available in both room-bin and bedroom-bin ACS/MMS versions.

Source:
`code/data/mms_center_periphery/output_family_size_supply/ACS_MMS_FAMILY_SIZE_SUPPLY_PACKET.md`

Bedroom source:
`code/data/mms_center_periphery/output_couillard_bedroom_supply/ACS_MMS_COUILLARD_BEDROOM_SUPPLY_PACKET.md`

Current room-bin facts:

- Size bins are `S_1_4`, `M_5_6`, and `L_7plus` rooms.
- Center housing is more tilted toward small units than the periphery. The
  center has a substantially higher small/family stock ratio relative to the
  periphery.
- Central family-capable rental share among central rentals is 0.344.
- The mean central M-size rent premium relative to small units is -0.027 log
  points in this room-bin construction, so the cleanest current evidence is
  more about physical stock/access than about a positive central room-price
  premium.

AHS absolute stock facts:

- The 2023 AHS has 90.7 million 3+ bedroom units.
- 70.6 million of those are owner units, 13.4 million are renter units, and
  4.4 million are vacant.
- 5.8 million renter households with children live in 0-2 bedroom units.

Bedroom ACS/MMS facts:

- Bedroom bins are `B_0_1`, `B_2`, and `B_3plus`.
- Weighted mean central 0-1/3+ bedroom stock scarcity is 0.914 log points.
- Parent households have 0.991 more bedrooms than childless households.
- Parents are 21.2 percentage points more likely to own than childless
  households in the young ACS/MMS sample.
- The mean central 3+ bedroom rent premium relative to 0-1 bedrooms is -0.020
  log points, so the bedroom evidence is again stronger on stock/access than
  on a positive renter-price premium.

Interpretation:

The aggregate claim cannot be "there are no family-sized units." The sharper
claim is that family-sized units are unevenly distributed across location,
tenure, and structure. Buying is empirically central: parenthood is strongly
associated with ownership, and family-sized units are often accessed through a
different ownership bundle rather than through a same-location rental upgrade.

## Fact 3: Parents Are Less Central Where Central Family Space Is Scarce

Status: newly stronger, but still descriptive.

Updated source:
`code/data/mms_center_periphery/analyze_family_size_supply_menu.R`

New output:
`code/data/mms_center_periphery/output_family_size_supply/acs_childless_parent_gap_vs_stock_scarcity.png`

Define the parent centrality gap as

$$
G_m = \Pr(C \mid \text{childless},m) - \Pr(C \mid \text{parent},m).
$$

Define central family-space stock scarcity as

$$
A_m =
\log\left(\frac{Q_{C,S,m}}{Q_{C,F,m}}\right)
-
\log\left(\frac{Q_{P,S,m}}{Q_{P,F,m}}\right),
$$

where \(S\) is 1-4 rooms and \(F\) is 5+ rooms in the current ACS/MMS run.
Larger \(A_m\) means the center is more tilted toward small units relative to
family-capable units than the periphery.

Current result:

- Weighted childless-minus-parent center-share gap: 0.139.
- Room-bin regression slope of \(G_m\) on \(A_m\): 0.2112, SE 0.0251.
- Room-bin regression slope of \(G_m\) on the central M-size price premium: 0.4220,
  SE 0.1349.
- Bedroom-bin regression slope of \(G_m\) on 3+ bedroom stock scarcity:
  0.1467, SE 0.0256.
- Bedroom-bin regression slope of \(G_m\) on the 3+ bedroom price premium:
  -0.0029, SE 0.0688.

Interpretation:

This is the strongest current bridge from Couillard to the paper. It says
parents are less central in exactly the metros where the central stock is most
tilted away from family-capable units. The price-premium evidence is weaker,
especially in the bedroom construction, which may reflect rent measurement,
selection into observed rental units, or the fact that family-space access is
often resolved through buying rather than renting. The evidence is not causal,
because family types and housing supply are jointly determined.

## Supply-Side Content

Status: partly available; causal supply evidence is missing.

What we have:

- AHS absolute stock counts by bedroom, tenure, and structure.
- ACS/MMS center/periphery stock shares by room bin and tenure.
- ACS/MMS center/periphery stock shares by bedroom bin and tenure.
- ACS/MMS metro-level association between parent centrality gaps and
  central family-space stock scarcity.
- A descriptive ACS/MMS early-to-late central 3+ bedroom stock-share change.

What this supports:

- The paper can say that the supply menu is consistent with a family-space
  access constraint.
- The paper cannot yet say how much of the fertility or sorting pattern is
  caused by exogenous housing supply restrictions.

What is missing:

- A sharper rent-gradient panel at finer geography or with Census summary
  tables, because current ACS/IPUMS bedroom rent premia are noisy and flat.
- A cleaner stock-growth panel with new construction or vintage, not just
  repeated cross-sectional stock shares.
- A causal supply-response test: high lagged family-space price premia should
  predict more family-sized stock growth in elastic markets, and weak growth
  in constrained markets.
- External supply-constraint shifters, such as Saiz land-unavailability or
  regulation/zoning measures, joined consistently to metro geographies.
- Location-specific PSID to test whether births cause center-to-periphery
  moves differentially in high family-space-scarcity metros.

## Immediate Empirical Priority

The next data object should sharpen the supply side rather than duplicate the
current cross-section:

1. Add new-construction or vintage cells by bedroom and location to distinguish
   inherited stock from active supply response.
2. Join external supply constraints, especially Saiz-style land unavailability
   and regulation/zoning proxies.
3. Build a finer rent-gradient panel from Census summary tables if a Census API
   key or cached API download is available.
4. Use location-specific PSID only for the dynamic claim:
   birth \( \rightarrow \) center-to-periphery movement by family-space
   scarcity.

## Bottom Line

The current facts support the paper spine:

1. Births raise rooms.
2. Family-sized housing is the relevant physical object.
3. Parents are less central where central family-capable stock is relatively
   scarce, in both room and bedroom versions.

The current facts do not yet identify the causal supply channel. The next
step is not more calibration. It is to measure whether family-capable stock is
slow to respond in constrained metros, and whether buying is the main observed
access technology when households move into family-sized space.
