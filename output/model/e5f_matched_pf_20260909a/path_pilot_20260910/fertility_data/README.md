# Official fertility history for the timing pilot

Annual2007–2023 period TFR and live-birth counts come from the2015 and2023 NCHS final reports, with URLs, publication dates and downloaded-file hashes in source_manifest.json. The overlapping2010–2015 observations agree exactly. Revised2007 period TFR is2.1200; the author-selected initial model benchmark2.1 remains a separate assumption.

empirical_blocks.csv sums live births for2008–11,2012–15,2016–19 and2020–23. These align with the model decisions2007,2011,2015 and2019 under the explicit continuation-calendar convention. Its TFR field is an equal-year average of published rates, not an exposure-pooled rate. The builder checks four unique annual observations per block. Published age-specific rates are retained, but female exposure counts were not extracted and must not be reverse-engineered from rounded rates.

HH-3 provides annual March/CPS ASEC household stocks. Revised rows are selected and alternatives retained. These include all head ages; the extra birth-to-household columns are stock-proxy diagnostics, not literal annual exposures or counterparts to model starting-head rates. The primary pilot comparison uses aggregate birth-count indices. No empirical female-TFR/model-household-rate equivalence or new target contract is asserted.

Reproduce using build_receipt.py with pypdf and xlrd, placing the three official source files at the names listed in source_manifest.json beneath sources/. Source URLs provide retrieval locations. Downloaded source files remain local; the reader runtime is in task scratch rather than vendored in this packet. The lead independently checked all17annual dates, allfourblock totals/TFR averages and the three source hashes.
