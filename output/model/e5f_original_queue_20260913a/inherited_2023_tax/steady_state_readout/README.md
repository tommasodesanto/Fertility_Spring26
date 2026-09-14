# Property tax: verified stationary comparison

Slide: `output/pdf/property_tax_steady_states.pdf`.
Editable Beamer source: `latex/appendix_property_tax_steady_states.tex`.
Builder: `code/model/tools/build_e5f_stationary_tax_slide.py`.

The exercise compares the two already verified stationary endpoints at the
same post-decline fertility preference, with annual property taxes of 1% and 2%,
all property-tax revenue equally rebated to household heads, fixed payroll
tax 0.179 and balanced pensions. Both use the retained housing-supply curve and
original household-entry queues, with no immigration or population rescaling.
It does not require or assert that the currently computed transition has
converged, nor does it describe the reform's immediate 2023 effect. Selection
and uniqueness of the stationary equilibrium are not established by this slide.

`comparison.csv` contains raw levels, displayed indices and exact effects.
`verification.json` pins both terminal receipts and the fresh baseline audit.
The builder checks stationary convergence, population/queue stationarity,
household/fiscal gates, housing clearing and payroll consistency, then derives
every table cell directly from their saved `endpoint_reference` dictionaries.
No model solve is performed.

Higher tax gives 1.8051% more households, 0.7263% more aggregate housing,
1.0597% less housing per household, 1.1553% higher house prices, 25.5725% higher
unit rents, and 0.6265 percentage points lower ownership. These are differences
between the computed steady states, not cumulative or fixed-date effects.

Regenerate source and numbers with:

```sh
python code/model/tools/build_e5f_stationary_tax_slide.py
```

Compile the resulting TeX twice with `pdflatex`, directing build outputs to
this folder's `build/` subdirectory. Copy the PDF to the stated output path
and render it with `pdftoppm` for visual inspection. The PDF skill's artifact
marker script was unavailable on this host; generation and PDF verification
used local tools.
