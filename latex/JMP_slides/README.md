# JMP Slides

The single continuing deck is `JMP_slides.tex`; its reader PDF is
`../../output/pdf/JMP_slides.pdf`.

The October 3 model exposition matches the canonical stationary production
engine and its post-interest chain-13 input snapshot: nonlinear child benefit,
physical parenthood room floor before the owner premium, zero-estate-normalized
bequest utility, income-inclusive purchase financing, sale and stayer debt rules,
stationary user-cost supply, zero household property-tax rebate and delayed
household entry. Estate-A and larger intended-birth menus are excluded experiments.
The author manuscript and mock are unchanged; remaining discrepancies are
recorded in `../README.md`.

The existing calibration, empirical and policy-result frames retain their
historical sources. They have not been refreshed for the current working model;
this model-section edit does not validate those results under the revised inputs.

Compile twice from the parent `latex/` directory so figure paths remain valid:

```sh
mkdir -p ../tmp/jmp_model_update/build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/jmp_model_update/build JMP_slides/JMP_slides.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/jmp_model_update/build JMP_slides/JMP_slides.tex
cp ../tmp/jmp_model_update/build/JMP_slides.pdf ../output/pdf/JMP_slides.pdf
```

Four missing appendix graphs were recovered as original embedded raster image
objects from the retained September PDF, with no changes to the plotted results.
Their provenance is in `assets/README.md`.

Verification: two final builds, no errors, undefined references or overfull boxes;
all changed model frames and both recovered-graph frames inspected at 150 dpi.
No model solve or recalibration was needed for this source-based update.
