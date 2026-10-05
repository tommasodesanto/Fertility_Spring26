# Household mechanisms

`frames.tex` contains the five mechanism frames included immediately before
`Policy` by the continuing `../JMP_slides.tex`. `mechanism_excerpt.tex` builds
the same frames by themselves; it is an excerpt, not a competing presentation.

Reader PDFs: `../../../output/pdf/JMP_mechanisms.pdf` and
`../../../output/pdf/JMP_slides.pdf`.

The editable plotting scripts, plotted CSVs, source identities, experiment
definitions and limitations are in
[`../../../output/model/credit_mechanism_20261004/mechanism_slides/README.md`](../../../output/model/credit_mechanism_20261004/mechanism_slides/README.md).
All five numerical figures use the same post-interest chain-13 reference
without Estate A. The pre-existing deck's historical calibration, transition
and policy frames retain their original evidence; these additions do not
refresh or reconcile their quantitative sources.

From the repository root, reproduce figures and both PDFs without a model solve:

```sh
bash output/model/credit_mechanism_20261004/mechanism_slides/regenerate.sh
```

Only the five selected figure PDFs are copied into `figures/`. Build products
and rendered verification images remain in the analysis packet's `build/`.
No manuscript, mock, frozen September 14 reference, or model source is changed.
