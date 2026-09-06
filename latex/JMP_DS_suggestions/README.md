# JMP Draft Suggestions

This directory is the agent-writable staging area for material proposed for
`../JMP_DS_draft/`. Nothing here is part of the paper unless Tommaso copies it
into the protected draft by hand.

Use clear, task-specific filenames. Do not overwrite or synchronize files in
the protected draft directory.

## Compact theory results (5 September 2026)

[`simplified_olg_paper_core.tex`](simplified_olg_paper_core.tex) is a four-page
proposal: three pages of results to follow the existing household setup, with
the two existing figures, and one separate page assessing the transition
assumptions and stationary welfare. The short allocation proof, original
notation and conditional fertility condition are retained. The general welfare
derivative is derived and checked in the supporting proof's section 11.

PDF: [`output/pdf/simplified_olg_paper_core.pdf`](../../output/pdf/simplified_olg_paper_core.pdf).
The longer ten-page assessment is
[`simplified_olg_simple_assessment.tex`](simplified_olg_simple_assessment.tex),
with its [PDF](../../output/pdf/simplified_olg_simple_assessment.pdf).
Evidence and reproduction commands are indexed in
[`output/model/simplified_olg_amendments/README.md`](../../output/model/simplified_olg_amendments/README.md).
These are suggestions for discussion; no author decision or manuscript edit
is implied. Compile the compact source twice from the repository root:

```sh
mkdir -p tmp/pdfs/simplified_olg_paper_core output/pdf
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_paper_core latex/JMP_DS_suggestions/simplified_olg_paper_core.tex
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_paper_core latex/JMP_DS_suggestions/simplified_olg_paper_core.tex
cp tmp/pdfs/simplified_olg_paper_core/simplified_olg_paper_core.pdf output/pdf/simplified_olg_paper_core.pdf
```

## Independent theory review (4 September 2026)

[`simplified_olg_independent_review.tex`](simplified_olg_independent_review.tex)
preserves the full independent review of the simplified OLG theory: the model
map, formal claim audit, welfare claims A--C, transition assessment, economic
assessment, exposition, and proposed decision priorities. It records proposals
for discussion, not accepted author decisions or model amendments. Source labels
refer to the working-tree versions inspected on the review date.

The rendered report is
[`output/pdf/simplified_olg_independent_review.pdf`](../../output/pdf/simplified_olg_independent_review.pdf).
It includes a clickable contents page, landscape pages for the claim audit, and
a source-file guide. The standalone TeX source is sufficient to rebuild it.

From the repository root:

```sh
mkdir -p tmp/pdfs/simplified_olg_independent_review output/pdf
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_independent_review latex/JMP_DS_suggestions/simplified_olg_independent_review.tex
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_independent_review latex/JMP_DS_suggestions/simplified_olg_independent_review.tex
cp tmp/pdfs/simplified_olg_independent_review/simplified_olg_independent_review.pdf output/pdf/simplified_olg_independent_review.pdf
```

This review does not modify the author-controlled manuscript.
