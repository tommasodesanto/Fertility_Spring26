# JMP Draft Suggestions

This directory is the agent-writable staging area for material proposed for
`../JMP_DS_draft/`. Nothing here is part of the paper unless Tommaso copies it
into the protected draft by hand.

Use clear, task-specific filenames. Do not overwrite or synchronize files in
the protected draft directory.

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
