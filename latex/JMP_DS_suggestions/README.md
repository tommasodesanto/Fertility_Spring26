# JMP Draft Suggestions

This directory is the agent-writable staging area for material proposed for
`../JMP_DS_draft/`. Nothing here is part of the paper unless Tommaso copies it
into the protected draft by hand.

Use clear, task-specific filenames. Do not overwrite or synchronize files in
the protected draft directory.

## Current full-note status (6 September 2026)

Tommaso has placed the new mixed-tenure transition and stationary welfare
analysis in the appendix as secondary material for now. The main priority is
the simple housing-allocation argument.

The four-page core below is a results extract, not the complete rewritten
note. The earlier full-model proposal is
[`simplified_olg_amendment_proposal.tex`](simplified_olg_amendment_proposal.tex),
with its [PDF](../../output/pdf/simplified_olg_amendment_proposal.pdf). It contains
the environment, household problems and equilibrium, but has not yet been
updated to incorporate all the later discussion. A consolidated revision of
that note remains to be done, preserving the author's notation and keeping
the new transition analysis in the appendix. The existing PDFs are unchanged
by this status clarification.

## Compact theory results (6 September 2026)

[`simplified_olg_paper_core.tex`](simplified_olg_paper_core.tex) is a four-page
proposal: three pages of results to follow the existing household setup, with
the allocation and analytical transition figures, and one separate page
assessing the transition assumptions and stationary welfare. The allocation
and conditional-fertility pages are unchanged. Figure 2 now uses a proved local
mixed-tenure transition with substantial renting and positive child goods
costs. The old limiting figure and its proof remain in the longer assessment.

[`simplified_olg_mixed_transition_proof.tex`](simplified_olg_mixed_transition_proof.tex)
is the new eight-page supporting appendix, with its
[PDF](../../output/pdf/simplified_olg_mixed_transition_proof.pdf). It derives the
exact original-equation map and initial-old boundary, proves local transitions
for a nondegenerate mixed family, and distinguishes demographic and welfare
signs. The general root and stationary signs hold at every positive taste
scale in that family; initial fertility and the boundary are certified on
the whole scale interval `[1,4]`. Finite-reform oscillations are proved near
scale `1`. Reform sizes remain local and unquantified.

PDF: [`output/pdf/simplified_olg_paper_core.pdf`](../../output/pdf/simplified_olg_paper_core.pdf).
The longer ten-page assessment is
[`simplified_olg_simple_assessment.tex`](simplified_olg_simple_assessment.tex),
with its [PDF](../../output/pdf/simplified_olg_simple_assessment.pdf).
Evidence and reproduction commands are indexed in
[`output/model/simplified_olg_amendments/README.md`](../../output/model/simplified_olg_amendments/README.md).
These are proposed results for discussion; their development did not change
planner permissions or the manuscript. The subsequent author placement
decision is recorded above. Compile the compact source twice from the repository root:

```sh
mkdir -p tmp/pdfs/simplified_olg_paper_core output/pdf
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_paper_core latex/JMP_DS_suggestions/simplified_olg_paper_core.tex
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=tmp/pdfs/simplified_olg_paper_core latex/JMP_DS_suggestions/simplified_olg_paper_core.tex
cp tmp/pdfs/simplified_olg_paper_core/simplified_olg_paper_core.pdf output/pdf/simplified_olg_paper_core.pdf
```

Compile the appendix similarly, using
`tmp/pdfs/simplified_olg_mixed_transition` for all build products and
`simplified_olg_mixed_transition_proof.tex` as the source. Compile twice and
copy only the final PDF to `output/pdf/`. Neither build writes into the
protected author draft.

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
