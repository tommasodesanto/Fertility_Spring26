# Mechanism slides: review preview

Five-frame preview of Astra's October 4 economics proposal, intended for author
review before integration immediately before Policy in the continuing
`latex/JMP_slides/JMP_slides.tex`. This is not a replacement presentation.

`mechanism_preview.tex` is the editable source; `mechanism_preview.pdf` is the
reader artifact. The first three frames use working chain-13 fixed-price
diagnostics, not the September historical results in the surrounding main deck
or the experimental Estate-A 2023 transition. Frames four and five define
proposed comparisons and contain no estimated results. No model solve was run
to produce this preview.

Sources: `../census/README.md`, `../DECISION_MEMO.md`,
`../diagnostics/README.md`, and `../diagnostics/joint_budget/README.md`.
Slide one describes space admissibility, not budget feasibility. The proposed
credit benchmark still requires an explicit debt/solvency contract. The proposed
fixed-stock transition still requires a quantity and terminal-closure contract.

Rebuild from this directory with `pdflatex -interaction=nonstopmode
-halt-on-error mechanism_preview.tex` twice. The native Codex LaTeX editor also
compiles this self-contained source. All five frames were rendered and reviewed.
