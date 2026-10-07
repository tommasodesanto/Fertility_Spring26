# Tenure segmentation theory note and slides (October 7, 2026)

Why the rental cap plus the parenthood space floor make the family-space option depend on
ownership access, and what that does to the cost of income risk for prospective parents.

- `tenure_segmentation_theory.tex` / `.pdf`: two-page note. Proposition 1 (exact within-period
  rental demand with children and the binding interval), Proposition 2 (envelope condition for the
  rental cap in the value gap of a birth), Proposition 3 (local convexity of the segmentation loss and
  the added risk penalty), then the full-model check and its verdict.
- `tenure_segmentation_slides.tex` / `.pdf`: two Beamer frames in the `latex/JMP_slides` style, one
  diagram each (desired space against expenditure with the cap and the floor; the loss-chord picture).
- `chatgpt_inputs/`: the three ChatGPT outputs (note PDF, slides PDF, slides pptx) copied from the
  Desktop on October 6, unchanged. The sources zip was not on the Desktop at copy time.

Base and evidence: chain 17 saved 2007 base (`output/model/production/2007/`, price 0.779), fixed
price throughout. Full-model check: `output/model/precautionary_2x2_20261006/round2_tables_segmentation_nakakuni.md`
(segmentation.py, cap8_readout.py), chain 17, fixed price, rebate on. October 5 segmentation test at
chain 11: `output/model/credit_mechanism_20261004/thread_credit_null/estate_a_chain11_seg/README.md`.

Compile either file with `pdflatex` twice; keep the `.aux`/`.log`/`.nav`/`.out`/`.snm`/`.toc`
byproducts out of the folder.
