# Economics Presentation Style Guide

Use this guide for paper talks, seminar decks, advisor updates, and quantitative-result
slides in this project. The intended reader is technically sophisticated but does not know
the codebase or the sequence of internal experiments.

## 1. Default look

- Follow the visual language of a conventional economics job-market or production-paper
  presentation: restrained, flat, and easy to scan.
- Use a single full-width composition for explanatory text and tables. Two columns are
  appropriate for two graphs or a genuine side-by-side comparison, not for prose plus a
  second wall of prose.
- Use conventional noun-phrase titles: `Toy economy`, `Competitive equilibrium`,
  `2023 calibration`, `Population dynamics`. Put the argument in the frame body.
- Use academic language even in advisor updates. Avoid deictic titles such as `Today`,
  `Now`, `Current status`, `What we established`, or `Where this leaves us` unless the
  date or sequence is itself substantive.
- Use `booktabs` tables without vertical rules, boxes, cards, dashboard elements, or
  decorative labels.
- Keep the main deck readable at presentation distance. Move robustness, derivations,
  implementation details, and alternative diagrams to the appendix.

## 2. Frame architecture

- One idea per frame.
- Keep bullets to one line whenever possible. If a bullet needs a paragraph, shorten it or
  split the argument across frames.
- A standard frame contains at most: one short framing sentence, one display/table/figure,
  and one interpretation sentence.
- Define each parameter or variable in plain English where it first appears. Do not use a
  generic state such as `$x$` as a substitute for naming the economically relevant states.
- Describe what the model does. State limitations only when they change the interpretation
  of the result; do not fill slides with lists of what the exercise does not do.
- Never expose workflow language such as `current reading`, `what we tried`, `the model
  shown to X`, `implementation status`, or internal experiment names. Translate it into the
  economic object an audience needs.
- Keep issue-ledger language out of outward-facing slides: `provisional`, `pending`,
  `not final`, `needs review`, run readiness, and implementation status belong in project
  notes unless they are necessary to interpret the economics shown on that frame.

## 3. Calibration and quantitative tables

- The default calibration table is exactly:

  `Moment | Target | Model`

- Do not add `Gap`, `Weight`, `Loss`, `Loss contribution`, standardized residuals, parameter
  counts, or objective diagnostics unless the author explicitly asks for them.
- Do not replace a conventional calibration table with a target-fit bar chart.
- If two or three misses deserve attention, bold those rows and name them in one sentence
  below the table. Do not add a diagnostic side panel.
- A parameter table defaults to `Parameter | Economic interpretation | Value`. Bounds,
  transforms, Jacobians, and near-bound diagnostics belong in the appendix unless requested.
- A policy table may show `Baseline | Reform | Change` when the change is the result being
  presented. Use percent changes for levels and percentage points for rates.
- Preserve units in row labels or column headers. Round for reading, while retaining enough
  precision to distinguish the compared objects.

## 4. Model and equilibrium slides

- Introduce the simplified economy before presenting its results: agents, timing, choices,
  prices, and market clearing.
- Define a steady state as a constant allocation, population/cohort structure, and price
  satisfying household optimality, demographic reproduction, and market clearing.
- When presenting a shock, separate three objects when they matter:
  1. the initial steady state;
  2. the impact or partial-equilibrium response;
  3. the demographic adjustment and new equilibrium.
- Use equations selectively. Every display gets one sentence before it and one economic
  interpretation after it. Derivations and proofs go to the appendix.

## 5. Figures and dynamics

- Use a figure only when it carries a mechanism or result more clearly than a table or one
  equation. Avoid diagrams of workflow, pipelines, or internal architecture.
- Label axes and equilibrium points economically. The primitive shock and the direction of
  movement should be visible without narration.
- For dynamics, show the impact point, any population/price trough, and the long-run point
  when those are the economic content. Choose illustrative parameters that make distinct
  stages visually distinct, without changing the qualitative economics.
- Keep the preferred mechanism graph in the main deck. Put alternative visualizations in
  the appendix and link to them only when useful during discussion.

## 6. Quantitative claims and caveats

- Distinguish targets, imposed inputs, untargeted validation series, and model outcomes.
- A path fixed to observed totals or age shares is an input, not a model fit. Say this once,
  plainly, where the path is introduced.
- Do not call a diagnostic closure a forecast, a finite transition a new steady state, or an
  impact calculation a complete policy experiment.
- Keep caveats short and economically substantive. Source notes belong in a small footer;
  code hashes, run identifiers, residual tolerances, and solver details do not belong in the
  visible main deck.
- Source footers should be standard citations or dataset names, not explanatory prose. If a
  measurement mismatch changes the economics, explain it once in the frame body.

## 7. Final check

Before delivering a deck:

1. Check every main frame against this guide.
2. Compile twice and require no errors, undefined references, or material overfull boxes.
3. Render the changed frames at full size and inspect them visually.
4. Confirm tables reproduce their source data and use consistent rounding.
5. Confirm the main deck contains economic content rather than a technical progress report.
