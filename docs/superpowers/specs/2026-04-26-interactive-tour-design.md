# Interactive Tour of the DT Center–Periphery Model — Design Spec

**Date:** 2026-04-26
**Author:** Claude (with Tommaso)
**Status:** Draft for user review
**Output artifact:** `code/model/docs/tour.html`

## 1. Goal

Build a single, self-contained, Distill-style interactive guide to the discrete-time center–periphery lifecycle model. The page is rigorous (real math from the analytical notes, real code excerpts from the Python port, real model output from a calibrated solver run) and *fun* (interactive widgets that genuinely teach the model). It complements — does not replace — the existing reference doc at `code/model/docs/architecture.html`.

## 2. Audience and reading mode

Primary audience: the author (Tommaso) and any economics colleague who would otherwise have to read three .tex files plus 4400 lines of Python to understand the model. Secondary audience: future Claude/Codex sessions that want a single high-rigor pointer to "what the project is."

Reading mode: top-to-bottom narrative. Right-rail sticky nav for jumping. Each interactive widget is also rewarding when *not* touched — i.e., the page reads as a paper even if every slider is left at its default.

## 3. Scope

### In scope

- A single HTML file `tour.html` (~5,000–8,000 lines including embedded JSON and JS).
- Math rendered with KaTeX; code with Prism.js; charts with hand-rolled SVG + vanilla JS.
- Seven narrative chapters (see §5).
- Eight math↔code dual-view panels.
- Five interactive widgets:
  1. Period utility heatmap with parameter sliders.
  2. Tenure-choice decision tree with $b$ slider.
  3. Logit temperature visualizer (location and fertility).
  4. Backward-induction sweep animation.
  5. Policy explorer (real solver output, multiple slicing controls).
- One GE replay console.
- A Python data-prep script `tools/dump_explorer_data.py` and assembly script `tools/assemble_tour.py`.

### Out of scope

- In-browser solving of the model (no Pyodide). Any "live" computation is interpolation over pre-computed slices.
- Editing or replacing `architecture.html`.
- Writing new economics content beyond what is faithfully extracted from existing project sources.
- Authentication, server, deployment. The artifact is a static file.
- Mobile responsiveness beyond "doesn't visually break at narrow widths." The target is a desktop reading experience.

## 4. Visual and stylistic direction

- **Vibe:** Distill.pub × Bret Victor's explorable explanations, with rigor as the non-negotiable. No whimsy or decorative animation. Every interactive element must illuminate the model.
- **Typography:** serif body (e.g., *Source Serif 4* or *Crimson Pro* via Google Fonts), sans-serif headings (e.g., *Inter*), monospace code (*JetBrains Mono* or system mono).
- **Color:** family-resemblance to existing `architecture.html`. Accent `#2b5780`. Cream/off-white background `#fafaf7`. Code background `#f0eee8`. New additions: a soft callout green `#3a8060` for "interactive — try this" affordances, and a muted red `#a84040` for warnings/gotchas.
- **Layout:** wide main column (~720 px) with a left margin gutter (~220 px) for sidenotes, equation references, and small figures (Distill-style). Right-side sticky nav (~200 px). On narrow screens, the gutter collapses into inline notes.
- **Motion:** `prefers-reduced-motion` respected. Animations are short (≤2 s) and skippable.

## 5. Chapter map

### Chapter 0 — Cold open
Title, one-paragraph hook. Animated state-space sketch as the hero (six axes appear in sequence). Reading-time estimate. "What's interactive" badge listing the five widgets.

### Chapter 1 — The economic question
Why housing × fertility. 2–3 paragraphs adapted from `latex/draft_march.tex` introduction and `latex/Slides_AM_Oct8.tex`. One embedded empirical figure pulled from `code/model/benchmarks/diag11_lifecycle_python_x0.png` or a project-root event-study log if a richer one exists. Inline citations to Becker & Lewis (1973), Doepke–Kindermann (2019), Sommer–Sullivan–Verbrugge (2013), Gyourko–Mayer–Sinai (2013).

### Chapter 2 — The model
Math-heavy. Each subsection: rigorous derivation (faithful KaTeX from `Analytical_Notes.tex` / `Model_Draft_Feb.tex` / `draft_march.tex`), with one interactive widget where it adds rigor.

1. **State space** $(b, \mathrm{ten}, i, n, \mathrm{cs}, j)$ — interactive layered diagram; hover each axis for units, grid sizes, and source line in `parameters.py`.
2. **Period utility**
$$u(c, h, n) = \frac{\big[(c - \bar c_0 - \bar c_n n)^\alpha (h - \bar h_0 - \bar h_n n - \bar h_{\mathrm{jump}}\mathbf{1}[n\geq 1])^{1-\alpha}\big]^{1-\sigma} - 1}{1-\sigma}.$$
**Widget 1:** sliders for $\alpha, \sigma, \bar c_n, \bar h_{\mathrm{jump}}$ → SVG heatmap of $u(c, h)$ at chosen $n$.
3. **Bellman recursion** — full LaTeX, with the four-stage decomposition (savings → tenure → location → fertility) shown as a stack diagram.
4. **Tenure choice** — three branches (stay, sell-and-rent, buy/rebuy) with exact wealth bookkeeping. **Widget 2:** decision tree at sample $(j, n, \mathrm{cs})$ with $b$ slider; argmax flips visibly.
5. **Location logit** — $\Pr(i' \mid \cdot) = \exp((V_{i'} + E_{i'} - mc)/\kappa_{\mathrm{loc}}) / Z$. **Widget 3a:** $\kappa_{\mathrm{loc}}$ slider with probability bars.
6. **Fertility logit** — same primitive, parity domain. **Widget 3b:** $\kappa_{\mathrm{fert}}$ slider.
7. **Bequest** — $V_{\mathrm{bq}}(b, \mathrm{ten}, i, n) = \theta_0(n)\big((b + p_i H_{\mathrm{own}}[\mathrm{ten}-1] + \theta_1)^{1-\sigma}-1\big)/(1-\sigma)$. Small chart of $V_{\mathrm{bq}}$ vs $b$.

### Chapter 3 — From math to code (dual view)
The signature widget. Eight panels, each: equation on left, real Python excerpt on right (with file:line label), hoverable identifiers cross-highlight, expandable "why this code matches" footnote.

| Panel | Math object | Code location |
|-------|-------------|---------------|
| 3.1 | Period utility | `kernels.py:eval_renter_scalar`, `kernels.py:eval_owner_scalar` |
| 3.2 | Renter Bellman | `kernels.py:full_renter_block_kernel` |
| 3.3 | Owner Bellman | `kernels.py:full_owner_block_kernel` |
| 3.4 | Tenure max | `kernels.py:tenure_choice_kernel` |
| 3.5 | Location logsumexp | `kernels.py:location_logit_kernel` |
| 3.6 | Fertility logit | `solver.py` fertility block |
| 3.7 | Forward distribution | `kernels.py:forward_distribution_fast_kernel` |
| 3.8 | GE update | `solver.py:run_model_cp_dt` price/share update |

### Chapter 4 — How we solve it
- **4.1 Backward induction sweep** — **Widget 4:** SVG ticking ages 60→0; side panel shows which kernel is firing.
- **4.2 Howard policy iteration** — when `do_full` vs eval; cost-per-iteration drop chart from saved benchmark JSONs.
- **4.3 Forward distribution** — mass flowing forward; small static SVG sketch (no widget — the policy explorer in Ch. 5 covers the interactive part).
- **4.4 GE fixed-point dance** — **GE replay console:** real recorded GE run, $(p_C, p_P)$ trajectory in 2-D, equilibrium error chart, scrolling iteration log.

### Chapter 5 — What the model says (policy explorer)
**Widget 5:** sliders for $(j, n, \mathrm{cs}, \mathrm{tenure\ origin}, \mathrm{location})$. Charts (all from the inlined JSON dump):
- Savings policy $b'(b)$.
- Consumption policy $c(b)$.
- Tenure choice $\mathrm{tn}^*(b)$ — shaded regions.
- Location probabilities.
- Fertility probabilities (only at fertile ages).
- Lifecycle simulated paths from the saved forward distribution.

### Chapter 6 — How we calibrate
- SMM target table from `parameters.py` (target / model fit / weight).
- Identification arguments — which moment ties down which parameter — adapted from `docs/archive/calibration_docs_2026-05-07/CALIBRATION_PLAN_MERGED.md` and `docs/archive/root_notes_2026-05-07/CALIBRATION_LIT_REVIEW_MARCH2026.md`.
- Geography inversion explanation.
- Calibrated $\theta^*$ table.

### Chapter 7 — Reference
Compact: file-by-file (lifted from existing `architecture.html`), gotchas, symbol glossary linked back to Chapter 2, "where to look next" pointers.

## 6. Interactive widget designs

All widgets are vanilla JS + SVG. No frameworks, no D3. Each widget is a single `<div data-widget="...">` initialized by a small `widgets.js` block embedded in `<script>` at the bottom of `tour.html`.

### 6.1 Widget 1: utility heatmap
- Inputs: four range sliders ($\alpha, \sigma, \bar c_n, \bar h_{\mathrm{jump}}$) and a parity selector $n \in \{0, 1, 2, 3\}$.
- Output: 80×80 SVG heatmap of $u(c, h)$ over $(c, h) \in [\bar c_0+\bar c_n n + \epsilon, c_{\max}] \times [\bar h_0 + \bar h_n n + \bar h_{\mathrm{jump}}\mathbf{1}[n\geq 1] + \epsilon, h_{\max}]$. Iso-utility contours overlaid.
- Pure JS computation — no JSON dependence.

### 6.2 Widget 2: tenure-choice tree
- Inputs: range slider for $b$; selectors for $(j, n, \mathrm{cs}, i)$.
- Reads policy slice from JSON dump.
- Output: SVG tree showing three branches with their continuation values; the argmax branch glows. Budget constraints render alongside as text annotations.

### 6.3 Widget 3a / 3b: logit visualizers
- Single slider for $\kappa$.
- Output: probability bars over destinations (locations or parities). Pure analytical — uses representative $V$ values from the JSON dump.

### 6.4 Widget 4: backward induction sweep
- Play / pause / step buttons.
- Animation: an "age cursor" moves from $j=J-1$ down to $j=0$. At each frame, four boxes labeled "Vd / VH / VI / V" pulse in sequence. Side panel shows current age, time elapsed (mock), and which kernel ran.
- Pure cosmetic, but synced to actual kernel order.

### 6.5 Widget 5: policy explorer
- The flagship widget. Sliders/dropdowns for $(j, n, \mathrm{cs}, \mathrm{ten\_origin}, i)$.
- Six SVG line charts, all reading from the JSON.
- Tooltip on hover shows numeric values and the corresponding $b$ on the wealth grid.

### 6.6 GE replay console
- Plays back a recorded GE iteration trace.
- Three synced views: 2-D scatter of $(p_C, p_P)$ trajectory, line chart of equilibrium error vs iteration, scrolling text log of iteration messages.
- Play/pause/restart, no scrubbing (it's short).

## 7. Data pipeline

### 7.0 The canonical $\theta^\*$ for the run

Use the four calibrated values stored in user memory (`calibration_dt_april`):
$\beta = 0.947, \chi = 1.08, \kappa_{\mathrm{loc}} = 1.92, \bar h_{\mathrm{jump}} = 1.43.$
The remaining nine entries of `theta` (full 13-vector defined in `theta.py`) take their **defaults from `parameters.P_base`** at `setup_mode="benchmark"`. The data-dump script logs the full resolved $\theta$ to `meta` so the page can show "what was actually run."

If the .py defaults differ from what produced the memory's reported moments (own=65.2%, TFR=1.93, etc.), the page reports the *actual* moments from the dump, not the memorized ones. Ground truth is the JSON, not memory.

### 7.1 `tools/dump_explorer_data.py`

**New file.** Dumps a curated slice of one solver run to JSON.

```python
def dump_tour_data(theta, theta_names, out_path, setup_mode="benchmark"):
    sol, P, p_eq = solve_theta(theta, setup_mode=setup_mode, ge_trace=True)
    ...
    data = {
        "meta": {...theta, theta_names, setup_mode, p_eq, timestamp},
        "grids": {b_grid, H_own, age_years, parities, child_stages},
        "policies": {b_prime, c, tenure, loc_probs, fert_probs}, # sliced
        "ge_trace": [{iter, p_C, p_P, err, msg}, ...],
        "lifecycle": {...moments},
        "moments": {targets, model, weights},
    }
    json.dump(data, open(out_path, "w"))
```

The slicing strategy: do **not** dump the full 6-D arrays (~MB-scale). Instead dump a representative cross-section: all $(j, n, \mathrm{cs})$ tuples that the explorer reaches, restricted to the displayable subset (one renter origin + one owner origin per tenure). Target output size: 100–300 KB JSON.

### 7.2 `tools/assemble_tour.py`

**New file.** Reads `tour_template.html` (with a placeholder marker) and the JSON, inlines the JSON into a `<script id="tour-data" type="application/json">` block, writes `tour.html`.

This split lets us re-run the data dump without touching the HTML structure.

### 7.3 GE iteration trace

`run_model_cp_dt` already iterates internally. Either:
- (a) Add a lightweight `ge_trace` argument that appends `(iter, p, err)` per iteration to a list and returns it, or
- (b) Capture stdout from a verbose run and parse it.

**Choose (a).** Cleaner, deterministic. Adds ~5 lines to `solver.py`. Default off.

## 8. Source extraction strategy

- **Math:** open `Analytical_Notes.tex`, `Model_Draft_Feb.tex`, `draft_march.tex`. For each equation we need, copy the .tex source and convert to KaTeX-compatible. Keep the original .tex line in an HTML comment for re-sync (`<!-- src: latex/Analytical_Notes.tex:142 -->`).
- **Code:** read directly from `dt_cp_model/*.py`. Excerpt with file:start–end line labels. Truncate long blocks with `...` and a "[show full]" expander.
- **Empirical figure:** start with `code/model/benchmarks/diag11_lifecycle_python_x0.png` as the embedded image (base64-inline). If a better event-study figure exists in the project root (e.g., from `moving_first_birth_*` logs, but those are .log not .png — we'd need a real PNG), substitute later.

## 9. File organization

```
code/model/
├── docs/
│   ├── architecture.html         # existing — untouched
│   └── tour.html                 # NEW — output artifact
├── tools/
│   ├── dump_explorer_data.py     # NEW
│   ├── assemble_tour.py          # NEW
│   └── tour_template.html        # NEW — pre-inlining version
└── dt_cp_model/
    └── solver.py                 # +5 lines for ge_trace argument
```

The `tour_template.html` is the source of truth — readable, editable. `tour.html` is generated, gitignored or committed-as-built. (Recommendation: commit `tour.html` so colleagues can `git pull` and double-click.)

## 10. Acceptance criteria

The page is "done" when:

- [ ] `tour.html` opens by double-clicking from `file://`, no server required.
- [ ] All KaTeX equations render. No `[Math Processing Error]`.
- [ ] All Python code excerpts visible with syntax highlighting.
- [ ] Each of the eight Chapter 3 dual-view panels renders both math and code, and at least one identifier cross-highlights on hover.
- [ ] All five interactive widgets respond to input within 100 ms.
- [ ] The policy explorer's $b'(b)$ chart at $(j, n, \mathrm{cs}) = (10, 0, 0)$ matches the policy returned by `solve_theta` at the calibrated $\theta^*$ to within 1e-6 (i.e., we are not lying about what we display).
- [ ] The GE replay console plays through at least 10 iterations of an actual recorded run.
- [ ] The page passes a manual read-through end-to-end without broken anchors, missing figures, or unrendered math.
- [ ] File size of `tour.html` ≤ 5 MB (sanity bound).
- [ ] No errors in browser console on Safari and Chrome.

## 11. Risks and tradeoffs

- **JSON size growth.** If we want richer slicing in Widget 5, the JSON could balloon. Mitigation: sample the wealth grid (every other point) and document the downsampling.
- **KaTeX coverage.** A few macros from `Analytical_Notes.tex` may not render in KaTeX. Mitigation: rewrite to KaTeX-compatible during extraction; flag any non-trivial loss.
- **Math-source drift.** If `latex/*.tex` evolves later, the page becomes stale. Mitigation: HTML comments `<!-- src: ... -->` make re-sync straightforward; not auto-tracked.
- **CDN dependency.** KaTeX/Prism via CDN means offline opens won't render. Mitigation: optional follow-up to inline them; not in initial scope.
- **Scope creep.** The "expansive" framing risks an open-ended doc. Mitigation: this spec's chapter map and acceptance criteria are the contract; anything beyond goes in a follow-up.

## 12. Out of scope (explicitly)

- Live MATLAB-side comparison (the page is about the Python port).
- Calibration sliders that *re-solve* the model. Pre-computed only.
- Any modification to active solver or calibration code.
- New economics content beyond what is in existing project sources.
