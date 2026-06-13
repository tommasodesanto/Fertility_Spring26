# Handoff — bring Part 1 up to speed, then think (2026-06-13)

You are picking up an academic housing–fertility theory project mid-stream. Over the
last few days the **short note** (`latex/intergen_housing_fertility_short_note.tex`,
just renamed from `..._advisor_note...`) was revised heavily in a line-by-line pass with
the author. The **long version**, Part 1
(`latex/intergenerational_housing_fertility_part1.tex`), is now *behind* the short note
on a specific set of readability/framing decisions made today, and needs to be brought
into line — **carefully, without sweeping changes.**

There have been *many* small decisions, several of them author hand-edits that are final.
Your job is to adapt those decisions into Part 1, not to re-litigate, re-derive, or
rewrite anything beyond them. Read this whole prompt before touching a file.

---

## Step 0 — Load context first (do not skip)

1. `memory/AGENT_MEMORY.md`
2. `CALIBRATION_STATUS.md` (model status; orientation only — this task is paper-facing)
3. **`memory/project_part1_policy_stat_audit.md`** — this is the detailed running log.
   Read the dated UPDATE blocks from 2026-06-12 onward in full. The last few entries
   contain the exact "tomorrow's Part 1 mirror list" and every cross-document divergence
   you must preserve.
4. `CLAUDE.md` / `AGENTS.md` — project standards (LaTeX, minimal-change editing rules).
5. `docs/style/econ_writing_style_guide.md` — the house style for any prose you write.

Then read **both documents end to end**: the short note first (it reflects the agreed
decisions), then Part 1 (the target). Do not start editing until you can state, for each
item below, what the short note now says and what Part 1 currently says.

---

## Task 1 (primary) — sync Part 1 to today's decisions

The short note is the reference for *these specific passages only*. For each, bring
Part 1's parallel passage in line **with the substance of the decision**, adapted to
Part 1's fuller register — do **not** paste the note's terse phrasing into the long
version, and do not shorten Part 1's deeper treatment.

Today's note changes to mirror (all in §§ Equilibrium policies / Efficiency / Fertility):

1. **Child-price display.** The fertility FOC divided by the marginal utility of
   consumption is now read as an object: lead-in "the solved policies give the level of
   fertility; the first-order condition behind them shows where housing enters it,"
   then the display with two **underbraces** — left = "willingness to pay for a marginal
   child," right = "full price of a child" — followed by the clause that a binding access
   constraint "acts as a tax of $\kappa\,\zeta_i^m$ per child." Mirror the underbraces,
   the lead-in, and the per-child-tax clause into Part 1's version of this display.

2. **$F$ gloss.** At the first use of $\zeta_i^{O,F}$, the note now says "($F$ for
   financing)." Add the same one-time gloss in Part 1.

3. **Surplus display — term-by-term reading.** After the efficiency proposition the note
   adds a paragraph reading the expression: buyer's marginal value (effective price plus
   gap) minus incumbent's (net-of-tax sale price); the market price cancels; the two tax
   terms that remain are transfers the planner recovers, not resource costs. Part 1 has
   a parallel two-sided paragraph — make sure it carries this term-by-term reading. The
   proposition display also gained **underbraces** ("buyer's marginal value" /
   "incumbent's marginal value") and an explicit "$>0$"; mirror if Part 1's display
   doesn't already.

4. **"Tenure fixed but space moves" clause.** The note's proposition now says the
   variation "resizes two owner-occupied houses and converts no tenures." Add this
   clarifying clause to Part 1's proposition statement.

5. **Worked-example figure paragraph — compensated trade.** The note derives "frees X in
   goods" as a compensated trade: pay the incumbent at $q-\bar\ell$, charge the buyer at
   $q^O+\zeta^{O,F}$, both as well off, residual left over. Part 1's worked-example prose
   should phrase the surplus the same way (it currently may assert it). Keep all of
   Part 1's numbers.

6. **$\mathrm dn/\mathrm dH$ — two-force reading.** The note now explains the derivative
   as a race: extra space makes the child's space floor $\kappa$ easier to meet (fertility
   up) vs paying $q^m$ tightens the budget around the goods floor $\chi$ (fertility down);
   the savings margin weakens the second force because only $1/(1+k)$ of the housing bill
   falls on current spending, so the threshold scales up by $(1+k)$. Mirror this reading
   into Part 1 where the threshold $\chi\le(1+k)\kappa q^m/\alpha$ is discussed.

7. **$H$-setup / experiment framing.** The note's Fertility section now opens by naming
   the object — housing "pinned at the effective family-housing cap of Section 2" (rental
   size limit for a renter; shorter of owner size limit and collateral bound for a buyer)
   — writes $H$ for that cap's level, and states the experiment ("we relax the cap by a
   unit of space and study fertility; differentiating with respect to the cap level").
   Make sure Part 1's analogous opening names $H$ and the comparative static the same way.

8. **Policy levers as acts on parameters.** Before the $\mathrm dN/\mathrm dz$ display the
   note frames the policy: levers act on named primitives — a down-payment grant raises
   $a_i$, financial deregulation raises $\phi$ (both loosen the collateral bound
   $a_i\rho/((1-\phi)q)$), making family-sized units available to rent raises $\bar h^R$;
   "none of this is construction: the stock $\bar H$ is unchanged… access, not supply."
   The author also wrapped this in a **fertility-neutral framing** (footnote: not a planner
   targeting a desired fertility level, just one relaxing the housing constraint to ease
   misallocation and observing the fertility impact). Mirror the lever-by-parameter framing
   and the access-not-supply point into Part 1; carry the fertility-neutral framing if
   Part 1 lacks it.

### Divergences you must PRESERVE (do not "sync" these away)

- **$r^F$ stays in Part 1.** The note dropped the implementation-cost term $r_i^F$ from its
  surplus display and section intro (keeping it as one prose sentence). Part 1 **keeps the
  full $r^F$ treatment** — display, the $0/$positive$/\infty$ trichotomy, and the proof.
  Do **not** remove it from Part 1.
- **The full switching proposition stays in Part 1.** The note compressed tenure-switching
  to a two-sentence deferral and **removed the boundary numbers**. Part 1 keeps the full
  proposition (the $\Delta=\vartheta\gamma\log(q/(q-\bar\ell))$ continuation-gap
  trichotomy: price dominance, wedge dominance, fertility-jump decomposition), its proof,
  and the worked boundary numbers. Do **not** cut or shorten it. A light plain-English
  polish of its prose is fine; deleting content is not.
- **The efficiency proof stays in Part 1.** The note cut its proof and two cautions; Part 1
  keeps the full proof.
- **Part 1's worked example keeps all its numbers.** The note carries none of the boundary
  switching numbers; do not import that omission into Part 1.
- **Figure:** the note now uses two **standalone** panels
  (`example_misallocation_only.pdf` in Efficiency, `example_fertility_cap.pdf` in Fertility).
  Part 1 keeps the **combined** two-panel `example_misallocation.pdf`. The figure script
  `code/model/tools/make_part1_example_figure.py` already emits all three from one run; you
  should not need to touch it. Whether Part 1 *also* splits its figure is an open question —
  **ask the author before changing Part 1's figure**, don't decide it yourself.

### How to do the sync (guardrails)

- **Minimal change.** Edit only the named passages. Do not refactor surrounding prose,
  opening lines, footnotes, or author notes. Preserve existing wording unless the change
  is one of the items above. If a passage seems obsolete, **flag it; do not delete or
  rewrite it.**
- **Author hand-edits are final.** Both documents contain wording the author wrote by hand
  (e.g. the note's "If any constraint holds…", "We relax the cap…", the fertility-neutral
  footnote, the informal intro, several footnotes). Match the *decision*, not necessarily
  the exact words, and never "improve" the author's sentences in passing.
- **Notation is settled — do not change it.** $\beta$ = discount factor, $\vartheta$ =
  fertility weight; $k=\beta(1+\gamma+b)$; explicit composites $(c-\chi n)$, $(h-\kappa n)$;
  $q^m$ = the mode's effective price ($q^O=q+\bar\ell/(1+r)$, $q^R=q$); $\bar\ell$ = per-unit
  gains tax; $\varpi$ = inflation; $\pi$ = entry; deterministic tenure argmax; $\phi=0.80$.
  No balanced-growth path. If something looks inconsistent, check the audit memo before
  assuming it's a bug.
- **Verify before claiming done.** Compile Part 1 with `latexmk -pdf` from `latex/`; confirm
  zero new `Overfull` warnings and that the page count / cross-references are intact. Don't
  report success without the compile output.
- **Commit discipline.** `git status -sb` first. If the author left uncommitted edits,
  commit those *separately and first* with no Claude trailer, then make your changes in a
  clean commit. End commit messages with the standard `Co-Authored-By: Claude …` trailer.
  Push to GitHub `main` when the sync compiles. Never force-push.
- When the Part 1 sync is done, update `memory/project_part1_policy_stat_audit.md` with a
  dated entry recording what you mirrored and any item you deferred to the author.

---

## Task 2 (secondary) — think about the open threads

Breadth is welcome **in thinking**, not in editing. Surface and reason about the items
below; propose concretely; **implement only what the author explicitly approves.** Several
are flagged in the audit memo.

- **Paper file is stale.** The archived
  `latex/intergenerational_housing_fertility_paper.tex` (and any slide decks) still carry
  the *pre-inflation* capital-gains framing and older notation. It does not yet reflect the
  inflation device, the savings margin, the $q^m$ effective-price machinery, or the
  $\beta\leftrightarrow\vartheta$ swap. Scope what a sync would involve; do not start it
  without sign-off.
- **Deferred note item.** The author deferred an explanation of the old-incumbent
  state/intergenerational-books in the note ("complicated, we'll fix it next time"). Worth
  drafting a clean version for review.
- **Measurement.** How $\zeta^{O,F}$ is read from bundles and $\bar\ell$ from deeds /
  assessment records — the empirical bridge. There is empirical scaffolding under
  `code/data/`; think about what a measurement section would claim and what identifies it.
- **Quantitative model.** Calibration is the live workhorse (`code/model/`,
  `CALIBRATION_STATUS.md`). The theory note's frictions (down-payment constraint, gains-tax
  lock-in) need to map onto the quant model's wedges. Orient before proposing anything.
- Anything else you notice that strengthens the argument — but raise it, don't build it.

The author's standing mode of work: **they read and edit; you answer questions and make
changes only on an explicit go, after they've closed the file.** Default to discussion.
When in doubt about scope, ask.
