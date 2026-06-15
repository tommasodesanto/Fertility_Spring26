# Part 1 Audit Disposition, June 15 Draft

This note records the June 15 audit pass on the analytical housing-fertility drafts.

## Files touched

- `latex/intergenerational_housing_fertility_part1.tex`
- `latex/intergen_housing_fertility_short_note.tex`
- `code/model/tools/make_part1_example_figure.py`
- regenerated `latex/figures/example_misallocation.pdf`, `example_misallocation_only.pdf`, and `example_fertility_cap.pdf`

Pre-audit copies are in `latex/archive/part1_june15_pre_audit/`. A ChatGPT Pro attachment bundle was refreshed on the Desktop at `/Users/tommasodesanto/Desktop/part1_chatgpt_pro_review_2026-06-15`.

## Fixed in the June 15 draft

- Added the maintained interior-saver/interior-old-retention qualification behind \(q^O=q+\bar\ell/(1+r)\). For constrained savers, the text now points to the type-specific continuation price \(q_i^O=q+\beta\bar\ell/(\lambda_i c_{j'}^O)\).
- Qualified the owner financing-gap measurement so \(\zeta_i^{O,F}\) is tied to slack owner size caps plus the same interior-margin assumptions.
- Clarified the planner entry object: the displayed \(V_i^P\) is a current-period object; the local results hold entry fixed, while a lifetime planner-entry problem would add discounted old-age value.
- Repaired the efficiency proposition language in both drafts. A positive welfare statement now requires a feasible local reallocation away from real size-limit corners, and the implementation cost \(r_i^F\) is explicit in the short note.
- Replaced “planner recovers transfers” language with the cleaner compensated-allocation reading: tax terms are not real resource costs and can offset transfers.
- Relabeled the worked example as a local/cross-sectional allocation rather than a stationary example. The old incumbent state is taken as given; the buyer anticipates her own future tax through \(q^O\). All numerical values were preserved.
- Tightened the policy-access paragraph: down-payment grants raise \(a_i\), deregulation raises \(\phi\), and rental family-space access raises \(\bar h^R\). Holding \(\bar H\) fixed, the rental lever means reclassification, conversion, or reallocation of existing units, with real conversion costs in \(\Xi^\Pi\).
- Expanded \(\Theta_i^m\) at definition to \(\partial W_i^m/\partial H_i^m=\lambda_i^m\zeta_i^m=\zeta_i^m/(c_i^m-\chi n_i^m)\), and described the last term as the normalized access gap.
- Added the tenure-switching orientation convention: the normal points toward the owner region \(\{\Delta_i>0\}\).
- Qualified the rent-own switching prose: \(q^O>q^R\) alone does not imply an upsize or fertility jump without the proposition's dominance and part (iii) conditions.
- Replaced “baseline” for \(q^O=q^R\) with “no-tax/no-tenure-price-wedge benchmark.”
- Cleaned the figure script notation: the fertility weight is `vartheta`, the total bundle gap is `gap_over_q`, the financing component is `zeta_OF`, and output paths are repo-relative.

## Follow-up Audit Patch

- Strengthened the old/saving regularity condition: the text now imposes \(0\le\bar\ell<q\), positive old liquidation wealth, \(0<h^O_{j'}<H^0_{j'}\), interior saving, and interior old-stage choices wherever old-renter FOCs or the switching constant are used.
- Added the balance-sheet clarification that the down-payment inequality is a collateral/access test, not a second goods expenditure: home equity and mortgage debt net out in \(a_{j'}=(1+r)a_i'-qH^0_{j'}\).
- Split the log warm-glow case from the no-bequest case: \(b>0\) uses \(\mathcal B(e)=b\log e\) with an interior bequest margin; \(b=0\) uses \(\mathcal B\equiv0\) and allows \(e=0\).
- Restricted the \(q^O<q^R\) tenure-switching discussion and corollaries: those signs are not implied unless the alternative environment also satisfies the proposition's dear-mode-constant condition, or a separate lemma allows \(\Delta<0\).
- Recovered and documented the switching-boundary calculation behind the \(0.007\) fertility jump: \(w=0.75304\), \(H^A=0.15\), \(H^B=0.20\), \(q^A=1\), \(q^B=1.14815\), \(k=0.26445\), \(\Delta=0.01967\), \(n_A=0.14025\), \(n_B=0.14769\).
- Added the short-note nonhousing-estate sentence, baseline outside-fertility normalization, access-only grant caveat, and precise normalized-gap policy ranking.
- Made the figure script robust to a flat copied bundle: it searches upward for `latex/figures` and otherwise writes to `Path.cwd()/figures`.

## Final June 15 Cleanup

- Corrected the competitive-equilibrium/planner fertility comparison to use the mode effective price \(q^{CE,m}\): \(q^{CE,R}=q^{CE}\) and, under Assumption 1, \(q^{CE,O}=q^{CE}+\bar\ell/(1+r)\). This avoids silently dropping the owner anticipated-tax wedge.
- Changed the Stone--Geary equal-scaling discussion in Section 9 to the mode-price condition \(\chi=\kappa q^m/\alpha\), while preserving the earlier market-price benchmark as the renter/no-tax version.
- Made Proposition 3 part (iii) self-contained by displaying the within-level-set versus level-shift inequality in the proposition statement.
- Added fixed-\(\bar\ell\) qualifications to the property-tax capitalization channel: if \(\bar\ell=\tau^gP\varpi/(1+\varpi)\) moves with \(P\), the same policy also changes \(q^O\) and the old-retention wedge.
- Removed the short note's unconditional \(>0\) from the surplus display; positive surplus now holds iff the real implementation cost is below the value gap.
- Tightened the short-note policy derivative to the fixed-price, fixed-tenure, fixed-active-set object \(\left.\dd N/\dd z\right|_{q,o,\mathcal A}\), and aligned the tenure-switching closing sentence with Part 1.
- Set the figure script's savings load from the displayed paper parameters, `beta_disc = 0.1763` and `k = beta_disc * (1.0 + gamma)`, leaving the printed numerical check as validation rather than reverse engineering.
- Verification after this pass regenerated the three example figures, recompiled Part 1 and the short note with `latexmk -pdf -interaction=nonstopmode -halt-on-error`, and found no Overfull boxes, Underfull boxes, undefined references, citation warnings, or rerun-needed warnings. The only warning-pattern log hits are the `rerunfilecheck` package-identification lines.

## Preserved or deferred

- Part 1 still keeps the full \(r_i^F\) trichotomy, the full switching proposition/proof/numbers, and the efficiency proof.
- Part 1 still uses the combined figure `example_misallocation.pdf`; the short note still uses two standalone panels. No decision was made to split Part 1's figure.
- The switching-boundary numbers are retained and now documented in a Part 1 footnote. They are a separate boundary calculation, not the worked local allocation above.
- The old-incumbent state / intergenerational-books explanation remains a deeper exposition item, not solved by this pass.
- The stale full-paper file `latex/intergenerational_housing_fertility_paper.tex` and slide decks still need a separate sync. They likely carry the pre-inflation gains framing, old notation, no savings margin, and no \(q^m\).
- The quantitative model still needs a mapping pass from theory frictions to code wedges: capital-gains/step-up, estate-to-entry adding-up, PTI or debt-service constraints, old retention, and fertility timing remain the main alignment gaps.
- Measurement bridge remains to be scaffolded: \(\zeta^{O,F}\) requires bundles, preference parameters, cap slackness, and interior-margin assumptions; \(\bar\ell\) requires basis, current value, exclusions, marginal tax rates, and sale-versus-bequest treatment.

## Verification

- Ran `code/model/.venv/bin/python code/model/tools/make_part1_example_figure.py`.
- Ran `latexmk -pdf -interaction=nonstopmode -halt-on-error intergenerational_housing_fertility_part1.tex` from `latex/`.
- Ran `latexmk -pdf -interaction=nonstopmode -halt-on-error intergen_housing_fertility_short_note.tex` from `latex/`.
- Final logs have no Overfull boxes, Underfull boxes, undefined references, citation warnings, or rerun-needed warnings. The only `rg` hit on `Warning|Rerun|Overfull|Underfull|undefined` is the package-identification line for `rerunfilecheck`.
