# Round-2 review prompt — Part I (2026-06-12, post-fix verification)

Attach the updated `intergenerational_housing_fertility_part1.tex` (+ PDF) and
`intergen_housing_fertility_advisor_note.tex` (+ PDF), then paste everything below the line.

---

This is the revised version of the two documents you reviewed. Your findings were dispositioned as follows; your job now is to verify the fixes are mathematically correct as implemented, not merely present, and to sweep for anything the fixes themselves broke. Same standard as before: derive, don't pattern-match.

**Implemented — verify each one.**

1. (your A1/B3/D4) Bequest adding-up now reads $\bar M\int[\cdots]\,\mathrm dG_Y=\int e_j\,\mathrm d(F_O+F_R)(j)$ and $\Gamma_Y$ conditions on $(F_O,F_R,e)$ at all six occurrence sites. Check the accounting is now closed against goods feasibility, which integrates $(c_j^O+e_j-a_j)$ over $F_O+F_R$.
2. (A2/A5) The assumption is now three-part: (i) kernel draws of $(g_{j'},\ell_{j'},\mathcal B_{j'})$ invariant in $h_i$ beyond the displayed books; (ii) interior old retention $h^O_{j'}<H^0_{j'}$; (iii) $\mathbb E[(g_{j'}-\ell_{j'})/c^O_{j'}]=0$. Verify (i)+(ii) are exactly sufficient for $\partial V^O/\partial H^0=(q-\ell_{j'})/c^O_{j'}$ and that no further condition is needed (e.g., differentiability of the levy at the margin, or the old budget's dependence on $H^0$ through the levy term $-\ell_j(H^0_j-h^O_j)$ — note $\partial/\partial H^0$ of that term is $-\ell_j$, which is already inside the stated envelope; confirm).
3. (A3/C3) Switching proof: ray prose now $Y=\alpha X/(1+k)$ and region $(1+k)Y>\alpha X$ throughout; the reduction paragraph adds that the concentrated continuation contributes a common constant across modes when induced old choices are interior in both, the continuation being $(1+\gamma+\mathbb E b_{j'})\log((1+r)a')$ plus a draw constant. Check this constant really is mode-invariant under the stated hypotheses (including that the $\mathcal B_{j'}$ draw may differ across $\Gamma_O$ and $\Gamma_R$ — is the stated hypothesis strong enough, or does it implicitly require the warm-glow draw distribution to be tenure-invariant?).
4. (A4/D1/D2/D3) Advisor note: savings rule now carries "interior old choices, repricing degenerate at zero, $k$ common across households"; neutrality stated as the exact MU-weighted condition; surplus proposition adds "owner menu slack." Check consistency with the canonical statements.
5. (A6/D5) $\vartheta=0.2164$ so $k=\vartheta(1+\gamma)=0.325$ to displayed precision. Recheck the example chain at this precision.
6. (B1) Units convention added at the books: old housing wealth in per-period goods units, a carried unit worth one period of services $q$, the shell recirculating through the stock, $P$ confined to the down-payment requirement. Does this convention make the old budget, the books, and the down-payment constraint mutually consistent, or do you see a remaining dimensional objection?
7. (B2, B4–B8, C2, C4, E1–E3) $\Gamma_R$ in the primitives lists; density renamed $f(i)$; old-stage superscript note; $\mathcal A$ defined; financing domains in the note; $B_i$ as gross assigned bequest taxed on receipt; "savings" in the converse conditioning; normalized-gap targeting sentence; prose trims.

**Not implemented — do not resubmit these; rebut only with a concrete derivation if you believe the rejection is wrong.**

- C1: the three owner constraints are linear upper bounds on $h$, so the feasible set is exactly $h\le\min\{\bar h^O,A_i\rho/((1-\phi)q),\psi y_i/q\}$ — the collapse is exact, not local. Derivative caveats (unique binding component, zero-measure kinks) are stated where derivatives are taken.
- C5: the sentence already contains the fixed-$q$ qualifier "at a given user cost."
- E4, E5: register choices ("Three things come out"; the three "old-side amplifies" statements scope different sections and are not verbatim).

**Standing do-not-relitigate list (unchanged):** planner-FOC entry terms cancel via a common factor; no balanced growth path; deterministic tenure argmax; $\phi=0.80$; explicit composites notation; single floorspace market with $q^R=q^O=q$ baseline; entry/tenure/transition margins outside the planner comparison; existence/uniqueness not claimed.

**Deliverables.** (a) For items 1–7: a one-line verdict each — correctly fixed, or incorrectly fixed with the derivation showing why. (b) Any new errors introduced by the revision (changed displays, broken cross-references, inconsistencies between the two documents created by the edits). (c) A final judgment: is the pair of documents ready to send to a senior advisor, and if not, the shortest list of blockers. No wholesale rewrites; no restating the model.
