# Part I handoff — Claude session 2026-06-11

Scope: the Part I analytical note `latex/intergenerational_housing_fertility_part1.tex`
(+ slides `latex/intergenerational_housing_fertility_part1_slides_v2.tex`). A parallel
ChatGPT handoff and audit exist and will be pasted alongside this; I have NOT seen them.
The new session should adjudicate that audit against this state, not start fresh.

## Where the draft stands (HEAD = a5508ce, working copy clean)

The note is now a Coven-style single-floorspace economy with Menzio-style
environment-first exposition:

- **§1 Environment**: labeled primitive paragraphs (Agents, Preferences, Housing,
  Financing, Property-tax assessments, Entry). Stone–Geary child inputs
  ($c_i = C_i - \chi_i n_i$, $s_i = h_i - \kappa n_i$, $u = \log c + \alpha\log s + \beta\log n$);
  full-price-of-a-child identity $c + qs + (\chi+\kappa q)n = y$. Tenure is a size-menu
  distinction only: $h \le \bar h^R$ (rent) vs $h \le \bar h^O$ (own), $\bar h^R < \bar h^O$;
  these are occupancy menus, not separate stocks. One market: all floorspace clears at
  user cost $q = \rho_p P$, $\rho_p = r + \delta^O + \tau^p$; fixed stock $\bar H$;
  $H_{\text{pass}} = \bar H - \int H_j^0\,dF_O \ge 0$ held by passive landlords.
  Construction / tenure-specific supply deferred to a footnote as extensions.
- **§2–4**: young problem (owner constraints $(1-\phi)Ph \le A_i(x)$, $qh \le \psi y_i$;
  effective cap $H_i^O$; $A_i(x) = a_i + xB_i + T_i^b(x)$ pledgeable-only), old problem
  (incumbent assessment discount $L_j^\tau = \tau^p P - \tilde\tau_j\tilde P_j \in [0,q)$,
  "subsidy to staying put"), entry (logit with scale $\kappa_E$).
- **§5 CE**: single demand $H^D(q,\tau^p,x) = \bar H$. Note $\tau^p$ is now an explicit
  argument of $W^O, W^Y, \pi^E, H^D$ (added in a5508ce): given $q$, $\tau^p$ moves
  $P = q/\rho_p$ and hence the DP constraint.
- **§6–8 Planner and efficiency**: single feasibility with multiplier $\Lambda Q$;
  financial constraints are *implementation* constraints, not resource constraints.
  a5508ce REPLACED the old efficiency proposition with "Planner–equilibrium comparison,
  local": pairwise reallocation surplus
  $(q + \zeta_i^{O,F} - r_i^F) - (q - L_j^\tau) = \zeta_i^{O,F} + L_j^\tau - r_i^F$,
  where $r_i^F \ge 0$ is the marginal real resource cost of relaxing $i$'s financing
  (pure-finance case $r_i^F = 0$); inefficiency iff some feasible pair has positive
  surplus; converse direction stated. **I have not audited this new proposition** —
  it likely interacts with the ChatGPT audit.
- **§9 Fertility and policy**: $dn/dH$ comparative static, sign condition
  $\chi \le \kappa q/\alpha$ (equal scaling $\chi = \kappa q/\alpha$ as benchmark);
  local policy decomposition
  $dN/dz|_{q,o,A} = \bar M \sum_m \int_{B_m} [\pi^E \Psi^m + (n^m - \bar n^E)\pi^E(1-\pi^E)/\kappa_E \cdot \Theta^m] P^m\,dG$
  with $\Theta = \zeta/c$; corollary $\Psi = (\kappa/\alpha)(\zeta/c)\,\Omega$ under
  equal scaling ($\Psi \propto \zeta$); tenure-switching subsection (proposition in
  branch labels A/B; boundary flux $\bar M \int_{I(z)} \pi^E (n^O - n^R)\,\partial_z\Delta/\|\nabla\Delta\|\,g\,dS$;
  trichotomy in $q^O$ vs $q^R$, with parity — the single-market baseline — giving a
  zero switching jump, so the no-switching formula is exact at the baseline).

## Math verified in this session (do not re-derive, do not let an audit "fix" these)

- Switching machinery: $X = \chi/c$, $Y = \alpha\kappa/s$, $n = \beta/(X+Y)$,
  $U = K - G(X,Y)$ with $G = \log X + \alpha\log Y + \beta\log(X+Y)$; level sets are
  globally decreasing graphs; $G_{XY}G_Y - G_{YY}G_X = \alpha/(XY^2) + \beta/(XT^2) + \alpha\beta X/(Y^2T^2) > 0$
  ($T = X+Y$); single crossing of $Y = \alpha X$; fertility single-peaked along
  indifference curves at $s/c = \kappa/\chi$.
- Switching proposition: $q^A < q^B$ + indifference $\Rightarrow$ the cheap-side cap
  binds strictly and $h_B > h_A$; under $\chi \le \kappa q^A/\alpha$ also $c_B < c_A$,
  $s_B > s_A$, $n_B > n_A$. (The $h_B > h_A$ part is ES-free; see open item 3.)
- Planner FOC entry-response terms cancel via a common factor (an external "Patch D"
  claiming otherwise was proved wrong and rejected).
- $\Psi \propto \zeta$ factorization under equal scaling, with
  $\Omega = (\alpha/s + q/c)/(\beta/n^2 + \chi^2/c^2 + \alpha\kappa^2/s^2) > 0$.
- A proposal to drop $q^O > q^R$ and re-found switching on $h_O > h_R$ as a primitive
  was proved incoherent; the trichotomy was adopted instead.

## Session log (compressed)

Full audit of the note → adjudicated vs an external audit → switching result derived
from scratch, reviewed twice, trichotomy added → readability rewrite (de-wedging:
primitives/observables before derived objects, remarks folded into prose) → slides v2
distilled to a conventional deck → structural simplification to Coven-style single
floorspace clearing + Menzio-style environment-first §1 (commit 9f357f7) → Coven PDF
actually read (toy model, setup, Definition 1, government; pp. 1–4, 13, 17–25) →
a5508ce external patch (described above, unaudited) → created the
`writing-econ-papers` skill (commit 8d5d90a).

## Files

- `latex/intergenerational_housing_fertility_part1.tex` (+pdf) — THE draft; clean at HEAD.
- `latex/intergenerational_housing_fertility_part1_slides_v2.tex` — my deck; has
  UNCOMMITTED user edits (mid-sync to single market; last I checked, leftovers: $q^m$
  in the branch-problem display, $P^O$ in the $\zeta$ collateral term). Don't revert.
- `latex/intergenerational_housing_fertility_model_summary_2page.tex` — user's two-page
  summary; UNCOMMITTED edits in flight; may still describe the old two-stock model.
- `latex/intergenerational_housing_fertility_part1_slides.tex` — ChatGPT's old deck,
  comparison only. `latex/part1_policy_sufficient_stat_figures_20260610.tex` — superseded.
- `docs/style/econ_writing_style_guide.md` — writing style guide (mirror of the
  `writing-econ-papers` personal skill); CLAUDE.md/AGENTS.md point to it. Invoke the
  skill for ALL paper-facing text; attach the docs mirror to ChatGPT/Codex prompts.
- `docs/prompts/part1_switching_*_review_2026-06-11.md` — the two switching review prompts.
- Coven, Golder, Gupta, Ndiaye (2025) PDF: `~/Desktop/Property_Tax.pdf`. Verified facts:
  single elastic supply $H^{supply} = c_i P^{\rho_i}$, $R = (\tau + r + \delta - \gamma)P$,
  ownership premium $\Xi^O$, min owned size $\underline H$ whose stated role is making
  the DP bind, and (p. 4) they explicitly ABSTRACT from tenure-based assessment lock-in —
  our $L^\tau$ is exactly that channel (positioning sentence not yet added).

## Open items

1. Audit a5508ce's new local no-reallocation proposition (and adjudicate the incoming
   ChatGPT audit against the verified-math list above).
2. Finish slides v2 single-market sync once the user's edits settle; then the two-page
   summary.
3. Switching proposition: the bundle ordering ($c_B < c_A$, $s_B > s_A$) is stated under
   equal scaling but the proof is ES-free — restore the unconditional statement.
4. Possible positioning sentence: $L^\tau$ = the assessment channel Coven et al. abstract from.
5. Quant model: measure benchmark per-unit owner user cost vs rent to select the
   operative branch of the switching trichotomy under wedges.

## Constraints and gotchas

- $\phi = 0.80$ non-negotiable ($\phi$ = financed share; DP uses $1-\phi$). No structural
  model changes unless explicitly asked. Tenure choice stays deterministic argmax.
  $(n,s)$ notation for parity/child-state.
- Never describe a paper's model from memory — open the PDF. poppler/pdftotext and the
  Read-PDF path are BROKEN on this machine (missing libintl dylib); render with
  `gs -dNOPAUSE -dBATCH -dQUIET -sDEVICE=png16m -r70 -sOutputFile=/tmp/p%d.png file.pdf`
  and Read the PNGs.
- Compile LaTeX twice (refs); after structural edits grep for orphaned `\ref`/`\eqref`
  and leftover old-structure symbols. `output/` is gitignored. Transient pdflatex
  `fflush ... Operation timed out` errors: just recompile.
- The repo has many uncommitted user files (slides v2, 2-page summary, model code,
  data scripts) — never stage, clean, or revert them wholesale.
