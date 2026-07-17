# Task: literature survey — how structural models reconcile subsistence floors with realistic income risk

You are surveying literature for a quantitative housing–fertility lifecycle
model. Deliver a decision memo, not an essay. Every claim about a paper must
carry an exact quote or equation with page/table number, verified against the
paper itself (not an abstract, not a summary site). If you cannot verify a
number, write NOT VERIFIED. Do not pad; we act on this tomorrow.

## Our situation (verified facts, do not re-derive)

- Discrete-time OLG lifecycle model, 4-year periods, entry at 18. Households
  face Stone–Geary subsistence requirements every period: consumption floor
  c̄₀ (calibrated, currently ≈ 28–31% of mean annual household income) plus a
  minimum rental dwelling h̄₀ (rent ≈ 2–4% of mean income). Renters cannot
  borrow (unsecured line is zero); owner debt is collateralized; no default.
- Current income process: 5-state Rouwenhorst AR(1), annual rho 0.9602,
  innovation sd 0.0645 — a provisional construction (reverse-engineered from
  an old grid), NOT defensible. Stationary log variance 0.053, roughly 10–25%
  of literature values.
- When we substitute a literature process (rho = 0.90, sd = 0.12–0.20 annual,
  e.g. Sommer–Sullivan), the model becomes infeasible: households persistently
  at the bottom income state (income ≈ 21% of mean at young ages) cannot
  afford c̄₀ + rent(h̄₀) ≈ 32% of mean income from income plus any borrowing
  or asset sales. Verified by direct probes: ~0.16% of simulated mass at age
  22 has an empty budget set even when everyone enters with zero debt. The
  binding households are young childless RENTERS at the lowest persistent
  income state — not leveraged owners, not debtors.
- Known standard devices, already verified by us from the papers:
  (a) De Nardi–French–Jones (2010 JPE, eq. 10, p. 44): means-tested transfers
      b = max{0, c̲ + m − resources}, asset-tested (if transfers > 0 then
      c = c̲ and a' = 0), floor ESTIMATED at $2,700/yr (1998$), following
      Hubbard–Skinner–Zeldes (1994, 1995). Unfunded (partial equilibrium).
  (b) De Nardi (2004 REStud, Table A.2): no explicit floor; the income
      process itself is earnings PLUS social-insurance transfers, and the
      model has no subsistence floors, so empty budget sets cannot arise.
- Our objection to (a) at our parameters: pinning the floor to OUR bundle
  (~32% of mean income) would make transfers cover a large share of young
  households — a mass welfare program, not a safety net — because our
  calibrated c̄₀ is far above real-world program generosity (~15–20%).
  Setting a realistic floor instead imposes an implicit upper bound on c̄₀
  whose consequences for the fertility/housing calibration are unknown.

## Questions (answer each; rank candidate solutions at the end)

1. **Floor levels in practice.** For every structural savings/housing model
   you can verify that uses an HSZ-style consumption floor (working-age or
   retirement): what floor level, in dollars and as a share of mean/median
   household income? At minimum: Hubbard–Skinner–Zeldes (1994, 1995);
   French (2005); French–Jones (2011); Kaplan–Violante (2014, if floor);
   Ameriks–Briggs–Caplin–Shapiro–Tonetti (2020); Scholz–Seshadri–Khitatrakun
   (2006). Exact quotes.
2. **Stone–Geary + income risk.** Find structural papers that combine
   Stone–Geary/subsistence preferences with persistent income risk. How do
   they avoid empty budget sets — floor transfers, truncated income support,
   subsistence levels kept below minimum income, borrowing against future
   income, or something else? Fertility papers with child-cost floors
   (Barro–Becker quantitative implementations, Daruich, Doepke–Kindermann,
   Sommer 2016 "Fertility choice in a life cycle model with idiosyncratic
   uninsurable earnings risk") are especially relevant — Sommer 2016 is
   probably the closest cousin: fertility + earnings risk. What does she do
   at the bottom?
3. **Housing minimums + income risk.** How do housing models with a minimum
   dwelling size or discrete rental floor (Sommer–Sullivan 2018; Kaplan–
   Mitman–Violante 2020; Boar–Gorea–Midrigan; Favilukis–Ludvigson–Van
   Nieuwerburgh; Greaney–Parkhomenko–Van Nieuwerburgh if applicable) keep the
   poorest households feasible? Do any allow doubling up / moving in with
   family / homelessness states, or is the rental floor simply set low enough
   relative to the income support?
4. **Income concept.** Which of the standard income-process estimates are for
   PRE-transfer earnings vs POST-tax-and-transfer income? For each value
   commonly imported (Sommer–Sullivan's 0.90/0.20 in particular — trace what
   they estimated it on, PSID concept and sample), state the concept. Is
   there a citable post-transfer household income process (rho, sigma) we
   could import directly (e.g., from Blundell–Pistaferri–Preston, Krueger–
   Perri, Heathcote–Storesletten–Violante calibrations)?
5. **Funding.** Among papers with transfer floors, which fund them within
   the model (budget balance, payroll/income tax) and which leave them
   unfunded? One line each.
6. **The c̄₀ tension.** Any precedent for subsistence consumption calibrated
   near 30% of mean income? What magnitudes do calibrated Stone–Geary papers
   use? If ours is an outlier, say so plainly.

## Deliverable

A ranked menu of implementable solutions for our model (each with: exact
mechanism, citations with verified quotes, expected side effects on wealth
distribution / fertility margins / ownership, one-line implementation cost),
plus your single recommendation and the two or three sentences we would write
in the paper. Flag every claim you could not verify.
