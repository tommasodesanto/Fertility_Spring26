# Direct-period feasibility diagnosis

The first `new_balanced` stationary attempt stopped at `forward_age_30` with
dead mass (3.16347282998\times 10^{-12}), just above the unchanged
(10^{-12}) gate.  The same income-only control stopped at
(3.1634703\times 10^{-12}).  The difference is (2.5\times10^{-18}), so
the current-income purchase adapter is not the primary cause.

The failure census is informative. Every displayed state is a childless renter
at age 30, with (b=-3.046511627906977), no transfer, and one of the cold
persistent income states (z\in[0.0435,0.1216]).  Current income is positive.
Several rows have positive reported current cash slack (for example 0.0086,
0.0316, 0.0945, and 0.1252), while others are negative.  This is not a
contradiction: the census computes current-flow slack in
`_dead_mass_census_at_age` (frozen solver lines 3873–3914), whereas the dead
flag is `V <= -1e9` (lines 3865–3868).  A state can afford current flow and
still have no admissible continuation for the full next-period income support.

The backward timing confirms that interpretation. At each current income
state, continuation value averages over every `Pi_z[zz, znext] > 0` at the
next age (solver lines 2499–2512).  Forward propagation then carries the saved
wealth unchanged across the income transition (lines 5401–5448).  Thus a
negative renter asset can be selected at one age/income state and subsequently
meet a sufficiently cold persistent draw at the same wealth.  The debt rule
allows unsecured debt to roll according to
\[
\underline b_{j+1}=\min\{s_{j+1}\min(b_j,0),-D_{j+1}\},
\]
(`parameters.py:622` and `solver.py:116–126`); before the taper this can leave
negative debt at a level whose low-income continuation is infeasible.

There is also a numerical risk, but the receipt does not establish it as the
cause. The forward saving scatter clips `bp_pol` to `[b_grid[0],b_grid[-1]]`
before interpolation (solver lines 5391–5399). The purchase adapter’s ordinary
transaction map likewise clips exact transaction wealth in
`e5f_earnings_wealth_contract.py`. If the failing node is a clipped image of
an out-of-grid branch, it is interpolation-induced; if it is an on-grid
`bp_pol` node, it is a genuine natural-debt/support failure. The value
(-3.0465116279) is an interior grid-looking node, not evidence by itself of
entry projection or transaction clipping. The preserved entry receipt also
shows candidate and reference wealth marginals agree to (1.1\times10^{-16}),
so rank coupling is not changing the wealth marginal.

Falsifiable diagnostic for the next controlled run: record, for every mass
entering the failing `(age,b,z,tenure)` cell, the pre-clipping `bp_pol`, the
clipping indicator, and the predecessor `(b,z)` state. Separately recompute
the next-age value at the same (b) with (i) the native `Pi_z`, (ii) identity
income transition, and (iii) the same transition with the coldest destination
states removed. If (i) is dead but (ii)/(iii) are alive, the failure is
structural natural-debt support under the full persistent-income support. If
the clipping indicator is positive, it is a grid-support artifact and must be
fixed by support/accounting treatment rather than a gate relaxation.

An age-18 zero-liquid-asset entrant rule could be a separately labelled
Sommer/De Nardi-style sensitivity, but it changes the entry wealth contract
and cannot be used as a numerical repair. No gate relaxation or production
model change is justified by this receipt.
