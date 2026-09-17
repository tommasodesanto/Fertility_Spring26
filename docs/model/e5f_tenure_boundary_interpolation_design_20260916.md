# Design note: interpolating the tenure boundary inside a wealth-grid cell

September 16, 2026. Proposal only; nothing implemented. Written for the author
decision recorded as active point 11 in `docs/model/ACTIVE_DECISION_LEDGER.md`
(status WAIT until the utility-function re-evaluation is settled).

## Problem

With the tenure taste scale at its frozen value \(\kappa_H=0.005\), the
own-versus-rent choice is numerically an argmax. Every household on one
wealth node \(b_k\) (same age, children, income state) makes the same choice,
so when the price path moves the indifference point across \(b_k\) the whole
node's mass flips tenure at once. Along the 104-date announced transition
this produces ownership jumps of 1.2--1.7 percentage points between adjacent
dates at two or three dates, housing-demand jumps of about 0.2 percent, and
rebate-budget jumps of about 0.1 percent. No smooth root exists across such a
date at the retained gates, and the August 27--29 conditional-share
(complementarity) route showed that treating each flipping atom as an
equilibrium object is exact but combinatorial. Raising \(\kappa_H\) removes
the flips but destroys the family ownership gap (0.161 at 0.005, 0.121 at
0.01, 0.062 at 0.02, 0.029 at 0.05; target 0.168), so it is not available.

## Proposal

Keep every household's choice deterministic and keep the Bellman value as the
max. Change only how a node's mass is allocated in the forward step near a
tenure boundary:

1. In the Bellman loop, store alongside the argmax the value gap between the
   chosen tenure and the best alternative, \(\Delta_k = V^{\text{chosen}}_k -
   V^{\text{alt}}_k \ge 0\), signed by which tenure is chosen (one float per
   node and state, same shape as `tenure_choice`).
2. In the forward step, whenever neighbouring nodes \(b_k, b_{k+1}\) choose
   different tenures, locate the crossing by linear interpolation of the
   signed gap, \(b^* = b_k + \frac{\Delta_k}{\Delta_k + \Delta_{k+1}}(b_{k+1}-b_k)\),
   interpret node \(k\)'s mass as spread over its cell
   \([b_k - h/2,\, b_k + h/2]\), and send the fraction of that cell lying
   beyond \(b^*\) to the neighbour's tenure. Symmetrically for node \(k+1\).
   Away from crossings nothing changes.
3. For the switched fraction, use the neighbour node's stored policy under
   its own tenure (consumption, saving, housing product). This is an
   \(O(h)\) approximation to that fraction's true conditional policy; storing
   conditional policies for both tenures at crossing nodes would remove it at
   the cost of memory, and can be a second step.

The mapping from prices to housing demand and tax revenue then becomes
continuous: as the price moves, \(b^*\) slides through the cell and the mass
moves in proportion, instead of jumping when \(b^*\) crosses a node. No taste
noise is added, so the family ownership gap is unaffected up to the cell
approximation, which shrinks with the grid.

## Where it lives

Live strand `code/model/intergen_eqscale_seq_optimized/`: the deterministic
argmax in `kernels.py` (`tenure_choice_kernel`, around line 596, and the
in-loop `np.argmax(Vopt, axis=3)` in `solver.py`), and the forward scatter
that reads `tenure_choice[b_loc, ...]` in `kernels.py` (lines 273--478). The
change is confined to (i) one extra stored array from the Bellman loop and
(ii) the branch of the forward scatter that assigns a node's mass by tenure.
It would be built in a copied kernel first, never in production.

## Verification before any use

1. Nesting: with the split disabled the copied kernel must reproduce the
   current solution byte for byte (policies, distribution, moments).
2. Grid convergence: the split kernel and the current kernel must agree on
   stationary moments as the wealth grid is refined; the gap between them
   must shrink with \(h\).
3. Continuity: the two-date and ten-date derivative smokes (seven native
   mappings) must show the ownership response to a price perturbation is
   smooth; the 104-date fixed-path probe must show the date 44 and 70--71
   jumps become ramps.
4. Economics: the stationary family ownership gap, mean rooms and the
   ownership rate must move by less than the grid-refinement gap in step 2.
5. Root: the ten-period root certifies at the retained gates, then the
   104-date root with the measured-Jacobian start.

## Effort

Prototype in a copied kernel with unit tests, about one day; the derivative
smokes and ten-period root, a few hours on Torch; the 104-date root, one
overnight. It should follow, not precede, the utility-function decision,
because a larger owner premium for parents changes the stakes at the
boundary and therefore the size of the cells that matter.
