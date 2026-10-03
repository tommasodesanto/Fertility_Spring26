# Python model playground

The saved soft case is available as arrays, and the selected soft reference can
be rebuilt as native Python inputs. Loading either does not solve the model.
The native interface uses the selected soft parameter vector, the checked
target and weight contract, nonnegative-mean entry, zero unsecured renter credit,
soft purchase financing, and the derived \(H_0\) and exact reference price.
The direct interface was checked on October 2: one 6.14-second native solve
reproduced all 11 checked policy/distribution arrays exactly. This is a replay
check, not a grid-convergence certificate.

Start the read-only browser explorer by double-clicking
`code/model/tools/start_model_explorer.command`. It serves saved cases at
<http://127.0.0.1:8765>; it does not run solves or overwrite output. Its command
window stays open while the server runs. If that URL already responds, the
launcher reports it and leaves the existing process alone.

For direct Python control, double-click
`code/model/tools/start_model_playground.command`. It opens the interactive
session below with all six numerical thread limits set to one; loading the
saved arrays at startup does not solve the model.

From the repository root, start an interactive Python session:

```bash
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  code/model/.venv/bin/python -i code/model/tools/model_playground.py
```

`sol = load_saved_solution()` loads the hash-checked original-timing arrays into
a `SimpleNamespace`, including `sol.b_grid`, `sol.V`, `sol.c_pol`, policy
arrays, and `sol.price`. For example, plot consumption by wealth at age 30 for
income state five, childless, conditional on renting (not averaged over tenure
choices). Axes are wealth, tenure, location, age, income, children ever born,
and children at home. Age index 3 is age 30; income index 4 is state five.

```python
import matplotlib.pyplot as plt
plt.plot(sol.b_grid, sol.c_pol[:, 0, 0, 3, 4, 0, 0], ".-", markersize=2)
plt.xlabel("Financial wealth / mean annual gross earnings")
plt.ylabel("Four-year consumption / mean annual gross earnings")
plt.show()
```

Build the corresponding native inputs without solving:

```python
P, b_grid, solver = load_reference_model()
P.reference_price
P.H0
```

To change a preference parameter and inspect the resulting fixed-price
partial-equilibrium solution, edit the parameter, rebuild shared objects, and
explicitly request one solve:

```python
P.hbar_first_child_jump -= 0.1
SD = solver.precompute_shared(P, b_grid)
changed = solver.solve_markov_income_at_prices(
    [P.reference_price], P, b_grid, SD=SD, fast_stats=False
)
from small_credit_lab.engine import diagnostics
from pathlib import Path
diagnostics.write_diagnostics(changed, P, Path("output/model/my_fixed_price_diagnostic"))
```

For another price with the same parameters, call
`solve_at_price(P, b_grid, solver, price)`. A changed price does not solve the
renewal-price root or re-derive \(H_0\); this is a fixed-price exercise. A changed
parameter likewise does not recalibrate or clear markets. Use a copy of `P` if
you want to preserve the original initialized inputs.

The saved grid has 120 nodes, from -12 to 3000. The upper endpoint has zero
mass. A nonuniform grid is not a household distribution: inspect
`sol.g_beginning_distribution` for mass after fertility and before tenure,
and `sol.g` for the realized cross-section. `sol.c_pol` is a conditional policy;
owner stayers have their own `sol.c_pol_stay` and `sol.bp_pol_stay` arrays.
The browser defaults to raw conditional policies on saved grid nodes. Its
separate average-over-tenure mode applies the transaction maps and choice
probabilities. The range selector only changes the displayed interval.

## Model source map

The active playground imports the hash-checked `small_credit_lab` solver. Its
model stages are byte-identical to the canonical extracted stages documented
in `code/model/refactor_lab/README.md`:

| File | Role |
|---|---|
| `code/model/tools/model_playground.py` | Python entry point and saved/reference loading |
| `code/model/refactor_lab/engine/household.py` | Household choices and Bellman solution |
| `code/model/refactor_lab/engine/distribution.py` | Population distribution and moments |
| `code/model/refactor_lab/engine/equilibrium.py` | Fixed-price solve and market equilibrium |
| `code/model/refactor_lab/engine/parameters.py` | Primitive helpers and defaults |

The initialized `P` comes from the selected soft input bundle and checked
reference contract. Module defaults are not a substitute for those loaded
parameters.
