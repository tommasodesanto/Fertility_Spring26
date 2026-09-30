# Credit-rule quick v1 — prepared, not submitted

This packet compares exactly two fixed-price household/cohort cases from **2007 stationary reference — block0506, September 28 verified export**: `ours` has `lambda_d=0`, zero debt caps, and the direct (J+1) saving-floor weights that implement \(b'\geq\min(b,0)\) only when death is impossible; `author` has zero taper weights and caps plus the private compiled strict tenure kernels.  The latter makes an owner-to-renter sale infeasible when raw \(b+(1-\psi)pH<0\), before the tenure maximum/logsum.  It does not alter owner stayers, the estate floor, purchases/LTV, the 160-node grid, prices, fiscal inputs, preferences, entry assets, or source identity.

Each fresh child writes `preliminary_numbers.json` immediately after its one lifecycle solve, then evaluates the inherited authenticated `g_pre` with the supplied policy without another lifecycle solve.  It subsequently attempts the maintained gates and writes the 14-row fit and 31-parameter tables.  Gate failures remain explicit and make the preliminary numbers uncertified.  The packet never clears markets, re-normalizes fertility, recalibrates, or runs GE.

Pins: driver `192ab1e6022e4c8a7328647132841fab82dddabaad2c68bc685e6c45a40d0ca8`; strict module `c64f98bf72484facc36be5eb4d38d581c1f813b419feae9ab28eee3e96117677`; original helper `96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44`.

Stage these four packet files unchanged to `/scratch/td2248/projects/fixed_reference_credit_rule_quick_20260929/source_credit_rule_quick_v1/`, verify the three hashes, then (before the hard deadline only) submit:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/credit_rule_quick_v1/launch_credit_rule_quick.sh
```

The array launches one CPU/24 GiB task per case in parallel and refuses duplicate outputs.  The supplied boundary fixture tests negative, zero, and positive sale equity; run it only in a Torch authenticated environment because it imports the frozen package.
