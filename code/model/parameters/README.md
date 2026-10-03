# Editable stationary-model inputs

`best_params.py` holds the adopted October 3 post-interest chain-13 working
continuation inputs under the retained old wealth target. It is a verified
selected point, not a global optimum or a certified paper baseline.
`toy_params.py` is an independent personal copy with annual beta lower by 0.001.
The other inputs match. Edit the toy file or copy a complete file for another
experiment. A stationary solve uses the supplied parameters; it does not search
for a better parameter vector. Calibration is a separate search followed by a
fresh native selected-point check.

The runner selects the file with `PARAMETER_FILE = "best_params.py"`; change that
line to `"toy_params.py"` for the personal example. Relative file names always
resolve in this directory, irrespective of the shell's working directory.
Absolute paths can select a run-local calibration export. The files contain
ordinary dictionaries, comments, literals and arithmetic (including list
repetition). Executable imports, calls or control flow are rejected. There is no
module cache; the selected file is reread and copied each time.

- `PARAMETERS` contains all ten calibration coordinates. Annual beta is converted
  to the four-year model period; costs retain their native units.
- `EXTERNAL_INPUTS` holds explicit fixed primitives, including the financed share
  `phi`, renter debt capacity `unsecured_credit_limit`, period rates, gross
  earnings, payroll tax and the physical housing supply coefficient `H0`.
  `phi = [0.8] * 4` gives a 20% down-payment threshold under the retained soft purchase rule;
  financed shares must be uniform. It is not a strict liquid-cash requirement.
  Disposable earnings and balanced pension income are derived from earnings and
  payroll primitives. Entry-law and structural-grid changes require a new bundle.
- `NATIVE_OVERRIDES` is for supported advanced native fields. Unknown fields,
  incompatible shapes, contradictory controls and unsupported derived or
  structural edits fail explicitly.
- `PRICE_GUESS`, `CLOSURE` and `BUDGET_SECONDS` control one stationary solve.
  `fixed_h0` holds the supply coefficient and reports the implied population
  scale; `population_one` fixes normalized population and derives `H0`.
- `PROVENANCE` is descriptive. It does not authorize adopting an experiment.

Only the exact canonical `best_params.py` path writes production cases under
`output/model/local_solution`. Other files write under
`output/model/experiments/<file-stem>`. Run-local exports named
`best_params.py` receive a parent-label and source-path digest suffix so distinct
calibration exports remain separate. An output source marker rejects another
file with the same stem from taking over an existing experiment directory. A
parameter file cannot set its own output destination.

Future canonical post-interest calibration runs export their own `best_params.py`
only after the final fresh native acceptance and exact repeat, with pinned target
and weight fingerprints. They include effective fixed primitives, verified price
and derived `H0`, and must reproduce every input field and the grid when reloaded.
The fresh native receipt also pins the exact candidate-bound caller inputs and
grid before solver mutations. Export requires that identity before substituting
the verified derived `H0`; a changed caller or an older receipt without the
fingerprint is rejected. Roundtrip equality alone is not native authentication.
Unsupported mappings fail instead of restoring defaults. Provisional search
checkpoints, historical original-timing runs and different target contracts do
not export or overwrite this canonical default. Global promotion remains an
explicit author decision.
