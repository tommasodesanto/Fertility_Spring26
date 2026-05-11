# Empirical Roundup Implementation

This turn adds actual implementation artifacts rather than just a prompt.

## Files

- [compile_empirical_roundup_v1.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/compile_empirical_roundup_v1.py)
- [empirical_roundup_first_birth_by_wealth_v1.do](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/empirical_roundup_first_birth_by_wealth_v1.do)
- [empirical_roundup_income_elasticity_v1.do](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/empirical_roundup_income_elasticity_v1.do)

## What Each File Does

### `compile_empirical_roundup_v1.py`

Reads the existing ACS and PSID output files already on disk and writes:

- `output/empirical_roundup_v1/empirical_roundup_v1.csv`
- `output/empirical_roundup_v1/empirical_roundup_v1.md`

This is the consolidated audit artifact.

### `empirical_roundup_first_birth_by_wealth_v1.do`

Implements the most feasible missing extension with the current PSID data:

- first-birth event studies by pre-birth liquid-wealth quintile
- outcomes:
  - `rooms`
  - `own`
  - `moved_for_size`

Outputs are written to:

- `output/empirical_roundup_first_birth_by_wealth_v1/`

### `empirical_roundup_income_elasticity_v1.do`

Implements the feasible parent-vs-childless income-elasticity exercise with the
current PSID variables:

- household FE regressions
- outcomes:
  - `ln_rooms`
  - `own`
- regressor:
  - `ln_income` interacted with `any_kid`

Outputs are written to:

- `output/empirical_roundup_income_elasticity_v1/`

## Run Commands

Roundup compiler:

```bash
python3 /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/compile_empirical_roundup_v1.py
```

Stata scripts:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/empirical_roundup_first_birth_by_wealth_v1.do
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/empirical/roundup/empirical_roundup_income_elasticity_v1.do
```

## Execution Status In This Environment

- The Python roundup compiler can be run here.
- Stata is not installed in this environment, so the `.do` files were implemented
  but not executed here.
