# Supplemental slide inputs

**2007 stationary reference — block0506, September 28 verified export.**

Ready outputs: [housing-supply figure PDF](rendered_output_v2/supplemental_housing_supply.pdf),
[PNG preview](rendered_output_v2/supplemental_housing_supply.png),
[borrowing table TeX](rendered_output_v2/borrowing_comparison.tex), and
[full-precision table CSV](rendered_output_v2/borrowing_comparison.csv).
Torch render job **18801894** completed in six seconds with zero model solves.
The lead visually inspected the PNG and checked the table against the source
receipts. Labels distinguish stationary household population, physical housing,
prices and mapped rents; the figure explicitly retains lifetime repayment.
The [manifest](rendered_output_v2/manifest.json) pins the actual inputs and
frozen render source. `output_v2.sha256` authenticates transferred output files.
The original render (18801617) remains in `rendered_output_v1/`.
No new deck, manuscript edit or assembled PDF report has been made. The local
price-elasticity figure awaits verified numerical results.

`build_slide_inputs.py` makes compact booktabs tables from passed, source-pinned receipts and result CSVs. It adds the supplemental housing-supply figure when its verified result pair is supplied, and adds the prescribed-price figure and local-elasticity table when their passed results are supplied. It uses no model code. Run rendering on Torch after lead authorization; this Mac-side preparation does not import the builder, calculate results, or render figures.

The borrowing table requires four receipts, their four existing fit CSVs, and an output path. Supply both `--supply-comparison` and `--supply-verification` to add the housing-supply figure. To add the prescribed-price figure and local-elasticity table, also pass all three of `--price-comparison`, `--elasticities`, and `--elasticity-completed`; each is checked against the passed completion receipt's SHA-256 and frozen reference label before output is produced. The comparison CSV must have the runner schema (`regime`, `scope`, `price_factor`, `outcome`, `value`); result paths are explicit CLI arguments, so later versioned result folders can be used without changing this builder.

Example after both result packets are available (the bounded Torch job script in this folder supplies these paths):

```sh
python build_slide_inputs.py \
  --reference-receipt /path/to/control/receipt.json \
  --grid-receipt /path/to/grid_control/receipt.json \
  --credit-receipt /path/to/credit/receipt.json \
  --ge-receipt /path/to/selected_repeat/receipt.json \
  --reference-fit /path/to/control/target_fit.csv \
  --grid-fit /path/to/grid_control/target_fit.csv \
  --credit-fit /path/to/credit/target_fit.csv \
  --ge-fit /path/to/selected_repeat/target_fit.csv \
  --price-comparison ../elasticity_v1/results_v2/comparison.csv \
  --elasticities ../elasticity_v1/results_v2/elasticities.csv \
  --elasticity-completed ../elasticity_v1/results_v2/completed.json \
  --supply-comparison ../supply_v1/results_v1/comparison.csv \
  --supply-verification ../supply_v1/results_v1/verification.json \
  --output generated_v2
```

The existing fit paths supply `nchs_mean_age` because the borrowing receipts omit mean age at first birth. Output directories are immutable: an existing destination causes an error. The manifest lists the exact reference label and SHA-256 hashes for every actual input and this builder.

For the already-passed borrowing and supply results named in `run_existing_render.sh`, stage that script beside the builder on Torch and submit it with `sbatch -o /scratch/td2248/projects/fixed_reference_economic_figures_20260929/results/render_%j.out /scratch/td2248/projects/fixed_reference_economic_figures_20260929/source_v1/run_existing_render.sh`. It runs one CPU with 4 GiB for at most five minutes and performs zero model solves. The script does not include the optional price-elasticity inputs.

The figures use a restrained, flat academic style and carry the frozen reference label. The prescribed-price figure uses changes relative to each regime's own baseline price; its impact panel shows immediate births and its cohort panel shows completed fertility. The supply figure uses endpoint values relative to the frozen baseline and labels population in household units. Neither figure claims a transition result. TeX tables round to three decimals; CSV tables retain full precision. Blank or em-dash population entries indicate partial-equilibrium columns.

For experiment definitions and status, see [the fixed-reference economics README](../README.md). The underlying 14-moment, 31-parameter fit remains in the source run packet and is not duplicated in these slide tables; see [the fertility-identification README](../../fertility_identification_20260928/README.md) for that index.
