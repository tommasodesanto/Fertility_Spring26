# Part A tables — capped renter parents (pre-tenure weights)

Population: origin-tenure renters (to=0) with children at home m>=1 whose
conditional renter policy hR_pol is at the 6-room cap, weighted by saved
pre-tenure mass g_beginning_distribution. Per-cell CSVs: partA_<arm>.csv;
pooled: partA_pooled.csv. Choice-specific tenure values are NOT saved in
stage/solution_arrays.npz (kernels return VH/tcj/probs only), so no
rent-vs-best-rung value gap is reported.

Operative screen (native_purchase_income=True): renter buying rung tn is
feasible iff b >= (1-phi)*p*H - y/Rg AND b - p*H >= max(-phi*p*H - y/Rg, b_grid[0]).

z_grid: 0.1035, 0.1713, 0.2836, 0.4696, 0.7776, 1.287, 2.132, 3.529, 5.843
Rg: 1.08243
