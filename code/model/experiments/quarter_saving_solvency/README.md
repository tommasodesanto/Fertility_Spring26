# Quarter-saving buyer estate floor

This source copies the quarter-saving v1 engines and changes only buyer saving floors in `engine/household.py` and `engine/kernels.py`. It enforces the existing no-negative-net-estate rule inside the Bellman saving choice when death is possible: \(b'\geq-(1-\psi)Q\), where \(\psi\) is the 6% selling cost. The quarter-saving floor remains active, and owner-stayer saving rules and the forward transaction map are unchanged. The one-case \(\phi=1\) retry driver and provenance pins are in `output/model/fixed_reference_economics_20260928/quarter_saving_solvency_v2/`.

At \(\phi=0.8\), the ordinary buyer collateral floor \(-0.8Q\) is already stricter than the estate floor \(-0.94Q\), so this addition cannot change the completed v1 80% case. At \(\phi=1\), the estate floor prevents a zero-down-payment buyer from choosing debt that exceeds the home's net liquidation value at a possible death.
