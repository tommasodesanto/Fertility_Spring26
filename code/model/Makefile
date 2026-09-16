.PHONY: tour-data tour clean-tour sandbox sandbox-check sandbox-test

# Run the solver once at the calibrated theta and dump curated JSON
tour-data:
	PYTHONPATH=. python tools/dump_explorer_data.py fast

# Inline JSON + code excerpts + widget JS into tour.html
tour:
	PYTHONPATH=. python tools/assemble_tour.py

# Re-build everything from scratch
tour-rebuild: tour-data tour

clean-tour:
	rm -f docs/tour.html docs/tour_data.json

# Solve the 2007 stationary state for one spec: make sandbox SPEC=kappa_h_zero
sandbox:
	PYTHONPATH=. python sandbox/run_ss.py --spec $(SPEC)

# Unit tests for the sandbox's default-off mechanism switches
sandbox-test:
	PYTHONPATH=. python -m pytest sandbox/tests -q

# Regression gate: baseline spec must reproduce the retained September 13
# 2007 calibration (output/model/e5f_final_night_20260913/corrected_initial/)
sandbox-check: sandbox-test
	PYTHONPATH=. python sandbox/check_regression.py
