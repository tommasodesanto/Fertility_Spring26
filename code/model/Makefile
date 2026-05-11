.PHONY: tour-data tour clean-tour

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
