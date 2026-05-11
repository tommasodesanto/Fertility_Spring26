# `code/model/docs/`

Two HTML guides to the Python model code:

- **[`architecture.html`](architecture.html)** — dry reference doc. Top-to-bottom file-by-file
  walkthrough of `dt_cp_model/`. Use this when you need to navigate the code.

- **[`tour.html`](tour.html)** — Distill-style interactive tour. Math from the analytical
  notes, code excerpts from the model, and a full policy / GE-replay explorer
  driven by one inlined solver run. Use this when you want to *understand* the
  model, not just the directory layout.

## Rebuilding `tour.html`

```bash
cd code/model
make tour-rebuild         # re-runs the solver dump + reassembles the HTML
# or, finer-grained:
make tour-data            # only dump tour_data.json
make tour                 # only re-inline data + JS + code excerpts
```

The build steps:

1. `tools/dump_explorer_data.py` runs the Python solver at `setup.x0` (in `fast`
   mode) with `P.collect_ge_trace = True` and writes `tour_data.json`. This
   captures policy slices, GE iteration trace, moment targets/weights/model
   values, and a summary block. Takes ~20 s.
2. `tools/assemble_tour.py` reads `tools/tour_template.html`, extracts eight
   code excerpts from `kernels.py` and `solver.py` for the chapter-3 dual-view
   panels, inlines `tools/tour_widgets.js` and `tour_data.json`, and writes
   `tour.html`. Takes ~50 ms.

The resulting `tour.html` is self-contained (~4.3 MB, mostly the JSON) and opens
from `file://` — no server required. KaTeX, Prism, and Google Fonts come from
CDN, so an internet connection is needed for first render (after which the
browser cache covers offline opens).
