#!/usr/bin/env python3
"""Build the compact ACS--PSID first-birth comparison figure.

The script consumes only the small estimator receipts and the saved PSID
reference CSV.  It never opens the national panel or an RDS checkpoint.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[4]
ACS_DIR = ROOT / "output/national_acs_comparison/national_continuation_20260921b"
PSID_FILE = ROOT / "output/psid_fullsample_staging_20260921/psid_reference_for_acs.csv"
OUT = ROOT / "output/national_acs_comparison/national_acs_primary_vs_psid.png"


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    acs = rows(ACS_DIR / "national_contrasts.csv")
    psid = rows(PSID_FILE)
    acs_map = {(r["outcome"], r["specification"]): r for r in acs}
    psid_map = {r["arm"]: r for r in psid}

    # Ownership is plotted in percentage points; rooms remain in room units.
    panels = [
        {
            "title": "Rooms (rooms)",
            "outcome": "rooms9",
            "acs_full": ("rooms9", "full"),
            "psid": ["first_birth_rooms_corrected"],
            "scale": 1.0,
            "psid_labels": ["PSID rooms"],
        },
        {
            "title": "Ownership (percentage points)",
            "outcome": "ownership_lw",
            "acs_full": ("ownership_lw", "full"),
            "psid": [
                "first_birth_ownership_original",
                "first_birth_ownership_corrected_hh_year",
            ],
            "scale": 100.0,
            "psid_labels": ["PSID original", "PSID aligned sensitivity"],
        },
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.7), constrained_layout=True)
    for ax, panel in zip(axes, panels):
        outcome = panel["outcome"]
        scale = panel["scale"]
        # Show all ACS specifications so the full specification is visibly a
        # prespecified primary result rather than a selected point.
        specs = [("event_only", "event only"), ("age_only", "+ age FE"),
                 ("state_year", "+ state + year FE"), ("full", "full (primary)")]
        y_labels = []
        y_positions = []
        for j, (spec, label) in enumerate(specs):
            r = acs_map[(outcome, spec)]
            x = float(r["estimate"]) * scale
            lo = float(r["conf.low"]) * scale
            hi = float(r["conf.high"]) * scale
            is_full = spec == "full"
            ax.errorbar(
                x,
                3.0 - j * 0.42,
                xerr=[[x - lo], [hi - x]],
                fmt="o" if is_full else "o",
                color="#1b5e9e" if is_full else "#7b8794",
                markerfacecolor="#1b5e9e" if is_full else "white",
                markeredgecolor="#1b5e9e" if is_full else "#7b8794",
                markersize=6 if is_full else 5,
                capsize=3,
                lw=1.5 if is_full else 1.1,
                zorder=3,
            )
            y_positions.append(3.0 - j * 0.42)
            y_labels.append(f"ACS {label}")

        psid_y = [0.95, 0.53] if len(panel["psid"]) == 2 else [0.74]
        for y, arm, label in zip(psid_y, panel["psid"], panel["psid_labels"]):
            r = psid_map[arm]
            x = float(r["estimate"]) * scale
            lo = float(r["ci_lo"]) * scale
            hi = float(r["ci_hi"]) * scale
            ax.errorbar(x, y, xerr=[[x - lo], [hi - x]], fmt="D",
                        color="#b24a2b", markerfacecolor="#b24a2b",
                        markersize=5.5, capsize=3, lw=1.4, zorder=4)
            y_positions.append(y)
            y_labels.append(label)

        ax.axvline(0, color="#9aa5ad", lw=0.8, zorder=0)
        ax.set_title(panel["title"], fontsize=11, pad=10)
        ax.set_yticks(y_positions, y_labels)
        ax.tick_params(axis="y", labelsize=8.5, length=0, pad=5)
        ax.set_xlabel("+3 minus −1 estimate with 95% CI", fontsize=9)
        ax.grid(axis="x", color="#e5e7eb", lw=0.7)
        ax.set_ylim(0.25, 3.35)
        ax.tick_params(axis="x", labelsize=8.5)

    fig.suptitle("First-birth housing response: ACS matched pseudo-panel and PSID reference",
                 fontsize=12.5, y=1.02)
    fig.text(0.5, -0.035,
             "ACS full: state + age + year FE, source-household clustered SE; "
             "PSID estimates use the saved reference definitions.",
             ha="center", fontsize=8.2, color="#4b5563")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
