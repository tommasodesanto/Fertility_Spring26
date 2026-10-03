"""Plot saved household children-ever-born distributions by lifecycle age.

This is a read-only reduction of two saved ``solution_arrays.npz`` files; it
does not call the model solver. Use --phi-080, --phi-100 and --output to apply
the same measurement to another matched pair of solutions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
COUNT_LABELS = ("0 children", "1 child", "2 children", "3+ children")
COLORS = ("#456990", "#e08d3c", "#4c956c", "#9b5de5")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_executed_fields(path: Path) -> dict:
    """Read only top-level scalar fields; the snapshot contains large arrays."""
    required = {"J", "age_start", "da", "n_parity", "use_postdecision_current_distribution", "phi"}
    fields = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith('  "'):
                continue
            key = line.split('"', 2)[1]
            if key not in required:
                continue
            value = line.split(":", 1)[1].strip().rstrip(",")
            if key == "phi" and value == "[":
                items = []
                for subline in stream:
                    if subline.strip() == "],":
                        break
                    items.append(float(subline.strip().rstrip(",")))
                fields[key] = items
            else:
                fields[key] = json.loads(value)
            if required <= fields.keys():
                return fields
    raise ValueError(f"Missing executed parameter fields {required - fields.keys()} in {path}")


def extract_arm(npz_path: Path, expected_phi: float) -> tuple[dict, dict]:
    p_path = npz_path.with_name("executed_P.json")
    p = read_executed_fields(p_path)
    if int(p["n_parity"]) != 4 or not bool(p["use_postdecision_current_distribution"]):
        raise ValueError(f"Unexpected count or timing contract: {p}")
    if not np.allclose(p["phi"], expected_phi, rtol=0, atol=1e-12):
        raise ValueError(f"Executed phi does not match requested arm: {p['phi']}")

    with np.load(npz_path) as saved:
        # g: assets, tenure, location, lifecycle age, income, n ever born, m at home.
        g = saved["g"]
        beginning = saved["g_beginning_distribution"]
        if g.ndim != 7 or g.shape[3] != int(p["J"]) or g.shape[5] != 4:
            raise ValueError(f"Unexpected g shape {g.shape} and P.J={p['J']}")
        if beginning.shape != g.shape:
            raise ValueError("Beginning and current distributions have different shapes")
        if not np.all(np.isfinite(g)) or not np.all(np.isfinite(beginning)):
            raise ValueError("A saved distribution has nonfinite cells")
        min_cell = float(g.min())
        if min_cell < -1e-14:
            raise ValueError(f"Negative distribution cell {min_cell}")
        age_count_mass = np.sum(g, axis=(0, 1, 2, 4, 6))
        beginning_age_count_mass = np.sum(beginning, axis=(0, 1, 2, 4, 6))

    age_mass = age_count_mass.sum(axis=1)
    if not np.isclose(age_mass.sum(), 1.0, rtol=0, atol=1e-10):
        raise ValueError(f"Saved g has unexpected total mass {age_mass.sum()}")
    if np.any(age_mass <= 0):
        raise ValueError("At least one age cell has zero household mass")
    shares = age_count_mass / age_mass[:, None]
    cdf = np.cumsum(shares, axis=1)
    if not np.allclose(shares.sum(axis=1), 1, rtol=0, atol=2e-14):
        raise ValueError("Age-conditional count shares do not sum to one")
    if np.min(np.diff(cdf, axis=1)) < -1e-14 or np.min(cdf) < -1e-14:
        raise ValueError("CDF is not monotone/nonnegative")
    if not np.allclose(cdf[:, -1], 1, rtol=0, atol=2e-14):
        raise ValueError("CDF terminal value is not one")
    if not np.allclose(age_count_mass, beginning_age_count_mass, rtol=0, atol=2e-11):
        raise ValueError("Current choices changed age/count margins unexpectedly")

    ages = float(p["age_start"]) + np.arange(int(p["J"])) * float(p["da"])
    record = {
        "ages": ages,
        "age_count_mass": age_count_mass,
        "age_mass": age_mass,
        "shares": shares,
        "cdf": cdf,
    }
    checks = {
        "solution_arrays": str(npz_path.resolve()),
        "solution_arrays_sha256": sha256(npz_path),
        "executed_P": str(p_path.resolve()),
        "executed_P_sha256": sha256(p_path),
        "executed_fields": p,
        "g_shape": list(g.shape),
        "g_mass": float(age_mass.sum()),
        "minimum_g_cell": min_cell,
        "minimum_age_mass": float(age_mass.min()),
        "maximum_age_share_sum_error": float(np.abs(shares.sum(axis=1) - 1).max()),
        "minimum_cdf_step": float(np.diff(cdf, axis=1).min()),
        "maximum_cdf_terminal_error": float(np.abs(cdf[:, -1] - 1).max()),
        "maximum_current_beginning_age_count_mass_gap": float(np.abs(age_count_mass - beginning_age_count_mass).max()),
    }
    return record, checks


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_outputs(phi080: Path, phi100: Path, output: Path, title: str) -> None:
    arms = {}
    checks = {}
    for label, path, phi in (("phi_080", phi080, 0.8), ("phi_100", phi100, 1.0)):
        arms[label], checks[label] = extract_arm(path, phi)
    base, relaxed = arms["phi_080"], arms["phi_100"]
    if not np.array_equal(base["ages"], relaxed["ages"]):
        raise ValueError("The two arms have different lifecycle age cells")
    output.mkdir(parents=True, exist_ok=True)

    rows = []
    diff_rows = []
    for j, age in enumerate(base["ages"]):
        for label in ("phi_080", "phi_100"):
            arm = arms[label]
            row = {"arm": label, "age_start": age, "age_end_exclusive": age + float(checks[label]["executed_fields"]["da"]),
                   "age_cell_index": j, "household_mass": arm["age_mass"][j]}
            for n, suffix in enumerate(("0", "1", "2", "3plus")):
                row[f"mass_n_{suffix}"] = arm["age_count_mass"][j, n]
                row[f"share_n_{suffix}"] = arm["shares"][j, n]
            for n in range(3):
                row[f"cdf_n_le_{n}"] = arm["cdf"][j, n]
            rows.append(row)
        diff = {"age_start": age, "age_end_exclusive": age + float(checks["phi_080"]["executed_fields"]["da"]),
                "age_cell_index": j}
        for n, suffix in enumerate(("0", "1", "2", "3plus")):
            diff[f"share_n_{suffix}_difference_pp"] = 100 * (relaxed["shares"][j, n] - base["shares"][j, n])
        for n in range(3):
            diff[f"cdf_n_le_{n}_difference_pp"] = 100 * (relaxed["cdf"][j, n] - base["cdf"][j, n])
        diff_rows.append(diff)
    write_csv(output / "children_by_age.csv", list(rows[0]), rows)
    write_csv(output / "differences_by_age.csv", list(diff_rows[0]), diff_rows)

    ages = base["ages"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True, layout="constrained")
    for n, ax in enumerate(axes.flat):
        for label, style in (("phi_080", "-"), ("phi_100", "--")):
            ax.plot(ages, arms[label]["shares"][:, n], style, linewidth=2.2,
                    color=COLORS[n], label=r"$\phi=0.8$" if label == "phi_080" else r"$\phi=1.0$")
        ax.set_title(f"{COUNT_LABELS[n]} ever born")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.2)
        if n >= 2:
            ax.set_xlabel("Age cell start (years)")
        if n % 2 == 0:
            ax.set_ylabel("Share of households at this age")
    axes[0, 0].legend(frameon=False)
    fig.suptitle(title + "\nDistribution includes childless households; four-year age cells", fontsize=13)
    fig.savefig(output / "child_count_shares_by_age.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True, layout="constrained",
                             gridspec_kw={"width_ratios": (1.2, 1)})
    for n in range(3):
        left, right = axes[n]
        for label, style in (("phi_080", "-"), ("phi_100", "--")):
            left.plot(ages, arms[label]["cdf"][:, n], style, linewidth=2.2, color=COLORS[n],
                      label=r"$\phi=0.8$" if label == "phi_080" else r"$\phi=1.0$")
        right.plot(ages, 100 * (relaxed["cdf"][:, n] - base["cdf"][:, n]),
                   color=COLORS[n], linewidth=2.2)
        right.axhline(0, color="#666666", linewidth=0.8)
        left.set_ylim(0, 1)
        left.set_ylabel(rf"$P(N\leq {n}\mid\mathrm{{age}})$")
        right.set_ylabel("Difference (pp)")
        for ax in (left, right):
            ax.grid(alpha=0.2)
        if n == 0:
            left.set_title("CDF of children ever born")
            right.set_title(r"$\phi=1.0$ minus $\phi=0.8$")
            left.legend(frameon=False)
        if n == 2:
            left.set_xlabel("Age cell start (years)")
            right.set_xlabel("Age cell start (years)")
    fig.suptitle(title + "\nCDF includes all households, including childless", fontsize=13)
    fig.savefig(output / "child_count_cdf_by_age.png", dpi=180)
    plt.close(fig)

    checks["definition"] = {
        "age": "age of the reproductive household member; cell starts at age_start + j*da",
        "distribution": "saved post-birth current-period household cross-section g",
        "g_axes": ["asset", "tenure", "location", "age", "income", "children_ever_born_n", "children_at_home_m"],
        "conditioning": "within each arm and each age cell, divide mass in n by all household mass at that age",
        "count_states": ["0", "1", "2", "3 or more"],
        "cdf": "P(N <= k | age) for k=0,1,2",
        "difference": "phi=1.0 minus phi=0.8 in percentage points",
        "price_closure": "input saved solutions; consult their source receipt for fixed-price status",
    }
    (output / "checks_and_provenance.json").write_text(json.dumps(checks, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phi-080", type=Path, default=PARENT / "phi_080" / "solution_arrays.npz")
    parser.add_argument("--phi-100", type=Path, default=PARENT / "phi_100" / "solution_arrays.npz")
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--title", default="Fixed-price credit relaxation")
    args = parser.parse_args()
    make_outputs(args.phi_080, args.phi_100, args.output, args.title)


if __name__ == "__main__":
    main()
