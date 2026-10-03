"""Run one explicit fixed-price experiment from editable values below.

Importing this module is inert. Execute the file only when you want the solve.
"""
from __future__ import annotations

OVERRIDES = {"beta_annual": 0.98}
PRICE = None  # None uses the authenticated soft reference price.


def main() -> None:
    import matplotlib.pyplot as plt

    from model_playground import ModelPlayground, load_saved_solution, ModelResult

    model = ModelPlayground()
    baseline_sol = load_saved_solution("soft")
    baseline = ModelResult(
        baseline_sol, model.P, model.params, float(baseline_sol.price),
        label="saved selected soft baseline",
    )
    changed = model.solve(price=PRICE, overrides=OVERRIDES)
    comparison = model.compare(baseline, changed)

    print("Fixed-price experiment; no calibration or market-price root was run.")
    print(f"Baseline price: {baseline.price:.12g}; changed price: {changed.price:.12g}")
    print("Parameter changes:")
    for key, value in OVERRIDES.items():
        print(f"  {key}: {baseline.parameters[key]:.12g} -> {value:.12g}")
    print("Overall aggregate levels and changed-minus-baseline differences:")
    for name, row in comparison["overall"].items():
        print(f"  {name}: {row['baseline']} -> {row['changed']} (change {row['difference']})")

    fields = ("mean_consumption", "mean_rooms", "ownership_rate")
    ages = [row["age"] for row in comparison["by_age"]]
    fig, axes = plt.subplots(1, len(fields), figsize=(13, 3.8), constrained_layout=True)
    for ax, field in zip(axes, fields):
        old = [row[field]["baseline"] for row in comparison["by_age"]]
        new = [row[field]["changed"] for row in comparison["by_age"]]
        ax.plot(ages, old, "o-", label="saved soft baseline")
        ax.plot(ages, new, "o-", label="fixed-price experiment")
        ax.set(title=field.replace("_", " "), xlabel="Age")
    axes[0].legend(frameon=False)
    plt.show()


if __name__ == "__main__":
    main()
