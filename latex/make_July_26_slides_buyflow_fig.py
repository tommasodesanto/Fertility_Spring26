from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SOURCE_CSV = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/"
    "Fertility_Spring26_fable_size_mapping_audit_20260701/"
    "output/m1_buy_timing/buy_flows_by_age.csv"
)
OUTPUT_PNG = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/"
    "latex/July_26_slides_buy_flow_by_age_children.png"
)

TITLE_TEXT = "Home purchases by age and number of children"
LEGEND_TEXT = ["no children", "one child", "two children"]


def main() -> None:
    data = pd.read_csv(SOURCE_CSV)
    if data.shape[1] <= 12:
        raise ValueError("Input file must contain at least 13 columns.")

    ages = data["age"]
    total_flow = data["buy_flow"]
    shares = [data.iloc[:, idx] for idx in (10, 11, 12)]
    parts = [total_flow * share for share in shares]

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)

    lower = pd.Series(0.0, index=data.index)
    for values, label in zip(parts, LEGEND_TEXT):
        ax.bar(ages, values, bottom=lower, width=0.82, label=label)
        lower = lower + values

    ax.set_title(TITLE_TEXT, fontsize=15, pad=12)
    ax.set_xlabel("Age", fontsize=12)
    ax.set_ylabel("Purchase flow", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(frameon=False, fontsize=10)
    ax.margins(x=0.01)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
