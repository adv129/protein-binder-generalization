"""Test whether generator/provenance metadata explains cross-target performance."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from protein_binder_generalization.features import build_feature_table
from protein_binder_generalization.provenance import (
    FEATURE_GROUPS,
    add_provenance,
    paired_ap_bootstrap,
    run_provenance_evaluation,
    summarize_provenance,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
TABLES = ROOT / "reports" / "tables"
FIGURES = ROOT / "reports" / "figures"


def plot_summary(summary: pd.DataFrame) -> None:
    labels = {
        "consensus": "Structural confidence summaries",
        "generator_only": "Generator identity only",
        "provenance_only": "All design provenance",
        "consensus_plus_generator": "Confidence + generator",
        "consensus_plus_provenance": "Confidence + all provenance",
    }
    plot_data = summary.reset_index().copy()
    plot_data["label"] = plot_data["feature_group"].map(labels)
    plot_data = plot_data.sort_values("mean_average_precision")

    figure, axis = plt.subplots(figsize=(9.2, 5.3))
    colors = [
        "#4C78A8" if group == "consensus" else "#A7A7A7"
        for group in plot_data["feature_group"]
    ]
    axis.barh(plot_data["label"], plot_data["mean_average_precision"], color=colors)
    axis.set_xlabel("Mean average precision across held-out targets")
    axis.set_title("Design provenance does not improve cross-target ranking", loc="left")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(axis="x", alpha=0.25)
    for bar, value in zip(
        axis.patches, plot_data["mean_average_precision"], strict=True
    ):
        axis.text(
            value + 0.006,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.3f}",
            va="center",
            fontsize=9,
        )
    axis.set_xlim(0, max(0.65, plot_data["mean_average_precision"].max() + 0.07))
    figure.tight_layout()
    figure.savefig(
        FIGURES / "design_provenance_comparison.png", dpi=200, bbox_inches="tight"
    )
    plt.close(figure)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    designs = pd.read_parquet(DATA)
    features = add_provenance(build_feature_table(designs), designs)
    folds, predictions = run_provenance_evaluation(features)
    summary = summarize_provenance(folds)
    sensitivity = summarize_provenance(folds[folds["n_binders"].ge(5)])
    bootstrap = pd.DataFrame(
        [
            paired_ap_bootstrap(folds, comparison_group=group)
            for group in FEATURE_GROUPS
            if group != "consensus"
        ]
    )

    folds.to_csv(TABLES / "loto_provenance_fold_metrics.csv", index=False)
    predictions.to_csv(TABLES / "loto_provenance_predictions.csv", index=False)
    summary.to_csv(TABLES / "loto_provenance_summary.csv")
    sensitivity.to_csv(TABLES / "loto_provenance_summary_min5_binders.csv")
    bootstrap.to_csv(TABLES / "loto_provenance_bootstrap.csv", index=False)
    plot_summary(sensitivity)

    print("Sensitivity analysis: targets with at least five binders")
    print(sensitivity.round(3).to_string())
    print("\nPaired target bootstrap relative to structural confidence")
    print(bootstrap.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
