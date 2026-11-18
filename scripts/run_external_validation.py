"""Run independent target-held-out and frozen-model validation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from protein_binder_generalization.calibration import (
    calibration_bins,
    expected_calibration_error,
    target_calibration_summary,
)
from protein_binder_generalization.external_validation import (
    evaluate_direct_rankings,
    evaluate_external_loto,
    evaluate_frozen_primary_model,
    normalize_target_name,
    paired_method_bootstrap,
    prepare_external_features,
    summarize_external_metrics,
)
from protein_binder_generalization.features import build_feature_table

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DATA = (
    ROOT
    / "data"
    / "raw"
    / "external"
    / "overath_2025"
    / "final_dataset.csv"
)
PRIMARY_DATA = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
TABLES = ROOT / "reports" / "tables"
FIGURES = ROOT / "reports" / "figures"


def plot_method_comparison(summary: pd.DataFrame) -> None:
    data = summary[summary["cohort"].eq("all_targets")].copy()
    data = data.sort_values("mean_average_precision")
    colors = [
        "#4C78A8" if method == "AF3 ipSAE" else "#A7A7A7"
        for method in data["method"]
    ]
    figure, axis = plt.subplots(figsize=(8.8, 6.0))
    axis.barh(data["method"], data["mean_average_precision"], color=colors)
    axis.set_xlabel("Mean average precision across 15 targets")
    axis.set_title("AF3 ipSAE is the strongest external ranking baseline", loc="left")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(axis="x", alpha=0.25)
    for bar, value in zip(
        axis.patches, data["mean_average_precision"], strict=True
    ):
        axis.text(
            value + 0.006,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.3f}",
            va="center",
            fontsize=9,
        )
    axis.set_xlim(0, max(0.58, data["mean_average_precision"].max() + 0.06))
    figure.tight_layout()
    figure.savefig(
        FIGURES / "external_validation_method_comparison.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_external_calibration(calibration: pd.DataFrame) -> None:
    data = calibration[calibration["method"].eq("External LOTO logistic")]
    figure, axis = plt.subplots(figsize=(8.2, 6.8))
    axis.scatter(
        data["mean_predicted_probability"],
        data["observed_rate"],
        s=70,
        color="#4C78A8",
        edgecolor="white",
        linewidth=0.7,
    )
    axis.plot([0, 0.65], [0, 0.65], linestyle="--", color="#555555")
    labeled = data.nlargest(8, "absolute_calibration_bias")
    for row in labeled.itertuples():
        axis.annotate(
            row.target,
            (row.mean_predicted_probability, row.observed_rate),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )
    axis.set_xlim(-0.01, 0.65)
    axis.set_ylim(-0.01, 0.65)
    axis.set_xlabel("Mean predicted binding probability")
    axis.set_ylabel("Observed binder rate")
    axis.set_title(
        "Independent data reproduce target-level calibration errors", loc="left"
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(
        FIGURES / "external_target_calibration.png", dpi=200, bbox_inches="tight"
    )
    plt.close(figure)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    external = prepare_external_features(
        pd.read_csv(EXTERNAL_DATA, low_memory=False)
    )
    primary = build_feature_table(pd.read_parquet(PRIMARY_DATA))
    primary_target_names = {
        normalize_target_name(target) for target in primary["target"].unique()
    }
    external["seen_in_primary"] = external["normalized_target"].isin(
        primary_target_names
    )

    direct_folds = evaluate_direct_rankings(external)
    loto_folds, loto_predictions = evaluate_external_loto(external)
    frozen_folds, frozen_predictions = evaluate_frozen_primary_model(
        primary, external
    )
    folds = pd.concat([direct_folds, loto_folds, frozen_folds], ignore_index=True)
    predictions = pd.concat(
        [loto_predictions, frozen_predictions], ignore_index=True
    )
    summary = summarize_external_metrics(folds, external)
    bootstrap = pd.DataFrame(
        [
            paired_method_bootstrap(
                folds, "AF3 ipSAE", "Frozen primary logistic"
            ),
            paired_method_bootstrap(
                folds, "AF3 ipSAE", "External LOTO logistic"
            ),
            paired_method_bootstrap(
                folds, "Frozen primary logistic", "External LOTO logistic"
            ),
            paired_method_bootstrap(folds, "Consensus mean", "AF3 ipSAE"),
        ]
    )

    calibration_rows = []
    bin_rows = []
    for method, group in predictions.groupby("method"):
        target_summary = target_calibration_summary(group)
        target_summary.insert(0, "method", method)
        calibration_rows.append(target_summary)
        bins = calibration_bins(
            group["binder_final"], group["predicted_probability"], n_bins=10
        )
        bins.insert(0, "method", method)
        bins["expected_calibration_error"] = expected_calibration_error(bins)
        bin_rows.append(bins)
    calibration = pd.concat(calibration_rows, ignore_index=True)
    calibration_bin_table = pd.concat(bin_rows, ignore_index=True)

    folds.to_csv(TABLES / "external_validation_fold_metrics.csv", index=False)
    predictions.to_csv(TABLES / "external_validation_predictions.csv", index=False)
    summary.to_csv(TABLES / "external_validation_summary.csv", index=False)
    bootstrap.to_csv(TABLES / "external_validation_bootstrap.csv", index=False)
    calibration.to_csv(TABLES / "external_target_calibration.csv", index=False)
    calibration_bin_table.to_csv(
        TABLES / "external_calibration_bins.csv", index=False
    )
    plot_method_comparison(summary)
    plot_external_calibration(calibration)

    print("External validation summary")
    print(summary.round(3).to_string(index=False))
    print("\nCalibration summary")
    print(
        calibration.groupby("method")
        .agg(
            mean_absolute_target_bias=("absolute_calibration_bias", "mean"),
            median_absolute_target_bias=("absolute_calibration_bias", "median"),
        )
        .round(3)
        .to_string()
    )
    print("\nPooled expected calibration error")
    print(
        calibration_bin_table.groupby("method")["expected_calibration_error"]
        .first()
        .round(3)
        .to_string()
    )
    print("\nPaired target bootstrap")
    print(bootstrap.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
