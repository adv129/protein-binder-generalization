"""Analyze target-level calibration and candidate-ranking failure cases."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from protein_binder_generalization.calibration import (
    calibration_bins,
    expected_calibration_error,
    select_failure_cases,
    target_calibration_summary,
    top_k_summary,
)
from protein_binder_generalization.features import build_feature_table

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
PREDICTIONS = ROOT / "reports" / "tables" / "loto_predictions.csv"
TABLES = ROOT / "reports" / "tables"
FIGURES = ROOT / "reports" / "figures"
FEATURE_GROUP = "consensus_level"
MODEL = "logistic_regression"


def primary_predictions(all_predictions: pd.DataFrame) -> pd.DataFrame:
    selected = all_predictions[
        all_predictions["feature_group"].eq(FEATURE_GROUP)
        & all_predictions["model"].eq(MODEL)
    ].copy()
    if selected.empty:
        raise ValueError("Primary model predictions were not found.")
    return selected


def plot_target_calibration(summary: pd.DataFrame) -> None:
    label_offsets = {
        "15-PGDH": (-48, -5),
        "BBF-14": (6, 18),
        "MBP": (6, 5),
        "TNFa": (6, -10),
        "Cas9": (6, 7),
        "EGFR": (6, 7),
    }
    figure, axis = plt.subplots(figsize=(8.2, 7.0))
    axis.scatter(
        summary["mean_predicted_probability"],
        summary["observed_rate"],
        s=70,
        color="#4C78A8",
        edgecolor="white",
        linewidth=0.7,
    )
    axis.plot([0, 1], [0, 1], linestyle="--", color="#555555", linewidth=1.3)
    for row in summary.itertuples():
        offset = label_offsets.get(row.target, (5, 4))
        axis.annotate(
            row.target,
            (row.mean_predicted_probability, row.observed_rate),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
        )
    axis.set_xlim(-0.02, 1.0)
    axis.set_ylim(-0.02, 1.0)
    axis.set_xlabel("Mean predicted binding probability")
    axis.set_ylabel("Observed binder rate")
    axis.set_title("Calibration varies substantially across held-out targets", loc="left")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(FIGURES / "target_calibration.png", dpi=200, bbox_inches="tight")
    plt.close(figure)


def plot_reliability(bins: pd.DataFrame) -> None:
    figure, axis = plt.subplots(figsize=(7.5, 6.2))
    axis.plot(
        bins["mean_probability"],
        bins["observed_rate"],
        marker="o",
        linewidth=2,
        color="#4C78A8",
    )
    axis.plot([0, 1], [0, 1], linestyle="--", color="#555555", linewidth=1.3)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Mean predicted probability")
    axis.set_ylabel("Observed binder frequency")
    axis.set_title("Out-of-fold reliability across all targets", loc="left")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(FIGURES / "pooled_reliability_curve.png", dpi=200, bbox_inches="tight")
    plt.close(figure)


def plot_top_k(summary: pd.DataFrame) -> None:
    ordered = summary.sort_values("top_k_hit_rate")
    y = np.arange(len(ordered))
    figure, axis = plt.subplots(figsize=(9.2, 6.8))
    axis.hlines(
        y,
        ordered["target_prevalence"],
        ordered["top_k_hit_rate"],
        color="#BBBBBB",
        linewidth=2,
    )
    axis.scatter(
        ordered["target_prevalence"],
        y,
        color="#777777",
        label="Overall target binder rate",
        zorder=3,
    )
    axis.scatter(
        ordered["top_k_hit_rate"],
        y,
        color="#4C78A8",
        label="Top-10 hit rate",
        zorder=3,
    )
    axis.set_yticks(y, ordered["target"])
    axis.set_xlim(-0.02, 1.02)
    axis.set_xlabel("Fraction experimentally confirmed as binders")
    axis.set_title("Top-10 selection improves some targets and fails on others", loc="left")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(axis="x", alpha=0.2)
    axis.legend(frameon=False, loc="lower right")
    figure.tight_layout()
    figure.savefig(FIGURES / "top10_hit_rate_by_target.png", dpi=200, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    designs = pd.read_parquet(DATA)
    predictions = primary_predictions(pd.read_csv(PREDICTIONS))

    target_summary = target_calibration_summary(predictions)
    top10 = top_k_summary(predictions, k=10)
    bins = calibration_bins(
        predictions["binder_final"], predictions["predicted_probability"], n_bins=10
    )
    cases = select_failure_cases(predictions, cases_per_target=3)

    feature_table = build_feature_table(designs)
    descriptive = designs[["uuid", "full_name", "generator", "rank"]]
    case_details = cases.merge(feature_table, on=["uuid", "target", "binder_final"])
    case_details = case_details.merge(descriptive, on="uuid")
    case_details = case_details.sort_values(
        ["target", "case_type", "predicted_probability"],
        ascending=[True, True, False],
    )

    bins.to_csv(TABLES / "calibration_bins.csv", index=False)
    target_summary.to_csv(TABLES / "target_calibration.csv", index=False)
    top10.to_csv(TABLES / "top10_target_summary.csv", index=False)
    case_details.to_csv(TABLES / "failure_cases.csv", index=False)

    plot_target_calibration(target_summary)
    plot_reliability(bins)
    plot_top_k(top10)

    print(f"Pooled expected calibration error: {expected_calibration_error(bins):.3f}")
    print(
        "Mean absolute target-level calibration bias: "
        f"{target_summary['absolute_calibration_bias'].mean():.3f}"
    )
    print("\nLargest target-level calibration errors")
    print(
        target_summary[
            [
                "target",
                "observed_rate",
                "mean_predicted_probability",
                "calibration_bias",
                "brier_score",
            ]
        ]
        .head(10)
        .round(3)
        .to_string(index=False)
    )
    print("\nTop-10 results")
    print(top10.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
