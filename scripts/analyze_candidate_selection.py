"""Evaluate held-out-target rankings at fixed experimental budgets."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from protein_binder_generalization.selection import (
    bootstrap_gain_over_random,
    candidate_selection_table,
    summarize_candidate_selection,
)

ROOT = Path(__file__).resolve().parents[1]
PREDICTIONS = ROOT / "reports" / "tables" / "loto_predictions.csv"
TABLES = ROOT / "reports" / "tables"
FIGURES = ROOT / "reports" / "figures"
FEATURE_GROUP = "consensus_level"
MODEL = "logistic_regression"
BUDGETS = (1, 3, 5, 10, 20, 30)


def primary_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    selected = predictions[
        predictions["feature_group"].eq(FEATURE_GROUP)
        & predictions["model"].eq(MODEL)
    ].copy()
    if selected.empty:
        raise ValueError("Primary model predictions were not found.")
    return selected


def plot_budget_curve(summary: pd.DataFrame) -> None:
    data = summary[summary["cohort"].eq("all_measured_targets")].sort_values(
        "requested_budget"
    )
    budgets = data["requested_budget"].to_numpy()

    figure, axis = plt.subplots(figsize=(9.2, 6.3))
    axis.fill_between(
        budgets,
        data["q25_model_hits"].to_numpy(),
        data["q75_model_hits"].to_numpy(),
        color="#4C78A8",
        alpha=0.18,
        label="Middle 50% of targets",
    )
    axis.plot(
        budgets,
        data["mean_model_hits"],
        marker="o",
        linewidth=2.4,
        color="#4C78A8",
        label="Model-ranked candidates",
    )
    axis.plot(
        budgets,
        data["mean_random_expected_hits"],
        marker="o",
        linestyle="--",
        color="#E45756",
        label="Random-selection expectation",
    )
    axis.plot(
        budgets,
        data["mean_oracle_hits"],
        linestyle=":",
        linewidth=2.2,
        color="#555555",
        label="Oracle ceiling",
    )
    axis.set_xticks(budgets)
    axis.set_xlabel("Number of candidates tested per target")
    axis.set_ylabel("Mean experimentally confirmed binders")
    axis.set_title(
        "Structure-based ranking improves yield, but performance varies by target",
        loc="left",
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, loc="upper left")
    figure.tight_layout()
    figure.savefig(
        FIGURES / "candidate_selection_budget_curve.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    predictions = primary_predictions(pd.read_csv(PREDICTIONS))
    selection = candidate_selection_table(predictions, budgets=BUDGETS)
    summary = summarize_candidate_selection(selection)
    bootstrap = pd.concat(
        [
            bootstrap_gain_over_random(selection),
            bootstrap_gain_over_random(selection, minimum_binders=5),
        ],
        ignore_index=True,
    )

    selection.to_csv(TABLES / "candidate_selection_by_target.csv", index=False)
    summary.to_csv(TABLES / "candidate_selection_summary.csv", index=False)
    bootstrap.to_csv(TABLES / "candidate_selection_bootstrap.csv", index=False)
    plot_budget_curve(summary)

    print("All measured targets")
    print(
        summary[summary["cohort"].eq("all_measured_targets")]
        .round(3)
        .to_string(index=False)
    )
    print("\nTargets with at least five binders")
    print(
        summary[summary["cohort"].eq("targets_with_at_least_5_binders")]
        .round(3)
        .to_string(index=False)
    )
    print("\nBootstrap gain over random selection")
    print(bootstrap.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
