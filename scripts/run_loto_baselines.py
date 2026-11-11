"""Run fixed leave-one-target-out baseline models without test-set tuning."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from protein_binder_generalization.evaluation import fold_metrics
from protein_binder_generalization.features import FEATURE_GROUPS, build_feature_table

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
TABLES = ROOT / "reports" / "tables"
FIGURES = ROOT / "reports" / "figures"
RANDOM_STATE = 398


def logistic_model() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(C=1.0, max_iter=5000, random_state=RANDOM_STATE),
            ),
        ]
    )


def boosted_tree_model() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
            (
                "classifier",
                HistGradientBoostingClassifier(
                    learning_rate=0.05,
                    max_iter=200,
                    max_leaf_nodes=15,
                    l2_regularization=1.0,
                    random_state=RANDOM_STATE,
                ),
            ),
        ]
    )


MODELS = {
    "logistic_regression": logistic_model,
    "gradient_boosted_trees": boosted_tree_model,
}


def run_evaluation(feature_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    labeled = feature_table[feature_table["binder_final"].notna()].copy()
    labeled["binder_final"] = labeled["binder_final"].astype(int)
    fold_rows: list[dict] = []
    prediction_rows: list[pd.DataFrame] = []

    for held_out_target in sorted(labeled["target"].unique()):
        train = labeled[labeled["target"] != held_out_target]
        test = labeled[labeled["target"] == held_out_target]
        if set(train["target"]).intersection(test["target"]):
            raise RuntimeError("Target leakage detected in a held-out fold.")

        for feature_group, feature_names in FEATURE_GROUPS.items():
            x_train = train[feature_names]
            y_train = train["binder_final"]
            x_test = test[feature_names]
            y_test = test["binder_final"]

            for model_name, model_factory in MODELS.items():
                model = model_factory()
                model.fit(x_train, y_train)
                probability = model.predict_proba(x_test)[:, 1]
                metrics = fold_metrics(y_test, probability)
                fold_rows.append(
                    {
                        "held_out_target": held_out_target,
                        "feature_group": feature_group,
                        "model": model_name,
                        **metrics,
                    }
                )
                prediction_rows.append(
                    pd.DataFrame(
                        {
                            "uuid": test["uuid"].to_numpy(),
                            "target": held_out_target,
                            "binder_final": y_test.to_numpy(),
                            "feature_group": feature_group,
                            "model": model_name,
                            "predicted_probability": probability,
                        }
                    )
                )

    return pd.DataFrame(fold_rows), pd.concat(prediction_rows, ignore_index=True)


def summarize(folds: pd.DataFrame) -> pd.DataFrame:
    return (
        folds.groupby(["feature_group", "model"])
        .agg(
            targets=("held_out_target", "nunique"),
            evaluable_ap_targets=("average_precision", "count"),
            mean_average_precision=("average_precision", "mean"),
            median_average_precision=("average_precision", "median"),
            mean_normalized_ap_lift=("normalized_ap_lift", "mean"),
            mean_brier_score=("brier_score", "mean"),
            mean_top_10_hit_rate=("top_10_hit_rate", "mean"),
            mean_top_10_enrichment=("top_10_enrichment", "mean"),
        )
        .sort_values("mean_average_precision", ascending=False)
    )


def paired_target_bootstrap(
    folds: pd.DataFrame,
    feature_group_a: str = "combined",
    feature_group_b: str = "structure",
    model: str = "logistic_regression",
    minimum_binders: int = 5,
    n_bootstrap: int = 20_000,
) -> pd.DataFrame:
    """Bootstrap the paired target-level AP difference between two feature sets."""
    eligible = folds[
        folds["model"].eq(model) & folds["n_binders"].ge(minimum_binders)
    ]
    pivot = eligible.pivot(
        index="held_out_target", columns="feature_group", values="average_precision"
    )
    differences = (pivot[feature_group_a] - pivot[feature_group_b]).dropna().to_numpy()
    rng = np.random.default_rng(RANDOM_STATE)
    bootstrap_means = np.array(
        [
            rng.choice(differences, size=len(differences), replace=True).mean()
            for _ in range(n_bootstrap)
        ]
    )
    return pd.DataFrame(
        [
            {
                "model": model,
                "feature_group_a": feature_group_a,
                "feature_group_b": feature_group_b,
                "minimum_binders": minimum_binders,
                "targets": len(differences),
                "mean_ap_difference": differences.mean(),
                "ci_2_5": np.quantile(bootstrap_means, 0.025),
                "ci_97_5": np.quantile(bootstrap_means, 0.975),
                "bootstrap_probability_a_better": (bootstrap_means > 0).mean(),
            }
        ]
    )


def plot_model_summary(summary: pd.DataFrame) -> None:
    plot_data = summary.reset_index().copy()
    plot_data["label"] = (
        plot_data["model"].str.replace("_", " ")
        + " — "
        + plot_data["feature_group"]
    )
    plot_data = plot_data.sort_values("mean_average_precision")

    figure, axis = plt.subplots(figsize=(9, 5.2))
    axis.barh(
        plot_data["label"],
        plot_data["mean_average_precision"],
        color="#4C78A8",
    )
    axis.set_xlabel("Mean average precision across held-out targets")
    axis.set_ylabel("")
    axis.set_title(
        "Leave-one-target-out comparison (targets with at least 5 binders)",
        loc="left",
    )
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(axis="x", alpha=0.25)
    figure.tight_layout()
    figure.savefig(
        FIGURES / "loto_model_comparison.png", dpi=200, bbox_inches="tight"
    )
    plt.close(figure)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    designs = pd.read_parquet(DATA)
    feature_table = build_feature_table(designs)
    folds, predictions = run_evaluation(feature_table)
    summary = summarize(folds)
    sensitivity = summarize(folds[folds["n_binders"].ge(5)])
    bootstrap = paired_target_bootstrap(folds)

    schema_rows = [
        {"feature_group": group, "feature": feature}
        for group, features in FEATURE_GROUPS.items()
        for feature in features
    ]
    pd.DataFrame(schema_rows).to_csv(TABLES / "feature_matrix_schema.csv", index=False)
    folds.to_csv(TABLES / "loto_fold_metrics.csv", index=False)
    predictions.to_csv(TABLES / "loto_predictions.csv", index=False)
    summary.to_csv(TABLES / "loto_model_summary.csv")
    sensitivity.to_csv(TABLES / "loto_model_summary_min5_binders.csv")
    bootstrap.to_csv(TABLES / "loto_feature_group_bootstrap.csv", index=False)
    plot_model_summary(sensitivity)

    print("All labeled targets")
    print(summary.round(3).to_string())
    print("\nSensitivity analysis: targets with at least five binders")
    print(sensitivity.round(3).to_string())
    print("\nPaired target bootstrap")
    print(bootstrap.round(3).to_string(index=False))
    print(f"\nEvaluated {folds['held_out_target'].nunique()} held-out targets.")


if __name__ == "__main__":
    main()
