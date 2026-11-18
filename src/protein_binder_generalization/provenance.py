"""Evaluate whether design-pipeline provenance explains held-out performance."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from protein_binder_generalization.evaluation import fold_metrics
from protein_binder_generalization.features import CONSENSUS_LEVEL_FEATURES

RANDOM_STATE = 398

PROVENANCE_COLUMNS = (
    "generator",
    "design_model",
    "sequence_design_method",
    "campaign",
    "rank",
)

FEATURE_GROUPS: dict[str, tuple[list[str], list[str]]] = {
    "consensus": (list(CONSENSUS_LEVEL_FEATURES), []),
    "generator_only": ([], ["generator"]),
    "provenance_only": (
        ["rank"],
        ["generator", "design_model", "sequence_design_method", "campaign"],
    ),
    "consensus_plus_generator": (
        list(CONSENSUS_LEVEL_FEATURES),
        ["generator"],
    ),
    "consensus_plus_provenance": (
        [*CONSENSUS_LEVEL_FEATURES, "rank"],
        ["generator", "design_model", "sequence_design_method", "campaign"],
    ),
}


def add_provenance(
    feature_table: pd.DataFrame, designs: pd.DataFrame
) -> pd.DataFrame:
    """Join pre-experimental design-pipeline metadata to the feature table."""
    required = {"uuid", *PROVENANCE_COLUMNS}
    missing = required.difference(designs.columns)
    if missing:
        raise KeyError(f"Required provenance columns are absent: {sorted(missing)}")
    if designs["uuid"].duplicated().any():
        raise ValueError("Design UUIDs must be unique before joining provenance.")

    metadata = designs[["uuid", *PROVENANCE_COLUMNS]].copy()
    return feature_table.merge(metadata, on="uuid", how="left", validate="one_to_one")


def provenance_model(
    numeric_features: Sequence[str], categorical_features: Sequence[str]
) -> Pipeline:
    """Construct a fixed logistic model for mixed numeric/categorical inputs."""
    transformers = []
    if numeric_features:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    [
                        (
                            "imputer",
                            SimpleImputer(strategy="median", add_indicator=True),
                        ),
                        ("scaler", StandardScaler()),
                    ]
                ),
                list(numeric_features),
            )
        )
    if categorical_features:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("one_hot", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                list(categorical_features),
            )
        )
    if not transformers:
        raise ValueError("At least one numeric or categorical feature is required.")

    return Pipeline(
        [
            ("preprocessor", ColumnTransformer(transformers)),
            (
                "classifier",
                LogisticRegression(C=1.0, max_iter=5000, random_state=RANDOM_STATE),
            ),
        ]
    )


def run_provenance_evaluation(
    feature_table: pd.DataFrame,
    feature_groups: Mapping[str, tuple[list[str], list[str]]] = FEATURE_GROUPS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate provenance ablations with leave-one-target-out folds."""
    labeled = feature_table[feature_table["binder_final"].notna()].copy()
    labeled["binder_final"] = labeled["binder_final"].astype(int)
    fold_rows: list[dict] = []
    prediction_rows: list[pd.DataFrame] = []

    for held_out_target in sorted(labeled["target"].unique()):
        train = labeled[labeled["target"].ne(held_out_target)]
        test = labeled[labeled["target"].eq(held_out_target)]

        for group_name, (numeric, categorical) in feature_groups.items():
            feature_names = [*numeric, *categorical]
            model = provenance_model(numeric, categorical)
            model.fit(train[feature_names], train["binder_final"])
            probabilities = model.predict_proba(test[feature_names])[:, 1]
            metrics = fold_metrics(test["binder_final"], probabilities)
            fold_rows.append(
                {
                    "held_out_target": held_out_target,
                    "feature_group": group_name,
                    **metrics,
                }
            )
            prediction_rows.append(
                pd.DataFrame(
                    {
                        "uuid": test["uuid"].to_numpy(),
                        "target": held_out_target,
                        "binder_final": test["binder_final"].to_numpy(),
                        "feature_group": group_name,
                        "predicted_probability": probabilities,
                    }
                )
            )

    return pd.DataFrame(fold_rows), pd.concat(prediction_rows, ignore_index=True)


def summarize_provenance(folds: pd.DataFrame) -> pd.DataFrame:
    """Return macro-averaged held-out-target metrics for each feature group."""
    return (
        folds.groupby("feature_group")
        .agg(
            targets=("held_out_target", "nunique"),
            evaluable_ap_targets=("average_precision", "count"),
            mean_average_precision=("average_precision", "mean"),
            median_average_precision=("average_precision", "median"),
            mean_normalized_ap_lift=("normalized_ap_lift", "mean"),
            mean_brier_score=("brier_score", "mean"),
            mean_top_10_hit_rate=("top_10_hit_rate", "mean"),
        )
        .sort_values("mean_average_precision", ascending=False)
    )


def paired_ap_bootstrap(
    folds: pd.DataFrame,
    comparison_group: str,
    reference_group: str = "consensus",
    minimum_binders: int = 5,
    n_bootstrap: int = 20_000,
) -> dict[str, float | int | str]:
    """Bootstrap paired target AP differences for one provenance ablation."""
    eligible = folds[folds["n_binders"].ge(minimum_binders)]
    pivot = eligible.pivot(
        index="held_out_target", columns="feature_group", values="average_precision"
    )
    differences = (
        pivot[comparison_group] - pivot[reference_group]
    ).dropna().to_numpy()
    rng = np.random.default_rng(RANDOM_STATE)
    sampled_indices = rng.integers(
        0, len(differences), size=(n_bootstrap, len(differences))
    )
    bootstrap_means = differences[sampled_indices].mean(axis=1)
    return {
        "comparison_group": comparison_group,
        "reference_group": reference_group,
        "minimum_binders": minimum_binders,
        "targets": len(differences),
        "mean_ap_difference": float(differences.mean()),
        "ci_2_5": float(np.quantile(bootstrap_means, 0.025)),
        "ci_97_5": float(np.quantile(bootstrap_means, 0.975)),
        "bootstrap_probability_comparison_better": float(
            (bootstrap_means > 0).mean()
        ),
    }
