"""Utilities for independent validation on the Overath et al. release."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from protein_binder_generalization.evaluation import fold_metrics
from protein_binder_generalization.features import CONSENSUS_LEVEL_FEATURES

RANDOM_STATE = 398

EXTERNAL_IPSAE_COLUMNS = (
    "af3_ipSAE_min",
    "boltz1_ipSAE_min",
    "colab_ipSAE_min",
)

DIRECT_RANKINGS = {
    "AF3 ipSAE": "af3_ipSAE_min",
    "Boltz-1 ipSAE": "boltz1_ipSAE_min",
    "ColabFold ipSAE": "colab_ipSAE_min",
    "Consensus mean": "ipsae_mean",
    "Consensus median": "ipsae_median",
    "Consensus minimum": "ipsae_min_across_models",
}


def normalize_target_name(name: str) -> str:
    """Normalize superficial punctuation/case differences in target names."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def prepare_external_features(designs: pd.DataFrame) -> pd.DataFrame:
    """Create comparable three-predictor confidence summaries and labels."""
    required = {
        "binder_id",
        "target_id",
        "source",
        "binder",
        *EXTERNAL_IPSAE_COLUMNS,
    }
    missing = required.difference(designs.columns)
    if missing:
        raise KeyError(f"Required external columns are absent: {sorted(missing)}")

    labeled = designs[designs["binder"].notna()].copy()
    labeled["binder_final"] = labeled["binder"].astype(bool).astype(int)
    labeled["target"] = labeled["target_id"]
    labeled["uuid"] = labeled["binder_id"]
    labeled["normalized_target"] = labeled["target"].map(normalize_target_name)

    scores = labeled.loc[:, EXTERNAL_IPSAE_COLUMNS].astype(float)
    labeled["ipsae_mean"] = scores.mean(axis=1)
    labeled["ipsae_median"] = scores.median(axis=1)
    labeled["ipsae_min_across_models"] = scores.min(axis=1)
    labeled["ipsae_max_across_models"] = scores.max(axis=1)
    labeled["ipsae_std"] = scores.std(axis=1)
    labeled["ipsae_range"] = scores.max(axis=1) - scores.min(axis=1)
    return labeled


def logistic_model() -> Pipeline:
    """Return the same fixed regularized model used in the primary analysis."""
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


def evaluate_direct_rankings(features: pd.DataFrame) -> pd.DataFrame:
    """Evaluate released scores directly without fitting a model."""
    rows = []
    for target, group in features.groupby("target"):
        for method, column in DIRECT_RANKINGS.items():
            ranking_score = group[column].fillna(0.0).to_numpy()
            metrics = fold_metrics(group["binder_final"], ranking_score)
            metrics["brier_score"] = np.nan
            rows.append({"held_out_target": target, "method": method, **metrics})
    return pd.DataFrame(rows)


def evaluate_external_loto(
    features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit the fixed model on all external targets except the test target."""
    feature_names = list(CONSENSUS_LEVEL_FEATURES)
    fold_rows = []
    predictions = []
    for held_out_target in sorted(features["target"].unique()):
        train = features[features["target"].ne(held_out_target)]
        test = features[features["target"].eq(held_out_target)]
        model = logistic_model()
        model.fit(train[feature_names], train["binder_final"])
        probabilities = model.predict_proba(test[feature_names])[:, 1]
        fold_rows.append(
            {
                "held_out_target": held_out_target,
                "method": "External LOTO logistic",
                **fold_metrics(test["binder_final"], probabilities),
            }
        )
        predictions.append(
            pd.DataFrame(
                {
                    "uuid": test["uuid"].to_numpy(),
                    "target": held_out_target,
                    "binder_final": test["binder_final"].to_numpy(),
                    "method": "External LOTO logistic",
                    "predicted_probability": probabilities,
                }
            )
        )
    return pd.DataFrame(fold_rows), pd.concat(predictions, ignore_index=True)


def evaluate_frozen_primary_model(
    primary_features: pd.DataFrame,
    external_features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Train once on the primary release and apply unchanged to external data."""
    feature_names = list(CONSENSUS_LEVEL_FEATURES)
    train = primary_features[primary_features["binder_final"].notna()].copy()
    train["binder_final"] = train["binder_final"].astype(int)
    model = logistic_model()
    model.fit(train[feature_names], train["binder_final"])
    probabilities = model.predict_proba(external_features[feature_names])[:, 1]
    predictions = external_features[
        ["uuid", "target", "binder_final"]
    ].copy()
    predictions["method"] = "Frozen primary logistic"
    predictions["predicted_probability"] = probabilities

    fold_rows = []
    for target, group in predictions.groupby("target"):
        fold_rows.append(
            {
                "held_out_target": target,
                "method": "Frozen primary logistic",
                **fold_metrics(
                    group["binder_final"], group["predicted_probability"]
                ),
            }
        )
    return pd.DataFrame(fold_rows), predictions


def summarize_external_metrics(
    folds: pd.DataFrame,
    target_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Macro-average methods over all, sufficiently positive, and novel targets."""
    metadata = target_metadata.drop_duplicates("target")[
        ["target", "normalized_target", "seen_in_primary"]
    ]
    annotated = folds.merge(
        metadata, left_on="held_out_target", right_on="target", validate="many_to_one"
    )
    cohorts = {
        "all_targets": annotated,
        "targets_with_at_least_5_binders": annotated[annotated["n_binders"].ge(5)],
        "targets_not_seen_in_primary": annotated[
            annotated["seen_in_primary"].eq(False)
        ],
    }
    rows = []
    for cohort_name, cohort in cohorts.items():
        summary = (
            cohort.groupby("method", as_index=False)
            .agg(
                targets=("held_out_target", "nunique"),
                mean_average_precision=("average_precision", "mean"),
                median_average_precision=("average_precision", "median"),
                mean_normalized_ap_lift=("normalized_ap_lift", "mean"),
                mean_brier_score=("brier_score", "mean"),
                mean_top_10_hit_rate=("top_10_hit_rate", "mean"),
            )
            .sort_values("mean_average_precision", ascending=False)
        )
        summary.insert(0, "cohort", cohort_name)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def paired_method_bootstrap(
    folds: pd.DataFrame,
    method_a: str,
    method_b: str,
    n_bootstrap: int = 20_000,
) -> dict[str, float | int | str]:
    """Bootstrap the paired target-level AP difference between two methods."""
    pivot = folds.pivot(
        index="held_out_target", columns="method", values="average_precision"
    )
    differences = (pivot[method_a] - pivot[method_b]).dropna().to_numpy()
    rng = np.random.default_rng(RANDOM_STATE)
    sampled_indices = rng.integers(
        0, len(differences), size=(n_bootstrap, len(differences))
    )
    bootstrap_means = differences[sampled_indices].mean(axis=1)
    return {
        "method_a": method_a,
        "method_b": method_b,
        "targets": len(differences),
        "mean_ap_difference": float(differences.mean()),
        "ci_2_5": float(np.quantile(bootstrap_means, 0.025)),
        "ci_97_5": float(np.quantile(bootstrap_means, 0.975)),
        "target_wins_for_a": int((differences > 0).sum()),
        "bootstrap_probability_a_better": float((bootstrap_means > 0).mean()),
    }
