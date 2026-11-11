"""Calibration and failure-case summaries for held-out predictions."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss


def calibration_bins(
    labels: pd.Series | np.ndarray,
    probabilities: pd.Series | np.ndarray,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Create approximately equal-frequency bins for a reliability diagram."""
    frame = pd.DataFrame(
        {
            "label": np.asarray(labels, dtype=int),
            "probability": np.asarray(probabilities, dtype=float),
        }
    ).sort_values("probability")
    if frame.empty:
        raise ValueError("Calibration requires at least one prediction.")
    bin_count = min(n_bins, len(frame))
    frame["bin"] = pd.qcut(
        frame["probability"].rank(method="first"),
        q=bin_count,
        labels=False,
    )
    summary = (
        frame.groupby("bin", as_index=False)
        .agg(
            n=("label", "size"),
            mean_probability=("probability", "mean"),
            observed_rate=("label", "mean"),
            minimum_probability=("probability", "min"),
            maximum_probability=("probability", "max"),
        )
        .sort_values("mean_probability")
    )
    summary["calibration_gap"] = (
        summary["observed_rate"] - summary["mean_probability"]
    )
    summary["absolute_calibration_gap"] = summary["calibration_gap"].abs()
    return summary


def expected_calibration_error(bins: pd.DataFrame) -> float:
    """Return the sample-weighted absolute calibration gap."""
    return float(
        np.average(bins["absolute_calibration_gap"], weights=bins["n"])
    )


def target_calibration_summary(predictions: pd.DataFrame) -> pd.DataFrame:
    """Summarize ranking and probability behavior for every held-out target."""
    rows = []
    for target, group in predictions.groupby("target"):
        labels = group["binder_final"].astype(int)
        probabilities = group["predicted_probability"].astype(float)
        prevalence = labels.mean()
        rows.append(
            {
                "target": target,
                "n": len(group),
                "binders": int(labels.sum()),
                "observed_rate": prevalence,
                "mean_predicted_probability": probabilities.mean(),
                "calibration_bias": probabilities.mean() - prevalence,
                "absolute_calibration_bias": abs(probabilities.mean() - prevalence),
                "brier_score": brier_score_loss(labels, probabilities),
                "average_precision": (
                    average_precision_score(labels, probabilities)
                    if labels.nunique() == 2
                    else np.nan
                ),
                "maximum_predicted_probability": probabilities.max(),
            }
        )
    return pd.DataFrame(rows).sort_values("absolute_calibration_bias", ascending=False)


def top_k_summary(predictions: pd.DataFrame, k: int = 10) -> pd.DataFrame:
    """Summarize the experimental yield of the top-ranked candidates per target."""
    rows = []
    for target, group in predictions.groupby("target"):
        ranked = group.sort_values("predicted_probability", ascending=False)
        selected = ranked.head(min(k, len(ranked)))
        prevalence = group["binder_final"].mean()
        hit_rate = selected["binder_final"].mean()
        rows.append(
            {
                "target": target,
                "k": len(selected),
                "hits": int(selected["binder_final"].sum()),
                "false_positives": int((1 - selected["binder_final"]).sum()),
                "top_k_hit_rate": hit_rate,
                "target_prevalence": prevalence,
                "enrichment": hit_rate / prevalence if prevalence > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("top_k_hit_rate", ascending=False)


def select_failure_cases(
    predictions: pd.DataFrame,
    cases_per_target: int = 3,
) -> pd.DataFrame:
    """Select high-confidence failures and low-confidence successful binders."""
    selected = []
    for target, group in predictions.groupby("target"):
        false_positives = (
            group[group["binder_final"].eq(0)]
            .nlargest(cases_per_target, "predicted_probability")
            .assign(case_type="high_confidence_nonbinder")
        )
        low_confidence_binders = (
            group[group["binder_final"].eq(1)]
            .nsmallest(cases_per_target, "predicted_probability")
            .assign(case_type="low_confidence_binder")
        )
        selected.extend([false_positives, low_confidence_binders])
    return pd.concat(selected, ignore_index=True)
