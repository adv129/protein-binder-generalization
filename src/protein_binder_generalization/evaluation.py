"""Target-held-out evaluation utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss


def fold_metrics(y_true: pd.Series, probability: np.ndarray, top_k: int = 10) -> dict:
    """Calculate decision-focused metrics for one held-out target."""
    y = y_true.astype(int).to_numpy()
    probability = np.asarray(probability, dtype=float)
    order = np.argsort(-probability)
    top = y[order[: min(top_k, len(y))]]
    prevalence = float(y.mean())
    has_both_classes = np.unique(y).size == 2

    return {
        "n_test": len(y),
        "n_binders": int(y.sum()),
        "prevalence": prevalence,
        "average_precision": (
            float(average_precision_score(y, probability))
            if has_both_classes
            else np.nan
        ),
        "random_ap_reference": prevalence if has_both_classes else np.nan,
        "normalized_ap_lift": (
            float((average_precision_score(y, probability) - prevalence) / (1 - prevalence))
            if has_both_classes and prevalence < 1
            else np.nan
        ),
        "brier_score": float(brier_score_loss(y, probability)),
        "top_10_hit_rate": float(top.mean()),
        "top_10_enrichment": (
            float(top.mean() / prevalence) if prevalence > 0 else np.nan
        ),
    }
