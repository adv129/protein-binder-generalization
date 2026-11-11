import math

import numpy as np
import pandas as pd

from protein_binder_generalization.evaluation import fold_metrics


def test_fold_metrics_use_prevalence_as_random_ap_reference() -> None:
    metrics = fold_metrics(pd.Series([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9]))
    assert metrics["random_ap_reference"] == 0.5
    assert metrics["average_precision"] == 1.0
    assert metrics["normalized_ap_lift"] == 1.0


def test_average_precision_is_undefined_without_positive_examples() -> None:
    metrics = fold_metrics(pd.Series([0, 0, 0]), np.array([0.2, 0.1, 0.3]))
    assert math.isnan(metrics["average_precision"])
    assert math.isnan(metrics["normalized_ap_lift"])
