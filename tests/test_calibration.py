import numpy as np
import pandas as pd

from protein_binder_generalization.calibration import (
    calibration_bins,
    expected_calibration_error,
    select_failure_cases,
)


def test_perfectly_calibrated_two_bin_example() -> None:
    bins = calibration_bins([0, 0, 1, 1], [0.0, 0.0, 1.0, 1.0], n_bins=2)
    assert expected_calibration_error(bins) == 0.0


def test_failure_cases_select_highest_nonbinder_and_lowest_binder() -> None:
    predictions = pd.DataFrame(
        {
            "uuid": ["a", "b", "c", "d"],
            "target": ["T", "T", "T", "T"],
            "binder_final": [0, 0, 1, 1],
            "predicted_probability": [0.8, 0.4, 0.2, 0.9],
        }
    )
    selected = select_failure_cases(predictions, cases_per_target=1)
    assert set(selected["uuid"]) == {"a", "c"}
    assert np.isclose(selected["predicted_probability"].sum(), 1.0)
