import numpy as np
import pandas as pd

from protein_binder_generalization.selection import (
    candidate_selection_table,
    summarize_candidate_selection,
)


def test_candidate_selection_uses_probability_ranking() -> None:
    predictions = pd.DataFrame(
        {
            "target": ["T", "T", "T", "T"],
            "binder_final": [0, 1, 0, 1],
            "predicted_probability": [0.1, 0.9, 0.2, 0.8],
        }
    )
    result = candidate_selection_table(predictions, budgets=[1, 2])
    assert result["model_hits"].tolist() == [1, 2]
    assert result["oracle_hits"].tolist() == [1, 2]
    assert np.allclose(result["random_expected_hits"], [0.5, 1.0])


def test_selection_summary_keeps_all_targets_and_sensitivity_cohort() -> None:
    predictions = pd.DataFrame(
        {
            "target": ["T1"] * 6 + ["T2"] * 6,
            "binder_final": [1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
            "predicted_probability": np.linspace(1.0, 0.1, 12),
        }
    )
    selection = candidate_selection_table(predictions, budgets=[1])
    summary = summarize_candidate_selection(selection)
    target_counts = summary.set_index("cohort")["targets"].to_dict()
    assert target_counts["all_measured_targets"] == 2
    assert target_counts["targets_with_at_least_5_binders"] == 1
