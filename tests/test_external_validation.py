import numpy as np
import pandas as pd

from protein_binder_generalization.external_validation import (
    normalize_target_name,
    prepare_external_features,
)


def test_normalize_target_name_ignores_case_and_punctuation() -> None:
    assert normalize_target_name("IL-7Ra") == normalize_target_name("IL7RA")
    assert normalize_target_name("PD-L1") == normalize_target_name("Pdl1")


def test_prepare_external_features_drops_unlabeled_and_builds_consensus() -> None:
    raw = pd.DataFrame(
        {
            "binder_id": ["a", "b"],
            "target_id": ["T", "T"],
            "source": ["study", "study"],
            "binder": [True, np.nan],
            "af3_ipSAE_min": [0.9, 0.1],
            "boltz1_ipSAE_min": [0.6, 0.2],
            "colab_ipSAE_min": [0.3, 0.3],
        }
    )
    result = prepare_external_features(raw)
    assert result["uuid"].tolist() == ["a"]
    assert result["binder_final"].tolist() == [1]
    assert np.isclose(result.loc[result.index[0], "ipsae_mean"], 0.6)
    assert np.isclose(result.loc[result.index[0], "ipsae_range"], 0.6)
