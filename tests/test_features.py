import pandas as pd
import pytest

from protein_binder_generalization.features import (
    assert_no_leakage,
    sequence_features,
)


def test_sequence_features_for_simple_sequence() -> None:
    features = sequence_features("AKDE")
    assert features["seq_length"] == 4
    assert features["seq_positive_fraction"] == 0.25
    assert features["seq_negative_fraction"] == 0.50
    assert features["seq_net_charge_proxy"] == -0.25


def test_experimental_variables_are_rejected() -> None:
    with pytest.raises(ValueError, match="Post-experimental"):
        assert_no_leakage(pd.Index(["ipsae_min_chai1", "kd_nM_final"]))
