import pandas as pd

from protein_binder_generalization.provenance import add_provenance, provenance_model


def test_add_provenance_preserves_one_row_per_design() -> None:
    features = pd.DataFrame(
        {
            "uuid": ["a", "b"],
            "target": ["T1", "T2"],
            "binder_final": [0, 1],
        }
    )
    designs = pd.DataFrame(
        {
            "uuid": ["a", "b"],
            "generator": ["G1", "G2"],
            "design_model": ["M1", "M2"],
            "sequence_design_method": ["S1", "S2"],
            "campaign": ["C1", "C2"],
            "rank": [1, 2],
        }
    )
    result = add_provenance(features, designs)
    assert len(result) == len(features)
    assert result["generator"].tolist() == ["G1", "G2"]


def test_provenance_model_handles_unseen_generator() -> None:
    training = pd.DataFrame(
        {
            "ipsae_mean": [0.1, 0.2, 0.8, 0.9],
            "generator": ["G1", "G1", "G2", "G2"],
        }
    )
    model = provenance_model(["ipsae_mean"], ["generator"])
    model.fit(training, [0, 0, 1, 1])
    probabilities = model.predict_proba(
        pd.DataFrame({"ipsae_mean": [0.5], "generator": ["new_generator"]})
    )
    assert probabilities.shape == (1, 2)
