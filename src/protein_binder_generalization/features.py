"""Leakage-safe features available before wet-lab testing."""

from __future__ import annotations

from collections import Counter
from itertools import pairwise
from math import log2

import numpy as np
import pandas as pd

PREDICTORS = (
    "ef2fast",
    "ef2full",
    "ptxv2",
    "odde",
    "afm3",
    "boltz2",
    "chai1",
    "of3",
    "rf3",
    "af3of3",
)

IPSAE_FEATURES = tuple(f"ipsae_min_{predictor}" for predictor in PREDICTORS)

HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

HYDROPHOBIC = frozenset("AVILMFWY")
POSITIVE = frozenset("KR")
NEGATIVE = frozenset("DE")
AROMATIC = frozenset("FWY")

SEQUENCE_FEATURES = (
    "seq_length",
    "seq_mean_hydropathy",
    "seq_hydrophobic_fraction",
    "seq_positive_fraction",
    "seq_negative_fraction",
    "seq_net_charge_proxy",
    "seq_aromatic_fraction",
    "seq_cysteine_fraction",
    "seq_gly_pro_fraction",
    "seq_entropy",
    "seq_max_run_fraction",
)

DISAGREEMENT_FEATURES = (
    "ipsae_mean",
    "ipsae_median",
    "ipsae_std",
    "ipsae_min_across_models",
    "ipsae_max_across_models",
    "ipsae_range",
    "ipsae_models_available",
)

DESIGN_FEATURES = (
    "binder_length",
    "epitope_n_residues",
    "n_optimization_rounds",
)

LEAKAGE_COLUMNS = frozenset(
    {
        "binder_final",
        "kd_nM_final",
        "adaptyv_binding",
        "adaptyv_kd_nM",
        "twist_binding",
        "twist_kd_nM",
        "vendor_agreement",
        "mouse_binding_final",
        "cyno_binding_final",
    }
)


def _fraction(sequence: str, residues: frozenset[str]) -> float:
    return sum(residue in residues for residue in sequence) / len(sequence)


def _max_run(sequence: str) -> int:
    longest = current = 1
    for previous, residue in pairwise(sequence):
        current = current + 1 if residue == previous else 1
        longest = max(longest, current)
    return longest


def sequence_features(sequence: str) -> dict[str, float]:
    """Calculate simple interpretable properties from an amino-acid sequence."""
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("Protein sequence cannot be empty.")
    unknown = set(sequence).difference(HYDROPATHY)
    if unknown:
        raise ValueError(f"Unsupported amino-acid symbols: {sorted(unknown)}")

    counts = Counter(sequence)
    length = len(sequence)
    entropy = -sum(
        (count / length) * log2(count / length) for count in counts.values()
    )
    positive = _fraction(sequence, POSITIVE)
    negative = _fraction(sequence, NEGATIVE)

    return {
        "seq_length": float(length),
        "seq_mean_hydropathy": float(np.mean([HYDROPATHY[aa] for aa in sequence])),
        "seq_hydrophobic_fraction": _fraction(sequence, HYDROPHOBIC),
        "seq_positive_fraction": positive,
        "seq_negative_fraction": negative,
        "seq_net_charge_proxy": positive - negative,
        "seq_aromatic_fraction": _fraction(sequence, AROMATIC),
        "seq_cysteine_fraction": counts.get("C", 0) / length,
        "seq_gly_pro_fraction": (counts.get("G", 0) + counts.get("P", 0)) / length,
        "seq_entropy": entropy,
        "seq_max_run_fraction": _max_run(sequence) / length,
    }


def build_feature_table(designs: pd.DataFrame) -> pd.DataFrame:
    """Return identifiers, outcomes, and features known before experimentation."""
    required = {"uuid", "target", "sequence", "binder_final", *IPSAE_FEATURES}
    missing = required.difference(designs.columns)
    if missing:
        raise KeyError(f"Required columns are absent: {sorted(missing)}")

    sequence_table = pd.DataFrame(
        [sequence_features(sequence) for sequence in designs["sequence"]],
        index=designs.index,
    )
    ipsae = designs.loc[:, IPSAE_FEATURES].astype(float)
    disagreement = pd.DataFrame(
        {
            "ipsae_mean": ipsae.mean(axis=1),
            "ipsae_median": ipsae.median(axis=1),
            "ipsae_std": ipsae.std(axis=1),
            "ipsae_min_across_models": ipsae.min(axis=1),
            "ipsae_max_across_models": ipsae.max(axis=1),
            "ipsae_range": ipsae.max(axis=1) - ipsae.min(axis=1),
            "ipsae_models_available": ipsae.notna().sum(axis=1).astype(float),
        },
        index=designs.index,
    )

    metadata = designs[["uuid", "target", "binder_final"]].copy()
    design_features = designs.loc[:, DESIGN_FEATURES].astype(float)
    features = pd.concat(
        [metadata, ipsae, disagreement, sequence_table, design_features], axis=1
    )
    assert_no_leakage(features.columns.difference(metadata.columns))
    return features


def assert_no_leakage(feature_names: pd.Index | list[str] | tuple[str, ...]) -> None:
    """Fail loudly if a post-experimental variable enters a feature set."""
    leaked = LEAKAGE_COLUMNS.intersection(feature_names)
    if leaked:
        raise ValueError(f"Post-experimental variables found in features: {sorted(leaked)}")


FEATURE_GROUPS = {
    "structure": list(IPSAE_FEATURES + DISAGREEMENT_FEATURES),
    "sequence": list(SEQUENCE_FEATURES + DESIGN_FEATURES),
    "combined": list(
        IPSAE_FEATURES + DISAGREEMENT_FEATURES + SEQUENCE_FEATURES + DESIGN_FEATURES
    ),
}
