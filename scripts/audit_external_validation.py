"""Audit the Overath et al. table before external validation modeling."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = (
    ROOT
    / "data"
    / "raw"
    / "external"
    / "overath_2025"
    / "final_dataset.csv"
)
ANTHROPIC = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
TABLES = ROOT / "reports" / "tables"

VALIDATION_FEATURES = (
    "af3_ipSAE_min",
    "boltz1_ipSAE_min",
    "colab_ipSAE_min",
    "af3_ipSAE_avg",
    "boltz1_ipSAE_avg",
    "colab_ipSAE_avg",
    "af3_iptm_avg",
    "boltz1_iptm_avg",
    "colab_iptm_avg",
    "af3_rosetta_interface_sc",
    "af3_rosetta_interface_dG",
    "af3_rosetta_interface_dSASA",
    "af2_binder_aligned_rmsd",
)


def normalized_sequences(values: pd.Series) -> set[str]:
    return set(values.dropna().str.strip().str.upper())


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    designs = pd.read_csv(EXTERNAL, low_memory=False)
    required = {"binder_id", "target_id", "binder", "source", "A_seq"}
    missing = required.difference(designs.columns)
    if missing:
        raise KeyError(f"Required external columns are absent: {sorted(missing)}")

    labeled = designs[designs["binder"].notna()].copy()
    labeled["binder"] = labeled["binder"].astype(bool).astype(int)
    duplicate_labeled_ids = int(labeled["binder_id"].duplicated().sum())
    if duplicate_labeled_ids:
        raise ValueError("Labeled external binder IDs are not unique.")

    target_summary = (
        labeled.groupby("target_id", as_index=False)
        .agg(
            designs=("binder", "size"),
            binders=("binder", "sum"),
            success_rate=("binder", "mean"),
            source_campaigns=("source", "nunique"),
        )
        .sort_values("success_rate", ascending=False)
    )
    source_summary = (
        labeled.groupby("source", as_index=False)
        .agg(
            designs=("binder", "size"),
            targets=("target_id", "nunique"),
            binders=("binder", "sum"),
            success_rate=("binder", "mean"),
        )
        .sort_values("designs", ascending=False)
    )
    feature_missingness = pd.DataFrame(
        {
            "feature": VALIDATION_FEATURES,
            "missing_fraction": [
                labeled[feature].isna().mean() for feature in VALIDATION_FEATURES
            ],
        }
    ).sort_values("missing_fraction", ascending=False)

    external_sequences = normalized_sequences(labeled["A_seq"])
    sequence_overlap = None
    if ANTHROPIC.exists():
        anthropic = pd.read_parquet(ANTHROPIC)
        anthropic_sequences = normalized_sequences(anthropic["sequence"])
        sequence_overlap = len(external_sequences.intersection(anthropic_sequences))

    target_summary.to_csv(TABLES / "external_target_outcomes.csv", index=False)
    source_summary.to_csv(TABLES / "external_source_outcomes.csv", index=False)
    feature_missingness.to_csv(
        TABLES / "external_feature_missingness.csv", index=False
    )

    print(f"Raw table: {len(designs):,} rows x {len(designs.columns):,} columns")
    print(
        f"Labeled table: {len(labeled):,} designs, "
        f"{labeled['binder'].sum():,} binders"
    )
    print(
        f"Targets: {labeled['target_id'].nunique()}; "
        f"source campaigns: {labeled['source'].nunique()}"
    )
    print(f"Rows without experimental labels: {designs['binder'].isna().sum()}")
    print(f"Duplicate labeled binder IDs: {duplicate_labeled_ids}")
    if sequence_overlap is not None:
        print(f"Exact sequence overlap with primary dataset: {sequence_overlap}")
    print("\nOutcome by target")
    print(target_summary.round(3).to_string(index=False))
    print("\nValidation-feature missingness")
    print(feature_missingness.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
