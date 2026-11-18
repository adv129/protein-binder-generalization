# External validation plan

## Selected public dataset

The independent validation dataset is the Overath et al. meta-analysis,
**Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of
3,766 Experimentally Characterised Binders**.

- Paper: <https://doi.org/10.1101/2025.08.14.670059>
- Dataset: <https://doi.org/10.5281/zenodo.15722219>
- License: CC BY 4.0 for the dataset
- Scale: 3,766 experimentally tested designs across 15 targets
- Labels: binary experimental binding outcome
- Features: more than 200 released structural, energetic, and confidence
  measurements from AF2, AF3, Boltz-1, ColabFold, and Rosetta

Only the 82 MB `final_dataset.csv` table is required initially. The multi-GB
prediction archives will not be downloaded unless a later analysis requires
raw structures or confidence matrices.

## Why this dataset

This release has both binders and nonbinders rather than only publication-selected
successes. It was assembled independently from the Anthropic campaigns, includes
multiple source studies, and provides target identifiers suitable for grouped
evaluation. It is therefore a meaningful test of whether the current conclusions
generalize beyond one design release.

## Validation questions

1. Do interface-confidence scores rank experimental binders above nonbinders on
   independently collected targets?
2. Does apparent pooled calibration again hide large target-level probability
   errors?
3. Do simple models outperform single-score rankings when evaluated with
   source- and target-held-out splits?
4. Does fixed-budget selection improve expected wet-lab yield while retaining
   large worst-target failures?

## Leakage and comparability checks

Before modeling, the audit will:

- identify target and source-study overlap;
- confirm label definitions and missingness;
- distinguish features calculated for every design from source-specific fields;
- avoid treating designs from the same source/target campaign as independent;
- report results separately when an exact feature match to the Anthropic table
  is unavailable.

The first analysis will reproduce single-score and grouped cross-validation
baselines using only released table columns. Re-running AF3 or Boltz on thousands
of structures is outside the initial scope and is unnecessary for this audit.

## Reproduce the download

```bash
python scripts/download_external_validation.py
```

The script verifies the MD5 checksum published in the Zenodo record. Downloaded
data remain under the ignored `data/raw/` directory and are not committed.
