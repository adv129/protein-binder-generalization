# Protein Binder Generalization

Can computational confidence metrics rank experimentally successful de novo
protein binders for a previously unseen protein target?

This project evaluates whether structure confidence, interface measurements,
sequence properties, and model disagreement transfer across protein targets.
The initial benchmark uses Anthropic's open protein-binder design release,
which links computational predictions to measurements from two contract
research organizations.

## Research goals

1. Audit the experimental labels and quantify agreement between assays.
2. Establish random-ranking and single-score baselines.
3. Evaluate models with leave-one-target-out cross-validation.
4. Measure ranking quality, calibration, and target-to-target variability.
5. Identify conditions in which computational rankings should not be trusted.

The project treats the Anthropic release as a benchmark dataset, not as an
assumption that the accompanying scientific claims are correct.

## Quick start

Requires Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python scripts/download_data.py
python scripts/audit_dataset.py
```

The downloader retrieves only the tabular data and documentation required for
the initial audit. It does not download the full structural archive.

## Repository structure

```text
config/                 Analysis configuration
data/                   Local data directories; downloaded data are ignored
docs/                   Research plan and decisions
notebooks/              Exploratory notebooks
reports/figures/        Generated figures
scripts/                Reproducible command-line workflows
src/                    Reusable Python package
```

## Evaluation rule

All primary results will use leave-one-target-out evaluation: models train on
all but one protein target and are evaluated on the held-out target. Randomly
splitting individual designs is not an acceptable primary evaluation because
closely related designs for the same target could appear in both sets.

## Dataset

- [Dataset card](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design)
- [Technical report](https://systemcards.org/docs/anthropic-claude-opus-4-8-other)

The dataset is licensed CC BY 4.0. Confirm its citation requirements before
publishing derived results.

## Current milestone

Complete the dataset viability audit described in
[`docs/research-plan.md`](docs/research-plan.md) before fitting predictive
models.
