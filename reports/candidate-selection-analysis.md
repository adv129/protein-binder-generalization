# Candidate selection under fixed experimental budgets

## Decision question

If a lab can test only a limited number of designs for a previously unseen
protein target, how many experimentally confirmed binders does the model recover?

This analysis ranks each held-out target using predictions from the structural
confidence logistic model. Every ranking is produced by a model trained without
any designs from that target. The top 1, 3, 5, 10, 20, and 30 candidates are
then compared with two retrospective references:

- random-selection expectation: budget multiplied by the target's observed
  binder prevalence;
- oracle ceiling: the maximum possible hits if the experimental labels were
  known before selection.

The random expectation uses the observed prevalence only for retrospective
evaluation. That prevalence would not be known for a genuinely new target.

## Results across all measured targets

All 15 targets with experimental measurements are included here, including MBP
with no confirmed binders and the two targets with fewer than five binders.

| Candidates tested | Mean model hits | Random expectation | Oracle ceiling | Median hits | Middle 50% |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.73 | 0.27 | 0.93 | 1.0 | 0.5–1.0 |
| 3 | 1.73 | 0.80 | 2.67 | 2.0 | 1.0–3.0 |
| 5 | 2.80 | 1.33 | 4.27 | 3.0 | 1.0–4.5 |
| 10 | 4.87 | 2.65 | 8.27 | 5.0 | 1.0–8.0 |
| 20 | 8.60 | 5.31 | 13.93 | 10.0 | 2.5–13.0 |
| 30 | 12.27 | 7.96 | 17.33 | 12.0 | 5.0–18.5 |

At a budget of ten experiments, ranking gains 2.21 binders per target over the
random-selection expectation. A paired target bootstrap gives an interval of
[1.12, 3.35]. The average gain is therefore consistent across resampled sets of
targets, but it does not guarantee success on any individual target.

## Sensitivity analysis

For the 12 targets with at least five confirmed binders, the top ten contain an
average of 6.00 binders versus 3.26 expected from random selection. The paired
gain is 2.74 binders, with a bootstrap interval of [1.48, 3.95].

The target-level range remains large:

- TREM2: 10 of the top 10 are binders;
- PD-L1: 9 of 10;
- Latent GDF-8, Nipah-G, RBX1, and VEGF-A: 8 of 10;
- BHRF1 and Cas9: only 1 of 10.

Thus, the model improves expected experimental efficiency while still failing
badly on specific protein targets.

## Interpretation

The strongest supported practical claim is about **candidate prioritization**,
not absolute probability prediction:

> On average across held-out targets, structural-confidence ranking identifies
> more binders per experiment than random selection, but target identity causes
> large and currently unpredictable variation in yield.

The shaded region in the budget figure is the middle 50% of target outcomes,
not a confidence interval around the mean. It visualizes deployment variability.
The bootstrap intervals quantify uncertainty in the mean paired gain across
targets.

## Reproduce

```bash
PYTHONPATH=src python scripts/analyze_candidate_selection.py
```

The script writes target-level selection outcomes, cohort summaries, bootstrap
results, and `reports/figures/candidate_selection_budget_curve.png`.
