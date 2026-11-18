# Design-provenance confounder analysis

## Question

Do the cross-target predictions work because structural confidence generalizes,
or because the model can exploit which public design pipeline generated a
candidate?

Generator identity is a plausible confounder because the ten generators have
different overall experimental success rates. For example, the observed binder
rate is 0.423 for FreeBindCraft designs and 0.217 for RFdiffusion designs.
Those pooled rates also mix together different protein targets, so they cannot
by themselves establish that one generator is better.

## Evaluation

The comparison uses the same leave-one-target-out folds and fixed regularized
logistic regression as the primary analysis. All preprocessing is fit on the
training targets only. The evaluated feature groups are:

1. four structural confidence summaries: mean, median, minimum, and maximum
   ipSAE across the ten structure predictors;
2. generator identity alone;
3. all available pre-experimental provenance: generator, design model,
   sequence-design method, campaign, and submitted rank;
4. structural confidence plus generator identity;
5. structural confidence plus all provenance variables.

The primary sensitivity analysis includes the 12 targets with at least five
confirmed binders. Categorical variables are one-hot encoded with unknown
categories allowed at inference time.

## Results

| Feature group | Mean AP | Mean top-10 hit rate |
|---|---:|---:|
| Structural confidence summaries | 0.591 | 0.600 |
| Confidence + all provenance | 0.549 | 0.575 |
| Confidence + generator | 0.534 | 0.558 |
| All design provenance | 0.505 | 0.533 |
| Generator identity only | 0.361 | 0.350 |

Relative to structural confidence alone, adding all provenance variables
changes mean AP by -0.042. A paired target bootstrap gives an interval of
[-0.084, -0.006]. Generator identity alone is substantially worse than the
structural baseline: difference -0.230, interval [-0.317, -0.145].

## Interpretation

The structural-confidence result is not explained by a simple generator-identity
shortcut. Provenance variables contain some signal, but they rank candidates on
unseen targets less effectively than the four structural confidence summaries.
Adding provenance to the fixed structural model does not improve ranking and
slightly reduces average performance.

This is an ablation result, not a causal comparison of design generators. The
dataset was not randomized across generators and targets, generator sample sizes
are unequal, and the fixed logistic model may not represent every useful
interaction. The correct claim is therefore that provenance does not improve
this predefined cross-target model—not that generator choice never matters.

## Reproduce

```bash
PYTHONPATH=src python scripts/analyze_design_provenance.py
```

The script writes fold metrics, out-of-fold predictions, summary tables,
paired bootstrap results, and `reports/figures/design_provenance_comparison.png`.
