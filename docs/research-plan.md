# Research plan

## Primary question

Do computational confidence metrics reliably rank experimentally successful
de novo protein binders for unseen protein targets, and can they identify when
their predictions are unreliable?

## Hypotheses

1. No single confidence score will perform consistently across all targets.
2. Combining structural, interface, sequence, and disagreement features will
   improve ranking relative to individual scores.
3. Model disagreement and poor calibration will identify targets for which a
   ranking should not be trusted.

## Week-one viability gate

Continue to modeling only if the audit confirms that:

- experimental outcomes can be joined unambiguously to computational records;
- most targets contain both positive and negative designs;
- missingness does not eliminate essential predictors for most designs;
- the label definition can be reproduced from the documentation;
- assay disagreement can be quantified or isolated in a sensitivity analysis.

If these conditions fail, narrow the outcome, select a reliable subset, or
evaluate an alternative public binder benchmark before training models.

## Primary evaluation

Use leave-one-target-out cross-validation. For each fold:

1. Hold out every design for one protein target.
2. Fit preprocessing and the model using only the remaining targets.
3. Rank candidates for the held-out target.
4. Record ranking performance and probability calibration.
5. Repeat until every eligible target has been held out.

Report the distribution across targets in addition to a pooled value.

## Baselines

1. Random candidate ordering.
2. Original design rank.
3. Best individual confidence score.
4. Regularized logistic regression.
5. Gradient-boosted trees.

## Primary outcomes

- Precision-recall area under the curve
- Top-10 and top-20 experimental hit rate
- Enrichment over random selection
- Brier score and calibration plots
- Variation in all metrics across held-out targets

## Required controls

- Designs on which both experimental providers agree
- Alternative treatment of inconclusive measurements
- One representative per backbone or closely related design family
- Models with and without assay or target metadata
- Feature-group ablations

## Claims the project can support

- Whether computational scores enrich for experimental binders in this release
- Whether ranking performance transfers across the included targets
- Whether confidence estimates are calibrated
- Which target or design conditions are associated with failure

## Claims the project cannot support

- Therapeutic efficacy, safety, or clinical usefulness
- Biological function beyond the measured binding assays
- Generalization to every protein family or experimental platform
- The superiority of an AI agent over a matched human expert
