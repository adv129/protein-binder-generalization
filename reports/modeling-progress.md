# Modeling progress

## Question evaluated

Can features learned from other protein targets rank experimentally successful
binders for a completely held-out target?

## Feature matrix

The current leakage-safe feature table contains:

- ten model-specific `ipSAE_min` interface-confidence values;
- their mean, median, standard deviation, minimum, maximum, and range;
- simple sequence properties such as hydropathy, charge proxy, composition, and
  sequence entropy;
- binder length, epitope size, and number of design-optimization rounds.

Target identity, experimental binding calls, K_D measurements, vendor agreement,
and other post-experimental variables are excluded.

## Evaluation design

For each of the 15 targets with experimental labels:

1. Hold out every design for that target.
2. Fit preprocessing and the model using the other 14 targets.
3. Predict probabilities for the held-out target.
4. Calculate target-specific ranking and calibration metrics.

MBP is retained for Brier score and false-confidence analysis but has no positive
examples, so average precision is undefined for its fold.

## Initial results

Across all 14 targets for which average precision is defined, regularized
logistic regression using summary measures of cross-model confidence has the
highest mean average precision (0.580 to 0.581). The raw ten-score logistic
model reaches 0.579, while disagreement alone reaches 0.432.

These macro averages are sensitive to 15-PGDH and BBF-14, which contain only one
and three binders. In the predefined sensitivity analysis restricted to the 12
targets with at least five binders:

- confidence-level summaries: mean AP 0.591;
- confidence plus disagreement summaries: mean AP 0.590;
- ten individual scores: mean AP 0.589;
- disagreement alone: mean AP 0.459;
- confidence plus disagreement summaries: mean top-10 hit rate 0.617.

The near-identical average precision of confidence-level and confidence-plus-
disagreement summaries indicates that disagreement has not yet demonstrated
incremental ranking value. Its small increase in top-10 hit rate requires a
paired uncertainty analysis before it can be interpreted. Sequence features
also do not show a consistent improvement over the strongest structural
summary baseline.

The paired target bootstrap supports that cautious interpretation:

- adding disagreement to confidence-level summaries: mean AP difference -0.001,
  95% interval [-0.015, 0.011];
- confidence-level summaries versus the ten raw scores: difference 0.002,
  interval [-0.043, 0.062];
- adding sequence and design features to the summary baseline: difference
  -0.017, interval [-0.076, 0.036].

None of these comparisons currently establishes a reliable difference in
cross-target average precision.

## Target-specific failures

The structure-only logistic model performs poorly on BBF-14, Cas9, TNF-alpha,
and BHRF1. It performs substantially better on TREM2, Latent GDF-8, PD-L1,
VEGF-A, Nipah-G, and RBX1. This variation is central to the project: aggregate
performance alone does not show when a candidate ranking can be trusted.

## Next analysis

1. Compare model confidence and calibration target by target.
2. Test whether disagreement among structure predictors identifies failures.
3. Add feature-group ablations without selecting models on held-out targets.
4. Quantify uncertainty across targets with bootstrap intervals.
5. Inspect high-confidence false positives and low-confidence true binders.
