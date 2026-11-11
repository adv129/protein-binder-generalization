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
logistic regression using structural confidence and disagreement features has
the highest mean average precision (0.557). The combined structural and sequence
model has a lower mean average precision (0.518) but a slightly better mean
Brier score and top-10 hit rate.

These macro averages are sensitive to 15-PGDH and BBF-14, which contain only one
and three binders. In the predefined sensitivity analysis restricted to the 12
targets with at least five binders:

- combined logistic regression: mean AP 0.573;
- structure-only logistic regression: mean AP 0.564;
- combined logistic regression: mean top-10 hit rate 0.600;
- structure-only logistic regression: mean top-10 hit rate 0.567.

The paired target bootstrap estimates a combined-minus-structure mean AP
difference of 0.009, with a 95% interval from approximately -0.043 to 0.059.
The present evidence therefore does not establish that sequence features improve
cross-target ranking, although they may help on particular targets and improve
some decision metrics.

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
