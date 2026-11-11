# Calibration and failure-case analysis

## Model analyzed

The analysis uses out-of-fold predictions from regularized logistic regression
with four cross-model confidence-level features: mean, median, minimum, and
maximum ipSAE. Every design is predicted by a model that was trained without any
designs from that target.

## Aggregate calibration hides target-level failure

The pooled equal-frequency reliability analysis produces an expected calibration
error of 0.034. Considered alone, that would suggest good calibration.

Target-level analysis gives a different conclusion:

- mean absolute target calibration bias: 0.164;
- median absolute target calibration bias: 0.122;
- largest absolute target calibration bias: 0.439.

Large positive and negative errors cancel when designs from all targets are
pooled. The low pooled error therefore does not mean that predicted probabilities
are trustworthy for a particular unseen target.

## Largest systematic errors

| Target | Observed binder rate | Mean prediction | Bias |
|---|---:|---:|---:|
| VEGF-A | 0.600 | 0.161 | -0.439 |
| TREM2 | 0.800 | 0.365 | -0.435 |
| BHRF1 | 0.256 | 0.508 | +0.253 |
| TrkA | 0.222 | 0.445 | +0.223 |
| MBP | 0.000 | 0.192 | +0.192 |

The model substantially underestimates success on VEGF-A and TREM2 and
overestimates it on BHRF1, TrkA, and MBP.

## Candidate-selection behavior

The top ten ranked designs contain:

- 10 binders for TREM2;
- 9 for PD-L1;
- 8 each for VEGF-A, RBX1, Nipah-G, and Latent GDF-8;
- 1 each for Cas9, BHRF1, and 15-PGDH;
- 0 for BBF-14 and MBP.

BHRF1 is a particularly important failure. Its top-10 hit rate is 0.10 even
though its overall binder rate is 0.256, so the model ranking is worse than
selecting an average BHRF1 design.

## What the individual failures show

The three highest-confidence nonbinders are BHRF1 designs with predicted
probabilities from 0.690 to 0.706. Their mean ipSAE values are approximately
0.90, and the standard deviation across predictors is only 0.05. Multiple
structure predictors therefore agree confidently on designs that fail the
experiment.

Several successful binders show the opposite pattern. Low-confidence binders
from Cas9, VEGF-A, EGFR, BBF-14, and other targets have a minimum ipSAE of zero
and high cross-model variation. A single or small subset of predictors can fail
on a design that still binds experimentally.

These examples explain why model disagreement is not a universal uncertainty
measure in this dataset:

- agreement can be confidently wrong;
- disagreement can occur for genuine binders;
- predictor failure may depend systematically on the protein target.

## Current conclusion

The model is useful for ranking candidates on several targets, but its absolute
probabilities should not be interpreted as target-independent probabilities of
experimental binding. Target-level distribution shift is a larger problem than
the pooled calibration curve suggests.

## Next questions

1. Which target properties are associated with systematic overprediction or
   underprediction?
2. Are particular structure predictors responsible for low-confidence true
   binders on Cas9 and other difficult targets?
3. Can a target-shift warning be produced without using experimental labels from
   the new target?
4. Does a conservative candidate-selection rule improve worst-target behavior
   while preserving top-10 yield?
