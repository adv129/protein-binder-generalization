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

## Design-provenance confounder check

The public dataset contains candidates from ten generation pipelines with
different pooled experimental success rates. A separate leave-one-target-out
ablation therefore tested whether the structural model was merely exploiting
generator identity or other design provenance.

Among the 12 targets with at least five binders:

- structural confidence summaries: mean AP 0.591;
- all provenance variables without structure confidence: mean AP 0.505;
- generator identity alone: mean AP 0.361;
- structural confidence plus all provenance: mean AP 0.549.

Adding all provenance variables changes paired target AP by -0.042, with a
bootstrap interval of [-0.084, -0.006]. Generator identity therefore does not
explain the structural model's ranking performance and does not improve this
fixed cross-target model. This is a predictive ablation, not a causal ranking
of design generators; see `reports/design-provenance-analysis.md`.

## Target-specific failures

The structure-only logistic model performs poorly on BBF-14, Cas9, TNF-alpha,
and BHRF1. It performs substantially better on TREM2, Latent GDF-8, PD-L1,
VEGF-A, Nipah-G, and RBX1. This variation is central to the project: aggregate
performance alone does not show when a candidate ranking can be trusted.

## Fixed experimental budgets

Decision-focused evaluation now measures how many binders are recovered when a
lab tests only the top 1, 3, 5, 10, 20, or 30 candidates for each held-out
target. Across all 15 measured targets, the top ten contain an average of 4.87
binders, compared with 2.65 expected from random selection. The paired gain is
2.21 binders per target, with a bootstrap interval of [1.12, 3.35].

Among the 12 targets with at least five binders, the top-ten average is 6.00
versus a random expectation of 3.26. However, individual targets range from 1
to 10 recovered binders. The model therefore improves average experimental
yield without providing a reliable per-target guarantee. Full results are in
`reports/candidate-selection-analysis.md`.

## Independent external validation

The main conclusion was tested on the Overath et al. public meta-analysis after
excluding seven unlabeled rows: 3,669 designs, 394 binders, and 15 targets. No
external binder sequence exactly matches a sequence in the primary dataset.

Direct AF3 ipSAE obtains mean target AP 0.513. A model trained once on the
Anthropic release and applied unchanged to the external table obtains mean AP
0.477 across all targets and 0.501 across the 11 targets not seen by normalized
name in training. The paired target bootstrap does not establish that AF3,
external LOTO logistic regression, and the frozen primary model differ reliably
in ranking quality.

Calibration failure independently reproduces. External LOTO predictions have
pooled ECE 0.018 but mean absolute target bias 0.102. The frozen cross-dataset
model has pooled ECE 0.078 and mean absolute target bias 0.128. Full details are
in `reports/external-validation-results.md`.

## Next analysis

1. Determine whether a model can abstain on targets where rankings are unsafe.
2. Test whether candidate diversification across predictors or generators
   improves the worst-target yield.
3. Consolidate the primary and external results into the final class report and
   presentation, clearly separating reproduced findings from new analyses.
