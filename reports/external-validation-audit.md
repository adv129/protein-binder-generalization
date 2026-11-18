# External validation dataset audit

## Dataset accepted for validation

The Overath et al. public release is suitable for independent validation because
it contains experimentally tested successes and failures, not only selected
successful structures. The downloaded `final_dataset.csv` matches the MD5
checksum published by Zenodo.

After excluding seven rows without an experimental label, the usable table has:

- 3,669 labeled designs;
- 394 binders and 3,275 nonbinders;
- 15 protein targets;
- 20 source campaigns;
- 312 total columns, including independently calculated AF2, AF3, Boltz-1,
  ColabFold, Rosetta, and interface-confidence features.

There are no duplicate binder identifiers among labeled rows and no exact binder
sequence overlap with the 1,440-design Anthropic dataset used for the primary
analysis.

## Important imbalance

FGFR2 contributes 2,123 of the 3,669 labeled designs and 193 of the 394 binders.
A pooled random split would therefore be dominated by one target and would not
measure transfer to a new protein. Macro-averaged leave-one-target-out results
remain necessary.

Observed target success rates range from 0.021 for pMHC-SILSY1 to 0.573 for
Mdm2. This large variation provides an independent test of the target-level
calibration problem found in the primary dataset.

## Missingness and comparable scores

The three primary interface-confidence scores have low missingness among labeled
designs:

- AF3 `ipSAE_min`: 0.0% after removal of the seven unlabeled rows;
- Boltz-1 `ipSAE_min`: 0.1%;
- ColabFold `ipSAE_min`: 0.5%.

The release does not contain the same ten predictor-specific ipSAE columns as
the Anthropic dataset. External validation will therefore test the same
scientific hypothesis—cross-target ranking by multi-predictor interface
confidence—using the independently released AF3, Boltz-1, and ColabFold scores.
It will not claim to be an identical feature-by-feature replication.

## Source and target structure

Each source campaign contains only one target, although some targets occur in
multiple source campaigns. Consequently, source/laboratory effects are partly
confounded with target identity. Holding out the entire target is still safer
than a random split, but the analysis cannot cleanly estimate an independent
laboratory effect.

## Modeling decision

The first external evaluation will include:

1. direct per-target ranking by AF3, Boltz-1, and ColabFold `ipSAE_min`;
2. simple mean and median consensus rankings across the three predictors;
3. a fixed logistic model trained with leave-one-target-out validation;
4. target-level calibration and fixed-budget yield summaries.

All primary metrics will be macro-averaged across targets so FGFR2 cannot
dominate the conclusion.
