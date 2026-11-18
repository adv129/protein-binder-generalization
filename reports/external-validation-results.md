# Independent external validation results

## What was tested

The primary analysis used 1,440 designs from the Anthropic public release. To
test whether its conclusions depend on that one dataset, the same core question
was evaluated on the Overath et al. public meta-analysis:

- 3,669 designs with experimental labels;
- 394 binders across 15 protein targets;
- no exact binder-sequence overlap with the primary dataset;
- independently calculated AF3, Boltz-1, and ColabFold interface-confidence
  scores.

The external predictor panel is not identical to the ten predictors in the
primary release. This is a validation of the scientific hypothesis—whether
interface confidence supports cross-target prioritization—not an identical
feature-by-feature replication.

## Ranking results

Macro-averaged across all 15 external targets:

| Method | Mean AP | Top-10 hit rate |
|---|---:|---:|
| AF3 ipSAE | 0.513 | 0.460 |
| Three-predictor consensus mean | 0.477 | 0.440 |
| Frozen model trained on the primary dataset | 0.477 | 0.407 |
| External leave-one-target-out logistic model | 0.472 | 0.447 |
| Boltz-1 ipSAE | 0.402 | 0.347 |
| ColabFold ipSAE | 0.332 | 0.333 |

AF3 ipSAE has the highest observed mean AP. However, its paired difference from
the frozen primary model is 0.035 with a target-bootstrap interval of
[-0.082, 0.152]. Its difference from the external LOTO model is 0.041 with an
interval of [-0.017, 0.117]. The dataset therefore does not establish that
these top methods differ reliably across targets.

## Cross-dataset transfer

The frozen model was trained once on all labeled Anthropic designs and applied
to the Overath table without retraining or recalibration. It reaches mean AP
0.477 across all 15 external targets.

Four normalized target names occur in both releases: EGFR, TrkA, PD-L1, and
IL-7Ra. On the 11 external targets not seen by name in the primary release, the
frozen model reaches mean AP 0.501. This supports genuine ranking transfer and
is not explained solely by the shared target names. The external designs also
have no exact sequence overlap with the training release.

The frozen model does not outperform direct AF3 ipSAE. The correct conclusion
is that the simple structural-confidence relationship transfers across public
releases, not that the learned model is universally superior.

## Calibration reproduces the primary failure

The external leave-one-target-out model has pooled expected calibration error
of 0.018, which appears excellent. Target-level analysis again contradicts that
pooled impression:

- mean absolute target calibration bias: 0.102;
- Mdm2 observed binder rate: 0.573, mean prediction: 0.305;
- pMHC-SILSY1 observed rate: 0.021, mean prediction: 0.245;
- TrkA observed rate: 0.070, mean prediction: 0.276.

The frozen primary model has pooled ECE 0.078 and mean absolute target bias
0.128. Its mean prediction across all external designs is 0.176, while the
observed binder fraction is 0.107, demonstrating cross-release overprediction.

Thus, the central primary finding independently reproduces: **ranking is useful,
but global probability calibration does not transfer reliably to a particular
new target.**

## Relation to the original external paper

Overath et al. already established that interface-focused metrics, especially
AF3 ipSAE, are strong retrospective filters. This project does not claim that
observation as new. The added contribution here is a consistent cross-release
evaluation emphasizing:

- a frozen model transferred from a separate public design release;
- separation of shared-name and novel-name external targets;
- target-level versus pooled calibration;
- the practical distinction between ranking candidates and interpreting scores
  as probabilities of experimental binding.

## Reproduce

```bash
python scripts/download_external_validation.py
python scripts/audit_external_validation.py
PYTHONPATH=src python scripts/run_external_validation.py
```
