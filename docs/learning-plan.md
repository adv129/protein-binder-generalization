# Learning plan

The goal is to understand every scientific and computational decision well
enough to explain it, implement it, and defend it under questioning.

## How each module works

Every module follows the same cycle:

1. **Orient:** place the concept in the complete binder-design workflow.
2. **Predict:** state what you expect before inspecting the result.
3. **Learn:** study the biological, statistical, and computational ideas.
4. **Apply:** implement or analyze one part of the project.
5. **Retrieve:** answer a short closed-notes quiz.
6. **Explain:** teach the concept back in plain language.
7. **Defend:** answer a skeptical research question.
8. **Reflect:** record what changed in your understanding.

A module is complete when you can answer the checkpoint without relying on
memorized terminology alone.

## Conceptual map

```text
Target protein and desired binding surface
                  |
                  v
Backbone generation (for example, RFdiffusion)
                  |
                  v
Sequence design (for example, ProteinMPNN)
                  |
                  v
Complex prediction and scoring (AlphaFold-like models, Boltz, others)
                  |
                  v
Candidate filtering and experimental selection
                  |
                  v
Expression and binding assays (SPR/BLI, affinity, quality control)
                  |
                  v
Our project: test which computational evidence predicts wet-lab success
                  |
                  v
Decision: which candidates should be tested for a new target?
```

Our project does not treat a high-confidence predicted structure as proof of
binding. It measures how much that evidence changes the probability of an
experimental success and whether the relationship transfers to another target.

## Module 1: From binder design to experimental evidence

### Learn

- What a protein binder is
- The difference between backbone generation, sequence design, and structure
  prediction
- What RFdiffusion, ProteinMPNN, and complex predictors contribute
- Why computational plausibility is not experimental validation
- Where this project sits in the workflow

### Apply

Draw the workflow from memory and annotate the input and output of each model.

### Checkpoint

- Explain why RFdiffusion and AlphaFold solve different problems.
- Explain why a predicted complex cannot establish experimental binding.
- Identify the precise decision our project is intended to improve.

## Module 2: Binding experiments and label construction

### Learn

- Binding affinity and the dissociation constant, K_D
- Association and dissociation rates
- SPR and BLI at a conceptual level
- Expression failure, nonspecific binding, and inconclusive measurements
- Agreement and disagreement between experimental providers

### Apply

Read the wet-lab documentation and write an explicit decision rule for positive,
negative, excluded, and uncertain outcomes.

### Checkpoint

- Explain why lower K_D usually represents tighter binding.
- Distinguish a nonbinder from an unmeasurable design.
- Defend the chosen label rule against a skeptical experimentalist.

## Module 3: Dataset provenance and scientific auditing

### Learn

- Unit of observation
- Data provenance and table joins
- Missingness, duplicates, repeated backbones, and selection bias
- The distinction between an externally measured result and independent
  replication
- Why a recent technical release requires careful claim boundaries

### Apply

Run the viability audit, construct a data dictionary, and produce an outcome
table by target.

### Checkpoint

- State what one row represents in every table used.
- Identify three ways leakage or duplication could inflate performance.
- Explain which conclusions the dataset can and cannot support.

## Module 4: Generalization and data leakage

### Learn

- Training, validation, and test data
- Grouped splitting
- Leave-one-target-out evaluation
- Distribution shift between protein targets
- Why random row splitting answers the wrong scientific question

### Apply

Implement the target-level folds and verify that no target occurs in both the
training and test portions of a fold.

### Checkpoint

- Compare the questions answered by random splitting and target-level splitting.
- Predict which evaluation will produce the higher score and explain why.
- Defend the split to a reviewer concerned about target similarity.

## Module 5: Ranking metrics and calibration

### Learn

- Class imbalance
- Precision, recall, and the precision-recall curve
- Average precision
- Top-k hit rate and enrichment over random selection
- Probability calibration and the Brier score
- Ranking quality versus trustworthy probabilities

### Apply

Evaluate random ranking, original candidate rank, and individual confidence
scores for every held-out target.

### Checkpoint

- Calculate precision and recall from a small confusion matrix.
- Explain why ROC AUC may look favorable on an imbalanced dataset.
- Distinguish a model that ranks well from one that is calibrated.
- Explain top-k performance in terms of a laboratory testing budget.

## Module 6: Interpretable predictive modeling

### Learn

- Logistic regression and regularization
- Tree-based models
- Preprocessing inside a cross-validation fold
- Overfitting with small grouped datasets
- Feature importance, coefficients, and their limitations

### Apply

Train regularized logistic regression and gradient-boosted tree baselines using
identical target-level folds.

### Checkpoint

- Interpret a logistic-regression coefficient without claiming causality.
- Explain why preprocessing must be fitted using training data only.
- Defend the use of small models instead of a neural network.

## Module 7: Structural features and model disagreement

### Learn

- Predicted aligned error and complex confidence
- Interface contacts, clashes, and buried surface area
- Variation across seeds and structural models
- Epistemic uncertainty and out-of-distribution targets
- Why confidence values from different models may not be directly comparable

### Apply

Compare feature groups and test whether prediction disagreement is associated
with experimental failure or poor calibration.

### Checkpoint

- Explain what an interface-confidence score measures and what it does not.
- Give two reasons that different structure predictors might agree and still be
  wrong.
- Explain how model disagreement could be converted into a decision rule.

## Module 8: Scientific interpretation and defense

### Learn

- Association versus causation
- Sensitivity analyses and ablations
- Target-level uncertainty
- Reproducibility and claim boundaries
- Translating a metric into an experimental decision

### Apply

Write the results as a decision-oriented report and prepare a five-minute oral
defense.

### Checkpoint

- State the strongest supported conclusion in one sentence.
- Identify the most serious remaining alternative explanation.
- Describe the next experiment that would reduce the largest uncertainty.
- Explain why a successful binder is not automatically a therapeutic candidate.

## Testing standard

Quizzes will mix four types of questions:

- **Recall:** define a concept without notes.
- **Transfer:** apply it to a new protein or dataset.
- **Diagnosis:** identify an error in an analysis plan.
- **Defense:** respond to a difficult question from a professor or reviewer.

If an answer is incomplete, revisit the concept immediately and answer a new
transfer question. Scores are less important than being able to reason from the
workflow and evidence.
