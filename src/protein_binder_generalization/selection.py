"""Decision-focused evaluation under fixed experimental budgets."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

RANDOM_STATE = 398


def candidate_selection_table(
    predictions: pd.DataFrame,
    budgets: Sequence[int] = (1, 3, 5, 10, 20, 30),
) -> pd.DataFrame:
    """Calculate target-level yield for top-ranked experimental candidates."""
    required = {"target", "binder_final", "predicted_probability"}
    missing = required.difference(predictions.columns)
    if missing:
        raise KeyError(f"Required prediction columns are absent: {sorted(missing)}")
    if not budgets or any(budget <= 0 for budget in budgets):
        raise ValueError("Experimental budgets must be positive integers.")

    rows = []
    for target, group in predictions.groupby("target"):
        ranked = group.sort_values("predicted_probability", ascending=False)
        n_designs = len(ranked)
        n_binders = int(ranked["binder_final"].sum())
        prevalence = n_binders / n_designs
        for requested_budget in budgets:
            budget = min(requested_budget, n_designs)
            selected = ranked.head(budget)
            hits = int(selected["binder_final"].sum())
            oracle_hits = min(budget, n_binders)
            random_expected_hits = budget * prevalence
            rows.append(
                {
                    "target": target,
                    "requested_budget": requested_budget,
                    "budget": budget,
                    "n_designs": n_designs,
                    "n_binders": n_binders,
                    "target_prevalence": prevalence,
                    "model_hits": hits,
                    "model_hit_rate": hits / budget,
                    "random_expected_hits": random_expected_hits,
                    "oracle_hits": oracle_hits,
                    "hit_gain_over_random": hits - random_expected_hits,
                    "recall_of_all_binders": (
                        hits / n_binders if n_binders > 0 else np.nan
                    ),
                    "fraction_of_oracle": (
                        hits / oracle_hits if oracle_hits > 0 else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def summarize_candidate_selection(selection: pd.DataFrame) -> pd.DataFrame:
    """Summarize yield for all measured targets and a minimum-binder cohort."""
    cohorts = {
        "all_measured_targets": selection,
        "targets_with_at_least_5_binders": selection[
            selection["n_binders"].ge(5)
        ],
    }
    summaries = []
    for cohort_name, cohort in cohorts.items():
        summary = (
            cohort.groupby("requested_budget", as_index=False)
            .agg(
                targets=("target", "nunique"),
                mean_model_hits=("model_hits", "mean"),
                median_model_hits=("model_hits", "median"),
                minimum_model_hits=("model_hits", "min"),
                q25_model_hits=("model_hits", lambda values: values.quantile(0.25)),
                q75_model_hits=("model_hits", lambda values: values.quantile(0.75)),
                mean_random_expected_hits=("random_expected_hits", "mean"),
                mean_oracle_hits=("oracle_hits", "mean"),
                mean_model_hit_rate=("model_hit_rate", "mean"),
                mean_fraction_of_oracle=("fraction_of_oracle", "mean"),
            )
        )
        summary.insert(0, "cohort", cohort_name)
        summaries.append(summary)
    return pd.concat(summaries, ignore_index=True)


def bootstrap_gain_over_random(
    selection: pd.DataFrame,
    minimum_binders: int | None = None,
    n_bootstrap: int = 20_000,
) -> pd.DataFrame:
    """Bootstrap mean target-level hit gains over random selection."""
    if minimum_binders is not None:
        selection = selection[selection["n_binders"].ge(minimum_binders)]
    rng = np.random.default_rng(RANDOM_STATE)
    rows = []
    for budget, group in selection.groupby("requested_budget"):
        differences = group["hit_gain_over_random"].to_numpy(dtype=float)
        sampled_indices = rng.integers(
            0, len(differences), size=(n_bootstrap, len(differences))
        )
        bootstrap_means = differences[sampled_indices].mean(axis=1)
        rows.append(
            {
                "requested_budget": budget,
                "minimum_binders": minimum_binders,
                "targets": len(differences),
                "mean_hit_gain_over_random": differences.mean(),
                "ci_2_5": np.quantile(bootstrap_means, 0.025),
                "ci_97_5": np.quantile(bootstrap_means, 0.975),
                "bootstrap_probability_gain_positive": (
                    bootstrap_means > 0
                ).mean(),
            }
        )
    return pd.DataFrame(rows)
