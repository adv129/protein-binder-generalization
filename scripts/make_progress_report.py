"""Create the dataset-audit figures and preliminary single-score baselines."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import average_precision_score


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "data" / "tables" / "design_summary.parquet"
FIGURES = ROOT / "reports" / "figures"
TABLES = ROOT / "reports" / "tables"


def target_summary(data: pd.DataFrame) -> pd.DataFrame:
    """Summarize label availability and binder rate for every target."""
    return (
        data.groupby("target", dropna=False)["binder_final"]
        .agg(designs="size", measured="count", binders="sum", success_rate="mean")
        .sort_values("success_rate", ascending=True, na_position="first")
    )


def plot_target_success(summary: pd.DataFrame) -> None:
    """Plot the large target-to-target variation in experimental success."""
    measured = summary.dropna(subset=["success_rate"]).copy()
    sns.set_theme(style="whitegrid")
    figure, axis = plt.subplots(figsize=(10, 6.4))
    bars = axis.barh(
        measured.index,
        measured["success_rate"] * 100,
        color="#4C78A8",
        edgecolor="white",
    )
    axis.set_xlabel("Experimentally confirmed binders (%)")
    axis.set_ylabel("")
    axis.set_xlim(0, 88)
    axis.set_title("Binder success varies strongly across protein targets", loc="left")
    axis.spines[["top", "right", "left"]].set_visible(False)
    for bar, (_, row) in zip(bars, measured.iterrows(), strict=True):
        axis.text(
            bar.get_width() + 1,
            bar.get_y() + bar.get_height() / 2,
            f"{int(row['binders'])}/{int(row['measured'])}",
            va="center",
            fontsize=9,
        )
    figure.tight_layout()
    figure.savefig(FIGURES / "target_success_rates.png", dpi=200, bbox_inches="tight")
    plt.close(figure)


def single_score_baselines(data: pd.DataFrame) -> pd.DataFrame:
    """Evaluate each precomputed score separately within each eligible target."""
    score_columns = sorted(column for column in data if column.startswith("ipsae_min_"))
    scores = {column: 1 for column in score_columns}
    scores["original_rank"] = -1

    records: list[dict[str, object]] = []
    labeled = data[data["binder_final"].notna()].copy()
    labeled["binder_final"] = labeled["binder_final"].astype(int)

    for target, target_data in labeled.groupby("target"):
        if target_data["binder_final"].nunique() < 2:
            continue
        prevalence = target_data["binder_final"].mean()
        for score_name, direction in scores.items():
            source_column = "rank" if score_name == "original_rank" else score_name
            subset = target_data[[source_column, "binder_final"]].dropna()
            if subset["binder_final"].nunique() < 2:
                continue
            score = subset[source_column] * direction
            order = score.sort_values(ascending=False).index
            top = subset.loc[order].head(10)
            top_rate = top["binder_final"].mean()
            records.append(
                {
                    "target": target,
                    "score": score_name,
                    "n": len(subset),
                    "positives": int(subset["binder_final"].sum()),
                    "prevalence": prevalence,
                    "average_precision": average_precision_score(
                        subset["binder_final"], score
                    ),
                    "top_10_hit_rate": top_rate,
                    "top_10_enrichment": top_rate / prevalence,
                }
            )
    return pd.DataFrame(records)


def summarize_baselines(results: pd.DataFrame) -> pd.DataFrame:
    """Aggregate single-score results without pooling designs across targets."""
    summary = (
        results.groupby("score")
        .agg(
            targets_evaluated=("target", "nunique"),
            mean_average_precision=("average_precision", "mean"),
            median_average_precision=("average_precision", "median"),
            mean_top_10_hit_rate=("top_10_hit_rate", "mean"),
            mean_top_10_enrichment=("top_10_enrichment", "mean"),
        )
        .sort_values("mean_average_precision", ascending=False)
    )
    return summary


def plot_baselines(summary: pd.DataFrame, random_baseline: float) -> None:
    """Plot preliminary target-wise average precision for each single score."""
    plot_data = summary.sort_values("mean_average_precision", ascending=True)
    figure, axis = plt.subplots(figsize=(10, 6.4))
    axis.barh(
        plot_data.index.str.replace("ipsae_min_", "", regex=False),
        plot_data["mean_average_precision"],
        color="#72A0C1",
        edgecolor="white",
    )
    axis.set_xlabel("Mean average precision across eligible targets")
    axis.set_ylabel("Single computational score")
    axis.set_xlim(0, max(0.65, plot_data["mean_average_precision"].max() + 0.05))
    axis.axvline(
        random_baseline,
        color="#555555",
        linestyle="--",
        linewidth=1.5,
        label=f"Mean random baseline ({random_baseline:.2f})",
    )
    axis.legend(frameon=False, loc="lower right")
    axis.set_title("Preliminary target-wise ranking performance", loc="left")
    axis.spines[["top", "right", "left"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(
        FIGURES / "single_score_average_precision.png", dpi=200, bbox_inches="tight"
    )
    plt.close(figure)


def write_feature_provenance(data: pd.DataFrame) -> None:
    """Record which column groups can be used at candidate-selection time."""
    rows = [
        ("binder_final", "Wet-lab outcome", "Outcome", False),
        ("kd_nM_final", "Wet-lab affinity", "Outcome characterization", False),
        ("adaptyv_* / twist_*", "Wet-lab assays", "Label construction and QC", False),
        ("vendor_agreement", "Wet-lab comparison", "Sensitivity analysis", False),
        ("ipsae_min_*", "Structure prediction", "Primary predictor", True),
        ("sc_dockq_*", "Structure prediction", "Secondary predictor", True),
        ("sequence", "Designed candidate", "Derived sequence features", True),
        ("binder_length", "Designed candidate", "Predictor", True),
        ("rank", "Candidate selection", "Baseline", True),
        ("generator", "Design provenance", "Metadata/ablation", True),
        ("campaign / design_model", "Design provenance", "Metadata/ablation", True),
    ]
    provenance = pd.DataFrame(
        rows,
        columns=["column_group", "created_by", "analysis_role", "available_before_lab"],
    )
    provenance["present_in_release"] = provenance["column_group"].map(
        lambda value: value in data.columns or "*" in value or "/" in value
    )
    provenance.to_csv(TABLES / "feature_provenance.csv", index=False)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    data = pd.read_parquet(DATA)

    targets = target_summary(data)
    targets.to_csv(TABLES / "target_outcomes.csv")
    plot_target_success(targets)

    baseline_results = single_score_baselines(data)
    baseline_results.to_csv(TABLES / "single_score_by_target.csv", index=False)
    baseline_summary = summarize_baselines(baseline_results)
    baseline_summary.to_csv(TABLES / "single_score_summary.csv")
    eligible_prevalence = (
        baseline_results.loc[
            baseline_results["score"].eq("ipsae_min_odde"), "prevalence"
        ].mean()
    )
    plot_baselines(baseline_summary, eligible_prevalence)

    write_feature_provenance(data)
    print(baseline_summary.round(3).to_string())
    print(f"\nWrote figures to {FIGURES.relative_to(ROOT)}/")
    print(f"Wrote tables to {TABLES.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
