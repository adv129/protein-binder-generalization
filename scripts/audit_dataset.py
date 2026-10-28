"""Run the first-pass dataset viability audit."""

from pathlib import Path

import pandas as pd
import yaml


def print_heading(title: str) -> None:
    print(f"\n{title}\n{'=' * len(title)}")


def summarize_missingness(frame: pd.DataFrame, limit: int = 20) -> pd.DataFrame:
    missing = frame.isna().mean().sort_values(ascending=False)
    return missing[missing.gt(0)].head(limit).rename("missing_fraction").to_frame()


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "config" / "default.yaml").read_text())
    design_path = root / config["dataset"]["design_table"]
    wetlab_path = root / config["dataset"]["wetlab_table"]

    if not design_path.exists() or not wetlab_path.exists():
        raise SystemExit("Data not found. Run `python scripts/download_data.py` first.")

    designs = pd.read_parquet(design_path)
    wetlab = pd.read_parquet(wetlab_path)
    columns = config["columns"]

    print_heading("Table dimensions")
    print(f"design_summary: {designs.shape[0]:,} rows × {designs.shape[1]:,} columns")
    print(f"wetlab_summary: {wetlab.shape[0]:,} rows × {wetlab.shape[1]:,} columns")

    print_heading("Identifier checks")
    for name, frame in (("design_summary", designs), ("wetlab_summary", wetlab)):
        identifier = columns["id"]
        if identifier in frame:
            print(
                f"{name}: {frame[identifier].nunique():,} unique IDs; "
                f"{frame[identifier].duplicated().sum():,} duplicate rows"
            )
        else:
            print(f"{name}: expected identifier {identifier!r} is absent")

    target = columns["target"]
    label = columns["label"]
    if {target, label}.issubset(designs.columns):
        print_heading("Outcome by target")
        outcome = (
            designs.groupby(target, dropna=False)[label]
            .agg(designs="size", measured="count", binders="sum", success_rate="mean")
            .sort_values("success_rate", ascending=False)
        )
        print(outcome.to_string(float_format=lambda value: f"{value:.3f}"))
    else:
        print_heading("Outcome by target")
        print(f"Expected columns were not found: {target!r}, {label!r}")

    print_heading("Highest missingness in design summary")
    missingness = summarize_missingness(designs)
    print(missingness.to_string(float_format=lambda value: f"{value:.1%}"))

    print_heading("Columns available for feature discovery")
    print("\n".join(sorted(designs.columns)))


if __name__ == "__main__":
    main()
