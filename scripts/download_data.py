"""Download the tabular portion of the public binder-design release."""

from pathlib import Path

from huggingface_hub import hf_hub_download


REPOSITORY = "Anthropic/claude-protein-binder-design"
FILES = (
    "data/tables/design_summary.parquet",
    "data/tables/wetlab/summary.parquet",
    "data/tables/column_dictionary.parquet",
    "data/docs/DATA_NOTES.md",
    "data/docs/WETLAB.md",
    "data/docs/INSILICO.md",
    "data/README.md",
)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    destination = root / "data" / "raw"
    destination.mkdir(parents=True, exist_ok=True)

    for filename in FILES:
        local_path = hf_hub_download(
            repo_id=REPOSITORY,
            filename=filename,
            repo_type="dataset",
            local_dir=destination,
        )
        print(f"Downloaded {Path(local_path).relative_to(root)}")


if __name__ == "__main__":
    main()
