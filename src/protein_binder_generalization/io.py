"""Dataset loading helpers."""

from pathlib import Path

import pandas as pd


def read_parquet(path: str | Path) -> pd.DataFrame:
    """Read a required Parquet table with a helpful missing-file error."""
    table_path = Path(path)
    if not table_path.exists():
        raise FileNotFoundError(
            f"Missing {table_path}. Run `python scripts/download_data.py` first."
        )
    return pd.read_parquet(table_path)
