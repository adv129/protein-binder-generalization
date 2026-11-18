"""Download the compact Overath et al. external-validation table from Zenodo."""

from __future__ import annotations

import hashlib
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "data" / "raw" / "external" / "overath_2025"
FILENAME = "final_dataset.csv"
URL = (
    "https://zenodo.org/api/records/15722219/files/"
    "final_dataset.csv/content"
)
EXPECTED_MD5 = "3a69ee9b0fecf53924a8c6479bac146e"


def md5sum(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    DESTINATION.mkdir(parents=True, exist_ok=True)
    output = DESTINATION / FILENAME
    if output.exists() and md5sum(output) == EXPECTED_MD5:
        print(f"Already downloaded and verified: {output.relative_to(ROOT)}")
        return

    temporary = output.with_suffix(".csv.partial")
    if temporary.exists() and md5sum(temporary) == EXPECTED_MD5:
        temporary.replace(output)
        print(f"Downloaded and verified: {output.relative_to(ROOT)}")
        return

    downloaded_bytes = temporary.stat().st_size if temporary.exists() else 0
    request = urllib.request.Request(URL)
    if downloaded_bytes:
        request.add_header("Range", f"bytes={downloaded_bytes}-")

    with urllib.request.urlopen(request) as response:
        is_resumed = response.status == 206
        mode = "ab" if downloaded_bytes and is_resumed else "wb"
        with temporary.open(mode) as handle:
            shutil.copyfileobj(response, handle)

    actual_md5 = md5sum(temporary)
    if actual_md5 != EXPECTED_MD5:
        raise ValueError(
            f"Checksum mismatch: expected {EXPECTED_MD5}, received {actual_md5}"
        )
    temporary.replace(output)

    print(f"Downloaded and verified: {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
