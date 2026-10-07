from __future__ import annotations

import os
import tarfile
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

from .utils import sha256_file


def safe_extract_tar(archive: str | Path, destination: str | Path) -> list[Path]:
    archive = Path(archive)
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    with tarfile.open(archive) as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if destination not in target.parents and target != destination:
                raise ValueError(f"Unsafe path in tar archive: {member.name}")
        handle.extractall(destination, filter="data")
        extracted = [destination / member.name for member in handle.getmembers() if member.isfile()]
    return extracted


def download_file(url: str, destination: str | Path, timeout: int = 120) -> Path:
    import requests

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    temporary = destination.with_suffix(destination.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    os.replace(temporary, destination)
    return destination


def download_manifest_files(manifest: pd.DataFrame, raw_dir: str | Path, include_large: bool = False) -> pd.DataFrame:
    raw_dir = Path(raw_dir)
    records: list[dict[str, object]] = []
    large_accessions = {"GSE207890", "GSE208629"}
    for row in manifest.loc[manifest["enabled"]].itertuples(index=False):
        if row.accession in large_accessions and not include_large:
            records.append({"accession": row.accession, "kind": "processed", "status": "skipped_large", "path": "", "sha256": ""})
            continue
        for kind, url in (("processed", row.processed_url), ("splicing", row.splicing_url)):
            if not url:
                continue
            filename = Path(urlparse(url).path).name
            destination = raw_dir / row.accession / filename
            try:
                path = download_file(url, destination)
                records.append({"accession": row.accession, "kind": kind, "status": "downloaded", "path": str(path), "sha256": sha256_file(path)})
            except Exception as exc:  # retain the full audit instead of hiding partial failures
                records.append({"accession": row.accession, "kind": kind, "status": "failed", "path": str(destination), "sha256": "", "message": str(exc)})
    return pd.DataFrame(records)


def file_inventory(root: str | Path) -> pd.DataFrame:
    root = Path(root)
    records = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        records.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return pd.DataFrame(records)
