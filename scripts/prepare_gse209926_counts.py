#!/usr/bin/env python3
"""Create a gene-by-sample count matrix from the GSE209926 XLS archive."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sma_xspecies.io import safe_extract_tar


GENE_ALIASES = ("gene_id", "geneid", "ensembl_gene_id", "id", "gene")
COUNT_ALIASES = ("count", "counts", "read_count", "readcount", "raw_count", "raw counts")


def choose_column(frame: pd.DataFrame, aliases: tuple[str, ...], kind: str) -> str:
    lookup = {str(column).strip().lower(): str(column) for column in frame.columns}
    for alias in aliases:
        if alias in lookup:
            return lookup[alias]
    if kind == "count":
        numeric = [column for column in frame if pd.to_numeric(frame[column], errors="coerce").notna().mean() > 0.9]
        if len(numeric) == 1:
            return numeric[0]
    raise ValueError(f"Could not unambiguously identify {kind} column in {list(frame.columns)}")


def read_excel_table(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path)
    if frame.shape[1] < 2:
        frame = pd.read_excel(path, header=1)
    return frame


def build_matrix(archive: Path, sample_map_path: Path, work_dir: Path) -> pd.DataFrame:
    files = safe_extract_tar(archive, work_dir)
    sample_map = pd.read_csv(sample_map_path)
    columns = []
    for sample_id in sample_map["sample_id"]:
        matches = [path for path in files if sample_id in path.name and path.suffix.lower() in {".xls", ".xlsx"}]
        if len(matches) != 1:
            raise ValueError(f"Expected one Excel file for {sample_id}; found {len(matches)}")
        frame = read_excel_table(matches[0])
        gene_column = choose_column(frame, GENE_ALIASES, "gene")
        count_column = choose_column(frame, COUNT_ALIASES, "count")
        values = pd.to_numeric(frame[count_column], errors="coerce")
        if (values.dropna() < 0).any() or not np.allclose(values.dropna(), np.round(values.dropna())):
            raise ValueError(f"{matches[0].name} does not contain non-negative integer counts")
        columns.append(pd.Series(values.to_numpy(), index=frame[gene_column].astype(str), name=sample_id))
    matrix = pd.concat(columns, axis=1, join="outer").fillna(0)
    matrix.index.name = "gene_id"
    return matrix.reset_index()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("sample_map", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--work-dir", type=Path, default=Path("data/interim/GSE209926"))
    args = parser.parse_args()
    result = build_matrix(args.archive, args.sample_map, args.work_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {result.shape[0]:,} genes × {result.shape[1] - 1} samples to {args.output}")


if __name__ == "__main__":
    main()

