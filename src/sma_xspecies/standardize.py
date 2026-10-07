from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm

from .utils import bh_adjust, strip_ensembl_version


ALIASES = {
    "gene_id": ["gene_id", "test_id", "GeneID", "ensembl_gene_id", "id"],
    "gene_symbol": ["gene_symbol", "gene", "gene_name", "symbol", "GeneSymbol"],
    "effect": ["effect", "log2(fold_change)", "log2FoldChange", "logFC", "estimate"],
    "se": ["se", "lfcSE", "std_error", "standard_error"],
    "p_value": ["p_value", "pvalue", "P.Value", "PValue", "p"],
    "fdr": ["fdr", "q_value", "padj", "adj.P.Val", "FDR", "qvalue"],
    "z": ["z", "stat", "t"],
}


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    lower = path.name.lower()
    if lower.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    compression = "gzip" if lower.endswith(".gz") else "infer"
    try:
        return pd.read_csv(path, sep=None, engine="python", compression=compression)
    except Exception:
        return pd.read_csv(path, sep="\t", compression=compression)


def _find_column(columns: Iterable[str], canonical: str, required: bool = False) -> str | None:
    lookup = {str(column).lower(): str(column) for column in columns}
    for alias in ALIASES[canonical]:
        if alias.lower() in lookup:
            return lookup[alias.lower()]
    if required:
        raise ValueError(f"Could not identify required {canonical} column; found {list(columns)}")
    return None


def standardize_de_frame(
    frame: pd.DataFrame,
    *,
    accession: str,
    species: str,
    tissue: str,
    role: str,
    eligibility: str,
    contrast: str = "disease-control",
) -> pd.DataFrame:
    gene_id_col = _find_column(frame.columns, "gene_id", required=True)
    symbol_col = _find_column(frame.columns, "gene_symbol")
    effect_col = _find_column(frame.columns, "effect", required=True)
    se_col = _find_column(frame.columns, "se")
    p_col = _find_column(frame.columns, "p_value")
    fdr_col = _find_column(frame.columns, "fdr")
    z_col = _find_column(frame.columns, "z")

    out = pd.DataFrame({
        "gene_id": strip_ensembl_version(frame[gene_id_col]),
        "gene_symbol": frame[symbol_col].astype("string") if symbol_col else pd.Series(pd.NA, index=frame.index, dtype="string"),
        "effect": pd.to_numeric(frame[effect_col], errors="coerce"),
    })
    out["se"] = pd.to_numeric(frame[se_col], errors="coerce") if se_col else np.nan
    out["p_value"] = pd.to_numeric(frame[p_col], errors="coerce") if p_col else np.nan
    out["fdr"] = pd.to_numeric(frame[fdr_col], errors="coerce") if fdr_col else bh_adjust(out["p_value"])
    if z_col:
        out["z"] = pd.to_numeric(frame[z_col], errors="coerce")
    else:
        out["z"] = out["effect"] / out["se"]
        missing_z = ~np.isfinite(out["z"]) & np.isfinite(out["p_value"]) & np.isfinite(out["effect"])
        clipped = out.loc[missing_z, "p_value"].clip(lower=np.finfo(float).tiny, upper=1.0)
        out.loc[missing_z, "z"] = np.sign(out.loc[missing_z, "effect"]) * norm.isf(clipped / 2.0)

    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["gene_id", "effect"])
    out = out.loc[out["gene_id"].str.len().fillna(0) > 0].copy()
    out.insert(0, "accession", accession)
    out["species"] = species
    out["tissue"] = tissue
    out["role"] = role
    out["eligibility"] = eligibility
    out["contrast"] = contrast
    out["effect_unit"] = "log2_fold_change"
    return out.drop_duplicates(["accession", "contrast", "gene_id"], keep="first")


def standardize_published_summary(path: str | Path, manifest_row: pd.Series) -> pd.DataFrame:
    return standardize_de_frame(
        read_table(path),
        accession=manifest_row["accession"],
        species=manifest_row["species"],
        tissue=manifest_row["tissue"],
        role=manifest_row["role"],
        eligibility=manifest_row["eligibility"],
    )


def prepare_expression_matrix(frame: pd.DataFrame, sample_annotations: pd.DataFrame) -> pd.DataFrame:
    """Normalize a processed expression table to gene_id/gene_symbol/GSM columns.

    This function does not calculate DE. It prepares a matrix for the blocked limma
    engine and fails when the sample mapping is ambiguous.
    """
    gene_id_col = _find_column(frame.columns, "gene_id", required=True)
    symbol_col = _find_column(frame.columns, "gene_symbol")
    rename: dict[str, str] = {}
    missing = []
    for row in sample_annotations.itertuples(index=False):
        candidates = [row.sample_id]
        if hasattr(row, "sample_name") and pd.notna(row.sample_name):
            candidates.append(str(row.sample_name))
        match = next((column for column in frame.columns if str(column) in candidates), None)
        if match is None:
            missing.append(row.sample_id)
        else:
            rename[str(match)] = row.sample_id
    if missing:
        raise ValueError(f"Processed matrix is missing annotated samples: {missing}")
    out = pd.DataFrame({"gene_id": strip_ensembl_version(frame[gene_id_col])})
    out["gene_symbol"] = frame[symbol_col].astype("string") if symbol_col else pd.Series(pd.NA, index=frame.index, dtype="string")
    for source, target in rename.items():
        out[target] = pd.to_numeric(frame[source], errors="coerce")
    out = out.dropna(subset=["gene_id"])
    numeric = list(rename.values())
    collapsed = out.groupby("gene_id", as_index=False).agg(
        {"gene_symbol": "first", **{column: "mean" for column in numeric}}
    )
    if collapsed[numeric].notna().mean().min() < 0.8:
        raise ValueError("More than 20% of expression values are nonnumeric or missing")
    return collapsed
