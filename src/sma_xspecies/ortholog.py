from __future__ import annotations

from pathlib import Path

import pandas as pd

from .utils import strip_ensembl_version


REQUIRED = {
    "human_ensembl_gene_id", "human_symbol", "mouse_ensembl_gene_id",
    "mouse_symbol", "orthology_type", "ensembl_release",
}


def read_one_to_one_mapping(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing = REQUIRED.difference(frame.columns)
    if missing:
        raise ValueError(f"Ortholog mapping is missing columns: {sorted(missing)}")
    frame = frame.loc[frame["orthology_type"] == "ortholog_one2one"].copy()
    frame["human_ensembl_gene_id"] = strip_ensembl_version(frame["human_ensembl_gene_id"])
    frame["mouse_ensembl_gene_id"] = strip_ensembl_version(frame["mouse_ensembl_gene_id"])
    if frame["ensembl_release"].nunique() != 1:
        raise ValueError("Ortholog mapping must contain exactly one Ensembl release")
    if frame.duplicated("human_ensembl_gene_id").any() or frame.duplicated("mouse_ensembl_gene_id").any():
        raise ValueError("Primary mapping must be one-to-one after filtering")
    return frame


def attach_ortholog_ids(effects: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    def attach_species(source: pd.DataFrame, id_column: str, symbol_column: str) -> pd.DataFrame:
        source = source.copy().reset_index(drop=True)
        source["_row_id"] = source.index
        by_id = source.merge(mapping, left_on="gene_id", right_on=id_column, how="inner")
        by_id["mapping_key"] = "ensembl_gene_id"
        unmatched = source.loc[~source["_row_id"].isin(by_id["_row_id"])].copy()
        # A symbol fallback is allowed only as a lookup into the same versioned,
        # one-to-one Compara table and only when that species symbol is unique.
        unique_symbols = mapping.loc[
            mapping[symbol_column].notna() & ~mapping[symbol_column].duplicated(keep=False)
        ]
        by_symbol = unmatched.merge(unique_symbols, left_on="gene_symbol", right_on=symbol_column, how="inner")
        by_symbol["mapping_key"] = "versioned_unique_symbol"
        return pd.concat([by_id, by_symbol], ignore_index=True, sort=False)

    human = attach_species(
        effects.loc[effects["species"] == "Homo_sapiens"], "human_ensembl_gene_id", "human_symbol"
    )
    mouse = attach_species(
        effects.loc[effects["species"] == "Mus_musculus"], "mouse_ensembl_gene_id", "mouse_symbol"
    )
    combined = pd.concat([human, mouse], ignore_index=True, sort=False)
    combined["ortholog_id"] = combined["human_ensembl_gene_id"]
    combined["ortholog_symbol"] = combined["human_symbol"].fillna(combined["ortholog_id"])
    return combined.drop(columns=["_row_id"], errors="ignore")
