from pathlib import Path
import sys

import numpy as np
import pandas as pd

from sma_xspecies.meta import (
    read_effect_files,
    collapse_within_study_orthologs,
)
from sma_xspecies.ortholog import (
    read_one_to_one_mapping,
    attach_ortholog_ids,
)

project = Path(sys.argv[1]).resolve()
output_root = Path(sys.argv[2]).resolve()
output_root.mkdir(parents=True, exist_ok=True)

target_accessions = [
    "GSE108094",
    "GSE209926",
    "GSE56284",
]

effects = read_effect_files(
    project / "results" / "de"
)

neural = effects.loc[
    effects["accession"].isin(target_accessions)
    & effects["contrast"].eq("disease-control")
    & effects["tissue"].isin([
        "motor_neuron",
        "spinal_cord",
    ])
    & effects["eligibility"].eq("eligible")
].copy()

observed = sorted(
    neural["accession"].unique()
)

if observed != sorted(target_accessions):
    raise RuntimeError(
        "Unexpected neural accessions: "
        f"{observed}"
    )

mapping = read_one_to_one_mapping(
    project
    / "resources"
    / "ensembl_human_mouse_one_to_one.tsv"
)

mapped = attach_ortholog_ids(
    neural,
    mapping,
)

required_mapped = {
    "accession",
    "ortholog_id",
    "ortholog_symbol",
    "species",
    "tissue",
    "effect",
    "se",
    "z",
}

missing = required_mapped.difference(
    mapped.columns
)

if missing:
    raise RuntimeError(
        "Mapped effects missing columns: "
        f"{sorted(missing)}"
    )

audit_records = []

for accession in target_accessions:
    input_subset = neural.loc[
        neural["accession"].eq(accession)
    ]

    mapped_subset = mapped.loc[
        mapped["accession"].eq(accession)
    ]

    audit_records.append({
        "accession": accession,
        "input_rows": len(input_subset),
        "mapped_rows": int(
            mapped_subset["ortholog_id"]
            .notna()
            .sum()
        ),
        "mapping_fraction": float(
            mapped_subset["ortholog_id"]
            .notna()
            .mean()
        ),
    })

collapsed = collapse_within_study_orthologs(
    mapped
)

collapsed["gene_symbol"] = (
    collapsed["ortholog_symbol"]
    .astype("string")
    .fillna(collapsed["ortholog_id"])
)

for column in [
    "effect",
    "se",
    "z",
]:
    collapsed[column] = pd.to_numeric(
        collapsed[column],
        errors="coerce",
    )

valid = collapsed.loc[
    collapsed["ortholog_id"].notna()
    & np.isfinite(collapsed["effect"])
    & np.isfinite(collapsed["se"])
    & (collapsed["se"] > 0)
].copy()

coverage = (
    valid.groupby(
        "ortholog_id",
        as_index=False,
    )
    .agg(
        gene_symbol=(
            "gene_symbol",
            "first",
        ),
        k=(
            "accession",
            "nunique",
        ),
        n_species=(
            "species",
            "nunique",
        ),
        datasets=(
            "accession",
            lambda values: ";".join(
                sorted(set(values))
            ),
        ),
        species=(
            "species",
            lambda values: ";".join(
                sorted(set(values))
            ),
        ),
    )
)

eligible = coverage.loc[
    (coverage["k"] >= 2)
    & (coverage["n_species"] >= 2)
].copy()

analysis = valid.loc[
    valid["ortholog_id"].isin(
        eligible["ortholog_id"]
    )
].copy()

analysis = analysis[[
    "ortholog_id",
    "gene_symbol",
    "accession",
    "species",
    "tissue",
    "effect",
    "se",
    "z",
    "n_source_rows",
]].sort_values(
    [
        "ortholog_id",
        "accession",
    ]
)

existing_meta = pd.read_csv(
    project
    / "results"
    / "meta"
    / "gene_meta_neural_joint.csv"
)

existing_ids = set(
    existing_meta["ortholog_id"]
    .dropna()
    .astype(str)
)

prepared_ids = set(
    eligible["ortholog_id"]
    .dropna()
    .astype(str)
)

missing_from_prepared = sorted(
    existing_ids - prepared_ids
)

extra_in_prepared = sorted(
    prepared_ids - existing_ids
)

comparison = pd.DataFrame({
    "category": [
        "existing_meta_ids",
        "prepared_meta_ids",
        "missing_from_prepared",
        "extra_in_prepared",
    ],
    "count": [
        len(existing_ids),
        len(prepared_ids),
        len(missing_from_prepared),
        len(extra_in_prepared),
    ],
})

audit = pd.DataFrame(audit_records)

collapsed_counts = (
    analysis.groupby(
        "accession"
    )["ortholog_id"]
    .nunique()
    .rename(
        "analysis_orthologs"
    )
    .reset_index()
)

audit = audit.merge(
    collapsed_counts,
    on="accession",
    how="left",
)

analysis.to_csv(
    output_root
    / "study_level_neural_effects.csv",
    index=False,
)

coverage.to_csv(
    output_root
    / "ortholog_study_coverage.csv",
    index=False,
)

audit.to_csv(
    output_root
    / "study_effect_preparation_audit.csv",
    index=False,
)

comparison.to_csv(
    output_root
    / "prepared_vs_existing_meta_ids.csv",
    index=False,
)

pd.DataFrame({
    "missing_from_prepared": pd.Series(
        missing_from_prepared,
        dtype="string",
    ),
}).to_csv(
    output_root
    / "orthologs_missing_from_prepared.csv",
    index=False,
)

pd.DataFrame({
    "extra_in_prepared": pd.Series(
        extra_in_prepared,
        dtype="string",
    ),
}).to_csv(
    output_root
    / "orthologs_extra_in_prepared.csv",
    index=False,
)

print("=" * 90)
print("STUDY-EFFECT PREPARATION AUDIT")
print("=" * 90)
print(audit.to_string(index=False))

print()
print("=" * 90)
print("ORTHOLOG COVERAGE")
print("=" * 90)

print(
    eligible["k"]
    .value_counts()
    .sort_index()
    .rename_axis("k")
    .rename("orthologs")
    .to_string()
)

print()
print(comparison.to_string(index=False))

print()
print(
    "Study-effect rows:",
    len(analysis),
)

print(
    "Eligible cross-species orthologs:",
    len(eligible),
)

if missing_from_prepared or extra_in_prepared:
    raise RuntimeError(
        "Prepared ortholog universe does not match "
        "the existing neural meta-analysis. "
        "Review the comparison files before fitting."
    )

print()
print(
    "PASS: study-level effects exactly match "
    "the existing neural meta-analysis universe."
)
