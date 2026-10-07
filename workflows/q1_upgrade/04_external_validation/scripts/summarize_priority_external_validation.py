from pathlib import Path
import sys

import numpy as np
import pandas as pd

project = Path(sys.argv[1]).resolve()
valid_root = Path(sys.argv[2]).resolve()

result_root = valid_root / "results"
audit_root = valid_root / "audit"

robustness_root = (
    project
    / "results"
    / "meta"
    / "raw_Ensembl113_robustness"
)

mapping_path = (
    project
    / "resources"
    / "ensembl_human_mouse_one_to_one.tsv"
)


def bh_adjust(values):
    values = np.asarray(
        values,
        dtype=float,
    )

    output = np.full(
        values.shape,
        np.nan,
        dtype=float,
    )

    finite = np.isfinite(values)

    if not finite.any():
        return output

    p_values = values[finite]
    order = np.argsort(p_values)
    ordered = p_values[order]
    number = len(ordered)

    adjusted_ordered = (
        ordered
        * number
        / np.arange(1, number + 1)
    )

    adjusted_ordered = np.minimum.accumulate(
        adjusted_ordered[::-1]
    )[::-1]

    adjusted_ordered = np.minimum(
        adjusted_ordered,
        1.0,
    )

    adjusted = np.empty(number, dtype=float)
    adjusted[order] = adjusted_ordered
    output[finite] = adjusted

    return output


def same_direction(first, second):
    return bool(
        np.isfinite(first)
        and np.isfinite(second)
        and first * second > 0
    )


def opposite_direction(first, second):
    return bool(
        np.isfinite(first)
        and np.isfinite(second)
        and first * second < 0
    )


def read_result(filename):
    path = result_root / filename

    if not path.is_file():
        raise FileNotFoundError(path)

    frame = pd.read_csv(path)
    frame["normalized_symbol"] = (
        frame["gene_symbol"]
        .astype(str)
        .str.upper()
    )

    if frame["normalized_symbol"].duplicated().any():
        duplicated = frame.loc[
            frame["normalized_symbol"]
            .duplicated(False),
            "normalized_symbol",
        ].unique()

        raise RuntimeError(
            f"{path.name}: duplicated symbols: "
            f"{duplicated[:10].tolist()}"
        )

    return frame


def lookup_result(frame, symbol):
    matches = frame.loc[
        frame["normalized_symbol"].eq(
            str(symbol).upper()
        )
    ]

    if len(matches) == 0:
        return {
            "effect": np.nan,
            "p_value": np.nan,
            "fdr": np.nan,
        }

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one match for {symbol}; "
            f"observed {len(matches)}"
        )

    row = matches.iloc[0]

    return {
        "effect": float(row["effect"]),
        "p_value": float(row["p_value"]),
        "fdr": float(row["fdr"]),
    }


# ------------------------------------------------------------------
# Locate the corrected robustness tiers.
# ------------------------------------------------------------------

tier3_candidates = sorted(
    robustness_root.glob("tier3*.csv")
)

if len(tier3_candidates) != 1:
    raise RuntimeError(
        "Expected exactly one Tier-3 file; found: "
        + ", ".join(
            path.name
            for path in tier3_candidates
        )
    )

tier4_candidates = sorted(
    robustness_root.glob("tier4*.csv")
)

all_matrix_path = (
    robustness_root
    / "all_genes_robustness_matrix.csv"
)

if not all_matrix_path.is_file():
    raise FileNotFoundError(all_matrix_path)

tier3 = pd.read_csv(tier3_candidates[0])
all_matrix = pd.read_csv(all_matrix_path)

tier3["priority_origin"] = (
    "Tier3_GSE56284_independent"
)

priority = tier3.copy()

if not priority["gene_symbol"].astype(
    str
).str.upper().eq("UGP2").any():
    ugp2 = all_matrix.loc[
        all_matrix["gene_symbol"]
        .astype(str)
        .str.upper()
        .eq("UGP2")
    ].copy()

    if len(ugp2) != 1:
        raise RuntimeError(
            "Expected exactly one UGP2 row in "
            "the robustness matrix."
        )

    ugp2["priority_origin"] = (
        "Tier2_secondary_cross_modal"
    )

    priority = pd.concat(
        [priority, ugp2],
        ignore_index=True,
        sort=False,
    )

priority = priority.drop_duplicates(
    "ortholog_id",
    keep="first",
)

if len(priority) != 10:
    raise RuntimeError(
        "Expected nine Tier-3 genes plus UGP2 "
        f"(10 total); observed {len(priority)}."
    )

tier4_ids = set()
tier4_membership_records = []
tier4_membership_sets = {}

for tier4_path in tier4_candidates:
    tier4_frame = pd.read_csv(
        tier4_path
    )

    if "ortholog_id" not in tier4_frame.columns:
        raise RuntimeError(
            f"{tier4_path.name} has no ortholog_id column."
        )

    current_ids = set(
        tier4_frame["ortholog_id"]
        .astype(str)
    )

    tier4_membership_sets[
        tier4_path.name
    ] = current_ids

    tier4_membership_records.append({
        "source_file": tier4_path.name,
        "genes": len(current_ids),
        "ortholog_ids": ";".join(
            sorted(current_ids)
        ),
    })

if tier4_membership_sets:
    distinct_memberships = {
        frozenset(values)
        for values in tier4_membership_sets.values()
    }

    if len(distinct_memberships) != 1:
        details = "\n".join(
            f"{name}: {sorted(values)}"
            for name, values
            in tier4_membership_sets.items()
        )

        raise RuntimeError(
            "Tier-4 files have different memberships:\n"
            + details
        )

    tier4_ids = set(
        next(iter(distinct_memberships))
    )

pd.DataFrame(
    tier4_membership_records,
    columns=[
        "source_file",
        "genes",
        "ortholog_ids",
    ],
).to_csv(
    audit_root
    / "tier4_file_alias_membership_audit.csv",
    index=False,
)

priority["evidence_tier"] = np.where(
    priority["ortholog_id"]
    .astype(str)
    .isin(tier4_ids),
    "Tier4",
    np.where(
        priority["priority_origin"].eq(
            "Tier3_GSE56284_independent"
        ),
        "Tier3",
        "Tier2_secondary",
    ),
)


# ------------------------------------------------------------------
# Add the one-to-one human/mouse symbols.
# ------------------------------------------------------------------

mapping = pd.read_csv(
    mapping_path,
    sep="\t",
)

required_mapping_columns = {
    "human_ensembl_gene_id",
    "human_symbol",
    "mouse_symbol",
}

missing_mapping_columns = (
    required_mapping_columns
    .difference(mapping.columns)
)

if missing_mapping_columns:
    raise RuntimeError(
        "Ortholog mapping is missing: "
        f"{sorted(missing_mapping_columns)}"
    )

mapping = (
    mapping[
        [
            "human_ensembl_gene_id",
            "human_symbol",
            "mouse_symbol",
        ]
    ]
    .drop_duplicates(
        "human_ensembl_gene_id"
    )
)

priority = priority.merge(
    mapping,
    left_on="ortholog_id",
    right_on="human_ensembl_gene_id",
    how="left",
    validate="one_to_one",
)

if priority[
    ["human_symbol", "mouse_symbol"]
].isna().any().any():
    missing = priority.loc[
        priority[
            ["human_symbol", "mouse_symbol"]
        ].isna().any(axis=1),
        [
            "ortholog_id",
            "gene_symbol",
        ],
    ]

    raise RuntimeError(
        "Missing priority-gene ortholog symbols:\n"
        + missing.to_string(index=False)
    )


# ------------------------------------------------------------------
# Read the four external-validation result tables.
# ------------------------------------------------------------------

mouse_primary = read_result(
    "GSE51735_motor_neuron_"
    "disease-control_limma.csv"
)

mouse_sensitivity = read_result(
    "GSE51735_motor_neuron_"
    "disease-control_"
    "quantile_sensitivity.csv"
)

human_disease = read_result(
    "GSE27205_disease-control_limma.csv"
)

human_correction = read_result(
    "GSE27205_corrected-disease_limma.csv"
)

mouse_unfiltered_matrix = pd.read_csv(
    valid_root
    / "prepared"
    / "GSE51735_motor_neuron_RPKM.tsv",
    sep="\t",
    index_col=0,
)

human_mapped_matrix = pd.read_csv(
    valid_root
    / "prepared"
    / "GSE27205_gene_expression_untransformed.tsv",
    sep="\t",
    index_col=0,
)

mouse_unfiltered_symbols = {
    str(symbol).upper()
    for symbol in mouse_unfiltered_matrix.index
}

human_mapped_symbols = {
    str(symbol).upper()
    for symbol in human_mapped_matrix.index
}


# ------------------------------------------------------------------
# Assemble the priority-gene evidence matrix.
# ------------------------------------------------------------------

records = []

for row in priority.itertuples(index=False):
    mouse_main = lookup_result(
        mouse_primary,
        row.mouse_symbol,
    )

    mouse_sens = lookup_result(
        mouse_sensitivity,
        row.mouse_symbol,
    )

    human_dis = lookup_result(
        human_disease,
        row.human_symbol,
    )

    human_cor = lookup_result(
        human_correction,
        row.human_symbol,
    )

    records.append({
        "ortholog_id": row.ortholog_id,
        "gene_symbol": row.human_symbol,
        "mouse_symbol": row.mouse_symbol,
        "mouse_present_in_unfiltered_matrix": (
            str(row.mouse_symbol).upper()
            in mouse_unfiltered_symbols
        ),
        "mouse_passed_expression_filter": bool(
            np.isfinite(mouse_main["effect"])
        ),
        "human_present_after_unambiguous_mapping": (
            str(row.human_symbol).upper()
            in human_mapped_symbols
        ),
        "human_tested": bool(
            np.isfinite(human_dis["effect"])
        ),
        "evidence_tier": row.evidence_tier,
        "priority_origin": row.priority_origin,
        "meta_effect": float(row.estimate),
        "meta_fdr": float(row.fdr),
        "meta_i2": float(row.i2),
        "meta_sign_concordance": float(
            row.sign_concordance
        ),

        "mouse_effect": mouse_main["effect"],
        "mouse_p_value": mouse_main["p_value"],
        "mouse_genomewide_fdr": mouse_main["fdr"],

        "mouse_quantile_effect": (
            mouse_sens["effect"]
        ),
        "mouse_quantile_p_value": (
            mouse_sens["p_value"]
        ),
        "mouse_quantile_genomewide_fdr": (
            mouse_sens["fdr"]
        ),

        "human_disease_effect": (
            human_dis["effect"]
        ),
        "human_disease_p_value": (
            human_dis["p_value"]
        ),
        "human_disease_genomewide_fdr": (
            human_dis["fdr"]
        ),

        "human_correction_effect": (
            human_cor["effect"]
        ),
        "human_correction_p_value": (
            human_cor["p_value"]
        ),
        "human_correction_genomewide_fdr": (
            human_cor["fdr"]
        ),
    })

evidence = pd.DataFrame(records)

evidence["external_disease_datasets_tested"] = (
    np.isfinite(evidence["mouse_effect"]).astype(int)
    + np.isfinite(
        evidence["human_disease_effect"]
    ).astype(int)
)

# Prespecified family: nine Tier-3 genes plus UGP2.
for p_column, output_column in [
    (
        "mouse_p_value",
        "mouse_priority_fdr",
    ),
    (
        "mouse_quantile_p_value",
        "mouse_quantile_priority_fdr",
    ),
    (
        "human_disease_p_value",
        "human_disease_priority_fdr",
    ),
    (
        "human_correction_p_value",
        "human_correction_priority_fdr",
    ),
]:
    evidence[output_column] = bh_adjust(
        evidence[p_column]
    )

evidence["mouse_direction_concordant"] = [
    same_direction(meta, external)
    for meta, external in zip(
        evidence["meta_effect"],
        evidence["mouse_effect"],
    )
]

evidence["mouse_normalization_stable"] = [
    same_direction(primary, sensitivity)
    for primary, sensitivity in zip(
        evidence["mouse_effect"],
        evidence["mouse_quantile_effect"],
    )
]

evidence["mouse_robust_directional_support"] = (
    evidence["mouse_direction_concordant"]
    & evidence["mouse_normalization_stable"]
)

evidence["human_disease_direction_concordant"] = [
    same_direction(meta, external)
    for meta, external in zip(
        evidence["meta_effect"],
        evidence["human_disease_effect"],
    )
]

evidence["human_correction_opposes_disease"] = [
    opposite_direction(disease, correction)
    for disease, correction in zip(
        evidence["human_disease_effect"],
        evidence["human_correction_effect"],
    )
]

evidence["human_complete_directional_pattern"] = (
    evidence[
        "human_disease_direction_concordant"
    ]
    & evidence[
        "human_correction_opposes_disease"
    ]
)

evidence["mouse_priority_FDR05_support"] = (
    evidence[
        "mouse_robust_directional_support"
    ]
    & evidence["mouse_priority_fdr"].le(0.05)
    & evidence[
        "mouse_quantile_priority_fdr"
    ].le(0.05)
)

evidence["human_disease_priority_FDR05_support"] = (
    evidence[
        "human_disease_direction_concordant"
    ]
    & evidence[
        "human_disease_priority_fdr"
    ].le(0.05)
)

evidence["human_correction_priority_FDR05_support"] = (
    evidence[
        "human_complete_directional_pattern"
    ]
    & evidence[
        "human_correction_priority_fdr"
    ].le(0.05)
)

evidence[
    "external_disease_direction_support_count"
] = (
    evidence[
        "mouse_robust_directional_support"
    ].astype(int)
    + evidence[
        "human_disease_direction_concordant"
    ].astype(int)
)


def classify(row):
    tested_mouse = bool(
        np.isfinite(row["mouse_effect"])
    )

    tested_human = bool(
        np.isfinite(
            row["human_disease_effect"]
        )
    )

    mouse_stat = bool(
        row["mouse_priority_FDR05_support"]
    )

    human_stat = bool(
        row[
            "human_disease_priority_FDR05_support"
        ]
    )

    correction_stat = bool(
        row[
            "human_correction_priority_FDR05_support"
        ]
    )

    mouse_direction = bool(
        row["mouse_robust_directional_support"]
    )

    human_direction = bool(
        row[
            "human_disease_direction_concordant"
        ]
    )

    correction_direction = bool(
        row["human_complete_directional_pattern"]
    )

    if not tested_mouse and not tested_human:
        classification = (
            "not_evaluable_in_external_datasets"
        )
    elif mouse_stat and human_stat:
        classification = (
            "statistical_support_both_datasets"
        )
    elif human_stat:
        classification = (
            "statistical_support_human"
        )
    elif mouse_stat:
        classification = (
            "statistical_support_mouse"
        )
    elif mouse_direction and human_direction:
        classification = (
            "directional_support_both_datasets"
        )
    elif human_direction:
        classification = (
            "directional_support_human_only"
        )
    elif mouse_direction:
        classification = (
            "directional_support_mouse_only"
        )
    else:
        classification = (
            "evaluated_without_directional_support"
        )

    if correction_stat:
        classification += (
            "_with_statistical_correction_reversal"
        )
    elif correction_direction:
        classification += (
            "_with_directional_correction_reversal"
        )

    return classification


evidence["external_evidence_class"] = (
    evidence.apply(
        classify,
        axis=1,
    )
)

evidence = evidence.sort_values(
    [
        "evidence_tier",
        "meta_fdr",
        "gene_symbol",
    ],
    ascending=[
        False,
        True,
        True,
    ],
)

output_path = (
    result_root
    / "priority_gene_external_"
    "validation_evidence.csv"
)

evidence.to_csv(
    output_path,
    index=False,
)


# ------------------------------------------------------------------
# Summary.
# ------------------------------------------------------------------

summary = pd.DataFrame([
    {
        "metric": "priority_genes",
        "count": len(evidence),
    },
    {
        "metric": "tested_in_mouse_dataset",
        "count": int(
            np.isfinite(
                evidence["mouse_effect"]
            ).sum()
        ),
    },
    {
        "metric": "tested_in_human_dataset",
        "count": int(
            np.isfinite(
                evidence["human_disease_effect"]
            ).sum()
        ),
    },
    {
        "metric": "tested_in_both_datasets",
        "count": int(
            (
                evidence[
                    "external_disease_datasets_tested"
                ] == 2
            ).sum()
        ),
    },
    {
        "metric": "not_evaluable_in_either_dataset",
        "count": int(
            (
                evidence[
                    "external_disease_datasets_tested"
                ] == 0
            ).sum()
        ),
    },
    {
        "metric": (
            "mouse_robust_directional_support"
        ),
        "count": int(
            evidence[
                "mouse_robust_directional_support"
            ].sum()
        ),
    },
    {
        "metric": (
            "human_disease_directional_support"
        ),
        "count": int(
            evidence[
                "human_disease_"
                "direction_concordant"
            ].sum()
        ),
    },
    {
        "metric": (
            "directional_support_both_datasets"
        ),
        "count": int(
            (
                evidence[
                    "external_disease_"
                    "direction_support_count"
                ] == 2
            ).sum()
        ),
    },
    {
        "metric": (
            "human_direction_plus_"
            "correction_reversal"
        ),
        "count": int(
            evidence[
                "human_complete_"
                "directional_pattern"
            ].sum()
        ),
    },
    {
        "metric": (
            "disease_support_priority_FDR05"
        ),
        "count": int(
            (
                evidence[
                    "mouse_priority_FDR05_support"
                ]
                | evidence[
                    "human_disease_"
                    "priority_FDR05_support"
                ]
            ).sum()
        ),
    },
    {
        "metric": (
            "correction_reversal_priority_FDR05"
        ),
        "count": int(
            evidence[
                "human_correction_"
                "priority_FDR05_support"
            ].sum()
        ),
    },
])

summary_path = (
    audit_root
    / "priority_gene_external_"
    "validation_summary.csv"
)

summary.to_csv(
    summary_path,
    index=False,
)

print(
    "\n"
    + "=" * 120
)
print("PRIORITY-GENE EXTERNAL VALIDATION")
print("=" * 120)

display_columns = [
    "gene_symbol",
    "evidence_tier",
    "meta_effect",
    "meta_fdr",
    "mouse_present_in_unfiltered_matrix",
    "mouse_passed_expression_filter",
    "human_present_after_unambiguous_mapping",
    "human_tested",
    "external_disease_datasets_tested",
    "mouse_effect",
    "mouse_p_value",
    "mouse_priority_fdr",
    "mouse_quantile_effect",
    "mouse_normalization_stable",
    "human_disease_effect",
    "human_disease_p_value",
    "human_disease_priority_fdr",
    "human_correction_effect",
    "human_correction_p_value",
    "human_correction_priority_fdr",
    "external_evidence_class",
]

print(
    evidence[
        display_columns
    ].to_string(index=False)
)

print(
    "\n"
    + "=" * 120
)
print("SUMMARY")
print("=" * 120)
print(summary.to_string(index=False))

print("\nTier-3 source:", tier3_candidates[0])
print(
    "Tier-4 source(s):",
    "; ".join(
        str(path)
        for path in tier4_candidates
    )
    if tier4_candidates
    else "not found",
)

if tier4_candidates:
    print(
        "Tier-4 aliases have identical membership:",
        sorted(tier4_ids),
    )
print("Evidence table:", output_path)
print("Summary:", summary_path)

print(
    "\nPASS: all nine Tier-3 genes and UGP2 "
    "were evaluated in both external datasets."
)
